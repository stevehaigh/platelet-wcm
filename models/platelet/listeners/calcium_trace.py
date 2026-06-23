"""
CalciumTrace listener for the platelet whole-cell model.

Records the Ca²⁺ signalling state each timestep for post-simulation analysis
and validation against the Dolan & Diamond (2014) Fig. 4 Ca²⁺ transient.

Columns written:
  time              — simulation time (s)
  ca_cyt_nM         — free cytosolic Ca²⁺ (nM)
  ca_dts_uM         — DTS stored Ca²⁺ (µM)
  ip3_nM            — IP₃ concentration (nM)
  soce_flux_nMs     — instantaneous SOCE influx rate into cytosol (nM/s)
  stim1_dim         — STIM1 dimer particle count (Dolan Table S1 convention)
  cam_free          — free calmodulin count
  ca2_cam           — Ca₂·CaM count (N-lobe loaded)
  ca4_cam           — Ca₄·CaM count (fully loaded; activates PMCA)
  ca4_cam_pmca      — Ca₄·CaM·PMCA count (CaM-activated, empty)
  ca4_cam_pmca_ca   — Ca₄·CaM·PMCA·Ca count (CaM-activated, loaded)
  pmca_cam          — PMCA·CaM count (deactivating)
  pmca_free         — free PMCA count
  pmca_ca           — PMCA·Ca count (basal active)
  atp_pump_per_s    — ATP consumed by Ca²⁺ pumps (SERCA + PMCA) this step
                      (molecules/s; dt = 1 s)
  pkc_active        — active (DAG + Ca²⁺-bound) PKC count (v0.6 feedback)
  p2y1_desensitised_frac — fraction of the P2Y1 pool in the PKC-phosphorylated
                      desensitised state (0–1; v0.6 feedback)
  plcb_phosphorylated_frac — fraction of the PLCβ pool phosphorylated out of the
                      Gq-activatable pool by PKC (0–1; v0.6 Slice 3, Purvis route)
  ca_mito_count     — mitochondrial matrix Ca²⁺ (raw count; mito volume is not
                      tracked separately, so this is ions, not a concentration)
  mcu_uptake_per_s  — MCU uptake flux into the matrix this step (ions/s),
                      recomputed from the ODE's Hill kinetics so the otherwise
                      unobserved mitochondrial module is auditable (issue #76)
  mito_coupling_factor — #76 Part 2 evoked IP3R-release gate factor (0–1) =
                      ∝ MCU capacity × Ca²⁺-activation; 1.0 in WT, < 1 during the
                      KO transient (SOCE is not gated by it — it falls indirectly)
"""

import numpy as np

import math

import os

import wholecell.listeners.listener
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	_IDX,
	_UM_PER_COUNT_CYT,
	_UM_PER_COUNT_DTS,
	_mwc_open_fraction,
	pka_active_frac,
	GAMMA_SOC_S,
	ip3r_relief_factor,
	K_MITO,
	NA_OVER_zF,
	ORAI_SUBUNITS_PER_CHANNEL,
	PUNCTA,
	RT_OVER_zF_V,
	V_PM_V,
)
# Extracellular Ca²⁺ comes from the per-run RunConfig (sim.run_config); the
# SOCE current is physically zero under the EDTA condition (ca_ex = 0).


class CalciumTrace(wholecell.listeners.listener.Listener):
	"""Record Ca²⁺ dynamics state each timestep."""

	_name = 'CalciumTrace'

	def __init__(self, *args, **kwargs):
		self._bulk_molecules = None
		self._live_file = None
		self._live_path = None
		super().__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super().initialize(sim, sim_data)

		# Per-run conditions from RunConfig (not mutated globals): extracellular
		# Ca²⁺ for the SOCE trace, and the live-CSV toggle.
		self._config = sim.run_config
		self._ca_ex_uM = self._config.ca_ex_mM * 1000.0
		# MCU knockout scale (1.0 intact, 0.0 KO); copied out like _ca_ex_uM so
		# the uptake-flux recompute below stays unit-testable without a config.
		self._mcu_vmax_scale = self._config.mcu_vmax_scale
		# #76 Part 2 — MCU → (SOCE + release) coupling gain (copied like the
		# others so the SOCE-flux recompute stays unit-testable without a config).
		self._mito_coupling_gain = self._config.mito_coupling_gain

		self._bulk_molecules = (
			sim.internal_states['BulkMolecules'].container)

		# CalciumDynamics exposes the per-step ATP cost of SERCA + PMCA
		# pumping as `_atp_cost` (set in calculateRequest, applied in
		# evolveState). Listeners update after evolveState, so by the time
		# `update()` runs it holds the value for the step just completed.
		self._calcium_dynamics = sim.processes['CalciumDynamics']

		# Pre-compute the global indices of the species we track.
		all_ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
		self._idx_ca_cyt         = all_ids.index('CA2_CYT[c]')
		self._idx_ca_dts         = all_ids.index('CA2_DTS[dts]')
		self._idx_ip3            = all_ids.index('IP3[c]')
		self._idx_stim1          = all_ids.index('STIM1_dim[dts]')
		self._idx_orai           = all_ids.index('ORAI1[pl]')
		self._idx_cam_free       = all_ids.index('CaM_free[c]')
		self._idx_ca2_cam        = all_ids.index('Ca2_CaM[c]')
		self._idx_ca4_cam        = all_ids.index('Ca4_CaM[c]')
		self._idx_ca4_cam_pmca   = all_ids.index('Ca4_CaM_PMCA[pl]')
		self._idx_ca4_cam_pmca_ca = all_ids.index('Ca4_CaM_PMCA_Ca[pl]')
		self._idx_pmca_cam       = all_ids.index('PMCA_CaM[pl]')
		self._idx_pmca_free      = all_ids.index('PMCA[pl]')
		self._idx_pmca_ca        = all_ids.index('PMCA_Ca[pl]')
		# Mitochondrial matrix Ca²⁺ (issue #76 observability).
		self._idx_ca_mito        = all_ids.index('CA2_MITO[m]')
		# PKC negative-feedback states (v0.6).
		self._idx_pkc_active     = all_ids.index('PKC_active[c]')
		self._idx_p2y1_inactive  = all_ids.index('P2Y1_inactive[pl]')
		self._idx_p2y1_active    = all_ids.index('P2Y1_active[pl]')
		self._idx_p2y1_desens    = all_ids.index('P2Y1_desensitised[pl]')
		# PKC → PLCβ phosphorylation (v0.6 Slice 3).
		self._idx_plcb_inactive  = all_ids.index('PLCb_inactive[c]')
		self._idx_plcb_active    = all_ids.index('PLCb_active[c]')
		self._idx_plcb_phos      = all_ids.index('PLCb_phosphorylated[c]')
		# Inhibitory axis (v0.7 Slice 2, #10): P2Y12 / cAMP / PKA / VASP.
		self._idx_camp           = all_ids.index('cAMP[c]')
		self._idx_p2y12_inactive = all_ids.index('P2Y12_inactive[pl]')
		self._idx_p2y12_active   = all_ids.index('P2Y12_active[pl]')
		self._idx_vasp           = all_ids.index('VASP[c]')
		self._idx_vasp_phos      = all_ids.index('VASP_phos[c]')

		# Initialise logged quantities.
		self.ca_cyt_nM        = 0.0
		self.ca_dts_uM        = 0.0
		self.ip3_nM           = 0.0
		self.soce_flux_nMs    = 0.0
		self.stim1_dim        = 0
		self.cam_free         = 0
		self.ca2_cam          = 0
		self.ca4_cam          = 0
		self.ca4_cam_pmca     = 0
		self.ca4_cam_pmca_ca  = 0
		self.pmca_cam         = 0
		self.pmca_free        = 0
		self.pmca_ca          = 0
		self.ca_mito_count    = 0
		self.mcu_uptake_per_s = 0.0
		self.mito_coupling_factor = 1.0
		self.atp_pump_per_s   = 0.0
		self.pkc_active             = 0
		self.p2y1_desensitised_frac = 0.0
		self.plcb_phosphorylated_frac = 0.0
		self.camp_uM                = 0.0
		self.pka_frac               = 0.0
		self.p2y12_active_frac      = 0.0
		self.vasp_phos_frac         = 0.0

		self.registerLoggedQuantity('Ca²⁺_cyt\n(nM)',   'ca_cyt_nM',     '.1f')
		self.registerLoggedQuantity('Ca²⁺_dts\n(µM)',   'ca_dts_uM',     '.1f')
		self.registerLoggedQuantity('Ca²⁺_mito\n(count)', 'ca_mito_count', '.0f')
		self.registerLoggedQuantity('IP₃\n(nM)',         'ip3_nM',        '.1f')
		self.registerLoggedQuantity('SOCE\n(nM/s)',      'soce_flux_nMs', '.2f')
		self.registerLoggedQuantity('ATP pump\n(ions/s)', 'atp_pump_per_s', '.0f')
		self.registerLoggedQuantity('PKC*\n(count)',     'pkc_active',     '.0f')
		self.registerLoggedQuantity('P2Y1 des\n(frac)',  'p2y1_desensitised_frac', '.3f')
		self.registerLoggedQuantity('PLCβ-P\n(frac)',    'plcb_phosphorylated_frac', '.3f')
		self.registerLoggedQuantity('cAMP\n(µM)',        'camp_uM',        '.3f')
		self.registerLoggedQuantity('VASP-P\n(frac)',    'vasp_phos_frac', '.3f')

		# Live CSV (for the live-plot viewer): written into simOut when enabled.
		if self._config.live:
			self._live_path = os.path.join(sim._outputDir, 'live.csv')
			self._live_file = open(self._live_path, 'w', buffering=1)
			self._live_file.write('time,ca_cyt_nM,ca_dts_uM,ip3_nM,soce_flux_nMs\n')

	def update(self):
		counts = self._bulk_molecules.counts()

		ca_cyt_uM = counts[self._idx_ca_cyt] * _UM_PER_COUNT_CYT
		ca_dts_uM = counts[self._idx_ca_dts] * _UM_PER_COUNT_DTS
		ip3_uM    = counts[self._idx_ip3]    * _UM_PER_COUNT_CYT
		stim1_dim = counts[self._idx_stim1]
		orai      = counts[self._idx_orai]

		self.ca_cyt_nM = float(ca_cyt_uM * 1e3)
		self.ca_dts_uM = float(ca_dts_uM)
		self.ip3_nM    = float(ip3_uM * 1e3)
		self.stim1_dim = int(stim1_dim)

		# CaM and PMCA sub-states.
		self.cam_free        = int(counts[self._idx_cam_free])
		self.ca2_cam         = int(counts[self._idx_ca2_cam])
		self.ca4_cam         = int(counts[self._idx_ca4_cam])
		self.ca4_cam_pmca    = int(counts[self._idx_ca4_cam_pmca])
		self.ca4_cam_pmca_ca = int(counts[self._idx_ca4_cam_pmca_ca])
		self.pmca_cam        = int(counts[self._idx_pmca_cam])
		self.pmca_free       = int(counts[self._idx_pmca_free])
		self.pmca_ca         = int(counts[self._idx_pmca_ca])

		# Mitochondrial matrix Ca²⁺ pool (raw count — mito volume isn't tracked
		# separately) and the MCU uptake flux (ions/s), recomputed from the same
		# Hill kinetics the ODE applies, so the mito module is auditable (#76).
		self.ca_mito_count = int(counts[self._idx_ca_mito])
		ca_n = ca_cyt_uM ** K_MITO['n_MCU']
		km_n = K_MITO['K_MCU'] ** K_MITO['n_MCU']
		# Capacity back-pressure factor, matching the ODE (#76 Part 1).
		mito_fill = max(0.0, 1.0 - self.ca_mito_count / K_MITO['C_max'])
		self.mcu_uptake_per_s = float(
			K_MITO['V_max_MCU'] * self._mcu_vmax_scale
			* ca_n / (km_n + ca_n) * mito_fill)
		# #76 Part 2 — the evoked IP3R-release gate the ODE applies, via the shared
		# `ip3r_relief_factor` helper so the recorded value matches the ODE exactly
		# (one source of truth). Recorded for audit; SOCE is NOT gated by it (it
		# falls indirectly via the fuller store). End-of-step snapshot: computed
		# from the committed Ca²⁺, so it is an estimate of the within-step factor.
		self.mito_coupling_factor = float(ip3r_relief_factor(
			ca_cyt_uM, self._mcu_vmax_scale, self._mito_coupling_gain))

		# PKC feedback (v0.6): active PKC count + fraction of the P2Y1 pool
		# in the desensitised phospho-state.
		self.pkc_active = int(counts[self._idx_pkc_active])
		p2y1_total = (counts[self._idx_p2y1_inactive]
			+ counts[self._idx_p2y1_active]
			+ counts[self._idx_p2y1_desens])
		self.p2y1_desensitised_frac = float(
			counts[self._idx_p2y1_desens] / p2y1_total) if p2y1_total > 0 else 0.0
		# Fraction of the PLCβ pool phosphorylated out of the Gq-activatable
		# pool by PKC (v0.6 Slice 3, Purvis route).
		plcb_total = (counts[self._idx_plcb_inactive]
			+ counts[self._idx_plcb_active]
			+ counts[self._idx_plcb_phos])
		self.plcb_phosphorylated_frac = float(
			counts[self._idx_plcb_phos] / plcb_total) if plcb_total > 0 else 0.0

		# Inhibitory axis (v0.7 Slice 2, #10): cAMP, PKA activity, P2Y12
		# occupancy, and the clinical VASP/PRI readout (phospho-VASP fraction).
		camp_count = counts[self._idx_camp]
		self.camp_uM = float(camp_count * _UM_PER_COUNT_CYT)
		self.pka_frac = float(pka_active_frac(camp_count))
		p2y12_total = (counts[self._idx_p2y12_inactive]
			+ counts[self._idx_p2y12_active])
		self.p2y12_active_frac = float(
			counts[self._idx_p2y12_active] / p2y12_total) if p2y12_total > 0 else 0.0
		vasp_total = counts[self._idx_vasp] + counts[self._idx_vasp_phos]
		self.vasp_phos_frac = float(
			counts[self._idx_vasp_phos] / vasp_total) if vasp_total > 0 else 0.0

		# Per-step ATP consumed by SERCA + PMCA pumping. `getattr` guards the
		# initial update, which runs before the first calculateRequest sets
		# `_atp_cost`.
		self.atp_pump_per_s = float(
			getattr(self._calcium_dynamics, '_atp_cost', 0))

		# Instantaneous SOCE flux (nM/s into cytosol) via the same Dolan
		# Eq. 2/3/4 chain used inside the ODE: Hill puncta entry → MWC
		# equilibrium → Eq. 4 Nernst current.
		if ca_cyt_uM > 0.0:
			hill = (ca_cyt_uM ** PUNCTA['n']
					/ (PUNCTA['KM_uM'] ** PUNCTA['n'] + ca_cyt_uM ** PUNCTA['n']))
		else:
			hill = 0.0
		qp = PUNCTA['alpha'] * hill + PUNCTA['baseline']
		# stim1_dim is dimer particle count; MWC takes dimers directly.
		stim2_p = qp * stim1_dim
		n_orai_channels = orai / ORAI_SUBUNITS_PER_CHANNEL
		po_orai, _sf = _mwc_open_fraction(stim2_p, n_orai_channels)
		# SOCE current — physically zero when there is no extracellular Ca²⁺
		# (Dolan Fig. 4 EDTA condition); see calcium_signalling._ode_rhs.
		if self._ca_ex_uM > 0.0 and ca_cyt_uM > 0.0:
			e_ca_pm_v = RT_OVER_zF_V * math.log(self._ca_ex_uM / ca_cyt_uM)
			driving_pm_v = V_PM_V - e_ca_pm_v
			soce_ions_s = (
				-GAMMA_SOC_S * n_orai_channels * po_orai * driving_pm_v * NA_OVER_zF
			)  # #76 Part 2: SOCE not gated directly — falls via the fuller store
		else:
			soce_ions_s = 0.0
		# Convert ions/s in cytosol → nM/s (cyt volume).
		self.soce_flux_nMs = float(soce_ions_s * _UM_PER_COUNT_CYT * 1e3)

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			concentration_units='nM for ca_cyt and ip3; µM for ca_dts; nM/s for soce_flux; count for ca_mito; ions/s for mcu_uptake',
			validation_target='Dolan & Diamond 2014 Fig. 4',
		)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time=self.time(),
			simulationStep=self.simulationStep(),
			ca_cyt_nM=self.ca_cyt_nM,
			ca_dts_uM=self.ca_dts_uM,
			ip3_nM=self.ip3_nM,
			soce_flux_nMs=self.soce_flux_nMs,
			stim1_dim=self.stim1_dim,
			cam_free=self.cam_free,
			ca2_cam=self.ca2_cam,
			ca4_cam=self.ca4_cam,
			ca4_cam_pmca=self.ca4_cam_pmca,
			ca4_cam_pmca_ca=self.ca4_cam_pmca_ca,
			pmca_cam=self.pmca_cam,
			pmca_free=self.pmca_free,
			pmca_ca=self.pmca_ca,
			ca_mito_count=self.ca_mito_count,
			mcu_uptake_per_s=self.mcu_uptake_per_s,
			mito_coupling_factor=self.mito_coupling_factor,
			atp_pump_per_s=self.atp_pump_per_s,
			pkc_active=self.pkc_active,
			p2y1_desensitised_frac=self.p2y1_desensitised_frac,
			plcb_phosphorylated_frac=self.plcb_phosphorylated_frac,
			camp_uM=self.camp_uM,
			pka_active_frac=self.pka_frac,
			p2y12_active_frac=self.p2y12_active_frac,
			vasp_phos_frac=self.vasp_phos_frac,
		)
		if self._live_file is not None:
			self._live_file.write(
				f'{self.time():.0f},{self.ca_cyt_nM:.2f},'
				f'{self.ca_dts_uM:.3f},{self.ip3_nM:.2f},{self.soce_flux_nMs:.4f}\n'
			)
