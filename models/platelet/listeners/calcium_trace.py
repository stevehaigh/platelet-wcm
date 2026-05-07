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
"""

import numpy as np

import math

import wholecell.listeners.listener
from reconstruction.platelet.dataclasses.process import calcium_signalling as cs_mod
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	_IDX,
	_UM_PER_COUNT_CYT,
	_UM_PER_COUNT_DTS,
	_mwc_open_fraction,
	GAMMA_SOC_S,
	NA_OVER_zF,
	ORAI_SUBUNITS_PER_CHANNEL,
	PUNCTA,
	RT_OVER_zF_V,
	V_PM_V,
)
# CA_EX_UM is read via cs_mod.CA_EX_UM rather than imported by value,
# so a runscript that overrides it (e.g. for the Phase 3 EDTA condition)
# is reflected in the listener's SOCE trace.


class CalciumTrace(wholecell.listeners.listener.Listener):
	"""Record Ca²⁺ dynamics state each timestep."""

	_name = 'CalciumTrace'

	# If set by a runscript, a lightweight CSV is written here each timestep
	# with an immediate flush so a live-plot viewer can tail it in real time.
	# Format: time,ca_cyt_nM,ca_dts_uM,ip3_nM,soce_flux_nMs
	_live_path = None

	def __init__(self, *args, **kwargs):
		self._bulk_molecules = None
		self._live_file = None
		super().__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super().initialize(sim, sim_data)

		self._bulk_molecules = (
			sim.internal_states['BulkMolecules'].container)

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

		self.registerLoggedQuantity('Ca²⁺_cyt\n(nM)',   'ca_cyt_nM',     '.1f')
		self.registerLoggedQuantity('Ca²⁺_dts\n(µM)',   'ca_dts_uM',     '.1f')
		self.registerLoggedQuantity('IP₃\n(nM)',         'ip3_nM',        '.1f')
		self.registerLoggedQuantity('SOCE\n(nM/s)',      'soce_flux_nMs', '.2f')

		if self._live_path is not None:
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
		if cs_mod.CA_EX_UM > 0.0 and ca_cyt_uM > 0.0:
			e_ca_pm_v = RT_OVER_zF_V * math.log(cs_mod.CA_EX_UM / ca_cyt_uM)
			driving_pm_v = V_PM_V - e_ca_pm_v
			soce_ions_s = (
				-GAMMA_SOC_S * n_orai_channels * po_orai * driving_pm_v * NA_OVER_zF
			)
		else:
			soce_ions_s = 0.0
		# Convert ions/s in cytosol → nM/s (cyt volume).
		self.soce_flux_nMs = float(soce_ions_s * _UM_PER_COUNT_CYT * 1e3)

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			concentration_units='nM for ca_cyt and ip3; µM for ca_dts; nM/s for soce_flux',
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
		)
		if self._live_file is not None:
			self._live_file.write(
				f'{self.time():.0f},{self.ca_cyt_nM:.2f},'
				f'{self.ca_dts_uM:.3f},{self.ip3_nM:.2f},{self.soce_flux_nMs:.4f}\n'
			)
