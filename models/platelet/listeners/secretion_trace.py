"""
SecretionTrace listener for the platelet whole-cell model (v0.61 Slice 1).

Records granule-release state each timestep: secreted dense-granule cargo
(ADP, serotonin) and α-granule cargo (fibrinogen) in the extracellular /
open-canalicular space ``[e]``, surface P-selectin, and the released /
surface-exposed fractions — the standard single-platelet activation readouts
(P-selectin exposure by flow cytometry; serotonin / ATP release by
lumi-aggregometry).

Columns written:
  time                       — simulation time (s)
  adp_e                      — secreted ADP count ([e])
  adp_e_uM                   — pericellular secreted-ADP concentration (µM);
                               the autocrine P2Y1 drive (v0.61 Slice 2)
  amp_e                      — ecto-NTPDase product AMP count ([e])
  serotonin_e                — secreted serotonin (5-HT) count ([e])
  fibrinogen_e               — secreted fibrinogen (FGA) count ([e])
  pselectin_surface          — surface P-selectin count ([pl])
  adp_released_frac          — fraction of the dense-granule ADP pool released
  serotonin_released_frac    — fraction of the dense-granule 5-HT pool released
  fibrinogen_released_frac   — fraction of the α-granule FGA pool released
  pselectin_surface_frac     — fraction of the P-selectin pool on the surface
  secretion_gate             — PKC* × Ca²⁺ secretion gate value (0–1)
"""

import wholecell.listeners.listener
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	_UM_PER_COUNT_EX,
)


class SecretionTrace(wholecell.listeners.listener.Listener):
	"""Record granule-secretion state each timestep."""

	_name = 'SecretionTrace'

	def __init__(self, *args, **kwargs):
		self._bulk_molecules = None
		super().__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super().initialize(sim, sim_data)

		self._bulk_molecules = sim.internal_states['BulkMolecules'].container
		self._secretion = sim.processes['GranuleSecretion']

		all_ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
		self._idx_adp_dg   = all_ids.index('ADP[dg]')
		self._idx_adp_e    = all_ids.index('ADP[e]')
		self._idx_amp_e    = all_ids.index('AMP[e]')
		self._idx_5ht_dg   = all_ids.index('5HT[dg]')
		self._idx_5ht_e    = all_ids.index('5HT[e]')
		self._idx_fga_ag   = all_ids.index('FGA[ag]')
		self._idx_fga_e    = all_ids.index('FGA[e]')
		self._idx_selp_ag  = all_ids.index('SELP[ag]')
		self._idx_selp_sfc = all_ids.index('SELP_surface[pl]')

		self.adp_e = 0
		self.adp_e_uM = 0.0
		self.amp_e = 0
		self.serotonin_e = 0
		self.fibrinogen_e = 0
		self.pselectin_surface = 0
		self.adp_released_frac = 0.0
		self.serotonin_released_frac = 0.0
		self.fibrinogen_released_frac = 0.0
		self.pselectin_surface_frac = 0.0
		self.secretion_gate = 0.0

		self.registerLoggedQuantity(
			'ADP rel\n(frac)', 'adp_released_frac', '.3f')
		self.registerLoggedQuantity(
			'P-sel\n(frac)', 'pselectin_surface_frac', '.3f')

	@staticmethod
	def _frac(released, source):
		total = released + source
		return float(released / total) if total > 0 else 0.0

	def update(self):
		counts = self._bulk_molecules.counts()

		self.adp_e = int(counts[self._idx_adp_e])
		self.adp_e_uM = float(counts[self._idx_adp_e] * _UM_PER_COUNT_EX)
		self.amp_e = int(counts[self._idx_amp_e])
		self.serotonin_e = int(counts[self._idx_5ht_e])
		self.fibrinogen_e = int(counts[self._idx_fga_e])
		self.pselectin_surface = int(counts[self._idx_selp_sfc])

		self.adp_released_frac = self._frac(
			counts[self._idx_adp_e], counts[self._idx_adp_dg])
		self.serotonin_released_frac = self._frac(
			counts[self._idx_5ht_e], counts[self._idx_5ht_dg])
		self.fibrinogen_released_frac = self._frac(
			counts[self._idx_fga_e], counts[self._idx_fga_ag])
		self.pselectin_surface_frac = self._frac(
			counts[self._idx_selp_sfc], counts[self._idx_selp_ag])

		self.secretion_gate = float(getattr(self._secretion, '_gate', 0.0))

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			units='counts for [e]/surface species; fractions are dimensionless',
		)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time=self.time(),
			simulationStep=self.simulationStep(),
			adp_e=self.adp_e,
			adp_e_uM=self.adp_e_uM,
			amp_e=self.amp_e,
			serotonin_e=self.serotonin_e,
			fibrinogen_e=self.fibrinogen_e,
			pselectin_surface=self.pselectin_surface,
			adp_released_frac=self.adp_released_frac,
			serotonin_released_frac=self.serotonin_released_frac,
			fibrinogen_released_frac=self.fibrinogen_released_frac,
			pselectin_surface_frac=self.pselectin_surface_frac,
			secretion_gate=self.secretion_gate,
		)
