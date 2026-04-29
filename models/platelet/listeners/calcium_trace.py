"""
CalciumTrace listener for the platelet whole-cell model.

Records the Ca²⁺ signalling state each timestep for post-simulation analysis
and validation against the Dolan & Diamond (2014) Fig. 4 Ca²⁺ transient.

Columns written:
  time            — simulation time (s)
  ca_cyt_nM       — free cytosolic Ca²⁺ (nM)
  ca_dts_uM       — DTS stored Ca²⁺ (µM)
  ip3_nM          — IP₃ concentration (nM)
  soce_flux_nMs   — instantaneous SOCE influx rate into cytosol (nM/s)
  stim1_dim       — STIM1 dimer count (monomer-equivalent units)
"""

import numpy as np

import wholecell.listeners.listener
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	_IDX,
	_UM_PER_COUNT_CYT,
	_UM_PER_COUNT_DTS,
	CA_EX_UM,
	K_SOCE,
)


class CalciumTrace(wholecell.listeners.listener.Listener):
	"""Record Ca²⁺ dynamics state each timestep."""

	_name = 'CalciumTrace'

	def __init__(self, *args, **kwargs):
		self._bulk_molecules = None
		super().__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super().initialize(sim, sim_data)

		self._bulk_molecules = (
			sim.internal_states['BulkMolecules'].container)

		# Pre-compute the global indices of the species we track.
		all_ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
		self._idx_ca_cyt  = all_ids.index('CA2_CYT[c]')
		self._idx_ca_dts  = all_ids.index('CA2_DTS[dts]')
		self._idx_ip3     = all_ids.index('IP3[c]')
		self._idx_stim1   = all_ids.index('STIM1_dim[dts]')

		# Initialise logged quantities for shell display.
		self.ca_cyt_nM     = 0.0
		self.ca_dts_uM     = 0.0
		self.ip3_nM        = 0.0
		self.soce_flux_nMs = 0.0
		self.stim1_dim     = 0

		self.registerLoggedQuantity('Ca²⁺_cyt\n(nM)',   'ca_cyt_nM',     '.1f')
		self.registerLoggedQuantity('Ca²⁺_dts\n(µM)',   'ca_dts_uM',     '.1f')
		self.registerLoggedQuantity('IP₃\n(nM)',         'ip3_nM',        '.1f')
		self.registerLoggedQuantity('SOCE\n(nM/s)',      'soce_flux_nMs', '.2f')

	def update(self):
		counts = self._bulk_molecules.counts()

		ca_cyt_uM = counts[self._idx_ca_cyt] * _UM_PER_COUNT_CYT
		ca_dts_uM = counts[self._idx_ca_dts] * _UM_PER_COUNT_DTS
		ip3_uM    = counts[self._idx_ip3]    * _UM_PER_COUNT_CYT
		stim1_dim = counts[self._idx_stim1]

		self.ca_cyt_nM     = float(ca_cyt_uM * 1e3)
		self.ca_dts_uM     = float(ca_dts_uM)
		self.ip3_nM        = float(ip3_uM * 1e3)
		self.stim1_dim     = int(stim1_dim)

		# Instantaneous SOCE flux (nM/s into cytosol).
		soce_uMs           = K_SOCE['k_orai'] * stim1_dim * (CA_EX_UM - ca_cyt_uM)
		self.soce_flux_nMs = float(soce_uMs * 1e3)

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
		)
