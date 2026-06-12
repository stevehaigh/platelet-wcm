"""
ThromboxaneTrace listener for the platelet whole-cell model (v0.61 Slice A).

Records thromboxane synthesis each timestep — the TXA₂/TXB₂ generation that
flow-cytometry / ELISA assays measure, and the aspirin (COX-1 knockout)
perturbation target.

Columns written:
  time            — simulation time (s)
  txa2            — thromboxane A₂ count ([e]; short-lived)
  txa2_uM         — pericellular TXA₂ concentration (µM) — the Slice-B TP drive
  txb2            — stable TXB₂ metabolite count ([e]; the cumulative ELISA readout)
  txa2_synth_gate — PKC* × Ca²⁺ synthesis gate value (0–1)
  tp_active_frac  — fraction of the TP receptor pool TXA₂-activated (Slice B loop)
"""

import wholecell.listeners.listener
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	_UM_PER_COUNT_EX,
)


class ThromboxaneTrace(wholecell.listeners.listener.Listener):
	"""Record thromboxane A₂ synthesis state each timestep."""

	_name = 'ThromboxaneTrace'

	def __init__(self, *args, **kwargs):
		self._bulk_molecules = None
		super().__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super().initialize(sim, sim_data)

		self._bulk_molecules = sim.internal_states['BulkMolecules'].container
		self._synthesis = sim.processes['ThromboxaneSynthesis']

		all_ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
		self._idx_txa2 = all_ids.index('TXA2[e]')
		self._idx_txb2 = all_ids.index('TXB2[e]')
		self._idx_tp_i = all_ids.index('TP_inactive[pl]')
		self._idx_tp_a = all_ids.index('TP_active[pl]')

		self.txa2 = 0
		self.txa2_uM = 0.0
		self.txb2 = 0
		self.txa2_synth_gate = 0.0
		self.tp_active_frac = 0.0

		self.registerLoggedQuantity('TXA₂\n(µM)', 'txa2_uM', '.3f')
		self.registerLoggedQuantity('TXB₂\n(count)', 'txb2', '.0f')
		self.registerLoggedQuantity('TP act\n(frac)', 'tp_active_frac', '.3f')

	def update(self):
		counts = self._bulk_molecules.counts()
		self.txa2 = int(counts[self._idx_txa2])
		self.txa2_uM = float(counts[self._idx_txa2] * _UM_PER_COUNT_EX)
		self.txb2 = int(counts[self._idx_txb2])
		self.txa2_synth_gate = float(getattr(self._synthesis, '_gate', 0.0))
		tp_total = counts[self._idx_tp_i] + counts[self._idx_tp_a]
		self.tp_active_frac = float(
			counts[self._idx_tp_a] / tp_total) if tp_total > 0 else 0.0

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			units='count for txa2/txb2; µM for txa2_uM; gate dimensionless',
		)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time=self.time(),
			simulationStep=self.simulationStep(),
			txa2=self.txa2,
			txa2_uM=self.txa2_uM,
			txb2=self.txb2,
			txa2_synth_gate=self.txa2_synth_gate,
			tp_active_frac=self.tp_active_frac,
		)
