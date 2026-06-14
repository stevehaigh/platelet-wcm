"""
ThromboxaneSynthesis process — platelet whole-cell model (v0.61 Slice A).

Active PKC + Ca²⁺ activate cPLA₂ → arachidonic acid → COX-1 → thromboxane
synthase → thromboxane A₂. The chain is lumped into a single Ca²⁺ × PKC-gated
production term, scaled by ``cox1_factor`` (1 = intact COX-1; 0 = aspirin
knockout). TXA₂ is short-lived and decays first-order to the stable metabolite
TXB₂ (the ELISA readout).

The gate keys off PKC activation *above* a resting-tone floor, so resting TXA₂
production is exactly zero (resting-quiescence invariant).

This slice is **production only**: TXA₂ does not yet drive Gq (the autocrine
TP → Gq loop is Slice B), so the calcium ODE is untouched and the Dolan
goldens are preserved. Rate constants live in
``reconstruction/platelet/dataclasses/process/thromboxane_synthesis.py``.
"""

import numpy as np

import wholecell.processes.process as process
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	pkc_ca_gate,
)


class ThromboxaneSynthesis(process.Process):
	"""PKC + Ca²⁺-gated thromboxane A₂ synthesis and decay."""

	_name = 'ThromboxaneSynthesis'

	def __init__(self):
		super(ThromboxaneSynthesis, self).__init__()

	def initialize(self, sim, sim_data):
		super(ThromboxaneSynthesis, self).initialize(sim, sim_data)

		# COX-1 availability (aspirin knob) — per-run, read once from RunConfig.
		self._cox1_factor = sim.run_config.cox1_factor

		tx = sim_data.process.thromboxane_synthesis
		self._k_prod = tx.k_prod
		self._k_decay = tx.k_decay
		self._PKC_floor_uM = tx.PKC_floor_uM
		self._K_pkc_uM = tx.K_pkc_uM
		self._K_ca_uM = tx.K_ca_uM
		self._n_ca = tx.n_ca

		self._drivers = self.bulkMoleculesView(
			np.array([tx.pkc_active_id, tx.ca_cyt_id], dtype='U30'))
		self._txa2 = self.bulkMoleculesView(
			np.array([tx.txa2_id], dtype='U30'))
		self._txb2 = self.bulkMoleculesView(
			np.array([tx.txb2_id], dtype='U30'))

		# Exposed for the listener (diagnostic only).
		self._gate = 0.0
		self._produced = 0
		self._decayed = 0

	def _compute_gate(self, pkc_count, ca_count):
		"""PKC_active × Ca²⁺ coincidence gate in [0, 1); 0 at resting tone.

		Shared definition: ``calcium_signalling.pkc_ca_gate`` (same gate used by
		GranuleSecretion).
		"""
		return pkc_ca_gate(pkc_count, ca_count, self._PKC_floor_uM,
			self._K_pkc_uM, self._K_ca_uM, self._n_ca)

	def calculateRequest(self):
		dt = self.timeStepSec()
		pkc_count, ca_count = self._drivers.total_counts()
		self._gate = self._compute_gate(pkc_count, ca_count)

		# De-novo TXA₂ production (COX-1-gated by the per-run aspirin knob);
		# first-order decay → TXB₂.
		self._produced = int(round(
			self._k_prod * self._cox1_factor * self._gate * dt))
		txa2_count = self._txa2.total_counts()[0]
		self._decayed = int(np.floor(
			txa2_count * (1.0 - np.exp(-self._k_decay * dt))))
		# We consume the decayed TXA₂; production is an increment (no request).
		self._txa2.requestIs(np.array([self._decayed], dtype=np.int64))

	def evolveState(self):
		# Net TXA₂ change: produced (new) − decayed (→ TXB₂).
		self._txa2.countsInc(
			np.array([self._produced - self._decayed], dtype=np.int64))
		self._txb2.countsInc(np.array([self._decayed], dtype=np.int64))
