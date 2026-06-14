"""
GranuleSecretion process — platelet whole-cell model (v0.61 Slice 1).

Active PKC together with cytosolic Ca²⁺ drives SNARE-mediated dense- and
α-granule fusion, releasing their cargo. This process relocates the cargo that
already exists as resting scaffolding from the granule lumen to the
extracellular / open-canalicular space ``[e]`` — and, for the membrane marker
P-selectin, to a surface state ``SELP_surface[pl]``.

Per granule pool the fusion rate is first-order in the remaining cargo, scaled
by a PKC_active × Ca²⁺ coincidence gate in [0, 1):

    gate = PKC*/(PKC* + K_pkc) · Ca^n/(Ca^n + K_ca^n)

The gate is multiplicative in PKC_active, which is zero in the resting platelet,
so resting secretion is exactly zero (gate 1, the resting-quiescence invariant).

The Ca²⁺ ODE is untouched: secreted ADP is not yet fed back onto P2Y1, so the
Dolan & Diamond Ca²⁺ validation is unaffected (the autocrine loop is a later
slice). Rate constants live in
``reconstruction/platelet/dataclasses/process/granule_secretion.py``.
"""

import numpy as np

import wholecell.processes.process as process
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	pkc_ca_gate,
)


class GranuleSecretion(process.Process):
	"""PKC + Ca²⁺-gated dense- and α-granule release."""

	_name = 'GranuleSecretion'

	def __init__(self):
		super(GranuleSecretion, self).__init__()

	def initialize(self, sim, sim_data):
		super(GranuleSecretion, self).initialize(sim, sim_data)

		gs = sim_data.process.granule_secretion
		self._PKC_floor_uM = gs.PKC_floor_uM
		self._K_pkc_uM = gs.K_pkc_uM
		self._K_ca_uM = gs.K_ca_uM
		self._n_ca = gs.n_ca

		source_ids = [c[0] for c in gs.cargo]
		dest_ids = [c[1] for c in gs.cargo]
		self._k_sec = np.array(
			[gs.k_sec[c[2]] for c in gs.cargo], dtype=np.float64)

		self._sources = self.bulkMoleculesView(
			np.array(source_ids, dtype='U30'))
		self._dests = self.bulkMoleculesView(
			np.array(dest_ids, dtype='U30'))
		self._drivers = self.bulkMoleculesView(
			np.array([gs.pkc_active_id, gs.ca_cyt_id], dtype='U30'))

		# Ecto-NTPDase clearance: ADP[e] (a secretion destination) → AMP[e].
		self._k_ntpdase = gs.k_ntpdase
		self._adp_dest_idx = dest_ids.index(gs.adp_ex_id)
		self._amp = self.bulkMoleculesView(
			np.array([gs.amp_ex_id], dtype='U30'))

		# Exposed for the SecretionTrace listener (diagnostic only).
		self._gate = 0.0
		self._released = np.zeros(len(gs.cargo), dtype=np.int64)
		self._cleared = 0

	def _compute_gate(self, pkc_count, ca_count):
		"""PKC_active × Ca²⁺ coincidence gate in [0, 1); 0 at resting tone.

		Keys off PKC activation *above* the resting floor, so the small tonic
		resting PKC_active does not drive secretion (gate exactly 0 at rest).
		Shared definition: ``calcium_signalling.pkc_ca_gate`` (same gate used by
		ThromboxaneSynthesis).
		"""
		return pkc_ca_gate(pkc_count, ca_count, self._PKC_floor_uM,
			self._K_pkc_uM, self._K_ca_uM, self._n_ca)

	def calculateRequest(self):
		dt = self.timeStepSec()
		pkc_count, ca_count = self._drivers.total_counts()
		self._gate = self._compute_gate(pkc_count, ca_count)

		src_counts = self._sources.total_counts()
		released_frac = 1.0 - np.exp(-self._k_sec * self._gate * dt)
		self._released = np.floor(src_counts * released_frac).astype(np.int64)
		self._sources.requestIs(self._released)

		# Ecto-NTPDase clearance of the existing ADP[e] pool (first-order).
		adp_e_count = self._dests.total_counts()[self._adp_dest_idx]
		self._cleared = int(np.floor(
			adp_e_count * (1.0 - np.exp(-self._k_ntpdase * dt))))
		# Request the cleared ADP[e] (the only destination species we consume).
		dest_request = np.zeros(len(self._released), dtype=np.int64)
		dest_request[self._adp_dest_idx] = self._cleared
		self._dests.requestIs(dest_request)

	def evolveState(self):
		self._sources.countsDec(self._released)
		# Net destination change: secretion adds cargo; clearance removes ADP[e].
		dest_delta = self._released.copy()
		dest_delta[self._adp_dest_idx] -= self._cleared
		self._dests.countsInc(dest_delta)
		self._amp.countsInc(np.array([self._cleared], dtype=np.int64))
