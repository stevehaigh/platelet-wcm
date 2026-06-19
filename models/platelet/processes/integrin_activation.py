"""
IntegrinActivation process — platelet whole-cell model (v0.61 §3).

PKC together with CalDAG-GEFI/Ca²⁺ → Rap1b → talin/kindlin switches αIIbβ3
(GPIIb-IIIa) from its resting (low-affinity) to its active (high-affinity)
conformation — *inside-out* activation. The high-affinity state is the one the
activation-specific antibody PAC-1 reports, so the **active fraction is the
per-cell PAC-1 readout**.

This is the 2-state minimal model: resting ⇌ active on the assembled surface
heterodimer. The forward (activating) rate is first-order in the resting pool,
scaled by the shared PKC_active × Ca²⁺ coincidence gate (the lumped inside-out
trigger) and by the per-run ``integrin_act_scale`` knob (1 = intact;
0 = αIIbβ3 antagonist / Glanzmann thrombasthenia). A slow gate-independent
reverse rate relaxes unstimulated integrin back to resting.

The gate keys off PKC activation *above* a resting-tone floor, so resting
activation is exactly zero (resting-quiescence invariant). αIIbβ3 is a terminal
output (no Gq feedback), so the calcium ODE is untouched and the Dolan goldens
are preserved. Rate constants live in
``reconstruction/platelet/dataclasses/process/integrin_activation.py``.
"""

import numpy as np

import wholecell.processes.process as process
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	pkc_ca_gate,
	pka_brake_factor,
	K_PKA,
)


class IntegrinActivation(process.Process):
	"""PKC + Ca²⁺-gated αIIbβ3 inside-out (resting ⇌ active) switch."""

	_name = 'IntegrinActivation'

	def __init__(self):
		super(IntegrinActivation, self).__init__()

	def initialize(self, sim, sim_data):
		super(IntegrinActivation, self).initialize(sim, sim_data)

		# αIIbβ3 antagonist / Glanzmann knob — per-run, read once from RunConfig.
		self._act_scale = sim.run_config.integrin_act_scale

		ig = sim_data.process.integrin_activation
		self._k_act = ig.k_act
		self._k_inact = ig.k_inact
		self._PKC_floor_uM = ig.PKC_floor_uM
		self._K_pkc_uM = ig.K_pkc_uM
		self._K_ca_uM = ig.K_ca_uM
		self._n_ca = ig.n_ca

		self._drivers = self.bulkMoleculesView(
			np.array([ig.pkc_active_id, ig.ca_cyt_id], dtype='U30'))
		self._resting = self.bulkMoleculesView(
			np.array([ig.resting_id], dtype='U30'))
		self._active = self.bulkMoleculesView(
			np.array([ig.active_id], dtype='U30'))

		# PKA inhibitory brake on inside-out activation (v0.7 Slice 2, #10).
		# cAMP/PKA brakes the Rap1→integrin step; ADP→P2Y12→Gi lowers cAMP/PKA
		# → dis-inhibits αIIbβ3 (the dominant, clinically-relevant P2Y12
		# effect). Normalised to 1.0 at resting cAMP → resting/Dolan unchanged.
		self._camp = self.bulkMoleculesView(np.array(['cAMP[c]'], dtype='U30'))
		self._pka_brake_gain = K_PKA['integrin_brake_gain']
		self._pka_b_max = K_PKA['integrin_b_max']

		# Exposed for the IntegrinTrace listener (diagnostic only).
		self._gate = 0.0
		self._pka_brake = 1.0
		self._activated = 0
		self._reverted = 0

	def _compute_gate(self, pkc_count, ca_count):
		"""PKC_active × Ca²⁺ coincidence gate in [0, 1); 0 at resting tone.

		Shared definition: ``calcium_signalling.pkc_ca_gate`` (same gate used by
		GranuleSecretion and ThromboxaneSynthesis).
		"""
		return pkc_ca_gate(pkc_count, ca_count, self._PKC_floor_uM,
			self._K_pkc_uM, self._K_ca_uM, self._n_ca)

	def calculateRequest(self):
		dt = self.timeStepSec()
		pkc_count, ca_count = self._drivers.total_counts()
		self._gate = self._compute_gate(pkc_count, ca_count)

		# PKA dis-inhibition brake (≥1): 1.0 at resting cAMP, >1 once ADP →
		# P2Y12 → Gi lowers cAMP/PKA, releasing the brake on integrin.
		camp_count = self._camp.total_counts()[0]
		self._pka_brake = pka_brake_factor(
			camp_count, self._pka_brake_gain, self._pka_b_max)

		# Forward: resting → active (gated by PKC×Ca, PKA brake, knockout knob).
		resting_count = self._resting.total_counts()[0]
		act_frac = 1.0 - np.exp(
			-self._k_act * self._act_scale * self._gate * self._pka_brake * dt)
		self._activated = int(np.floor(resting_count * act_frac))

		# Reverse: active → resting (slow, gate-independent relaxation).
		active_count = self._active.total_counts()[0]
		self._reverted = int(np.floor(
			active_count * (1.0 - np.exp(-self._k_inact * dt))))

		# Reserve the counts each direction consumes.
		self._resting.requestIs(np.array([self._activated], dtype=np.int64))
		self._active.requestIs(np.array([self._reverted], dtype=np.int64))

	def evolveState(self):
		# Net conformational transfer (mass-neutral: same molecule, two states).
		self._resting.countsInc(
			np.array([self._reverted - self._activated], dtype=np.int64))
		self._active.countsInc(
			np.array([self._activated - self._reverted], dtype=np.int64))
