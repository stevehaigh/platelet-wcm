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
	N_P2Y12_TOTAL,
)


def akt_rap_step(akt, rap, gate, p2y12_frac, dt,
		k_akt_on, k_akt_off, k_rap_form, k_rap_gap, f_akt_gap):
	"""One Euler-exact 1 s step of the PI3K/Akt → Rap1b-GTP arm (#73).

	Pure function (no state) so it can be unit-tested in isolation and stays the
	single source of truth for the mechanism (cf. ``ip3r_relief_factor`` for the
	MCU coupling). Akt (lumped PI3K/Akt) tracks P2Y12 occupancy — the
	autocrine-ADP / Gi sustain arm — and dis-inhibits the Rap1b-GAP Rasa3.
	Rap1b-GTP FORMS via the inside-out gate (already 0 at rest, so the
	resting-quiescence invariant holds) and is removed by the GAP; raised Akt
	slows the GAP, so Rap stays high *while P2Y12 is driven* and decays once ADP
	clears. Each species relaxes toward its instantaneous steady state (exact
	over a constant-rate step), keeping both in [0, 1). The GAP term is clamped
	≥ 0 so full Akt inhibition cannot drive a negative removal rate. Akt is
	stepped first; Rap's GAP uses the updated Akt. Returns ``(akt_next,
	rap_next)``.
	"""
	akt_in = k_akt_on * p2y12_frac
	akt_tot = akt_in + k_akt_off
	akt_ss = akt_in / akt_tot if akt_tot > 0.0 else 0.0
	akt_next = akt_ss + (akt - akt_ss) * np.exp(-akt_tot * dt)

	a_form = k_rap_form * gate
	b_gap = k_rap_gap * max(0.0, 1.0 - f_akt_gap * akt_next)
	rap_tot = a_form + b_gap
	rap_ss = a_form / rap_tot if rap_tot > 0.0 else 0.0
	rap_next = rap_ss + (rap - rap_ss) * np.exp(-rap_tot * dt)
	return akt_next, rap_next


class IntegrinActivation(process.Process):
	"""PKC + Ca²⁺-gated αIIbβ3 inside-out (resting ⇌ active) switch."""

	_name = 'IntegrinActivation'

	def __init__(self):
		super(IntegrinActivation, self).__init__()

	def initialize(self, sim, sim_data):
		super(IntegrinActivation, self).initialize(sim, sim_data)

		# αIIbβ3 antagonist / Glanzmann knob — per-run, read once from RunConfig.
		self._act_scale = sim.run_config.integrin_act_scale
		# #73 — Rap1b knockout knob (0 = no inside-out activation).
		self._rap1b_scale = sim.run_config.rap1b_scale

		ig = sim_data.process.integrin_activation
		self._k_act = ig.k_act
		self._k_inact = ig.k_inact
		self._PKC_floor_uM = ig.PKC_floor_uM
		self._K_pkc_uM = ig.K_pkc_uM
		self._K_ca_uM = ig.K_ca_uM
		self._n_ca = ig.n_ca

		# #73 — PI3K/Akt → Rap1b arm params (the P2Y12-sustained route).
		self._k_rap_form = ig.k_rap_form
		self._k_rap_gap = ig.k_rap_gap
		self._f_akt_gap = ig.f_akt_gap
		self._k_akt_on = ig.k_akt_on
		self._k_akt_off = ig.k_akt_off

		self._drivers = self.bulkMoleculesView(
			np.array([ig.pkc_active_id, ig.ca_cyt_id], dtype='U30'))
		# #73 — P2Y12-active drives Akt (the autocrine-ADP / Gi sustain arm).
		self._p2y12 = self.bulkMoleculesView(
			np.array([ig.p2y12_active_id], dtype='U30'))
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
		# #73 — Akt (lumped PI3K/Akt) and Rap1b-GTP activity (0–1), Euler-exact
		# stepped each 1 s. Both 0 at rest → resting integrin activation exactly 0.
		self._akt = 0.0
		self._rap = 0.0

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

		# #73 — PI3K/Akt → Rap1b arm. Akt tracks P2Y12 occupancy (the autocrine-
		# ADP / Gi sustain arm) and dis-inhibits the Rap1b-GAP Rasa3. Rap1b-GTP
		# FORMS via the gate (fast, 0 at rest) and is removed by the GAP (slow).
		# Pure Euler-exact step (single source of truth, also unit-tested).
		p12 = (self._p2y12.total_counts()[0] / N_P2Y12_TOTAL
			if N_P2Y12_TOTAL else 0.0)
		self._akt, self._rap = akt_rap_step(
			self._akt, self._rap, self._gate, p12, dt,
			self._k_akt_on, self._k_akt_off,
			self._k_rap_form, self._k_rap_gap, self._f_akt_gap)

		# Forward: resting → active, driven by Rap1b-GTP (× PKA brake × knockout
		# knobs). Rap1b is now the proximal driver (was the PKC×Ca gate directly).
		resting_count = self._resting.total_counts()[0]
		act_frac = 1.0 - np.exp(
			-self._k_act * self._act_scale * self._rap1b_scale
			* self._rap * self._pka_brake * dt)
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
