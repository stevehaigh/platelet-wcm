"""
Integrin αIIbβ3 inside-out activation dataclass (v0.61 §3 — minimal model).

Active PKC together with the Ca²⁺-sensing GEF CalDAG-GEFI activates Rap1b,
which recruits talin/kindlin to the αIIbβ3 (GPIIb-IIIa) cytoplasmic tail and
switches the integrin from a low-affinity (resting) to a high-affinity (active)
conformation — *inside-out* activation. The high-affinity state binds
fibrinogen/VWF and is the conformation reported by the activation-specific
antibody **PAC-1** (flow cytometry).

This is the **2-state minimal model** from the design doc: a single
resting ⇌ active conformational switch on the assembled surface heterodimer.
The forward (activating) rate is gated by the same PKC_active × Ca²⁺
coincidence term used by GranuleSecretion / ThromboxaneSynthesis — which is
exactly the inside-out trigger (PKC + CalDAG-GEFI/Ca²⁺) lumped into one term —
so resting activation is exactly zero (the resting-quiescence invariant). A
slow first-order reversion returns unstimulated integrin to the resting state.

The **active fraction is the per-cell PAC-1 readout**. Aggregation itself is
inter-cellular and out of scope for a single-platelet model (see the design
doc); only the affinity *state* is represented here. The activation rate is
scaled by ``RunConfig.integrin_act_scale`` (1 = intact; 0 = αIIbβ3 antagonist /
Glanzmann thrombasthenia), read by the process in initialize().

This slice adds no calcium-ODE state — αIIbβ3 is a terminal output, not a Gq
feedback loop — so the Dolan goldens are preserved (byte-identical). Design:
``reports/design/pkc-downstream-effects-2026-06-12.qmd`` §3.
"""

# αIIbβ3 conformational states on the plasmalemma (the assembled ITGA2B/ITGB3
# surface heterodimer). Both states carry the same mass — activation is a
# conformational switch, not synthesis — so the resting ⇌ active transfer is
# mass-neutral. RestingDecay excludes these (their counts are managed here).
#
# Mass accounting (species TSV): these states carry the β3 (ITGB3, ~90 kDa)
# partner's mass only. The αIIb subunit mass is already in the inventory as the
# pre-existing ITGA2B[c] (~129 kDa, 80 000 copies), so the assembled αIIbβ3
# receptor is mass-counted exactly once across the two entries (no double-count).
INTEGRIN_MANAGED_SPECIES = (
	'aIIbb3_resting[pl]',
	'aIIbb3_active[pl]',
)


class IntegrinActivation:
	"""Parameters for the IntegrinActivation process."""

	def __init__(self, sim_data):
		# First-order conformational switch rates (s⁻¹). The forward rate is the
		# value at full gate; it is scaled each step by the PKC×Ca gate (and by
		# the per-run RunConfig.integrin_act_scale knockout knob). The reverse
		# (reversion) rate is gate-independent — unstimulated integrin relaxes
		# back to low affinity. k_act ≫ k_inact so a sustained agonist drives a
		# high active (PAC-1⁺) fraction; the ratio sets the steady-state level.
		self.k_act = 0.05      # resting → active, at full gate
		self.k_inact = 0.005   # active → resting reversion (slow)

		# Coincidence-detector gate (shared form with GranuleSecretion /
		# ThromboxaneSynthesis): gate = pkc_term · ca_term, each in [0, 1).
		# The floor makes resting activation exactly zero.
		self.PKC_floor_uM = 0.05
		self.K_pkc_uM = 0.30
		self.K_ca_uM  = 0.40
		self.n_ca     = 2

		# ── #73 — PI3K/Akt → Rap1b arm (P2Y12-sustained activation) ──────────
		# Zou 2022 (IJMS 23:12512) / Stolla 2011 (Blood): Rap1b-GTP is the
		# integrin's *proximal* driver. It FORMS via the fast inside-out trigger
		# (the existing PKC×Ca gate — already 0 at rest, so the resting-quiescence
		# invariant holds) and is DEACTIVATED by the Rap1b-GAP Rasa3. Akt (lumped
		# PI3K/Akt, driven by P2Y12 occupancy = the autocrine-ADP / Gi arm)
		# dis-inhibits Rasa3, so Rap stays high *while P2Y12 is driven* and decays
		# when ADP clears → P2Y12 off → Akt off → GAP re-engages → integrin
		# REVERSES (the persistence/reversibility axis the lumped gate could not
		# capture). Rap1b-GTP replaces the gate as the integrin's forward driver:
		#   k_act_eff = k_act · Rap1b_GTP · pka_brake · act_scale · rap1b_scale
		# Two-timescale (Stolla): fast formation (gate), slow Akt-modulated GAP.
		# All model-choice (no measured constants); calcium ODE untouched →
		# Dolan goldens preserved. Design: pi3k-akt-rap1b-arm-2026-06-22.qmd.
		self.k_rap_form = 0.30   # Rap1b-GTP formation at full gate (s⁻¹, fast rise)
		self.k_rap_gap  = 0.08   # Rasa3 GAP deactivation of Rap1b-GTP (s⁻¹, slow)
		self.f_akt_gap  = 0.90   # max fractional Rasa3-GAP inhibition by Akt (0–1)
		self.k_akt_on   = 0.04   # Akt activation by P2Y12 occupancy (s⁻¹, slow)
		self.k_akt_off  = 0.02   # Akt deactivation when P2Y12 clears (s⁻¹)

		# Species the process reads (gate drivers / P2Y12) / writes (the states).
		self.pkc_active_id   = 'PKC_active[c]'
		self.ca_cyt_id       = 'CA2_CYT[c]'
		self.p2y12_active_id = 'P2Y12_active[pl]'
		self.resting_id      = 'aIIbb3_resting[pl]'
		self.active_id       = 'aIIbb3_active[pl]'
