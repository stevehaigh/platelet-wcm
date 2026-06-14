"""
Thromboxane A₂ synthesis dataclass (v0.61 Slice A — production only).

Active PKC + Ca²⁺ drive cytosolic phospholipase A₂ (cPLA₂), liberating
arachidonic acid; COX-1 (PTGS1) and thromboxane synthase (TBXAS1) convert it
to thromboxane A₂. The whole enzymatic chain is lumped into a single
Ca²⁺ × PKC-gated production term scaled by ``cox1_factor`` — the aspirin knob
(1 = intact COX-1, 0 = irreversibly acetylated / aspirinised). TXA₂ is
short-lived (t½ ≈ 30 s) and decays to the stable, ELISA-measurable metabolite
TXB₂.

The gate keys off PKC activation *above* a resting-tone floor (as in
GranuleSecretion), so resting TXA₂ production is exactly zero.

This slice is **production only** — TXA₂ does not yet feed back onto Gq (the
autocrine TP → Gq loop is Slice B), so the calcium ODE is untouched. Design:
``reports/design/pkc-downstream-effects-2026-06-12.qmd`` §2.
"""

class ThromboxaneSynthesis:
	"""Parameters for the ThromboxaneSynthesis process."""

	def __init__(self, sim_data):
		# Maximum TXA₂ production rate at full gate (molecules·s⁻¹). The COX-1
		# availability factor (aspirin knob) is the per-run RunConfig field
		# cox1_factor, read by the process in initialize().
		self.k_prod = 2500.0

		# TXA₂ → TXB₂ first-order decay (t½ ≈ 30 s).
		self.k_decay = 0.0231        # s⁻¹  (ln 2 / 30 s)

		# Coincidence-detector gate (shared form with GranuleSecretion):
		#   pkc_drive = max(PKC* − PKC_floor, 0); pkc_term = drive/(drive+K_pkc)
		#   ca_term   = Ca^n / (Ca^n + K_ca^n);   gate = pkc_term · ca_term
		self.PKC_floor_uM = 0.05
		self.K_pkc_uM = 0.30
		self.K_ca_uM  = 0.40
		self.n_ca     = 2

		# Species the process reads / writes.
		self.pkc_active_id = 'PKC_active[c]'
		self.ca_cyt_id     = 'CA2_CYT[c]'
		self.txa2_id       = 'TXA2[e]'
		self.txb2_id       = 'TXB2[e]'
