"""
Granule-secretion dataclass for the platelet whole-cell model (v0.61 Slice 1).

PKC + Ca²⁺ drive SNARE-mediated dense- and α-granule release. Cargo that
already exists as resting scaffolding — ADP / serotonin in dense granules,
fibrinogen and P-selectin in α-granules — is relocated from the granule
lumen to the extracellular / open-canalicular space ``[e]``; the membrane
marker P-selectin instead moves to a surface state ``SELP_surface[pl]``.

Release is gated on a PKC_active × Ca²⁺ coincidence term with a hard resting
floor of exactly zero (the gate is multiplicative in PKC_active, which is 0 at
rest), so the un-stimulated platelet does not secrete. The rate is first-order
in the remaining granule pool, scaled by that gate.

Design: ``reports/design/pkc-downstream-effects-2026-06-12.qmd`` §1.
"""

# Proteins relocated by GranuleSecretion. Excluded from RestingDecay so the
# two processes don't both lay claim to the same counts in one timestep.
SECRETION_MANAGED_PROTEINS = (
	'FGA[ag]', 'SELP[ag]', 'FGA[e]', 'SELP_surface[pl]',
)


class GranuleSecretion:
	"""Parameters for the GranuleSecretion process."""

	def __init__(self, sim_data):
		# Cargo relocation map: (source_id, destination_id, granule_type).
		# Dense granules → [e]; α-granule fibrinogen → [e]; P-selectin → the
		# plasmalemma surface state (the canonical flow-cytometry marker).
		self.cargo = (
			('ADP[dg]',  'ADP[e]',            'dense'),
			('5HT[dg]',  '5HT[e]',            'dense'),
			('FGA[ag]',  'FGA[e]',            'alpha'),
			('SELP[ag]', 'SELP_surface[pl]',  'alpha'),
		)

		# First-order fusion rate per granule pool (s⁻¹), scaled by the gate.
		# α-granule release is slower than dense-granule release.
		self.k_sec = {'dense': 0.10, 'alpha': 0.05}

		# Coincidence-detector gate: gate = pkc_term · ca_term, each in [0, 1).
		#   pkc_drive = max(PKC* − PKC_floor, 0)   (activation above resting tone)
		#   pkc_term  = pkc_drive / (pkc_drive + K_pkc)
		#   ca_term   = Ca^n / (Ca^n + K_ca^n)
		# The floor makes resting secretion exactly zero: the un-stimulated
		# platelet carries a small tonic PKC_active (~0.6 % of the pool), which
		# sits below PKC_floor, so the gate is 0 at rest (resting-quiescence
		# invariant) yet ≈unchanged during activation (PKC* ≫ floor).
		self.PKC_floor_uM = 0.05   # resting-tone floor on active PKC (µM)
		self.K_pkc_uM = 0.30   # active-PKC half-saturation above floor (µM)
		self.K_ca_uM  = 0.40   # cytosolic-Ca²⁺ half-saturation (µM)
		self.n_ca     = 2      # Ca²⁺ Hill coefficient

		# Driver species the gate reads (cytosolic counts → µM in the process).
		self.pkc_active_id = 'PKC_active[c]'
		self.ca_cyt_id     = 'CA2_CYT[c]'

		# Ecto-NTPDase (CD39) clearance of secreted ADP: ADP[e] → AMP[e],
		# first-order (v0.61 Slice 2). Makes the autocrine ADP → P2Y1 loop
		# self-limiting — the cell recovers as secreted ADP is hydrolysed
		# (the loop is also limited by P2Y1 desensitisation and finite
		# granule cargo). CD39 clears ADP over tens of seconds.
		self.k_ntpdase = 0.05        # s⁻¹ (τ ≈ 20 s)
		self.adp_ex_id = 'ADP[e]'    # cleared species (a secretion destination)
		self.amp_ex_id = 'AMP[e]'    # hydrolysis product (sink)
