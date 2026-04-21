"""
Resting-state protein decay dataclass for the platelet whole-cell model.

Anucleate platelets cannot replenish proteins via transcription, so all
proteins decay with a characteristic half-life (~7 days, Burkhart 2012).
"""

# Protein half-life in seconds (7 days).
# Source: Burkhart 2012 (platelet proteome turnover study).
_PROTEIN_HALF_LIFE_SEC = 7 * 24 * 3600  # 604800 s


class RestingDecay:
	"""Parameters for the RestingDecay process."""

	def __init__(self, sim_data):
		self.molecule_names = sim_data.internal_state.bulk_molecules.protein_molecule_names
		self.protein_half_life = _PROTEIN_HALF_LIFE_SEC
