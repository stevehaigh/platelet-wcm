"""
Molecule groups for the platelet whole-cell model.

Platelets do not divide, so all division-related lists are empty.
The unique_molecules division groups mirror the E. coli names exactly because
the engine references them by those names.
"""


class MoleculeGroups:
	"""Categorised molecule lists consumed by the simulation engine."""

	def __init__(self):
		# BulkMolecules division — empty for a non-dividing cell
		self.bulk_molecules_binomial_division = []
		self.bulk_molecules_equal_division = []

		# UniqueMolecules division — E. coli names kept verbatim (engine
		# references these by exact attribute name); all empty for platelets
		self.unique_molecules_active_ribosome_division = []
		self.unique_molecules_RNA_division = []
		self.unique_molecules_domain_index_division = []
		self.unique_molecules_chromosomal_segment_division = []
