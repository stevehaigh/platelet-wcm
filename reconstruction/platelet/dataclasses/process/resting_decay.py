"""
Resting-state protein decay dataclass for the platelet whole-cell model.

Anucleate platelets cannot replenish proteins via transcription, so all
proteins decay with a characteristic half-life (~7 days, Burkhart 2012).

Molecules managed by the CalciumDynamics ODE are excluded — their counts
evolve under the ODE each timestep and must not be independently decayed.
Species relocated by GranuleSecretion (its cargo sources + destinations) are
excluded for the same reason; the non-protein members are harmless no-ops here.
"""

import numpy as np

from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	MOLECULE_NAMES as _CALCIUM_MOLECULE_NAMES,
)
from reconstruction.platelet.dataclasses.process.granule_secretion import (
	SECRETION_MANAGED_SPECIES as _SECRETION_SPECIES,
)

# Protein half-life in seconds (7 days).
# Source: Burkhart 2012 (platelet proteome turnover study).
_PROTEIN_HALF_LIFE_SEC = 7 * 24 * 3600  # 604800 s

_EXCLUDED_SET = frozenset(_CALCIUM_MOLECULE_NAMES) | frozenset(_SECRETION_SPECIES)


class RestingDecay:
	"""Parameters for the RestingDecay process."""

	def __init__(self, sim_data):
		all_proteins = sim_data.internal_state.bulk_molecules.protein_molecule_names
		self.molecule_names = np.array(
			[m for m in all_proteins if m not in _EXCLUDED_SET],
			dtype=all_proteins.dtype,
		)
		self.protein_half_life = _PROTEIN_HALF_LIFE_SEC
