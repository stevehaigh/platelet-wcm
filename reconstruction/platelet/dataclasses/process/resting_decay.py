"""
Resting-state protein decay dataclass for the platelet whole-cell model.

Anucleate platelets cannot replenish proteins via transcription, so all
proteins decay with a characteristic half-life (~7 days, Burkhart 2012).

Molecules managed by the CalciumDynamics ODE are excluded — their counts
evolve under the ODE each timestep and must not be independently decayed.
"""

import numpy as np

from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	MOLECULE_NAMES as _CALCIUM_MOLECULE_NAMES,
)

# Protein half-life in seconds (7 days).
# Source: Burkhart 2012 (platelet proteome turnover study).
_PROTEIN_HALF_LIFE_SEC = 7 * 24 * 3600  # 604800 s

_CALCIUM_SET = frozenset(_CALCIUM_MOLECULE_NAMES)


class RestingDecay:
	"""Parameters for the RestingDecay process."""

	def __init__(self, sim_data):
		all_proteins = sim_data.internal_state.bulk_molecules.protein_molecule_names
		self.molecule_names = np.array(
			[m for m in all_proteins if m not in _CALCIUM_SET],
			dtype=all_proteins.dtype,
		)
		self.protein_half_life = _PROTEIN_HALF_LIFE_SEC
