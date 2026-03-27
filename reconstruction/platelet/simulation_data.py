"""
SimulationDataPlatelet

Minimal parameter object for the platelet whole-cell model.  Satisfies the
engine's interface contract so that BulkMolecules, UniqueMolecules, and
LocalEnvironment can all call initialize() without error.

Real biology is added incrementally:
  - Issue #18: Burkhart proteome — populates internal_state.bulk_molecules
  - Issue #24: Purvis calcium parameters
  - Issue #28: granule UniqueMolecule definitions
  - Issue #34: agonist-addition timeline replaces ExternalState stub
"""

from reconstruction.platelet.dataclasses.constants import Constants
from reconstruction.platelet.dataclasses.molecule_groups import MoleculeGroups
from reconstruction.platelet.dataclasses.internal_state import InternalState
from reconstruction.platelet.dataclasses.external_state import ExternalState


class SimulationDataPlatelet:
	"""Parameter object consumed by the platelet simulation engine."""

	def __init__(self):
		# ── Submass name → index mapping ──────────────────────────────────────
		# Subset of E. coli submass types relevant to platelets.
		# The engine uses this to construct per-submass mass arrays.
		self.submass_name_to_index = {
			'protein':    0,
			'metabolite': 1,
		}

		# ── Compartment abbreviation → index mapping ───────────────────────────
		# Platelet compartments:
		#   c  = cytoplasm
		#   dg = dense granule lumen
		#   ag = alpha-granule lumen
		#   m  = mitochondrial matrix
		#   e  = extracellular / open canalicular system
		self.compartment_abbrev_to_index = {
			'c':  0,
			'dg': 1,
			'ag': 2,
			'm':  3,
			'e':  4,
		}

		# ── Sub-objects read by engine states ──────────────────────────────────
		self.constants = Constants()
		self.molecule_groups = MoleculeGroups()
		self.internal_state = InternalState()
		self.external_state = ExternalState()
