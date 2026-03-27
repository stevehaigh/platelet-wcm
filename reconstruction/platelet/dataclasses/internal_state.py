"""
Internal state specifications for the platelet whole-cell model.

Provides the bulk_molecules and unique_molecule sub-objects that the engine's
BulkMolecules and UniqueMolecules states read during initialize().
"""

import numpy as np

from wholecell.utils import units


# Placeholder molecule inventory.  Real copy numbers will be populated in
# issue #18 (Burkhart proteome curation).  For now, one dummy molecule per
# compartment is enough to prove the engine initialises without error.
_PLACEHOLDER_MOLECULES = [
	('DUMMY_PROTEIN[c]',  1.0),   # cytoplasm
	('DUMMY_LIPID[e]',    1.0),   # extracellular
]


class _BulkMoleculesSpec:
	"""Mirrors the interface of reconstruction.ecoli.dataclasses.state.bulk_molecules."""

	def __init__(self):
		ids = np.array([m[0] for m in _PLACEHOLDER_MOLECULES], dtype=str)
		# Mass array shape: (n_molecules, n_submass_types).
		# Submass indices match SimulationDataPlatelet.submass_name_to_index.
		# All mass attributed to 'protein' (index 0) for the stub.
		n_submass = 2  # protein, metabolite (see SimulationDataPlatelet)
		masses_raw = np.zeros((len(ids), n_submass))
		for i, (_, mass_fg) in enumerate(_PLACEHOLDER_MOLECULES):
			masses_raw[i, 0] = mass_fg  # protein column

		# Engine reads: bulk_data['id'] and bulk_data['mass'].asNumber(fg/mol)
		# so mass must carry units of fg/mol (multiply per-molecule fg by N_A).
		n_avogadro = 6.022e23
		self.bulk_data = {
			'id':   ids,
			'mass': masses_raw * n_avogadro * units.fg / units.mol,
		}


class _UniqueMoleculeSpec:
	"""Mirrors the interface of reconstruction.ecoli.dataclasses.state.unique_molecules."""

	def __init__(self):
		# No unique molecules in the v0.1 stub.
		# Granules will be added as UniqueMolecules in issue #28.
		self.unique_molecule_definitions = {}
		self.unique_molecule_masses = {
			'id':   np.array([], dtype=str),
			'mass': np.array([]) * units.fg / units.mol,
		}


class InternalState:
	"""Top-level internal_state namespace read by the engine states."""

	def __init__(self):
		self.bulk_molecules = _BulkMoleculesSpec()
		self.unique_molecule = _UniqueMoleculeSpec()
