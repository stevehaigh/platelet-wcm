"""
Internal state specifications for the platelet whole-cell model.

Provides the bulk_molecules and unique_molecule sub-objects that the engine's
BulkMolecules and UniqueMolecules states read during initialize().
"""

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray


# Placeholder molecule inventory. Real copy numbers will be populated in
# issue #18 (Burkhart proteome curation). For now, a tiny inventory is enough
# to prove the simulation runtime can initialize and step cleanly.
_PLACEHOLDER_MOLECULES = [
	('DUMMY_PROTEIN[c]', 1.0, 1),  # cytoplasm
	('DUMMY_LIPID[e]', 1.0, 1),  # extracellular
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
		for i, (_, mass_fg, _) in enumerate(_PLACEHOLDER_MOLECULES):
			masses_raw[i, 0] = mass_fg  # protein column

		bulk_data = np.zeros(
			len(ids),
			dtype=[
				('id', 'U50'),
				('mass', '{}f8'.format(n_submass)),
				],
			)

		bulk_data['id'] = ids

		# Engine reads: bulk_data['mass'].asNumber(fg/mol), so mass must carry
		# units of fg/mol (multiply per-molecule fg by N_A).
		n_avogadro = 6.022e23
		bulk_data['mass'] = masses_raw * n_avogadro
		self.bulk_data = UnitStructArray(bulk_data, {
			'id': None,
			'mass': units.fg / units.mol,
			})
		self.initial_counts = np.array(
			[count for _, _, count in _PLACEHOLDER_MOLECULES], dtype=np.int64)


class _UniqueMoleculeSpec:
	"""Mirrors the interface of reconstruction.ecoli.dataclasses.state.unique_molecules."""

	def __init__(self):
		# No unique molecules in the v0.1 stub.
		# Granules will be added as UniqueMolecules in issue #28.
		self.unique_molecule_definitions = {}
		unique_molecule_masses = np.zeros(
			0,
			dtype=[
				('id', 'U50'),
				('mass', '2f8'),
				],
			)
		self.unique_molecule_masses = UnitStructArray(unique_molecule_masses, {
			'id': None,
			'mass': units.fg / units.mol,
			})


class InternalState:
	"""Top-level internal_state namespace read by the engine states."""

	def __init__(self):
		self.bulk_molecules = _BulkMoleculesSpec()
		self.unique_molecule = _UniqueMoleculeSpec()
