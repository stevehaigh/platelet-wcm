"""
Internal state specifications for the platelet whole-cell model.

Provides the bulk_molecules and unique_molecule sub-objects that the engine's
BulkMolecules and UniqueMolecules states read during initialize().
"""

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray


# Placeholder molecule inventory. Real values from issue #18 (Burkhart proteome
# curation). Each entry: (molecule_id, mass_fg, initial_count, molecule_class)
# molecule_class: 'protein' | 'lipid' — used to route molecules to the correct
# decay process (RestingDecay targets proteins only; lipids get their own
# process once turnover data is available).
#
# Lipid count: platelet plasma membrane contains ~3×10^8 phospholipid molecules
# (van Meer 2008; Purvis 2008 surface area estimate). Using 3e8 as placeholder.
# Lipid half-life: membrane phospholipid turnover in resting anucleate platelets
# is poorly characterised — no active synthesis pathway running, so turnover is
# negligible on a 7-day timescale. Lipid decay deferred to a future process.
_PLACEHOLDER_MOLECULES = [
	('DUMMY_PROTEIN[c]', 1.0, 5000,    'protein'),  # Burkhart 2012 median
	('DUMMY_LIPID[e]',   1.0, 300_000_000, 'lipid'),  # van Meer 2008 order-of-magnitude
]

# Molecule IDs that RestingDecay should act on (proteins + mRNAs).
# Lipids, metabolites, and ions are excluded — they have different turnover
# mechanisms and will be handled by dedicated processes.
# Replace this list with a proper classification query when #18 lands.
PROTEIN_MOLECULE_IDS = np.array(
	[m[0] for m in _PLACEHOLDER_MOLECULES if m[3] == 'protein'],
	dtype='U50',
)


class _BulkMoleculesSpec:
	"""Mirrors the interface of reconstruction.ecoli.dataclasses.state.bulk_molecules."""

	def __init__(self):
		ids = np.array([m[0] for m in _PLACEHOLDER_MOLECULES], dtype=str)
		n_submass = 2  # protein (0), metabolite (1) — see SimulationDataPlatelet
		masses_raw = np.zeros((len(ids), n_submass))
		for i, (_, mass_fg, _, _) in enumerate(_PLACEHOLDER_MOLECULES):
			masses_raw[i, 0] = mass_fg  # protein submass column for all stubs

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
			[m[2] for m in _PLACEHOLDER_MOLECULES], dtype=np.int64)

		# Subset of molecule IDs that undergo protein-class decay.
		# Used by RestingDecay; replaces the all-molecules view.
		self.protein_molecule_names = PROTEIN_MOLECULE_IDS


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
