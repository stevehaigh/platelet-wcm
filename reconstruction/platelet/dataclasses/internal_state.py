"""
Internal state specifications for the platelet whole-cell model.

Provides the bulk_molecules and unique_molecule sub-objects that the engine's
BulkMolecules and UniqueMolecules states read during initialize().
"""

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray


# Minimal molecule inventory for the v0.1 platelet stub.
# Raw data and citations: reconstruction/platelet/raw_data/molecules.tsv
#
# Each entry: (molecule_id, mass_fg, initial_count, molecule_class)
#   molecule_class: 'protein' | 'metabolite'
#   mass_fg: per-molecule mass in femtograms = mw_da × 1.661e-9
#   initial_count: resting-state copy number per platelet (see TSV for sources)
#
# Submass routing (matches SimulationDataPlatelet.submass_name_to_index):
#   'protein'    → submass column 0
#   'metabolite' → submass column 1
#
# Compartments:
#   [c]   = cytoplasm     [dts] = dense tubular system (Ca2+ store)
#   [dg]  = dense granule [ag]  = alpha granule
_MOLECULES = [
	# id              mass_fg      initial_count  molecule_class
	# ── metabolites (concentrations derived from Purvis 2008 / Dolan 2014 / Sveshnikova 2025) ──
	('CA2_CYT[c]',   6.660e-8,    361,           'metabolite'),  # 100 nM × 6 fL
	('CA2_DTS[dts]', 6.660e-8,    38842,         'metabolite'),  # 250 µM × 4.3% × 6 fL
	('ATP[c]',       8.424e-7,    3_613_200,     'metabolite'),  # 1 mM × 6 fL
	('ADP[c]',       7.096e-7,    361_320,       'metabolite'),  # 0.1 mM × 6 fL
	('5HT[dg]',      2.927e-7,    3_500_000,     'metabolite'),  # serotonin; dense granule
	('ADP[dg]',      7.096e-7,    400_000,       'metabolite'),  # ADP; dense granule
	('IP3[c]',       6.977e-7,    181,           'metabolite'),  # 50 nM × 6 fL
	# ── proteins (copy numbers from Burkhart 2012 unless noted) ──
	('GP1BA[c]',     1.378e-4,    25_000,        'protein'),   # GpIbα; surface receptor
	('ITGA2B[c]',    2.149e-4,    80_000,        'protein'),   # αIIb integrin
	('ACTB[c]',      6.933e-5,    2_000_000,     'protein'),   # β-actin
	('FGA[ag]',      5.647e-4,    30_000,        'protein'),   # fibrinogen hexamer; alpha granule
	('SELP[ag]',     1.493e-4,    30_000,        'protein'),   # P-selectin; alpha granule
	('ITPR2[c]',     5.110e-4,    1_700,         'protein'),   # IP3 receptor type 2
	('ATP2A3[c]',    1.814e-4,    16_300,        'protein'),   # SERCA3 Ca2+-ATPase
	('STIM1[c]',     1.285e-4,    7_400,         'protein'),   # STIM1 Ca2+ sensor
]

# Submass column index for each class (mirrors SimulationDataPlatelet).
_SUBMASS_COL = {'protein': 0, 'metabolite': 1}

# IDs of molecules that undergo protein-class decay (RestingDecay).
# Metabolites and ions are excluded — different turnover mechanisms.
PROTEIN_MOLECULE_IDS = np.array(
	[m[0] for m in _MOLECULES if m[3] == 'protein'],
	dtype='U50',
)


class _BulkMoleculesSpec:
	"""Mirrors the interface of reconstruction.ecoli.dataclasses.state.bulk_molecules."""

	def __init__(self):
		ids = np.array([m[0] for m in _MOLECULES], dtype=str)
		n_submass = 2  # protein (0), metabolite (1) — see SimulationDataPlatelet
		masses_raw = np.zeros((len(ids), n_submass))
		for i, (_, mass_fg, _, molecule_class) in enumerate(_MOLECULES):
			masses_raw[i, _SUBMASS_COL[molecule_class]] = mass_fg

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
			[m[2] for m in _MOLECULES], dtype=np.int64)

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
