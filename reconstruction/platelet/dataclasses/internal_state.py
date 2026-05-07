"""
Internal state specifications for the platelet whole-cell model.

Provides the bulk_molecules and unique_molecule sub-objects that the engine's
BulkMolecules and UniqueMolecules states read during initialize().
"""

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray


# Molecule inventory for the v0.2 platelet model (calcium core).
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
#   [pl]  = plasmalemma (plasma membrane)
#
# IP3R, SERCA, PMCA, STIM1, and CaM are split into their kinetic sub-states so
# the CalciumDynamics ODE solver can integrate the 6-state Sneyd IP3R, the SERCA
# E1/E2 cycle, the 5-state Caride 2007 CaM-coupled PMCA scheme, the STIM1
# DTS-bound / free / dimer sensor cycle, and the 3-state CaM Ca²⁺-binding ladder.
# Sub-state initial counts are the Dolan & Diamond 2014 Table S1 representative
# configuration; protein totals match the per-protein totals in that table
# (1,328 IP3R; 11,892 SERCA; 769 PMCA; 4,265 STIM1; 1,447 Orai1; 20,481 CaM).
_MOLECULES = [
	# id                   mass_fg      initial_count  molecule_class
	# ── metabolites (concentrations derived from Purvis 2008 / Dolan 2014 / Sveshnikova 2025) ──
	('CA2_CYT[c]',        6.660e-8,    361,           'metabolite'),  # 100 nM × 6 fL
	('CA2_DTS[dts]',      6.660e-8,    38842,         'metabolite'),  # 250 µM × 4.3% × 6 fL
	('ATP[c]',            8.424e-7,    10_839_600,    'metabolite'),  # 3 mM × 6 fL (Holmsen 1979/1981)
	('ADP[c]',            7.096e-7,    1_083_960,     'metabolite'),  # 0.3 mM × 6 fL (ATP:ADP = 10:1)
	('PI[c]',             1.580e-7,    361_320,       'metabolite'),  # 100 µM × 6 fL inorganic phosphate
	('5HT[dg]',           2.927e-7,    3_500_000,     'metabolite'),  # serotonin; dense granule
	('ADP[dg]',           7.096e-7,    400_000,       'metabolite'),  # ADP; dense granule
	('IP3[c]',            6.977e-7,    181,           'metabolite'),  # 50 nM × 6 fL
	# ── proteins (copy numbers from Burkhart 2012 unless noted) ──
	('GP1BA[c]',          1.378e-4,    25_000,        'protein'),     # GpIbα; surface receptor
	('ITGA2B[c]',         2.149e-4,    80_000,        'protein'),     # αIIb integrin
	('ACTB[c]',           6.933e-5,    2_000_000,     'protein'),     # β-actin
	('FGA[ag]',           5.647e-4,    30_000,        'protein'),     # fibrinogen hexamer; alpha granule
	('SELP[ag]',          1.493e-4,    30_000,        'protein'),     # P-selectin; alpha granule
	# ── IP3R sub-states (type 2; mass = ITPR2 monomer × Burkhart total/Dolan total) ──
	('IP3R_n[dts]',       5.110e-4,    809,           'protein'),     # neutral
	('IP3R_o[dts]',       5.110e-4,    261,           'protein'),     # open
	('IP3R_a[dts]',       5.110e-4,    65,            'protein'),     # active (Ca²⁺-bound, conducting)
	('IP3R_i1[dts]',      5.110e-4,    167,           'protein'),     # inhibited-1 (Ca²⁺ at inhibitory site)
	('IP3R_i2[dts]',      5.110e-4,    25,            'protein'),     # inhibited-2
	('IP3R_s[dts]',       5.110e-4,    1,             'protein'),     # shut
	# ── SERCA3b sub-states (E1/E2 cycle; mass = ATP2A3 monomer) ──
	# E1 ↔ E1·Ca pre-equilibrated at cyt=100 nM (lab-book 2026-05-05 fix iii):
	# E1·Ca / E1 = k_bind_f · cyt² / k_bind_r = 1000 · 0.01 / 10 = 1.0
	# so the (E1 + E1·Ca = 5926) Dolan total splits ~2963 each. The Dolan
	# Table S1 values (5920, 6) are not at binding equilibrium for our ODE
	# and produced a spurious 118 k ions/s cyt → E1·Ca pulse on t=0.
	('SERCA_E1[dts]',     1.814e-4,    2_963,         'protein'),     # E1 (cytosol-facing, empty)
	('SERCA_E2[dts]',     1.814e-4,    5_927,         'protein'),     # E2 (DTS-facing, empty)
	('SERCA_E1Ca[dts]',   1.814e-4,    2_963,         'protein'),     # E1·2Ca²⁺
	('SERCA_E1PCa[dts]',  1.814e-4,    7,             'protein'),     # E1P·2Ca²⁺ (phosphorylated)
	('SERCA_E2PCa[dts]',  1.814e-4,    4,             'protein'),     # E2P·2Ca²⁺
	('SERCA_E2P[dts]',    1.814e-4,    28,            'protein'),     # E2P (Ca²⁺ released to DTS)
	# ── PMCA4b sub-states (Caride 2007 Table 3 5-state CaM-coupled scheme) ──
	# Basal path (steps 4–5): PMCA ⇌ PMCA·Ca → Ca²⁺_ex
	# CaM-activated path (steps 8–11): PMCA + Ca₄·CaM ⇌ Ca₄·CaM·PMCA
	#   → Ca₄·CaM·PMCA·Ca → Ca²⁺_ex; Ca₄·CaM·PMCA ⇌ PMCA·CaM + 4 Ca²⁺
	# Complex sub-states carry combined mass: PMCA (2.114e-4) + CaM (2.775e-5) = 2.391e-4 fg
	('PMCA[pl]',              2.114e-4,    765,     'protein'),  # PMCA free
	('PMCA_Ca[pl]',           2.114e-4,    4,       'protein'),  # PMCA·Ca²⁺ (basal)
	('Ca4_CaM_PMCA[pl]',      2.391e-4,    0,       'protein'),  # Ca₄·CaM·PMCA (CaM-activated, empty)
	('Ca4_CaM_PMCA_Ca[pl]',   2.391e-4,    0,       'protein'),  # Ca₄·CaM·PMCA·Ca²⁺
	('PMCA_CaM[pl]',          2.391e-4,    0,       'protein'),  # PMCA·CaM (deactivating)
	# ── Calmodulin sub-states (Caride 2007 Table 3 steps 6–7; mass = CaM monomer 16,706 Da) ──
	# CaM that is bound to PMCA is accounted for in the PMCA complex sub-states above;
	# these three states cover all free (unbound) CaM.
	# ICs are set at thermodynamic equilibrium with 100 nM free Ca²⁺ (Caride k6/k6r, k7/k7r).
	# Dolan Table S1 gave total CaM = 20,481 with sub-state breakdown (15, 1) that reflects
	# their original model without explicit CaM Ca²⁺-binding kinetics; those values are NOT
	# at equilibrium with 100 nM Ca²⁺ in our ODE and cause a spurious CaM-loading burst in
	# the first timestep.  The values below satisfy detailed balance exactly.
	('CaM_free[c]',       2.775e-5,    20_062,  'protein'),  # free CaM (equilibrated)
	('Ca2_CaM[c]',        2.775e-5,    200,     'protein'),  # Ca₂·CaM (N-lobe loaded; ~0.010 × free)
	('Ca4_CaM[c]',        2.775e-5,    219,     'protein'),  # Ca₄·CaM (fully loaded; ~1.10 × Ca₂·CaM)
	# ── STIM1 sub-states (sensor cycle; mass = STIM1 monomer) ──
	('STIM1_free[dts]',   1.285e-4,    438,           'protein'),     # free monomer (active sensor pool)
	('STIM1_Ca[dts]',     1.285e-4,    3_805,         'protein'),     # DTS-bound (inactive)
	('STIM1_dim[dts]',    2.570e-4,    11,            'protein'),     # 11 STIM1 dimer particles (Dolan Table S1)
	# ── Orai1 (CRAC channel pore-forming subunit) ──
	('ORAI1[pl]',         5.108e-5,    1_447,         'protein'),     # 30,768 Da; tetramerises to form ~360 channels
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
