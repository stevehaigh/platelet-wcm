"""
Internal state specifications for the platelet whole-cell model.

Provides the bulk_molecules and unique_molecule sub-objects that the engine's
BulkMolecules and UniqueMolecules states read during initialize().

Species inventory
-----------------
The molecule inventory (id, mass, initial count, class for all 66 species)
lives in ``reports/params/species-v0.6.tsv`` (issue #32 Phase 2, species
side). It is loaded at import time via ``_species_loader.load_species()``
and exposed below as ``_MOLECULES`` — bit-equal to the pre-refactor
in-source list.

Each row of the TSV is ``(id, mass_fg, initial_count, molecule_class)``:
  - id              : str — molecule identifier with compartment tag
  - mass_fg         : float — per-molecule mass in femtograms (mw_da × 1.661e-9)
  - initial_count   : int — resting-state copy number per platelet
  - molecule_class  : 'protein' | 'metabolite'

Submass routing (matches SimulationDataPlatelet.submass_name_to_index):
  'protein'    → submass column 0
  'metabolite' → submass column 1

Compartments:
  [c]   = cytoplasm     [dts] = dense tubular system (Ca²⁺ store)
  [dg]  = dense granule [ag]  = alpha granule
  [pl]  = plasmalemma (plasma membrane)
  [m]   = mitochondrial matrix

Sub-state breakdowns
--------------------
IP3R, SERCA, PMCA, STIM1, and CaM are split into kinetic sub-states so
CalciumDynamics can integrate the Li-Rinzel IP3R inactivation, the SERCA
E1/E2 cycle, the Caride 2007 CaM-coupled PMCA scheme, the STIM1 sensor
cycle (free / DTS-bound / dimer), and the CaM Ca²⁺-binding ladder. Sub-
state initial counts are pre-equilibrated at the Dolan resting state
(cyt = 100 nM, DTS = 250 µM).

Derivations for the rate constants that pin those equilibria are
documented in two complementary places (the kinetics-as-data refactor,
issue #32, externalised most rate constants into TOML):

  - ``reports/params/calcium-v0.6.toml`` — source-of-truth for the
    rate constants themselves, with literature attribution and
    calibration-coupling notes inline.
  - ``reports/design/kinetics-v0.5-review.pdf`` — clickable PDF
    rendering of the above; regenerate via
    ``runscripts/manual/buildKineticsReview.py``.
  - ``reconstruction/platelet/dataclasses/process/calcium_signalling.py``
    — Python module-level names (``K_SERCA`` etc.) that the engine
    reads; today these are populated from the TOML at import time.

DTS luminal buffers (CALR, HSP90B1, BiP, CREC) carry ~95–99 % of the
DTS Ca²⁺ at rest; their free/bound sub-state counts are equilibrium-
derived from the corresponding ``K_*`` Kd at DTS [Ca²⁺] = 250 µM.

PI cycle / GPCR cascade sub-states (Phase 4 / issue #31 + Phase v0.4 /
issue #9) are at the tonic resting equilibrium of their respective
rate constants — see calcium_signalling.py for the Gq tonic-floor /
PLCβ activation fraction / IP3 baseline calibration chain.
"""

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.platelet.dataclasses._species_loader import load_species


# Externalised species inventory — see module docstring above and the
# v0.6 TSV for the source-of-truth values. `_MOLECULES` retains the
# legacy 4-tuple shape so all downstream call sites (initial_counts
# array, mass matrix, PROTEIN_MOLECULE_IDS) are unchanged.
_MOLECULES = load_species()

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
