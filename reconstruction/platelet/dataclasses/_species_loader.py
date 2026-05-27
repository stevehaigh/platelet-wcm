"""TSV-backed species loader for the platelet model.

Reads the molecule inventory from ``reports/params/species-v<N>.tsv`` and
returns a list of 4-tuples ``(id, mass_fg, initial_count, molecule_class)``
shaped to drop in for the legacy in-source ``_MOLECULES`` list in
``internal_state.py`` (issue #32 Phase 2 — species side).

Versioning matches the kinetics loader: the filename suffix (``v0.5``,
``v0.6``, …) bumps when the column schema changes; value-only updates
stay on the current file.

Float semantics
---------------
``repr(x)`` in Python guarantees float round-trip — ``float(repr(x)) is x``
at the bit level for IEEE 754 binary64. The TSV dump used ``repr()`` so
reload via ``float(...)`` is bit-equal to the original Python literals.
``int(...)`` round-trips trivially.
"""

import csv
import os
from typing import List, Tuple


_DEFAULT_VERSION = 'v0.5'

# This file lives at:
#   reconstruction/platelet/dataclasses/_species_loader.py
# Repo root is three levels up.
_PARAMS_DIR = os.path.normpath(os.path.join(
	os.path.dirname(os.path.abspath(__file__)),
	'..', '..', '..',
	'reports', 'params',
))


SpeciesRow = Tuple[str, float, int, str]


def load_species(version: str = _DEFAULT_VERSION) -> List[SpeciesRow]:
	"""Read ``reports/params/species-<version>.tsv`` and return its rows
	as ``(id, mass_fg, initial_count, molecule_class)`` 4-tuples in file
	order — bit-equal to the pre-refactor in-source ``_MOLECULES`` list.
	"""
	path = os.path.join(_PARAMS_DIR, f'species-{version}.tsv')
	rows: List[SpeciesRow] = []
	with open(path, newline='', encoding='utf-8') as f:
		reader = csv.DictReader(f, delimiter='\t')
		for r in reader:
			rows.append((
				r['id'],
				float(r['mass_fg']),
				int(r['initial_count']),
				r['molecule_class'],
			))
	return rows
