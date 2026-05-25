"""TOML-backed parameter loader for the platelet calcium-signalling
module.

Reads rate constants from ``reports/params/calcium-v<N>.toml`` and
returns a nested dict. The constants in ``calcium_signalling.py`` are
populated from this loader at import time, so no call site downstream
of the module-level ``K_*`` dicts changes (issue #32, Phase 2).

Versioning
----------
The TOML filename suffix (``v0.5``, ``v0.6``, …) bumps when the schema
shape changes; value-only updates stay on the current file. Pin the
default in ``_DEFAULT_VERSION`` so a stale Python branch can still
locate the matching TOML.

Float semantics
---------------
TOML 1.0 floats are IEEE 754 binary64 — the same representation Python
uses for ``float`` literals. ``0.13`` parsed from TOML is bit-equal to
``0.13`` written in Python source, so externalising a value is
byte-identical at the simulation-output level.
"""

import os
import tomllib
from typing import Any


_DEFAULT_VERSION = 'v0.5'

# `_params_loader.py` lives at
#   reconstruction/platelet/dataclasses/process/_params_loader.py
# so the repo root is four directories up.
_PARAMS_DIR = os.path.normpath(os.path.join(
	os.path.dirname(os.path.abspath(__file__)),
	'..', '..', '..', '..',
	'reports', 'params',
))


def load_calcium_kinetics(version: str = _DEFAULT_VERSION) -> dict[str, Any]:
	"""Read ``reports/params/calcium-<version>.toml`` and return the
	parsed structure as a nested dict.

	Caller is expected to extract the section it needs, e.g.
	``load_calcium_kinetics()['ip3r']['k_dyk']`` for the IP3R deYoung-
	Keizer parameters.
	"""
	path = os.path.join(_PARAMS_DIR, f'calcium-{version}.toml')
	with open(path, 'rb') as f:
		return tomllib.load(f)
