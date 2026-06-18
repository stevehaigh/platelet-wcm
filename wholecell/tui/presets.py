"""
Built-in and user-saved presets for the platelet TUI (phase P1).

A preset captures the *form state* — field values, knockout-checkbox states,
and the at-rest switch — not a run spec; loading one repopulates the form
deterministically (overrides applied over the schema defaults). Built-ins live
here; user presets are JSON files under a presets directory (overridable via
`$PLATELET_TUI_PRESETS_DIR`, default `~/.config/platelet-wcm/tui-presets/`).

The built-ins mirror the webapp's four conditions (Agonist transient / 60 s
settle / EDTA / Resting) plus two curated knockouts (Aspirin / Glanzmann).
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Tuple

from wholecell.tui import runspec

# name -> partial overrides; any field not named falls back to its schema
# default. `values` are field strings, `knockouts` flips KO checkboxes.
BUILTINS: Dict[str, dict] = {
	'Defaults': {},
	'Agonist transient': {'values': {'length_sec': '200'}},
	'Agonist (60 s settle)': {
		'values': {'length_sec': '260', 'agonist_delay_s': '60'}},
	'EDTA (no Ca_ex)': {'values': {'length_sec': '200', 'ca_ex_mM': '0.0'}},
	'Resting': {'values': {'length_sec': '300'}, 'at_rest': True},
	'Aspirin (COX-1 KO)': {
		'values': {'length_sec': '200'}, 'knockouts': {'cox1_factor': True}},
	'Glanzmann (aIIbb3 KO)': {
		'values': {'length_sec': '200'},
		'knockouts': {'integrin_act_scale': True}},
}

# (values, knockouts, at_rest) — a full, ready-to-apply form state.
PresetState = Tuple[Dict[str, str], Dict[str, bool], bool]


def _from_overrides(overrides: dict) -> PresetState:
	"""Expand partial overrides into a full form state over the defaults."""
	values = {key: runspec.field_default(key) for key in runspec.ALL_INPUT_KEYS}
	values.update({k: str(v) for k, v in overrides.get('values', {}).items()})
	knockouts = {key: False for key in runspec.KNOCKOUT_KEYS}
	knockouts.update(overrides.get('knockouts', {}))
	return values, knockouts, bool(overrides.get('at_rest', False))


def presets_dir() -> str:
	"""Directory holding user presets (env-overridable for tests)."""
	return os.environ.get(
		'PLATELET_TUI_PRESETS_DIR',
		os.path.expanduser('~/.config/platelet-wcm/tui-presets'))


def list_user_presets() -> List[str]:
	"""Saved user-preset names (filenames without `.json`), sorted."""
	directory = presets_dir()
	if not os.path.isdir(directory):
		return []
	return sorted(
		name[:-5] for name in os.listdir(directory) if name.endswith('.json'))


def all_names() -> List[str]:
	"""Built-in names followed by user-preset names."""
	return list(BUILTINS) + list_user_presets()


def resolve(name: str) -> PresetState:
	"""Form state for a preset name (a user preset shadows a built-in)."""
	user_path = os.path.join(presets_dir(), name + '.json')
	if os.path.isfile(user_path):
		with open(user_path) as handle:
			return _from_overrides(json.load(handle))
	return _from_overrides(BUILTINS.get(name, {}))


def _safe_name(name: str) -> str:
	"""Filesystem-safe preset name (also its canonical display name)."""
	safe = ''.join(
		c if c.isalnum() or c in ' -_' else '_' for c in name).strip()
	return safe or 'preset'


def save(name: str, values: Dict[str, str], knockouts: Dict[str, bool],
		at_rest: bool) -> str:
	"""Write the form state as a user preset; return its canonical name."""
	safe = _safe_name(name)
	directory = presets_dir()
	os.makedirs(directory, exist_ok=True)
	payload = {
		'values': values,
		'knockouts': {k: v for k, v in knockouts.items() if v},
		'at_rest': at_rest,
	}
	with open(os.path.join(directory, safe + '.json'), 'w') as handle:
		json.dump(payload, handle, indent=2)
	return safe
