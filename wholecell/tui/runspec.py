"""
Form schema and RunConfig assembly for the platelet TUI (phase P1).

`GROUPS` defines the editable run conditions, grouped by subsystem and mapped
1:1 onto `RunConfig` fields (see `reconstruction/platelet/run_config.py`).
Fields flagged `knockout=True` get an inline "KO" checkbox in the UI that, when
ticked, forces that knob to 0 (a knockout) at build time. `RUN_PARAMS` are the
two run-level arguments that are *not* RunConfig fields (length, seed).

`build_spec` turns the raw form strings into the JSON spec consumed by
`runscripts/manual/runFromConfig.py`.
"""

from __future__ import annotations

import json
from typing import Dict, Iterable, List, Tuple

from reconstruction.platelet.knockouts import ENTITIES as _ENTITIES
from reconstruction.platelet.knockouts import expand as _expand_entities

# Logical entities offered as expression (copy-number) knockouts in the UI.
KNOCKOUT_ENTITIES: List[str] = list(_ENTITIES)

# Run-level params (not RunConfig fields): (key, label, default, kind).
RUN_PARAMS: Tuple[Tuple[str, str, str, str], ...] = (
	('length_sec', 'Length (s)', '60', 'int'),
	('seed', 'Seed', '0', 'int'),
)

# Subsystem groups. Each field: (key, label, default, kind, knockout).
# `key` is a RunConfig field; `knockout` True → render an inline KO checkbox.
GROUPS: Tuple[Tuple[str, Tuple[Tuple[str, str, str, str, bool], ...]], ...] = (
	('Stimulus', (
		('thrombin_peak_nM',  'Thrombin (nM)',     '1.0',  'float', False),
		('adp_peak_uM',       'ADP (uM)',          '10.0', 'float', False),
		('atp_ex_peak_uM',    'ATP_ex (uM)',       '10.0', 'float', False),
		('ca_ex_mM',          'Ca2+_ex (mM)',      '1.2',  'float', False),
		('agonist_delay_s',   'Agonist delay (s)', '0',    'float', False),
	)),
	('Feedback loops', (
		('autocrine_adp_gain',  'Autocrine ADP gain',  '1.0', 'float', True),
		('autocrine_txa2_gain', 'Autocrine TXA2 gain', '1.0', 'float', True),
		('cox1_factor',         'COX-1 (aspirin)',     '1.0', 'float', True),
		('p2y12_block',         'P2Y12 block (clopidogrel)', '0.0', 'float', False),
	)),
	('Inhibitory / drugs (cAMP-raising)', (
		('pgi2_nM',    'PGI2 / iloprost (nM)',     '0.0', 'float', False),
		('forskolin',  'Forskolin (AC fold)',      '0.0', 'float', False),
		('pde3_block', 'PDE3 block (cilostazol)',  '0.0', 'float', False),
	)),
	('Pumps / brakes (scale)', (
		('pmca_kcat_scale',    'PMCA k_cat',         '1.0', 'float', True),
		('mcu_vmax_scale',     'MCU V_max',          '1.0', 'float', True),
		('mito_coupling_gain', 'MCU→IP3R coupling',  '1.0', 'float', True),
		('k_des_scale',        'P2Y1 desens',        '1.0', 'float', True),
		('k_plcb_phos_scale',  'PLCb phos',          '1.0', 'float', True),
		('integrin_act_scale', 'aIIbb3 activation',  '1.0', 'float', True),
		('rap1b_scale',        'Rap1b (CalDAG-GEFI)', '1.0', 'float', True),
	)),
)

# The three agonist peaks zeroed by the "run at rest" shorthand.
AGONIST_KEYS = ('thrombin_peak_nM', 'adp_peak_uM', 'atp_ex_peak_uM')

# key -> (label, kind, knockoutable); RunConfig field keys (excludes run params).
_LOOKUP: Dict[str, Tuple[str, str, bool]] = {}
for _key, _label, _default, _kind in RUN_PARAMS:
	_LOOKUP[_key] = (_label, _kind, False)
for _title, _fields in GROUPS:
	for _key, _label, _default, _kind, _ko in _fields:
		_LOOKUP[_key] = (_label, _kind, _ko)

RUN_CONFIG_KEYS: List[str] = [
	key for _title, fields in GROUPS for (key, *_rest) in fields]

# Knockoutable RunConfig fields (those with an inline KO checkbox).
KNOCKOUT_KEYS: List[str] = [
	key for _title, fields in GROUPS
	for (key, _l, _d, _k, ko) in fields if ko]

# Every editable input key (run params + RunConfig fields), for collection.
ALL_INPUT_KEYS: List[str] = [k for k, *_ in RUN_PARAMS] + RUN_CONFIG_KEYS


def field_default(key: str) -> str:
	"""Default string value for a field (for resetting the form)."""
	for key2, _label, default, _kind in RUN_PARAMS:
		if key2 == key:
			return default
	for _title, fields in GROUPS:
		for key2, _label, default, _kind, _ko in fields:
			if key2 == key:
				return default
	raise KeyError(key)


def diff_from_defaults(values: Dict[str, str], knockouts: Dict[str, bool],
		at_rest: bool) -> List[str]:
	"""Labels of knobs that differ from their schema defaults (for an indicator)."""
	changed = []
	for key in ALL_INPUT_KEYS:
		if str(values.get(key, '')).strip() != field_default(key):
			changed.append(_LOOKUP[key][0])
	for key in KNOCKOUT_KEYS:
		if knockouts.get(key):
			changed.append(_LOOKUP[key][0] + ' KO')
	if at_rest:
		changed.append('at rest')
	return changed


def build_spec(values: Dict[str, str], knockouts: Dict[str, bool],
		at_rest: bool, entity_knockouts: Iterable[str] = ()) -> dict:
	"""Assemble a {length_sec, seed, run_config} spec from raw form strings.

	`values` maps field key → raw string; `knockouts` maps a knockoutable key →
	whether its KO checkbox is ticked; `entity_knockouts` are logical-entity
	names (from `KNOCKOUT_ENTITIES`) to zero out as expression knockouts (Tier
	2 — fed to `RunConfig.count_overrides`). Raises ValueError (naming the
	field) on a missing/non-numeric/negative entry.
	"""

	def parse(key: str):
		label, kind, _ko = _LOOKUP[key]
		raw = str(values.get(key, '')).strip()
		if raw == '':
			raise ValueError(f'{label} is empty')
		try:
			val = float(raw)
		except ValueError:
			raise ValueError(f'{label} is not a number: {raw!r}')
		if val < 0:
			raise ValueError(f'{label} must be ≥ 0')
		return int(val) if kind == 'int' else val

	run_config: Dict[str, object] = {}
	for key in RUN_CONFIG_KEYS:
		_label, _kind, knockoutable = _LOOKUP[key]
		if knockoutable and knockouts.get(key):
			run_config[key] = 0.0
		else:
			run_config[key] = parse(key)
	if at_rest:
		for key in AGONIST_KEYS:
			run_config[key] = 0.0
	run_config['live'] = True
	run_config['count_overrides'] = _expand_entities(entity_knockouts)

	return {
		'length_sec': parse('length_sec'),
		'seed': parse('seed'),
		'run_config': run_config,
	}


def write_spec(spec: dict, path: str) -> None:
	"""Write a run spec to JSON (also serves as a saved preset)."""
	with open(path, 'w') as handle:
		json.dump(spec, handle, indent=2)
