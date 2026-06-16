"""Tests for the expression-knockout entity map (reconstruction.platelet.knockouts)."""

import pytest

from reconstruction.platelet import knockouts
from reconstruction.platelet.dataclasses._species_loader import load_species


def test_every_entity_id_is_a_real_species():
	real = {m[0] for m in load_species()}
	for name, ids in knockouts.ENTITIES.items():
		for species_id in ids:
			assert species_id in real, f'{name}: {species_id} not a species'


def test_expand_zeroes_all_sub_states():
	out = knockouts.expand(['PAR1 (thrombin)'])
	assert out == {
		'PAR1_inactive[pl]': 0,
		'PAR1_active[pl]': 0,
		'PAR1_internalized[pl]': 0,
	}


def test_expand_merges_multiple_entities():
	out = knockouts.expand(['P2Y1 (ADP receptor)', 'TP (TXA2 receptor)'])
	assert out['P2Y1_inactive[pl]'] == 0
	assert out['TP_active[pl]'] == 0
	assert all(v == 0 for v in out.values())


def test_expand_unknown_entity_raises():
	with pytest.raises(KeyError):
		knockouts.expand(['NotAnEntity'])
