"""Tests for Tier-2 expression knockouts via RunConfig.count_overrides."""

import numpy as np
import pytest

from reconstruction.platelet import initialization as init
from reconstruction.platelet.run_config import RunConfig
from reconstruction.platelet.simulation_data import SimulationDataPlatelet
from wholecell.containers.bulk_objects_container import BulkObjectsContainer


def _sim_data_and_container():
	sim_data = SimulationDataPlatelet()
	ids = sim_data.internal_state.bulk_molecules.bulk_data['id']
	return sim_data, BulkObjectsContainer(ids), list(ids)


def test_run_config_default_count_overrides_is_empty():
	rc = RunConfig()
	assert rc.count_overrides == {}
	assert rc.to_metadata()['count_overrides'] == {}


def test_no_overrides_seeds_baseline_byte_identical():
	sim_data, container, _ids = _sim_data_and_container()
	init.initialize_bulk_molecules(container, sim_data)
	expected = np.array(
		sim_data.internal_state.bulk_molecules.initial_counts, dtype=np.int64)
	assert np.array_equal(container.counts(), expected)


def test_override_sets_a_species_count():
	sim_data, container, ids = _sim_data_and_container()
	target = ids[0]
	init.initialize_bulk_molecules(
		container, sim_data, count_overrides={target: 0})
	assert container.counts()[ids.index(target)] == 0
	# every other species is untouched
	baseline = np.array(
		sim_data.internal_state.bulk_molecules.initial_counts, dtype=np.int64)
	got = container.counts()
	assert np.array_equal(
		np.delete(got, ids.index(target)),
		np.delete(baseline, ids.index(target)))


def test_unknown_id_raises():
	sim_data, container, _ids = _sim_data_and_container()
	with pytest.raises(ValueError, match='unknown molecule id'):
		init.initialize_bulk_molecules(
			container, sim_data, count_overrides={'NOPE[c]': 0})
