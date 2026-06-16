"""
Initialization helpers for the platelet simulation scaffold.
"""

import numpy as np


def initialize_bulk_molecules(bulk_container, sim_data, count_overrides=None):
	"""Seed the bulk-molecule counts from the sim_data resting baseline.

	``count_overrides`` (a per-run ``{molecule_id: count}`` mapping from
	``RunConfig``) overrides individual resting copy numbers — e.g. zeroing an
	entity's sub-states is an expression knockout. Empty/None → the unmodified
	baseline (byte-identical). An unknown id raises (surfaces typos early).
	"""
	counts = np.array(
		sim_data.internal_state.bulk_molecules.initial_counts, dtype=np.int64)
	if count_overrides:
		ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
		index = {mol_id: i for i, mol_id in enumerate(ids)}
		for mol_id, count in count_overrides.items():
			if mol_id not in index:
				raise ValueError(
					f'count_overrides: unknown molecule id {mol_id!r}')
			counts[index[mol_id]] = int(count)
	bulk_container.countsIs(counts)


def initialize_unique_molecules(unique_container, sim_data):
	"""Keep the unique-molecule state explicit even while the scaffold is empty."""
	del unique_container, sim_data
