"""
Initialization helpers for the platelet simulation scaffold.
"""

import numpy as np


def initialize_bulk_molecules(bulk_container, sim_data):
	"""Load the placeholder platelet bulk counts into the simulation state."""
	bulk_container.countsIs(np.array(
		sim_data.internal_state.bulk_molecules.initial_counts, dtype=np.int64))


def initialize_unique_molecules(unique_container, sim_data):
	"""Keep the unique-molecule state explicit even while the scaffold is empty."""
	del unique_container, sim_data
