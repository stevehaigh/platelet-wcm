"""
Minimal platelet process data for runtime scaffold wiring.
"""

import numpy as np


class PlateletStub:
	"""Expose the smallest process-facing molecule list for the scaffold."""

	def __init__(self, sim_data):
		self.molecule_names = np.array(
			sim_data.internal_state.bulk_molecules.bulk_data['id'], dtype='U')
