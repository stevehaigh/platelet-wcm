"""
Minimal platelet process used to prove the runtime contract.
"""

import numpy as np

import wholecell.processes.process


class PlateletStub(wholecell.processes.process.Process):
	"""No-op platelet process that exercises request and evolve hooks."""

	_name = 'PlateletStub'

	def __init__(self):
		super(PlateletStub, self).__init__()
		self.calculate_request_calls = 0
		self.evolve_state_calls = 0
		self.initialized = False

	def initialize(self, sim, sim_data):
		super(PlateletStub, self).initialize(sim, sim_data)
		self.molecule_names = sim_data.process.platelet_stub.molecule_names
		self.molecules = self.bulkMoleculesView(self.molecule_names)
		self.initialized = True

	def calculateRequest(self):
		self.calculate_request_calls += 1
		self.molecules.requestIs(
			np.zeros(self.molecule_names.size, dtype=np.int64))

	def evolveState(self):
		self.evolve_state_calls += 1
