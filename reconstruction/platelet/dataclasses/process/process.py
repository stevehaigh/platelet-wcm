"""
SimulationData platelet-process namespace.
"""

from .resting_decay import RestingDecay
from .stub import PlateletStub


class Process:
	"""Minimal process namespace for platelet simulation data."""

	def __init__(self, sim_data):
		self.platelet_stub = PlateletStub(sim_data)
		self.resting_decay = RestingDecay(sim_data)
