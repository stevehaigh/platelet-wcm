"""
SimulationData platelet-process namespace.
"""

from .stub import PlateletStub


class Process:
	"""Minimal process namespace for platelet simulation data."""

	def __init__(self, sim_data):
		self.platelet_stub = PlateletStub(sim_data)
