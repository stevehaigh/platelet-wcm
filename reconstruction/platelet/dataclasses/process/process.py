"""
SimulationData platelet-process namespace.
"""

from .calcium_signalling import CalciumSignalling
from .granule_secretion import GranuleSecretion
from .integrin_activation import IntegrinActivation
from .resting_decay import RestingDecay
from .thromboxane_synthesis import ThromboxaneSynthesis


class Process:
	"""Process namespace for platelet simulation data."""

	def __init__(self, sim_data):
		self.resting_decay = RestingDecay(sim_data)
		self.calcium_signalling = CalciumSignalling(sim_data)
		self.granule_secretion = GranuleSecretion(sim_data)
		self.thromboxane_synthesis = ThromboxaneSynthesis(sim_data)
		self.integrin_activation = IntegrinActivation(sim_data)
