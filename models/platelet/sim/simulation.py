from wholecell.sim.simulation import Simulation

from reconstruction.platelet.run_config import RunConfig
from wholecell.states.bulk_molecules import BulkMolecules
from wholecell.states.unique_molecules import UniqueMolecules
from wholecell.states.local_environment import LocalEnvironment

from models.platelet.listeners.calcium_trace import CalciumTrace
from models.platelet.listeners.mass import Mass
from models.platelet.listeners.secretion_trace import SecretionTrace
from models.platelet.listeners.thromboxane_trace import ThromboxaneTrace
from models.platelet.processes.calcium_dynamics import CalciumDynamics
from models.platelet.processes.granule_secretion import GranuleSecretion
from models.platelet.processes.resting_decay import RestingDecay
from models.platelet.processes.thromboxane_synthesis import ThromboxaneSynthesis
from models.platelet.sim.initial_conditions import calcInitialConditions


class PlateletSimulation(Simulation):
	_internalStateClasses = (
		BulkMolecules,
		UniqueMolecules,
		)

	_externalStateClasses = (
		LocalEnvironment,
		)

	_processClasses = (
		(
			RestingDecay,
			CalciumDynamics,
			GranuleSecretion,
			ThromboxaneSynthesis,
		),
	)

	_listenerClasses = (Mass, CalciumTrace, SecretionTrace, ThromboxaneTrace)

	_hookClasses = ()

	_initialConditionsFunction = calcInitialConditions

	_shellColumnHeaders = (
		'Time (s)',
		)

	def __init__(self, run_config=None, **kwargs):
		# Per-run conditions (extracellular Ca²⁺, agonist peaks, feedback gains,
		# perturbation scales). Stored before super().__init__() runs
		# _initialize(), so processes/listeners can read sim.run_config in
		# their own initialize(). Replaces the old module-global monkeypatching.
		self.run_config = run_config if run_config is not None else RunConfig()
		super().__init__(**kwargs)


def platelet_simulation(**options):
	"""Instantiate a platelet simulation using the scaffold runtime."""
	return PlateletSimulation(**options)
