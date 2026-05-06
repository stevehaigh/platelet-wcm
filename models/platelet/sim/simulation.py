from wholecell.sim.simulation import Simulation

from wholecell.states.bulk_molecules import BulkMolecules
from wholecell.states.unique_molecules import UniqueMolecules
from wholecell.states.local_environment import LocalEnvironment

from models.platelet.listeners.calcium_trace import CalciumTrace
from models.platelet.listeners.mass import Mass
from models.platelet.processes.calcium_dynamics import CalciumDynamics
from models.platelet.processes.resting_decay import RestingDecay
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
		),
	)

	_listenerClasses = (Mass, CalciumTrace)

	_hookClasses = ()

	_initialConditionsFunction = calcInitialConditions

	_shellColumnHeaders = (
		'Time (s)',
		)


def platelet_simulation(**options):
	"""Instantiate a platelet simulation using the scaffold runtime."""
	return PlateletSimulation(**options)
