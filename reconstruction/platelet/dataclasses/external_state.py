"""
External state stub for the platelet whole-cell model.

LocalEnvironment expects media/timeline machinery from sim_data.external_state.
For the v0.1 stub, platelets have no nutrient environment — this will be
replaced with an agonist-addition timeline in issue #34.

The stub provides the minimum interface to keep LocalEnvironment from raising
an AttributeError, while doing nothing.
"""


_STUB_MEDIA_ID = 'platelet_resting'
_STUB_MEDIA = {_STUB_MEDIA_ID: {}}   # empty environment — no exchange molecules


def empty_exchange_data_from_concentrations(concentrations):
	"""Return an empty exchange-data payload for the platelet stub."""
	del concentrations
	return {
		'importUnconstrainedExchangeMolecules': [],
		'importConstrainedExchangeMolecules': {},
		}


class _MakeMedia:
	"""Stub for the make_media object LocalEnvironment calls."""

	def make_timeline(self, timeline):
		"""Return a single-entry timeline pointing at the stub resting media."""
		return [(0, _STUB_MEDIA_ID)]


class ExternalState:
	"""Stub external_state namespace consumed by LocalEnvironment.initialize()."""

	def __init__(self):
		self.make_media = _MakeMedia()

		# No variant-specified timeline override
		self.current_timeline_id = None

		# Single stub media entry so LocalEnvironment can look up the initial media
		self.saved_timelines = {}
		self.saved_media = _STUB_MEDIA

		# No exchange reactions for a resting platelet stub
		self.exchange_data_from_concentrations = (
			empty_exchange_data_from_concentrations)
		self.exchange_to_env_map = {}
		self.import_constraint_threshold = 0.0
