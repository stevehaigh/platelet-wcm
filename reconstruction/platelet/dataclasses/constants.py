"""
Constants for the platelet whole-cell model.
"""

from wholecell.utils import units


class Constants:
	"""Physical and biological constants used by the platelet simulation."""

	def __init__(self):
		self.n_avogadro = 6.022e23 / units.mol
