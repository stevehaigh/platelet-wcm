"""Common code for platelet single-run analysis plots."""

from models.platelet.analysis import analysisPlot
from runscripts.manual import analysisPlatelet


class SingleAnalysisPlot(analysisPlot.AnalysisPlot):
	"""Abstract base class for platelet single-run analysis plots."""

	def cli(self):
		"""Run this platelet plot from the command line."""
		script = analysisPlatelet.AnalysisPlatelet(analysis_plotter=self)
		script.cli()
