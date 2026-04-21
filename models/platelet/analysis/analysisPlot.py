"""Minimal matplotlib-based analysis base for platelet plots."""

import abc
import os

import matplotlib as mp
import numpy as np

from wholecell.utils import filepath as fp


class AnalysisPlot(metaclass=abc.ABCMeta):
	"""Small analysis base that avoids importing E. coli-specific modules."""

	_suppress_numpy_warnings = False

	def __init__(self, cpus=0):
		self.cpus = cpus
		self._axeses = {}
		self.ap = None

	@abc.abstractmethod
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		raise NotImplementedError

	def plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		def do_plot():
			self.do_plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata)

		if not os.path.isdir(inputDir):
			raise RuntimeError(
				'Input directory ({}) does not currently exist.'.format(inputDir))
		fp.makedirs(plotOutDir)

		with mp.rc_context():
			if self._suppress_numpy_warnings:
				with np.errstate(divide='ignore'), np.errstate(invalid='ignore'):
					do_plot()
			else:
				do_plot()

		self._axeses = {}

	@classmethod
	def main(cls, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile=None, metadata=None, cpus=0, analysis_paths=None):
		instance = cls(cpus)
		instance.ap = analysis_paths
		instance.plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)
