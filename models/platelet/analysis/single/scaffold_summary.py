"""Summarize the current platelet scaffold run."""

import os

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from models.platelet.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""Plot placeholder bulk counts and runtime timings for the scaffold."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile

		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		bulk_reader = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
		eval_reader = TableReader(os.path.join(simOutDir, 'EvaluationTime'))

		time = main_reader.readColumn('time').flatten()
		molecule_ids = bulk_reader.readAttribute('objectNames')
		counts = bulk_reader.readColumn('counts')
		calculate_request = eval_reader.readColumn('calculate_request_total').flatten()
		evolve_state = eval_reader.readColumn('evolve_state_total').flatten()

		fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

		for index, molecule_id in enumerate(molecule_ids):
			axes[0].plot(time, counts[:, index], marker='o', label=molecule_id)
		axes[0].set_ylabel('Count')
		axes[0].set_title('Platelet scaffold bulk molecules')
		axes[0].legend(loc='best')
		axes[0].grid(True, alpha=0.3)

		axes[1].plot(time, 1000 * calculate_request, marker='o',
			label='calculateRequest')
		axes[1].plot(time, 1000 * evolve_state, marker='o',
			label='evolveState')
		axes[1].set_xlabel('Time (s)')
		axes[1].set_ylabel('Time (ms)')
		axes[1].set_title('Platelet scaffold process timings')
		axes[1].legend(loc='best')
		axes[1].grid(True, alpha=0.3)

		if np.allclose(counts, counts[0]):
			fig.suptitle(
				'Platelet scaffold summary: placeholder counts stay constant by design',
				fontsize=12)

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
