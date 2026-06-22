"""
Reversibility demo figure (TUI) — reversible αIIbβ3 activation.

Focused figure for the integrin-reversibility demo (TUI demo 5), mapping to Zou
et al. 2022 (Int J Mol Sci 23:12512): αIIbβ3 activation is an intrinsically
reversible inside-out switch whose persistence depends on sustained autocrine
ADP. Under a *transient* agonist the PAC-1 active fraction rises then relaxes,
and the relaxation tracks the autocrine ADP[e] being cleared by ecto-NTPDase.
With a baseline overlay (PLATELET_BASELINE_SIMOUT), a weak (ADP) run reverses
visibly more than a strong (thrombin) one.

Usage: analysisPlatelet.py --plot demo_reversibility [sim_dir]
"""

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt

from models.platelet.analysis import singleAnalysisPlot
from models.platelet.analysis.single._demo_common import (
	make_draw, resolve_baseline, suptitle)
from wholecell.analysis.analysis_tools import exportFigure


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""PAC-1 active fraction (rise-then-fall) and its autocrine-ADP driver."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile
		base = resolve_baseline(simOutDir)
		draw = make_draw(simOutDir, base)

		fig, axes = plt.subplots(1, 2, figsize=(13, 5))
		fig.subplots_adjust(wspace=0.26, top=0.85, bottom=0.14)

		# PAC-1 active fraction — the reversible activation
		ax = axes[0]
		draw(ax, 'IntegrinTrace', 'active_frac', 'tab:purple', 'PAC-1 active')
		ax.set_ylabel('PAC-1 active fraction')
		ax.set_ylim(0, 1)
		ax.set_title('Reversible αIIbβ3 activation  (PAC-1 rises then falls)',
			fontsize=10, fontweight='bold')
		ax.legend(loc='best', fontsize=8)

		# Autocrine ADP[e] — the reversal driver (cleared by ecto-NTPDase)
		ax = axes[1]
		draw(ax, 'SecretionTrace', 'adp_e_uM', 'tab:cyan',
			r'autocrine $[\mathrm{ADP}]_e$')
		ax.set_ylabel(r'secreted $[\mathrm{ADP}]_e$ ($\mu$M)')
		ax.set_ylim(bottom=0)
		ax.set_title('Autocrine $[\\mathrm{ADP}]_e$  (the driver; cleared by ecto-NTPDase)',
			fontsize=10, fontweight='bold')
		ax.legend(loc='best', fontsize=8)

		for a in axes:
			a.set_xlabel('Time (s)')
			a.grid(True, alpha=0.3)
		suptitle(fig, 'Integrin reversibility — PAC-1 relaxes as the ADP drive '
			'is withdrawn', base)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
