"""
Integrin demo figure (TUI) — PAC-1 activation and its cAMP/PKA brake.

Focused figure for the integrin demo (demo 3, clopidogrel) and the Glanzmann
foil. Two stacked panels sharing a time axis (no twin axis — that conflated the
fraction and brake scales):

  Top    — PAC-1 active fraction (the flow-cytometry readout) + the PKC×Ca²⁺
           activation gate. 0–1.
  Bottom — the PKA brake factor (≥1) that the cAMP arm sets: clopidogrel keeps
           it engaged (≈1), so PAC-1 falls with the integrin fully intact.

Grey baseline overlay via PLATELET_BASELINE_SIMOUT (e.g. control vs clopidogrel).
For a clean overlay, run the baseline and current sims at the SAME length.

Usage: analysisPlatelet.py --plot demo_integrin [sim_dir]
"""

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt

from models.platelet.analysis import singleAnalysisPlot
from models.platelet.analysis.single._demo_common import (
	make_draw, resolve_baseline, suptitle)
from wholecell.analysis.analysis_tools import exportFigure


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""PAC-1 active fraction (top) and the PKA brake (bottom), baseline overlay."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile
		base = resolve_baseline(simOutDir)
		draw = make_draw(simOutDir, base)

		fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 8),
			sharex=True)
		fig.subplots_adjust(hspace=0.22, top=0.9, bottom=0.09, right=0.97)

		# ── Top: PAC-1 active fraction (the readout) ─────────────────────
		draw(ax_top, 'IntegrinTrace', 'active_frac', 'tab:purple',
			r'PAC-1$^{+}$ active fraction')
		# Gate is current-run context only (no baseline) to keep it readable.
		draw(ax_top, 'IntegrinTrace', 'integrin_gate', 'tab:olive',
			r'PKC$\times\mathrm{Ca}^{2+}$ activation gate', lw=1.3, ls=':',
			baseline=False)
		ax_top.set_ylabel('active fraction  (0–1)')
		ax_top.set_ylim(0, 1)
		ax_top.set_title('αIIbβ3 activation — PAC-1 readout',
			fontsize=11, fontweight='bold')
		ax_top.legend(loc='upper left', fontsize=8)
		ax_top.grid(True, alpha=0.3)

		# ── Bottom: PKA brake (the cAMP mechanism) ───────────────────────
		draw(ax_bot, 'IntegrinTrace', 'pka_brake', 'tab:brown',
			'PKA brake factor')
		ax_bot.axhline(1.0, color='grey', lw=0.8, ls=':',
			label='1.0 = brake fully engaged (resting cAMP)')
		ax_bot.set_ylabel('PKA brake factor\n(≥1; higher → more activation)')
		ax_bot.set_ylim(bottom=0.95)
		ax_bot.set_xlabel('Time (s)')
		ax_bot.set_title(
			'cAMP/PKA brake on activation  (P2Y12 → Gi → cAMP → PKA)',
			fontsize=11, fontweight='bold')
		ax_bot.legend(loc='best', fontsize=8)
		ax_bot.grid(True, alpha=0.3)

		suptitle(fig,
			'Integrin αIIbβ3 — PAC-1 activation and its cAMP/PKA brake', base)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
