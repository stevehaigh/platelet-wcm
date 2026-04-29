"""
Ca²⁺ transient analysis plot for the platelet whole-cell model.

Reads the CalciumTrace listener output and produces a 3-panel figure:

  Panel 1 — Cytosolic Ca²⁺ (nM) vs time, overlaid with a schematic
             reference curve approximating the Dolan & Diamond (2014)
             Fig. 4 transient shape. The peak and plateau values are
             taken from the Dolan 2014 text (peak ~400 nM, SOCE-sustained
             plateau ~200 nM, baseline ~100 nM).

  Panel 2 — DTS store Ca²⁺ (µM) vs time, showing store depletion and
             partial SOCE-driven refill.

  Panel 3 — IP₃ (nM, left axis) and SOCE flux (nM/s, right axis) vs time.

Validation target: Dolan & Diamond (2014) Biophys J 106:2049-60, Fig. 4.

Usage
-----
python runscripts/manual/analysisPlatelet.py --plot calcium_trace [sim_dir]
"""

import os

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from models.platelet.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


# ── Dolan 2014 Fig. 4 schematic reference (approximate analytical fit) ────
# The published trace shows:
#   baseline ~100 nM, peak ~400 nM at t≈5 s, SOCE plateau ~200 nM,
#   return toward baseline over ~120 s.
# These values are taken from the Dolan 2014 text and used as a schematic
# guide only; the exact digitised trace is not available.
_DOLAN_CA_BASELINE_NM   = 100.0
_DOLAN_CA_PEAK_NM       = 400.0
_DOLAN_CA_PLATEAU_NM    = 200.0
_DOLAN_T_PEAK_S         = 5.0
_DOLAN_TAU_RISE_S       = 2.0
_DOLAN_TAU_DECAY_S      = 30.0
_DOLAN_TAU_PLATEAU_S    = 120.0


def _dolan_reference_nM(t):
	"""Schematic Ca²⁺ transient matching Dolan & Diamond (2014) Fig. 4 shape.

	Returns cytosolic Ca²⁺ in nM for scalar or array `t` (seconds from
	stimulus onset).  This is an analytical approximation; it is not
	digitised data from the original figure.
	"""
	t = np.asarray(t, dtype=float)
	rise   = 1.0 - np.exp(-np.maximum(t, 0.0) / _DOLAN_TAU_RISE_S)
	decay  = np.exp(-np.maximum(t - _DOLAN_T_PEAK_S, 0.0) / _DOLAN_TAU_DECAY_S)
	plateau_rise = 1.0 - np.exp(-np.maximum(t, 0.0) / _DOLAN_TAU_PLATEAU_S)
	transient = (_DOLAN_CA_PEAK_NM - _DOLAN_CA_BASELINE_NM) * rise * decay
	plateau   = (_DOLAN_CA_PLATEAU_NM - _DOLAN_CA_BASELINE_NM) * plateau_rise * decay
	return _DOLAN_CA_BASELINE_NM + transient + plateau


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""Ca²⁺ transient plot with Dolan 2014 reference overlay."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile

		reader = TableReader(os.path.join(simOutDir, 'CalciumTrace'))
		time          = reader.readColumn('time').flatten()
		ca_cyt_nM     = reader.readColumn('ca_cyt_nM').flatten()
		ca_dts_uM     = reader.readColumn('ca_dts_uM').flatten()
		ip3_nM        = reader.readColumn('ip3_nM').flatten()
		soce_flux_nMs = reader.readColumn('soce_flux_nMs').flatten()

		# Reference curve evaluated over the same time window.
		t_ref = np.linspace(0.0, time[-1] if len(time) else 200.0, 400)
		ca_ref_nM = _dolan_reference_nM(t_ref)

		fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

		# ── Panel 1: cytosolic Ca²⁺ ──────────────────────────────────────
		ax1 = axes[0]
		ax1.plot(time, ca_cyt_nM, color='tab:blue', linewidth=1.5,
			label='simulated')
		ax1.plot(t_ref, ca_ref_nM, color='tab:orange', linewidth=1.2,
			linestyle='--', label='Dolan 2014 (schematic)')
		ax1.axhline(_DOLAN_CA_BASELINE_NM, color='grey', linewidth=0.7,
			linestyle=':', label='100 nM baseline')
		ax1.axvline(0.0, color='k', linewidth=0.8, linestyle='--',
			alpha=0.6, label='IP3 stimulus (t = 0)')
		ax1.set_ylabel('Cytosolic Ca2+ (nM)')
		ax1.set_title('Platelet Ca2+ transient — Dolan & Diamond 2014 validation')
		ax1.legend(loc='upper right', fontsize=8)
		ax1.grid(True, alpha=0.3)
		ax1.set_ylim(bottom=0)

		# ── Panel 2: DTS store Ca²⁺ ──────────────────────────────────────
		ax2 = axes[1]
		ax2.plot(time, ca_dts_uM, color='tab:red', linewidth=1.5)
		ax2.set_ylabel('DTS Ca2+ (uM)')
		ax2.set_title('DTS store Ca2+')
		ax2.grid(True, alpha=0.3)
		ax2.set_ylim(bottom=0)

		# ── Panel 3: IP₃ and SOCE flux ────────────────────────────────────
		ax3 = axes[2]
		color_ip3  = 'tab:green'
		color_soce = 'tab:purple'
		ax3.plot(time, ip3_nM, color=color_ip3, linewidth=1.5, label='IP3 (nM)')
		ax3.set_ylabel('IP3 (nM)', color=color_ip3)
		ax3.tick_params(axis='y', labelcolor=color_ip3)

		ax3r = ax3.twinx()
		ax3r.plot(time, soce_flux_nMs, color=color_soce, linewidth=1.2,
			linestyle='--', label='SOCE flux (nM/s)')
		ax3r.set_ylabel('SOCE influx (nM/s)', color=color_soce)
		ax3r.tick_params(axis='y', labelcolor=color_soce)

		# Combined legend for both axes.
		lines1, labels1 = ax3.get_legend_handles_labels()
		lines2, labels2 = ax3r.get_legend_handles_labels()
		ax3.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=8)
		ax3.set_xlabel('Time (s)')
		ax3.set_title('IP3 stimulus and SOCE influx')
		ax3.grid(True, alpha=0.3)

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
