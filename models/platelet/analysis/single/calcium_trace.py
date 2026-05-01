"""
Ca²⁺ transient analysis plot for the platelet whole-cell model.

Reads the CalciumTrace listener and BulkMolecules outputs and produces a
5-panel figure:

  Panel 1 — Cytosolic Ca²⁺ (nM) vs time, overlaid with a schematic
             Dolan & Diamond (2014) Fig. 4 reference curve.

  Panel 2 — Calmodulin sub-state dynamics (counts): CaM_free, Ca₂·CaM,
             Ca₄·CaM, and PMCA·CaM, shown as a stacked area plot.
             Illustrates CaM acting as a cytosolic Ca²⁺ buffer.

  Panel 3 — PMCA activation state (counts): free PMCA, basal-occupied
             PMCA·Ca, CaM-activated forms (Ca₄·CaM·PMCA and
             Ca₄·CaM·PMCA·Ca), and deactivating PMCA·CaM.

  Panel 4 — DTS store Ca²⁺ (µM) with STIM1 dimer count on the right axis,
             showing store depletion driving SOCE activation.

  Panel 5 — IP₃ (nM, left) and SOCE flux (nM/s, right) vs time.

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


# ── Dolan 2014 Fig. 4 schematic reference ────────────────────────────────
# Baseline ~100 nM, peak ~400 nM at t≈5 s, SOCE plateau ~200 nM,
# return over ~120 s. Analytical approximation, not digitised data.
_DOLAN_CA_BASELINE_NM  = 100.0
_DOLAN_CA_PEAK_NM      = 400.0
_DOLAN_CA_PLATEAU_NM   = 200.0
_DOLAN_T_PEAK_S        = 5.0
_DOLAN_TAU_RISE_S      = 2.0
_DOLAN_TAU_DECAY_S     = 30.0
_DOLAN_TAU_PLATEAU_S   = 120.0


def _dolan_reference_nM(t):
	"""Schematic Ca²⁺ transient matching Dolan & Diamond (2014) Fig. 4 shape."""
	t = np.asarray(t, dtype=float)
	rise          = 1.0 - np.exp(-np.maximum(t, 0.0) / _DOLAN_TAU_RISE_S)
	decay         = np.exp(-np.maximum(t - _DOLAN_T_PEAK_S, 0.0) / _DOLAN_TAU_DECAY_S)
	plateau_rise  = 1.0 - np.exp(-np.maximum(t, 0.0) / _DOLAN_TAU_PLATEAU_S)
	transient = (_DOLAN_CA_PEAK_NM - _DOLAN_CA_BASELINE_NM) * rise * decay
	plateau   = (_DOLAN_CA_PLATEAU_NM - _DOLAN_CA_BASELINE_NM) * plateau_rise * decay
	return _DOLAN_CA_BASELINE_NM + transient + plateau


def _dolan_reference_band_nM(t):
	"""±30% band around the Dolan reference to visualise tolerance range."""
	ref = _dolan_reference_nM(t)
	return ref * 0.7, ref * 1.3


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""5-panel Ca²⁺ transient plot — cytosol, CaM, PMCA, DTS/STIM, IP3/SOCE."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile

		# ── CalciumTrace listener data ────────────────────────────────────
		ct = TableReader(os.path.join(simOutDir, 'CalciumTrace'))
		time          = ct.readColumn('time').flatten()
		ca_cyt_nM     = ct.readColumn('ca_cyt_nM').flatten()
		ca_dts_uM     = ct.readColumn('ca_dts_uM').flatten()
		ip3_nM        = ct.readColumn('ip3_nM').flatten()
		soce_flux_nMs = ct.readColumn('soce_flux_nMs').flatten()
		stim1_dim     = ct.readColumn('stim1_dim').flatten()

		# ── BulkMolecules: CaM and PMCA sub-states ───────────────────────
		rb = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
		mol_ids    = list(rb.readAttribute('objectNames'))
		all_counts = rb.readColumn('counts')
		# BulkMolecules may have an extra initial row; align to CalciumTrace length.
		n = min(len(time), all_counts.shape[0])
		all_counts = all_counts[:n]

		def get_counts(name):
			try:
				return all_counts[:, mol_ids.index(name)].astype(float)
			except ValueError:
				return np.zeros(n)

		cam_free         = get_counts('CaM_free[c]')
		ca2_cam          = get_counts('Ca2_CaM[c]')
		ca4_cam          = get_counts('Ca4_CaM[c]')
		pmca_cam         = get_counts('PMCA_CaM[pl]')
		pmca_free        = get_counts('PMCA[pl]')
		pmca_ca          = get_counts('PMCA_Ca[pl]')
		ca4_cam_pmca     = get_counts('Ca4_CaM_PMCA[pl]')
		ca4_cam_pmca_ca  = get_counts('Ca4_CaM_PMCA_Ca[pl]')

		t = time[:n]

		# ── Reference curve ───────────────────────────────────────────────
		t_ref = np.linspace(0.0, t[-1] if len(t) else 200.0, 500)
		ca_ref_nM = _dolan_reference_nM(t_ref)
		ca_lo, ca_hi = _dolan_reference_band_nM(t_ref)

		# ── Figure layout: 5 stacked panels ──────────────────────────────
		fig, axes = plt.subplots(5, 1, figsize=(10, 18), sharex=True)
		fig.subplots_adjust(hspace=0.45)

		# ── Panel 1: Cytosolic Ca²⁺ ──────────────────────────────────────
		ax1 = axes[0]
		ax1.fill_between(t_ref, ca_lo, ca_hi, alpha=0.15, color='tab:orange',
			label='Dolan 2014 ±30%')
		ax1.plot(t_ref, ca_ref_nM, color='tab:orange', linewidth=1.2,
			linestyle='--', label='Dolan 2014 (schematic)')
		ax1.plot(t, ca_cyt_nM, color='tab:blue', linewidth=2.0,
			label='simulated (Phase 1)')
		ax1.axhline(100.0, color='grey', linewidth=0.8, linestyle=':',
			label='100 nM resting target')
		ax1.set_ylabel('[Ca²⁺]_cyt (nM)', fontsize=10)
		ax1.set_title('Panel 1 — Cytosolic Ca²⁺ transient  (Phase 1: CaM + 5-state PMCA)',
			fontsize=10, fontweight='bold')
		ax1.legend(loc='upper right', fontsize=8)
		ax1.grid(True, alpha=0.3)
		ax1.set_ylim(bottom=0)
		# Annotate early transient peak (within first 20 s — biological peak,
		# not the slow drift maximum that can occur later).
		early = np.where(t <= 20.0)[0]
		early_peak_idx = int(early[np.argmax(ca_cyt_nM[early])]) if len(early) else 0
		ax1.annotate(
			f'transient peak = {ca_cyt_nM[early_peak_idx]:.0f} nM @ t={t[early_peak_idx]:.0f} s',
			xy=(t[early_peak_idx], ca_cyt_nM[early_peak_idx]),
			xytext=(t[early_peak_idx] + 12, ca_cyt_nM[early_peak_idx] * 1.08),
			fontsize=8, color='tab:blue',
			arrowprops=dict(arrowstyle='->', color='tab:blue', lw=0.8),
		)
		# Annotate steady-state value at t=200 s.
		ax1.annotate(
			f'steady state ≈ {ca_cyt_nM[-1]:.0f} nM',
			xy=(t[-1], ca_cyt_nM[-1]),
			xytext=(t[-1] - 60, ca_cyt_nM[-1] + 40),
			fontsize=8, color='grey',
			arrowprops=dict(arrowstyle='->', color='grey', lw=0.7),
		)

		# ── Panel 2: CaM sub-state dynamics (stacked area) ───────────────
		ax2 = axes[1]
		# Stacked from bottom: free → Ca2 → Ca4 → bound to PMCA
		y_free = cam_free
		y_ca2  = ca2_cam
		y_ca4  = ca4_cam
		y_pmca_cam_col = pmca_cam   # CaM stuck in PMCA·CaM deactivating state
		total_cam = y_free + y_ca2 + y_ca4 + y_pmca_cam_col

		ax2.stackplot(
			t,
			y_free, y_ca2, y_ca4, y_pmca_cam_col,
			labels=['CaM free', 'Ca₂·CaM', 'Ca₄·CaM', 'PMCA·CaM (deact.)'],
			colors=['#d4e6f1', '#5dade2', '#1a5276', '#e59866'],
			alpha=0.85,
		)
		ax2.set_ylabel('CaM sub-state (count)', fontsize=10)
		ax2.set_title('Panel 2 — Calmodulin Ca²⁺-binding ladder  (buffering effect)',
			fontsize=10, fontweight='bold')
		ax2.legend(loc='center right', fontsize=8)
		ax2.grid(True, alpha=0.3)
		# Annotate total CaM conserved.
		ax2.axhline(total_cam[0], color='k', linewidth=0.7, linestyle=':',
			label='total CaM (conserved)')

		# ── Panel 3: PMCA activation state ───────────────────────────────
		ax3 = axes[2]
		pmca_cam_act = ca4_cam_pmca + ca4_cam_pmca_ca   # CaM-activated pool
		ax3.stackplot(
			t,
			pmca_free, pmca_ca, pmca_cam_act, pmca_cam,
			labels=[
				'PMCA free',
				'PMCA·Ca (basal active)',
				'Ca₄·CaM·PMCA(·Ca) (CaM-activated)',
				'PMCA·CaM (deactivating)',
			],
			colors=['#aed6f1', '#2e86c1', '#1b4f72', '#e59866'],
			alpha=0.85,
		)
		ax3.set_ylabel('PMCA sub-state (count)', fontsize=10)
		ax3.set_title('Panel 3 — PMCA activation state  (basal vs CaM-activated)',
			fontsize=10, fontweight='bold')
		ax3.legend(loc='center right', fontsize=8)
		ax3.grid(True, alpha=0.3)
		# Note: CaM-activated pool (dark blue) is near-zero at any instant because
		# k10 = 30 s⁻¹ cycles it faster than the 1-s outer timestep can capture.
		# Its contribution to extrusion is real but invisible in this integer snapshot.
		ax3.text(0.02, 0.97,
			'Ca₄·CaM·PMCA(·Ca) cycles at k₁₀=30 s⁻¹\n'
			'→ near-zero instantaneous count,\n'
			'but ~5× higher turnover rate than basal',
			transform=ax3.transAxes, fontsize=7, verticalalignment='top',
			bbox=dict(boxstyle='round', facecolor='#d6eaf8', alpha=0.7))

		# ── Panel 4: DTS Ca²⁺ + STIM1 dimer count ───────────────────────
		ax4 = axes[3]
		color_dts  = 'tab:red'
		color_stim = 'tab:purple'
		ax4.plot(t, ca_dts_uM, color=color_dts, linewidth=1.8, label='DTS Ca²⁺')
		ax4.set_ylabel('[Ca²⁺]_DTS (µM)', color=color_dts, fontsize=10)
		ax4.tick_params(axis='y', labelcolor=color_dts)
		ax4.set_ylim(bottom=0)

		ax4r = ax4.twinx()
		ax4r.plot(t, stim1_dim, color=color_stim, linewidth=1.4,
			linestyle='--', label='STIM1 dimers')
		ax4r.set_ylabel('STIM1 dimer count', color=color_stim, fontsize=10)
		ax4r.tick_params(axis='y', labelcolor=color_stim)

		lines4a, labs4a = ax4.get_legend_handles_labels()
		lines4b, labs4b = ax4r.get_legend_handles_labels()
		ax4.legend(lines4a + lines4b, labs4a + labs4b, loc='center right', fontsize=8)
		ax4.set_title('Panel 4 — DTS store depletion + STIM1 mobilisation',
			fontsize=10, fontweight='bold')
		ax4.grid(True, alpha=0.3)

		# ── Panel 5: IP₃ and SOCE flux ────────────────────────────────────
		ax5 = axes[4]
		color_ip3  = 'tab:green'
		color_soce = 'darkorange'
		ax5.plot(t, ip3_nM, color=color_ip3, linewidth=1.8, label='IP₃ (nM)')
		ax5.set_ylabel('IP₃ (nM)', color=color_ip3, fontsize=10)
		ax5.tick_params(axis='y', labelcolor=color_ip3)

		ax5r = ax5.twinx()
		ax5r.plot(t, soce_flux_nMs, color=color_soce, linewidth=1.4,
			linestyle='--', label='SOCE flux (nM/s)')
		ax5r.set_ylabel('SOCE influx (nM/s)', color=color_soce, fontsize=10)
		ax5r.tick_params(axis='y', labelcolor=color_soce)

		lines5a, labs5a = ax5.get_legend_handles_labels()
		lines5b, labs5b = ax5r.get_legend_handles_labels()
		ax5.legend(lines5a + lines5b, labs5a + labs5b, loc='upper right', fontsize=8)
		ax5.set_xlabel('Time (s)', fontsize=10)
		ax5.set_title('Panel 5 — IP₃ stimulus and SOCE influx',
			fontsize=10, fontweight='bold')
		ax5.grid(True, alpha=0.3)

		# ── Acceptance-criteria annotation box ────────────────────────────
		# Peak is searched within the first 20 s (biological transient peak,
		# not the slow SOCE plateau or any long-term drift).
		early_mask = t[:n] <= 20.0
		if np.any(early_mask):
			peak_nM = float(np.max(ca_cyt_nM[early_mask]))
		else:
			peak_nM = float(np.max(ca_cyt_nM))
		at50_nM   = float(ca_cyt_nM[min(50, len(ca_cyt_nM) - 1)])
		at200_nM  = float(ca_cyt_nM[min(200, len(ca_cyt_nM) - 1)])
		dts_min   = float(np.min(ca_dts_uM))
		# DTS minimum > 0: Phase 1 uses Dolan Table S1 ICs which are not at
		# mathematical steady state — IP3R flux drains DTS in <1 s.  This
		# criterion requires Phase 2 (true steady-state IC derivation).
		dts_note = '(Phase 2 reqd)' if dts_min <= 0 else ''
		lines = [
			'Phase 1 acceptance criteria:',
			f"  peak [Ca²⁺]_cyt (t≤20s): {peak_nM:.0f} nM"
				f"  {'✓ PASS' if 200 <= peak_nM <= 800 else '✗ FAIL'}"
				f"  (target 200–800 nM)",
			f"  at t=50 s:  {at50_nM:.0f} nM"
				f"  {'✓ PASS' if at50_nM < 1000 else '✗ FAIL'}"
				f"  (target <1000 nM)",
			f"  at t=200 s: {at200_nM:.0f} nM"
				f"  {'✓ PASS' if at200_nM < 1000 else '✗ FAIL'}"
				f"  (target <1000 nM)",
			f"  DTS min: {dts_min:.4f} µM"
				f"  {'✓ PASS' if dts_min > 0 else '✗ FAIL'}"
				f"  (target >0 µM) {dts_note}",
		]
		fig.text(0.01, 0.005, '\n'.join(lines), fontsize=7.5,
			verticalalignment='bottom',
			bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
