"""
Granule-secretion kinetics analysis plot for the platelet whole-cell model.

Reads the SecretionTrace listener (granule cargo release) and the CalciumTrace
listener (cytosolic Ca²⁺ trigger) from a finished run and produces a 3-panel
figure — the single-platelet activation readout comparable to lumi-aggregometry
(ATP / serotonin release) and flow cytometry (surface P-selectin):

  Panel 1 — Cytosolic Ca²⁺ (nM) vs time. The trigger; shown so the release
             curves below can be read against the Ca²⁺ transient.

  Panel 2 — Cargo released / surface-exposed (%) vs time: dense-granule ADP and
             serotonin vs α-granule fibrinogen and surface P-selectin. The
             PKC × Ca²⁺ secretion gate (0–1) is overlaid on the right axis.
             This is the headline observable — sigmoidal: ~zero at rest, a rapid
             ramp once the gate opens, then a plateau as the pools deplete.

  Panel 3 — Cumulative secreted signalling cargo in the pericellular space:
             ADP (µM, left axis) and serotonin (count, right axis) — the
             autocrine feedback that drives the second wave.

Granule representation note: granules are bulk cargo pools in the current model
(``ADP[dg]`` → ``ADP[e]`` etc.), so release is shown as the fraction of each
conserved cargo pool that has reached the extracellular space. Per-instance
``loaded → fusing → released`` granule state awaits the UniqueMolecule refactor
(issue #5).

Usage
-----
python runscripts/manual/analysisPlatelet.py --plot granule_secretion [sim_dir]
"""

import os

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from models.platelet.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


# Release is considered to have "onset" when the fastest cargo (dense-granule
# ADP) first crosses this fraction released. Threshold-crossing, not argmax —
# the release curves are flat-topped once the pool depletes.
_ONSET_FRAC = 0.01


def _onset_time(t, frac):
	"""First time the released fraction crosses _ONSET_FRAC, or None."""
	crossed = np.where(frac >= _ONSET_FRAC)[0]
	return float(t[crossed[0]]) if len(crossed) else None


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""3-panel granule-secretion kinetics plot — Ca²⁺ trigger, release, autocrine cargo."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile

		# ── SecretionTrace listener data ─────────────────────────────────
		st = TableReader(os.path.join(simOutDir, 'SecretionTrace'))
		t          = st.readColumn('time').flatten()
		gate       = st.readColumn('secretion_gate').flatten()
		adp_frac   = 100 * st.readColumn('adp_released_frac').flatten()
		sht_frac   = 100 * st.readColumn('serotonin_released_frac').flatten()
		fga_frac   = 100 * st.readColumn('fibrinogen_released_frac').flatten()
		psel_frac  = 100 * st.readColumn('pselectin_surface_frac').flatten()
		adp_e_uM   = st.readColumn('adp_e_uM').flatten()
		sht_e      = st.readColumn('serotonin_e').flatten()

		# ── CalciumTrace: cytosolic Ca²⁺ trigger (align lengths) ─────────
		ct = TableReader(os.path.join(simOutDir, 'CalciumTrace'))
		ca_cyt_nM = ct.readColumn('ca_cyt_nM').flatten()
		n = min(len(t), len(ca_cyt_nM))
		t, ca_cyt_nM = t[:n], ca_cyt_nM[:n]
		gate, adp_frac, sht_frac, fga_frac, psel_frac, adp_e_uM, sht_e = (
			a[:n] for a in
			(gate, adp_frac, sht_frac, fga_frac, psel_frac, adp_e_uM, sht_e))

		# ── Figure layout: 3 stacked panels ──────────────────────────────
		fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
		fig.subplots_adjust(hspace=0.35)

		# ── Panel 1: Cytosolic Ca²⁺ trigger ──────────────────────────────
		ax1 = axes[0]
		ax1.plot(t, ca_cyt_nM, color='tab:blue', linewidth=2.0,
			label=r'cytosolic $\mathrm{Ca^{2+}}$')
		ax1.axhline(100.0, color='grey', linewidth=0.8, linestyle=':',
			label='100 nM resting')
		ax1.set_ylabel(r'$[\mathrm{Ca^{2+}}]_\mathrm{cyt}$ (nM)', fontsize=10)
		ax1.set_title('Panel 1 — Cytosolic $\\mathrm{Ca^{2+}}$ trigger',
			fontsize=10, fontweight='bold')
		ax1.legend(loc='upper right', fontsize=8)
		ax1.grid(True, alpha=0.3)
		ax1.set_ylim(bottom=0)

		# ── Panel 2: cargo released / surface-exposed + secretion gate ───
		ax2 = axes[1]
		ax2.plot(t, adp_frac, lw=2.4, color='#c0392b',
			label='ADP (dense granule)')
		ax2.plot(t, sht_frac, lw=1.4, ls='--', color='#e67e22',
			label='serotonin / 5-HT (dense granule)')
		ax2.plot(t, fga_frac, lw=2.4, color='#2980b9',
			label='fibrinogen (alpha granule)')
		ax2.plot(t, psel_frac, lw=1.4, ls='--', color='#27ae60',
			label='P-selectin surface (alpha granule)')
		ax2.set_ylabel('cargo released / surface-exposed (%)', fontsize=10)
		ax2.set_title('Panel 2 — Granule cargo release  (dense leads, alpha follows)',
			fontsize=10, fontweight='bold')
		ax2.grid(True, alpha=0.3)
		ax2.set_ylim(bottom=0)

		ax2g = ax2.twinx()
		ax2g.plot(t, gate, lw=1.6, color='#8e44ad', alpha=0.7,
			label=r'secretion gate (PKC$^*$ $\times$ $\mathrm{Ca^{2+}}$)')
		ax2g.fill_between(t, gate, 0, color='#8e44ad', alpha=0.08)
		ax2g.set_ylabel('secretion gate (0–1)', color='#8e44ad', fontsize=10)
		ax2g.tick_params(axis='y', labelcolor='#8e44ad')
		ax2g.set_ylim(bottom=0)

		l2a, lab2a = ax2.get_legend_handles_labels()
		l2g, lab2g = ax2g.get_legend_handles_labels()
		ax2.legend(l2a + l2g, lab2a + lab2g, loc='center right', fontsize=8)

		# ── Panel 3: cumulative secreted signalling cargo ────────────────
		ax3 = axes[2]
		color_adp = '#c0392b'
		color_5ht = '#e67e22'
		ax3.plot(t, adp_e_uM, color=color_adp, linewidth=2.0,
			label='secreted ADP (pericellular)')
		ax3.set_ylabel(r'$[\mathrm{ADP}]_\mathrm{e}$ (µM)',
			color=color_adp, fontsize=10)
		ax3.tick_params(axis='y', labelcolor=color_adp)
		ax3.set_ylim(bottom=0)

		ax3r = ax3.twinx()
		ax3r.plot(t, sht_e, color=color_5ht, linewidth=1.4, linestyle='--',
			label='secreted serotonin')
		ax3r.set_ylabel('secreted 5-HT (count)', color=color_5ht, fontsize=10)
		ax3r.tick_params(axis='y', labelcolor=color_5ht)
		ax3r.set_ylim(bottom=0)

		l3a, lab3a = ax3.get_legend_handles_labels()
		l3r, lab3r = ax3r.get_legend_handles_labels()
		ax3.legend(l3a + l3r, lab3a + lab3r, loc='lower right', fontsize=8)
		ax3.set_xlabel('time (s)', fontsize=10)
		ax3.set_title('Panel 3 — Cumulative secreted autocrine cargo  '
			'(drives the second wave)', fontsize=10, fontweight='bold')
		ax3.grid(True, alpha=0.3)

		# ── Sigmoidal-shape / acceptance annotation ──────────────────────
		onset = _onset_time(t, adp_frac / 100.0)
		plateau = float(adp_frac[-1])
		is_sigmoidal = (
			float(adp_frac[0]) < 1.0          # ~zero at rest
			and onset is not None             # rises
			and plateau > float(adp_frac[0])) # to a higher plateau
		lines = [
			'Secretion kinetics (dense-granule ADP):',
			f"  at rest (t=0):  {adp_frac[0]:.1f} %",
			f"  release onset:  "
				+ (f"t = {onset:.0f} s" if onset is not None else 'none'),
			f"  plateau (t={t[-1]:.0f} s):  {plateau:.0f} %",
			f"  sigmoidal (rest→ramp→plateau):  "
				f"{'✓' if is_sigmoidal else '✗'}",
		]
		fig.text(0.01, 0.005, '\n'.join(lines), fontsize=7.5,
			verticalalignment='bottom',
			bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

		fig.suptitle('Granule secretion kinetics', fontsize=13, fontweight='bold')
		fig.text(0.5, 0.945,
			'PKC + $\\mathrm{Ca^{2+}}$-gated SNARE secretion relocates pre-existing '
			'granule cargo to the extracellular space (P-selectin to the surface). '
			'Release is ~zero at rest and dense-granule cargo leads alpha-granule cargo.',
			ha='center', va='top', fontsize=8, wrap=True)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
