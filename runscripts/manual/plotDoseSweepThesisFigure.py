#!/usr/bin/env python3
"""Render the two-panel dose-sweep thesis figure from committed sweep data.

Panel A: 3-D surface of peak cytosolic Ca²⁺ over the ADP × thrombin grid
         — the binary store-release plateau with a sub-threshold cliff.
Panel B: heatmap of Ca²⁺ AUC above resting baseline — the graded
         integrated response.

Both panels use viridis (one colourblind-safe colormap) so the figure
reads as a single unit. Input is the committed transition-sweep NPZ; no
simulation is rerun. Run from the repo root:

    PYTHONPATH=$PWD python runscripts/manual/plotDoseSweepThesisFigure.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 — registers 3D projection

NPZ = 'reports/figures/v0.6/dose-sweep-9x9-transition.npz'
OUT = 'reports/figures/v0.6/dose-sweep-transition-surface-auc-v0.6.png'


def _log_edges(centers: np.ndarray) -> np.ndarray:
	"""Cell edges in log space so pcolormesh quads centre on grid points."""
	log_c = np.log10(centers)
	mid = 0.5 * (log_c[:-1] + log_c[1:])
	edges = np.concatenate([
		[2.0 * log_c[0] - mid[0]], mid, [2.0 * log_c[-1] - mid[-1]],
	])
	return 10.0 ** edges


def _decade_exponents(centers: np.ndarray) -> list[int]:
	"""Integer powers of ten spanning the data range — clean log ticks
	decoupled from the (non-decade) grid points."""
	lo = int(np.floor(np.log10(centers.min())))
	hi = int(np.ceil(np.log10(centers.max())))
	return list(range(lo, hi + 1))


def main() -> None:
	d = np.load(NPZ)  # numeric arrays only — no allow_pickle needed
	adp = d['adp_grid']
	thr = d['thr_grid']
	peak_ca = d['peak_ca_nM']
	auc = d['auc_ca_nMs']

	fig = plt.figure(figsize=(13, 5.5))

	# Panel A — 3-D surface of peak cytosolic Ca²⁺ (log10 stimulus axes).
	axA = fig.add_subplot(1, 2, 1, projection='3d')
	A, T = np.meshgrid(np.log10(adp), np.log10(thr))
	surf = axA.plot_surface(A, T, peak_ca, cmap='viridis',
		edgecolor='gray', linewidth=0.3, antialiased=True)
	adp_exp = _decade_exponents(adp)
	thr_exp = _decade_exponents(thr)
	axA.set_xticks(adp_exp)
	axA.set_xticklabels([rf'$10^{{{e}}}$' for e in adp_exp], fontsize=8)
	axA.set_yticks(thr_exp)
	axA.set_yticklabels([rf'$10^{{{e}}}$' for e in thr_exp], fontsize=8)
	axA.set_xlabel('ADP peak (µM)', fontsize=9)
	axA.set_ylabel('Thrombin peak (nM)', fontsize=9)
	axA.set_zlabel(r'Peak cytosolic Ca$^{2+}$ (nM)', fontsize=9)
	axA.set_title(r'A  Peak cytosolic Ca$^{2+}$ — binary store release',
		fontsize=11, loc='left')
	fig.colorbar(surf, ax=axA, shrink=0.55, pad=0.1, label='nM')

	# Panel B — AUC heatmap (graded integrated response).
	axB = fig.add_subplot(1, 2, 2)
	pcm = axB.pcolormesh(_log_edges(adp), _log_edges(thr), auc,
		shading='flat', cmap='viridis')
	axB.set_xscale('log')
	axB.set_yscale('log')
	axB.set_xlabel('ADP peak (µM)', fontsize=9)
	axB.set_ylabel('Thrombin peak (nM)', fontsize=9)
	axB.set_title(r'B  Ca$^{2+}$ AUC above rest — graded response',
		fontsize=11, loc='left')
	fig.colorbar(pcm, ax=axB, label='nM·s')

	fig.savefig(OUT, dpi=200, bbox_inches='tight')
	plt.close(fig)

	print(f'wrote {OUT}')
	print(f'peak Ca²⁺  min/max: {peak_ca.min():.0f} / {peak_ca.max():.0f} nM')
	print(f'AUC        min/max: {auc.min():.0f} / {auc.max():.0f} nM·s')
	print(f'ADP range:      {adp.min():g} – {adp.max():g} µM')
	print(f'thrombin range: {thr.min():g} – {thr.max():g} nM')


if __name__ == '__main__':
	main()
