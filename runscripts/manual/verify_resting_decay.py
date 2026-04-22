"""
Verify and plot RestingDecay process: run a 1-day platelet sim, report
molecule counts, and save a chart of simulated protein decay vs theory.

Usage:
    PYTHONPATH="$PWD" python runscripts/manual/verify_resting_decay.py [out_dir]

If out_dir is omitted, defaults to out/decay-verify.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.io.tablereader import TableReader
from reconstruction.platelet.dataclasses.internal_state import (
	PROTEIN_MOLECULE_IDS,
	_MOLECULES,
)

_HALF_LIFE = 7 * 24 * 3600  # 7 days in seconds

_CLASS_BY_ID = {m[0]: m[3] for m in _MOLECULES}


def theoretical_fraction(t):
	"""Fraction remaining under exponential decay: exp(-ln2 * t / t_half)."""
	return np.exp(-np.log(2) * t / _HALF_LIFE)


def main(out_dir='out/decay-verify'):
	print("Running 1-day platelet simulation...")
	result = run_platelet_sim(out_dir, length_sec=86400, seed=0)
	sim_out_dir = result['sim_out_dir']

	reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
	counts = reader.readColumn('counts')   # shape: (timesteps, n_molecules)
	ids = list(reader.readAttribute('objectNames'))
	n_steps = counts.shape[0]
	t = np.arange(n_steps)

	# Console report
	p_per_step = 1 - np.exp(-np.log(2) / _HALF_LIFE)
	expected_pct = 100.0 * (1 - (1 - p_per_step) ** 86400)

	print(f"\n=== RestingDecay results (1 simulated day) ===")
	print(f"  {'Molecule':<25}  {'Class':<12}  {'Start':>12}  {'End':>12}  {'Lost':>8}  {'Decay %':>8}")
	print("  " + "-" * 90)
	for i, mol_id in enumerate(ids):
		c0, cf = int(counts[0, i]), int(counts[-1, i])
		lost = c0 - cf
		pct = 100.0 * lost / max(c0, 1)
		mol_class = _CLASS_BY_ID.get(mol_id, '?')
		print(f"  {mol_id:<25}  {mol_class:<12}  {c0:>12,}  {cf:>12,}  {lost:>8,}  {pct:>7.2f}%")

	print(f"\n  Theoretical protein decay over 1 day: {expected_pct:.2f}%")
	print(f"  (half-life = 7 days; p per 1s step = {p_per_step:.2e})")

	# Plot: all proteins normalised to fraction remaining
	fig, ax = plt.subplots(figsize=(11, 6))
	fig.suptitle('RestingDecay — protein decay vs theory (1 day, normalised)', fontsize=13)

	protein_indices = [i for i, mol_id in enumerate(ids) if mol_id in PROTEIN_MOLECULE_IDS]
	cmap = plt.get_cmap('tab10')

	for colour_i, i in enumerate(protein_indices):
		n0 = counts[0, i]
		if n0 == 0:
			continue
		frac = counts[:, i] / n0
		ax.step(t / 3600, frac, where='post', alpha=0.75,
			color=cmap(colour_i),
			label=f'Sim: {ids[i]}')

	# single theory curve (same for all proteins)
	t_theory = np.linspace(0, n_steps - 1, 500)
	ax.plot(t_theory / 3600, theoretical_fraction(t_theory),
		color='black', linestyle='--', linewidth=2,
		label=r'Theory: $e^{-\ln 2 \cdot t\;/\;t_{1/2}}$')

	ax.set_xlabel('Simulated time (hours)')
	ax.set_ylabel('Fraction remaining  ($n\\ /\\ n_0$)')
	ax.set_xlim(0, (n_steps - 1) / 3600)
	ax.set_ylim(bottom=0)

	annotation = (
		'Theory (dashed):\n'
		r'  $N(t)/N_0 = e^{-\ln 2 \cdot t \;/\; t_{1/2}}$'
		'\n\n'
		'Sim (solid lines):\n'
		r'  $p = 1 - e^{-\ln 2 \cdot \Delta t \;/\; t_{1/2}}$'  '\n'
		r'  $\Delta n \sim \mathrm{Binomial}(n,\; p)$'  '\n'
		f'  $t_{{1/2}}$ = 7 days,  $\\Delta t$ = 1 s\n'
		f'  $p$ per step = {p_per_step:.2e}'
	)
	ax.text(
		0.97, 0.97, annotation,
		transform=ax.transAxes,
		fontsize=9,
		verticalalignment='top',
		horizontalalignment='right',
		bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='grey', alpha=0.9),
		)

	ax.legend(loc='lower left', fontsize=8)
	plt.tight_layout()

	plot_out = os.path.join(os.path.dirname(sim_out_dir), 'plotOut', 'resting_decay.png')
	os.makedirs(os.path.dirname(plot_out), exist_ok=True)
	fig.savefig(plot_out, dpi=150)
	print(f"\nPlot saved to: {plot_out}")


if __name__ == '__main__':
	main(sys.argv[1] if len(sys.argv) > 1 else 'out/decay-verify')
