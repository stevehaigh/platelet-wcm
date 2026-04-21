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
from reconstruction.platelet.dataclasses.internal_state import PROTEIN_MOLECULE_IDS


_HALF_LIFE = 7 * 24 * 3600  # 7 days in seconds


def theoretical_decay(t, n0):
	"""Exponential decay: N(t) = N0 * exp(-ln2 * t / t_half)."""
	return n0 * np.exp(-np.log(2) * t / _HALF_LIFE)


def main(out_dir='out/decay-verify'):
	print("Running 1-day platelet simulation...")
	result = run_platelet_sim(out_dir, length_sec=86400, seed=0)
	sim_out_dir = result['sim_out_dir']

	reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
	counts = reader.readColumn('counts')   # shape: (timesteps, n_molecules)
	ids = list(reader.readAttribute('objectNames'))
	n_steps = counts.shape[0]
	t = np.arange(n_steps)

	# ── Console report ────────────────────────────────────────────────────────
	p_per_step = 1 - np.exp(-np.log(2) / _HALF_LIFE)
	expected_pct = 100.0 * (1 - (1 - p_per_step) ** 86400)

	print(f"\n=== RestingDecay results (1 simulated day) ===")
	print(f"  {'Molecule':<25}  {'Start':>12}  {'End':>12}  {'Lost':>8}  {'Decay %':>8}  {'Tracked':>8}")
	print("  " + "-" * 82)
	for i, mol_id in enumerate(ids):
		c0, cf = int(counts[0, i]), int(counts[-1, i])
		lost = c0 - cf
		pct = 100.0 * lost / max(c0, 1)
		tracked = 'protein' if mol_id in PROTEIN_MOLECULE_IDS else 'lipid'
		print(f"  {mol_id:<25}  {c0:>12,}  {cf:>12,}  {lost:>8,}  {pct:>7.2f}%  {tracked:>8}")

	print(f"\n  Theoretical protein decay over 1 day: {expected_pct:.2f}%")
	print(f"  (half-life = 7 days; p per 1s step = {p_per_step:.2e})")

	# ── Plot (proteins only) ──────────────────────────────────────────────────
	fig, ax = plt.subplots(figsize=(9, 5))
	fig.suptitle('RestingDecay — simulated protein decay vs theory (1 day)', fontsize=13)

	protein_indices = [i for i, mol_id in enumerate(ids) if mol_id in PROTEIN_MOLECULE_IDS]

	for i in protein_indices:
		n0 = int(counts[0, i])
		ax.step(t / 3600, counts[:, i], where='post', alpha=0.8, label=f'Sim: {ids[i]}')
		t_theory = np.linspace(0, n_steps - 1, 500)
		ax.plot(t_theory / 3600, theoretical_decay(t_theory, n0=n0),
			color='black', linestyle='--', linewidth=1.5,
			label=f'Theory (N0={n0:,}, t½=7 days)')

	ax.set_xlabel('Simulated time (hours)')
	ax.set_ylabel('Molecule count')
	ax.legend()
	ax.set_xlim(0, (n_steps - 1) / 3600)
	plt.tight_layout()

	plot_out = os.path.join(os.path.dirname(sim_out_dir), 'plotOut', 'resting_decay.png')
	os.makedirs(os.path.dirname(plot_out), exist_ok=True)
	fig.savefig(plot_out, dpi=150)
	print(f"\nPlot saved to: {plot_out}")


if __name__ == '__main__':
	main(sys.argv[1] if len(sys.argv) > 1 else 'out/decay-verify')
