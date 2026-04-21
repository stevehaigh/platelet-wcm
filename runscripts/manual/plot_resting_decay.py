"""
Plot RestingDecay verification results.

Reads the simulation output produced by verify_resting_decay.py and generates:
  1. Actual simulated molecule counts over time.
  2. Theoretical exponential decay curve (scaled to a realistic copy number)
     to illustrate the expected 7-day half-life.

Usage:
    PYTHONPATH="$PWD" python runscripts/manual/plot_resting_decay.py [sim_out_dir]

If sim_out_dir is omitted, defaults to out/decay-verify/simOut.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader


# 7-day half-life in seconds
_HALF_LIFE = 7 * 24 * 3600

# Representative platelet protein copy number for the theoretical curve.
# Based on Burkhart 2012 median protein abundance (~5000 copies / platelet).
_REPRESENTATIVE_COUNT = 5000


def theoretical_decay(t, n0=_REPRESENTATIVE_COUNT):
	"""Exponential decay: N(t) = N0 * 2^(-t / t_half)."""
	return n0 * np.exp(-np.log(2) * t / _HALF_LIFE)


def main(sim_out_dir):
	reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
	counts = reader.readColumn('counts')   # shape: (timesteps, n_molecules)
	ids = reader.readAttribute('objectNames')
	n_steps = counts.shape[0]
	t = np.arange(n_steps)  # seconds

	fig, ax = plt.subplots(figsize=(9, 5))
	fig.suptitle('RestingDecay process — 1-day simulation vs theory', fontsize=13)

	# Simulated counts
	for i, mol_id in enumerate(ids):
		ax.step(t / 3600, counts[:, i], where='post', alpha=0.8, label=f'Sim: {mol_id}')

	# Theoretical decay over the same 1-day window, scaled to starting count
	n0 = counts[0].mean()
	t_theory = np.linspace(0, n_steps - 1, 500)
	ax.plot(t_theory / 3600, theoretical_decay(t_theory, n0=n0),
		color='black', linestyle='--', linewidth=1.5, label=f'Theory (N0={int(n0):,}, t½=7 days)')

	ax.set_xlabel('Simulated time (hours)')
	ax.set_ylabel('Molecule count')
	ax.legend()
	ax.set_xlim(0, (n_steps - 1) / 3600)

	plt.tight_layout()
	out_path = os.path.join(os.path.dirname(sim_out_dir), 'plotOut', 'resting_decay.png')
	os.makedirs(os.path.dirname(out_path), exist_ok=True)
	fig.savefig(out_path, dpi=150)
	print(f"Plot saved to: {out_path}")
	return out_path


if __name__ == '__main__':
	sim_out = sys.argv[1] if len(sys.argv) > 1 else 'out/decay-verify/simOut'
	main(sim_out)
