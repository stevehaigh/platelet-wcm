"""
Verify RestingDecay process: run a 1-day platelet sim and report molecule counts.

Usage:
    PYTHONPATH="$PWD" python runscripts/manual/verify_resting_decay.py
"""

import os
import numpy as np

from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.io.tablereader import TableReader


def main():
	print("Running 1-day platelet simulation...")
	result = run_platelet_sim('out/decay-verify', length_sec=86400, seed=0)
	sim_out_dir = result['sim_out_dir']

	reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
	counts = reader.readColumn('counts')
	ids = reader.readAttribute('objectNames')

	print("\n=== RestingDecay results (1 simulated day) ===")
	print(f"  {'Molecule':<25}  {'Start':>6}  {'End':>6}  {'Lost':>6}  {'Decay %':>8}")
	print("  " + "-" * 58)
	for i, mol_id in enumerate(ids):
		c0, cf = int(counts[0, i]), int(counts[-1, i])
		lost = c0 - cf
		pct = 100.0 * lost / max(c0, 1)
		print(f"  {mol_id:<25}  {c0:>6}  {cf:>6}  {lost:>6}  {pct:>7.2f}%")

	# Theoretical expectation
	t_half = 7 * 24 * 3600   # 7 days in seconds
	p_per_step = 1 - np.exp(-np.log(2) / t_half)
	expected_pct = 100.0 * (1 - (1 - p_per_step) ** 86400)
	print(f"\n  Theoretical decay over 1 day: {expected_pct:.2f}%")
	print(f"  (half-life = 7 days; p per 1s step = {p_per_step:.2e})")


if __name__ == '__main__':
	main()
