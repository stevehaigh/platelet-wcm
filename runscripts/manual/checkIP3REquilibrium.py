"""
Compute the equilibrium distribution of the 6-state Sneyd-Dufour IP3R
sub-state Markov chain at IP3 = 50 nM, Ca²⁺_cyt = 100 nM, with our
Purvis 2008 rate constants. Compare against the Dolan Table S1
initial conditions to test whether Dolan IC is a fixed point of our
ODE.

Phase 0 / Step 4 of lab-book-2026-05-07-phase-0-biology-audit.md.

Usage:
    PYTHONPATH=$PWD pyenv exec python runscripts/manual/checkIP3REquilibrium.py

Saves results to reports/data/ip3r-equilibrium-2026-05-07.json.
"""

import json
import os

import numpy as np

from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	K_IP3R,
	_phi_n_i1_fwd, _phi_n_i1_rev,
	_phi_n_o_fwd,  _phi_n_o_rev,
	_phi_o_a_fwd,  _phi_o_a_rev,
	_phi_a_i2_fwd, _phi_a_i2_rev,
	_phi_o_s_fwd,  _phi_o_s_rev,
)


# State indices (matches MOLECULE_NAMES ordering).
N, O, A, I1, I2, S = 0, 1, 2, 3, 4, 5
STATE_NAMES = ['n', 'o', 'a', 'i1', 'i2', 's']

# Conditions for the equilibrium check.
IP3_UM   = 0.05    # 50 nM (Dolan baseline)
CA_CYT_UM = 0.10   # 100 nM (Dolan rest target)

# Total IP3R count, Dolan Table S1 (sum across sub-states).
N_TOTAL_IP3R = 1327  # 809 + 167 + 261 + 65 + 25 + 0 ≈ Dolan reports 1328

# Dolan Table S1 sub-state counts (from
# reports/data/calcium-data-provenance.md).
DOLAN_TABLE_S1 = {
	'n':  809,
	'o':  261,
	'a':   65,
	'i1': 167,
	'i2':  25,
	's':    0,
}


def build_rate_matrix(ip3, ca):
	"""Build the 6x6 transition rate matrix Q.

	Q[i, j] = rate (s⁻¹) from state j to state i (for i ≠ j).
	Q[i, i] = -Σ_{j ≠ i} Q[j, i]  (diagonal balance).

	Equilibrium distribution p satisfies Q · p = 0 with Σp = 1.
	"""
	Q = np.zeros((6, 6))

	# Forward rates (j → i)
	pairs = [
		(N,  I1, _phi_n_i1_fwd(ca),       _phi_n_i1_rev()),
		(N,  O,  _phi_n_o_fwd(ip3, ca),   _phi_n_o_rev(ca)),
		(O,  A,  _phi_o_a_fwd(ca),        _phi_o_a_rev(ca)),
		(A,  I2, _phi_a_i2_fwd(ca),       _phi_a_i2_rev()),
		(O,  S,  _phi_o_s_fwd(ca),        _phi_o_s_rev()),
	]
	for j, i, fwd, rev in pairs:
		# j → i at rate fwd; i → j at rate rev
		Q[i, j] += fwd
		Q[j, i] += rev

	# Diagonal: column sums must be zero (mass conservation).
	for k in range(6):
		Q[k, k] = -(np.sum(Q[:, k]) - Q[k, k])

	return Q


def solve_equilibrium(Q):
	"""Find the stationary distribution p with Q·p = 0, Σp = 1.

	Replace one row of Q with the normalisation constraint and solve.
	"""
	A = Q.copy()
	A[5, :] = 1.0   # Σ p = 1
	b = np.zeros(6)
	b[5] = 1.0
	p = np.linalg.solve(A, b)
	return p


def po(p):
	"""Channel open probability per the Dolan / Purvis weighting:
	Po = 0.9·a/total + 0.1·o/total. Here p is already normalised so
	a/total = p[A], o/total = p[O]."""
	return 0.9 * p[A] + 0.1 * p[O]


def main():
	print(f'IP3R sub-state equilibrium at IP3 = {IP3_UM*1000:.0f} nM, '
		f'Ca²⁺_cyt = {CA_CYT_UM*1000:.0f} nM')
	print('=' * 70)

	# Print the rate constants for transparency.
	Q = build_rate_matrix(IP3_UM, CA_CYT_UM)
	print('\nTransition rate constants (s⁻¹) at given conditions:')
	print(f'  n → i1: {_phi_n_i1_fwd(CA_CYT_UM):8.4f}    i1 → n: {_phi_n_i1_rev():8.4f}')
	print(f'  n → o : {_phi_n_o_fwd(IP3_UM, CA_CYT_UM):8.4f}    o  → n: {_phi_n_o_rev(CA_CYT_UM):8.4f}')
	print(f'  o → a : {_phi_o_a_fwd(CA_CYT_UM):8.4f}    a  → o: {_phi_o_a_rev(CA_CYT_UM):8.4f}')
	print(f'  a → i2: {_phi_a_i2_fwd(CA_CYT_UM):8.4f}    i2 → a: {_phi_a_i2_rev():8.4f}')
	print(f'  o → s : {_phi_o_s_fwd(CA_CYT_UM):8.4f}    s  → o: {_phi_o_s_rev():8.4f}')

	# Solve for equilibrium.
	p_eq = solve_equilibrium(Q)
	dolan_total = sum(DOLAN_TABLE_S1.values())
	p_dolan = np.array([
		DOLAN_TABLE_S1['n']  / dolan_total,
		DOLAN_TABLE_S1['o']  / dolan_total,
		DOLAN_TABLE_S1['a']  / dolan_total,
		DOLAN_TABLE_S1['i1'] / dolan_total,
		DOLAN_TABLE_S1['i2'] / dolan_total,
		DOLAN_TABLE_S1['s']  / dolan_total,
	])

	print('\nSub-state distributions:')
	print(f'  {"state":<6}{"Dolan IC":>12}{"our equil.":>12}{"Δ (rel.)":>12}')
	for i, name in enumerate(STATE_NAMES):
		dolan_pct = p_dolan[i] * 100
		eq_pct = p_eq[i] * 100
		if dolan_pct > 0.01:
			delta = (eq_pct - dolan_pct) / dolan_pct
			delta_str = f'{delta:+8.2%}'
		else:
			delta_str = '  N/A   '
		print(f'  {name:<6}{dolan_pct:>11.2f}%{eq_pct:>11.2f}%   {delta_str}')

	po_dolan = po(p_dolan)
	po_eq = po(p_eq)
	print(f'\nPo (channel open prob):')
	print(f'  Po at Dolan IC      = {po_dolan:.4f}')
	print(f'  Po at our equil.    = {po_eq:.4f}')
	print(f'  Ratio (eq / Dolan)  = {po_eq / po_dolan:.2f}×')

	po4_dolan = po_dolan ** 4
	po4_eq = po_eq ** 4
	print(f'\nPo⁴ (4 IP3R-tetramer factor in flux):')
	print(f'  Po⁴ at Dolan IC     = {po4_dolan:.3e}')
	print(f'  Po⁴ at our equil.   = {po4_eq:.3e}')
	print(f'  Ratio (eq / Dolan)  = {po4_eq / po4_dolan:.1f}×')

	# Check whether dy/dt is meaningfully non-zero at Dolan IC, in counts/s.
	# This tells us the magnitude of the initial sub-state racing.
	dy_dolan = Q @ (p_dolan * N_TOTAL_IP3R)
	print(f'\nd(sub-state count)/dt at Dolan IC (count/s):')
	for i, name in enumerate(STATE_NAMES):
		print(f'  d({name})/dt = {dy_dolan[i]:+8.2f}')
	max_drift = np.max(np.abs(dy_dolan))
	print(f'  Max |drift| = {max_drift:.2f} count/s '
		f'({max_drift / N_TOTAL_IP3R * 100:.2f}%/s of total IP3R)')

	# Save JSON for the lab book.
	out = {
		'conditions': {'ip3_uM': IP3_UM, 'ca_cyt_uM': CA_CYT_UM},
		'states': STATE_NAMES,
		'dolan_table_s1_counts': DOLAN_TABLE_S1,
		'dolan_fractions': p_dolan.tolist(),
		'our_equilibrium_fractions': p_eq.tolist(),
		'po_dolan': float(po_dolan),
		'po_our_equilibrium': float(po_eq),
		'po4_dolan': float(po4_dolan),
		'po4_our_equilibrium': float(po4_eq),
		'po4_ratio_eq_over_dolan': float(po4_eq / po4_dolan),
		'd_substate_dt_at_dolan_ic_counts_per_s': dy_dolan.tolist(),
	}
	repo_root = os.path.join(os.path.dirname(__file__), '..', '..')
	out_path = os.path.normpath(os.path.join(
		repo_root, 'reports', 'data', 'ip3r-equilibrium-2026-05-07.json'))
	with open(out_path, 'w') as f:
		json.dump(out, f, indent=2)
	print(f'\nSaved: {out_path}')


if __name__ == '__main__':
	main()
