"""
Re-derive the resting initial conditions by integrating the calcium ODE
with no IP3 forcing for 600 s, starting from the current
``internal_state.py`` initial counts. The converged state is a true
fixed point of the model — no more sub-state racing in the first
timestep, no more startup spike from a non-equilibrium IP3R Markov
chain.

This is the 19a deliverable (issue #19): re-derive the sub-state IC
without changing biology. Output JSON is intended to be pasted
into ``reconstruction/platelet/dataclasses/internal_state.py``.

Usage:
    PYTHONPATH=$PWD pyenv exec python runscripts/manual/restConvergence.py \\
        [--length 600] [--out reports/data/rest-converged-2026-05-07.json]
"""

import argparse
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	_ode_rhs,
	MOLECULE_NAMES,
	N_SPECIES,
	_UM_PER_COUNT_CYT,
	_UM_PER_COUNT_DTS,
)
from reconstruction.platelet.simulation_data import SimulationDataPlatelet


DEFAULT_LENGTH_SEC = 600.0
DEFAULT_OUT = 'reports/data/rest-converged-2026-05-07.json'

# Species the ODE cares about (subset of full sim_data molecule list).
ODE_SPECIES = list(MOLECULE_NAMES)


def initial_counts_vector(sim_data):
	"""Pull the 27 ODE species' starting counts from sim_data, keeping ODE order."""
	all_ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
	all_counts = sim_data.internal_state.bulk_molecules.initial_counts
	y0 = np.zeros(N_SPECIES)
	for i, name in enumerate(ODE_SPECIES):
		y0[i] = float(all_counts[all_ids.index(name)])
	return y0


def main(argv=None):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--length', type=float, default=DEFAULT_LENGTH_SEC,
		help=f'Convergence integration length in seconds (default {DEFAULT_LENGTH_SEC:.0f}).')
	parser.add_argument('--out', default=DEFAULT_OUT,
		help=f'Output JSON path (default {DEFAULT_OUT}).')
	args = parser.parse_args(argv)

	sim_data = SimulationDataPlatelet()
	y0 = initial_counts_vector(sim_data)

	# Integrate the ODE with no IP3 forcing — IP3 stays at its initial value.
	# t_sim_start is unused when ip3_forced=False.
	print(f'Integrating calcium ODE for {args.length:.0f} s with no IP3 forcing…')
	sol = solve_ivp(
		_ode_rhs, (0.0, args.length), y0,
		method='BDF',
		args=(0.0, False, 0.0),     # t_sim_start, ip3_forced, ip3_delay
		atol=1e-3, rtol=1e-6,
		max_step=10.0,
	)
	if not sol.success:
		raise RuntimeError(f'Convergence integration failed: {sol.message}')

	y_final = np.maximum(sol.y[:, -1], 0.0)

	# Convergence check: dy/dt at the final state should be ~ 0.
	dy_final = _ode_rhs(args.length, y_final, 0.0, False, 0.0)
	max_drift = float(np.max(np.abs(dy_final)))

	# Per-species comparison.
	print('\nConverged state vs. initial:')
	print(f'  {"species":<28}{"initial":>10}{"converged":>12}{"|dy/dt|":>12}')
	rows = []
	for i, name in enumerate(ODE_SPECIES):
		initial = int(round(y0[i]))
		converged = int(round(y_final[i]))
		drift = float(dy_final[i])
		print(f'  {name:<28}{initial:>10}{converged:>12}{drift:>+12.3e}')
		rows.append({
			'species': name,
			'initial_count': initial,
			'converged_count': converged,
			'd_count_per_s_at_converged': drift,
		})

	# Macroscale concentrations.
	idx = {name: i for i, name in enumerate(ODE_SPECIES)}
	ca_cyt_uM = y_final[idx['CA2_CYT[c]']] * _UM_PER_COUNT_CYT
	ca_dts_uM = y_final[idx['CA2_DTS[dts]']] * _UM_PER_COUNT_DTS
	ip3_uM    = y_final[idx['IP3[c]']]      * _UM_PER_COUNT_CYT
	print('\nMacroscale concentrations at convergence:')
	print(f'  Ca²⁺_cyt = {ca_cyt_uM*1e3:>8.3f} nM   (Dolan target: 100 nM)')
	print(f'  Ca²⁺_DTS = {ca_dts_uM:>8.3f} µM   (Dolan target: 200–300 µM)')
	print(f'  IP3      = {ip3_uM*1e3:>8.3f} nM   (Dolan baseline: 50 nM)')

	# Sub-state mass-conservation check (totals should be preserved).
	def total(prefix):
		return sum(y_final[idx[s]] for s in ODE_SPECIES if s.startswith(prefix))

	print('\nMass conservation check (total per protein):')
	totals = {
		'IP3R':    sum(y_final[idx[s]] for s in ODE_SPECIES if s.startswith('IP3R_')),
		'SERCA':   sum(y_final[idx[s]] for s in ODE_SPECIES if s.startswith('SERCA_')),
		'PMCA':    (y_final[idx['PMCA[pl]']] + y_final[idx['PMCA_Ca[pl]']]
		          + y_final[idx['Ca4_CaM_PMCA[pl]']]
		          + y_final[idx['Ca4_CaM_PMCA_Ca[pl]']]
		          + y_final[idx['PMCA_CaM[pl]']]),
		'CaM':     (y_final[idx['CaM_free[c]']] + y_final[idx['Ca2_CaM[c]']]
		          + y_final[idx['Ca4_CaM[c]']]
		          + y_final[idx['Ca4_CaM_PMCA[pl]']]
		          + y_final[idx['Ca4_CaM_PMCA_Ca[pl]']]
		          + y_final[idx['PMCA_CaM[pl]']]),
		# STIM1 total monomers = free + STIM1_Ca + 2·dimers
		'STIM1':   (y_final[idx['STIM1_free[dts]']] + y_final[idx['STIM1_Ca[dts]']]
		          + 2.0 * y_final[idx['STIM1_dim[dts]']]),
	}
	expected = {'IP3R': 1328, 'SERCA': 11892, 'PMCA': 769, 'CaM': 20481, 'STIM1': 4265}
	for k, v in totals.items():
		exp = expected[k]
		dev = (v - exp) / exp * 100
		print(f'  {k:<8} converged total = {v:>10.1f}  '
			f'(expected {exp}; Δ {dev:+.3f}%)')

	# Convergence verdict.
	print(f'\nMax |dy/dt| at converged state: {max_drift:.4f} count/s')
	if max_drift < 1.0:
		print('✓ Convergence target met (max drift < 1 count/s)')
	else:
		print('✗ Convergence target NOT met — increase --length or check ODE')

	repo_root = os.path.join(os.path.dirname(__file__), '..', '..')
	out_path = os.path.normpath(os.path.join(repo_root, args.out))
	out = {
		'integration_length_s': args.length,
		'max_drift_count_per_s': max_drift,
		'converged_concentrations': {
			'ca_cyt_nM': ca_cyt_uM * 1e3,
			'ca_dts_uM': ca_dts_uM,
			'ip3_nM':    ip3_uM * 1e3,
		},
		'totals': {k: float(v) for k, v in totals.items()},
		'species': rows,
	}
	with open(out_path, 'w') as f:
		json.dump(out, f, indent=2)
	print(f'\nSaved: {out_path}')


if __name__ == '__main__':
	main()
