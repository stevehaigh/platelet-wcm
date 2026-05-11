"""
Diagnostic: deYoung-Keizer 1992 IP3R gating properties.

Computes m∞, h∞, Po, and IP3R Ca²⁺ flux at a grid of IP3 and Ca²⁺_cyt
concentrations, to verify that the new model gates appropriately at
resting (IP3 = 50 nM) vs stimulated (IP3 = 1–10 µM) conditions.

Replaces the former Sneyd-Dufour Markov-chain equilibrium checker
(which was removed with the Sneyd-Dufour → deYoung-Keizer migration,
lab-book-2026-05-11-dyk-ip3r-design.md, issue #27).

Usage:
    PYTHONPATH=$PWD python runscripts/manual/checkIP3REquilibrium.py
"""

import json
import os

from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	K_DYK,
	N_IP3R,
	GAMMA_IP3R_S,
	V_IM_V,
	RT_OVER_zF_V,
	NA_OVER_zF,
)


def m_inf(ip3_uM, ca_uM):
	return (ip3_uM / (ip3_uM + K_DYK['d1'])) * (ca_uM / (ca_uM + K_DYK['d5']))


def h_inf(ca_uM):
	return K_DYK['d2'] / (ca_uM + K_DYK['d2'])


def tau_h(ca_uM):
	return 1.0 / (K_DYK['a2'] * (ca_uM + K_DYK['d2']))


def po(ip3_uM, ca_uM):
	return (m_inf(ip3_uM, ca_uM) ** 4) * h_inf(ca_uM)


def main():
	print('deYoung-Keizer 1992 IP3R gating diagnostic')
	print('=' * 60)

	print(f'\nParameters:')
	for k, v in K_DYK.items():
		print(f'  {k} = {v}')
	print(f'  N_IP3R = {N_IP3R}')

	# ── Resting conditions ─────────────────────────────────────────────
	ip3_rest = 0.05    # µM (50 nM)
	ca_rest  = 0.10    # µM (100 nM)
	ca_dts   = 250.0   # µM (Dolan resting DTS)

	import math
	e_ca = RT_OVER_zF_V * math.log(ca_dts / ca_rest)
	driving = V_IM_V - e_ca
	flux_rest = -GAMMA_IP3R_S * N_IP3R * po(ip3_rest, ca_rest) * driving * NA_OVER_zF

	print(f'\nResting conditions  (IP3 = {ip3_rest*1000:.0f} nM, '
		  f'Ca_cyt = {ca_rest*1000:.0f} nM, Ca_DTS = {ca_dts:.0f} µM):')
	print(f'  m∞        = {m_inf(ip3_rest, ca_rest):.4f}')
	print(f'  h∞        = {h_inf(ca_rest):.4f}')
	print(f'  τ_h       = {tau_h(ca_rest):.1f} s')
	print(f'  Po        = {po(ip3_rest, ca_rest):.3e}')
	print(f'  IP3R flux = {flux_rest:.0f} ions/s  (into cytosol)')

	# ── IP3 dose–response at resting Ca²⁺ ─────────────────────────────
	print(f'\nIP3 dose–response (Ca_cyt = {ca_rest*1000:.0f} nM):')
	print(f'  {"IP3 (nM)":>10}  {"m∞":>8}  {"Po":>10}  {"flux (ions/s)":>15}')
	results = []
	for ip3_nM in [10, 50, 100, 500, 1000, 5000, 10000]:
		ip3_uM = ip3_nM / 1000.0
		p = po(ip3_uM, ca_rest)
		flux = -GAMMA_IP3R_S * N_IP3R * p * driving * NA_OVER_zF
		print(f'  {ip3_nM:>10}  {m_inf(ip3_uM, ca_rest):>8.4f}  {p:>10.3e}  {flux:>15.0f}')
		results.append({'ip3_nM': ip3_nM, 'm_inf': m_inf(ip3_uM, ca_rest),
						'po': p, 'flux_ions_s': flux})

	# ── Ca²⁺ dependence at stimulated IP3 ─────────────────────────────
	ip3_stim = 1.0   # µM (mid-stimulation)
	print(f'\nCa²⁺ dependence (IP3 = {ip3_stim*1000:.0f} nM):')
	print(f'  {"Ca_cyt (nM)":>12}  {"m∞":>8}  {"h∞":>8}  {"Po":>10}')
	for ca_nM in [50, 100, 200, 500, 1000, 2000]:
		ca_uM = ca_nM / 1000.0
		print(f'  {ca_nM:>12}  {m_inf(ip3_stim, ca_uM):>8.4f}  '
			  f'{h_inf(ca_uM):>8.4f}  {po(ip3_stim, ca_uM):>10.3e}')

	# ── Stimulated peak flux estimate ──────────────────────────────────
	ip3_peak = 0.275   # µM (5.5× resting = Dolan IP3 fold × 0.05 µM)
	ca_peak  = 0.40    # µM (400 nM, mid-transient estimate)
	e_ca_peak = RT_OVER_zF_V * math.log(ca_dts * 0.5 / ca_peak)
	driving_peak = V_IM_V - e_ca_peak
	flux_peak = -GAMMA_IP3R_S * N_IP3R * po(ip3_peak, ca_peak) * driving_peak * NA_OVER_zF
	print(f'\nEstimated peak transient flux')
	print(f'  (IP3 = {ip3_peak*1000:.0f} nM, Ca_cyt = {ca_peak*1000:.0f} nM, '
		  f'Ca_DTS = {ca_dts*0.5:.0f} µM):')
	print(f'  Po   = {po(ip3_peak, ca_peak):.3e}')
	print(f'  flux = {flux_peak:.0f} ions/s')

	# ── Save JSON ──────────────────────────────────────────────────────
	out = {
		'model': 'deYoung-Keizer 1992 / Li-Rinzel 1994',
		'parameters': K_DYK,
		'N_IP3R': N_IP3R,
		'resting': {
			'ip3_nM': ip3_rest * 1000,
			'ca_cyt_nM': ca_rest * 1000,
			'm_inf': m_inf(ip3_rest, ca_rest),
			'h_inf': h_inf(ca_rest),
			'tau_h_s': tau_h(ca_rest),
			'po': po(ip3_rest, ca_rest),
			'flux_ions_s': flux_rest,
		},
		'ip3_dose_response': results,
	}
	repo_root = os.path.join(os.path.dirname(__file__), '..', '..')
	out_path = os.path.normpath(os.path.join(
		repo_root, 'reports', 'data', 'ip3r-dyk-gating.json'))
	with open(out_path, 'w') as f:
		json.dump(out, f, indent=2)
	print(f'\nSaved: {out_path}')


if __name__ == '__main__':
	main()
