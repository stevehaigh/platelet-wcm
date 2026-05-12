"""
Diagnostic plot for the post-transient recovery dynamics.

Shows cyt + DTS Ca²⁺, IP3 + PIP2, and PMCA state distribution over a
long simulation (designed for 3000s+ runs). Helps diagnose why the
system is or isn't returning to resting state after a transient.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotLongRecovery.py [sim_outdir]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from wholecell.io.tablereader import TableReader


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('simOutDir', nargs='?',
		default='out/long-recovery/platelet_stub_000000/000000/generation_000000/000000/simOut')
	parser.add_argument('--out', default='reports/figures/long-recovery-2026-05-12.png')
	parser.add_argument('--stim-onset', type=float, default=60.0)
	args = parser.parse_args()

	ct = TableReader(os.path.join(args.simOutDir, 'CalciumTrace'))
	t = ct.readColumn('time')
	ca_cyt = ct.readColumn('ca_cyt_nM').flatten()
	ca_dts = ct.readColumn('ca_dts_uM').flatten()
	ip3 = ct.readColumn('ip3_nM').flatten()

	rb = TableReader(os.path.join(args.simOutDir, 'BulkMolecules'))
	ids = list(rb.readAttribute('objectNames'))
	counts = rb.readColumn('counts')
	n = min(len(t), counts.shape[0])
	t = t[:n]; ca_cyt = ca_cyt[:n]; ca_dts = ca_dts[:n]; ip3 = ip3[:n]
	def get(name):
		try: return counts[:n, ids.index(name)].astype(float)
		except ValueError: return np.zeros(n)

	pip2 = get('PIP2[c]')
	plcb_i = get('PLCb_inactive[c]')
	plcb_a = get('PLCb_active[c]')
	pmca_free = get('PMCA[pl]')
	pmca_ca = get('PMCA_Ca[pl]')
	cam_pmca = get('Ca4_CaM_PMCA[pl]')
	cam_pmca_ca = get('Ca4_CaM_PMCA_Ca[pl]')
	pmca_cam = get('PMCA_CaM[pl]')

	fig, axes = plt.subplots(5, 1, figsize=(11, 14), sharex=True)
	fig.subplots_adjust(hspace=0.30)

	stim_kwargs = dict(color='red', linestyle='--', alpha=0.5)

	# Panel 1: cyt Ca$^{2+}$
	axes[0].plot(t, ca_cyt, color='tab:blue', linewidth=2)
	axes[0].axhline(100.0, color='gray', linewidth=0.8, linestyle=':', label='Resting target (100 nM)')
	axes[0].axvline(args.stim_onset, **stim_kwargs)
	axes[0].set_ylabel('[Ca$^{2+}$]_cyt (nM)')
	axes[0].set_title('Cytosolic free Ca$^{2+}$ — long recovery')
	axes[0].set_yscale('log')
	axes[0].legend(loc='upper right', fontsize=9); axes[0].grid(alpha=0.3)

	# Panel 2: DTS Ca²⁺
	axes[1].plot(t, ca_dts, color='tab:blue', linewidth=2)
	axes[1].axhline(250.0, color='gray', linewidth=0.8, linestyle=':', label='Resting target (250 µM)')
	axes[1].axvline(args.stim_onset, **stim_kwargs)
	axes[1].set_ylabel('[Ca$^{2+}$]_DTS (µM)')
	axes[1].set_title('DTS free Ca$^{2+}$ — long recovery (peak then decline)')
	axes[1].legend(loc='upper right', fontsize=9); axes[1].grid(alpha=0.3)

	# Panel 3: IP3 + PIP2
	axes[2].plot(t, ip3, color='tab:purple', linewidth=2, label='IP$_3$ (model output)')
	axes[2].axhline(50.0, color='gray', linewidth=0.8, linestyle=':', label='Resting target (50 nM)')
	axes[2].axvline(args.stim_onset, **stim_kwargs)
	axes[2].set_ylabel('[IP$_3$] (nM)')
	axes[2].set_title('IP$_3$ (left) and PIP$_2$ (right) — does the PI cycle re-equilibrate?')
	ax2b = axes[2].twinx()
	ax2b.plot(t, pip2 / 1000.0, color='tab:cyan', linewidth=1.5, label='PIP$_2$')
	ax2b.axhline(112.0, color='tab:cyan', linewidth=0.8, linestyle=':', alpha=0.6)
	ax2b.set_ylabel('PIP$_2$ (× 1 000 count)', color='tab:cyan')
	ax2b.tick_params(axis='y', labelcolor='tab:cyan')
	lines1, l1 = axes[2].get_legend_handles_labels()
	lines2, l2 = ax2b.get_legend_handles_labels()
	axes[2].legend(lines1 + lines2, l1 + l2, loc='upper right', fontsize=9)
	axes[2].grid(alpha=0.3)

	# Panel 4: PLCβ active/inactive
	axes[3].plot(t, plcb_a, color='tab:orange', linewidth=2, label='PLCβ active')
	axes[3].plot(t, plcb_i, color='tab:gray', linewidth=1.5, label='PLCβ inactive')
	axes[3].axhline(143.0, color='tab:orange', linewidth=0.8, linestyle=':',
		alpha=0.6, label='Resting plcb_a (143)')
	axes[3].axvline(args.stim_onset, **stim_kwargs)
	axes[3].set_ylabel('PLCβ count')
	axes[3].set_title('PLCβ activation — should return to ~143 at rest after stim')
	axes[3].legend(loc='upper right', fontsize=9); axes[3].grid(alpha=0.3)

	# Panel 5: PMCA state distribution
	axes[4].stackplot(t,
		pmca_free, pmca_ca, cam_pmca, cam_pmca_ca, pmca_cam,
		labels=['PMCA free', 'PMCA·Ca²⁺ (basal active)',
				'Ca₄·CaM·PMCA', 'Ca₄·CaM·PMCA·Ca²⁺ (pumping)',
				'PMCA·CaM (trap)'],
		colors=['tab:green', 'tab:olive', 'tab:cyan', 'tab:orange', 'tab:red'],
		alpha=0.8)
	axes[4].axvline(args.stim_onset, **stim_kwargs)
	axes[4].set_xlabel('time (s)')
	axes[4].set_ylabel('PMCA molecules (of 769 total)')
	axes[4].set_title('PMCA state distribution — is recovery rate-limited by the CaM trap?')
	axes[4].legend(loc='upper right', fontsize=8); axes[4].grid(alpha=0.3)

	os.makedirs(os.path.dirname(args.out), exist_ok=True)
	fig.savefig(args.out, dpi=140, bbox_inches='tight')
	print(f'Saved: {args.out}')

	# Numerical summary at key timepoints
	def at(tp): return int(np.argmin(np.abs(t - tp)))
	print(f'\n{"t":>6} {"cyt nM":>8} {"DTS µM":>8} {"IP3 nM":>8} {"PIP2 k":>8} {"PLCb_a":>8} {"PMCA":>6} {"PMCA·CaM":>9}')
	for tp in [59, 200, 400, 600, 1000, 1500, 2000, 2500, 3000]:
		if tp <= t[-1]:
			i = at(tp)
			print(f'{t[i]:>6.0f} {ca_cyt[i]:>8.1f} {ca_dts[i]:>8.1f} {ip3[i]:>8.1f} '
				  f'{pip2[i]/1000:>8.1f} {plcb_a[i]:>8.0f} {pmca_free[i]:>6.0f} {pmca_cam[i]:>9.0f}')


if __name__ == '__main__':
	main()
