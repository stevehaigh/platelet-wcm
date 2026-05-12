"""
Diagnostic plot: cytosolic and DTS Ca$^{2+}$ — free vs bound, over time.

Reads an existing simulation output (CalciumTrace + BulkMolecules) and
produces a 2-panel figure:
  Panel 1 — cytosol: free Ca$^{2+}$, CaM-bound Ca$^{2+}$ (2·Ca₂CaM + 4·Ca₄CaM +
            complex sub-states), GSN-bound Ca$^{2+}$.
  Panel 2 — DTS: free Ca$^{2+}$, STIM1-bound Ca$^{2+}$.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotCaBoundFree.py [sim_outdir]

If no sim dir is given, the most recent run under out/ is used.
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from wholecell.io.tablereader import TableReader


def find_latest_simout():
	candidates = sorted(
		glob.glob('out/*/platelet_stub_*/*/generation_000000/000000/simOut'),
		key=os.path.getmtime,
	)
	if not candidates:
		sys.exit('No sim output found under out/')
	return candidates[-1]


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('simOutDir', nargs='?', default=None,
		help='Path to simOut/ (default: most recent under out/)')
	parser.add_argument('--out', default='reports/figures/ca-bound-free.png')
	parser.add_argument('--ip3-stim-onset', type=float, default=60.0,
		help='Time of IP3 stimulus onset for annotation (s)')
	args = parser.parse_args()

	simOutDir = args.simOutDir or find_latest_simout()
	print(f'Reading: {simOutDir}')

	# CalciumTrace — free Ca$^{2+}$ concentrations
	ct = TableReader(os.path.join(simOutDir, 'CalciumTrace'))
	t = ct.readColumn('time')
	ca_cyt_nM = ct.readColumn('ca_cyt_nM').flatten()
	ca_dts_uM = ct.readColumn('ca_dts_uM').flatten()
	ip3_nM    = ct.readColumn('ip3_nM').flatten()

	# BulkMolecules — bound-state counts
	rb = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
	mol_ids = list(rb.readAttribute('objectNames'))
	all_counts = rb.readColumn('counts')
	n = min(len(t), all_counts.shape[0])
	t = t[:n]; ca_cyt_nM = ca_cyt_nM[:n]; ca_dts_uM = ca_dts_uM[:n]
	ip3_nM = ip3_nM[:n]; all_counts = all_counts[:n]

	def get(name):
		try:
			return all_counts[:, mol_ids.index(name)].astype(float)
		except ValueError:
			return np.zeros(n)

	# ── Cytosol bound Ca$^{2+}$ accounting ─────────────────────────────────
	# CaM ladder: each Ca₂·CaM carries 2 Ca$^{2+}$, each Ca₄·CaM carries 4.
	# PMCA-CaM complex sub-states are membrane-localised but carry Ca$^{2+}$
	# that originated from the cytosol; include them in the "bound" total
	# for the cyt-side mass balance.
	ca2_cam         = get('Ca2_CaM[c]')
	ca4_cam         = get('Ca4_CaM[c]')
	ca4_cam_pmca    = get('Ca4_CaM_PMCA[pl]')    # 4 Ca$^{2+}$ on the CaM
	ca4_cam_pmca_ca = get('Ca4_CaM_PMCA_Ca[pl]') # 4 on CaM + 1 on PMCA pumping site
	pmca_ca         = get('PMCA_Ca[pl]')         # basal PMCA·Ca$^{2+}$
	gsn_ca          = get('GSN_Ca[c]')           # gelsolin proxy

	cam_bound_ions = 2.0*ca2_cam + 4.0*ca4_cam + 4.0*ca4_cam_pmca \
	                 + 5.0*ca4_cam_pmca_ca + pmca_ca
	gsn_bound_ions = gsn_ca

	# Convert cyt counts to nM (same factor used in CalciumTrace).
	# CalciumTrace gives ca_cyt_nM as the *free* concentration; recover free
	# count, then express bound pools in the same nM scale.
	# Free count = ca_cyt_nM × (cyt_volume × Avogadro × 1e-9) ; we don't need
	# the constant explicitly — derive it from row 0.
	free_count_0 = ca_cyt_nM[0]  # at t=0
	# But ca_cyt_nM is already a concentration; we need the conversion factor.
	# From calcium_signalling.py: cyt count of 361 ions = 100 nM.
	NM_PER_COUNT_CYT = 100.0 / 361.0
	cam_bound_nM = cam_bound_ions * NM_PER_COUNT_CYT
	gsn_bound_nM = gsn_bound_ions * NM_PER_COUNT_CYT

	# ── DTS bound Ca$^{2+}$ accounting ─────────────────────────────────────
	stim1_ca = get('STIM1_Ca[dts]')     # 1 Ca$^{2+}$ per STIM1
	stim1_dim = get('STIM1_dim[dts]')   # dimer; Ca-bound form
	calr_ca = get('CALR_Ca[dts]')       # CALR C-domain low-affinity sites
	calr_p_ca = get('CALR_P_Ca[dts]')   # CALR P-domain high-affinity slow site
	hsp90_m_ca = get('HSP90B1_M_Ca[dts]')  # HSP90B1 medium-affinity (slow)
	hsp90_l_ca = get('HSP90B1_L_Ca[dts]')  # HSP90B1 low-affinity (fast)
	bip_ca   = get('BiP_Ca[dts]')          # BiP/HSPA5
	crec_ca  = get('CREC_Ca[dts]')         # CALU+RCN1+RCN2 pool

	# DTS volume = 4.3% × 6 fL = 0.258 fL.  From internal_state.py:
	#   38842 ions = 250 µM  →  µM per count = 250/38842
	UM_PER_COUNT_DTS = 250.0 / 38842.0
	stim_bound_uM   = (stim1_ca + 2.0*stim1_dim) * UM_PER_COUNT_DTS
	calr_bound_uM   = calr_ca   * UM_PER_COUNT_DTS
	calr_p_bound_uM = calr_p_ca * UM_PER_COUNT_DTS
	hsp90_m_uM      = hsp90_m_ca * UM_PER_COUNT_DTS
	hsp90_l_uM      = hsp90_l_ca * UM_PER_COUNT_DTS
	bip_uM          = bip_ca     * UM_PER_COUNT_DTS
	crec_uM         = crec_ca    * UM_PER_COUNT_DTS
	additional_dts_uM = hsp90_m_uM + hsp90_l_uM + bip_uM + crec_uM

	# ── Plot — 4 panels: cyt-free, cyt-bound, DTS, IP3 ─────────────────
	fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 13), sharex=True)
	fig.subplots_adjust(hspace=0.32)

	stim_kwargs = dict(color='red', linestyle='--', alpha=0.6)
	stim_label = f'IP3 onset (t = {args.ip3_stim_onset:.0f} s)'

	# Panel 1 — cytosol FREE only (true scale)
	ax1.plot(t, ca_cyt_nM, color='tab:blue', linewidth=2.0, label='Free Ca$^{2+}$_cyt')
	ax1.axhline(100.0, color='gray', linewidth=0.8, linestyle=':',
		label='Dolan resting target (100 nM)')
	ax1.axvline(args.ip3_stim_onset, **stim_kwargs)
	ax1.set_ylabel('[Ca$^{2+}$]_cyt (nM)')
	ax1.set_title('Cytosolic free Ca$^{2+}$' + f'   |   peak = {ca_cyt_nM.max():.0f} nM @ t = {t[np.argmax(ca_cyt_nM)]:.0f} s')
	ax1.legend(loc='upper right', fontsize=9)
	ax1.grid(alpha=0.3)

	# Panel 2 — cytosol BOUND (larger scale)
	total_bound = cam_bound_nM + gsn_bound_nM
	ax2.plot(t, cam_bound_nM,
		label=r'CaM-bound (2$\cdot$Ca$_2$CaM + 4$\cdot$Ca$_4$CaM + complexes)',
		color='tab:orange', linewidth=2)
	ax2.plot(t, gsn_bound_nM, label='GSN-bound (scaffold; 50× under-sized)',
		color='tab:green', linewidth=2)
	ax2.plot(t, total_bound, label='Total bound', color='black',
		linewidth=1.2, linestyle=':')
	ax2.axvline(args.ip3_stim_onset, **stim_kwargs)
	ax2.set_ylabel('Bound Ca$^{2+}$ (nM equiv)')
	ratio_rest = total_bound[int(args.ip3_stim_onset) - 1] / ca_cyt_nM[int(args.ip3_stim_onset) - 1]
	ratio_peak = total_bound[np.argmax(ca_cyt_nM)] / ca_cyt_nM.max()
	ax2.set_title(
		'Cytosolic bound Ca$^{2+}$   |   bound:free ratio  '
		+ f'rest ≈ {ratio_rest:.1f}:1,  peak ≈ {ratio_peak:.1f}:1  '
		+ '(literature ~50:1)'
	)
	ax2.legend(loc='upper right', fontsize=9)
	ax2.grid(alpha=0.3)

	# Panel 3 — DTS (free + bound on one axis; both are µM-scale)
	ax3.plot(t, ca_dts_uM, label='Free Ca$^{2+}$_DTS', color='tab:blue', linewidth=2)
	ax3.plot(t, calr_bound_uM, label='CALR C-domain (fast)', color='tab:red', linewidth=2)
	ax3.plot(t, hsp90_m_uM, label='HSP90B1 medium (slow)', color='tab:green', linewidth=2)
	ax3.plot(t, hsp90_l_uM + bip_uM + crec_uM,
		label='HSP90B1-L + BiP + CREC (fast)', color='tab:olive', linewidth=2)
	ax3.plot(t, calr_p_bound_uM, label='CALR P-domain (slow)', color='tab:purple', linewidth=2)
	ax3.plot(t, ca_dts_uM, label='Free Ca$^{2+}$_DTS', color='tab:blue', linewidth=2)
	ax3.plot(t, stim_bound_uM, label='STIM1-bound', color='tab:orange', linewidth=1.5)
	ax3.plot(t, ca_dts_uM + stim_bound_uM + calr_bound_uM + calr_p_bound_uM
				+ hsp90_m_uM + hsp90_l_uM + bip_uM + crec_uM,
		label='Total DTS Ca$^{2+}$', color='black', linewidth=1.2, linestyle=':')
	ax3.axvline(args.ip3_stim_onset, **stim_kwargs)
	ax3.set_ylabel('[Ca$^{2+}$]_DTS (µM)')
	ax3.set_title('DTS Ca$^{2+}$ — multi-buffer (CALR C+P, HSP90B1 M+L, BiP, CREC)')
	ax3.legend(loc='upper right', fontsize=9)
	ax3.grid(alpha=0.3)

	# Panel 4 — IP3 (drive)
	ax4.plot(t, ip3_nM, color='tab:purple', linewidth=2, label='IP$_3$')
	ax4.axvline(args.ip3_stim_onset, **stim_kwargs)
	ax4.set_xlabel('time (s)')
	ax4.set_ylabel('[IP$_3$] (nM)')
	ax4.set_title('IP$_3$ forcing curve (Dolan 2014 Fig. S2 fit, delayed)')
	ax4.legend(loc='upper right', fontsize=9)
	ax4.grid(alpha=0.3)

	os.makedirs(os.path.dirname(args.out), exist_ok=True)
	fig.savefig(args.out, dpi=140, bbox_inches='tight')
	print(f'Saved: {args.out}')

	# Numeric summary
	idx_pre  = np.searchsorted(t, args.ip3_stim_onset - 1)
	idx_peak = np.argmax(ca_cyt_nM)
	idx_end  = -1
	print(f'\nSummary:')
	print(f'  Pre-stim (t={t[idx_pre]:.0f} s):  '
		f'cyt free = {ca_cyt_nM[idx_pre]:.1f} nM, '
		f'CaM-bound = {cam_bound_nM[idx_pre]:.1f} nM equiv, '
		f'GSN-bound = {gsn_bound_nM[idx_pre]:.1f} nM equiv, '
		f'DTS free = {ca_dts_uM[idx_pre]:.1f} µM, '
		f'STIM-bound = {stim_bound_uM[idx_pre]:.1f} µM')
	print(f'  Peak     (t={t[idx_peak]:.0f} s):  '
		f'cyt free = {ca_cyt_nM[idx_peak]:.1f} nM, '
		f'CaM-bound = {cam_bound_nM[idx_peak]:.1f} nM equiv, '
		f'GSN-bound = {gsn_bound_nM[idx_peak]:.1f} nM equiv, '
		f'DTS free = {ca_dts_uM[idx_peak]:.1f} µM')
	print(f'  End      (t={t[idx_end]:.0f} s):   '
		f'cyt free = {ca_cyt_nM[idx_end]:.1f} nM, '
		f'CaM-bound = {cam_bound_nM[idx_end]:.1f} nM equiv, '
		f'GSN-bound = {gsn_bound_nM[idx_end]:.1f} nM equiv, '
		f'DTS free = {ca_dts_uM[idx_end]:.1f} µM')


if __name__ == '__main__':
	main()
