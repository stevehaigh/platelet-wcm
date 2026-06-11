"""
Regenerate ca2-transient-reference.png.

Run from the repo root:
    PYTHONPATH=$PWD python reports/figures/gen_ca2_transient_reference.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import os

# ── Curve parameters (match calcium_trace.py and calcium_signalling.py) ──────
_BASELINE_NM   = 100.0
_PEAK_NM       = 400.0
_PLATEAU_NM    = 200.0
_T_PEAK_S      = 5.0
_TAU_RISE_S    = 2.0
_TAU_DECAY_S   = 30.0
_TAU_PLATEAU_S = 120.0

_IP3_REST_NM   = 50.0
_IP3_FOLD      = 5.5
_IP3_T_PEAK    = 3.0
_IP3_TAU_RISE  = 3.0
_IP3_TAU_DECAY = 60.0


def _ca_transient_nM(t):
	t = np.asarray(t, dtype=float)
	rise         = 1.0 - np.exp(-np.maximum(t, 0.0) / _TAU_RISE_S)
	decay        = np.exp(-np.maximum(t - _T_PEAK_S, 0.0) / _TAU_DECAY_S)
	plateau_rise = 1.0 - np.exp(-np.maximum(t, 0.0) / _TAU_PLATEAU_S)
	transient = (_PEAK_NM - _BASELINE_NM) * rise * decay
	plateau   = (_PLATEAU_NM - _BASELINE_NM) * plateau_rise * decay
	return _BASELINE_NM + transient + plateau


def _ip3_nM(t):
	t = np.asarray(t, dtype=float)
	rise  = 1.0 - np.exp(-np.maximum(t, 0.0) / _IP3_TAU_RISE)
	decay = np.exp(-np.maximum(t - _IP3_T_PEAK, 0.0) / _IP3_TAU_DECAY)
	return _IP3_REST_NM * (1.0 + (_IP3_FOLD - 1.0) * rise * decay)


t = np.linspace(0, 300, 1000)
ca = _ca_transient_nM(t)
ip3 = _ip3_nM(t)

peak_t = t[np.argmax(ca)]
peak_ca = ca.max()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
fig.suptitle(
	r'Expected Ca$^{2+}$ Transient — Schematic (target: Dolan & Diamond 2014 Fig. 4)',
	fontsize=11,
)

# ── Upper panel: cytosolic Ca²⁺ ──────────────────────────────────────────────
ax1.plot(t, ca, color='steelblue', linewidth=2)
ax1.axhline(_BASELINE_NM, color='gray', linestyle='--', linewidth=1,
			label='Resting (100 nM)')
ax1.axhline(_PLATEAU_NM + _BASELINE_NM / 10, color='tomato', linestyle=':',
			linewidth=1.2, label=r'SOCE plateau (~210 nM)')
ax1.annotate(
	f'Peak\n~{peak_ca:.0f} nM',
	xy=(peak_t, peak_ca),
	xytext=(peak_t + 25, peak_ca - 70),
	arrowprops=dict(arrowstyle='->', color='black', lw=0.9),
	fontsize=8.5,
)
ax1.text(90, 222, 'SOCE-sustained plateau', color='tomato', fontsize=7.5)
ax1.set_ylabel(r'Cytosolic Ca$^{2+}$ (nM)')
ax1.set_ylim(bottom=0)
ax1.legend(fontsize=8, loc='upper right')

# ── Lower panel: IP3 ─────────────────────────────────────────────────────────
ax2.plot(t, ip3, color='forestgreen', linewidth=2)
ax2.axhline(_IP3_REST_NM, color='gray', linestyle='--', linewidth=1,
			label='Resting (50 nM)')
ax2.set_ylabel(r'IP3 (nM)')
ax2.set_xlabel('Time (s)')
ax2.text(0.01, 0.96, r'IP3 forcing function (§3.2)', transform=ax2.transAxes,
		 va='top', fontsize=9)
ax2.legend(fontsize=8, loc='upper right')

fig.tight_layout()

out = os.path.join(os.path.dirname(__file__), 'ca2-transient-reference.png')
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f'Saved {out}')
