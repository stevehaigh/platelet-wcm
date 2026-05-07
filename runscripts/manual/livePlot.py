"""
Live Ca²⁺ trace viewer — polls a CSV written by runPlateletSim --live.

Usage:
    # Watch a specific file:
    PYTHONPATH=$PWD pyenv exec python runscripts/manual/livePlot.py out/my_run/live.csv

    # Auto-discover the most recent live.csv under out/:
    PYTHONPATH=$PWD pyenv exec python runscripts/manual/livePlot.py

Start this before or during the simulation — it will display "Waiting..." until
data arrives and then update the plot once per second.
"""

import argparse
import csv
import glob
import os
import sys

import matplotlib

# Pick an interactive backend before importing pyplot. matplotlib otherwise
# defaults to Agg in headless contexts, which silently fails plt.show().
# Try macOS native first (works with stock framework Python), then Tk, then Qt.
def _select_interactive_backend():
	import os
	# Respect an explicit user choice via $MPLBACKEND only if it's interactive.
	current = os.environ.get('MPLBACKEND', '').lower()
	if current and current != 'agg':
		return matplotlib.get_backend()
	# Order matters: TkAgg first because the MacOSX backend has a known
	# issue where FuncAnimation timer callbacks fire but the canvas never
	# redraws, so the live plot looks frozen.
	for backend in ('TkAgg', 'QtAgg', 'Qt5Agg', 'MacOSX'):
		try:
			matplotlib.use(backend, force=True)
			return backend
		except (ImportError, ValueError):
			continue
	raise RuntimeError(
		'No interactive matplotlib backend available. Install tkinter or PyQt, '
		'or run on macOS with framework Python.')

_BACKEND = _select_interactive_backend()

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Dolan 2014 reference bands (Fig. 3B / Fig. 4).
DOLAN_CYT_MIN_NM  = 200.0
DOLAN_CYT_MAX_NM  = 800.0
DOLAN_DTS_MIN_UM  = 80.0
DOLAN_DTS_MAX_UM  = 250.0

POLL_INTERVAL_MS = 1000   # redraw interval


def find_latest_live_csv():
	"""Return the most recently modified live.csv under out/, or None."""
	repo_root = os.path.join(os.path.dirname(__file__), '..', '..')
	pattern = os.path.join(repo_root, 'out', '*', 'live.csv')
	files = glob.glob(pattern)
	if not files:
		return None
	return max(files, key=os.path.getmtime)


def read_csv(path):
	"""Read the live CSV, tolerating partially-written final rows."""
	t, ca_cyt, ca_dts, ip3, soce = [], [], [], [], []
	try:
		with open(path, 'r') as f:
			reader = csv.DictReader(f)
			for row in reader:
				try:
					t.append(float(row['time']))
					ca_cyt.append(float(row['ca_cyt_nM']))
					ca_dts.append(float(row['ca_dts_uM']))
					ip3.append(float(row['ip3_nM']))
					soce.append(float(row['soce_flux_nMs']))
				except (ValueError, KeyError):
					pass  # skip partial last row written mid-timestep
	except FileNotFoundError:
		pass
	return (np.asarray(t), np.asarray(ca_cyt), np.asarray(ca_dts),
			np.asarray(ip3), np.asarray(soce))


def build_figure():
	fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
	fig.subplots_adjust(hspace=0.08, top=0.93, bottom=0.08, left=0.1, right=0.97)

	ax_cyt, ax_dts, ax_soce = axes

	# Dolan reference bands.
	ax_cyt.axhspan(DOLAN_CYT_MIN_NM, DOLAN_CYT_MAX_NM,
				   color='steelblue', alpha=0.08, label='Dolan 2014 target band')
	ax_dts.axhspan(DOLAN_DTS_MIN_UM, DOLAN_DTS_MAX_UM,
				   color='green', alpha=0.08, label='Dolan 2014 DTS range')

	line_cyt, = ax_cyt.plot([], [], color='royalblue', lw=1.8)
	line_dts, = ax_dts.plot([], [], color='seagreen', lw=1.8)
	line_soce, = ax_soce.plot([], [], color='darkorange', lw=1.5)

	ax_cyt.set_ylabel('[Ca²⁺]_cyt  (nM)', fontsize=10)
	ax_dts.set_ylabel('[Ca²⁺]_DTS  (µM)', fontsize=10)
	ax_soce.set_ylabel('SOCE flux  (nM/s)', fontsize=10)
	ax_soce.set_xlabel('Simulation time  (s)', fontsize=10)

	ax_cyt.legend(fontsize=8, loc='upper right')
	ax_dts.legend(fontsize=8, loc='upper right')

	ax_cyt.set_ylim(bottom=0)
	ax_dts.set_ylim(bottom=0)

	# Title / status label.
	title = fig.suptitle('Waiting for simulation data…', fontsize=11,
						 fontweight='bold')

	return fig, axes, (line_cyt, line_dts, line_soce), title


def make_updater(live_path_ref, axes, lines, title):
	"""Return a FuncAnimation update function."""
	ax_cyt, ax_dts, ax_soce = axes
	line_cyt, line_dts, line_soce = lines

	def update(_frame):
		path = live_path_ref[0]
		if path is None:
			path = find_latest_live_csv()
			if path is None:
				return lines
			live_path_ref[0] = path

		t, ca_cyt, ca_dts, ip3, soce = read_csv(path)
		if len(t) == 0:
			return lines

		n = len(t)
		run_name = os.path.basename(os.path.dirname(path))
		title.set_text(f'{run_name}  —  t = {t[-1]:.0f} s  ({n} steps)')

		line_cyt.set_data(t, ca_cyt)
		line_dts.set_data(t, ca_dts)
		line_soce.set_data(t, soce)

		# Rescale x to data range.
		ax_soce.set_xlim(0, max(t[-1] + 5, 20))

		# Rescale y with 10 % headroom above the max, keep floor at 0.
		def rescale_y(ax, arr):
			ymax = float(np.max(arr)) if len(arr) else 1.0
			ax.set_ylim(0, max(ymax * 1.15, ax.get_ylim()[1]))

		rescale_y(ax_cyt, ca_cyt)
		rescale_y(ax_dts, ca_dts)
		# SOCE can be zero for long stretches — symmetric y axis looks better.
		soce_max = float(np.max(np.abs(soce))) if len(soce) else 1.0
		ax_soce.set_ylim(-soce_max * 0.1, max(soce_max * 1.15, 1.0))

		# Some backends (notably MacOSX) don't auto-redraw on set_data when
		# blit=False; force an idle redraw so the window actually updates.
		ax_cyt.figure.canvas.draw_idle()

		return lines

	return update


def main(argv=None):
	parser = argparse.ArgumentParser(
		description='Live Ca²⁺ trace viewer for platelet-wcm simulations.')
	parser.add_argument(
		'live_csv',
		nargs='?',
		default=None,
		help=('Path to the live.csv written by runPlateletSim --live. '
			  'If omitted, the most recently modified out/*/live.csv is used.'))
	args = parser.parse_args(argv)

	live_path = args.live_csv
	if live_path is None:
		live_path = find_latest_live_csv()
		if live_path is not None:
			print(f'Auto-discovered: {live_path}')
		else:
			print('No live.csv found yet; will watch for out/*/live.csv to appear.')

	fig, axes, lines, title = build_figure()

	# Use a mutable container so the closure can update the path on discovery.
	live_path_ref = [live_path]
	update = make_updater(live_path_ref, axes, lines, title)

	ani = animation.FuncAnimation(
		fig, update, interval=POLL_INTERVAL_MS, blit=False, cache_frame_data=False)

	plt.show()
	return ani  # keep reference so GC doesn't collect it


if __name__ == '__main__':
	ani = main()
