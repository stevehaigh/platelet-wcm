"""Standalone rich.Live diagnostic — no platelet code, no sim data.

Runs four progressively-more-complex rich.Live demos for 5 s each and
then exits cleanly. Use to bisect the replayTui blank-cell bug:

    1. Inline (screen=False), single Text   — simplest
    2. Inline,                  Panel(Text)  — adds a border
    3. Inline,                  Layout       — adds the layout split
    4. Alt-screen (screen=True), Layout      — same shape as replayTui

A counter ticks up each second so you can tell if Live is actually
refreshing. If any mode shows blank where the previous worked, the
bug is in the *delta* between those two modes.

Usage:
    PYTHONPATH=$PWD pyenv exec python runscripts/manual/richDemo.py
    PYTHONPATH=$PWD pyenv exec python runscripts/manual/richDemo.py --mode 1
"""

from __future__ import annotations

import argparse
import sys
import time

from rich.console import Console
from rich.layout import Layout
from rich.live import Live
from rich.panel import Panel
from rich.text import Text


FRAMES = 5
FRAME_S = 1.0


def mode1() -> None:
	"""Live + a bare Text, no Panel, no Layout, inline (no alt-screen)."""
	console = Console()
	console.print('[bold]Mode 1[/bold]: bare Text, inline (screen=False)')
	with Live(Text('frame 0'), console=console, screen=False,
			auto_refresh=False) as live:
		for i in range(1, FRAMES + 1):
			time.sleep(FRAME_S)
			live.update(Text(f'frame {i} — should tick up'))
			live.refresh()
	console.print('Mode 1 done.\n')


def mode2() -> None:
	"""Live + Panel(Text), inline."""
	console = Console()
	console.print('[bold]Mode 2[/bold]: Panel(Text), inline (screen=False)')
	with Live(Panel(Text('frame 0'), title='counter'), console=console,
			screen=False, auto_refresh=False) as live:
		for i in range(1, FRAMES + 1):
			time.sleep(FRAME_S)
			live.update(Panel(Text(f'frame {i}\nline 2\nline 3'),
				title='counter'))
			live.refresh()
	console.print('Mode 2 done.\n')


def mode3() -> None:
	"""Live + Layout(header / content / footer), inline."""
	console = Console()
	console.print('[bold]Mode 3[/bold]: Layout (header/content/footer), inline')
	layout = Layout()
	layout.split_column(
		Layout(Panel('header'), name='h', size=3),
		Layout(Panel('content (this is the one that goes blank in replayTui)'),
			name='c', ratio=1),
		Layout(Panel('footer'), name='f', size=3),
	)
	with Live(layout, console=console, screen=False,
			auto_refresh=False) as live:
		for i in range(1, FRAMES + 1):
			time.sleep(FRAME_S)
			layout['c'].update(Panel(
				f'frame {i}\nline 2\nline 3\nline 4\nline 5',
				title='content'))
			live.refresh()
	console.print('Mode 3 done.\n')


def mode4() -> None:
	"""Live + Layout, alt-screen (screen=True) — same shape as replayTui."""
	console = Console()
	console.print('[bold]Mode 4[/bold]: Layout, alt-screen (screen=True)')
	console.print('(takes over the terminal for ~5 s, then restores)')
	time.sleep(1.0)
	layout = Layout()
	layout.split_column(
		Layout(Panel('header'), name='h', size=3),
		Layout(Panel('content'), name='c', ratio=1),
		Layout(Panel('footer'), name='f', size=3),
	)
	with Live(layout, console=console, screen=True,
			auto_refresh=False) as live:
		for i in range(1, FRAMES + 1):
			time.sleep(FRAME_S)
			layout['c'].update(Panel(
				f'frame {i}\nline 2\nline 3\nline 4\nline 5\n\n'
				'If you see this updating, alt-screen mode is fine.\n'
				"If this is blank, the bug is rich.Live's screen=True path "
				'in your terminal.',
				title='content (alt-screen)'))
			live.refresh()
	console.print('Mode 4 done. Did the content update each second?')


MODES = {1: mode1, 2: mode2, 3: mode3, 4: mode4}


def main(argv: list[str] | None = None) -> int:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--mode', type=int, default=None,
		help='Run one mode (1-4). Default = run all in sequence.')
	args = parser.parse_args(argv)

	if args.mode is not None:
		if args.mode not in MODES:
			print(f'unknown mode {args.mode}; choose 1..4', file=sys.stderr)
			return 2
		MODES[args.mode]()
	else:
		for n in (1, 2, 3, 4):
			MODES[n]()
	return 0


if __name__ == '__main__':
	sys.exit(main())
