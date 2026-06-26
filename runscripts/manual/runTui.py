"""
Launch the platelet TUI experiment bench (Textual).

    OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" uv run python runscripts/manual/runTui.py

or simply:

    make tui

Both run under the uv-managed 3.11.5 interpreter (.venv), which has the sim
deps plus `textual` / `textual-plotext`. Plain `python3` may resolve to a
different interpreter that lacks them.

Edit the run conditions, press `r` (or click Run) to launch a simulation as
a subprocess, and watch the Ca²⁺ trace stream live. See
`reports/design/tui-tinkering-dashboard-2026-06-15.qmd`.
"""

from wholecell.tui.app import PlateletBenchApp


def main() -> None:
	"""Run the platelet TUI."""
	PlateletBenchApp().run()


if __name__ == '__main__':
	main()
