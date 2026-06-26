#!/usr/bin/env bash
# Quick launcher for the platelet WCM terminal replay TUI.
#
# Defaults to the Phase 3 +Ca²⁺ simOut from `out/phase3_issue44_final/` at
# 0.2× speed (slow-mo). Both can be overridden positionally:
#
#   ./runscripts/manual/launch-tui.sh                     # defaults
#   ./runscripts/manual/launch-tui.sh out/.../simOut       # custom sim
#   ./runscripts/manual/launch-tui.sh out/.../simOut 1.0   # custom sim + speed
#
# Designed to be paste-and-go in iTerm or any other real terminal.
# Requires the uv-managed Python (pinned to 3.11.5 via .python-version)
# with deps installed (`uv sync --all-extras`).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

DEFAULT_SIMOUT="out/phase3_issue44_final/with_ca/platelet_stub_000000/000000/generation_000000/000000/simOut"
SIMOUT="${1:-$DEFAULT_SIMOUT}"
SPEED="${2:-0.2}"

if [[ ! -d "$SIMOUT" ]]; then
	echo "error: simOut not found at $SIMOUT" >&2
	echo "  Run 'PYTHONPATH=\$PWD uv run python runscripts/manual/runPhase3.py phase3_issue44_final --length 200' first," >&2
	echo "  or pass a different simOut directory as the first argument." >&2
	exit 1
fi

# Ensure the viz extras (rich + textual) are importable. One-shot; installs
# the viz extra so versions match the lock.
if ! uv run python -c "import textual, rich" 2>/dev/null; then
	echo "viz extras (rich/textual) not installed; installing now ..." >&2
	uv sync --extra viz >/dev/null
fi

echo "Launching replay TUI in your terminal..."
echo "  simOut : $SIMOUT"
echo "  speed  : ${SPEED}× (real-time = 1.0; 0.2 = slow-mo)"
echo
echo "Keys inside the TUI:"
echo "  q          quit"
echo "  space      pause / resume"
echo "  + / -      speed up / slow down"
echo "  ← / →      step back / forward 1 s (pauses)"
echo "  r          restart"
echo
echo "(Footer at the bottom of the TUI auto-renders the same hint.)"
echo

exec env PYTHONPATH="$REPO_ROOT" uv run python \
	"$REPO_ROOT/runscripts/manual/replayTui.py" "$SIMOUT" --speed "$SPEED"
