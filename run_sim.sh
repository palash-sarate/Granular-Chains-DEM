#!/usr/bin/env bash
set -Eeuo pipefail

# Paths and files
BIN="${HOME}/LIGGGHTS-PUBLIC/src/lmp_auto"
INPUT_FILE="in.ball_drop"
LOG_FILE="log.liggghts"

# Pre-flight checks
if [[ ! -x "$BIN" ]]; then
	echo "Error: LIGGGHTS binary not found or not executable at: $BIN" >&2
	exit 1
fi
if [[ ! -f "$INPUT_FILE" ]]; then
	echo "Error: Input file not found: $INPUT_FILE" >&2
	exit 1
fi

# Start streaming the log to the terminal while the simulation runs.
# We don't truncate the existing log; we follow only new lines.
touch "$LOG_FILE"
stdbuf -oL -eL tail -n 0 -f "$LOG_FILE" &
TAIL_PID=$!

cleanup() {
	# Stop tailing when the simulation finishes or on interruption
	kill "$TAIL_PID" >/dev/null 2>&1 || true
}
trap cleanup EXIT INT TERM

"$BIN" < "$INPUT_FILE"
