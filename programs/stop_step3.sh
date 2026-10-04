#!/bin/bash
# Safely stop a running step3_do_inversions.sh: its work runs inside an awk
# that calls system(), so killing the driver shell orphans the awk and its
# mpirun/finv children. This targets the actual process tree, verifying each
# PID's working directory against the given IRIS dir first (never kill by
# name alone - see AUDIT pkill incidents).
# usage: stop_step3.sh /abs/path/to/IRIS-dir
IRIS=$(readlink -f "${1:?usage: stop_step3.sh <IRIS-dir>}")
killed=0
for pid in $(pgrep -x awk; pgrep -x mpirun; pgrep -x mpiexec; pgrep -x finv); do
    cwd=$(readlink -f "/proc/$pid/cwd" 2>/dev/null) || continue
    case "$cwd" in
        "$IRIS"|"$IRIS"/*)
            echo "killing $(ps -o comm= -p $pid) pid $pid (cwd $cwd)"
            kill "$pid" 2>/dev/null && killed=$((killed+1))
            ;;
    esac
done
sleep 2
for pid in $(pgrep -x mpirun; pgrep -x finv); do
    cwd=$(readlink -f "/proc/$pid/cwd" 2>/dev/null) || continue
    case "$cwd" in
        "$IRIS"|"$IRIS"/*) kill -9 "$pid" 2>/dev/null ;;
    esac
done
echo "stop_step3: $killed process(es) signalled under $IRIS"
