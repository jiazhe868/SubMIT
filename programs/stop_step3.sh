#!/bin/bash
# stop_step3.sh IRIS_DIR - stop a running step 3 (the driver loop and its mpirun/finv/ffwd children).
#
# Killing only the driver shell leaves mpirun and its ranks running, and vice versa the driver would
# start the next inversion. Every candidate PID is checked against the IRIS folder through
# /proc/PID/cwd before it is signalled: never kill by name alone (another event may be running).
#
# usage: sh ../programs/stop_step3.sh "$PWD"      (or ./submit-level: sh programs/stop_step3.sh work/<name>/IRIS)
iris_dir=$(readlink -f "${1:?usage: stop_step3.sh <IRIS-dir>}")
signalled_count=0
signal_processes_under_iris() {
    local signal=$1 pid working_dir
    for pid in $(pgrep -f "step3_do_inversions|calibrate_from_1sub"; pgrep -x awk; pgrep -x mpirun; pgrep -x mpiexec; pgrep -x finv; pgrep -x ffwd); do
        working_dir=$(readlink -f "/proc/$pid/cwd" 2>/dev/null) || continue
        case "$working_dir" in
            "$iris_dir"|"$iris_dir"/*)
                [ "$signal" = TERM ] && echo "stopping $(ps -o comm= -p "$pid") pid $pid (in $working_dir)"
                kill -s "$signal" "$pid" 2>/dev/null && signalled_count=$((signalled_count + 1)) ;;
        esac
    done
}
signal_processes_under_iris TERM
sleep 2
signal_processes_under_iris KILL
echo "stop_step3: $signalled_count signal(s) sent under $iris_dir"
