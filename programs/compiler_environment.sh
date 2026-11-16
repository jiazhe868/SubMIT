# compiler_environment.sh - MPI rank count and compilers for steps 3 and 4 (sourced, not run).
#
# Why a shared file: the same 40 lines were copied into step3_do_inversions.sh,
# step3_do_inversions_tempnsub.sh and step4_fwd_tempnsub.sh, and drifted apart.
#
# Sets: submit_rank_count (MPI ranks = chains), SUBMIT_PYTHON (the user's python, remembered BEFORE
# the Intel environment can put intelpython first on PATH), I_MPI_WAIT_MODE, SUBMIT_SHM_GF, and the
# Intel oneAPI compilers when they are installed and work (SUBMIT_NO_INTEL=1 skips them: an Intel
# mpirun cannot launch an OpenMPI/MPICH binary). Without Intel the Makefile uses mpicc/mpif90.
#
# usage (in a step script):  . ../programs/compiler_environment.sh

# physical cores only, capped at 32: 32 chains suffice, and rank counts at the core count let Intel
# MPI's busy-wait spinning starve straggler chains (5x wall-time blow-up observed at 64 of 64)
threads_per_core=$(lscpu 2>/dev/null | awk -F: '/^Thread\(s\) per core/{gsub(/ /,"",$2); print $2}')
[ -n "$threads_per_core" ] && [ "$threads_per_core" -ge 1 ] 2>/dev/null || threads_per_core=1
physical_core_count=$(( $(nproc) / threads_per_core ))
[ "$physical_core_count" -gt 32 ] && physical_core_count=32
submit_rank_count=${SUBMIT_NCORE:-$physical_core_count}
[ "$submit_rank_count" -ge 1 ] 2>/dev/null || submit_rank_count=1

export SUBMIT_PYTHON=${SUBMIT_PYTHON:-$(command -v python3 || command -v python)}
# finished ranks sleep instead of spinning while the slowest chain finishes
export I_MPI_WAIT_MODE=1
# one node-shared copy of the Green's functions (forward_model.c); SUBMIT_SHM_GF=0 opts out
export SUBMIT_SHM_GF=${SUBMIT_SHM_GF:-1}

# Intel's setvars.sh and vars.sh read unset variables, and the step scripts run under set -u, which
# would stop them inside the sourced file: nounset is lifted while the Intel scripts run
case $- in *u*) restore_nounset=1 ;; *) restore_nounset=0 ;; esac
set +u
if [ "${SUBMIT_NO_INTEL:-0}" != "1" ]; then
    for intel_setvars in ${INTEL_SETVARS:-} /opt/intel/oneapi/setvars.sh $HOME/intel/oneapi/setvars.sh; do
        [ -f "$intel_setvars" ] || continue
        . "$intel_setvars" --force > /dev/null 2>&1
        # some installs do not export the compiler through setvars: source its own vars.sh
        if ! command -v ifx > /dev/null 2>&1; then
            for compiler_vars in $(ls -1d "$(dirname "$intel_setvars")"/compiler/*/env/vars.sh 2>/dev/null | sort -rV); do
                . "$compiler_vars" > /dev/null 2>&1 && command -v ifx > /dev/null 2>&1 && break
            done
        fi
        if command -v ifx > /dev/null 2>&1; then
            echo "Intel environment: $intel_setvars (ifx: $(command -v ifx))"
            break
        fi
        echo "Intel environment $intel_setvars has no working ifx, trying the next one"
    done
fi
[ "$restore_nounset" = 1 ] && set -u
command -v mpirun > /dev/null || { echo "ERROR: mpirun not in PATH"; exit 1; }
