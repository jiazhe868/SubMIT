#!/bin/bash

original_dir=$(pwd)
# cores: override with SUBMIT_NCORE; use half the cores on very large nodes
# physical cores only (nproc/threads-per-core), capped at 32 by default:
# 32 chains are statistically sufficient, and rank counts at/near the core
# count let Intel MPI's busy-wait spinning starve straggler chains (observed
# 5x wall-time blowup at 64/64 under load). Override with SUBMIT_NCORE.
tpc=$(lscpu 2>/dev/null | awk -F: '/^Thread\(s\) per core/{gsub(/ /,"",$2); print $2}')
[ -n "$tpc" ] && [ "$tpc" -ge 1 ] 2>/dev/null || tpc=1
phys=$(( $(nproc) / tpc ))
[ "$phys" -gt 32 ] && phys=32
ncore=${SUBMIT_NCORE:-$phys}
[ "$ncore" -ge 1 ] || ncore=1
# finished ranks sleep instead of burning a core while waiting for stragglers
export I_MPI_WAIT_MODE=1
if [ -z "$1" ]; then
    echo "No argument provided. Exiting."
    exit 1
else
    nn=$1
fi
nn=$1

# Source the Intel environment
# find an Intel setvars if one exists (env INTEL_SETVARS wins); otherwise rely
# on generic MPI wrappers already in PATH (Makefile auto-detects the compiler)
# try Intel env candidates until one provides a WORKING ifx (stale installs may
# lack it); if none does, the Makefile falls back to generic MPI wrappers
# remember the user's python BEFORE the Intel env can put intelpython first on
# PATH; the 1-subevent hooks (model-based screening) inherit it
export SUBMIT_PYTHON=${SUBMIT_PYTHON:-$(command -v python3 || command -v python)}
# SUBMIT_NO_INTEL=1 skips this (e.g. when building with GNU mpicc/mpif90:
# an Intel mpirun cannot launch an OpenMPI/MPICH-built binary)
[ "${SUBMIT_NO_INTEL:-0}" = "1" ] && INTEL_CANDIDATES="" || \
    INTEL_CANDIDATES="$INTEL_SETVARS /opt/intel/oneapi/setvars.sh $HOME/intel/oneapi/setvars.sh"
for sv in $INTEL_CANDIDATES; do
    [ -n "$sv" ] && [ -f "$sv" ] || continue
    . "$sv" --force > /dev/null 2>&1
    if ! command -v ifx > /dev/null 2>&1; then
        # some installs do not export the compiler via setvars: source it directly
        for cv in $(ls -1d "$(dirname "$sv")"/compiler/*/env/vars.sh 2>/dev/null | sort -rV); do
            [ -f "$cv" ] && . "$cv" > /dev/null 2>&1 && command -v ifx > /dev/null 2>&1 && break
        done
    fi
    if command -v ifx > /dev/null 2>&1; then
        echo "Intel env OK: $sv (ifx: $(command -v ifx))"
        break
    fi
    echo "Intel env $sv has no working ifx, trying next candidate"
done
command -v mpirun > /dev/null || { echo "ERROR: mpirun not in PATH"; exit 1; }

# Read and process the file using awk
awk -v original_dir="$original_dir" -v nn="$nn" -v ncore="$ncore" '
{
    event_name = $1
    mw = $4
    nsub = nn
    while (nsub <= nn) {
        target_dir = "inv_" event_name "_" nsub "sub"
        if (system("[ -d " target_dir " ]") == 0) {
            # Create the runmpi.sh script in the target directory
            run_script = target_dir "/runmpi.sh"
            print "#!/bin/bash" > run_script
            print "make clean" >> run_script
            print "make || { echo COMPILE FAILED in " target_dir "; exit 1; }" >> run_script
            print "[ -x ./finv ] || { echo finv missing in " target_dir "; exit 1; }" >> run_script
            print "rm -f *best.dat best_model_*.dat *chain.dat mpi.out" >> run_script
            print "mpirun -np " ncore " ./finv > mpi.out 2>&1" >> run_script
            print "ls *best.dat > /dev/null 2>&1 || { echo INVERSION PRODUCED NO OUTPUT in " target_dir "; exit 1; }" >> run_script
            print "mode=exploration; [ \"$SUBMIT_RESEED\" = 1 ] && mode=ensemble; [ \"$SUBMIT_RESEED\" = 2 ] && mode=hybrid" >> run_script
            print "sort -gk3 *best.dat | head -1 > best_model_${mode}.dat" >> run_script

            # Make the runmpi.sh script executable
            system("chmod +x " run_script)

            # Run the runmpi.sh script
            system("cd " target_dir " && bash runmpi.sh")

            print "Executed runmpi.sh in", target_dir
        }
        nsub++
    }
}' event_list.dat

