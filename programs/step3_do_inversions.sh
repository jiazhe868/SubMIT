#!/bin/bash
# INTERRUPT WARNING: the per-nsub work below runs inside an awk process that
# calls system() - killing this driver script (or its shell) ORPHANS the awk,
# which keeps orchestrating mpirun jobs. To stop a run, target the awk (and
# its mpirun/finv children), not the driver: use ../programs/stop_step3.sh
# <IRIS-dir>, which matches processes by /proc/PID/cwd before killing.

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
# node-shared Green's functions (POSIX shm): one physical copy + one disk
# read per node instead of per rank (validated bit-identical vs per-rank,
# 39% faster at 8 ranks on Calama; ~20 GB RAM saved at 32 ranks). Opt out
# with SUBMIT_SHM_GF=0.
export SUBMIT_SHM_GF=${SUBMIT_SHM_GF:-1}

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

# GF-coverage preflight: sub_init's brsac read of a MISSING Green's function
# is silent (garbage kernel for that station biases every inversion - 2024
# Chile AU.MAW, fix 31). Verify every station distance in the 1sub info
# files has GF files before burning compute; fail loudly with the list.
gfmiss=0
while read -r event_id rest; do
    i1="inv_${event_id}_1sub"
    gdir=$(ls -d gf_*/greenFuncDir_disp 2>/dev/null | head -1)
    [ -d "$i1" ] && [ -n "$gdir" ] || continue
    vd=$(ls -d "$gdir"/vmodel_* 2>/dev/null | head -1)
    [ -n "$vd" ] || continue
    for f in stations.info stationsSH.info; do
        [ -s "$i1/$f" ] || continue
        for dist in $(gawk '{print $NF}' "$i1/$f" | sort -n | uniq); do
            ls "$vd/$dist".grn.* > /dev/null 2>&1 || {
                echo "ERROR: no Green's functions at distance $dist (listed in $i1/$f)"
                echo "       regenerate GFs (gf_SubMIT/dotel.sh now unions .z+.t distances)"
                echo "       or remove the station from the info files + Par num_sta counts"
                gfmiss=1
            }
        done
    done
done < event_list.dat
[ "$gfmiss" -eq 0 ] || { echo "step3 ABORTED: GF coverage incomplete"; exit 1; }

# Read and process the file using awk
awk -v original_dir="$original_dir" -v ncore="$ncore" '
{
    event_name = $1
    mw = $4
    nsub = 1
    while (nsub <= 9) {
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
            # after the 1-sub scan, tighten the n>=2 search bounds to the
            # empirically observed total source time (cen1+dura1): kills late
            # noise-absorber subevents that a magnitude-scaled cen_max admits
            # (2024 Chile M7.4 deep event: actual ~16 s vs Mw-scaled 31 s).
            # No-op when the 1-sub duration railed its cap (giant events).
            if (nsub == 1) {
                system("sh ../programs/prepare_SubMIT/tighten_bounds_from_1sub.sh " target_dir " " event_name)
                system("sh ../programs/prepare_SubMIT/balance_weights_from_1sub.sh " target_dir " " event_name)
                system("sh ../programs/prepare_SubMIT/stage2_screen_from_1sub.sh " target_dir " " event_name)
                # the hooks retune windows and weights for N>=2 only; re-run
                # N=1 under that SAME configuration so every point of the
                # L-curve is computed on identical data (otherwise N=1 is
                # scored on longer windows / other weights and the 5% rule
                # compares incomparable misfits). The first pass is kept.
                n2 = "inv_" event_name "_2sub"
                if (system("[ -f " n2 "/Par.file ]") == 0) {
                    system("f=$(ls -t " target_dir "/best_model_*.dat 2>/dev/null | head -1); " \
                           "[ -n \"$f\" ] && cp \"$f\" " target_dir "/calibration_pass_best.dat; " \
                           "cp " n2 "/Par.file " target_dir "/Par.file && cd " target_dir " && bash runmpi.sh")
                    print "Re-ran " target_dir " under the calibrated N>=2 configuration"
                }
            }
        }
        nsub++
    }
}' event_list.dat

