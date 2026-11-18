#!/bin/bash
vmodel_folders=$(find . -maxdepth 1 -type d -name "vmodel_*")
if [ -n "$vmodel_folders" ]; then
    echo "Deleting vmodel_* folders..."
    rm -rf vmodel_*
fi

path=$(pwd)
nm=$(echo "$path" | awk -F'/' '{print $(NF-1)}' | sed 's/^gf_//')
# after step2 folder replication inv_<event> is renamed inv_<event>_1sub;
# fall back so out-of-band GF reruns keep working
invd=../../inv_${nm}; [ -d "$invd" ] || invd=../../inv_${nm}_1sub
max_distance=$(saclst dist f  ${invd}/dataloc/*.z | sort -nrk2 | head -1 | gawk '{print (int($2+200) < 950) ? int($2+200) : 950}')
min_distance=$(saclst dist f  ${invd}/dataloc/*.z | sort -nk2 | head -1 | gawk '{print (int($2-200) > 2) ? int($2-200) : 2}')
evlo_evla=$(saclst evlo evla f ${invd}/dataloc/*.z | sort -nrk2 | head -1 | gawk '{print $2, $3}')
evlo=$(echo $evlo_evla | awk '{print $1}')
evla=$(echo $evlo_evla | awk '{print $2}')

# Generate 1D velocity model from CRUST1.0
# run in place: make_model.py reads CRUST1.0 relative to itself and writes
# vmodel.txt here, so concurrent events never share a scratch file
python ../../../programs/crust1.0/make_model.py $evlo $evla
# Velocity-model split (user rule 2026-08-09):
# - TELESEISMIC GFs (../vmodel, built below from vmodel.txt) RETAIN the full
#   crust1.0 column including water and marine sediments: they are SOURCE-side
#   structure (pP/sP/pwP travel through them; the hand-tuned California and
#   Kamchatka references both kept them, and stripping them degraded fits
#   ~50%). The tel receiver side is the separate hard-rock block.
# - LOCAL fk GFs (single model for source AND receivers) strip leading
#   Vs<2.0 layers: stations sit on hard rock, and a soft veneer/water column
#   misrepresents them. Stripped thickness is absorbed into the first hard
#   layer so interface/source depths keep the sea-surface reference and stay
#   consistent with the tel depth grid. (Replaces the old $2>0.1 water-only
#   filter - soft sediments are now stripped too.)
cat vmodel.txt | gawk '{ if (!done && $2 < 2.0) { acc += $4; next }
                         if (!done) { print $4+acc, $2, $1, $3; done = 1; next }
                         print $4, $2, $1, $3 }' > vmodel
nlen=$(cat vmodel.txt | wc -l)
# deep (intraslab, evdp>60 km) events: the single uniform 8.2 km/s x 250 km
# mantle layer biases depth phases (pP/sP ~2 s early, ~12% weak vs a graded
# model on 2024 Chile M7.4 -> source depths ~10 km too deep); use ak135-graded
# upper-mantle layers instead. Shallow events keep the legacy template.
# evdp from TELESEISMIC headers (dataloc headers can carry unset/-12345
# evdp, and sort|head takes the minimum -> deep branch never fires)
evdp0=$(saclst evdp f ${invd}/data/*.z | sort -nk2 | head -1 | gawk '{print int($2)}')
if [ -n "$evdp0" ] && [ "$evdp0" -gt 60 ]; then
cat vmodel.txt | gawk -v nl="$nlen" 'BEGIN{print "#vmodel";printf("%d ",nl+7)} {if ($4>0) {print $1,$2,$3,$4} else {print $1,$2,$3,25}} END{print "8.05 4.50 3.32 45\n8.18 4.51 3.36 90\n8.40 4.56 3.42 55\n8.56 4.62 3.46 60\n8.8 4.77 3.51 100\n9.7 5.28 3.88 250\n10.9 6.1 4.4 0\n3 5.8 3.36 2.72 20\n6.5 3.75 2.92 15\n8.04 4.47 3.32 0"}' > ../vmodel
else
cat vmodel.txt | gawk -v nl="$nlen" 'BEGIN{print "#vmodel";printf("%d ",nl+4)} {if ($4>0) {print $1,$2,$3,$4} else {print $1,$2,$3,25}} END{print "8.2 4.5 3.4 250\n8.8 4.77 3.51 100\n9.7 5.28 3.88 250\n10.9 6.1 4.4 0\n3 5.8 3.36 2.72 20\n6.5 3.75 2.92 15\n8.04 4.47 3.32 0"}' > ../vmodel
fi

# Calculate the maximum possible interval that results in less than 200 elements
interval=$(( (max_distance - min_distance) / 150 ))
if [ "$interval" -lt 1 ]; then
    interval=1
fi

evdp=$(saclst evdp f  ${invd}/data/*.z | sort -nk2 | head -1| gawk '{print int($2)}')
# event-adaptive depth aperture (../dotel.sh reads the grid from distdep.dat): deep
# intraslab ruptures scale with the Wells-Coppersmith length; shallow keep 20
mw=$(gawk -v ev="$nm" '$1==ev{print $4; exit}' ../../event_list.dat)
halfw=20
if [ "$evdp" -gt 60 ] && [ -n "$mw" ]; then
    halfw=$(gawk -v mw="$mw" 'BEGIN{h=int(0.65*10^(0.59*mw-2.44)+0.5); if(h<20)h=20; if(h>70)h=70; print h}')
fi
depth_min=$((evdp - halfw))
if [ "$depth_min" -lt 3 ]; then
    depth_min=3
fi
# The depth grid starts in rock: at least 0.5 km below the water and the soft sediments (Vs < 2.0
# km/s) of the CRUST1.0 column, the layers that are stripped from the regional model above and kept
# in the teleseismic one. The grid used to start at 3 km whatever the column held. The 2024
# California example has 2.75 km of water over 0.63 km of sediment (Vs 0.34 km/s), so its 3-km node
# was a source in that sediment layer: its teleseismic P Green's functions are 11 (dd) and 29 (ds)
# times those of the 5-km node and ring for 100 s with a period near 7.5 s, as the water
# reverberations of the data do (4 x 2.75 km / 1.5 km/s = 7.3 s), which a strike-slip source in rock
# hardly excites. Every v3 model of that example had its shallow subevents at 3.0-3.1 km, with
# mechanisms on the 30-degree cone (PITFALLS B5).
soft_column_km=$(gawk '{ if ($2 < 2.0) soft += $4; else exit } END { printf "%.2f", soft + 0 }' vmodel.txt)
first_rock_depth=$(gawk -v soft="$soft_column_km" 'BEGIN { depth = int(soft + 0.5); if (depth < soft + 0.5) depth++; print depth }')
if [ "$depth_min" -lt "$first_rock_depth" ]; then
    echo "doloc: the depth grid starts at $first_rock_depth km (water and soft sediments: $soft_column_km km)"
    depth_min=$first_rock_depth
fi
depth_max=$((evdp + halfw))
# NOTE: no Moho cap on the GF grid - megathrust interface events rupture
# below the local (oceanic) Moho; the crustal-depth cap for strike-slip
# events is applied at search time in get_search_par.py instead
diff=$(expr $depth_max - $depth_min)
remainder=$(expr $diff % 2)
if [ "$remainder" -ne 0 ]; then
    depth_max=$(expr $depth_max - 1)
fi
# SUBMIT_GF_INPUTS=<dir>: rebuild the Green's functions of a PUBLISHED run
# exactly - its velocity models (vmodel_tel, vmodel_loc) and its grid
# (distdep.dat: dist_min dist_max dist_step depth_min depth_max depth_step) -
# instead of deriving them with the current rules
if [ -n "$SUBMIT_GF_INPUTS" ] && [ -f "$SUBMIT_GF_INPUTS/distdep.dat" ]; then
    cp "$SUBMIT_GF_INPUTS/vmodel_loc" vmodel
    cp "$SUBMIT_GF_INPUTS/vmodel_tel" ../vmodel
    read min_distance max_distance interval depth_min depth_max _ddp < "$SUBMIT_GF_INPUTS/distdep.dat"
    echo "doloc: using published GF inputs from $SUBMIT_GF_INPUTS"
fi
ddp=2
###first step
gawk -v dp1="$depth_min" -v dp2="$depth_max" -v dd="$ddp" 'BEGIN{for (i=dp1;i<=dp2;i=i+dd) print "mkdir vmodel_"i}' | sh
gawk -v dp1="$depth_min" -v dp2="$depth_max" -v dd="$ddp" -v ddist="$interval" -v mind="$min_distance" -v md="$max_distance" 'BEGIN{for (i=dp1;i<=dp2;i=i+dd) print "sh genfk.sh "i,ddist,mind,md}' | sh

# cap concurrent fk jobs (override with SUBMIT_NCORE); previously unbounded "&"
NCORE=${SUBMIT_NCORE:-$(( $(nproc) / 2 > 0 ? $(nproc) / 2 : 1 ))}
initial_dir=$(pwd)
for dir in vmodel_*; do
    while [ "$(jobs -rp | wc -l)" -ge "$NCORE" ]; do sleep 1; done
    if [ -d "$dir" ]; then
        (
            cd "$dir" || exit  # Enter the folder
            cp ../vmodel .      # Copy the vmodel file
            sh fk.cmd         # Run fk.cmd in the background (parallel)
        ) &
    fi
done
wait
cd "$initial_dir" || exit

# generate distlst file
distlst=""
dd=$interval
md=$max_distance
mind=$min_distance
if [ "$mind" -gt 0 ] 2>/dev/null && [ "$md" -gt 0 ] 2>/dev/null; then
    # Populate the string with values from dd to md, incremented by dd
    i=$mind
    while [ $i -le $md ]; do
        distlst="$distlst $i"
        i=$((i + dd))
    done
    distlst=$(echo $distlst)
    dist_min=$(echo "$distlst" | tr ' ' '\n' | sort -n | head -n 1)
    dist_max=$(echo "$distlst" | tr ' ' '\n' | sort -n | tail -n 1)
    echo "$distlst" | tr ' ' '\n' > distlst
else
    echo "Both min_d and max_d must be positive integers."
fi
echo "$dist_min $dist_max $dd $depth_min $depth_max $ddp" > distdep.dat

###second step
ls -1 -d vmodel_* | gawk '{print "mv "$1"/"$1"/* "$1"/"}' | sh
#ls -1 -d vmodel_* | gawk '{print "cd "$1;print "cp ../transfer.sh ../dotransfer.sh ../distlst .";print "sh dotransfer.sh";print "cd ..";}' | sh

# Function to run commands in each vmodel_* directory
run_transfer() {
    dir="$1"
    (
        cd "$dir" || exit
        cp ../transfer.sh ../dotransfer.sh ../distlst .
        sh dotransfer.sh
    )
}

# Run transfer script in parallel
for dir in vmodel_*; do
    if [ -d "$dir" ]; then
        while [ "$(jobs -rp | wc -l)" -ge "$NCORE" ]; do sleep 1; done
        run_transfer "$dir" &
    fi
done
# Wait for all transfer operations to finish
wait
cd $initial_dir
#ls -1 -d vmodel_* | gawk '{print "rm "$1"/*.grn.?"}' | sh
