#!/bin/bash
# Stage-2 model-based station screen orchestrator (step3 hook, after the
# balance hook): full ffwd on the 1-sub best -> polarity/gain screens ->
# if anything excluded, refresh station info + Par counts in every inv dir
# and re-run the weight balancing on the cleaned set.
# usage: stage2_screen_from_1sub.sh <inv_..._1sub dir> <event_name>
d1=$1
ev=$2
fwd="fwd_${ev}"
[ -x "$fwd/ffwd" ] || { echo "stage2_screen: no ffwd in $fwd - skipped"; exit 0; }
# an unrepresentative 1-sub scan (giant events: duration railed or << Mw
# scaling) is a bad amplitude/CC reference - skip, same rule as tighten_bounds
best=""
for f in "$d1"/best_model_hybrid.dat "$d1"/best_model_ensemble.dat "$d1"/best_model_exploration.dat; do
    [ -s "$f" ] && best=$f && break
done
if [ -n "$best" ] && [ -s "$d1/mainshock.dat" ]; then
    dur1=$(gawk 'NR==1{printf "%.2f", $7}' "$best")
    mw1=$(gawk '{print $4}' "$d1/mainshock.dat")
    skip=$(gawk -v d="$dur1" -v mw="$mw1" 'BEGIN{print (d < 0.3*5*10^(0.5*(mw-6))) ? 1 : 0}')
    [ "$skip" = "1" ] && { echo "stage2_screen: 1-sub unrepresentative (dur $dur1 << Mw scaling) - skipped"; exit 0; }
fi
if [ "$(cat "$d1"/extended.flag 2>/dev/null)" = "1" ]; then
    echo "stage2_screen: extended/doublet - POLARITY screen only (gain ref unreliable)"
    export STAGE2_POLARITY_ONLY=1
fi
rm -f excluded_stations_model.txt
(cd "$fwd" && ./ffwd > /dev/null 2>&1)
[ -s "$fwd/sta_residual.dat" ] || { echo "stage2_screen: no sta_residual.dat - skipped"; exit 0; }
PYEXE=${SUBMIT_PYTHON:-$(command -v python3 || command -v python)}
case "$PYEXE" in */intel*|*oneapi*) PYEXE=/usr/bin/python3 ;; esac
$PYEXE ../programs/prepare_SubMIT/screen_stations_model.py "$ev" || { echo "stage2_screen: script failed - skipped"; exit 0; }
if [ -s excluded_stations_model.txt ]; then
    for d in inv_${ev}_*sub; do
        [ -d "$d" ] || continue
        (cd "$d" && sh prepare_stationinfo.sh > /dev/null 2>&1
         np=$(grep -c . stations.info); nsh=$(grep -c . stationsSH.info); nloc=$(grep -c . stationsloc.info)
         sed -i "s/^num_sta_P=.*/num_sta_P= $np #Number of all stations (P)/" Par.file
         sed -i "s/^num_sta_SH=.*/num_sta_SH= $nsh #Number of all stations (SH)/" Par.file
         sed -i "s/^num_sta_Rayl=.*/num_sta_Rayl= $nloc #Number of all stations (P)/" Par.file) || true
    done
    echo "stage2_screen: station sets changed - re-balancing weights"
    sh ../programs/prepare_SubMIT/balance_weights_from_1sub.sh "$d1" "$ev"
fi
