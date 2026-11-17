#!/bin/bash
# calibrate_from_1sub.sh - calibrate the N >= 2 inversions on the 1-subevent result (step 3 hook).
#
# Three stages, in order (they were three scripts before v3: tighten_bounds_from_1sub.sh,
# balance_weights_from_1sub.sh, stage2_screen_from_1sub.sh):
#
# 1. bounds and windows from the 1-subevent source time. The magnitude-scaled duration
#    5*10^((Mw-6)/2) overestimates compact ruptures (2024 Chile M7.4 intraslab: ~16 s against 25 s,
#    so cen_max 31 admitted noise-absorbing subevents at 26-30 s and a non-monotone L-curve):
#      cen_max  <- min(cen_max,  ceil(cen1 + dur1 + 6))
#      dura_max <- min(dura_max, ceil(dur1 + 5))
#      P/SH window end <- cen_max + dura_max + 10 s + depth-phase allowance 2z/4 (deep events),
#      rounded UP to 5 s. The rounding is new in v3: without it a 1-subevent result differing by a
#      second between two runs changed every N >= 2 configuration, and the L-curve with it.
#    Skipped for extended/doublet sources, when the 1-subevent duration railed its cap (giant events),
#    and when it is below 30% of the magnitude scaling (a stuck scan: California M7.0 dur1 = 3.0
#    would have forbidden the real ~9-s subevents). Bounds only ever shrink.
# 2. wave-type weights balanced so the misfit shares are P:SH:regional = 1.2 : 1 : 1, iterating
#    w <- w * sqrt(target/measured) on the forward model of the 1-subevent best model (at most 4
#    times, to within 15%), renormalized to weightSH = 1. A regional set of <= 2 stations is a filler
#    and keeps its weight (P and SH balance to 1.2 : 1 alone).
# 3. station screen on the forward run of the 1-subevent model (screen_stations_model.py):
#    anti-correlation with the prediction (polarity: 2024 Chile GT.DBIC, CN.DRLN, GE.SNAA), amplitude
#    off by more than x4 from the band median of obs/syn (gain: CN.GAC), and a record whose rms in the
#    band and window of the inversion is more than 8 times the median of its group (2026 Vanuatu
#    IM.IL31, a drifting record that carried 68% of the P misfit). Only the polarity test and the
#    record test run for extended sources. When the 1-subevent result is unrepresentative (its
#    duration below 30% of the magnitude scaling) the two tests against the prediction are skipped;
#    the record test needs no model and still runs (until 2026-11-17 the whole screen was skipped
#    then, which left that record in). Excluded traces leave every inversion folder; station lists
#    and weights are then rebuilt.
#
# Exit status: 0 when the calibration ran or a stage was skipped for a stated reason, 1 when a stage
# failed half-way (step 3 then stops instead of running N >= 2 on inconsistent configurations).
#
# usage (in IRIS/, called by step3_do_inversions.sh): calibrate_from_1sub.sh inv_<event>_1sub <event>
set -u
one_subevent_dir=$1
event_name=$2
forward_dir="fwd_${event_name}"
best_model="$one_subevent_dir/best_model_hybrid.dat"
python_executable=${SUBMIT_PYTHON:-$(command -v python3 || command -v python)}
case "$python_executable" in */intel*|*oneapi*) python_executable=/usr/bin/python3 ;; esac
[ -s "$best_model" ] || { echo "calibrate: no 1-subevent best model - skipped"; exit 0; }

centroid_1=$(gawk 'NR==1{printf "%.2f", $4}' "$best_model")
duration_1=$(gawk 'NR==1{printf "%.2f", $7}' "$best_model")
duration_cap_1=$(gawk '$1=="dura_max"{print $2}' "$one_subevent_dir/search_par.file")
magnitude=$(gawk '{print $4}' "$one_subevent_dir/mainshock.dat" 2>/dev/null)
hypocentre_depth=$(gawk '{print $3}' "$one_subevent_dir/mainshock.dat" 2>/dev/null)
is_extended=$(cat "$one_subevent_dir/extended.flag" 2>/dev/null || echo 0)
duration_railed=$(gawk -v d="$duration_1" -v m="$duration_cap_1" 'BEGIN{print (d >= 0.95*m) ? 1 : 0}')
duration_too_short=$(gawk -v d="$duration_1" -v mw="${magnitude:-0}" 'BEGIN{print (mw > 0 && d < 0.3*5*10^(0.5*(mw-6))) ? 1 : 0}')

# ------------------------------------------------------------------------------------------------
# 1. bounds and windows
# ------------------------------------------------------------------------------------------------
if [ "$is_extended" = "1" ]; then
    echo "calibrate bounds: extended/doublet source - skipped"
elif [ "$duration_railed" -eq 1 ]; then
    echo "calibrate bounds: 1-subevent duration $duration_1 railed dura_max $duration_cap_1 - skipped"
elif [ "$duration_too_short" -eq 1 ]; then
    echo "calibrate bounds: 1-subevent duration $duration_1 << magnitude scaling (stuck scan?) - skipped"
else
    for inversion_dir in inv_"${event_name}"_*sub; do
        case "$inversion_dir" in *_1sub) continue ;; esac
        search_file="$inversion_dir/search_par.file"
        [ -s "$search_file" ] || continue
        gawk -v c1="$centroid_1" -v d1="$duration_1" '
            function ceil(x) { return (x == int(x)) ? x : int(x) + 1 }
            $1=="cen_max"  { limit = ceil(c1 + d1 + 6); if (limit < $2) $2 = limit }
            $1=="dura_max" { limit = ceil(d1 + 5);      if (limit < $2) $2 = limit }
            { print }' "$search_file" > "$search_file.tmp" && mv "$search_file.tmp" "$search_file"
        centroid_cap=$(gawk '$1=="cen_max"{print $2}' "$search_file")
        duration_cap=$(gawk '$1=="dura_max"{print $2}' "$search_file")
        # the window must hold the last subevent's depth phases, not only its direct pulse: sS-S ~ 2z/4
        # (~60 s at 120 km); shallow events keep a 5-s allowance
        window_end=$(gawk -v c="$centroid_cap" -v d="$duration_cap" -v z="${hypocentre_depth:-0}" 'BEGIN{
            depth_phase = (z > 60) ? 2*z/4.0 : 5; v = c + d + 10 + depth_phase
            v = 5 * int((v + 4.999) / 5); if (v > 300) v = 300; printf "%d", v }')
        for key in nd_timeP nd_timeSH; do
            sed -i "s/^${key}=.*/${key}= $window_end #ending time of inverse time window (cen_max+dura_max+10+depth-phase 2z\/4, rounded up to 5 s)/" "$inversion_dir/Par.file"
        done
        echo "calibrate bounds: $inversion_dir cen_max $centroid_cap dura_max $duration_cap windows $window_end s"
    done
fi

# ------------------------------------------------------------------------------------------------
# 2. wave-type weights (needs ffwd in the forward template folder)
# ------------------------------------------------------------------------------------------------
build_forward_program() {
    [ -d "$forward_dir" ] || { echo "calibrate: no $forward_dir - weights and screen skipped"; return 1; }
    cp ../programs/code_SubMIT/{forward_model.c,ffwd.f90,finv.f90,mtdcmp.f,sacsubc.c,sacsubc.h,Makefile} "$forward_dir/"
    ( cd "$forward_dir" && make clean > /dev/null 2>&1 && make ffwd > make_calibration.log 2>&1 )
    [ -x "$forward_dir/ffwd" ] || { echo "calibrate: ffwd build failed (see $forward_dir/make_calibration.log)"; return 1; }
    # the forward template is a step-2 copy: refresh its configuration from the 1-subevent folder
    for configuration_file in Par.file stations.info stationsSH.info stationsloc.info weights.dat distdep.dat search_par.file totalmt.dat; do
        if [ -s "$one_subevent_dir/$configuration_file" ]; then cp "$one_subevent_dir/$configuration_file" "$forward_dir/"; fi
    done
    # ... except the Par.file, which comes from the N = 2 folder when there is one: stage 1 shortened
    # the windows there (not in the 1-subevent folder), and the weights must be balanced on the
    # windows that the N >= 2 runs and the N = 1 rerun use (until v3 they were balanced on the
    # step-2 windows, e.g. 155 s instead of 115 s for Chile). ffwd takes the subevent count from
    # Input.model, so the 2-subevent Par.file serves the 1-subevent model.
    two_subevent_par="inv_${event_name}_2sub/Par.file"
    if [ -s "$two_subevent_par" ]; then cp "$two_subevent_par" "$forward_dir/Par.file"; fi
    rm -rf "$forward_dir/data" "$forward_dir/dataloc"
    cp -r "$one_subevent_dir/data" "$one_subevent_dir/dataloc" "$forward_dir/"
    gawk 'NR==1{printf "%.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", $4,$5,$6,$7,$8,$9,$10}' "$best_model" > "$forward_dir/Input.model"
    return 0
}

balance_weights() {
    local p_weight sh_weight regional_weight regional_count is_filler iteration line converged
    p_weight=$(gawk -F= '/^weightP=/{print $2}' "$forward_dir/Par.file" | gawk '{print $1}')
    sh_weight=$(gawk -F= '/^weightSH=/{print $2}' "$forward_dir/Par.file" | gawk '{print $1}')
    regional_weight=$(gawk -F= '/^weightRayl=/{print $2}' "$forward_dir/Par.file" | gawk '{print $1}')
    regional_count=$(grep -c . "$forward_dir/stationsloc.info" 2>/dev/null || echo 0)
    is_filler=$(gawk -v n="$regional_count" 'BEGIN{print (n <= 2) ? 1 : 0}')
    if [ "$is_filler" -eq 1 ]; then echo "calibrate weights: filler regional set ($regional_count stations) - balancing P:SH only"; fi
    cp "$forward_dir/Par.file" "$forward_dir/Par.file.before_balance"
    for iteration in 1 2 3 4; do
        line=$(cd "$forward_dir" && FFWD_COMPACT=1 ./ffwd 2>/dev/null | grep "^bandresid:" | tail -1)
        case "$line" in
            ""|*nan*|*inf*|*NaN*)
                cp "$forward_dir/Par.file.before_balance" "$forward_dir/Par.file"
                echo "calibrate weights: forward run failed or not finite ($line) - weights unchanged"
                return 0 ;;
        esac
        read -r p_weight sh_weight regional_weight converged <<EOF
$(echo "$line" | gawk -v wp="$p_weight" -v wsh="$sh_weight" -v wr="$regional_weight" -v filler="$is_filler" '{
    rp = $2; rsh = $3; rr = $4
    if (filler == 1) { total = rp + rsh; sp = rp/total; ssh = rsh/total; tp = 1.2/2.2; tsh = 1.0/2.2; sr = tp; tr = tp }
    else { total = rp + rsh + rr; sp = rp/total; ssh = rsh/total; sr = rr/total; tp = 1.2/3.2; tsh = 1.0/3.2; tr = 1.0/3.2 }
    fp = sqrt(tp/(sp + 1e-12)); fsh = sqrt(tsh/(ssh + 1e-12)); fr = (filler == 1) ? 1.0 : sqrt(tr/(sr + 1e-12))
    ok = (sp > 0.85*tp && sp < 1.15*tp && ssh > 0.85*tsh && ssh < 1.15*tsh)
    if (filler == 0) ok = ok && (sr > 0.85*tr && sr < 1.15*tr)
    new_wp = wp*fp; new_wsh = wsh*fsh; new_wr = (filler == 1) ? wr : wr*fr
    printf "%.3g %.3g %.3g %d", new_wp/new_wsh, 1.0, (filler == 1) ? wr : new_wr/new_wsh, ok }')
EOF
        if [ "$(gawk -v a="$p_weight" -v b="$sh_weight" -v c="$regional_weight" 'BEGIN{print (a+0>0 && b+0>0 && c+0>0) ? 1 : 0}')" != "1" ]; then
            cp "$forward_dir/Par.file.before_balance" "$forward_dir/Par.file"
            echo "calibrate weights: malformed update ($p_weight/$sh_weight/$regional_weight) - weights unchanged"
            return 0
        fi
        sed -i "s/^weightP=.*/weightP= $p_weight #/; s/^weightSH=.*/weightSH= $sh_weight #/" "$forward_dir/Par.file"
        if [ "$is_filler" -eq 0 ]; then sed -i "s/^weightRayl=.*/weightRayl= $regional_weight/" "$forward_dir/Par.file"; fi
        if [ "$converged" -eq 1 ]; then break; fi
    done
    echo "calibrate weights: P:SH:regional = $p_weight : $sh_weight : $regional_weight (iteration $iteration)"
    for inversion_dir in inv_"${event_name}"_*sub; do
        case "$inversion_dir" in *_1sub) continue ;; esac
        [ -s "$inversion_dir/Par.file" ] || continue
        sed -i "s/^weightP=.*/weightP= $p_weight #/; s/^weightSH=.*/weightSH= $sh_weight #/" "$inversion_dir/Par.file"
        if [ "$is_filler" -eq 0 ]; then sed -i "s/^weightRayl=.*/weightRayl= $regional_weight/" "$inversion_dir/Par.file"; fi
    done
    # explicit: the status of the last test above (1 for a filler set) once became the script's exit
    # status, and step 3 stopped on a calibration that had worked
    return 0
}

build_forward_program || exit 0
balance_weights

# ------------------------------------------------------------------------------------------------
# 3. model-based station screen, then re-balance on the cleaned station set
# ------------------------------------------------------------------------------------------------
if [ "$duration_too_short" -eq 1 ]; then
    export SCREEN_DATA_ONLY=1
    echo "calibrate screen: 1-subevent result unrepresentative (duration << magnitude scaling) - only the amplitude of the records is screened"
fi
if [ "$is_extended" = "1" ]; then export STAGE2_POLARITY_ONLY=1; fi
rm -f excluded_stations_model.txt
( cd "$forward_dir" && ./ffwd > ffwd_screen.log 2>&1 )
[ -s "$forward_dir/sta_residual.dat" ] || { echo "calibrate screen: no sta_residual.dat - skipped"; exit 0; }
"$python_executable" ../programs/prepare_SubMIT/screen_stations_model.py "$event_name" || { echo "calibrate screen: script failed - skipped"; exit 0; }
if [ -s excluded_stations_model.txt ]; then
    for inversion_dir in inv_"${event_name}"_*sub "$forward_dir"; do
        [ -d "$inversion_dir" ] || continue
        cp ../programs/code_SubMIT/prepare_stationinfo.py "$inversion_dir/"
        # a failure here would leave stale station lists next to a cleaned data folder: stop step 3
        ( cd "$inversion_dir" && "$python_executable" prepare_stationinfo.py > prepare_stationinfo.log 2>&1 ) \
            || { echo "calibrate screen: prepare_stationinfo.py failed in $inversion_dir (see its prepare_stationinfo.log)"; exit 1; }
        ( cd "$inversion_dir" \
          && sed -i "s/^num_sta_P=.*/num_sta_P= $(grep -c . stations.info) #Number of all stations (P)/" Par.file \
          && sed -i "s/^num_sta_SH=.*/num_sta_SH= $(grep -c . stationsSH.info) #Number of all stations (SH)/" Par.file \
          && sed -i "s/^num_sta_Rayl=.*/num_sta_Rayl= $(grep -c . stationsloc.info) #Number of all stations (regional)/" Par.file ) \
            || { echo "calibrate screen: station counts not updated in $inversion_dir/Par.file"; exit 1; }
    done
    echo "calibrate screen: station sets changed - re-balancing the weights"
    build_forward_program && balance_weights
fi
exit 0
