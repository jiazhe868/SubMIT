#!/bin/bash
# Tighten the n>=2 subevent search bounds using the 1-sub scan's empirical
# source time. The magnitude-scaled duration (5*10^((Mw-6)/2)) overestimates
# compact/deep ruptures (2024 Chile M7.4 deep intraslab: actual ~16 s vs
# scaled 25 s -> cen_max 31 admitted noise-absorber subevents at 26-30 s and
# a non-monotone L-curve). The 1-sub best (cen1, dura1) measures the actual
# total source time; bound the multi-subevent searches by it:
#   cen_max  <- min(cen_max,  cen1 + dura1 + 6)
#   dura_max <- min(dura_max, dura1 + 5)
# Safety: NO-OP when the 1-sub duration railed its cap (>=0.95*dura_max, i.e.
# the scan itself was clamped - giant events) or when the 1-sub output is
# missing. Bounds only ever shrink, never grow.
#
# usage: tighten_bounds_from_1sub.sh <inv_..._1sub dir> <event_name>
d1=$1
ev=$2
best=""
for f in "$d1"/best_model_hybrid.dat "$d1"/best_model_ensemble.dat \
         "$d1"/best_model_exploration.dat; do
    [ -s "$f" ] && best=$f && break
done
[ -n "$best" ] || { echo "tighten_bounds: no 1-sub best model - skipped"; exit 0; }
cen1=$(gawk 'NR==1{printf "%.2f", $4}' "$best")
dur1=$(gawk 'NR==1{printf "%.2f", $7}' "$best")
dmax1=$(gawk '$1=="dura_max"{print $2}' "$d1"/search_par.file)
[ -n "$cen1" ] && [ -n "$dur1" ] && [ -n "$dmax1" ] || { echo "tighten_bounds: parse failure - skipped"; exit 0; }
if [ "$(cat "$d1"/extended.flag 2>/dev/null)" = "1" ]; then
    echo "tighten_bounds: extended/doublet source - skipped"
    exit 0
fi
railed=$(gawk -v d="$dur1" -v m="$dmax1" 'BEGIN{print (d >= 0.95*m) ? 1 : 0}')
if [ "$railed" -eq 1 ]; then
    echo "tighten_bounds: 1-sub duration $dur1 railed dura_max $dmax1 (giant/clamped) - skipped"
    exit 0
fi
# a 1-sub duration far BELOW the magnitude scaling is equally pathological (a
# stuck scan railed at dura_min poisons every n>=2 search: California M7.0,
# dur1=3.0 -> dura_max 8 forbade the real ~9 s subevents). Skip below 30% of
# the Mw-scaled duration.
mw1=$(gawk '{print $4}' "$d1"/mainshock.dat 2>/dev/null)
if [ -n "$mw1" ]; then
    tooshort=$(gawk -v d="$dur1" -v mw="$mw1" 'BEGIN{print (d < 0.3*5*10^(0.5*(mw-6))) ? 1 : 0}')
    if [ "$tooshort" -eq 1 ]; then
        echo "tighten_bounds: 1-sub duration $dur1 << Mw-scaled duration (stuck scan?) - skipped"
        exit 0
    fi
fi
for n in 2 3 4 5 6; do
    sp="inv_${ev}_${n}sub/search_par.file"
    [ -s "$sp" ] || continue
    gawk -v c1="$cen1" -v d1="$dur1" '
        $1=="cen_max"  { nc = c1 + d1 + 6; if (nc < $2) $2 = sprintf("%.1f", nc) }
        $1=="dura_max" { nd = d1 + 5;      if (nd < $2) $2 = sprintf("%.1f", nd) }
        { print }' "$sp" > "$sp.tmp" && mv "$sp.tmp" "$sp"
    echo "tighten_bounds: $sp -> cen_max $(gawk '$1=="cen_max"{print $2}' "$sp") dura_max $(gawk '$1=="dura_max"{print $2}' "$sp")"
    # body windows follow the bounds (user rule 2026-08-12): fitting long
    # codas past the last reachable subevent buys nothing for compact events
    # (California M7.0: auto nd_timeP was 86 s vs the hand-tuned 50). Window
    # = cen_max + dura_max + 15 s margin (depth phases + alignment slack);
    # shrink-only, same skip guards as the bounds (this code is unreached
    # for giant/doublet events).
    pf="inv_${ev}_${n}sub/Par.file"
    if [ -s "$pf" ]; then
        cmx=$(gawk '$1=="cen_max"{print $2}' "$sp"); dmx=$(gawk '$1=="dura_max"{print $2}' "$sp")
        # window must hold the last subevent's DEPTH PHASES, not just its
        # direct pulse (user rule 2026-08-13): sS-S ~ 2z/4 (~60 s at 120 km),
        # pP/sP similar order. Shallow events keep a 5-s phase margin.
        edp=$(gawk '{print $3}' "$d1"/mainshock.dat 2>/dev/null)
        tdp=$(gawk -v z="${edp:-0}" 'BEGIN{printf "%d", (z>60) ? 2*z/4.0 : 5}')
        ndnew=$(gawk -v c="$cmx" -v d="$dmx" -v t="$tdp" 'BEGIN{v=c+d+10+t; if (v>300) v=300; printf "%d", v}')
        for key in nd_timeP nd_timeSH; do
            sed -i "s/^${key}=.*/${key}= $ndnew #ending time of inverse time window (cen_max+dura_max+10+depth-phase 2z\/4)/" "$pf"
        done
        echo "tighten_bounds: $pf windows -> $ndnew s (depth-phase term $tdp s)"
    fi
done
