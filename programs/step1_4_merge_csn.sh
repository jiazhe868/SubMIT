#!/bin/bash
# Optional: merge CSN (evtdb.csn.uchile.cl) strong-motion SAC data over the Wilber
# local data, matching the SubMIT local format (velocity m/s, delta 1 s, cut -50/400
# relative to origin, o=0, t1=0). CSN wins on station-name collisions.
#
# usage (from IRISloc/):  sh ../programs/step1_4_merge_csn.sh <wilber_event_name>
# CSN input: ../StrongMotionData*/<csn_event>/C1.*..[enz] produced by
# step1_3_process_array.sh + convert_to_sac.py (delta 0.1 s, o=0, b rel. origin).
# If no CSN data exist (e.g. 2026 Venezuela), this script is a clean no-op (exit 0).
ev=$1
dest="${ev}_loc/data"
[ -d "$dest" ] || { echo "step1_4: no $dest (run step1_2 first)"; exit 1; }

csn=$(ls -1 -d ../StrongMotionData*/[12]???-??-*_*/ 2>/dev/null | head -1)
if [ -z "$csn" ] || ! ls "$csn"/C1.*..[enz] > /dev/null 2>&1; then
    echo "step1_4: no CSN strong-motion data found - skipping (Wilber local data only)"
    exit 0
fi

n=0
: > "$dest/.csn_merged.lst"
for f in "$csn"/C1.*..[enz]; do
    base=$(basename "$f")
    cp "$f" "$dest/$base"
    echo "$dest/$base" >> "$dest/.csn_merged.lst"
    n=$((n + 1))
done
# decimate to delta = 1 s (adaptive to source delta), cut to the SubMIT local
# window, t1 = origin time (0)
# process ONLY the files just copied from CSN (Wilber stations may share the
# C1.* naming and must not be touched)
saclst delta f $(cat "$dest/.csn_merged.lst") | gawk '{
    print "cuterr fillz"; print "cut -50 400"; print "r "$1;
    if (sqrt(($2-0.1)^2) < 1e-4)      { print "decimate 5"; print "decimate 2"; }
    else if (sqrt(($2-0.2)^2) < 1e-4) { print "decimate 5"; }
    else if (sqrt(($2-0.5)^2) < 1e-4) { print "decimate 2"; }
    else if (sqrt(($2-1.0)^2) > 1e-4) { print "* WARNING unexpected delta "$2" for "$1; }
    print "ch t1 0"; print "w over"} END{print "q"}' | sac > /dev/null
echo "step1_4: merged $n CSN component file(s) from $csn into $dest (CSN wins on collisions)"
