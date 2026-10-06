#!/bin/sh
# add_arrival.sh P|S file ...
# Writes the first iasp91 P (or S) arrival time, computed by TauP from the SAC
# headers evdp (km) and gcarc (deg), into header t1 of every file. Arrival
# names P/p (or S/s) are accepted, so direct up-going phases count at short
# distances. Needs taup_time (bundled in programs/TauP-2.0/bin), saclst, sac.
ph=$1; shift
case $ph in P) lc=p ;; S) lc=s ;; *) echo "usage: add_arrival.sh P|S file ..." >&2; exit 1 ;; esac
tmp=arrivals_$ph.$$.tmp
: > "$tmp"
for f in "$@"; do
    [ -f "$f" ] || continue
    saclst evdp gcarc f "$f" | gawk '{print "h"; print $2; print $3; print "q"}' | taup_time \
        | gawk -v a="$ph" -v b="$lc" '$3 == a || $3 == b {print $4}' | sort -n | head -1 \
        | gawk -v f="$f" '{print f, $1}' >> "$tmp"
done
[ -s "$tmp" ] && gawk '{print "r", $1; print "ch t1", $2; print "wh"} END {print "quit"}' "$tmp" | sac > /dev/null
rm -f "$tmp"
