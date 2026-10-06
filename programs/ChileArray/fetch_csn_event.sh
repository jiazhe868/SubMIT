#!/bin/bash
# fetch_csn_event.sh EVENT_ID
# Optional, Chile only: download the strong-motion records of one earthquake from
# the Centro Sismologico Nacional database (https://evtdb.csn.uchile.cl) and
# convert them to SubMIT regional SAC files (velocity, C1.<STA>..[enz]).
# EVENT_ID is the 32-character id in the event page address,
#   https://evtdb.csn.uchile.cl/event/<EVENT_ID>
# Run from the event folder that holds IRIS/ and IRISloc/ ("./submit fetch-csn"
# does this). Writes StrongMotionData/<date-time>_<mag>_<EVENT_ID>/, which
# step1_4_merge_csn.sh ("./submit prepare") merges into the regional data.
# Needs wget, unzip, python (obspy), sac.
id=${1:?usage: fetch_csn_event.sh EVENT_ID}
here=$(cd "$(dirname "$0")" && pwd)
for t in wget unzip sac python; do
    command -v $t > /dev/null || { echo "fetch_csn_event: $t not found"; exit 1; }
done
mkdir -p StrongMotionData && cd StrongMotionData || exit 1
wget -q -O "event_$id.html" "https://evtdb.csn.uchile.cl/event/$id" \
    || { echo "fetch_csn_event: cannot download the event page for $id"; exit 1; }
# event header: "Evento del YYYY-MM-DD HH:MM:SS", Latitud, Longitud, Profundidad, Magnitud
read -r date clock lat lon dep mag <<< "$(gawk '
    /Evento del/ {for (i = 1; i <= NF; i++) if ($i == "del") {d = $(i+1); t = $(i+2)}}
    /Latitud:/ && /strong/ {gsub(/<[^>]*>/, ""); la = $2}
    /Longitud:/ && /strong/ {gsub(/<[^>]*>/, ""); lo = $2}
    /Profundidad:/ && /strong/ {gsub(/<[^>]*>/, ""); de = $2}
    /Magnitud:/ && /strong/ {gsub(/<[^>]*>/, ""); ma = $2}
    END {print d, t, la, lo, de, ma}' "event_$id.html")"
[ -n "$mag" ] || { echo "fetch_csn_event: could not read the event header of $id"; exit 1; }
tt="$date-$(echo "$clock" | tr -d :)"
dir="${tt}_${mag}_${id}"
mkdir -p "$dir"
echo "$lon $lat $dep ${date}T$(echo "$clock" | tr -d :)" > "$dir/loladep.dat"
echo "CSN event $id: $date $clock  lon $lon lat $lat depth $dep km  M $mag"
cp "$here/convert_to_sac.py" "$here/sub_fetchdata.sh" "$here/sub_processdata.sh" .
rm -f "event_$id.html" "$id"       # sub_fetchdata.sh downloads the page again as "$id"
sh sub_fetchdata.sh "$tt" "$mag" "$id" > "$dir/fetch.log" 2>&1
sh sub_processdata.sh "$tt" "$mag" "$id" > "$dir/process.log" 2>&1
n=$(ls "$dir"/C1.*..z 2>/dev/null | wc -l)
echo "fetch_csn_event: $n stations converted in StrongMotionData/$dir"
[ "$n" -gt 0 ]
