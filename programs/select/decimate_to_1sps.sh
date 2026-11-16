#!/bin/sh
# decimate_to_1sps.sh rtz|enz - bring every record of the current folder to 1 sample/s (step 1).
#
# SAC's anti-alias "decimate" (factors 2-7 only) chained per input rate: 500, 200, 100, 80, 50, 40,
# 25 and 20 samples/s. A record at any other rate would be read with the wrong time step by the
# inversion, so its whole station moves to others/ with a warning. Shared by MyExtractSeed.sh
# (teleseismic r, t, z) and MyExtractSeed_loc.sh (regional e, n, z); the chains are those of v2,
# so the processed records of the examples are reproduced exactly.
components=${1:?usage: decimate_to_1sps.sh rtz|enz}
saclst delta f *.*.*.[$components] 2>/dev/null | gawk '
function chain(file, factors,    n, i, f) {
    n = split(factors, f, " ")
    print "r", file
    for (i = 1; i <= n; i++) print "decimate", f[i]
    print "w over"
}
{
    if ((($2 - 0.002) ^ 2) ^ 0.5 < 1e-4)  chain($1, "2 5 5 2 5")      # 500 sps
    if ((($2 - 0.005) ^ 2) ^ 0.5 < 1e-4)  chain($1, "2 2 5 2 5")      # 200 sps
    if ((($2 - 0.01) ^ 2) ^ 0.5 < 1e-4)   chain($1, "2 5 2 5")        # 100 sps
    if ((($2 - 0.0125) ^ 2) ^ 0.5 < 1e-5) chain($1, "2 2 4 5")        # 80 sps
    if ((($2 - 0.02) ^ 2) ^ 0.5 < 1e-4)   chain($1, "5 5 2")          # 50 sps
    if ((($2 - 0.025) ^ 2) ^ 0.5 < 1e-4)  chain($1, "2 2 2 5")        # 40 sps
    if ((($2 - 0.04) ^ 2) ^ 0.5 < 1e-4)   chain($1, "5 5")            # 25 sps
    if ((($2 - 0.05) ^ 2) ^ 0.5 < 1e-4)   chain($1, "2 2 5")          # 20 sps
} END { print "quit" }' | sac
mkdir -p others
saclst delta f *.*.*.[$components] 2>/dev/null | gawk '{if (($2 - 1) ^ 2 > 1e-6) print $1}' | sed 's/\.[a-z]$//' | sort -u | while read base; do
    echo "WARNING: $base is not at 1 sample/s after decimation - moved to others/"
    mv "$base".[$components] others/
done
