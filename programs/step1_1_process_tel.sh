ls -1 *.tar | gawk '{print "tar -xvf "$1}' | sh
mkdir tar
mv *.tar tar
# normalize Wilber download suffixes: a tel tar may extract to
# <event>_tel-<region>/; the IRISloc pairing (step2 expects
# ../IRISloc/<event>_loc) requires the bare <event> dir name
ls -1 -d [12]???-??-??-*_tel* 2>/dev/null | while read d; do mv "$d" "${d%%_tel*}"; done
ls -1 -d [12]???-??-??-* | gawk '{print "cp ../programs/MyExtractSeed.sh "$1;print "cd "$1;print "sh MyExtractSeed.sh";print "cd .."}' | sh
