ls -1 *.tar | gawk '{print "tar -xvf "$1}' | sh
mkdir tar
mv *.tar tar
ls -1 -d [12]???-??-??-* | gawk '{print "cp ../programs/MyExtractSeed_loc.sh "$1;print "cd "$1;print "sh MyExtractSeed_loc.sh";print "cd .."}' | sh
