#!/bin/bash
dp=$1
dd=$2
mind=$3
md=$4
# Initialize an empty string for distlst
distlst=""
if [ "$dd" -gt 0 ] 2>/dev/null && [ "$md" -gt 0 ] 2>/dev/null; then
    # Populate the string with values from dd to md, incremented by dd
    i=$mind
    while [ $i -le $md ]; do
        distlst="$distlst $i"
        i=$((i + dd))
    done

    distlst=$(echo $distlst)

else
    echo "Both dd and md must be positive integers."
fi

DELTAT=$(gawk -F= '/^DELTAT=/{print $2}' ../../../programs/submit.conf 2>/dev/null); DELTAT=${DELTAT:-1}
gawk -v dist="$distlst" -v dep="$dp" -v dt="$DELTAT" 'BEGIN  { print " perl ../../../../programs/fk3.2/fk.pl -Mvmodel/"dep"/f -S2 -N1024/"dt"/1 "dist }' > vmodel_${dp}/fk.cmd
gawk -v dist="$distlst" -v dep="$dp" -v dt="$DELTAT" 'BEGIN  { print " perl ../../../../programs/fk3.2/fk.pl -Mvmodel/"dep"/f -S0 -N1024/"dt"/1 "dist }' >> vmodel_${dp}/fk.cmd

