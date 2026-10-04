#!/bin/bash

# select stations to keep in fwd_*_1sub

while read line; do
    event_name=$(echo $line | awk '{print $1}')
    mw=$(echo $line | awk '{print $4}')
    nsub=1
    while [ $nsub -le 5 ]; do
        cp -r "fwd_${event_name}/stations.info" "fwd_${event_name}/stationsSH.info" "fwd_${event_name}/stationsloc.info" "fwd_${event_name}/Par.file" "inv_${event_name}_${nsub}sub"
        nsub=$((nsub + 1))
    done
done < event_list.dat
