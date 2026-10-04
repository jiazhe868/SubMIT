#!/bin/bash

saclst gcarc az f data/*.z | grep -v '^data/vel_' > sta_dist_az

# Create a distance list from the station distance and azimuth data
cat sta_dist_az | gawk '{print "saclst dist f "$1}' | sh | gawk '{printf ("%d\n",$2)}' > dist_list

# Run the slowness and travel-time script
sh doslownesstaup.sh

# Combine slowness data with distance data for station info
paste allslowness dist_list | gawk '{print $1,$2,$3,$4,$5,$6}' > stations.info

# Combine slowness data with distance data, prefixing with "vel_" for velocity info
paste allslowness dist_list | gawk '{print $1,$2,$3,$4,$5,"vel_"$6}' > vel_stations.info

# Modify paths in the velocity stations info
sed -i 's/data\//data\/vel_/g' vel_stations.info

# Append the velocity stations info to the station info
cat vel_stations.info >> stations.info

# Process .t files similarly
saclst gcarc az f data/*.t > sta_dist_az
cat sta_dist_az | gawk '{print "saclst dist f "$1}' | sh | gawk '{printf ("%d\n",$2)}' > dist_list
sh doslownesstaup.sh
paste allslowness dist_list | gawk '{print $1,$2,$3,$4,$5,$6}' > stationsSH.info

# Process location data. Station cutoff follows the local GF coverage:
# the GF grid ends at nd_localdist (distdep.dat col 2) and carries a +200 km
# buffer for subevent search offsets, so usable stations are those within
# nd_localdist-200 (Calama: 650-200 = the old hardcoded 450)
ndloc=$(gawk '{print $2}' distdep.dat 2>/dev/null); ndloc=${ndloc:-650}
cutoff=$((ndloc - 200))
saclst stlo stla dist f dataloc/*.z | sort -nk4 | gawk -v c="$cutoff" '{if ($4<c) print $0}' > loc_data_sorted
cat loc_data_sorted | gawk 'BEGIN{FS="[//.]"} {print $2"."$3" 1 1 1"}' > loc_processed
paste loc_data_sorted loc_processed | gawk '{print $1,$2,$3,$5,"1 1 1 "}' > stationsloc.info

# Clean up the file extensions in the stationsloc.info file
sed -i 's/.z/./g' stationsloc.info

# Clean up temporary files
rm gcarc_az_data loc_data_sorted loc_processed dist_list sta_dist_az allslowness

