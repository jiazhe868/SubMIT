saclst gcarc az f data/*.z > aaa
nsta=$(cat aaa | wc -l | gawk '{print $1/2}')
cat aaa | gawk -v nn="$nsta" '{if (NR<=nn) print $0}' > sta_dist_az
cat sta_dist_az  | gawk '{print "saclst dist f "$1}'| sh | gawk '{printf ("%d\n",$2)}' > distlst
sh doslownesstaup.sh
paste allslowness distlst | gawk '{print $1,$2,$3,$4,$5,$6}' > stations.info
paste allslowness distlst | gawk '{print $1,$2,$3,$4,$5,"vel_"$6}' > bbb2
sed -i 's/data\//data\/vel_/g' bbb2
cat bbb2 >> stations.info

saclst gcarc az f data/*.t > sta_dist_az
cat sta_dist_az  | gawk '{print "saclst dist f "$1}'| sh | gawk '{printf ("%d\n",$2)}' > distlst
sh doslownesstaup.sh
paste allslowness distlst | gawk '{print $1,$2,$3,$4,$5,$6}' > stationsSH.info

saclst stlo stla dist f dataloc/*.z | sort -nk4 > bbb1 
cat bbb1 | gawk 'BEGIN{FS="[//.]"} {print $2"."$3" 1 1 1"}' > bbb2 
paste bbb1  bbb2   | gawk '{print $1,$2,$3,$5,"1 1 1 "}' > stationsloc.info
sed -i 's/.z/./g' stationsloc.info
rm aaa bbb1 bbb2 distlst sta_dist_az allslowness
