[ -f pslowness.tmp ] && rm pslowness.tmp
[ -f sslowness.tmp ] && rm sslowness.tmp
[ -f stationlist1 ] && rm stationlist1 > /dev/null
cat sta_dist_az| gawk '{print $1,$2,$3}' > stationlist1
cat sta_dist_az | gawk '{print "sh slownesstaup.sh "$1}'| sh
cat sta_dist_az | gawk '{print "sh slownesstaupS.sh "$1}'| sh
paste stationlist1 pslowness.tmp sslowness.tmp > allslowness
[ -f pslowness.tmp ] && rm pslowness.tmp
[ -f sslowness.tmp ] && rm sslowness.tmp
[ -f stationlist1 ] && rm stationlist1 > /dev/null
