NAM=$1
saclst evdp gcarc f $NAM | gawk ' { print "h";print $2;print $3;print "q"; } ' | ../../programs/TauP-2.0/bin/taup_time -mod iasp91 | gawk ' {if ($3=="P"||$3=="p") print $0}' | gawk '{if (NR==1) print 111.195/$5}' >> pslowness.tmp
#cat taup.tmp | gawk -v var=$NAM '{if (NR==1) print var,"P",$1;}' >> ptime.tmp
