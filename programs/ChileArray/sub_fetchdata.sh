tt=$1
mm=$2
nn=$3
link="https://evtdb.csn.uchile.cl/event/"${nn}
wget $link
grep "href=\"\/station\/" $nn | gawk 'BEGIN{FS="[//\"]"} {print $4}' > stalist_$nn
cat stalist_$nn | gawk -v nn="$nn" '{print "wget https://evtdb.csn.uchile.cl/write/"nn"/"$1" -O "$1".zip"}'  | sh
cat stalist_$nn | gawk -v nn="$nn" -v tt="$tt" -v mm="$mm" '{print "mv "$1".zip "tt"_"mm"_"nn"/"}'   | sh
mv $nn stalist_$nn ${tt}_${mm}_${nn}
cd ${tt}_${mm}_${nn}
ls -1 *.zip | gawk '{print "unzip "$1}' | sh
cd ..
