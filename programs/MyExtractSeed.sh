#!/bin/sh

#########################################################
# * copyright: 2014-, Earth Science School of USTC
# * File name: LeadData.cmd
# * Description: Prepare seismic data for earthquakes
# * Author:Jiazhe
# * Date:06/22/2014
# * History:07/05/2013 v0.1
# *	    11/20/2013 v0.2
# *	    06/22/2014 v1.0
#########################################################


### Instruction: sh MyExtractSeed.sh ***.seed ###
### Require SAC, TauP and rdseed being installed ###
### Arrival times: programs/select/add_arrival.sh (TauP iasp91, bundled) ###


###############    extract the seed, add info and change name     #####################
#seedname=$(ls -1 *.seed | gawk '{print $1}')
ls -1 *.SAC > 1SAC1
cat 1SAC1 | gawk '{split($1,m,".");print "mv "$1, m[1]"."m[2]"."m[3]"."m[4]}' | sh

###############    make directories    #################################
if [ -d RESP ]; then
rm RESP/*
echo "found RESP/"
else
mkdir RESP
fi
if [ -d data ]; then
echo "found Vel/";
else
mkdir data;
fi
if [ -d Vel ]; then
echo "found Vel/";
else
mkdir Vel;
fi

################     remove the instrumental response     ###############################
find . -name "SACPZ*" > 1PZs
cat 1PZs | gawk '{split($1,mm,"[_/.]"); if (mm[7]=="--") mm[7]="";print "mv "$1" RESP/"mm[5]"."mm[6]"."mm[7]"."mm[8]}'  | sh
sed -i '6,8s/Z//g' RESP/*
ls -1 *.BH? > BH1
saclst o b nzyear nzjday nzhour nzmin nzsec nzmsec gcarc f *.BH? > headerlist.dat
gawk '
function is_leap(year) { return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0 }
{
    file = $1; o_time = $2; b_time = $3; nzyear = $4; nzjday = $5
    seconds = $6 * 3600 + $7 * 60 + $8 + $9 / 1000 + o_time
    # crossing midnight moves the reference day, across New Year and leap years too (v2 assumed
    # 365-day years)
    while (seconds < 0) { seconds += 86400; nzjday -= 1; if (nzjday < 1) { nzyear -= 1; nzjday = is_leap(nzyear) ? 366 : 365 } }
    while (seconds >= 86400) { seconds -= 86400; nzjday += 1; if (nzjday > (is_leap(nzyear) ? 366 : 365)) { nzjday = 1; nzyear += 1 } }
    hour = int(seconds / 3600); seconds -= hour * 3600
    minute = int(seconds / 60); seconds -= minute * 60
    print "r " file
    print "ch o 0 b " b_time - o_time
    printf "ch nzyear %d nzjday %d nzhour %02d nzmin %02d nzsec %02d nzmsec %03d\n", nzyear, nzjday, hour, minute, int(seconds), int((seconds - int(seconds)) * 1000)
    # records prepared outside Wilber may lack gcarc/dist/az (README 7.3 promises they work): let SAC
    # compute them; Wilber records keep their own values, so the examples are reproduced exactly
    if ($10 == -12345) print "ch lcalda true"
    print "wh"
}
END { print "q" }' headerlist.dat | sac
cat BH1 |gawk '{
	print "r",$1;
	print "transfer from polezero subtype RESP/"$1,"to vel freq 0.001 0.005 10 20";
        print "rmean";
        print "rtr";
        print "taper";
	print "w Vel/"$1;
}
END{print "q"}' | sac
rm *.BH?
cd Vel/

#########     Remove stations that do not have fual 3 components  (In case)     ###########################
#ls -1 *.BH* | gawk '
#BEGIN{
#        FS=".";
#        num=0;
#        a="AA";b="AAA";c="AA";
#        }
#{
#        if ($1==a&&$2==b&&$3==c)
#                {num+=1;}
#        else if(num<=2&&(a!="AA"||b!="AAA"))
#                {       print "rm "a"."b"."c".*";
#                        num=1;a=$1;b=$2;c=$3}
#                else {  num=1;a=$1;b=$2;c=$3}
#}
#END{
#        if (num<=2) print "rm "a"."b"."c".*";
#}' | sh

##################   Add P or S arrival times   ###################
sh ../../../programs/select/add_arrival.sh P *.BHZ
sh ../../../programs/select/add_arrival.sh S *.BH[EN12]

########   change those stations which CMPAZs have minor error   ##########
ls -1 *BH1 | gawk 'BEGIN{FS="."} {print $1"."$2"."$3".BH"}' > BH1list
ls -1 *BH2 | gawk 'BEGIN{FS="."} {print $1"."$2"."$3".BH"}' > BH2list
comm -12 BH1list BH2list > BH12list.dat
ls -1 *BHE | gawk 'BEGIN{FS="."} {print $1"."$2"."$3".BH"}' > BHElist
ls -1 *BHN | gawk 'BEGIN{FS="."} {print $1"."$2"."$3".BH"}' > BHNlist
comm -12 BHElist BHNlist > BHENlist.dat
cat BH12list.dat | gawk '{print "saclst cmpaz f "$1"1"}' | sh > 1az.dat
cat BH12list.dat | gawk '{print "saclst cmpaz f "$1"2"}' | sh > 2az.dat
cat BHENlist.dat | gawk '{print "saclst cmpaz f "$1"E"}' | sh > eaz.dat
cat BHENlist.dat | gawk '{print "saclst cmpaz f "$1"N"}' | sh > naz.dat
paste 1az.dat 2az.dat > 12az.dat
paste eaz.dat naz.dat > enaz.dat
cat 12az.dat | gawk '{
        diff1=($2+360-$4)%360;
        pota4=($2+270)%360;
        diff2=($4+360-$2)%360;
        potb4=($2+90)%360;
        if (diff1<diff2) {
                if (diff1!=90) print $3,pota4;
        }
        else {
                if (diff2!=90) print $3,potb4;
        }
}' | awk '{
        print "r "$1;
        print "ch cmpaz "$2;
        print "wh";
}
END{ print "q";}' | sac
cat enaz.dat | gawk '{
        diff1=($2+360-$4)%360;
        pota4=($2+270)%360;
        diff2=($4+360-$2)%360;
        potb4=($2+90)%360;
        if (diff1<diff2) {
                if (diff1!=90) print $3,pota4;
        }
        else {
                if (diff2!=90) print $3,potb4;
        }
}' | awk '{
        print "r "$1;
        print "ch cmpaz "$2;
        print "wh";
}
END{ print "q";}' | sac
rm *az.dat
##############    Select good data using SNR   #############################
#perl ../etc_commands/select/selectData8snr.pl
#rm *BH*
#mv ./good/*BH* .
#rm -rf good

##############    Rot to gcp and cut (For BH[EN] and BH[12], respectively)   ###############
# SAC's rotate needs both horizontals with the same length: for a pair whose records differ in length
# (a short gap in one component) it fails with ERROR 2010, and the next write stored the unrotated east
# and north records as .r and .t (PITFALLS A1: 26 of the 56 SH records of the published Venezuela
# model). Cut such pairs to their common window first; pairs of equal length are untouched.
for pair in EN 12; do
    first=$(echo $pair | cut -c1); second=$(echo $pair | cut -c2)
    while read stem; do
        # only a length mismatch: Wilber pairs differ in b by microseconds, which rotate tolerates
        window=$(saclst b e npts f ${stem}${first} ${stem}${second} | gawk 'NR==1{b1=$2; e1=$3; n1=$4} NR==2{b2=$2; e2=$3; n2=$4} END{if (n1 != n2) print (b1 > b2 ? b1 : b2), (e1 < e2 ? e1 : e2)}')
        [ -n "$window" ] && printf "cut %s\nr %s %s\nw over\ncut off\nq\n" "$window" ${stem}${first} ${stem}${second} | sac > /dev/null
    done < BH${pair}list.dat
done
mkdir rtr/
cat BHENlist.dat | gawk '{print "saclst t1 f "$1"Z"}' | sh | gawk '{a=$1;print "r "$1 ;print "rtr" ; gsub("BHZ","z",a);print "w rtr/"a ;a=$1; b=$1; gsub("BHZ", "BHE",a); gsub("BHZ","BHN", b);print "r",a,b; print "rtr" ; print "rot to gcp"; gsub("BHE","r",a); gsub("BHN","t",b); print "w rtr/"a, "rtr/"b;} END {print "quit";}' | tee | sac
cat BH12list.dat | gawk '{print "saclst t1 f "$1"Z"}' | sh | gawk '{a=$1;print "r "$1 ;print "rtr" ; gsub("BHZ","z",a);print "w rtr/"a ;a=$1; b=$1; gsub("BHZ", "BH1",a); gsub("BHZ","BH2", b);print "r",a,b; print "rtr" ; print "rot to gcp"; gsub("BH1","r",a); gsub("BH2","t",b); print "w rtr/"a, "rtr/"b;} END {print "quit";}' | tee | sac

##############    Quality control     ########################
cd rtr/
script_path="../../../../programs/select/selectSNR.pl"
if [ -f "$script_path" ]; then
    echo "Script $script_path found. Running the script..."
    # Run the script (you can pass any required arguments if needed)
    perl "$script_path"
else
    echo "Script $script_path does not exist. Exiting..."
    # Exit or break depending on the context
    exit 1
fi
rm *.[rtz]
mv ./good_data/*.[rtz] .
rm -rf good_data
mkdir others
saclst depmax f *.[rtz] | gawk '{if (sqrt($2^2)>1e-3) print "mv "$1" others"}' | sh > /dev/null

for base in $(ls *.r *.t *.z 2>/dev/null | sed -e 's/\.[rtz]//' | sort | uniq); do
    # Check if all three components (r, t, z) exist
    if [ -f "${base}.r" ]  &&  [ -f "${base}.t" ] && [ -f "${base}.z" ]; then
        echo "Complete set found for $base"

        # Perform further quality control checks (e.g., SNR, amplitude, etc.)
        # Insert additional QC checks here
        # ...

    else
        echo "Incomplete set for $base, moving to others"
        # Move all existing components to the others folder
        [ -f "${base}.r" ] && mv "${base}.r" others/
        [ -f "${base}.t" ] && mv "${base}.t" others/
        [ -f "${base}.z" ] && mv "${base}.z" others/
    fi
done

ls -1 *.00.[rtz] | gawk ' BEGIN {FS="."} {print "mv "$1"."$2".[1-9]0."$4" others"}' | sh 2>/dev/null
ls -1 *.10.[rtz] | gawk ' BEGIN {FS="."} {print "mv "$1"."$2".[2-9]0."$4" others"}' | sh 2>/dev/null

sh ../../../../programs/select/decimate_to_1sps.sh rtz

##############    Remove RESP files    #############################
cp *.[rtz] ../../data
cd ../../
rm -rf RESP
