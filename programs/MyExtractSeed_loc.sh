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
### Require addp.cmd, adds.cmd, pickp.cmd, picks.cmd in ~/bin/; or you can change the path ###


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
# ---- channel auto-detect: strong-motion HN? (e.g. Calama) or broadband BH?
# (e.g. Kamchatka, where no usable strong motion exists). Runs AFTER the SEED
# extraction above, when the SAC components exist. The rest of the processing
# is identical up to the channel code, so hand off to a BH-substituted copy.
if ! ls *.HN? > /dev/null 2>&1 && ls *.BH? > /dev/null 2>&1; then
    # extract the processing tail (from the channel-listing line to EOF) and
    # substitute the channel code
    awk '/^ls -1 \*\.HN\? > HN1$/,0' "$0" | sed 's/HN/BH/g' > rest_BH.sh
    echo "MyExtractSeed_loc: no HN? channels, BH? found - broadband mode"
    sh rest_BH.sh
    exit $?
fi
ls -1 *.HN? > HN1
saclst o b nzyear nzjday nzhour nzmin nzsec nzmsec b f *.HN? > headerlist.dat
cat headerlist.dat | awk '
function to_seconds(hour, min, sec, msec) {
    return hour * 3600 + min * 60 + sec + msec / 1000;
}

function from_seconds(total_seconds) {
    hour = int(total_seconds / 3600);
    total_seconds -= hour * 3600;
    min = int(total_seconds / 60);
    sec = int(total_seconds % 60);
    msec = int((total_seconds - int(total_seconds)) * 1000);
    return sprintf("%02d %02d %02d %03d", hour, min, sec, msec);
}

{
    # Read current values
    o_time = $2;
    b_time = $3;
    nzyear = $4;
    nzjday = $5;
    nzhour = $6;
    nzmin = $7;
    nzsec = $8;
    nzmsec = $9;
    file = $1;

    # Set o to 0 and adjust b
    new_b_time = b_time - o_time;

    # Calculate new time based on o_time
    time_in_seconds = to_seconds(nzhour, nzmin, nzsec, nzmsec) + o_time;

    # Handle possible year overflow (ignoring leap year complexity for now)
    if (nzjday > 365) {
        nzjday = 1;
        nzyear += 1;
    }

    # Convert back to hours, minutes, seconds, and milliseconds
    new_time = from_seconds(time_in_seconds);
    split(new_time, time_parts);

    new_hour = time_parts[1];
    new_min = time_parts[2];
    new_sec = time_parts[3];
    new_msec = time_parts[4];

    # Print SAC commands
    print "r " file;
    print "ch o 0 b " new_b_time;
    print "ch nzyear " nzyear " nzjday " nzjday " nzhour " new_hour " nzmin " new_min " nzsec " new_sec " nzmsec " new_msec;
    print "wh";
}
END { print "q" }' | sac

cat HN1 |gawk '{
	print "r",$1;
	print "transfer from polezero subtype RESP/"$1,"to vel freq 0.001 0.005 10 20";
        print "rmean";
        print "rtr";
	print "hp c 0.02 n 4 p 2"
        print "taper";
	print "w Vel/"$1;
}
END{print "q"}' | sac
rm *.HN?
cd Vel/

##################   Add P or S arrival times   ###################
sh ~/bin/addp.cmd *.HNZ
sh ~/bin/adds.cmd *.HN[EN12]

########   change those stations which CMPAZs have minor error   ##########
ls -1 *HN1 | gawk 'BEGIN{FS="."} {print $1"."$2"."$3".HN"}' > HN1list
ls -1 *HN2 | gawk 'BEGIN{FS="."} {print $1"."$2"."$3".HN"}' > HN2list
comm -12 HN1list HN2list > HN12list.dat
ls -1 *HNE | gawk 'BEGIN{FS="."} {print $1"."$2"."$3".HN"}' > HNElist
ls -1 *HNN | gawk 'BEGIN{FS="."} {print $1"."$2"."$3".HN"}' > HNNlist
comm -12 HNElist HNNlist > HNENlist.dat
cat HN12list.dat | gawk '{print "saclst cmpaz f "$1"1"}' | sh > 1az.dat
cat HN12list.dat | gawk '{print "saclst cmpaz f "$1"2"}' | sh > 2az.dat
cat HNENlist.dat | gawk '{print "saclst cmpaz f "$1"E"}' | sh > eaz.dat
cat HNENlist.dat | gawk '{print "saclst cmpaz f "$1"N"}' | sh > naz.dat
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
#rm *HN*
#mv ./good/*HN* .
#rm -rf good

##############    cut (For HN[EN] and HN[12], respectively)   ###############
mkdir enz/
cat HNENlist.dat | gawk '{print "saclst t1 f "$1"Z"}' | sh | gawk '{a=$1;print "r "$1 ; gsub("HNZ","z",a);print "w enz/"a ;a=$1; b=$1; gsub("HNZ", "HNE",a); gsub("HNZ","HNN", b);print "r",a,b; gsub("HNE","e",a); gsub("HNN","n",b); print "w enz/"a, "enz/"b;} END {print "quit";}' | tee | sac
cat HN12list.dat | gawk '{print "saclst t1 f "$1"Z"}' | sh | gawk '{a=$1;print "r "$1 ; gsub("HNZ","z",a);print "w enz/"a ;a=$1; b=$1; gsub("HNZ", "HN1",a); gsub("HNZ","HN2", b);print "r",a,b;  print "rot to 0"; gsub("HN1","e",a); gsub("HN2","n",b); print "w enz/"b, "enz/"a;} END {print "quit";}' | tee | sac

##############    Quality control     ########################
cd enz/
script_path="../../../../programs/select/selectSNRloc.pl"
if [ -f "$script_path" ]; then
    echo "Script $script_path found. Running the script..."
    # Run the script (you can pass any required arguments if needed)
    perl "$script_path"
else
    echo "Script $script_path does not exist. Exiting..."
    # Exit or break depending on the context
    exit 1
fi
rm *.[enz]
mv ./good_data/*.[enz] .
rm -rf good_data
mkdir others
#saclst depmax f *.[enz] | gawk '{if (sqrt($2^2)>1e-3) print "mv "$1" others"}' | sh > /dev/null

for base in $(ls *.e *.n *.z 2>/dev/null | sed -e 's/\.[enz]//' | sort | uniq); do
    if [ -f "${base}.e" ]  &&  [ -f "${base}.n" ] && [ -f "${base}.z" ]; then
        echo "Complete set found for $base"
    else
        echo "Incomplete set for $base, moving to others"
        [ -f "${base}.e" ] && mv "${base}.e" others/
        [ -f "${base}.n" ] && mv "${base}.n" others/
        [ -f "${base}.z" ] && mv "${base}.z" others/
    fi
done

#ls -1 *.[enz] | gawk ' BEGIN {FS="."} {print "mv "$1"."$2".[1-9]0."$4" others"}' | sh 2>/dev/null
#ls -1 *.[enz] | gawk ' BEGIN {FS="."} {print "mv "$1"."$2".[2-9]0."$4" others"}' | sh 2>/dev/null

saclst delta f *.*.*.[enz] | gawk '{
if(sqrt(($2-0.025)^2)<1e-4){
                print "r",$1;
                print "decimate 2";
		print "decimate 2";
                print "decimate 2";
                print "decimate 5";
                print "w over"
        }
if(sqrt(($2-0.01)^2)<1e-4){
                print "r",$1;
		print "decimate 2";
                print "decimate 5";
                print "decimate 2";
                print "decimate 5";
                print "w over";
        }
if(sqrt(($2-0.005)^2)<1e-4){
                print "r",$1;
                print"decimate 2";
		print "decimate 2";
                print "decimate 5";
                print "decimate 2";
                print "decimate 5";
                print "w over";
        }
if(sqrt(($2-0.002)^2)<1e-4){
                print "r",$1;
		print "decimate 2";
                print"decimate 5";
                print "decimate 5";
                print "decimate 2";
                print "decimate 5";
                print "w over";
        }
if(sqrt(($2-0.05)^2)<1e-4){
                print "r",$1;
		print "decimate 2";
                print "decimate 2";
                print "decimate 5";
                print "w over";
        }
} END{print "quit"}' | sac

##############    Remove RESP files    #############################
cp *.[enz] ../../data
cd ../../data/
ls -1 *.[enz] | gawk '{print "cuterr fillz"; print "cut -50 600";print "r "$1;print "ch t1 0"; print "w over"} END{print "q"}' | sac
ls -1 *.[enz] | gawk 'BEGIN{FS="."} {if ($1=="C") print "mv "$0" C1."$2"."$3"."$4}' | sh 
cd ..
rm -rf RESP
