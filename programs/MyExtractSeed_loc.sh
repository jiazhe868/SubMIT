#!/bin/sh
# MyExtractSeed_loc.sh - step 1 for the regional records of one event (run inside IRISloc/<event>_loc/,
# called by ./submit prepare through step1_2_process_loc.sh).
#
# From the raw Wilber SAC files (counts) and their SACPZ responses to IRISloc/<event>_loc/data/:
# ground velocity (m/s) high-passed at 0.02 Hz, origin time as time reference (o = 0), P/S arrivals
# in t1, horizontals rotated to east/north, signal-to-noise selection, 1 sample/s, cut from 50 s
# before to 600 s after the origin, one instrument per station.
#
# Channel families: every strong-motion (HN?) and broadband (HH?, BH?) family in the download is
# processed, and select/choose_regional_instruments.py then keeps one record per station (broadband
# unless clipped or inconsistent with a co-located strong-motion record). Until v3 only HN? was used
# whenever any HN? channel existed, broadband only when none did, and HH? never (cohort report
# 2026-10-23, item 3; the 2026 Panama clean-environment test kept 4 of ~20 stations).
# Rejected records stay in Vel/enz_<family>/bad_data (low SNR), Vel/enz_<family>/others and
# Vel/not_selected.
#
# Author: Zhe Jia (2014-, v1 for SEED volumes, rewritten for Wilber SAC downloads).

#---------------------------------------------------------------------------------------------------
# names, responses, channel families
#---------------------------------------------------------------------------------------------------
ls -1 *.SAC > 1SAC1
gawk '{split($1, m, "."); print "mv "$1, m[1]"."m[2]"."m[3]"."m[4]}' 1SAC1 | sh
mkdir -p RESP data Vel
rm -f RESP/*
find . -name "SACPZ*" > 1PZs
gawk '{split($1, mm, "[_/.]"); if (mm[7] == "--") mm[7] = ""; print "mv "$1" RESP/"mm[5]"."mm[6]"."mm[7]"."mm[8]}' 1PZs | sh
# SAC's polezero reader does not parse the trailing Z of the ISO times on lines 6-8 (CREATED, START, END)
sed -i '6,8s/Z//g' RESP/*
families=""
for family in HN HH BH; do
    ls ./*.${family}? > /dev/null 2>&1 && families="$families $family"
done
[ -n "$families" ] || { echo "MyExtractSeed_loc: no HN?, HH? or BH? records here"; exit 1; }
echo "MyExtractSeed_loc: channel families:$families"
python ../../programs/select/choose_regional_instruments.py clip-check

#---------------------------------------------------------------------------------------------------
# origin time as reference time (o = 0), distances when the header lacks them
#---------------------------------------------------------------------------------------------------
saclst o b nzyear nzjday nzhour nzmin nzsec nzmsec gcarc f ./*.[HB][HN]? > headerlist.dat
gawk '
function is_leap(year) { return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0 }
{
    file = $1; o_time = $2; b_time = $3; nzyear = $4; nzjday = $5
    seconds = $6 * 3600 + $7 * 60 + $8 + $9 / 1000 + o_time
    # crossing midnight (an event a few minutes after 00:00 UTC, records start before it) moves the
    # reference day, including across New Year and leap years
    while (seconds < 0) { seconds += 86400; nzjday -= 1; if (nzjday < 1) { nzyear -= 1; nzjday = is_leap(nzyear) ? 366 : 365 } }
    while (seconds >= 86400) { seconds -= 86400; nzjday += 1; if (nzjday > (is_leap(nzyear) ? 366 : 365)) { nzjday = 1; nzyear += 1 } }
    hour = int(seconds / 3600); seconds -= hour * 3600
    minute = int(seconds / 60); seconds -= minute * 60
    print "r " file
    print "ch o 0 b " b_time - o_time
    printf "ch nzyear %d nzjday %d nzhour %02d nzmin %02d nzsec %02d nzmsec %03d\n", nzyear, nzjday, hour, minute, int(seconds), int((seconds - int(seconds)) * 1000)
    # records prepared outside Wilber may lack gcarc/dist/az: let SAC compute them
    if ($10 == -12345) print "ch lcalda true"
    print "wh"
}
END { print "q" }' headerlist.dat | sac

#---------------------------------------------------------------------------------------------------
# one channel family: response, arrivals, orientation, rotation, SNR, complete sets, 1 sample/s
#---------------------------------------------------------------------------------------------------
process_family() {
    family=$1
    work=Vel/enz_$family
    ls -1 ./*.${family}? | sed 's|^\./||' | gawk '{
        print "r", $1
        print "transfer from polezero subtype RESP/"$1, "to vel freq 0.001 0.005 10 20"
        print "rmean"
        print "rtr"
        print "hp c 0.02 n 4 p 2"
        print "taper"
        print "w Vel/"$1
    } END { print "q" }' | sac
    rm -f ./*.${family}?
    (
    cd Vel || exit 1
    sh ../../../programs/select/add_arrival.sh P *.${family}Z
    sh ../../../programs/select/add_arrival.sh S *.${family}[EN12]
    # horizontals whose metadata azimuths are not exactly 90 deg apart: the second one is set to
    # first + 90 (or - 90), whichever is closer
    for pair in 12 EN; do
        first=$(echo $pair | cut -c1); second=$(echo $pair | cut -c2)
        ls -1 *${family}${first} 2>/dev/null | gawk 'BEGIN{FS="."} {print $1"."$2"."$3"."}' > list_first
        ls -1 *${family}${second} 2>/dev/null | gawk 'BEGIN{FS="."} {print $1"."$2"."$3"."}' > list_second
        comm -12 list_first list_second > ${family}${pair}list.dat
        while read stem; do
            saclst cmpaz f ${stem}${family}${first} ${stem}${family}${second} | gawk 'NR==1{a=$2} NR==2{b=$2; f=$1} END{
                d1 = (a + 360 - b) % 360; d2 = (b + 360 - a) % 360
                if (d1 < d2) { if (d1 != 90) print "r "f"\nch cmpaz "(a + 270) % 360"\nwh" }
                else         { if (d2 != 90) print "r "f"\nch cmpaz "(a + 90) % 360"\nwh" } }'
        done < ${family}${pair}list.dat | { cat; echo q; } | sac > /dev/null
    done
    rm -f list_first list_second
    # SAC's rotate refuses two horizontals of different length (ERROR 2010), and the next write then
    # stores the unrotated 1 and 2 records as north and east (the teleseismic form of this defect is
    # PITFALLS A1): cut such 1/2 pairs to their common window first
    while read stem; do
        window=$(saclst b e npts f ${stem}${family}1 ${stem}${family}2 | gawk 'NR==1{b1=$2; e1=$3; n1=$4} NR==2{b2=$2; e2=$3; n2=$4} END{if (n1 != n2) print (b1 > b2 ? b1 : b2), (e1 < e2 ? e1 : e2)}')
        [ -n "$window" ] && printf "cut %s\nr %s %s\nw over\ncut off\nq\n" "$window" ${stem}${family}1 ${stem}${family}2 | sac > /dev/null
    done < ${family}12list.dat
    mkdir -p enz_$family
    # vertical as is; E/N copied, 1/2 rotated to north/east (rot to 0)
    gawk -v f=$family '{print "r "$1 f "Z"; print "w enz_" f "/"$1"z"; print "r "$1 f "E " $1 f "N"; print "w enz_" f "/"$1"e enz_" f "/"$1"n"} END {print "q"}' ${family}ENlist.dat | sac > /dev/null
    gawk -v f=$family '{print "r "$1 f "Z"; print "w enz_" f "/"$1"z"; print "r "$1 f "1 " $1 f "2"; print "rot to 0"; print "w enz_" f "/"$1"n enz_" f "/"$1"e"} END {print "q"}' ${family}12list.dat | sac > /dev/null
    # check the property the rotation must create: a rotated north record has cmpaz 0 (a pair that
    # still failed to rotate keeps the azimuth of channel 1); such stations are set aside, not fitted
    mkdir -p enz_$family/others
    while read stem; do
        [ -f enz_$family/${stem}n ] || continue
        north_azimuth=$(saclst cmpaz f enz_$family/${stem}n | gawk '{print $2}')
        if gawk -v a="$north_azimuth" 'BEGIN{exit !(a > 0.01 && a < 359.99)}'; then
            echo "MyExtractSeed_loc: ${stem}${family}1/2 not rotated (north cmpaz $north_azimuth) - moved to enz_$family/others"
            mv enz_$family/${stem}[enz] enz_$family/others/ 2>/dev/null
        fi
    done < ${family}12list.dat
    cd enz_$family || exit 1
    # the signal-to-noise test of every component starts at the P pick (t2 = t1 of the vertical): the
    # 50 s before the S pick are P coda for deep sources (2026 Tonga, 234 km at 5-7 deg: S-P ~63 s,
    # horizontal "SNR" 1.5-1.9 dB, every regional station rejected)
    saclst t1 f ./*.z 2>/dev/null | gawk '{stem = $1; sub(/z$/, "", stem); sub(/^\.\//, "", stem)
        print "r " stem "e " stem "n " stem "z"; print "ch t2 " $2; print "wh"} END {print "q"}' | sac > /dev/null
    perl ../../../../programs/select/selectSNRloc.pl
    rm -f ./*.[enz]
    mv ./good_data/*.[enz] . 2>/dev/null
    rm -rf good_data
    mkdir -p others
    for base in $(ls ./*.[enz] 2>/dev/null | sed -e 's|^\./||; s/\.[enz]$//' | sort -u); do
        if [ ! -f "$base.e" ] || [ ! -f "$base.n" ] || [ ! -f "$base.z" ]; then
            echo "Incomplete set for $base ($family), moving to others"
            mv "$base".[enz] others/ 2>/dev/null
        fi
    done
    sh ../../../../programs/select/decimate_to_1sps.sh enz
    ) || exit 1
}
for family in $families; do
    echo "MyExtractSeed_loc: processing $family? records"
    process_family $family
done

#---------------------------------------------------------------------------------------------------
# one instrument per station, then the common window and network renaming
#---------------------------------------------------------------------------------------------------
python ../../programs/select/choose_regional_instruments.py choose
cp Vel/enz_selected/*.[enz] data/ 2>/dev/null
cd data || exit 1
ls -1 ./*.[enz] 2>/dev/null | gawk '{print "cuterr fillz"; print "cut -50 600"; print "r "$1; print "ch t1 0"; print "w over"} END{print "q"}' | sac
# the Chilean network code C becomes C1 (the CSN records merged in step1_4_merge_csn.sh use C1)
ls -1 *.[enz] 2>/dev/null | gawk 'BEGIN{FS="."} {if ($1 == "C") print "mv "$0" C1."$2"."$3"."$4}' | sh
cd ..
rm -rf RESP
