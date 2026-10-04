#! /bin/bash
gmtset LABEL_FONT_SIZE 18p
gmtset ANNOT_FONT_SIZE_PRIMARY 16p
gmtset ANNOT_OFFSET_PRIMARY 0.15c
gmtset TICK_LENGTH -0.15c
gmtset GRID_PEN_PRIMARY 0.24.3_2:p
gmtset BASEMAP_TYPE plain
gmtset LABEL_OFFSET 0.2c
gmtset FRAME_PEN 1.6
gmtset COLOR_NAN white

cp ../Input.model ../stf/stf_*.sac ../fm.dat .

x1=$(minmax Input.model | gawk 'BEGIN{FS="[</>]"} {print $5-30}')
x2=$(minmax Input.model | gawk 'BEGIN{FS="[</>]"} {print $6+30}')
y1=$(minmax Input.model | gawk 'BEGIN{FS="[</>]"} {print $8-30}')
y2=$(minmax Input.model | gawk 'BEGIN{FS="[</>]"} {print $9+30}')
scale=$(echo $x1 $x2 $y1 $y2 | gawk '{a1=$2-$1;a2=$4-$3;len=(a1>a2)?a1:a2;print len/5}')
xt=$(echo $x1 | gawk -v ll="$scale" '{print $1+ll/1.5}')
yt=$(echo $y1 | gawk -v ll="$scale" '{print $1+ll/1.5}')

aa=$(echo $x1 $x2 | gawk -v ll="$scale" '{print ($2-$1)/ll"i"}')
bb=$(echo $y1 $y2 | gawk -v ll="$scale" '{print ($2-$1)/ll"i"}')
cc=$(echo $y1 $y2 | gawk -v ll="$scale" '{print ($2-$1)/ll-1"i"}')

psbasemap -R${x1}/${x2}/${y1}/${y2} -JX${aa}/${bb} -X1.5i -Y2i -P -K -Ba20f10:"X (km)":/a20f10:"Y (km)":WSen > fig.ps
cat fm.dat | gawk '{print sqrt(($2^2+2*$3^2+2*$4^2+$5^2+2*$6^2+$7^2)/2)"e27",$2,$3,$4,$5,$6,$7}' > FM.dat
cat fm.dat | gawk 'BEGIN{a2=0;a3=0;a4=0;a5=0;a6=0;a7=0} {a2=a2+$2;a3=a3+$3;a4=a4+$4;a5=a5+$5;a6=a6+$6;a7=a7+$7;} END{print sqrt(($2^2+2*a3^2+2*a4^2+a5^2+2*a6^2+a7^2)/2)"e27",a2,a3,a4,a5,a6,a7}' > sumFM.dat
cat FM.dat | gawk '{printf("%.1f\n", 2/3*log($1)/log(10)-10.7)}' > FM_mw.dat
cat sumFM.dat | gawk '{printf("%.1f", 2/3*log($1)/log(10)-10.7)}' > sumFM_mw.dat
cat Input.model | gawk '{print $1,$2+0,$3+0,$4,1,0,1}' > Input.model1
scalarm0=$(cat sumFM.dat  | gawk '{print log($1)/log(10)}')
sh gentimes.sh > times.dat
cd ../stf/
ls -1 stf_????.sac > ../subeveplot/stf.dat
cd ../subeveplot/
paste stf.dat times.dat  | gawk '{print "r "$1;print "ch delta 1 b -10 t1 0 dist 5000";print "mul "$2;print "w over"} END{print "q"}' | sac
mamp=$(saclst depmax f stf_000*.sac  | minmax | gawk 'BEGIN{FS="[</>]"} {print $8}')
paste Input.model1 FM.dat FM_mw.dat | gawk -v m0="$scalarm0" '{print $2,$3, 10,$14,$9,$12,$11,-$13,-$10,log($8)/log(10),$2,$3,"E"NR,"Mw "$15}' | gawk -v m0="$scalarm0" '{if (NR==1) cl="red";if (NR==2) cl="yellow";if (NR==3) cl="green";if (NR==4) cl="0/255/255";if (NR==5) cl="blue";if (NR==6) cl="purple";print "echo \""$0"\" |  psmeca -R -J -P -K -O -G"cl" -Sm"($10-m0+2.5)*0.6"/14 -T0 -p -t  >> fig.ps"}'  | tac| sh 
paste sumFM.dat sumFM_mw.dat | gawk -v x="$xt" -v y="$yt" -v m0="$scalarm0" '{print x,y,10,$7,$2,$5,$4,-$6,-$3,log($1)/log(10),x,y,"Total Mw "$8}'  | gawk -v m0="$scalarm0" '{print "echo \""$0"\" | psmeca -R -J -P -K -O -G100/100/100 -Sm"($10-m0+2.5)*0.6"/14 -T0 -p -t  >> fig.ps"}' | sh
psxy -R -J -P -K -O -Sa0.3i -G0/0/0 >> fig.ps << END
0 0
END
psbasemap -R-5/65/0/1 -JX2i/1i -Y${cc} -P -K -O -Ba20f10:"Time (s)":eS >> fig.ps
cat stf.dat| gawk -v mamp="$mamp" '{if (NR==1) cl="255/0/0";if (NR==2) cl="255/255/0";if (NR==3) cl="0/255/0";if (NR==4) cl="0/255/255";if (NR==5) cl="0/0/255";if (NR==6) cl="255/0/255";print "echo \""$1,0,0"\" | pssac2 -R -J -P -K -O -W1p/0/0/0 -G"cl"/0/-5/65 -C-5/65 -M"mamp"/0 -Ent1 >> fig.ps";}' | sh
