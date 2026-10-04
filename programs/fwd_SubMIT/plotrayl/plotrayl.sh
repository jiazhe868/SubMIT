rm rayl*sac
evlo=$(cat ../Par.file | gawk '{if (NR==29) print $2}')
evla=$(cat ../Par.file | gawk '{if (NR==30) print $2}')
nm=$(cat ../stationsloc.info | wc -l)
cat ../stationsloc.info  | gawk -v elo="$evlo" -v ela="$evla" -v nn="$nm" '{if (NR<=nn) print $1,sqrt((($2-elo)*cos(ela/180*3.1416))^2+($3-ela)^2),atan2($2-elo,$3-ela)*180/3.141592653}' | gawk '{if ($3<0) $3=$3+360;printf("%s %.1f %d %04d %d\n",$1,$2,$3,NR-1,1)}' | gawk 'BEGIN{FS="/"} {print $2}' > sta_dist_az_rev1
nm=$(cat sta_dist_az_rev1 | wc -l)
cp ../waveforms/rayl*sac .
cat sta_dist_az_rev1 | gawk '{print $1}' | gawk 'BEGIN{FS="[.]"} {print $1"."$2}' > sta1etemp
cat sta_dist_az_rev1 | gawk '{print $1}' | gawk 'BEGIN{FS="[.]"} {print $1"."$2}' > sta1ntemp
cat sta_dist_az_rev1 | gawk '{print $1}' | gawk 'BEGIN{FS="[.]"} {print $1"."$2}' > sta1ztemp
cat sta_dist_az_rev1 | gawk '{print "1 1 1 rayle_syn_"$4".sac 0 0"}' > liste1
cat sta_dist_az_rev1 | gawk '{print "rayle_obs_"$4".sac 0 0"}' > liste2
cat sta_dist_az_rev1 | gawk '{print "1 1 1 rayln_syn_"$4".sac 0 0"}' > listn1
cat sta_dist_az_rev1 | gawk '{print "rayln_obs_"$4".sac 0 0"}' > listn2
cat sta_dist_az_rev1 | gawk '{print "1 1 1 raylz_syn_"$4".sac 0 0"}' > listz1
cat sta_dist_az_rev1 | gawk '{print "raylz_obs_"$4".sac 0 0"}' > listz2
paste sta_dist_az_rev1 liste1 liste2 sta1etemp    > allinfoe
paste sta_dist_az_rev1 listn1 listn2 sta1ntemp   > allinfon
paste sta_dist_az_rev1 listz1 listz2 sta1ztemp   > allinfoz
cp allinfoe allinfo
cat allinfon >> allinfo
cat allinfoz >> allinfo
bt=$(cat ../Par.file | gawk '{if (NR==5) print $2}')
echo "${bt}"
echo "r *.sac\nch b ${bt}\nch t1 0\nch dist 5000\nch delta 1\nwh\nq\n"  | sac
#ls -1 *.sac | gawk '{print "r "$1;print "int"; print "w over"} END{print "q"}' | sac
#ls -1 *_000[0-7].sac | gawk '{print "r "$1;print "mul 0.3";print "w over"} END{print "q"}' | sac
#nrow1=$(cat allinfo| wc -l | gawk '{a=$1;b=$1%3;c=a-b+3;d=c/3;print d;}')
#nrow2=$(cat allinfo| wc -l | gawk '{a=$1;b=$1%3;c=a-b+3;d=c/3;print d;}')
#nrow3=$(cat allinfo| wc -l | gawk '{a=$1;b=$1%3;c=a-b+3;d=c/3;print a-d-d;}')
nrow1=$nm
nrow2=$nm
nrow3=$nm
dy=$(echo $nrow1 | gawk '{print 10/$1}')
tb=$(cat ../Par.file | gawk '{if (NR==5) print $2+10}')
te=$(cat ../Par.file | gawk '{if (NR==6) print $2-10}')

txtloc=$(echo $tb | gawk '{print $1-20}')
am=$(saclst depmax f *.sac | minmax | gawk 'BEGIN{FS="[</>]"} {print $8*3}')

cat allinfo | gawk '{print $9,$10,0}' | gawk -v txtb="$txtloc" -v nrow1="$nrow1" -v dy="$dy" -v tb="$tb" -v te="$te" -v am="$am" '{if (NR<=nrow1) {if (NR==1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow1+1" -JX2i/9.9i -X0.8i -Y1i -P -K  -Ent1 -W1.5p/255/0/0 -M"am"/0 > map.ps << EOF"; print $1,$2,nrow1-NR+1; print "EOF";} if (NR>1&&NR<nrow1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow1+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/255/0/0 -M"am"/0  >> map.ps << EOF"; print $1,$2,nrow1-NR+1; print "EOF";}  if (NR==nrow1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow1+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/255/0/0 -M"am"/0 -Ba100f50/S >> map.ps << EOF"; print $1,$2,nrow1-NR+1; print "EOF";}}}'  | sh  > /dev/null
cat allinfo | gawk '{print $12,$13,$15,$2,$3}' | gawk -v txtb="$txtloc" -v nrow1="$nrow1" -v dy="$dy" -v tb="$tb" -v te="$te" -v am="$am" '{if (NR<=nrow1) {if (NR==1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow1+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0 >> map.ps << EOF"; print $1,$2,nrow1-NR+1; print "EOF"; print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1-NR+0.8,$5,$4);print "EOF";} if (NR>1&&NR<nrow1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow1+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0  >> map.ps << EOF"; print $1,$2,nrow1-NR+1; print "EOF";print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1-NR+0.8,$5,$4);print "EOF";} if (NR==nrow1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow1+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0  >> map.ps << EOF"; print $1,$2,nrow1-NR+1; print "EOF";print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1-NR+0.8,$5,$4);print "EOF";}}}'| sh  > /dev/null

cat allinfo | gawk '{print $9,$10,0}' | gawk -v txtb="$txtloc" -v nrow1="$nrow1" -v nrow2="$nrow2" -v dy="$dy" -v tb="$tb" -v te="$te" -v am="$am" '{if (NR>nrow1&&NR<=nrow1+nrow2) {if (NR==nrow1+1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow2+1" -JX2i/9.9i -X2.4i  -P -K -O -Ent1 -W1.5p/255/0/0 -M"am"/0 >> map.ps << EOF"; print $1,$2,nrow1+nrow2-NR+1; print "EOF";} if (NR>nrow1+1&&NR<nrow1+nrow2) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow2+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/255/0/0 -M"am"/0  >> map.ps << EOF"; print $1,$2,nrow1+nrow2-NR+1; print "EOF";} if (NR==nrow1+nrow2) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow2+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/255/0/0 -M"am"/0 -Ba100f50/S >> map.ps << EOF"; print $1,$2,nrow1+nrow2-NR+1; print "EOF";}}}' | sh  > /dev/null
cat allinfo | gawk '{print $12,$13,$15,$2,$3}' | gawk -v txtb="$txtloc" -v nrow1="$nrow1" -v nrow2="$nrow2" -v dy="$dy" -v tb="$tb" -v te="$te" -v am="$am" '{if (NR>nrow1&&NR<=nrow1+nrow2) {if (NR==nrow1+1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow2+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0 >> map.ps << EOF"; print $1,$2,nrow1+nrow2-NR+1; print "EOF";print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1+nrow2-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1+nrow2-NR+0.8,$5,$4);print "EOF";} if (NR>nrow1+1&&NR<nrow1+nrow2) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow2+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0  >> map.ps << EOF"; print $1,$2,nrow1+nrow2-NR+1; print "EOF";print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1+nrow2-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1+nrow2-NR+0.8,$5,$4);print "EOF";} if (NR==nrow1+nrow2) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow2+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0 >> map.ps << EOF"; print $1,$2,nrow1+nrow2-NR+1; print "EOF";print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1+nrow2-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1+nrow2-NR+0.8,$5,$4);print "EOF";}}}' | sh   > /dev/null
cat allinfo | gawk '{print $9,$10,0}' | gawk -v txtb="$txtloc" -v nrow1="$nrow1" -v nrow2="$nrow2"  -v nrow3="$nrow3" -v dy="$dy" -v tb="$tb" -v te="$te" -v am="$am" '{if (NR>nrow1+nrow2&&NR<=nrow1+nrow2+nrow3) {if (NR==nrow1+nrow2+1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow3+1" -JX2i/9.9i -X2.4i  -P -K -O -Ent1 -W1.5p/255/0/0 -M"am"/0 >> map.ps << EOF"; print $1,$2,nrow1+nrow2+nrow3-NR+1; print "EOF";} if (NR>nrow1+nrow2+1&&NR<nrow1+nrow2+nrow3) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow3+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/255/0/0 -M"am"/0  >> map.ps << EOF"; print $1,$2,nrow1+nrow2+nrow3-NR+1; print "EOF";}  if (NR==nrow1+nrow2+nrow3) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow3+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/255/0/0 -M"am"/0 -Ba100f50/S >> map.ps << EOF"; print $1,$2,nrow1+nrow2+nrow3-NR+1; print "EOF";}}}' | sh  > /dev/null
cat allinfo | gawk '{print $12,$13,$15,$2,$3}' | gawk -v txtb="$txtloc" -v nrow1="$nrow1" -v nrow2="$nrow2"  -v nrow3="$nrow3" -v dy="$dy" -v tb="$tb" -v te="$te" -v am="$am" '{if (NR>nrow1+nrow2&&NR<=nrow1+nrow2+nrow3) {if (NR==nrow1+nrow2+1) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow3+1" -JX2i/9.9i  -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0 >> map.ps << EOF"; print $1,$2,nrow1+nrow2+nrow3-NR+1; print "EOF";print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1+nrow2+nrow3-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1+nrow2+nrow3-NR+0.8,$5,$4);print "EOF";} if (NR>nrow1+nrow2+1&&NR<nrow1+nrow2+nrow3) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow3+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0  >> map.ps << EOF"; print $1,$2,nrow1+nrow2+nrow3-NR+1; print "EOF";print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1+nrow2+nrow3-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1+nrow2+nrow3-NR+0.8,$5,$4);print "EOF";} if (NR==nrow1+nrow2+nrow3) {print "pssac2 -C-200/2000 -R"tb"/"te"/0/"nrow3+1" -JX2i/9.9i -P -K -O -Ent1 -W1.5p/0/0/0 -M"am"/0 >> map.ps << EOF"; print $1,$2,nrow1+nrow2+nrow3-NR+1; print "EOF";print "pstext -N -N -N -R -J -P -K -O >> map.ps << EOF"; print txtb,nrow1+nrow2+nrow3-NR+1.2" 9 0 1 ML "$3;printf (txtb" %f 10 0 0 ML %.1f\260 %.1f\260\n",nrow1+nrow2+nrow3-NR+0.8,$5,$4);print "EOF";}}}' | sh  > /dev/null
