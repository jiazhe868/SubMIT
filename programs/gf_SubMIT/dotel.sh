#!/bin/bash
vmodel_folders=$(find greenFuncDir/ greenFuncDir_disp/ -maxdepth 1 -type d -name "vmodel_*")
if [ -n "$vmodel_folders" ]; then
    echo "Deleting vmodel_* folders..."
    rm -rf $vmodel_folders
fi

path=$(pwd)
nm=$(echo "$path" | awk -F'/' '{print $(NF)}' | sed 's/^gf_//')
invd=../inv_${nm}; [ -d "$invd" ] || invd=../inv_${nm}_1sub
evdp=$(saclst evdp f  ${invd}/data/*.z | sort -nk2 | head -1| gawk '{print int($2)}')
SUBMIT_DELTAT=$(gawk -F= '/^DELTAT=/{print $2}' ../programs/submit.conf 2>/dev/null)
export SUBMIT_DELTAT=${SUBMIT_DELTAT:-1}
# event-adaptive depth aperture: a deep intraslab rupture (evdp > 60 km) can
# run tens of km up-/down-dip inside the steeply dipping slab (2024 Chile M7.4:
# subevents to evdp+51 km), so the fixed +/-20 grid clips it; scale the
# half-width with the Wells-Coppersmith length. Shallow events keep +/-20
# (depth extent bounded by seismogenic width; crustal cap in get_search_par).
mw=$(gawk -v ev="$nm" '$1==ev{print $4; exit}' ../event_list.dat)
halfw=20
if [ "$evdp" -gt 60 ] && [ -n "$mw" ]; then
    halfw=$(gawk -v mw="$mw" 'BEGIN{h=int(0.65*10^(0.59*mw-2.44)+0.5); if(h<20)h=20; if(h>70)h=70; print h}')
fi
depth_min=$((evdp - halfw))
if [ "$depth_min" -lt 3 ]; then
    depth_min=3
fi
depth_max=$((evdp + halfw))
if [ -n "$SUBMIT_GF_INPUTS" ] && [ -f "$SUBMIT_GF_INPUTS/distdep.dat" ]; then
    read _a _b _c depth_min depth_max _d < "$SUBMIT_GF_INPUTS/distdep.dat"
    echo "dotel: using published depth grid $depth_min-$depth_max km"
fi
# UNION of P (.z) and SH (.t) distances: an SH-only station (P screened out,
# transverse kept) otherwise gets NO Green's functions and sub_init reads
# silent garbage for it (2024 Chile AU.MAW, fix 31)
gcdistlst=$(saclst dist f  ${invd}/data/vel_*.z ${invd}/data/vel_*.t 2>/dev/null | gawk '{printf("%d\n",$2)}' | sort -n | uniq | tr '\n' ' ')

depi=2
model=vmodel
model_ray=raysiasp

gawk -v dp1="$depth_min" -v dp2="$depth_max" -v depi="$depi" 'BEGIN{for (i=dp1;i<=dp2;i=i+depi) {print "mkdir greenFuncDir/vmodel_"i;print "cp raysiasp  vmodel paste.sh mtel3.pl greenFuncDir/vmodel_"i}}' | sh
gawk -v dp1="$depth_min" -v dp2="$depth_max" -v depi="$depi" -v gcdist="$gcdistlst" -v model="$model" -v model_ray="$model_ray" 'BEGIN {
    for (i=dp1; i<=dp2; i=i+depi) {
        dir="greenFuncDir/vmodel_"i;
        cmd = "time ./mtel3.pl -M" model "/" model_ray "/" i " -O. " gcdist;
        print cmd > dir "/tel3.sh";
    } 
}'
# cap concurrent GF jobs (override with SUBMIT_NCORE); previously unbounded "&"
NCORE=${SUBMIT_NCORE:-$(( $(nproc) / 2 > 0 ? $(nproc) / 2 : 1 ))}
gawk -v dp1="$depth_min" -v dp2="$depth_max" -v depi="$depi" -v model="$model" '
BEGIN {
    for (i = dp1; i <= dp2; i += depi) {
        print "cd greenFuncDir/vmodel_" i " && sh tel3.sh && sh paste.sh "model" "i" && mv vmodel_"i"/*.grn.???? ./ && rm -rf vmodel_"i"/ && cd ../../"
    }
}' | xargs -d '\n' -P "$NCORE" -I{} sh -c '{}' 

cp -r greenFuncDir/vmodel_* greenFuncDir_disp/

ls -1 -d greenFuncDir_disp/vmodel_* | while read dir; do
    while [ "$(jobs -rp | wc -l)" -ge "$NCORE" ]; do sleep 1; done
    (
        cd "$dir" || exit
        ls -1 ????*.grn.???? | awk '{print "r "$1;print "int";print "w over";} END{print "q"}' | sac
        cd - > /dev/null || exit
    ) &
done
wait

cd greenFuncDir/
ls -1 -d vmodel_*  | gawk '{print "cd "$1;print "ls -1 ????*.grn.???? | awk '\''{print \"cp \"$1\" ../../greenFuncDir_disp/"$1"/vel_\"$1}'\'' | sh"; print "cd ../"}' | sh

cd ../loc
ls -1 -d vmodel_* | gawk '{print "cp "$1"/*.grn.???? ../greenFuncDir_disp/"$1; print "rm "$1"/*.grn.*"}' | sh
cd ..
rm -rf greenFuncDir/vmodel_*
echo "finished calculating green's functions..."
