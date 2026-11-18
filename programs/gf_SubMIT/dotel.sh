#!/bin/bash
# both work folders are created here: git does not track empty folders (cohort report item 5)
mkdir -p greenFuncDir greenFuncDir_disp
vmodel_folders=$(find greenFuncDir/ greenFuncDir_disp/ -maxdepth 1 -type d -name "vmodel_*")
if [ -n "$vmodel_folders" ]; then
    echo "Deleting vmodel_* folders..."
    rm -rf $vmodel_folders
fi

path=$(pwd)
nm=$(echo "$path" | awk -F'/' '{print $(NF)}' | sed 's/^gf_//')
invd=../inv_${nm}; [ -d "$invd" ] || invd=../inv_${nm}_1sub
SUBMIT_DELTAT=$(gawk -F= '/^DELTAT=/{print $2}' ../programs/submit.conf 2>/dev/null)
export SUBMIT_DELTAT=${SUBMIT_DELTAT:-1}
# The depth grid is the one of the regional Green's functions (loc/distdep.dat, written by
# loc/doloc.sh, which runs first): event-adaptive aperture, start in rock, or the published grid of
# an example. Until 2026-11-17 this script derived the grid again with a copy of doloc.sh's rule.
read -r _dist_min _dist_max _dist_step depth_min depth_max _depth_step < loc/distdep.dat \
    || { echo "dotel: no loc/distdep.dat - run loc/doloc.sh first"; exit 1; }
echo "dotel: depth grid $depth_min-$depth_max km (loc/distdep.dat)"
# UNION of P (.z) and SH (.t) distances: an SH-only station (P screened out,
# transverse kept) otherwise gets NO Green's functions and sub_init reads
# silent garbage for it (2024 Chile AU.MAW); forward_model.c now stops on a missing file.
# Every .z and .t record counts: the list used to be built from vel_*.z and vel_*.t, but there are
# no velocity transverse records, so a rebuild after the model-based screen had removed the vertical
# records of six California stations left their SH records without Green's functions (2026-11-17).
gcdistlst=$(saclst dist f  ${invd}/data/*.z ${invd}/data/*.t 2>/dev/null | gawk '{printf("%d\n",$2)}' | sort -n | uniq | tr '\n' ' ')

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

# greenFuncDir_disp/ is empty in the repository, and git does not track empty folders: a fresh clone
# has none, the copy below failed, and the clean-up further down then deleted every computed Green's
# function (cohort report 2026-10-23, item 5). Create it, and stop instead of deleting on failure.
mkdir -p greenFuncDir_disp
cp -r greenFuncDir/vmodel_* greenFuncDir_disp/ || { echo "dotel: copying the teleseismic Green's functions failed"; exit 1; }

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
# regional Green's functions join the teleseismic ones of the same depth; delete the originals only
# after every copy succeeded
copy_failed=0
for depth_folder in vmodel_*; do
    cp "$depth_folder"/*.grn.???? "../greenFuncDir_disp/$depth_folder/" || copy_failed=1
done
[ "$copy_failed" -eq 0 ] || { echo "dotel: copying the regional Green's functions failed - originals kept in loc/"; exit 1; }
rm -f vmodel_*/*.grn.*
cd ..
rm -rf greenFuncDir/vmodel_*
echo "finished calculating green's functions..."
