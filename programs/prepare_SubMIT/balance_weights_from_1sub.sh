#!/bin/bash
# Balance band weights so each data type's FINAL MISFIT CONTRIBUTION hits the
# target shares P:SH:Rayl = 1.2 : 1 : 1 (user rule 2026-08-10; replaces the
# ad-hoc "trust factor" discounts). Contributions are measured by running the
# forward operator on the 1-sub scan's best model; since a band's share scales
# ~ w^2, weights update by w <- w*sqrt(target/measured), iterated <=4 times to
# within 15%. Weights are renormalized to weightSH=1 and written into the
# n>=2 inversion dirs' Par.file (the 1-sub scout keeps its band-true start).
# GUARD: a filler regional set (<=2 stations, weight <=0.002) is excluded -
# its weight stays put and P:SH balance to 1.2:1 alone.
#
# usage: balance_weights_from_1sub.sh <inv_..._1sub dir> <event_name>
d1=$1
ev=$2
fwd="fwd_${ev}"
best=""
for f in "$d1"/best_model_hybrid.dat "$d1"/best_model_ensemble.dat \
         "$d1"/best_model_exploration.dat; do
    [ -s "$f" ] && best=$f && break
done
[ -n "$best" ] || { echo "balance_weights: no 1-sub best model - skipped"; exit 0; }
[ -d "$fwd" ] || { echo "balance_weights: no $fwd dir - skipped"; exit 0; }

# build ffwd in the fwd dir (rebuild when the binary predates the
# bandresid instrumentation)
if [ -x "$fwd/ffwd" ] && ! strings "$fwd/ffwd" 2>/dev/null | grep -q "bandresid"; then
    rm -f "$fwd/ffwd"
fi
if [ ! -x "$fwd/ffwd" ]; then
    cp ../programs/fwd_SubMIT/*.c ../programs/fwd_SubMIT/*.f90 \
       ../programs/fwd_SubMIT/Makefile "$fwd/" 2>/dev/null
    cp ../programs/code_SubMIT/sub_header.h "$fwd/" 2>/dev/null
    (cd "$fwd" && make > make_balance.log 2>&1)
    [ -x "$fwd/ffwd" ] || { echo "balance_weights: ffwd build failed - skipped"; exit 0; }
fi

# refresh the fwd dir's CONFIG from the freshly regenerated 1-sub dir: the
# plain fwd dir is a stale step2-era copy (Venezuela: stale stationsloc vs
# rebuilt GFs -> residualrayl=nan -> corrupted Par via empty sed variables)
for f in Par.file stations.info stationsSH.info stationsloc.info weights.dat distdep.dat; do
    [ -s "$d1/$f" ] && cp "$d1/$f" "$fwd/$f"
done
# 1-sub model -> Input.model (cen x y dura vr theta z), counts set to 1
gawk 'NR==1{printf "%.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", $4,$5,$6,$7,$8,$9,$10}' \
    "$best" > "$fwd/Input.model"
sed -i 's/^num_subevent=.*/num_subevent= 1 #/' "$fwd/Par.file"
cp "$fwd/Par.file" "$fwd/Par.file.balance_snapshot"
# the plain fwd dir is copied from inv_<event> BEFORE get_search_par runs, so
# search_par.file may not exist: borrow the 1-sub one (ffwd reads neq from it)
if [ ! -s "$fwd/search_par.file" ]; then
    cp "$d1/search_par.file" "$fwd/search_par.file" 2>/dev/null
fi
sed -i 's/^neq_min.*/neq_min 1/; s/^neq_max.*/neq_max 1/' "$fwd/search_par.file" 2>/dev/null

wp=$(gawk -F= '/^weightP=/{print $2}' "$fwd/Par.file" | gawk '{print $1}')
wsh=$(gawk -F= '/^weightSH=/{print $2}' "$fwd/Par.file" | gawk '{print $1}')
wr=$(gawk -F= '/^weightRayl=/{print $2}' "$fwd/Par.file" | gawk '{print $1}')
nloc=$(grep -c . "$fwd/stationsloc.info" 2>/dev/null || echo 0)
# filler = station COUNT only: a tiny band-true weight on a real network
# (California: 86 stations at 0.0014) is exactly what balancing corrects
fillr=$(gawk -v n="$nloc" 'BEGIN{print (n<=2) ? 1 : 0}')
[ "$fillr" -eq 1 ] && echo "balance_weights: filler regional (n=$nloc, w=$wr) - balancing P:SH only"

restore_and_skip() {
    cp "$fwd/Par.file.balance_snapshot" "$fwd/Par.file" 2>/dev/null
    echo "balance_weights: $1 - skipped (Par restored)"
    exit 0
}
for it in 1 2 3 4; do
    line=$(cd "$fwd" && FFWD_COMPACT=1 ./ffwd 2>/dev/null | grep "^bandresid:" | tail -1)
    [ -n "$line" ] || restore_and_skip "no bandresid line (fwd run failed?)"
    case "$line" in *nan*|*inf*|*NaN*) restore_and_skip "non-finite band residual: $line" ;; esac
    read -r wp wsh wr done_flag <<EOF2
$(echo "$line" | gawk -v wp="$wp" -v wsh="$wsh" -v wr="$wr" -v fill="$fillr" '{
    rp=$2; rsh=$3; rr=$4
    if (fill==1) { tot=rp+rsh; sp=rp/tot; ssh=rsh/tot; tp=1.2/2.2; tsh=1.0/2.2; sr=tp; tr=tp }
    else { tot=rp+rsh+rr; sp=rp/tot; ssh=rsh/tot; sr=rr/tot; tp=1.2/3.2; tsh=1.0/3.2; tr=1.0/3.2 }
    fp=sqrt(tp/(sp+1e-12)); fsh=sqrt(tsh/(ssh+1e-12)); fr=(fill==1)?1.0:sqrt(tr/(sr+1e-12))
    # converged when every balanced share is within 15% of target
    ok = (sp>0.85*tp && sp<1.15*tp && ssh>0.85*tsh && ssh<1.15*tsh)
    if (fill==0) ok = ok && (sr>0.85*tr && sr<1.15*tr)
    nwp=wp*fp; nwsh=wsh*fsh; nwr=(fill==1)?wr:wr*fr
    # renormalize to SH=1
    printf "%.6g %.6g %.6g %d", nwp/nwsh, 1.0, (fill==1)?wr:nwr/nwsh, ok
}')
EOF2
    # never write non-numeric weights (the corruption path): validate first
    okw=$(gawk -v a="$wp" -v b="$wsh" -v c="$wr" 'BEGIN{print (a+0>0 && b+0>0 && c+0>0) ? 1 : 0}')
    [ "$okw" = "1" ] || restore_and_skip "malformed weight update ($wp/$wsh/$wr)"
    sed -i "s/^weightP=.*/weightP= $wp #/" "$fwd/Par.file"
    sed -i "s/^weightSH=.*/weightSH= $wsh #/" "$fwd/Par.file"
    [ "$fillr" -eq 0 ] && sed -i "s/^weightRayl=.*/weightRayl= $wr/" "$fwd/Par.file"
    [ "$done_flag" -eq 1 ] && break
done
echo "balance_weights: converged weights P:SH:Rayl = $wp : $wsh : $wr (iter $it)"
for n in 2 3 4 5 6; do
    pf="inv_${ev}_${n}sub/Par.file"
    [ -s "$pf" ] || continue
    sed -i "s/^weightP=.*/weightP= $wp #/" "$pf"
    sed -i "s/^weightSH=.*/weightSH= $wsh #/" "$pf"
    [ "$fillr" -eq 0 ] && sed -i "s/^weightRayl=.*/weightRayl= $wr/" "$pf"
done
echo "balance_weights: applied to n>=2 Par files"
