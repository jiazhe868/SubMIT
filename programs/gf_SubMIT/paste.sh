model=$1
i=$2
if [ -d ./${model}_$i ];then
    rm -r ${model}_$i
fi
mkdir ${model}_$i
saclst t3 f ${model}direct_$i/*.[0-8] | gawk -v model=$model '{split($1,aa,"_");print"cuterr fillz;cut "$2-50,$2+600 ;print"r "$1,model"core_"aa[2];print"w over";print"r "$1;print"addf "model"core_"aa[2];print"w "model"_"aa[2]}END{print"quit"}' | sac
saclst t3 f ${model}direct_$i/*.[abc] | gawk -v model=$model '{split($1,aa,"_");print"cuterr fillz;cut "$2-50,$2+600 ;print"r "$1,model"core_"aa[2];print"w over";print"r "$1;print"addf "model"core_"aa[2];print"w "model"_"aa[2]}END{print"quit"}' | sac
saclst t3 f ${model}direct_$i/*.???[rtz] | gawk -v model=$model '{split($1,aa,"_");print"cuterr fillz;cut "$2-50,$2+600 ;print"r "$1,model"core_"aa[2];print"w over";print"r "$1;print"addf "model"core_"aa[2];print"w "model"_"aa[2]}END{print"quit"}' | sac
rm -rf ${model}direct_$i ${model}core_$i
