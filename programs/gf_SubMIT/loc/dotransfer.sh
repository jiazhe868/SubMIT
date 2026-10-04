cat distlst | gawk '{print "r "$1".grn.0";print "div 1e10";print "w "$1".grn.zero";} END{print "q"}'  | sac
cat distlst | gawk '{print "sh transfer.sh "$1}' | sh
#rm *.[1247]
cat distlst | gawk '{print "r "$1".grn.????";print "ch t1 0";print "wh"} END{print "q"}' | sac
