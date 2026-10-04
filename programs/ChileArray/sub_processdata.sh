tt=$1
mm=$2
nn=$3
cp convert_to_sac.py ${tt}_${mm}_${nn}
cd ${tt}_${mm}_${nn}
ls -1 *[ENZ].txt > data_list.dat
python convert_to_sac.py
ls -1 *.[enz] | gawk '{print "cuterr fillz";print "cut -50 400";print "r "$1;print "w over"} END{print "q"}' | sac
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
if(sqrt(($2-0.1)^2)<1e-4){
                print "r",$1;
                print "decimate 2";
                print "decimate 5";
                print "w over";
        }
if(sqrt(($2-0.2)^2)<1e-4){
                print "r",$1;
                print "decimate 5";
                print "w over";
        }

} END{print "quit"}' | sac

cd ..
