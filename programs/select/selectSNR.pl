#!/usr/bin/perl 
$snr="../../../../programs/select/snr";
$zr=3.0;
$tr=3.0;
$rr=0.0;
`mkdir good_data`;
`mkdir bad_data`;
@HZ_data=<*.z>;
foreach $HZ (@HZ_data){
         ($tmp,$time,$khole,$com)=split(/\./,$HZ);
         $HT="$tmp"."\."."$time"."\."."$khole"."\."."t";
	 $HR="$tmp"."\."."$time"."\."."$khole"."\."."r";
if (-e $HZ) {

		($junk,$zsnr)=split('=',`$snr $HZ t1 0 50`);
		($junk,$tsnr)=split('=',`$snr $HT t1 0 50`);
		if (($zsnr>= $zr) && ($tsnr>=$tr)) {
			print "good z component: $HZ SNR=$zsnr\n";
			print "good t component: $HT SNR=$tsnr\n";
			`cp $HZ $HT $HR good_data`;
		}
		else {
			print "low SNR: $HZ SNR=$zsnr $HT SNR=$tsnr\n";
                        `cp $HZ $HT $HR bad_data`;
		}
	  }
}
