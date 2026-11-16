#!/usr/bin/perl 
# regional selection (step 1): every component needs >= 2 dB signal-to-noise in the 50 s after the P
# pick (header t2, written by MyExtractSeed_loc.sh) against the 50 s before it
$snr="../../../../programs/select/snr";
$zr=2;
$er=2;
$nr=2;
`mkdir good_data`;
`mkdir bad_data`;
@HZ_data=<*.z>;
foreach $HZ (@HZ_data){
         ($tmp,$time,$khole,$com)=split(/\./,$HZ);
         $HE="$tmp"."\."."$time"."\."."$khole"."\."."e";
	 $HN="$tmp"."\."."$time"."\."."$khole"."\."."n";
if (-e $HZ) {

		($junk,$zsnr)=split('=',`$snr $HZ t2 0 50`);
		($junk,$esnr)=split('=',`$snr $HE t2 0 50`);
		($junk,$nsnr)=split('=',`$snr $HN t2 0 50`);
		if (($zsnr>= $zr) && ($esnr>=$er) && ($nsnr>=$nr)) {
			print "good zen component: $HZ SNR= $zsnr,$esnr,$nsnr\n";
			`cp $HZ $HE $HN good_data`;
		}
		else {
			print "low SNR: $zsnr,$esnr,$nsnr\n";
                        `cp $HZ $HE $HN bad_data`;
		}
	  }
}
