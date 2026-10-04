/***************************************************************
 * * copyright: 2014-, Earth Science School of USTC
 * * File name: snr.c
 * * Description: program to calculate signal-to-noise ratio for
 * *              seismograms.
 * * Author:Jiazhe
 * * Date:08/07/2014
 * * History:08/07/2014 v0.1
 * *         03/01/2015 v1.0
 * ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sacio.h>

int main(int argc,char *argv[]){
	int NT=500000;
	if (argc<=4) {
		fprintf(stderr,"Usage: snr file_insac t_mark (t1/t2/t5) error_time time_window_length\n     error_time:Ignore a window of error_time close to t_mark. Has to be an integer\n     time_window_length:The length of time windows for calculation of SNR from each side of t_mark. Has to be an integer.)\n");
                exit(1);
	}
	char *filename=argv[1];
	char *tmark=argv[2];
	int error_time=atoi(argv[3]);
	int time_window_length=atoi(argv[4]);
//	printf ("%s %s %d\n",filename,tmark,time_window);
	float b0,t_bg,delta,u[NT],noise_mean=0,signal_mean=0;
	double noise_var=0,signal_var=0;
	double SNR;
//	nerr is the mark of error in reading SAC file.
	int i, npts, n_t1, n_error, n_length, nerr;
	rsac1(filename,u,&npts,&b0,&delta,&NT,&nerr,strlen(filename));
	if (nerr!=0) {
		fprintf(stderr,"Error reading the SAC file!\n");
		exit(1);
	}
//	read t1 or t2 or t5 mark in the sac file.
	if (strcmp(tmark,"t1")==0||strcmp(tmark,"T1")==0) {
		getfhv("T1",&t_bg,&nerr,strlen("T1"));
		if (nerr!=0) {
			fprintf(stderr,"Error reading %s in the SACHEAD!\n", tmark);
			exit(1);
		}
	}
	if (strcmp(tmark,"t2")==0||strcmp(tmark,"T2")==0) {
                getfhv("T2",&t_bg,&nerr,strlen("T2"));
                if (nerr!=0) {
                        fprintf(stderr,"Error reading %s in the SACHEAD!\n", tmark);
                        exit(1);
                }
        }
	if (strcmp(tmark,"t5")==0||strcmp(tmark,"T5")==0) {
                getfhv("T5",&t_bg,&nerr,strlen("T5"));
                if (nerr!=0) {
                        fprintf(stderr,"Error reading %s in the SACHEAD!\n", tmark);
                        exit(1);
                }
        }
//	read delta of the sac file.
	getfhv("DELTA",&delta,&nerr,strlen("DELTA")); 
	if (nerr!=0) {
		fprintf(stderr,"Error reading delta t in the SACHEAD!\n");
		exit(1);
	}
//	The number of data points for t1.
	n_t1=(int)((t_bg-b0)/delta);
	n_length=(int)(time_window_length/delta);
	n_error=(int)(error_time/delta);
//	Make sure that there can be time windows for SNR calculation.
	if (n_t1-n_error<=0||n_t1+n_error>=npts) {
		fprintf(stderr,"No time window for calculation of SNR\n");
		exit(1);
	}
//	If the time_window_length is too long for the data, then make it shorter.
	if (n_t1-n_error<n_length){
		n_length=n_t1-n_error;
	}
	for (i=1;i<=n_length;i++) {
		signal_mean=signal_mean+u[n_t1+n_error+i-1];
		noise_mean=noise_mean+u[n_t1-n_error-i+1];
	}
	signal_mean=signal_mean/n_length;
	noise_mean=noise_mean/n_length;
	for (i=1;i<=n_length;i++) {
		signal_var=signal_var+(u[n_t1+n_error+i-1]-signal_mean)*(u[n_t1+n_error+i-1]-signal_mean);
		noise_var=noise_var+(u[n_t1-n_error-i+1]-noise_mean)*(u[n_t1-n_error-i+1]-noise_mean);
	}
	signal_var=signal_var/n_length;
	noise_var=noise_var/n_length;
	SNR=10*log10(signal_var/noise_var);
	printf("SNR=%f\n",SNR);
}

