/***************************************************************
 * copyright: 2017-, GPS in Caltech
 * File name: sub_forward.c
 * Description: Source characterization with sub-events
 * Author:Jiazhe
 * Date:03/28/2017
 * History:03/28/2017 v0.1
 ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "sacsubc.h"
#include "sub_header.h"

void outputsac(int npts, float *arr, float dt, char *filename)
{
        int nerr;
        float b, e, depmax, depmin, depmen;
                scmxmn(arr,npts,&depmax,&depmin,&depmen);
                newhdr();
                setfhv("DEPMAX", depmax, &nerr);
                setfhv("DEPMIN", depmin, &nerr);
                setfhv("DEPMEN", depmen, &nerr);
                setnhv("NPTS    ",npts,&nerr);
                setfhv("DELTA   ",dt  ,&nerr);
                b = 0;
                setfhv("B       ",b  ,&nerr);
                setihv("IFTYPE  ","ITIME   ",&nerr);
                e = b + (npts -1 )*dt;
                setfhv("E       ",e     ,&nerr);
                setlhv("LEVEN   ",1,&nerr);
                setlhv("LOVROK  ",1,&nerr);
                setlhv("LCALDA  ",1,&nerr);
                setnhv("NZYEAR", 1970, &nerr);
                setnhv("NZJDAY", 1, &nerr);
                setnhv("NZHOUR", 0, &nerr);
                setnhv("NZMIN" , 0, &nerr);
                setnhv("NZSEC" , 0, &nerr);
                setnhv("NZMSEC", 0, &nerr);
                bwsac(npts,filename,arr);
}

void genstf(float gdeltat, float gsig, float gmu, int stfnpts, float *tdata1, float dura, int nn) {
        int k;
        float variance2 = gsig*gsig*2.0;
        float gaussfun[stfnpts];
        int bgpt=(int)((gmu-dura/1.5)/gdeltat);
        int ndpt=(int)((gmu+dura/1.5)/gdeltat);
        if (bgpt<0) {bgpt=0;}
        if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
        for (k=0;k<stfnpts;k++){
                gaussfun[k]=1*gdeltat/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2);
        }
        char BB[]=".sac"; char temps[]="stf_"; char tempp[]="_";
        char aaa[64];
        sprintf(aaa,"%s%04d%s",temps,nn,BB);
        outputsac(stfnpts, gaussfun, gdeltat, aaa);
}

void genstf1(int gnpts, float gdeltat, float gsig, float gmu, int stfnpts, float *tdata1, float dura, float times, float *submrf) {
        int k;
        float variance2 = gsig*gsig*2.0;
        float gaussfun[stfnpts];
        int bgpt=(int)((gmu-dura/1.5)/gdeltat);
        int ndpt=(int)((gmu+dura/1.5)/gdeltat);
        if (bgpt<0) {bgpt=0;}
        if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
        for (k=0;k<stfnpts;k++){
                submrf[k]=1/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2)*times;
        }
}

void output_forward(int num_subevent, int num_sta, modelformat modelvec[maxstanum], int npts, int stfnpts, float deltat, float *tdata, double *tempgxm, double *rclvd, double *GG, double *tempdx, float *CCC, float *l2res, float residual, float normresidual, float finalresidual, float *T_Cross){
	int i,j,k,l,ii,jj;
	float templog, realsig, ssig, realmu,tempresidual, xdd[npts], xss[npts],xssub[npts], ss, ssub, dd;
	system("mkdir -p stf waveforms waveforms_subevents");
	for (j=0;j<num_sta;j++) {
                for (i=0;i<num_subevent;i++) {
                        templog=(float)sqrt(2*log(10));
                        realsig=modelvec[i].L/modelvec[i].Vr;
                        ssig=realsig/(2*templog);
                        realmu=modelvec[i].ct;
                        genstf(deltat,ssig,realmu,stfnpts,tdata,realsig,i);
                }
	}
	for (i=0;i<num_subevent;i++) {
        	printf("1e27 %.4f %.4f %.4f %.4f %.4f %.4f \n", tempgxm[i*5], tempgxm[i*5+2], tempgxm[i*5+3], tempgxm[i*5+1], tempgxm[i*5+4], -tempgxm[i*5]-tempgxm[i*5+1]);
	}
	float mtimes[num_subevent];
	float submrf[num_subevent][stfnpts], summrf[stfnpts];
	for (i=0;i<num_subevent;i++) {
                mtimes[i]=sqrt(tempgxm[i*5]*tempgxm[i*5]+2*tempgxm[i*5+2]*tempgxm[i*5+2]+2*tempgxm[i*5+3]*tempgxm[i*5+3]+tempgxm[i*5+1]*tempgxm[i*5+1]+2*tempgxm[i*5+4]*tempgxm[i*5+4]+ (-tempgxm[i*5]-tempgxm[i*5+1])*(-tempgxm[i*5]-tempgxm[i*5+1])/2);
		printf("%.3f\n",mtimes[i]);
        }
	for (k=0;k<stfnpts;k++){
                summrf[k]=0;
        }
	for (i=0;i<num_subevent;i++) {
                        templog=(float)sqrt(2*log(10));
                        realsig=modelvec[i].L/modelvec[i].Vr;
                        ssig=realsig/(2*templog);
                        realmu=modelvec[i].ct;
                        genstf1(stfnpts, deltat,ssig,realmu,stfnpts,tdata,realsig,mtimes[i],submrf[i]);
                        for (k=0;k<stfnpts;k++){
                                summrf[k]=summrf[k]+submrf[i][k];
                        }
        }
	char cc2[]="mrf_sum.sac";
	outputsac(stfnpts, summrf, deltat, cc2);
	FILE *fp1;
        if ((fp1 = fopen("sta_residual.dat", "w+")) == NULL) {
                fprintf(stderr,"Cannot open station residual file.\n");
                exit(1);
        }
        char tempobs[]="obs_";char tempsyn[]="syn_"; char BB[]=".sac"; char tempsub[]="synsub_"; char tempp[]="_";
        char aaa[64],bbb[64],ccc[64];
	for (j=0;j<num_sta;j++) {
                tempresidual=0;
                sprintf(aaa,"%s%04d%s",tempobs,j,BB);
                sprintf(bbb,"%s%04d%s",tempsyn,j,BB);
                for (k=0;k<npts;k++) {
                        ss=0;
                        for (i=0;i<5*num_subevent;i++) {
                                jj=j*npts+k;
                                ss = ss +(double)tempgxm[i]*(double)GG[i*npts*num_sta+jj];
                        }
                        dd=tempdx[j*npts+k];
                        xdd[k]=dd; xss[k]=ss;
                }
		for (l=0;l<num_subevent;l++) {
			sprintf(ccc,"%s%04d%s%d%s",tempsub,j,tempp,l+1,BB);
			for (k=0;k<npts;k++) {
				ssub=0;
				for (i=0;i<5;i++) {
					ii=l*5+i;
					jj=j*npts+k;
					ssub=ssub+(double)tempgxm[ii]*(double)GG[ii*npts*num_sta+jj];
				}
				xssub[k]=ssub;
			}
			outputsac(npts, xssub, deltat,ccc);
		}
		tempresidual=l2res[j]/normresidual*exp(fabs(1-CCC[j]))*num_sta;
                outputsac(npts, xdd, deltat, aaa);
                outputsac(npts, xss, deltat, bbb);
                fprintf(fp1,"%d %f %f %f\n",j+1,tempresidual, CCC[j], T_Cross[j]);
        }
	fclose(fp1);
	system("rm stf/*.sac waveforms/*.sac waveforms_subevents/*.sac");
	system("mv stf_*.sac mrf*.sac stf");
	system("mv obs_*.sac syn_*.sac waveforms");
	system("mv synsub_*.sac waveforms_subevents");
	for (j=0;j<num_subevent;j++) {
		printf("clvd[%d] = %f\n", j+1, rclvd[j]);
	}
	printf("dataresidual= %f\nfinalresidual= %f\n",residual/normresidual,finalresidual);
}


void output_forward_joint(int num_subevent, int num_stap, int num_stash, int num_starayl, modelformat modelvec[maxstanum], int nptsp, int nptssh, int nptsrayl, int stfnpts, float deltat, float *tdata, double *tempgxm, double *rclvd, double *GG, double *tempdx, float *pCCC, float *shCCC, float *raylCCC, float *pT_Cross, float *shT_Cross, float *raylT_Cross, float *l2resp, float *l2ressh, float *l2resrayl, float residualp, float residualsh, float residualrayl, float normresidualp, float normresidualsh, float normresidualrayl, float finalresidual) {
        int i,j,k,l,ii,jj,kk;
        float templog, realsig, ssig, realmu,tempresidual, xddp[nptsp], xssp[nptsp], xddsh[nptssh], xsssh[nptssh], xssubp[nptsp], xssubsh[nptssh], ss, ssub, dd;
	float xddrayl[nptsrayl], xssrayl[nptsrayl], xssubrayl[nptsrayl];
        system("mkdir -p stf waveforms waveforms_subevents");
        for (i=0;i<num_subevent;i++) {
                templog=(float)sqrt(2*log(10));
                realsig=modelvec[i].L/modelvec[i].Vr;
                ssig=realsig/(2*templog);
                realmu=modelvec[i].ct;
                genstf(deltat,ssig,realmu,stfnpts,tdata,realsig,i);
        }
	float summt[6];
	for (i=0;i<6;i++) {summt[i]=0;};
	FILE *file = fopen("fm.dat", "w");
        if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
        }
        for (i=0;i<num_subevent;i++) {
                printf("1e27 %.4f %.4f %.4f %.4f %.4f %.4f \n", tempgxm[i*5], tempgxm[i*5+2], tempgxm[i*5+3], tempgxm[i*5+1], tempgxm[i*5+4], -tempgxm[i*5]-tempgxm[i*5+1]);
		fprintf(file, "1e27 %.4f %.4f %.4f %.4f %.4f %.4f\n", tempgxm[i*5], tempgxm[i*5+2], tempgxm[i*5+3], tempgxm[i*5+1], tempgxm[i*5+4], -tempgxm[i*5]-tempgxm[i*5+1]);
		summt[0]+=tempgxm[i*5];
		summt[1]+=tempgxm[i*5+2];
		summt[2]+=tempgxm[i*5+3];
		summt[3]+=tempgxm[i*5+1];
		summt[4]+=tempgxm[i*5+4];
		summt[5]+=-tempgxm[i*5]-tempgxm[i*5+1];
	}
	fclose(file);
	printf("\n1e27 %.4f %.4f %.4f %.4f %.4f %.4f \n",summt[0],summt[1],summt[2],summt[3],summt[4],summt[5]);
	float mtimes[num_subevent];
        float submrf[num_subevent][stfnpts], summrf[stfnpts], summrfsquare[stfnpts],diffsummrf[stfnpts];
        for (i=0;i<num_subevent;i++) {
                mtimes[i]=sqrt((tempgxm[i*5]*tempgxm[i*5]+2*tempgxm[i*5+2]*tempgxm[i*5+2]+2*tempgxm[i*5+3]*tempgxm[i*5+3]+tempgxm[i*5+1]*tempgxm[i*5+1]+2*tempgxm[i*5+4]*tempgxm[i*5+4]+ (-tempgxm[i*5]-tempgxm[i*5+1])*(-tempgxm[i*5]-tempgxm[i*5+1]))/2);
		printf("%.3f\n",mtimes[i]);
        }
        for (k=0;k<stfnpts;k++){
                summrf[k]=0;
		summrfsquare[k]=0;
        }
        for (i=0;i<num_subevent;i++) {
                        templog=(float)sqrt(2*log(10));
                        realsig=modelvec[i].L/modelvec[i].Vr;
                        ssig=realsig/(2*templog);
                        realmu=modelvec[i].ct;
                        genstf1(stfnpts, deltat,ssig,realmu,stfnpts,tdata,realsig,mtimes[i],submrf[i]);
                        for (k=0;k<stfnpts;k++){
                                summrf[k]=summrf[k]+submrf[i][k];
                        }
        }
	diffsummrf[0]=0;
	for (k=1;k<stfnpts;k++){
		diffsummrf[k]=(summrf[k]-summrf[k-1])/deltat;
	}
	for (k=0;k<stfnpts;k++){
		summrfsquare[k]=diffsummrf[k]*diffsummrf[k];
	}
        char cc2[]="mrf_sum.sac"; char cc3[]="mrf_sum_sq.sac"; char cc4[]="diffmrf_sum.sac";
        outputsac(stfnpts, summrf, deltat, cc2);
	outputsac(stfnpts, summrfsquare, deltat, cc3);
	outputsac(stfnpts, diffsummrf, deltat, cc4);
        FILE *fp1;
        if ((fp1 = fopen("sta_residual.dat", "w+")) == NULL) {
                fprintf(stderr,"Cannot open station residual file.\n");
                exit(1);
        }
	fprintf(fp1,"%s\n","# P waves");
        char tempobsp[]="P_obs_", tempsynp[]="P_syn_", tempobssh[]="SH_obs_", tempsynsh[]="SH_syn_", BB[]=".sac", tempsubp[]="P_synsub_", tempsubsh[]="SH_synsub_", tempp[]="_";
	char tempobsrayle[]="rayle_obs_", tempobsrayln[]="rayln_obs_", tempobsraylz[]="raylz_obs_", tempsynrayle[]="rayle_syn_", tempsynrayln[]="rayln_syn_", tempsynraylz[]="raylz_syn_", tempsubrayle[]="rayle_synsub_", tempsubrayln[]="rayln_synsub_", tempsubraylz[]="raylz_synsub_";
        char aaa[64],bbb[64],ccc[64];
        for (j=0;j<num_stap;j++) {
                tempresidual=0;
                sprintf(aaa,"%s%04d%s",tempobsp,j,BB);
                sprintf(bbb,"%s%04d%s",tempsynp,j,BB);
                for (k=0;k<nptsp;k++) {
                        ss=0;
                        for (i=0;i<5*num_subevent;i++) {
                                jj=j*nptsp+k;
                                ss = ss +(double)tempgxm[i]*(double)GG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj];
                        }
                        dd=tempdx[j*nptsp+k];
                        xddp[k]=dd; xssp[k]=ss;
                }
                for (l=0;l<num_subevent;l++) {
                        sprintf(ccc,"%s%04d%s%d%s",tempsubp,j,tempp,l+1,BB);
                        for (k=0;k<nptsp;k++) {
                                ssub=0;
                                for (i=0;i<5;i++) {
                                        ii=l*5+i;
                                        jj=j*nptsp+k;
                                        ssub=ssub+(double)tempgxm[ii]*(double)GG[ii*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj];
                                }
                                xssubp[k]=ssub;
                        }
                        outputsac(nptsp, xssubp, deltat,ccc);
                }
                tempresidual=l2resp[j]/normresidualp*exp(fabs(1-pCCC[j]))*num_stap;
                outputsac(nptsp, xddp, deltat, aaa);
                outputsac(nptsp, xssp, deltat, bbb);
                fprintf(fp1,"%d %f %f %f\n",j+1,tempresidual, pCCC[j], pT_Cross[j]);
        }
        fprintf(fp1,"%s\n","# SH waves");
	for (j=0;j<num_stash;j++) {
                tempresidual=0;
                sprintf(aaa,"%s%04d%s",tempobssh,j,BB);
                sprintf(bbb,"%s%04d%s",tempsynsh,j,BB);
                for (k=0;k<nptssh;k++) {
                        ss=0;
                        for (i=0;i<5*num_subevent;i++) {
                                jj=j*nptssh+nptsp*num_stap+k;
                                ss = ss +(double)tempgxm[i]*(double)GG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj];
                        }
                        dd=tempdx[j*nptssh+nptsp*num_stap+k];
                        xddsh[k]=dd; xsssh[k]=ss;
                }
                for (l=0;l<num_subevent;l++) {
                        sprintf(ccc,"%s%04d%s%d%s",tempsubsh,j,tempp,l+1,BB);
                        for (k=0;k<nptssh;k++) {
                                ssub=0;
                                for (i=0;i<5;i++) {
                                        ii=l*5+i;
                                        jj=j*nptssh+nptsp*num_stap+k;
                                        ssub=ssub+(double)tempgxm[ii]*(double)GG[ii*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj];
                                }
                                xssubsh[k]=ssub;
                        }
                        outputsac(nptssh, xssubsh, deltat,ccc);
                }
                tempresidual=l2ressh[j]/normresidualsh*exp(fabs(1-shCCC[j]))*num_stash;
                outputsac(nptssh, xddsh, deltat, aaa);
                outputsac(nptssh, xsssh, deltat, bbb);
                fprintf(fp1,"%d %f %f %f\n",j+1,tempresidual, shCCC[j], shT_Cross[j]);
        }
	fprintf(fp1,"%s\n","# Rayleigh waves");
	for (j=0;j<3*num_starayl;j++) {
                tempresidual=0;
		jj=j%3;
                kk=(int)(j/3);
		if (jj==0) {
			sprintf(aaa,"%s%04d%s",tempobsrayle,kk,BB);
			sprintf(bbb,"%s%04d%s",tempsynrayle,kk,BB);
                } else if (jj==1) {
			sprintf(aaa,"%s%04d%s",tempobsrayln,kk,BB);
                        sprintf(bbb,"%s%04d%s",tempsynrayln,kk,BB);
                } else if (jj==2) {
			sprintf(aaa,"%s%04d%s",tempobsraylz,kk,BB);
                        sprintf(bbb,"%s%04d%s",tempsynraylz,kk,BB);
                }
                for (k=0;k<nptsrayl;k++) {
                        ss=0;
                        for (i=0;i<5*num_subevent;i++) {
                                ss = ss +(double)tempgxm[i]*(double)GG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+j*nptsrayl+k+nptsp*num_stap+nptssh*num_stash];
                        }
                        dd=tempdx[j*nptsrayl+k+nptsp*num_stap+nptssh*num_stash];
                        xddrayl[k]=dd; xssrayl[k]=ss;
                }
                for (l=0;l<num_subevent;l++) {
			if (jj==0) {
				sprintf(ccc,"%s%04d%s%d%s",tempsubrayle,kk,tempp,l+1,BB);
	                } else if (jj==1) {
				sprintf(ccc,"%s%04d%s%d%s",tempsubrayln,kk,tempp,l+1,BB);
	                } else if (jj==2) {
				sprintf(ccc,"%s%04d%s%d%s",tempsubraylz,kk,tempp,l+1,BB);
	                }	
                        for (k=0;k<nptsrayl;k++) {
                                ssub=0;
                                for (i=0;i<5;i++) {
                                        ii=l*5+i;
                                        ssub=ssub+(double)tempgxm[ii]*(double)GG[ii*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+j*nptsrayl+k+nptsp*num_stap+nptssh*num_stash];
                                }
				xssubrayl[k]=ssub;
                        }
                        outputsac(nptsrayl, xssubrayl, deltat,ccc);
                }
                tempresidual=l2resrayl[j]/normresidualrayl*exp(fabs(1-raylCCC[j]))*num_starayl*3;
                outputsac(nptsrayl, xddrayl, deltat, aaa);
                outputsac(nptsrayl, xssrayl, deltat, bbb);
                fprintf(fp1,"%d %f %f %f\n",j+1,tempresidual, raylCCC[j], raylT_Cross[j]);
        }

        fclose(fp1);
        system("rm stf/*.sac waveforms/*.sac waveforms_subevents/*.sac");
        system("mv *mrf*sac stf_*.sac stf");
        system("mv *_obs_*.sac *_syn_*.sac waveforms");
        system("mv *_synsub_*.sac waveforms_subevents");
        for (j=0;j<num_subevent;j++) {
                printf("clvd[%d] = %f\n", j+1, rclvd[j]);
        }
	//printf("residualp= %f  residualsh= %f residualrayl= %f\n",residualp/normresidualp, residualsh/normresidualsh, residualrayl/normresidualrayl);
	printf("residualp= %f  residualsh= %f residualrayl= %f\n",residualp/(normresidualp+normresidualsh+normresidualrayl), residualsh/(normresidualp+normresidualsh+normresidualrayl), residualrayl/(normresidualp+normresidualsh+normresidualrayl));
        printf("dataresidual= %f  finalresidual= %f\n",(residualp+residualsh+residualrayl)/(normresidualp+normresidualsh+normresidualrayl),finalresidual);
}
