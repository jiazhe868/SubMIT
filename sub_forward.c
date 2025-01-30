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
#include "r8lib.h"
#include "sacsubc.h"
#include "sub_header.h"
#define DEG2RAD (M_PI / 180.0)

void apply(float    *data,
      int       nsamps,
      int       zp,
      float    *sn,
      float    *sd,
      int       nsects);
extern void mtdcmp_(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

void read_modelvec(float *modellist, int num_subevent, modelformat modelvec[maxstanum]) {
        int i,tempi;
        for (i=0;i<num_subevent;i++) {
                tempi=i*7;
                modelvec[i].ct=modellist[tempi];
		modelvec[i].cx=modellist[tempi+1];
		modelvec[i].cy=modellist[tempi+2];
                modelvec[i].L=modellist[tempi+3];
                modelvec[i].Vr=modellist[tempi+4];
                modelvec[i].theta_rupture=modellist[tempi+5];
		modelvec[i].dep=modellist[tempi+6];
        }
}

void convolve(float* kernel, float* kernels, float* in1, float* in2, float* out1, float* out2, int kernel_length, int length)
{
        int i,k;
	float *pout1, *pout2, *pin1, *pin2, *pkernel1, *pkernel2;
        for(i=0; i<length+kernel_length-1; i++){
	pout1=out1+i;
	pout2=out2+i;
        *pout1 = 0.;
	*pout2 = 0.;
        for(k=0; k<kernel_length; k++){
	pin1=in1+i+1-kernel_length+k;
	pin2=in2+i+1-kernel_length+k;
	pkernel1=kernel+kernel_length-k-1;
	pkernel2=kernels+kernel_length-k-1;
        if (i+1-kernel_length+k>=0&&i+1-kernel_length+k<length) {
		*pout1 += (*pin1++) * (*pkernel1++);
		*pout2 += (*pin2++) * (*pkernel2++);
	} else {
		pin1++;
		pin2++;
		pkernel1++;
		pkernel2++;
	}
        }
	}
	pout1=NULL;
	pout2=NULL;
	pin1=NULL;
	pin2=NULL;
	pkernel1=NULL;
	pkernel2=NULL;
}

void simpleconvolve(float* kernel, float* in1, float* out1, int kernel_length, int length){
        int i,k;
        float *pout1, *pin1, *pkernel1;
        for(i=0; i<length+kernel_length-1; i++){
        pout1=out1+i;
        *pout1 = 0.;
        for(k=0; k<kernel_length; k++){
        pin1=in1+i+1-kernel_length+k;
        pkernel1=kernel+kernel_length-k-1;
        if (i+1-kernel_length+k>=0&&i+1-kernel_length+k<length) {
                *pout1 += (*pin1++) * (*pkernel1++);
        } else {
                pin1++;
                pkernel1++;
        }
        }
        }
        pout1=NULL;
        pin1=NULL;
        pkernel1=NULL;
}

void genGaussP(int gnpts, float gdeltat, float gsig, float gmu, float gsigs, float gmus, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float *grnddpz, float *grnddsz, float *grndspz, float *grndssz, float *grnsspz, float *grnsssz, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, int nn, float dura) {
        int k;
        float variance2 = gsig*gsig*2.0;
        float variance2s = gsigs*gsigs*2.0;
        float gaussfun[stfnpts], gaussfuns[stfnpts], grntemp1[stfnpts+gnpts-1],grntemp2[stfnpts+gnpts-1],grntemp3[stfnpts+gnpts-1], grntemp4[stfnpts+gnpts-1],grntemp5[stfnpts+gnpts-1],grntemp6[stfnpts+gnpts-1];
        int bgpt=(int)((gmu-dura/1.5-stfbg)/gdeltat);
        int ndpt=(int)((gmu+dura/1.5-stfbg)/gdeltat);
        if (bgpt<0) {bgpt=0;}
        if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
        for (k=0;k<bgpt;k++){
                gaussfun[k]=0;
                gaussfuns[k]=0;
        }
        for (k=bgpt;k<ndpt;k++){
                gaussfun[k]=1*gdeltat/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2);
                gaussfuns[k]=1*gdeltat/(sqrt(2*3.1416)*gsigs)*exp(-((tdata1[k]-gmus)*(tdata1[k]-gmus))/variance2s);
        }
        for (k=ndpt;k<stfnpts;k++){
                gaussfun[k]=0;
                gaussfuns[k]=0;
        }
        float newbtime=databg+stfbg, newetime=datand+stfnd;
        int newnpt=(int)((newetime-newbtime)/gdeltat+1);
        int ib=(int)((databg-newbtime)/gdeltat+1), ie=(int)((datand-newbtime)/gdeltat+1);
        convolve(gaussfun, gaussfuns, grnddpz, grnddsz, grntemp1, grntemp2, stfnpts, gnpts);
        convolve(gaussfun, gaussfuns, grndspz, grndssz, grntemp3, grntemp4, stfnpts, gnpts);
        convolve(gaussfun, gaussfuns, grnsspz, grnsssz, grntemp5, grntemp6, stfnpts, gnpts);
        for (k=ib-1; k<ie; k++) {
                g1[k-ib+1]=-(grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2;
                g2[k-ib+1]=(grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2;
                g3[k-ib+1]=-(grntemp5[k]+grntemp6[k])*sin(az*2/rad2deg);
                g4[k-ib+1]=-(grntemp3[k]+grntemp4[k])*cos(az/rad2deg);
                g5[k-ib+1]=-(grntemp3[k]+grntemp4[k])*sin(az/rad2deg);
        }
        apply(g1,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g2,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g3,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g4,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g5,(long int) gnpts,1,flt_sn,flt_sd, nsects);
}

void genGaussP_interp(int gnpts, float gdeltat, float gsig, float gmu, float gsigs, float gmus, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float dz, float *grnddpz0, float *grnddpz1, float *grnddsz0, float *grnddsz1, float *grndspz0, float *grndspz1, float *grndssz0, float *grndssz1, float *grnsspz0, float *grnsspz1, float *grnsssz0, float *grnsssz1, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, int nn, float dura) {
        int k;
        float variance2 = gsig*gsig*2.0;
        float variance2s = gsigs*gsigs*2.0;
        float gaussfun[stfnpts], gaussfuns[stfnpts], grntemp1[stfnpts+gnpts-1],grntemp2[stfnpts+gnpts-1],grntemp3[stfnpts+gnpts-1], grntemp4[stfnpts+gnpts-1],grntemp5[stfnpts+gnpts-1],grntemp6[stfnpts+gnpts-1], grnddpz[gnpts], grnddsz[gnpts], grndspz[gnpts], grndssz[gnpts], grnsspz[gnpts], grnsssz[gnpts];
        for (k=0;k<gnpts;k++) {
                grnddpz[k]=grnddpz0[k]*(1-dz)+grnddpz1[k]*dz;
                grnddsz[k]=grnddsz0[k]*(1-dz)+grnddsz1[k]*dz;
                grndspz[k]=grndspz0[k]*(1-dz)+grndspz1[k]*dz;
                grndssz[k]=grndssz0[k]*(1-dz)+grndssz1[k]*dz;
                grnsspz[k]=grnsspz0[k]*(1-dz)+grnsspz1[k]*dz;
                grnsssz[k]=grnsssz0[k]*(1-dz)+grnsssz1[k]*dz;
        }
        int bgpt=(int)((gmu-dura/1.5-stfbg)/gdeltat);
        int ndpt=(int)((gmu+dura/1.5-stfbg)/gdeltat);
        if (bgpt<0) {bgpt=0;}
        if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
        for (k=0;k<bgpt;k++){
                gaussfun[k]=0;
                gaussfuns[k]=0;
        }
        for (k=bgpt;k<ndpt;k++){
                gaussfun[k]=1*gdeltat/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2);
                gaussfuns[k]=1*gdeltat/(sqrt(2*3.1416)*gsigs)*exp(-((tdata1[k]-gmus)*(tdata1[k]-gmus))/variance2s);
        }
        for (k=ndpt;k<stfnpts;k++){
                gaussfun[k]=0;
                gaussfuns[k]=0;
        }
        float newbtime=databg+stfbg, newetime=datand+stfnd;
        int newnpt=(int)((newetime-newbtime)/gdeltat+1);
        int ib=(int)((databg-newbtime)/gdeltat+1), ie=(int)((datand-newbtime)/gdeltat+1);
        convolve(gaussfun, gaussfuns, grnddpz, grnddsz, grntemp1, grntemp2, stfnpts, gnpts);
        convolve(gaussfun, gaussfuns, grndspz, grndssz, grntemp3, grntemp4, stfnpts, gnpts);
        convolve(gaussfun, gaussfuns, grnsspz, grnsssz, grntemp5, grntemp6, stfnpts, gnpts);
        for (k=ib-1; k<ie; k++) {
                g1[k-ib+1]=-(grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2;
                g2[k-ib+1]=(grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2;
                g3[k-ib+1]=-(grntemp5[k]+grntemp6[k])*sin(az*2/rad2deg);
                g4[k-ib+1]=-(grntemp3[k]+grntemp4[k])*cos(az/rad2deg);
                g5[k-ib+1]=-(grntemp3[k]+grntemp4[k])*sin(az/rad2deg);
        }
        apply(g1,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g2,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g3,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g4,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g5,(long int) gnpts,1,flt_sn,flt_sd, nsects);
}

void genGaussSH(int gnpts, float gdeltat, float gsig, float gmu, float gsigs, float gmus, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float *grndsst, float *grnssst, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, int nn, float dura) {
        int k;
        float variance2 = gsig*gsig*2.0;
	float variance2s = gsigs*gsigs*2.0;
        float gaussfun[stfnpts], gaussfuns[stfnpts], grntemp1[stfnpts+gnpts-1],grntemp2[stfnpts+gnpts-1];
	int bgpt=(int)((gmu-dura/1.5-stfbg)/gdeltat);
	int ndpt=(int)((gmu+dura/1.5-stfbg)/gdeltat);
	if (bgpt<0) {bgpt=0;}
	if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
	for (k=0;k<bgpt;k++){
                gaussfun[k]=0;
                gaussfuns[k]=0;
        }
        for (k=bgpt;k<ndpt;k++){
                gaussfun[k]=1*gdeltat/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2);
		gaussfuns[k]=1*gdeltat/(sqrt(2*3.1416)*gsigs)*exp(-((tdata1[k]-gmus)*(tdata1[k]-gmus))/variance2s);
        }
	for (k=ndpt;k<stfnpts;k++){
		gaussfun[k]=0;
		gaussfuns[k]=0;
	}
	float newbtime=databg+stfbg, newetime=datand+stfnd;
        int newnpt=(int)((newetime-newbtime)/gdeltat+1);
        int ib=(int)((databg-newbtime)/gdeltat+1), ie=(int)((datand-newbtime)/gdeltat+1);
	convolve(gaussfuns, gaussfuns, grndsst, grnssst, grntemp1, grntemp2, stfnpts, gnpts);
	for (k=ib-1; k<ie; k++) {
		g1[k-ib+1]=-grntemp2[k]/2*sin(az*2/rad2deg);
		g2[k-ib+1]=grntemp2[k]/2*sin(az*2/rad2deg);
		g3[k-ib+1]=grntemp2[k]*cos(az*2/rad2deg);
		g4[k-ib+1]=-grntemp1[k]*sin(az/rad2deg);
		g5[k-ib+1]=grntemp1[k]*cos(az/rad2deg);
	}
        apply(g1,(long int) gnpts,1,flt_sn,flt_sd, nsects);
	apply(g2,(long int) gnpts,1,flt_sn,flt_sd, nsects);
	apply(g3,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g4,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g5,(long int) gnpts,1,flt_sn,flt_sd, nsects);	
}

void genGaussSH_interp(int gnpts, float gdeltat, float gsig, float gmu, float gsigs, float gmus, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5,  float dz, float *grndsst0, float *grndsst1, float *grnssst0, float *grnssst1, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, int nn, float dura) {
        int k;
        float variance2 = gsig*gsig*2.0;
        float variance2s = gsigs*gsigs*2.0;
        float gaussfun[stfnpts], gaussfuns[stfnpts], grntemp1[stfnpts+gnpts-1],grntemp2[stfnpts+gnpts-1], grndsst[gnpts], grnssst[gnpts];
        for (k=0;k<gnpts;k++) {
                grndsst[k]=grndsst0[k]*(1-dz)+grndsst1[k]*dz;
                grnssst[k]=grnssst0[k]*(1-dz)+grnssst1[k]*dz;
        }
        int bgpt=(int)((gmu-dura/1.5-stfbg)/gdeltat);
        int ndpt=(int)((gmu+dura/1.5-stfbg)/gdeltat);
        if (bgpt<0) {bgpt=0;}
        if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
        for (k=0;k<bgpt;k++){
                gaussfun[k]=0;
                gaussfuns[k]=0;
        }
        for (k=bgpt;k<ndpt;k++){
                gaussfun[k]=1*gdeltat/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2);
                gaussfuns[k]=1*gdeltat/(sqrt(2*3.1416)*gsigs)*exp(-((tdata1[k]-gmus)*(tdata1[k]-gmus))/variance2s);
        }
        for (k=ndpt;k<stfnpts;k++){
                gaussfun[k]=0;
                gaussfuns[k]=0;
        }
        float newbtime=databg+stfbg, newetime=datand+stfnd;
        int newnpt=(int)((newetime-newbtime)/gdeltat+1);
        int ib=(int)((databg-newbtime)/gdeltat+1), ie=(int)((datand-newbtime)/gdeltat+1);
        convolve(gaussfuns, gaussfuns, grndsst, grnssst, grntemp1, grntemp2, stfnpts, gnpts);
        for (k=ib-1; k<ie; k++) {
                g1[k-ib+1]=-grntemp2[k]/2*sin(az*2/rad2deg);
                g2[k-ib+1]=grntemp2[k]/2*sin(az*2/rad2deg);
                g3[k-ib+1]=grntemp2[k]*cos(az*2/rad2deg);
                g4[k-ib+1]=-grntemp1[k]*sin(az/rad2deg);
                g5[k-ib+1]=grntemp1[k]*cos(az/rad2deg);
        }
        apply(g1,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g2,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g3,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g4,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g5,(long int) gnpts,1,flt_sn,flt_sd, nsects);
}

void genGaussRayle_interp(int gnpts, float gdeltat, float gsig, float gmu, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float dx, float dz, float *grnddpr000, float *grnddpr010, float *grnddpr001, float *grnddpr011, float *grnddsr000, float *grnddsr010, float *grnddsr001, float *grnddsr011, float *grndspr000, float *grndspr010, float *grndspr001, float *grndspr011, float *grndssr000, float *grndssr010, float *grndssr001, float *grndssr011, float *grnsspr000, float *grnsspr010, float *grnsspr001, float *grnsspr011, float *grnsssr000, float *grnsssr010, float *grnsssr001, float *grnsssr011, float *grndsst000, float *grndsst010, float *grndsst001, float *grndsst011, float *grnssst000, float *grnssst010, float *grnssst001, float *grnssst011, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, float wt, float dura, int nn) {
	int k;
        float variance2 = gsig*gsig*2.0;
        float gaussfun[stfnpts], grntemp1[stfnpts+gnpts-1],grntemp2[stfnpts+gnpts-1],grntemp3[stfnpts+gnpts-1], grntemp4[stfnpts+gnpts-1],grntemp5[stfnpts+gnpts-1],grntemp6[stfnpts+gnpts-1],grntemp7[stfnpts+gnpts-1],grntemp8[stfnpts+gnpts-1], grnddpr[gnpts], grnddsr[gnpts], grndspr[gnpts], grndssr[gnpts], grnsspr[gnpts], grnsssr[gnpts],grndsst[gnpts],grnssst[gnpts];
        for (k=0;k<gnpts;k++) {
                grnddpr[k]=grnddpr000[k]*(1-dx)*(1-dz)+grnddpr010[k]*dx*(1-dz)+grnddpr001[k]*(1-dx)*dz+grnddpr011[k]*dx*dz;
		grnddsr[k]=grnddsr000[k]*(1-dx)*(1-dz)+grnddsr010[k]*dx*(1-dz)+grnddsr001[k]*(1-dx)*dz+grnddsr011[k]*dx*dz;
		grndspr[k]=grndspr000[k]*(1-dx)*(1-dz)+grndspr010[k]*dx*(1-dz)+grndspr001[k]*(1-dx)*dz+grndspr011[k]*dx*dz;
		grndssr[k]=grndssr000[k]*(1-dx)*(1-dz)+grndssr010[k]*dx*(1-dz)+grndssr001[k]*(1-dx)*dz+grndssr011[k]*dx*dz;
		grnsspr[k]=grnsspr000[k]*(1-dx)*(1-dz)+grnsspr010[k]*dx*(1-dz)+grnsspr001[k]*(1-dx)*dz+grnsspr011[k]*dx*dz;
		grnsssr[k]=grnsssr000[k]*(1-dx)*(1-dz)+grnsssr010[k]*dx*(1-dz)+grnsssr001[k]*(1-dx)*dz+grnsssr011[k]*dx*dz;
		grndsst[k]=grndsst000[k]*(1-dx)*(1-dz)+grndsst010[k]*dx*(1-dz)+grndsst001[k]*(1-dx)*dz+grndsst011[k]*dx*dz;
		grnssst[k]=grnssst000[k]*(1-dx)*(1-dz)+grnssst010[k]*dx*(1-dz)+grnssst001[k]*(1-dx)*dz+grnssst011[k]*dx*dz;
		grnddpr[k]=grnddpr[k]*wt;
                grnddsr[k]=grnddsr[k]*wt;
		grndspr[k]=grndspr[k]*wt;
		grndssr[k]=grndssr[k]*wt;
		grnsspr[k]=grnsspr[k]*wt;
		grnsssr[k]=grnsssr[k]*wt;
		grndsst[k]=grndsst[k]*wt;
		grnssst[k]=grnssst[k]*wt;
        }
        int bgpt=(int)((gmu-dura/1.5-stfbg)/gdeltat);
        int ndpt=(int)((gmu+dura/1.5-stfbg)/gdeltat);
        if (bgpt<0) {bgpt=0;}
        if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
        for (k=0;k<bgpt;k++){
                gaussfun[k]=0;
        }
        for (k=bgpt;k<ndpt;k++){
                gaussfun[k]=1*gdeltat/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2);
        }
        for (k=ndpt;k<stfnpts;k++){
                gaussfun[k]=0;
        }
        float newbtime=databg+stfbg, newetime=datand+stfnd;
        int newnpt=(int)((newetime-newbtime)/gdeltat+1);
        int ib=(int)((databg-newbtime)/gdeltat+1), ie=(int)((datand-newbtime)/gdeltat+1);
	convolve(gaussfun, gaussfun, grnddpr, grnddsr, grntemp1, grntemp2, stfnpts, gnpts);
        convolve(gaussfun, gaussfun, grndspr, grndssr, grntemp3, grntemp4, stfnpts, gnpts);
        convolve(gaussfun, gaussfun, grnsspr, grnsssr, grntemp5, grntemp6, stfnpts, gnpts);
        convolve(gaussfun, gaussfun, grndsst, grnssst, grntemp7, grntemp8, stfnpts, gnpts);
        for (k=ib-1; k<ie; k++) {
		g1[k-ib+1]=(-(grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2)*sin(az/rad2deg)+(-grntemp8[k]/2*sin(az*2/rad2deg))*cos(az/rad2deg);
                g2[k-ib+1]=((grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2)*sin(az/rad2deg)+(grntemp8[k]/2*sin(az*2/rad2deg))*cos(az/rad2deg);
                g3[k-ib+1]=(-(grntemp5[k]+grntemp6[k])*sin(az*2/rad2deg))*sin(az/rad2deg)+(grntemp8[k]*cos(az*2/rad2deg))*cos(az/rad2deg);
                g4[k-ib+1]=(-(grntemp3[k]+grntemp4[k])*cos(az/rad2deg))*sin(az/rad2deg)+(-grntemp7[k]*sin(az/rad2deg))*cos(az/rad2deg);
                g5[k-ib+1]=(-(grntemp3[k]+grntemp4[k])*sin(az/rad2deg))*sin(az/rad2deg)+(grntemp7[k]*cos(az/rad2deg))*cos(az/rad2deg);
        }
	apply(g1,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g2,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g3,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g4,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g5,(long int) gnpts,1,flt_sn,flt_sd, nsects);
}

void genGaussRayln_interp(int gnpts, float gdeltat, float gsig, float gmu, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float dx, float dz, float *grnddpr000, float *grnddpr010, float *grnddpr001, float *grnddpr011, float *grnddsr000, float *grnddsr010, float *grnddsr001, float *grnddsr011, float *grndspr000, float *grndspr010, float *grndspr001, float *grndspr011, float *grndssr000, float *grndssr010, float *grndssr001, float *grndssr011, float *grnsspr000, float *grnsspr010, float *grnsspr001, float *grnsspr011, float *grnsssr000, float *grnsssr010, float *grnsssr001, float *grnsssr011, float *grndsst000, float *grndsst010, float *grndsst001, float *grndsst011, float *grnssst000, float *grnssst010, float *grnssst001, float *grnssst011, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, float wt, float dura, int nn) {
        int k;
        float variance2 = gsig*gsig*2.0;
        float gaussfun[stfnpts], grntemp1[stfnpts+gnpts-1],grntemp2[stfnpts+gnpts-1],grntemp3[stfnpts+gnpts-1], grntemp4[stfnpts+gnpts-1],grntemp5[stfnpts+gnpts-1],grntemp6[stfnpts+gnpts-1],grntemp7[stfnpts+gnpts-1],grntemp8[stfnpts+gnpts-1], grnddpr[gnpts], grnddsr[gnpts], grndspr[gnpts], grndssr[gnpts], grnsspr[gnpts], grnsssr[gnpts],grndsst[gnpts],grnssst[gnpts];
        for (k=0;k<gnpts;k++) {
		grnddpr[k]=grnddpr000[k]*(1-dx)*(1-dz)+grnddpr010[k]*dx*(1-dz)+grnddpr001[k]*(1-dx)*dz+grnddpr011[k]*dx*dz;
                grnddsr[k]=grnddsr000[k]*(1-dx)*(1-dz)+grnddsr010[k]*dx*(1-dz)+grnddsr001[k]*(1-dx)*dz+grnddsr011[k]*dx*dz;
                grndspr[k]=grndspr000[k]*(1-dx)*(1-dz)+grndspr010[k]*dx*(1-dz)+grndspr001[k]*(1-dx)*dz+grndspr011[k]*dx*dz;
                grndssr[k]=grndssr000[k]*(1-dx)*(1-dz)+grndssr010[k]*dx*(1-dz)+grndssr001[k]*(1-dx)*dz+grndssr011[k]*dx*dz;
                grnsspr[k]=grnsspr000[k]*(1-dx)*(1-dz)+grnsspr010[k]*dx*(1-dz)+grnsspr001[k]*(1-dx)*dz+grnsspr011[k]*dx*dz;
                grnsssr[k]=grnsssr000[k]*(1-dx)*(1-dz)+grnsssr010[k]*dx*(1-dz)+grnsssr001[k]*(1-dx)*dz+grnsssr011[k]*dx*dz;
                grndsst[k]=grndsst000[k]*(1-dx)*(1-dz)+grndsst010[k]*dx*(1-dz)+grndsst001[k]*(1-dx)*dz+grndsst011[k]*dx*dz;
                grnssst[k]=grnssst000[k]*(1-dx)*(1-dz)+grnssst010[k]*dx*(1-dz)+grnssst001[k]*(1-dx)*dz+grnssst011[k]*dx*dz;
		grnddpr[k]=grnddpr[k]*wt;
                grnddsr[k]=grnddsr[k]*wt;
                grndspr[k]=grndspr[k]*wt;
                grndssr[k]=grndssr[k]*wt;
                grnsspr[k]=grnsspr[k]*wt;
                grnsssr[k]=grnsssr[k]*wt;
                grndsst[k]=grndsst[k]*wt;
                grnssst[k]=grnssst[k]*wt;
        }
        int bgpt=(int)((gmu-dura/1.5-stfbg)/gdeltat);
        int ndpt=(int)((gmu+dura/1.5-stfbg)/gdeltat);
        if (bgpt<0) {bgpt=0;}
        if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
        for (k=0;k<bgpt;k++){
                gaussfun[k]=0;
        }
        for (k=bgpt;k<ndpt;k++){
                gaussfun[k]=1*gdeltat/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2);
        }
        for (k=ndpt;k<stfnpts;k++){
                gaussfun[k]=0;
        }
        float newbtime=databg+stfbg, newetime=datand+stfnd;
        int newnpt=(int)((newetime-newbtime)/gdeltat+1);
        int ib=(int)((databg-newbtime)/gdeltat+1), ie=(int)((datand-newbtime)/gdeltat+1);
        convolve(gaussfun, gaussfun, grnddpr, grnddsr, grntemp1, grntemp2, stfnpts, gnpts);
        convolve(gaussfun, gaussfun, grndspr, grndssr, grntemp3, grntemp4, stfnpts, gnpts);
        convolve(gaussfun, gaussfun, grnsspr, grnsssr, grntemp5, grntemp6, stfnpts, gnpts);
        convolve(gaussfun, gaussfun, grndsst, grnssst, grntemp7, grntemp8, stfnpts, gnpts);
        for (k=ib-1; k<ie; k++) {
                g1[k-ib+1]=(-(grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2)*cos(az/rad2deg)-(-grntemp8[k]/2*sin(az*2/rad2deg))*sin(az/rad2deg);
                g2[k-ib+1]=((grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2)*cos(az/rad2deg)-(grntemp8[k]/2*sin(az*2/rad2deg))*sin(az/rad2deg);
                g3[k-ib+1]=(-(grntemp5[k]+grntemp6[k])*sin(az*2/rad2deg))*cos(az/rad2deg)-(grntemp8[k]*cos(az*2/rad2deg))*sin(az/rad2deg);
		g4[k-ib+1]=(-(grntemp3[k]+grntemp4[k])*cos(az/rad2deg))*cos(az/rad2deg)-(-grntemp7[k]*sin(az/rad2deg))*sin(az/rad2deg);
                g5[k-ib+1]=(-(grntemp3[k]+grntemp4[k])*sin(az/rad2deg))*cos(az/rad2deg)-(grntemp7[k]*cos(az/rad2deg))*sin(az/rad2deg);
        }
        apply(g1,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g2,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g3,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g4,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g5,(long int) gnpts,1,flt_sn,flt_sd, nsects);
}

void genGaussRaylz_interp(int gnpts, float gdeltat, float gsig, float gmu, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float dx, float dz, float *grnddpz000, float *grnddpz010, float *grnddpz001, float *grnddpz011, float *grnddsz000, float *grnddsz010, float *grnddsz001, float *grnddsz011, float *grndspz000, float *grndspz010, float *grndspz001, float *grndspz011, float *grndssz000, float *grndssz010, float *grndssz001, float *grndssz011, float *grnsspz000, float *grnsspz010, float *grnsspz001, float *grnsspz011, float *grnsssz000, float *grnsssz010, float *grnsssz001, float *grnsssz011, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, float wt, float dura, int nn) {
        int k;
        float variance2 = gsig*gsig*2.0;
        float gaussfun[stfnpts], grntemp1[stfnpts+gnpts-1],grntemp2[stfnpts+gnpts-1],grntemp3[stfnpts+gnpts-1], grntemp4[stfnpts+gnpts-1],grntemp5[stfnpts+gnpts-1],grntemp6[stfnpts+gnpts-1], grnddpz[gnpts], grnddsz[gnpts], grndspz[gnpts], grndssz[gnpts], grnsspz[gnpts],grnsssz[gnpts];
        for (k=0;k<gnpts;k++) {
                grnddpz[k]=grnddpz000[k]*(1-dx)*(1-dz)+grnddpz010[k]*dx*(1-dz)+grnddpz001[k]*(1-dx)*dz+grnddpz011[k]*dx*dz;
                grnddsz[k]=grnddsz000[k]*(1-dx)*(1-dz)+grnddsz010[k]*dx*(1-dz)+grnddsz001[k]*(1-dx)*dz+grnddsz011[k]*dx*dz;
                grndspz[k]=grndspz000[k]*(1-dx)*(1-dz)+grndspz010[k]*dx*(1-dz)+grndspz001[k]*(1-dx)*dz+grndspz011[k]*dx*dz;
                grndssz[k]=grndssz000[k]*(1-dx)*(1-dz)+grndssz010[k]*dx*(1-dz)+grndssz001[k]*(1-dx)*dz+grndssz011[k]*dx*dz;
                grnsspz[k]=grnsspz000[k]*(1-dx)*(1-dz)+grnsspz010[k]*dx*(1-dz)+grnsspz001[k]*(1-dx)*dz+grnsspz011[k]*dx*dz;
                grnsssz[k]=grnsssz000[k]*(1-dx)*(1-dz)+grnsssz010[k]*dx*(1-dz)+grnsssz001[k]*(1-dx)*dz+grnsssz011[k]*dx*dz;
		grnddpz[k]=grnddpz[k]*wt;
                grnddsz[k]=grnddsz[k]*wt;
                grndspz[k]=grndspz[k]*wt;
                grndssz[k]=grndssz[k]*wt;
                grnsspz[k]=grnsspz[k]*wt;
                grnsssz[k]=grnsssz[k]*wt;
        }
        int bgpt=(int)((gmu-dura/1.5-stfbg)/gdeltat);
        int ndpt=(int)((gmu+dura/1.5-stfbg)/gdeltat);
        if (bgpt<0) {bgpt=0;}
        if (ndpt>stfnpts-1){ndpt=stfnpts-1;}
        for (k=0;k<bgpt;k++){
                gaussfun[k]=0;
        }
        for (k=bgpt;k<ndpt;k++){
                gaussfun[k]=1*gdeltat/(sqrt(2*3.1416)*gsig)*exp(-((tdata1[k]-gmu)*(tdata1[k]-gmu))/variance2);
        }
        for (k=ndpt;k<stfnpts;k++){
                gaussfun[k]=0;
        }
        float newbtime=databg+stfbg, newetime=datand+stfnd;
        int newnpt=(int)((newetime-newbtime)/gdeltat+1);
        int ib=(int)((databg-newbtime)/gdeltat+1), ie=(int)((datand-newbtime)/gdeltat+1);
        convolve(gaussfun, gaussfun, grnddpz, grnddsz, grntemp1, grntemp2, stfnpts, gnpts);
        convolve(gaussfun, gaussfun, grndspz, grndssz, grntemp3, grntemp4, stfnpts, gnpts);
        convolve(gaussfun, gaussfun, grnsspz, grnsssz, grntemp5, grntemp6, stfnpts, gnpts);
        for (k=ib-1; k<ie; k++) {
                g1[k-ib+1]=(-(grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2);
                g2[k-ib+1]=((grntemp5[k]+grntemp6[k])/2*cos(az*2/rad2deg)-(grntemp1[k]+grntemp2[k])/2);
                g3[k-ib+1]=(-(grntemp5[k]+grntemp6[k])*sin(az*2/rad2deg));
                g4[k-ib+1]=(-(grntemp3[k]+grntemp4[k])*cos(az/rad2deg));
                g5[k-ib+1]=(-(grntemp3[k]+grntemp4[k])*sin(az/rad2deg));
        }
        apply(g1,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g2,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g3,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g4,(long int) gnpts,1,flt_sn,flt_sd, nsects);
        apply(g5,(long int) gnpts,1,flt_sn,flt_sd, nsects);
}

void Search_Cross_Correlation(float *U0,float *U1,int N0,float Delta0,int N1,float T_shift,float dT_Search,float *T_Cross,float *Cross_Corre_Coefficient,float *Scale_Factor,float *Rms_Factor)
        {int dI;
         if((dI=(int)(dT_Search/Delta0))==0) {
                dI=1;
         }
         int I_Cross,I,J;
         float Cov_U0_U1,Cov_U0_U0,Cov_U1_U1;
         float R,R_max,R1,R2;
         float Corre_FACTOR,RMS_FACTOR;
         Cov_U0_U0=0.0;
         for(J=0;J<N0;J++)
               {Cov_U0_U0+=(*(U0+J))*(*(U0+J));
                }
         R_max=-1.0;
         for(I=0;I<N1;I=I+dI){
            Cov_U0_U1=0.0;
            Cov_U1_U1=0.0;
            R1=R;
            R=-1.0;
            for(J=I;J<I+N0;J++)
               {Cov_U0_U1+=(*(U0+J-I))*(*(U1+J));
                Cov_U1_U1+=(*(U1+J))*(*(U1+J));}
            R=Cov_U0_U1/sqrt((Cov_U0_U0*Cov_U1_U1));
            if(R>R_max)
               {R_max=R;
                I_Cross=I;
                Corre_FACTOR=Cov_U0_U1/Cov_U0_U0;
		RMS_FACTOR=sqrt(fabs(Cov_U1_U1)/(fabs(Cov_U0_U0)+1e-10));}

            if(I>dI&&R1>R&&R1>R2)
            R2=R1;
            }
         *Cross_Corre_Coefficient=R_max;
         *T_Cross=-T_shift+I_Cross*Delta0;
         *Scale_Factor=Corre_FACTOR;
         *Rms_Factor= RMS_FACTOR;
}

void calculate_p_t_axes(double strike, double dip, double rake, double p_axis[3], double t_axis[3]) {
    //strike *= DEG2RAD;
    //dip *= DEG2RAD;
    //rake *= DEG2RAD;
    double n[3] = {
        -sin(dip) * sin(strike),
        sin(dip) * cos(strike),
        -cos(dip)
    };
    double d[3] = {
        cos(rake) * cos(strike) + sin(rake) * cos(dip) * sin(strike),
        cos(rake) * sin(strike) - sin(rake) * cos(dip) * cos(strike),
        -sin(rake) * sin(dip)
    };
    for (int i = 0; i < 3; ++i) {
        p_axis[i] = (n[i] + d[i]) / 2;
    }
    double norm_p = sqrt(p_axis[0] * p_axis[0] + p_axis[1] * p_axis[1] + p_axis[2] * p_axis[2]);
    if (norm_p==0) {
            norm_p=norm_p+(1e-6);
    }
    for (int i = 0; i < 3; ++i) {
        p_axis[i] /= norm_p;
    }
    for (int i = 0; i < 3; ++i) {
        t_axis[i] = (n[i] - d[i]) / 2;
    }
    double norm_t = sqrt(t_axis[0] * t_axis[0] + t_axis[1] * t_axis[1] + t_axis[2] * t_axis[2]);
    if (norm_t==0) {
            norm_t=norm_t+(1e-6);
    }
    for (int i = 0; i < 3; ++i) {
        t_axis[i] /= norm_t;
    }
}

double calculate_angle(double axis1[3], double axis2[3]) {
    double dot_product = axis1[0] * axis2[0] + axis1[1] * axis2[1] + axis1[2] * axis2[2];
    if (dot_product>=1) {
	    dot_product=dot_product-(1e-6);
    }
    if (dot_product<=-1) {
            dot_product=dot_product+(1e-6);
    }
    double angle = acos(dot_product);
    angle = angle / DEG2RAD;
    if (angle > 90.0) {
        angle = 180.0 - angle;
    }
    return angle;
}

void sub_forward_(float *btimep, float *etimep, float *btimesh, float *etimesh, float *btimerayl, float *etimerayl, float *stfbtime, float *stfetime, int *ptsp, int *ptssh, int *ptsrayl, float *delta, int *nstap, int *nstash, int *nstarayl, int *nsubev, float *initmodel, char stainfop_stname[maxstanum][64], float *stainfop_stdist, float *stainfop_staz, float *stainfop_stcp, float *stainfop_stcs, char stainfosh_stname[maxstanum][64], float *stainfosh_stdist, float *stainfosh_staz, float *stainfosh_stcp, float *stainfosh_stcs, float *stainforayl_stlo, float *stainforayl_stla, float stainforayl_wt[maxlocnum][3], char seisdatap_stname[maxstanum][64], float *seisdatap_btime, float *seisdatap_t1, float seisdatap_xdata[maxstanum][maxnpts], char seisdatash_stname[maxstanum][64], float *seisdatash_btime, float *seisdatash_t1, float seisdatash_xdata[maxstanum][maxnpts], char seisdatarayl_stname[maxlocnum][64], float *seisdatarayl_btime, float *seisdatarayl_t1, float seisdatarayl_xdata[maxlocnum*3][maxnpts], float greendatap[maxndepth*maxstanum*6][maxnpts], float greendatash[maxndepth*maxstanum*2][maxnpts], float greendatarayl[14*maxndepth*maxdistnum][maxnpts], float *residualall, float *fltsnp, float *fltsdp, float *fltsnsh, float *fltsdsh, float *fltsnrayl, float *fltsdrayl, int *nsectp, int *nsectsh, int *nsectrayl, float *wtp, float *wtsh, float *wtrayl, float *CMTsolution, float *alphat, float *cmtscaling, float *dc, float *evlo, float *evla, float *evvp, float *evvs, float *begin_localdist, float *end_localdist, float *interval_localdist, int *n_distloc, float *begin_dep, float *end_dep, float *interval_dep, int *n_dep) {
	int i,j,k,l,snsectp, snsectsh, snsectrayl, ii,jj,kk, ii000,ii010,ii100,ii110, ii001,ii011,ii101,ii111;
	int num_subevent=*nsubev, ndist=*n_distloc;
	float wp=*wtp, wsh=*wtsh, wrayl=*wtrayl, alpha=*alphat, CMTscaling=*cmtscaling, elo=*evlo, ela=*evla, evp=*evvp, evs=*evvs, bg_dist=*begin_localdist, nd_dist=*end_localdist, interval_dist=*interval_localdist;
	double dcconstrain=(double)*dc;
	float sfltsnp[maxflen], sfltsdp[maxflen], sfltsnsh[maxflen], sfltsdsh[maxflen];
	float sfltsnrayl[maxflen], sfltsdrayl[maxflen];
	float bg_dep,nd_dep,dep_interval;
	int ndepth;
	char stastnamep[maxstanum][64], seisstnamep[maxstanum][64], stastnamesh[maxstanum][64], seisstnamesh[maxstanum][64];
	char seisstnamerayl[maxlocnum][64];
        float pstdist[maxstanum], pstaz[maxstanum], pstcp[maxstanum], pstcs[maxstanum], pbtime[maxstanum], pt1[maxstanum], shstdist[maxstanum], shstaz[maxstanum], shstcp[maxstanum], shstcs[maxstanum], shbtime[maxstanum], sht1[maxstanum], CMT[5], raylstwt[maxlocnum][3];
	float raylstlo[maxlocnum], raylstla[maxlocnum], raylbtime[maxlocnum], raylt1[maxlocnum];
	float stf_btime=*stfbtime;
        float stf_etime=*stfetime;
        float deltat=*delta;
        int stfnpts=(int)((stf_etime-stf_btime)/deltat+1);
        bg_dep=*begin_dep;
        nd_dep=*end_dep;
        dep_interval=*interval_dep;
        ndepth=*n_dep;
	float realmu,realsig,ssig;
	int indexdp[num_subevent],currentindex[num_subevent],indexdist,currentinddist,N0p,N1p,N0sh,N1sh,N0rayl,N1rayl;
	float tdata[stfnpts];
        for (i=0;i<stfnpts;i++) {
                tdata[i]=stf_btime+i*deltat;
        }
	float residualp=0, residualsh=0, residualrayl=0, normresidualp=0, normresidualsh=0, normresidualrayl=0, dd=0,ss=0,dataresidual=0, finalresidual=0, tempsig, tempsigs, templog, tempmu1, tempmu2, temprclvd=0;
	float resp[maxstanum],ressh[maxstanum],normresp1[maxstanum],normressh1[maxstanum],normresp2[maxstanum],normressh2[maxstanum];
	float resrayl[maxlocnum*3],normresrayl1[maxlocnum*3],normresrayl2[maxlocnum*3];
	*residualall=0;
	staformat stainfop[maxstanum], stainfosh[maxstanum];
	stalocformat stainforayl[maxlocnum];
        modelformat modelvec[maxsubeve];
	seisdataformat* seisdatap;
	seisdataformat* seisdatash;
	seisdatalocformat* seisdatarayl;
	snsectp=*nsectp;
	snsectsh=*nsectsh;
	snsectrayl=*nsectrayl;
	for (i=0;i<5;i++) {
                CMT[i]=CMTsolution[i];
        }
	float initmodel_rev[num_subevent*7];
        for (i=0;i<num_subevent;i++) {
		initmodel_rev[i*7+0]=initmodel[i*7+0];
                initmodel_rev[i*7+1]=initmodel[i*7+1];
                initmodel_rev[i*7+2]=initmodel[i*7+2];
                initmodel_rev[i*7+3]=initmodel[i*7+3]*initmodel[i*7+4];
                initmodel_rev[i*7+4]=initmodel[i*7+4];
                initmodel_rev[i*7+5]=initmodel[i*7+5];
                initmodel_rev[i*7+6]=initmodel[i*7+6];
        }
	read_modelvec(initmodel_rev, num_subevent, modelvec);
	double *ptempG, *ptempd, *ptempdx, tempgxm[num_subevent*5];
	double *rayltempG, *rayltempd, *rayltempdx;
	float **pmatrixG, **shmatrixG, **raylematrixG, **raylnmatrixG, **raylzmatrixG;
	double *pGG, *raylGG, GtG[num_subevent*5*num_subevent*5], AtA[num_subevent*5*num_subevent*5], GtG_AtA[num_subevent*5*num_subevent*5], *GtG_AtA_inv_Gt;
	double rclvd[num_subevent], mt[num_subevent], strike[num_subevent][2],dip[num_subevent][2],rake[num_subevent][2];
        double mt0,mtxx,mtxy,mtxz,mtyy,mtyz,mtzz, rclvd_total, strike_total[2], dip_total[2], rake_total[2];
        double tempsum=0, a1=1, a2=1;
	float maxshftp=2, maxshftsh=5, maxshftrayl=3;
        float Delta0=deltat, dT_Search=deltat; //Delta0: data delta. dT_Search: search step in shifting;
        float pT_shift=maxshftp,shT_shift=maxshftsh, raylT_shift=maxshftrayl; 
	for (i=0;i<num_subevent*5;i++) {
                for (j=0;j<num_subevent*5;j++) {
                        if (i==j){
                                AtA[i*num_subevent*5+j]=alpha*alpha;
                        }
                        else {
                                AtA[i*num_subevent*5+j]=0;
                        }
                }
        }
	if (wp<=0 || wsh<=0 || wrayl<=0) {
                fprintf(stderr,"At least one of weightp and weightsh should be positive! \n");
                exit(1);
        } else  {
		int nptsp=*ptsp, nptssh=*ptssh;
		int nptsrayl=*ptsrayl;
                int num_stap=*nstap, num_stash=*nstash;
                float b_timep=*btimep, b_timesh=*btimesh;
                float e_timep=*etimep, e_timesh=*etimesh;
                float sigp[num_subevent][num_stap], mup[num_subevent][num_stap], sigsp[num_subevent][num_stap], musp[num_subevent][num_stap], l2resp[num_stap];
		float sigsh[num_subevent][num_stash], mush[num_subevent][num_stash], sigssh[num_subevent][num_stash], mussh[num_subevent][num_stash], l2ressh[num_stash];
                float pT_Cross[num_stap], pCross_Corre_Coefficient[num_stap], pScale_Factor[num_stap], pRms_Factor[num_stap]; 
                float shT_Cross[num_stash], shCross_Corre_Coefficient[num_stash], shScale_Factor[num_stash], shRms_Factor[num_stash]; 
		seisdatap=(seisdataformat*)malloc(sizeof(seisdataformat)*maxstanum);
                seisdatash=(seisdataformat*)malloc(sizeof(seisdataformat)*maxstanum);
		int num_starayl=*nstarayl;
                float b_timerayl=*btimerayl;
                float e_timerayl=*etimerayl;
                float sigrayl[num_subevent][num_starayl], murayl[num_subevent][num_starayl], l2resrayl[num_starayl*3];
		float raylT_Cross[num_starayl*3], raylCross_Corre_Coefficient[num_starayl*3], raylScale_Factor[num_starayl*3], raylRms_Factor[num_starayl*3];
                seisdatarayl=(seisdatalocformat*)malloc(sizeof(seisdatalocformat)*maxlocnum);
		for (i=0;i<maxflen;i++) {
                        sfltsnp[i]=fltsnp[i];
                        sfltsdp[i]=fltsdp[i];
			sfltsnsh[i]=fltsnsh[i];
                        sfltsdsh[i]=fltsdsh[i];
			sfltsnrayl[i]=fltsnrayl[i];
                        sfltsdrayl[i]=fltsdrayl[i];
                }
                for (i=0;i<maxstanum;i++) {
                        strcpy(stastnamep[i],stainfop_stname[i]);
                        strcpy(seisstnamep[i],seisdatap_stname[i]);
			strcpy(stastnamesh[i],stainfosh_stname[i]);
                        strcpy(seisstnamesh[i],seisdatash_stname[i]);
                        pstdist[i]=stainfop_stdist[i];
                        pstaz[i]=stainfop_staz[i];
                        pstcp[i]=stainfop_stcp[i];
                        pstcs[i]=stainfop_stcs[i];
                        pbtime[i]=seisdatap_btime[i];
                        pt1[i]=seisdatap_t1[i];
                        shstdist[i]=stainfosh_stdist[i];
                        shstaz[i]=stainfosh_staz[i];
                        shstcp[i]=stainfosh_stcp[i];
                        shstcs[i]=stainfosh_stcs[i];
                        shbtime[i]=seisdatash_btime[i];
                        sht1[i]=seisdatash_t1[i];
                }
		for (i=0;i<maxlocnum;i++) {
			strcpy(seisstnamerayl[i],seisdatarayl_stname[i]);
                        raylstwt[i][0]=stainforayl_wt[i][0];
			raylstwt[i][1]=stainforayl_wt[i][1];
			raylstwt[i][2]=stainforayl_wt[i][2];
			raylstlo[i]=stainforayl_stlo[i];
                        raylstla[i]=stainforayl_stla[i];
                        raylbtime[i]=seisdatarayl_btime[i];
                        raylt1[i]=seisdatarayl_t1[i];
		}
		convert_arrays_staformat(stastnamep, pstdist, pstaz, pstcp, pstcs, stainfop);
                convert_arrays_seisformat(seisstnamep, pbtime, pt1, seisdatap_xdata, seisdatap);
		convert_arrays_staformat(stastnamesh, shstdist, shstaz, shstcp, shstcs, stainfosh);
                convert_arrays_seisformat(seisstnamesh, shbtime, sht1, seisdatash_xdata, seisdatash);
		convert_arrays_stalocformat(raylstlo, raylstla, raylstwt, stainforayl);
		convert_arrays_seislocformat(seisstnamerayl, raylbtime, raylt1, seisdatarayl_xdata, seisdatarayl);
                pmatrixG=alloc2d(num_subevent*5,nptsp);
		shmatrixG=alloc2d(num_subevent*5,nptssh);
		ptempG=(double *)malloc(num_subevent*5*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)*sizeof(double));
                ptempd=(double *)malloc((nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)*sizeof(double));
                ptempdx=(double *)malloc((nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)*sizeof(double));
                pGG=(double *)malloc(num_subevent*5*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)*sizeof(double));
                GtG_AtA_inv_Gt=(double *)malloc(num_subevent*5*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)*sizeof(double));
		raylematrixG=alloc2d(num_subevent*5,nptsrayl);
		raylnmatrixG=alloc2d(num_subevent*5,nptsrayl);
		raylzmatrixG=alloc2d(num_subevent*5,nptsrayl);
		templog=(float)sqrt(2*log(10));
		float x0,y0,x1,y1,tempaz[num_subevent],tempdist[num_subevent],dx[num_subevent],dz[num_subevent],gridx1[num_subevent],gridx2[num_subevent],griddep1[num_subevent],griddep2[num_subevent];
                int indx1[num_subevent],indx2[num_subevent],indz1[num_subevent],indz2[num_subevent];
                for (i=0;i<num_subevent;i++) {
			indz1[i]=(int)((modelvec[i].dep-bg_dep)/dep_interval);
                        indz2[i]=indz1[i]+1;
                        if (indz1[i]<0) {
                                indz1[i]=0; indz2[i]=0;
                        }
                        if (indz1[i]>=ndepth-1) {
                                indz1[i]=ndepth-1; indz2[i]=ndepth-1;
                        }
                        griddep1[i]=((float)indz1[i])*dep_interval+bg_dep;
                        griddep2[i]=((float)indz2[i])*dep_interval+bg_dep;
                        dz[i]=(modelvec[i].dep-griddep1[i])/dep_interval;
                        if (dz[i]<0) {
                                dz[i]=0;
                        }
                        if (dz[i]>1) {
                                dz[i]=1;
                        }
		}
		for (j=0;j<num_stap;j++) {
                        for (i=0;i<num_subevent;i++) {
                                tempsig=modelvec[i].L/modelvec[i].Vr-modelvec[i].L/stainfop[j].stcp*(cos(modelvec[i].theta_rupture/rad2deg)*cos(stainfop[j].staz/rad2deg)+sin(modelvec[i].theta_rupture/rad2deg)*sin(stainfop[j].staz/rad2deg));
                                tempsigs=modelvec[i].L/modelvec[i].Vr-modelvec[i].L/stainfop[j].stcp*(cos(modelvec[i].theta_rupture/rad2deg)*cos(stainfop[j].staz/rad2deg)+sin(modelvec[i].theta_rupture/rad2deg)*sin(stainfop[j].staz/rad2deg));
                                realsig=modelvec[i].L/modelvec[i].Vr;
                                ssig=realsig/(2*templog);
                                realmu=modelvec[i].ct;
                                sigp[i][j]=tempsig/(2*templog);
                                sigsp[i][j]=tempsigs/(2*templog);
                                tempmu1=(float)sqrt(modelvec[i].cx*modelvec[i].cx+modelvec[i].cy*modelvec[i].cy);
                                tempmu2=(float)atan2((double)modelvec[i].cx,(double)modelvec[i].cy);
                                mup[i][j]=modelvec[i].ct-tempmu1/fabs(stainfop[j].stcp)*cos(stainfop[j].staz/rad2deg-tempmu2)-(modelvec[i].dep-modelvec[0].dep)/sqrt(1/(1/(evp*evp)-(1/stainfop[j].stcp)*(1/stainfop[j].stcp)))*stainfop[j].stcp/fabs(stainfop[j].stcp);
                                musp[i][j]=modelvec[i].ct-tempmu1/fabs(stainfop[j].stcp)*cos(stainfop[j].staz/rad2deg-tempmu2)-(modelvec[i].dep-modelvec[0].dep)/sqrt(1/(1/(evs*evs)-(1/stainfop[j].stcs)*(1/stainfop[j].stcs)))*stainfop[j].stcp/fabs(stainfop[j].stcp);
				ii000=6*(indz1[i]*num_stap+j);
                                ii001=6*(indz2[i]*num_stap+j);
				//if (i==0) {printf("%f %f %f %f\n",sigp[i][j],mup[i][j],sigsp[i][j],musp[i][j]);}
				genGaussP_interp(nptsp, deltat, sigp[i][j], mup[i][j], sigsp[i][j], musp[i][j], stfnpts, stf_btime, stf_etime, b_timep, e_timep, tdata, stainfop[j].staz, pmatrixG[5*i], pmatrixG[5*i+1], pmatrixG[5*i+2], pmatrixG[5*i+3], pmatrixG[5*i+4], dz[i], greendatap[ii000], greendatap[ii001], greendatap[ii000+1], greendatap[ii001+1], greendatap[ii000+2], greendatap[ii001+2], greendatap[ii000+3], greendatap[ii001+3], greendatap[ii000+4], greendatap[ii001+4], greendatap[ii000+5], greendatap[ii001+5], sfltsnp, sfltsdp, snsectp, j, tempsig);
                        }
			for (i=0;i<nptsp;i++) {
                                ptempd[j*nptsp+i]=(double)seisdatap[j].xdata[i]*wp;
                        }
                        for (i=0;i<5*num_subevent;i++) {
                                for (k=0;k<nptsp;k++) {
                                        jj=j*nptsp+k;
                                        pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj]=(double)pmatrixG[i][k]*wp;
                                }
                        }
                }
		for (j=0;j<num_stash;j++) {
                        for (i=0;i<num_subevent;i++) {
                                tempsig=modelvec[i].L/modelvec[i].Vr-modelvec[i].L/stainfosh[j].stcs*(cos(modelvec[i].theta_rupture/rad2deg)*cos(stainfosh[j].staz/rad2deg)+sin(modelvec[i].theta_rupture/rad2deg)*sin(stainfosh[j].staz/rad2deg));
                                tempsigs=modelvec[i].L/modelvec[i].Vr-modelvec[i].L/stainfosh[j].stcs*(cos(modelvec[i].theta_rupture/rad2deg)*cos(stainfosh[j].staz/rad2deg)+sin(modelvec[i].theta_rupture/rad2deg)*sin(stainfosh[j].staz/rad2deg));
                                templog=(float)sqrt(2*log(10));
                                realsig=modelvec[i].L/modelvec[i].Vr;
                                ssig=realsig/(2*templog);
                                realmu=modelvec[i].ct;
                                sigsh[i][j]=tempsig/(2*templog);
                                sigssh[i][j]=tempsigs/(2*templog);
                                tempmu1=(float)sqrt(modelvec[i].cx*modelvec[i].cx+modelvec[i].cy*modelvec[i].cy);
                                tempmu2=(float)atan2((double)modelvec[i].cx,(double)modelvec[i].cy);
                                mush[i][j]=modelvec[i].ct-tempmu1/stainfosh[j].stcs*cos(stainfosh[j].staz/rad2deg-tempmu2)-(modelvec[i].dep-modelvec[0].dep)/sqrt(1/(1/(evs*evs)-(1/stainfosh[j].stcs)*(1/stainfosh[j].stcs)))*stainfosh[j].stcs/fabs(stainfosh[j].stcs);
                                mussh[i][j]=modelvec[i].ct-tempmu1/stainfosh[j].stcs*cos(stainfosh[j].staz/rad2deg-tempmu2)-(modelvec[i].dep-modelvec[0].dep)/sqrt(1/(1/(evs*evs)-(1/stainfosh[j].stcs)*(1/stainfosh[j].stcs)))*stainfosh[j].stcs/fabs(stainfosh[j].stcs);
				ii000=2*(indz1[i]*num_stash+j);
                                ii001=2*(indz2[i]*num_stash+j);
				genGaussSH_interp(nptssh, deltat, sigsh[i][j], mush[i][j], sigssh[i][j], mussh[i][j], stfnpts, stf_btime, stf_etime, b_timesh, e_timesh, tdata, stainfosh[j].staz, shmatrixG[5*i], shmatrixG[5*i+1], shmatrixG[5*i+2], shmatrixG[5*i+3], shmatrixG[5*i+4], dz[i], greendatash[ii000], greendatash[ii001], greendatash[ii000+1], greendatash[ii001+1], sfltsnsh, sfltsdsh, snsectsh, j, tempsig);
                        }
                        for (i=0;i<nptssh;i++) {
                                ptempd[nptsp*num_stap+j*nptssh+i]=(double)seisdatash[j].xdata[i]*wsh;
                        }
			for (i=0;i<5*num_subevent;i++) {
                                for (k=0;k<nptssh;k++) {
                                        jj=nptsp*num_stap+j*nptssh+k;
                                        pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj]=(double)shmatrixG[i][k]*wsh;
                                }
                        }
                }
		float newelo,newela,tempa,tempc;
		for (j=0;j<num_starayl;j++) {
                        for (i=0;i<num_subevent;i++) {
				newelo=elo+modelvec[i].cx/(111.32*cos(ela/rad2deg));
				newela=ela+modelvec[i].cy/110.574;
				tempa=pow(sin((stainforayl[j].stla-newela)/rad2deg/2),2)+cos(stainforayl[j].stla/rad2deg)*cos(newela/rad2deg)*pow(sin((stainforayl[j].stlo-newelo)/rad2deg/2),2);
				tempc=2*atan2(sqrt(tempa),sqrt(1-tempa));
				tempdist[i]=6371*tempc;
				tempaz[i]=atan2(sin((stainforayl[j].stlo-newelo)/rad2deg)*cos((stainforayl[j].stla)/rad2deg),cos((newela)/rad2deg)*sin((stainforayl[j].stla)/rad2deg)-sin((newela)/rad2deg)*cos((stainforayl[j].stla)/rad2deg)*cos((stainforayl[j].stlo-newelo)/rad2deg))*rad2deg;
				/*
				x0=(stainforayl[j].stlo-elo)*111.32*cos((stainforayl[j].stla+ela)/2*3.1416/180);
                                y0=(stainforayl[j].stla-ela)*110.574;
                                x1=x0-modelvec[i].cx;
                                y1=y0-modelvec[i].cy;
                                //azrayl[j]=90-atan2(y0,x0)*rad2deg;
                                if (sqrt(x1*x1)>0.01||sqrt(y1*y1)>0.01){
                                tempaz[i]=90-atan2(y1,x1)*rad2deg;} else {tempaz[i]=0;}
                                tempdist[i]=sqrt(x1*x1+y1*y1);
				*/
				indx1[i]=(int)((tempdist[i]-bg_dist)/interval_dist);
				indx2[i]=indx1[i]+1;
				if (indx1[i]<0) {
                        	        indx1[i]=0; indx2[i]=0;
                        	}
				if (indx1[i]>=ndist-1) {
                                        indx1[i]=ndist-1; indx2[i]=ndist-1;
                                }
				gridx1[i]=((float)indx1[i])*interval_dist+bg_dist;
				gridx2[i]=((float)indx2[i])*interval_dist+bg_dist;
				dx[i]=(tempdist[i]-gridx1[i])/interval_dist;
				if (dx[i]<0) {
                        	        dx[i]=0;
                        	}
				if (dx[i]>1) {
                        	        dx[i]=1;
				}
				
                                tempsig=modelvec[i].L/modelvec[i].Vr-modelvec[i].L/(evs*0.85)*(cos(modelvec[i].theta_rupture/rad2deg)*cos(tempaz[i]/rad2deg)+sin(modelvec[i].theta_rupture/rad2deg)*sin(tempaz[i]/rad2deg));
                                realsig=modelvec[i].L/modelvec[i].Vr;
                                ssig=realsig/(2*templog);
                                realmu=modelvec[i].ct;
                                sigrayl[i][j]=tempsig/(2*templog);
                                murayl[i][j]=modelvec[i].ct;
				ii000=14*(indz1[i]*ndist+indx1[i]);
				ii010=14*(indz1[i]*ndist+indx2[i]);
                                ii001=14*(indz2[i]*ndist+indx1[i]);
                                ii011=14*(indz2[i]*ndist+indx2[i]);
                                genGaussRayle_interp(nptsrayl, deltat, sigrayl[i][j], murayl[i][j], stfnpts, stf_btime, stf_etime, b_timerayl, e_timerayl, tdata, tempaz[i], raylematrixG[5*i], raylematrixG[5*i+1], raylematrixG[5*i+2], raylematrixG[5*i+3], raylematrixG[5*i+4], dx[i], dz[i],  greendatarayl[ii000+6], greendatarayl[ii010+6], greendatarayl[ii001+6], greendatarayl[ii011+6], greendatarayl[ii000+7], greendatarayl[ii010+7], greendatarayl[ii001+7], greendatarayl[ii011+7], greendatarayl[ii000+8], greendatarayl[ii010+8], greendatarayl[ii001+8], greendatarayl[ii011+8], greendatarayl[ii000+9], greendatarayl[ii010+9], greendatarayl[ii001+9], greendatarayl[ii011+9], greendatarayl[ii000+10], greendatarayl[ii010+10], greendatarayl[ii001+10], greendatarayl[ii011+10], greendatarayl[ii000+11], greendatarayl[ii010+11], greendatarayl[ii001+11], greendatarayl[ii011+11], greendatarayl[ii000+12], greendatarayl[ii010+12], greendatarayl[ii001+12], greendatarayl[ii011+12], greendatarayl[ii000+13], greendatarayl[ii010+13], greendatarayl[ii001+13], greendatarayl[ii011+13], sfltsnrayl, sfltsdrayl, snsectrayl, stainforayl[j].wte, tempsig, j);
				genGaussRayln_interp(nptsrayl, deltat, sigrayl[i][j], murayl[i][j], stfnpts, stf_btime, stf_etime, b_timerayl, e_timerayl, tdata, tempaz[i], raylnmatrixG[5*i], raylnmatrixG[5*i+1], raylnmatrixG[5*i+2], raylnmatrixG[5*i+3], raylnmatrixG[5*i+4], dx[i], dz[i],  greendatarayl[ii000+6], greendatarayl[ii010+6], greendatarayl[ii001+6], greendatarayl[ii011+6], greendatarayl[ii000+7], greendatarayl[ii010+7], greendatarayl[ii001+7], greendatarayl[ii011+7], greendatarayl[ii000+8], greendatarayl[ii010+8], greendatarayl[ii001+8], greendatarayl[ii011+8], greendatarayl[ii000+9], greendatarayl[ii010+9], greendatarayl[ii001+9], greendatarayl[ii011+9], greendatarayl[ii000+10], greendatarayl[ii010+10], greendatarayl[ii001+10], greendatarayl[ii011+10], greendatarayl[ii000+11], greendatarayl[ii010+11], greendatarayl[ii001+11], greendatarayl[ii011+11], greendatarayl[ii000+12], greendatarayl[ii010+12], greendatarayl[ii001+12], greendatarayl[ii011+12], greendatarayl[ii000+13], greendatarayl[ii010+13], greendatarayl[ii001+13], greendatarayl[ii011+13], sfltsnrayl, sfltsdrayl, snsectrayl, stainforayl[j].wtn, tempsig, j);
				genGaussRaylz_interp(nptsrayl, deltat, sigrayl[i][j], murayl[i][j], stfnpts, stf_btime, stf_etime, b_timerayl, e_timerayl, tdata, tempaz[i], raylzmatrixG[5*i], raylzmatrixG[5*i+1], raylzmatrixG[5*i+2], raylzmatrixG[5*i+3], raylzmatrixG[5*i+4], dx[i], dz[i],  greendatarayl[ii000], greendatarayl[ii010], greendatarayl[ii001], greendatarayl[ii011], greendatarayl[ii000+1], greendatarayl[ii010+1], greendatarayl[ii001+1], greendatarayl[ii011+1], greendatarayl[ii000+2], greendatarayl[ii010+2], greendatarayl[ii001+2], greendatarayl[ii011+2], greendatarayl[ii000+3], greendatarayl[ii010+3], greendatarayl[ii001+3], greendatarayl[ii011+3], greendatarayl[ii000+4], greendatarayl[ii010+4], greendatarayl[ii001+4], greendatarayl[ii011+4], greendatarayl[ii000+5], greendatarayl[ii010+5], greendatarayl[ii001+5], greendatarayl[ii011+5], sfltsnrayl, sfltsdrayl, snsectrayl, stainforayl[j].wtz, tempsig, j);
                        }
                        for (i=0;i<nptsrayl;i++) {
                                ptempd[nptsp*num_stap+nptssh*num_stash+3*j*nptsrayl+i]=(double)seisdatarayl[j].edata[i]*wrayl;
				ptempd[nptsp*num_stap+nptssh*num_stash+(3*j+1)*nptsrayl+i]=(double)seisdatarayl[j].ndata[i]*wrayl;
				ptempd[nptsp*num_stap+nptssh*num_stash+(3*j+2)*nptsrayl+i]=(double)seisdatarayl[j].zdata[i]*wrayl;
                        }
                        for (i=0;i<5*num_subevent;i++) {
                                for (k=0;k<nptsrayl;k++) {
                                        jj=nptsp*num_stap+nptssh*num_stash+3*j*nptsrayl+k;
                                        pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj]=(double)raylematrixG[i][k]*wrayl;
					pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj+nptsrayl]=(double)raylnmatrixG[i][k]*wrayl;
					pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj+2*nptsrayl]=(double)raylzmatrixG[i][k]*wrayl;
                                }
                        }
                }

		for (i=0;i<num_subevent*5;i++) {
                        for (j=0;j<nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3;j++) {
                                jj=i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+j;
                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=pGG[jj];
                        }
                        for (j=nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3;j<nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5;j++) {
                                if (i%5==0) {
                                        if (j==nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3) {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=(double)1*CMTscaling;
                                        }       else {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=0;
                                        }
                                }
                                if (i%5==1) {
                                        if (j==nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+1) {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=(double)1*CMTscaling;
                                        }       else {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=0;
                                        }
                                }
                                if (i%5==2) {
                                        if (j==nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+2) {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=(double)1*CMTscaling;
                                        }       else {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=0;
                                        }
                                }
                                if (i%5==3) {
                                        if (j==nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+3) {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=(double)1*CMTscaling;
                                        }       else {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=0;
                                        }
                                }
				if (i%5==4) {
                                        if (j==nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+4) {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=(double)1*CMTscaling;
                                        }       else {
                                                ptempG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5)+j]=0;
                                        }
                                }
                        }
                }
		for (i=num_stap*nptsp+nptssh*num_stash+nptsrayl*num_starayl*3;i<num_stap*nptsp+nptssh*num_stash+nptsrayl*num_starayl*3+5;i++) {
                        ptempd[i]=(double)(CMT[i-num_stap*nptsp-nptssh*num_stash-nptsrayl*num_starayl*3]*CMTscaling);
                }
                r8mat_mtm_new(num_subevent*5, nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5, num_subevent*5, ptempG, ptempG, GtG);
                r8mat_add(num_subevent*5, num_subevent*5, a1, GtG, a2, AtA, GtG_AtA);
                r8mat_cholesky_inverse(num_subevent*5, GtG_AtA);
                r8mat_mmt_new(num_subevent*5, num_subevent*5, nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5, GtG_AtA, ptempG, GtG_AtA_inv_Gt);
                r8mat_mm(num_subevent*5, nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5, 1, GtG_AtA_inv_Gt, ptempd, tempgxm);
                N0p=nptsp-(int)(maxshftp/deltat);
                N1p=(int)(maxshftp/deltat);
                N0sh=nptssh-(int)(maxshftsh/deltat);
                N1sh=(int)(maxshftsh/deltat);
		N0rayl=nptsrayl-(int)(maxshftrayl/deltat);
                N1rayl=(int)(maxshftrayl/deltat);

		float xddp[nptsp], xdd1p[nptsp], xssp[nptsp+N1p];
		float xddsh[nptssh], xdd1sh[nptssh], xsssh[nptssh+N1sh];
		float xddrayl[nptsrayl], xdd1rayl[nptsrayl], xssrayl[nptsrayl+N1rayl];
		int nptsclip=10;
                for (j=0;j<num_stap;j++) {
                        for (k=0;k<nptsp;k++) {
                                ss=0;
                                for (i=0;i<5*num_subevent;i++) {
                                        jj=j*nptsp+k;
                                        ss = ss +(double)tempgxm[i]*(double)pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj];
                                }
                                dd=seisdatap[j].xdata[k];
                                if (k>=N1p) {xdd1p[k-N1p]=dd;} else {xdd1p[nptsp-1-k]=0;}
                                xddp[k]=dd; xssp[k]=ss; //xss1[k]=ss;
                        }
			for (k=nptsp;k<nptsp+N1p;k++){
				xssp[k]=0;
			}
			pT_Cross[j]=0;
			pCross_Corre_Coefficient[j]=0;
			pScale_Factor[j]=0;
			pRms_Factor[j]=0;
                        Search_Cross_Correlation(xdd1p,xssp,N0p,Delta0,2*N1p,pT_shift,dT_Search, &pT_Cross[j],&pCross_Corre_Coefficient[j],&pScale_Factor[j],&pRms_Factor[j]);
                        if (pT_Cross[j]>0){
                                for (i=0;i<(int)(pT_Cross[j]/deltat);i++) {
                                        ptempdx[j*nptsp+i]=0;
                                }
                                for (i=(int)(pT_Cross[j]/deltat);i<nptsp;i++) {
                                        ptempdx[j*nptsp+i]=ptempd[j*nptsp+i-(int)(pT_Cross[j]/deltat)];
                                }
                        }
                        else {
                                for (i=0;i<nptsp+(int)(pT_Cross[j]/deltat);i++) {
                                        ptempdx[j*nptsp+i]=ptempd[j*nptsp+i-(int)(pT_Cross[j]/deltat)];
                                }
                                for (i=nptsp+(int)(pT_Cross[j]/deltat);i<nptsp;i++) {
                                        ptempdx[j*nptsp+i]=0;
                                }
                        }
                }
		for (j=0;j<num_stash;j++) {
                        for (k=0;k<nptssh;k++) {
                                ss=0;
                                for (i=0;i<5*num_subevent;i++) {
                                        jj=j*nptssh+nptsp*num_stap+k;
                                        ss = ss +(double)tempgxm[i]*(double)pGG[i*(nptssh*num_stash+nptsp*num_stap+nptsrayl*num_starayl*3)+jj];
                                }
                                dd=seisdatash[j].xdata[k];
                                if (k>=N1sh) {xdd1sh[k-N1sh]=dd;} else {xdd1sh[nptssh-1-k]=0;}
                                xddsh[k]=dd; xsssh[k]=ss; //xss1[k]=ss;
                        }
			for (k=nptssh;k<nptssh+N1sh;k++){
                                xsssh[k]=0;
                        }
                        Search_Cross_Correlation(xdd1sh,xsssh,N0sh,Delta0,2*N1sh,shT_shift,dT_Search, &shT_Cross[j],&shCross_Corre_Coefficient[j],&shScale_Factor[j],&shRms_Factor[j]);
                        if (shT_Cross[j]>0){
                                for (i=0;i<(int)(shT_Cross[j]/deltat);i++) {
                                        ptempdx[num_stap*nptsp+j*nptssh+i]=0;
                                }
                                for (i=(int)(shT_Cross[j]/deltat);i<nptssh;i++) {
                                        ptempdx[num_stap*nptsp+j*nptssh+i]=ptempd[num_stap*nptsp+j*nptssh+i-(int)(shT_Cross[j]/deltat)];
                                }
                        }
                        else {
                                for (i=0;i<nptssh+(int)(shT_Cross[j]/deltat);i++) {
                                        ptempdx[num_stap*nptsp+j*nptssh+i]=ptempd[num_stap*nptsp+j*nptssh+i-(int)(shT_Cross[j]/deltat)];
                                }
                                for (i=nptssh+(int)(shT_Cross[j]/deltat);i<nptssh;i++) {
                                        ptempdx[num_stap*nptsp+j*nptssh+i]=0;
                                }
                        }
                }

		for (j=0;j<num_starayl*3;j++) {
			jj=j%3;
			ii=(int)(j/3);
                        for (k=0;k<nptsrayl;k++) {
                                ss=0;
                                for (i=0;i<5*num_subevent;i++) {
                                        ss = ss +(double)tempgxm[i]*(double)pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+j*nptsrayl+nptsp*num_stap+nptssh*num_stash+k];
                                }
				if (jj==0) {
					dd=seisdatarayl[ii].edata[k];
				} else if (jj==1) {
					dd=seisdatarayl[ii].ndata[k];
				} else if (jj==2) {
                                        dd=seisdatarayl[ii].zdata[k];
                                }
                                if (k>=N1rayl) {xdd1rayl[k-N1rayl]=dd;} else {xdd1rayl[nptsrayl-1-k]=0;}
                                xssrayl[k]=ss; 
                        }
                        for (k=nptsrayl;k<nptsrayl+N1rayl;k++){
                                xssrayl[k]=0;
                        }
                        raylT_Cross[j]=0;
                        raylCross_Corre_Coefficient[j]=0;
                        raylScale_Factor[j]=0;
                        raylRms_Factor[j]=0;
                        Search_Cross_Correlation(xdd1rayl,xssrayl,N0rayl,Delta0,2*N1rayl,raylT_shift,dT_Search, &raylT_Cross[j],&raylCross_Corre_Coefficient[j],&raylScale_Factor[j],&raylRms_Factor[j]);
                        if (raylT_Cross[j]>0){
                                for (i=0;i<(int)(raylT_Cross[j]/deltat);i++) {
                                        ptempdx[nptsp*num_stap+nptssh*num_stash+j*nptsrayl+i]=0;
                                }
                                for (i=(int)(raylT_Cross[j]/deltat);i<nptsrayl;i++) {
                                        ptempdx[nptsp*num_stap+nptssh*num_stash+j*nptsrayl+i]=ptempd[nptsp*num_stap+nptssh*num_stash+j*nptsrayl+i-(int)(raylT_Cross[j]/deltat)];
                                }
                        }
                        else {
                                for (i=0;i<nptsrayl+(int)(raylT_Cross[j]/deltat);i++) {
                                        ptempdx[nptsp*num_stap+nptssh*num_stash+j*nptsrayl+i]=ptempd[nptsp*num_stap+nptssh*num_stash+j*nptsrayl+i-(int)(raylT_Cross[j]/deltat)];
                                }
                                for (i=nptsrayl+(int)(raylT_Cross[j]/deltat);i<nptsrayl;i++) {
                                        ptempdx[nptsp*num_stap+nptssh*num_stash+j*nptsrayl+i]=0;
                                }
                        }
                }

		for (i=num_stap*nptsp+nptssh*num_stash+nptsrayl*num_starayl*3;i<num_stap*nptsp+nptssh*num_stash+nptsrayl*num_starayl*3+5;i++) {
                        ptempdx[i]=(double)(CMT[i-num_stap*nptsp-nptssh*num_stash-nptsrayl*num_starayl*3]*CMTscaling);
                }
                r8mat_mm(num_subevent*5, nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3+5, 1, GtG_AtA_inv_Gt, ptempdx, tempgxm);

		for (i=0;i<num_subevent;i++) {
                        mtxx=tempgxm[i*5];
                        mtxy=tempgxm[i*5+2];
                        mtxz=tempgxm[i*5+3];
                        mtyy=tempgxm[i*5+1];
                        mtyz=tempgxm[i*5+4];
                        mtzz=-tempgxm[i*5]-tempgxm[i*5+1];
                        mt0=1e27;
                        mt[i]=sqrt((mtxx*mtxx+mtyy*mtyy+mtzz*mtzz+2*mtxy*mtxy+2*mtxz*mtxz+2*mtyz*mtyz)/2);
                        tempsum=tempsum+mt[i];
                        mtdcmp_(&mt0, &mtxx, &mtxy, &mtxz, &mtyy, &mtyz, &mtzz, &rclvd[i], strike[i], dip[i], rake[i]);
                }
		mtxx=0; mtxy=0; mtxz=0; mtyy=0; mtyz=0; mtzz=0;
		for (i=0;i<num_subevent;i++) {
                        mtxx+=tempgxm[i*5];
                        mtxy+=tempgxm[i*5+2];
                        mtxz+=tempgxm[i*5+3];
                        mtyy+=tempgxm[i*5+1];
                        mtyz+=tempgxm[i*5+4];
                        mtzz+=-tempgxm[i*5]-tempgxm[i*5+1];
		}
		//mtxx=CMT[0];
                //mtxy=CMT[2];
                //mtxz=CMT[3];
       		//mtyy=CMT[1];
                //mtyz=CMT[4];
                //mtzz=-CMT[0]-CMT[1];
		mtdcmp_(&mt0, &mtxx, &mtxy, &mtxz, &mtyy, &mtyz, &mtzz, &rclvd_total, strike_total, dip_total, rake_total);

		for (i=0;i<num_subevent;i++) {
			if (i==0) {temprclvd=temprclvd+(float)(sqrt(rclvd[i])*mt[i]/tempsum);}
                        temprclvd=temprclvd+(float)(sqrt(rclvd[i])*mt[i]/tempsum);
                }

		float dur_ratio[num_subevent],pendur,avg_dur_ratio;
                pendur=0;
                avg_dur_ratio=0;
                for (i=0;i<num_subevent;i++) {
                        dur_ratio[i]=sqrt(mt[i])/(modelvec[i].L/modelvec[i].Vr);
                        avg_dur_ratio+=dur_ratio[i];
                }
                avg_dur_ratio=avg_dur_ratio/num_subevent;
                for (i=0;i<num_subevent;i++) {
                        pendur+=max(dur_ratio[i]/avg_dur_ratio,avg_dur_ratio/dur_ratio[i])-1;
                }
                pendur=pendur/num_subevent;
		int npind[num_subevent];
		double pen_stress;
		double p_axis1[3], p_axis2[3];
		double t_axis1[3], t_axis2[3];
		double p_angle = 0;
	        double t_angle = 0;
		calculate_p_t_axes(strike_total[0], dip_total[0], rake_total[0], p_axis1, t_axis1);
		for (i=0;i<num_subevent;i++) {
			calculate_p_t_axes(strike[i][0], dip[i][0], rake[i][0], p_axis2, t_axis2);
			p_angle = max(p_angle, calculate_angle(p_axis1, p_axis2));
			t_angle = max(t_angle, calculate_angle(t_axis1, t_axis2));
		}
		//printf("%f %f %f %f %f %f\n",p_axis1[0],p_axis1[1],p_axis1[2],t_axis1[0],t_axis1[1],t_axis1[2]);
		pen_stress=max(0,p_angle/30-1)+max(0,t_angle/30-1);
                for (j=0;j<num_stap;j++) {
                        l2resp[j]=0;
			normresp1[j]=0;
			normresp2[j]=0;
                        for (k=0;k<nptsp;k++) {
                                ss=0;
                                for (i=0;i<5*num_subevent;i++) {
                                        jj=j*nptsp+k;
                                        ss = ss +(double)tempgxm[i]*(double)pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj];
                                }
                                dd=seisdatap[j].xdata[k]*wp;
                                if (k>=N1p) {
					xdd1p[k-N1p]=dd;
				} else {
					xdd1p[nptsp-1-k]=0;
				}
                                xddp[k]=dd; xssp[k]=ss;
                        }
			for (k=nptsp;k<nptsp+N1p;k++){
                                xssp[k]=0;
                        }
                        Search_Cross_Correlation(xdd1p,xssp,N0p,Delta0,2*N1p,pT_shift,dT_Search, &pT_Cross[j],&pCross_Corre_Coefficient[j],&pScale_Factor[j],&pRms_Factor[j]);
			for (k=nptsclip;k<nptsp-nptsclip;k++) {
				l2resp[j]=l2resp[j]+(xddp[k-(int)(pT_Cross[j]/deltat)]-xssp[k])*(xddp[k-(int)(pT_Cross[j]/deltat)]-xssp[k]);
                                normresp1[j]+=xddp[k-(int)(pT_Cross[j]/deltat)]*xddp[k-(int)(pT_Cross[j]/deltat)];
				normresp2[j]+=xssp[k]*xssp[k];
                        }	
			normresidualp+=sqrt(normresp1[j]*normresp2[j]);
                        residualp=residualp+l2resp[j]*exp(fabs(1-pCross_Corre_Coefficient[j]));
                }

		for (j=0;j<num_stash;j++) {
                        l2ressh[j]=0;
			normressh1[j]=0;
			normressh2[j]=0;
                        for (k=0;k<nptssh;k++) {
                                ss=0;
                                for (i=0;i<5*num_subevent;i++) {
                                        jj=j*nptssh+nptsp*num_stap+k;
                                        ss = ss +(double)tempgxm[i]*(double)pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+jj];
                                }
                                dd=seisdatash[j].xdata[k]*wsh;
                                if (k>=N1sh) {
					xdd1sh[k-N1sh]=dd;
				} else {
					xdd1sh[nptssh-1-k]=0;
				}
                                xddsh[k]=dd; xsssh[k]=ss; 
                        }
			for (k=nptssh;k<nptssh+N1sh;k++){
                                xsssh[k]=0;
                        }
                        Search_Cross_Correlation(xdd1sh,xsssh,N0sh,Delta0,2*N1sh,shT_shift,dT_Search, &shT_Cross[j],&shCross_Corre_Coefficient[j],&shScale_Factor[j],&shRms_Factor[j]);
			for (k=nptsclip;k<nptssh-nptsclip;k++) {
                                l2ressh[j]=l2ressh[j]+(xddsh[k-(int)(shT_Cross[j]/deltat)]-xsssh[k])*(xddsh[k-(int)(shT_Cross[j]/deltat)]-xsssh[k]);
				normressh1[j]+=(xddsh[k-(int)(shT_Cross[j]/deltat)]*xddsh[k-(int)(shT_Cross[j]/deltat)]);
				normressh2[j]+=(xsssh[k]*xsssh[k]);
                        }
			normresidualsh+=sqrt(normressh1[j]*normressh2[j]);
                        residualsh=residualsh+l2ressh[j]*exp(fabs(1-shCross_Corre_Coefficient[j]));
                }
		for (j=0;j<num_starayl*3;j++) {
                        l2resrayl[j]=0;
                        normresrayl1[j]=0;
                        normresrayl2[j]=0;
			jj=j%3;
                        ii=(int)(j/3);
                        for (k=0;k<nptsrayl;k++) {
                                ss=0;
                                for (i=0;i<5*num_subevent;i++) {
                                        ss = ss +(double)tempgxm[i]*(double)pGG[i*(nptsp*num_stap+nptssh*num_stash+nptsrayl*num_starayl*3)+j*nptsrayl+k+nptsp*num_stap+nptssh*num_stash];
                                }
				if (jj==0) {
                                        dd=seisdatarayl[ii].edata[k]*wrayl;
                                } else if (jj==1) {
                                        dd=seisdatarayl[ii].ndata[k]*wrayl;
                                } else if (jj==2) {
                                        dd=seisdatarayl[ii].zdata[k]*wrayl;
                                }
                                if (k>=N1rayl) {
                                        xdd1rayl[k-N1rayl]=dd;
                                } else {
                                        xdd1rayl[nptsrayl-1-k]=0;
                                }
                                xddrayl[k]=dd; xssrayl[k]=ss;
                        }
                        for (k=nptsrayl;k<nptsrayl+N1rayl;k++){
                                xssrayl[k]=0;
                        }
                        Search_Cross_Correlation(xdd1rayl,xssrayl,N0rayl,Delta0,2*N1rayl,raylT_shift,dT_Search, &raylT_Cross[j],&raylCross_Corre_Coefficient[j],&raylScale_Factor[j],&raylRms_Factor[j]);
                        for (k=nptsclip;k<nptsrayl-nptsclip;k++) {
                                l2resrayl[j]=l2resrayl[j]+(xddrayl[k-(int)(raylT_Cross[j]/deltat)]-xssrayl[k])*(xddrayl[k-(int)(raylT_Cross[j]/deltat)]-xssrayl[k]);
                                normresrayl1[j]+=xddrayl[k-(int)(raylT_Cross[j]/deltat)]*xddrayl[k-(int)(raylT_Cross[j]/deltat)];
                                normresrayl2[j]+=xssrayl[k]*xssrayl[k];
                        }
                        normresidualrayl+=sqrt(normresrayl1[j]*normresrayl2[j]);
                        residualrayl=residualrayl+l2resrayl[j]*exp(fabs(1-raylCross_Corre_Coefficient[j]));
                }
		//residualrayl=0; normresidualrayl=0;
		*residualall=(residualp+residualsh+residualrayl)/(normresidualp+normresidualsh+normresidualrayl);
                dataresidual=*residualall;
		*residualall=*residualall*exp(temprclvd/dcconstrain)*exp(pen_stress);
                finalresidual=*residualall;
		//printf("%f %f %f %f %f %f\n",dataresidual,exp(temprclvd/dcconstrain),exp(pen_stress),finalresidual,p_angle,t_angle);
                //output_forward_joint(num_subevent, num_stap, num_stash, num_starayl, modelvec, nptsp, nptssh, nptsrayl, stfnpts, deltat, tdata, tempgxm, rclvd, pGG, ptempdx, pCross_Corre_Coefficient, shCross_Corre_Coefficient, raylCross_Corre_Coefficient, pT_Cross, shT_Cross, raylT_Cross, l2resp, l2ressh, l2resrayl, residualp, residualsh, residualrayl, normresidualp, normresidualsh, normresidualrayl, finalresidual);
		free(ptempd);
                free(ptempdx);
                free(ptempG);
                free2d(pmatrixG);
                free(seisdatap);
		free2d(shmatrixG);
                free(seisdatash);
		free2d(raylematrixG);
		free2d(raylnmatrixG);
		free2d(raylzmatrixG);
                free(seisdatarayl);
                free(pGG);
                free(GtG_AtA_inv_Gt);
                ptempG=NULL;
                seisdatap=NULL;
                ptempd=NULL;
                ptempdx=NULL;
                pmatrixG=NULL;
                pGG=NULL;
                GtG_AtA_inv_Gt=NULL;
                seisdatash=NULL;
                shmatrixG=NULL;
		seisdatarayl=NULL;
                raylematrixG=NULL;
		raylnmatrixG=NULL;
		raylzmatrixG=NULL;
                GtG_AtA_inv_Gt=NULL;
	}
}
