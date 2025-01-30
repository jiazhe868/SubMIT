/***************************************************************
 * copyright: 2017-, GPS in Caltech
 * File name: sub_init.c
 * Description: Source characterization with sub-events
 * Author:Jiazhe
 * Date:03/28/2017
 * History:03/28/2017 v0.1
 ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sacsubc.h"
#include "sub_header.h"

void read_parameter(float *bg_time_P, float *nd_time_P, float *bg_time_SH, float *nd_time_SH, float *bg_time_Rayl, float *nd_time_Rayl, float *stf_btime, float *stf_etime, float *deltat, int *nptsP, int *nptsSH, int *nptsRayl, int *num_sta_P, int *num_sta_SH,  int *num_sta_Rayl, char stainfofileP[64], char stainfofileSH[64], char stainfofileRayl[64], float *weightP, float *weightSH, float *weightRayl, char inputmodelfile[64], char greenchar1[256], char greenchar2[64], char greenchar3[64], float *evlo, float *evla, float *evvp, float *evvs, float *bg_distloc, float *nd_distloc, float *distloc_interval, float *bg_dep, float *nd_dep, float *dep_interval, int *num_subevent, float *flp, float *fhp, float *fls, float *fhs, float *flr, float *fhr, float *alpha, float *cmtscaling, float *dcconstrain) {
	char junka[512],junkb[512];
        FILE *fp1;
        if ((fp1 = fopen("Par.file", "r")) == NULL) {
                fprintf(stderr,"Cannot open Par file.\n");
                exit(1);
        }

	fscanf(fp1,"%s%f%[^\n]s",junka,bg_time_P,junkb); //bg_time
	if(*bg_time_P==-99999) {
        fprintf(stderr,"Please check your Par.file for bg_time_P!\n");
        exit(1);
	}
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,nd_time_P,junkb); //nd_time
        if(*nd_time_P==-99999) {
        fprintf(stderr,"Please check your Par.file for nd_time_P!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,bg_time_SH,junkb); //bg_time
        if(*bg_time_SH==-99999) {
        fprintf(stderr,"Please check your Par.file for bg_time_SH!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,nd_time_SH,junkb); //nd_time
        if(*nd_time_SH==-99999) {
        fprintf(stderr,"Please check your Par.file for nd_time_SH!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%f%[^\n]s",junka,bg_time_Rayl,junkb); //bg_time
        if(*bg_time_Rayl==-99999) {
        fprintf(stderr,"Please check your Par.file for bg_time_P!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,nd_time_Rayl,junkb); //nd_time
        if(*nd_time_Rayl==-99999) {
        fprintf(stderr,"Please check your Par.file for nd_time_P!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%f%[^\n]s",junka,stf_btime,junkb); //stf bg_time
        if(*stf_btime==-99999) {
        fprintf(stderr,"Please check your Par.file for stf_btime!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);
        fscanf(fp1,"%s%f%[^\n]s",junka,stf_etime,junkb); //stf nd_time
        if(*stf_etime==-99999) {
        fprintf(stderr,"Please check your Par.file for stf_etime!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,deltat,junkb); //delta_t
        if(*deltat==-99999) {
        fprintf(stderr,"Please check your Par.file for deltat!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);
	float temp;
	temp=(*nd_time_P-*bg_time_P)/(*deltat)+1;
	*nptsP=(int)temp;
	temp=(*nd_time_SH-*bg_time_SH)/(*deltat)+1;
        *nptsSH=(int)temp;
	temp=(*nd_time_Rayl-*bg_time_Rayl)/(*deltat)+1;
        *nptsRayl=(int)temp;

        fscanf(fp1,"%s%d%[^\n]s",junka,num_sta_P,junkb); //num_sta
        if(*num_sta_P==-99999) {
        fprintf(stderr,"Please check your Par.file for num_sta_P!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%d%[^\n]s",junka,num_sta_SH,junkb); //num_sta
        if(*num_sta_SH==-99999) {
        fprintf(stderr,"Please check your Par.file for num_sta_SH!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%d%[^\n]s",junka,num_sta_Rayl,junkb); //num_sta
        if(*num_sta_Rayl==-99999) {
        fprintf(stderr,"Please check your Par.file for num_sta_P!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%s%[^\n]s",junka,stainfofileP,junkb); //char stainfofile
        if(stainfofileP[0]==0) {
        fprintf(stderr,"Please check your Par.file for stainfofileP!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%s%[^\n]s",junka,stainfofileSH,junkb); //char stainfofile
        if(stainfofileSH[0]==0) {
        fprintf(stderr,"Please check your Par.file for stainfofileSH!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%s%[^\n]s",junka,stainfofileRayl,junkb); //char stainfofile
        if(stainfofileRayl[0]==0) {
        fprintf(stderr,"Please check your Par.file for stainfofileP!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,weightP,junkb); // fl
        if(*weightP==-99999) {
        fprintf(stderr,"Please check your Par.file for weightP!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,weightSH,junkb); //fh
        if(*weightSH==-99999) {
        fprintf(stderr,"Please check your Par.file for weightSH!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%f%[^\n]s",junka,weightRayl,junkb); // fl
        if(*weightRayl==-99999) {
        fprintf(stderr,"Please check your Par.file for weightP!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%s%[^\n]s",junka,inputmodelfile,junkb); // char inputmodelfile
        if(inputmodelfile[0]==0) {
        fprintf(stderr,"Please check your Par.file for inputmodelfile!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%d%[^\n]s",junka,num_subevent,junkb); //num_subevent
        if(*num_subevent==-99999) {
        fprintf(stderr,"Please check your Par.file for num_subevent!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,flp,junkb); // fl
        if(*flp==-99999) {
        fprintf(stderr,"Please check your Par.file for P Low_frequency!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,fhp,junkb); //fh
        if(*fhp==-99999) {
        fprintf(stderr,"Please check your Par.file for P High_frequency!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,fls,junkb); // fl
        if(*fls==-99999) {
        fprintf(stderr,"Please check your Par.file for S Low_frequency!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,fhs,junkb); //fh
        if(*fhs==-99999) {
        fprintf(stderr,"Please check your Par.file for S High_frequency!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%f%[^\n]s",junka,flr,junkb); // fl
        if(*flr==-99999) {
        fprintf(stderr,"Please check your Par.file for P Low_frequency!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,fhr,junkb); //fh
        if(*fhr==-99999) {
        fprintf(stderr,"Please check your Par.file for P High_frequency!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,alpha,junkb); //fh
        if(*alpha==-99999) {
        fprintf(stderr,"Please check your Par.file for Alpha!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,cmtscaling,junkb); //fh
        if(*cmtscaling==-99999) {
        fprintf(stderr,"Please check your Par.file for CMTscaling!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%f%[^\n]s",junka,evlo,junkb); //fh
        if(*evlo==-99999) {
        fprintf(stderr,"Please check your Par.file for evlo!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,evla,junkb); //fh
        if(*evla==-99999) {
        fprintf(stderr,"Please check your Par.file for evla!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%f%[^\n]s",junka,evvp,junkb); //fh
        if(*evvp==-99999) {
        fprintf(stderr,"Please check your Par.file for evvp!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,evvs,junkb); //fh
        if(*evvs==-99999) {
        fprintf(stderr,"Please check your Par.file for evvs!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%f%[^\n]s",junka,dcconstrain,junkb); //fh
        if(*dcconstrain==-99999) {
        fprintf(stderr,"Please check your Par.file for DC_constrain!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fscanf(fp1,"%s%f%[^\n]s",junka,bg_distloc,junkb); //fh
        if(*bg_distloc==-99999) {
        fprintf(stderr,"Please check your Par.file for bg_localdist!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,nd_distloc,junkb); //fh
        if(*nd_distloc==-99999) {
        fprintf(stderr,"Please check your Par.file for nd_localdist!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,distloc_interval,junkb); //fh
        if(*distloc_interval==-99999) {
        fprintf(stderr,"Please check your Par.file for localdist_interval!\n");
        exit(1);
        }
	fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,bg_dep,junkb); //fh
        if(*bg_dep==-99999) {
        fprintf(stderr,"Please check your Par.file for bg_depth!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,nd_dep,junkb); //fh
        if(*nd_dep==-99999) {
        fprintf(stderr,"Please check your Par.file for nd_depth!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%f%[^\n]s",junka,dep_interval,junkb); //fh
        if(*dep_interval==-99999) {
        fprintf(stderr,"Please check your Par.file for depth_interval!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%s%[^\n]s",junka,greenchar1,junkb); // char inputmodelfile
        if(greenchar1[0]==0) {
        fprintf(stderr,"Please check your Par.file for greenchar1!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%s%[^\n]s",junka,greenchar2,junkb); // char inputmodelfile
        if(greenchar2[0]==0) {
        fprintf(stderr,"Please check your Par.file for greenchar2!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

        fscanf(fp1,"%s%s%[^\n]s",junka,greenchar3,junkb); // char inputmodelfile
        if(greenchar3[0]==0) {
        fprintf(stderr,"Please check your Par.file for greenchar3!\n");
        exit(1);
        }
        fscanf(fp1,"%c",junka);

	fclose(fp1);
}

void read_stainfo(staformat stainfo[maxstanum], char stainfofile[], int num_sta) {
	FILE *fp2;
	int i;
        if ((fp2 = fopen(stainfofile, "r")) == NULL) {
                fprintf(stderr,"Cannot open station info file.\n");
                exit(1);
        }
	for (i=0;i<num_sta;i++) {
		fscanf(fp2,"%s %f %f %f %f %s",stainfo[i].stname,&stainfo[i].stdist,&stainfo[i].staz,&stainfo[i].stcp, &stainfo[i].stcs, stainfo[i].greendist); //input: stainfofile; output: staformat stainfo
	}
	fclose(fp2);
}

void read_stalocinfo(stalocformat stainfo[maxlocnum], char stainfofile[], int num_sta) {
        FILE *fp2;
        int i;
	char tempe[]="e"; char tempn[]="n"; char tempz[]="z"; char tempname[64];
        if ((fp2 = fopen(stainfofile, "r")) == NULL) {
                fprintf(stderr,"Cannot open station info file.\n");
                exit(1);
        }
        for (i=0;i<num_sta;i++) {
                fscanf(fp2,"%s %f %f %s %f %f %f",tempname, &stainfo[i].stlo, &stainfo[i].stla, stainfo[i].stind, &stainfo[i].wte, &stainfo[i].wtn, &stainfo[i].wtz); 
		sprintf(stainfo[i].stename,"%s%s",tempname,tempe);
		sprintf(stainfo[i].stnname,"%s%s",tempname,tempn);
		sprintf(stainfo[i].stzname,"%s%s",tempname,tempz);
        }
        fclose(fp2);
}

void read_initialmodel(modelformat modelvec[maxstanum], float *modellist, char inputmodelfile[], int num_subevent) {
        FILE *fp3;
        if ((fp3 = fopen(inputmodelfile, "r")) == NULL) {
                fprintf(stderr,"Cannot open input model file %s. \n",inputmodelfile);
                exit(1);
        }
	int i,tempi;
        for (i=0;i<num_subevent;i++) {
		tempi=i*7;
		fscanf(fp3,"%f %f %f %f %f %f %f", &modellist[tempi], &modellist[tempi+1],&modellist[tempi+2],&modellist[tempi+3], &modellist[tempi+4], &modellist[tempi+5], &modellist[tempi+6]);
		modelvec[i].ct=modellist[tempi];
                modelvec[i].cx=modellist[tempi+1];
                modelvec[i].cy=modellist[tempi+2];
                modelvec[i].L=modellist[tempi+3];
                modelvec[i].Vr=modellist[tempi+4];
                modelvec[i].theta_rupture=modellist[tempi+5];
                modelvec[i].dep=modellist[tempi+6];
        }
        fclose(fp3);
}


void read_greenpath_P(greenpathformatP greenpath[maxndepth*maxstanum], staformat stainfo[maxstanum], char greenchar1[], char greenchar2[], char greenchar3[], float bg_dep, float nd_dep, float dep_interval, int *ndepth, int num_sta){
        int i,j,k,ii;
        float currentdepth,cdd;
        char aaa[128],bbb1[512],bbb2[512],bbb3[512],bbb4[512],bbb5[512],bbb6[512],tempdep[128];
        char tempddpz[]="ddpz",tempddsz[]="ddsz",tempdspz[]="dspz",tempdssz[]="dssz",tempsspz[]="sspz",tempsssz[]="sssz";
        float ndd=(nd_dep-bg_dep)/dep_interval*100;
        if (fabs(fmod(ndd,100))>0.01) {
                fprintf(stderr,"Depth range and interval are not correct.\n");
                exit(1);
        }
        *ndepth=(int)((nd_dep-bg_dep)/dep_interval+1);
        for (j=0;j<num_sta;j++) {
                strcpy(aaa,stainfo[j].greendist);
                for (i=0;i<*ndepth;i++) {
                        ii=i*num_sta+j;
                        currentdepth=bg_dep+i*dep_interval;
                        cdd=currentdepth*100+0.01;
                        if (fmod(currentdepth,1)<0.01) {
                                sprintf(tempdep,"%d",(int)currentdepth);
                        }
                        else if (fabs(fmod(cdd,10)/100-0)<0.01){
                                sprintf(tempdep,"%.1f",currentdepth);
                        }
                        else {
                                printf("%f %f\n",currentdepth,fmod(cdd,10)/100);
                                fprintf(stderr,"Only 1 digit after the decimal point is preserved for green's functions. Please change depth settings.\n");
                                exit(1);
                        }
                        sprintf(bbb1,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempddpz);
                        sprintf(bbb2,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempddsz);
                        sprintf(bbb3,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempdspz);
                        sprintf(bbb4,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempdssz);
                        sprintf(bbb5,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempsspz);
                        sprintf(bbb6,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempsssz);
                        strcpy(greenpath[ii].ddpzname,bbb1);
                        strcpy(greenpath[ii].ddszname,bbb2);
                        strcpy(greenpath[ii].dspzname,bbb3);
                        strcpy(greenpath[ii].dsszname,bbb4);
                        strcpy(greenpath[ii].sspzname,bbb5);
                        strcpy(greenpath[ii].ssszname,bbb6);
                }
        }
}

void read_greenpath_SH(greenpathformatSH greenpath[maxndepth*maxstanum], staformat stainfo[maxstanum], char greenchar1[], char greenchar2[], char greenchar3[], float bg_dep, float nd_dep, float dep_interval, int *ndepth, int num_sta){
        int i,j,k,ii;
	float currentdepth,cdd;
	char aaa[128],bbb1[512],bbb2[512],tempdep[128];
	char tempdsst[]="dsst",tempssst[]="ssst";
	float ndd=(nd_dep-bg_dep)/dep_interval*100;
        if (fabs(fmod(ndd,100))>0.01) {
		fprintf(stderr,"Depth range and interval are not correct.\n");
                exit(1);
	}
	*ndepth=(int)((nd_dep-bg_dep)/dep_interval+1);
	for (j=0;j<num_sta;j++) {
		strcpy(aaa,stainfo[j].greendist);
		for (i=0;i<*ndepth;i++) {
			ii=i*num_sta+j;
                	currentdepth=bg_dep+i*dep_interval;
			cdd=currentdepth*100+0.01;
                	if (fmod(currentdepth,1)<0.01) {
                        	sprintf(tempdep,"%d",(int)currentdepth);
                	}
                	else if (fabs(fmod(cdd,10)/100-0)<0.01){
				sprintf(tempdep,"%.1f",currentdepth);
                	}
                	else {
				printf("%f %f\n",currentdepth,fmod(cdd,10)/100);
                        	fprintf(stderr,"Only 1 digit after the decimal point is preserved for green's functions. Please change depth settings.\n");
                        	exit(1);
                	}
			sprintf(bbb1,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempdsst);
			sprintf(bbb2,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempssst);
			strcpy(greenpath[ii].dsstname,bbb1);
			strcpy(greenpath[ii].ssstname,bbb2);
		}
        }
}

void read_1dgreenpath(greenpathformat1d greenpath[maxndepth*maxdistnum], char greenchar1[], char greenchar2[], char greenchar3[],float bg_distloc, float nd_distloc, float distloc_interval, int *ndistloc, float bg_dep, float nd_dep, float dep_interval, int *ndepth){
        int i,j,k,ii;
        float currentdepth, cdd, currentdist, cld;
	char aaa[128],bbb1[512],bbb2[512],bbb3[512],bbb4[512],bbb5[512],bbb6[512],tempdep[128];
        char tempddpr[]="ddpr",tempddsr[]="ddsr",tempdspr[]="dspr",tempdssr[]="dssr",tempsspr[]="sspr",tempsssr[]="sssr",tempddpz[]="ddpz",tempddsz[]="ddsz",tempdspz[]="dspz",tempdssz[]="dssz",tempsspz[]="sspz",tempsssz[]="sssz",tempdsst[]="dsst",tempssst[]="ssst";
        float ndd=(nd_dep-bg_dep)/dep_interval*100;
        if (fabs(fmod(ndd,100))>0.01) {
                fprintf(stderr,"Depth range and interval are not correct.\n");
                exit(1);
        }
	if (fabs(fmod((nd_distloc-bg_distloc)/distloc_interval*100,100))>0.01) {
                fprintf(stderr,"Local distance range and interval are not correct.\n");
                exit(1);
        }
        *ndepth=(int)((nd_dep-bg_dep)/dep_interval+1);
	*ndistloc=(int)((nd_distloc-bg_distloc)/distloc_interval+1);
	for (j=0;j<*ndistloc;j++) {
                currentdist=bg_distloc+j*distloc_interval;
                cld=currentdist*100+0.01;
                if (fmod(currentdist,1)<0.01) {
                        sprintf(aaa,"%d",(int)currentdist);
                }
                else if (fabs(fmod(cld,10)/100-0)<0.01){
                        sprintf(aaa,"%.1f",currentdist);
                }
                else {
                        printf("%f %f\n",currentdist,fmod(cld,10)/100);
                        fprintf(stderr,"Only 1 digit after the decimal point is preserved for green's functions. Please change depth settings.\n");
                        exit(1);
                }
                for (i=0;i<*ndepth;i++) {
                        ii=i*(*ndistloc)+j;
                        currentdepth=bg_dep+i*dep_interval;
                        cdd=currentdepth*100+0.01;
                        if (fmod(currentdepth,1)<0.01) {
                                sprintf(tempdep,"%d",(int)currentdepth);
                        }
                        else if (fabs(fmod(cdd,10)/100-0)<0.01){
                                sprintf(tempdep,"%.1f",currentdepth);
                        }
                        else {
                                printf("%f %f\n",currentdepth,fmod(cdd,10)/100);
                                fprintf(stderr,"Only 1 digit after the decimal point is preserved for green's functions. Please change depth settings.\n");
                                exit(1);
                        }
			sprintf(bbb1,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempddpr);
                        sprintf(bbb2,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempddsr);
                        sprintf(bbb3,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempdspr);
                        sprintf(bbb4,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempdssr);
                        sprintf(bbb5,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempsspr);
                        sprintf(bbb6,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempsssr);
                        strcpy(greenpath[ii].ddprname,bbb1);
                        strcpy(greenpath[ii].ddsrname,bbb2);
                        strcpy(greenpath[ii].dsprname,bbb3);
                        strcpy(greenpath[ii].dssrname,bbb4);
                        strcpy(greenpath[ii].ssprname,bbb5);
                        strcpy(greenpath[ii].sssrname,bbb6);
                        sprintf(bbb1,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempddpz);
                        sprintf(bbb2,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempddsz);
                        sprintf(bbb3,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempdspz);
                        sprintf(bbb4,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempdssz);
                        sprintf(bbb5,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempsspz);
                        sprintf(bbb6,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempsssz);
                        strcpy(greenpath[ii].ddpzname,bbb1);
                        strcpy(greenpath[ii].ddszname,bbb2);
                        strcpy(greenpath[ii].dspzname,bbb3);
                        strcpy(greenpath[ii].dsszname,bbb4);
                        strcpy(greenpath[ii].sspzname,bbb5);
                        strcpy(greenpath[ii].ssszname,bbb6);
			sprintf(bbb1,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempdsst);
                        sprintf(bbb2,"%s%s%s%s%s%s",greenchar1,tempdep,greenchar2,aaa,greenchar3,tempssst);
                        strcpy(greenpath[ii].dsstname,bbb1);
                        strcpy(greenpath[ii].ssstname,bbb2);
                }
        }
}

void read_CMT(char CMTfile[], float *CMT){
	FILE *fp5;
        if ((fp5 = fopen(CMTfile, "r")) == NULL) {
                fprintf(stderr,"Cannot open input model file.\n");
                exit(1);
        }
	int i;
	for (i=0;i<5;i++) {
		fscanf(fp5,"%f",&CMT[i]);
	}
	fclose(fp5);
}

void rsac(char filename[], float bg_time, float nd_time, float deltat,  float *btime, float *t1, float xdata[maxnpts], int npts) { 
	float *x;
        int nerr,i;
        brsac(maxnpts*100,filename, &x, &nerr);
	getfhv("B",btime,&nerr);
	getfhv("T1",t1,&nerr);
	int nb=(*t1+bg_time-*btime)/deltat;
	for (i=0;i<npts;i++) {
		xdata[i]=x[nb+i];
	}
} 

void read_seismic_data(staformat stainfo[maxstanum], seisdataformat seisdata[maxstanum], float bg_time, float nd_time, float deltat, int gnpts, float gfl, float gfh, int num_sta, float *flt_sn, float *flt_sd, int *nsects, int velmul) { //Input: staformat stainfo, bg_time, nd_time, deltat, npts, deltat, fl, fh; Output: seisdataformat seisdata
	int i,j,k;
	long int order=4;
        char type[4]="BP";
	char proto[4]="BU";
	float pnl_sn[maxflen], pnl_sd[maxflen];
	for (j=0;j<num_sta;j++) {
		rsac(stainfo[j].stname, bg_time, nd_time, deltat, &seisdata[j].btime, &seisdata[j].t1, seisdata[j].xdata, gnpts);
		strcpy(seisdata[j].stname,stainfo[j].stname);
	}
	long int nsects1;
	design(order, type, proto, 1., 1., (double) gfl, (double) gfh, (double) deltat, pnl_sn, pnl_sd, &nsects1);
	for (i=0;i<maxflen;i++) {flt_sn[i]=pnl_sn[i]; flt_sd[i]=pnl_sd[i];}
	*nsects=(int) nsects1;
	int nsta2=num_sta/2;
	for (j=0;j<num_sta;j++) {
	apply(seisdata[j].xdata,(long int) gnpts,1,pnl_sn,pnl_sd,nsects1);	
		if (velmul>1) {
			if (j>=nsta2) {
				for (i=0;i<gnpts;i++) {
		                        seisdata[j].xdata[i]=seisdata[j].xdata[i]*velmul;
                		}
			}
		}
	}
}

void read_seismic_dataloc(stalocformat stainfo[maxlocnum], seisdatalocformat seisdata[maxlocnum], float bg_time, float nd_time, float deltat, int gnpts, float gfl, float gfh, int num_sta, float *flt_sn, float *flt_sd, int *nsects) {
        int i,j,k;
        long int order=4;
        char type[4]="BP";
        char proto[4]="BU";
        float pnl_sn[maxflen], pnl_sd[maxflen],pnl_sn2[maxflen], pnl_sd2[maxflen];
        for (j=0;j<num_sta;j++) {
                rsac(stainfo[j].stename, bg_time, nd_time, deltat, &seisdata[j].btime, &seisdata[j].t1, seisdata[j].edata, gnpts);
		rsac(stainfo[j].stnname, bg_time, nd_time, deltat, &seisdata[j].btime, &seisdata[j].t1, seisdata[j].ndata, gnpts);
		rsac(stainfo[j].stzname, bg_time, nd_time, deltat, &seisdata[j].btime, &seisdata[j].t1, seisdata[j].zdata, gnpts);
                strcpy(seisdata[j].stname,stainfo[j].stzname);
        }
        long int nsects1;
        design(order, type, proto, 1., 1., (double) gfl, (double) gfh, (double) deltat, pnl_sn, pnl_sd, &nsects1);
        for (i=0;i<maxflen;i++) {flt_sn[i]=pnl_sn[i]; flt_sd[i]=pnl_sd[i];}
        *nsects=(int) nsects1;
        for (j=0;j<num_sta;j++) {
	        apply(seisdata[j].edata,(long int) gnpts,1,pnl_sn,pnl_sd,nsects1);
		apply(seisdata[j].ndata,(long int) gnpts,1,pnl_sn,pnl_sd,nsects1);
		apply(seisdata[j].zdata,(long int) gnpts,1,pnl_sn,pnl_sd,nsects1);
		for (i=0;i<gnpts;i++) {
			seisdata[j].edata[i]=seisdata[j].edata[i]*stainfo[j].wte;
			seisdata[j].ndata[i]=seisdata[j].ndata[i]*stainfo[j].wtn;
			seisdata[j].zdata[i]=seisdata[j].zdata[i]*stainfo[j].wtz;
		}
        }
}

void read_green_data_P(greenpathformatP greenpath[maxndepth*maxstanum], float greendata[6*maxndepth*maxstanum][maxnpts], float greenbtime, float greenetime, float delta, int num_sta, int ndepth){
        int i,j,k,gfnpts,ii,l;
        float temp1,temp2;
        gfnpts=(int)((greenetime-greenbtime)/delta+1);
	int nsta2=num_sta/2;
        for (i=0;i<ndepth;i++) {
                for (l=0;l<num_sta;l++) {
                        ii=i*num_sta+l;
                        j=6*(i*num_sta+l);
                        rsac(greenpath[ii].ddpzname, greenbtime, greenetime, delta, &temp1, &temp2, greendata[j], gfnpts);
                        rsac(greenpath[ii].ddszname, greenbtime, greenetime, delta, &temp1, &temp2, greendata[j+1], gfnpts);
                        rsac(greenpath[ii].dspzname, greenbtime, greenetime, delta, &temp1, &temp2, greendata[j+2], gfnpts);
                        rsac(greenpath[ii].dsszname, greenbtime, greenetime, delta, &temp1, &temp2, greendata[j+3], gfnpts);
                        rsac(greenpath[ii].sspzname, greenbtime, greenetime, delta, &temp1, &temp2, greendata[j+4], gfnpts);
                        rsac(greenpath[ii].ssszname, greenbtime, greenetime, delta, &temp1, &temp2, greendata[j+5], gfnpts);
			if (l<nsta2) {
				for (k=0;k<gfnpts;k++) {
                                greendata[j][k]=greendata[j][k]*1e5;
                                greendata[j+1][k]=greendata[j+1][k]*1e5;
                                greendata[j+2][k]=greendata[j+2][k]*1e5;
                                greendata[j+3][k]=greendata[j+3][k]*1e5;
                                greendata[j+4][k]=greendata[j+4][k]*1e5;
                                greendata[j+5][k]=greendata[j+5][k]*1e5;
                        	}
			} else {
				for (k=0;k<gfnpts;k++) {
                                greendata[j][k]=greendata[j][k]*1e5*2;
                                greendata[j+1][k]=greendata[j+1][k]*1e5*2;
                                greendata[j+2][k]=greendata[j+2][k]*1e5*2;
                                greendata[j+3][k]=greendata[j+3][k]*1e5*2;
                                greendata[j+4][k]=greendata[j+4][k]*1e5*2;
                                greendata[j+5][k]=greendata[j+5][k]*1e5*2;
                        	}
			}
                }
        }
}

void read_green_data_SH(greenpathformatSH greenpath[maxndepth*maxstanum], float greendata[2*maxndepth*maxstanum][maxnpts], float greenbtime, float greenetime, float delta, int num_sta, int ndepth){
	int i,j,k,gfnpts,ii,l;
	float temp1,temp2;
	gfnpts=(int)((greenetime-greenbtime)/delta+1);
	for (i=0;i<ndepth;i++) {
		for (l=0;l<num_sta;l++) {
			ii=i*num_sta+l;
                	j=2*(i*num_sta+l);
			rsac(greenpath[ii].dsstname, greenbtime, greenetime, delta, &temp1, &temp2, greendata[j], gfnpts);
			rsac(greenpath[ii].ssstname, greenbtime, greenetime, delta, &temp1, &temp2, greendata[j+1], gfnpts);
                	for (k=0;k<gfnpts;k++) {
                        	greendata[j][k]=greendata[j][k]*1e5;
				greendata[j+1][k]=greendata[j+1][k]*1e5;
                	}
                }
	}
}


void read_green_data_Rayl(greenpathformat1d greenpath[maxndepth*maxdistnum], stalocformat stainfo[maxlocnum], float greendata[14*maxndepth*maxdistnum][maxnpts], float greenbtime, float greenetime, float delta, int ndistloc, int ndepth){
        int i,j,k,gfnpts,gfnptsbg,ii,l,m,n;
        float temp1,temp2;
        gfnpts=(int)((greenetime-greenbtime)/delta+1);
	gfnptsbg=(int)((0-greenbtime)/delta);
	float greendata1[14][maxnpts];
        for (i=0;i<ndepth;i++) {
		//printf("%s\n",greenpath[i*ndistloc+ndistloc-1].ssstname);
                for (l=0;l<ndistloc;l++) {
                        ii=i*ndistloc+l;
                        j=14*ii;	
                        rsac(greenpath[ii].ddpzname, 0, greenetime, delta, &temp1, &temp2, greendata1[0], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].ddszname, 0, greenetime, delta, &temp1, &temp2, greendata1[1], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].dspzname, 0, greenetime, delta, &temp1, &temp2, greendata1[2], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].dsszname, 0, greenetime, delta, &temp1, &temp2, greendata1[3], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].sspzname, 0, greenetime, delta, &temp1, &temp2, greendata1[4], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].ssszname, 0, greenetime, delta, &temp1, &temp2, greendata1[5], gfnpts-gfnptsbg);
			rsac(greenpath[ii].ddprname, 0, greenetime, delta, &temp1, &temp2, greendata1[6], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].ddsrname, 0, greenetime, delta, &temp1, &temp2, greendata1[7], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].dsprname, 0, greenetime, delta, &temp1, &temp2, greendata1[8], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].dssrname, 0, greenetime, delta, &temp1, &temp2, greendata1[9], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].ssprname, 0, greenetime, delta, &temp1, &temp2, greendata1[10], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].sssrname, 0, greenetime, delta, &temp1, &temp2, greendata1[11], gfnpts-gfnptsbg);
			rsac(greenpath[ii].dsstname, 0, greenetime, delta, &temp1, &temp2, greendata1[12], gfnpts-gfnptsbg);
                        rsac(greenpath[ii].ssstname, 0, greenetime, delta, &temp1, &temp2, greendata1[13], gfnpts-gfnptsbg);
			for (k=0;k<gfnptsbg;k++) {
				for (n=0;n<14;n++) {
					greendata[j+i][k]=0;
				}
			}
                        for (k=gfnptsbg;k<gfnpts;k++) {
                                greendata[j][k]=greendata1[0][k-gfnptsbg]*1e5;
                                greendata[j+1][k]=greendata1[1][k-gfnptsbg]*1e5;
                                greendata[j+2][k]=greendata1[2][k-gfnptsbg]*1e5;
                                greendata[j+3][k]=greendata1[3][k-gfnptsbg]*1e5;
                                greendata[j+4][k]=greendata1[4][k-gfnptsbg]*1e5;
                                greendata[j+5][k]=greendata1[5][k-gfnptsbg]*1e5;
				greendata[j+6][k]=greendata1[6][k-gfnptsbg]*1e5;
                                greendata[j+7][k]=greendata1[7][k-gfnptsbg]*1e5;
                                greendata[j+8][k]=greendata1[8][k-gfnptsbg]*1e5;
                                greendata[j+9][k]=greendata1[9][k-gfnptsbg]*1e5;
                                greendata[j+10][k]=greendata1[10][k-gfnptsbg]*1e5;
                                greendata[j+11][k]=greendata1[11][k-gfnptsbg]*1e5;
				greendata[j+12][k]=greendata1[12][k-gfnptsbg]*1e5;
                                greendata[j+13][k]=greendata1[13][k-gfnptsbg]*1e5;
                        }
			//printf("%f %f %f\n",stainfo[l].wte,stainfo[l].wtn,stainfo[l].wtz);
                }
        }
}


void sub_init_(float *btimep, float *etimep, float *btimesh, float *etimesh, float *btimerayl, float *etimerayl, float *stfbtime, float *stfetime, int *fnptsp, int *fnptssh, int *fnptsrayl, float *delta, int *nstap, int *nstash, int *nstarayl, int *nsubev, float *initmodel, char stainfop_stname[maxstanum][64], float *stainfop_stdist, float *stainfop_staz, float *stainfop_stcp, float *stainfop_stcs, char stainfosh_stname[maxstanum][64], float *stainfosh_stdist, float *stainfosh_staz, float *stainfosh_stcp, float *stainfosh_stcs, float *stainforayl_stlo, float *stainforayl_stla, float stainforayl_wt[maxlocnum][3], char seisdatap_stname[maxstanum][64], float *seisdatap_btime, float *seisdatap_t1, float seisdatap_xdata[maxstanum][maxnpts], char seisdatash_stname[maxstanum][64], float *seisdatash_btime, float *seisdatash_t1, float seisdatash_xdata[maxstanum][maxnpts], char seisdatarayl_stname[maxlocnum][64], float *seisdatarayl_btime, float *seisdatarayl_t1, float seisdatarayl_xdata[maxlocnum*3][maxnpts], float greendataP[6*maxndepth*maxstanum][maxnpts], float greendataSH[2*maxndepth*maxstanum][maxnpts], float greendatarayl[14*maxndepth*maxdistnum][maxnpts], float *fltsnp, float *fltsdp, float *fltsnsh, float *fltsdsh, float *fltsnrayl, float *fltsdrayl, int *nsectp, int *nsectsh, int *nsectrayl, float *wtp, float *wtsh, float *wtrayl, float *alpha, float *cmtscaling, float *dcconstrain, float *evlo, float *evla, float *evvp, float *evvs, float *begin_localdist, float *end_localdist, float *interval_localdist, int *n_distloc, float *begin_dep, float *end_dep, float *interval_dep, int *n_dep, int *nunknown) {
	float **sxdatap, **sxdatash, **sxdatarayl;
	float sfltsnp[maxflen], sfltsdp[maxflen], sfltsnsh[maxflen], sfltsdsh[maxflen], sfltsnrayl[maxflen], sfltsdrayl[maxflen];
	int i,j,k,snsectp,snsectsh, snsectrayl;
	staformat stainfoP[maxstanum],stainfoSH[maxstanum]; 
	stalocformat stainforayl[maxlocnum];
	modelformat modelvec[maxsubeve];
	greenpathformatP* greenpathP;
	greenpathformatSH* greenpathSH;
	greenpathformat1d* greenpathrayl;
	seisdataformat* seisdataP;
	seisdataformat* seisdataSH;
	seisdatalocformat* seisdatarayl;
	float bg_time_P=-99999, nd_time_P=-99999, bg_time_SH=-99999, nd_time_SH=-99999,  bg_time_rayl=-99999, nd_time_rayl=-99999, stf_btime=-99999, stf_etime=-99999, deltat=-99999, flp=-99999, fhp=-99999, fls=-99999, fhs=-99999, flrayl=-99999, fhrayl=-99999, weightp=-99999, weightsh=-99999, weightrayl=-99999, alphatemp=-99999, cmttemp=-99999, dctemp=-99999, bg_dep=-99999, nd_dep=-99999, dep_interval=-99999, fevlo=-99999, fevla=-99999, fevvp=-99999, fevvs=-99999, bg_distloc=-99999, nd_distloc=-99999, distloc_interval=-99999;
	int nptsP=-99999, nptsSH=-99999, nptsrayl=-99999, num_sta_P=-99999, num_sta_SH=-99999, num_sta_rayl=-99999, num_subevent=-99999, ndistloc=-99999, ndepth=-99999;
	char stainfofileP[64], stainfofileSH[64], stainfofilerayl[64], inputmodelfile[64], greenchar1[256], greenchar2[64], greenchar3[64];
	float mdllist[120];
	read_parameter(&bg_time_P, &nd_time_P, &bg_time_SH, &nd_time_SH, &bg_time_rayl, &nd_time_rayl, &stf_btime, &stf_etime, &deltat, &nptsP, &nptsSH,&nptsrayl, &num_sta_P, &num_sta_SH, &num_sta_rayl, stainfofileP, stainfofileSH, stainfofilerayl, &weightp, &weightsh, &weightrayl, inputmodelfile, greenchar1, greenchar2, greenchar3, &fevlo, &fevla, &fevvp, &fevvs, &bg_distloc, &nd_distloc, &distloc_interval, &bg_dep, &nd_dep, &dep_interval, &num_subevent, &flp, &fhp, &fls, &fhs, &flrayl, &fhrayl, &alphatemp, &cmttemp, &dctemp);
	*wtp=weightp;
	*wtsh=weightsh;
	*wtrayl=weightrayl;
	if (weightp<=0 || weightsh<=0 || weightrayl==0) {
		fprintf(stderr,"Please check your weighting settings! \n");
        	exit(1);
	} else  {
		sxdatap=alloc2d(maxstanum,maxnpts);
                greenpathP=(greenpathformatP*)malloc(sizeof(greenpathformatP)*maxndepth*maxstanum);
                seisdataP=(seisdataformat*)malloc(sizeof(seisdataformat)*maxstanum);
		sxdatash=alloc2d(maxstanum,maxnpts);
                greenpathSH=(greenpathformatSH*)malloc(sizeof(greenpathformatSH)*maxndepth*maxstanum);
                seisdataSH=(seisdataformat*)malloc(sizeof(seisdataformat)*maxstanum);
		sxdatarayl=alloc2d(maxlocnum*3,maxnpts);
		greenpathrayl=(greenpathformat1d*)malloc(sizeof(greenpathformat1d)*maxndepth*maxdistnum);
                seisdatarayl=(seisdatalocformat*)malloc(sizeof(seisdatalocformat)*maxlocnum);
		read_stainfo(stainfoP, stainfofileP, num_sta_P);
		read_stainfo(stainfoSH, stainfofileSH, num_sta_SH);
		read_stalocinfo(stainforayl, stainfofilerayl, num_sta_rayl);
		read_initialmodel(modelvec, mdllist, inputmodelfile, num_subevent);
		read_greenpath_P(greenpathP, stainfoP, greenchar1, greenchar2, greenchar3, bg_dep, nd_dep, dep_interval, &ndepth, num_sta_P);
		read_greenpath_SH(greenpathSH, stainfoSH, greenchar1, greenchar2, greenchar3, bg_dep, nd_dep, dep_interval, &ndepth, num_sta_SH);
		read_1dgreenpath(greenpathrayl, greenchar1, greenchar2, greenchar3, bg_distloc, nd_distloc, distloc_interval, &ndistloc, bg_dep, nd_dep, dep_interval, &ndepth);
		read_seismic_data(stainfoP, seisdataP, bg_time_P, nd_time_P, deltat, nptsP, flp, fhp, num_sta_P, sfltsnp, sfltsdp, &snsectp, 2);
		read_seismic_data(stainfoSH, seisdataSH, bg_time_SH, nd_time_SH, deltat, nptsSH, fls, fhs, num_sta_SH, sfltsnsh, sfltsdsh, &snsectsh, 1);
		read_seismic_dataloc(stainforayl, seisdatarayl, bg_time_rayl, nd_time_rayl, deltat, nptsrayl, flrayl, fhrayl, num_sta_rayl, sfltsnrayl, sfltsdrayl, &snsectrayl);
		read_green_data_P(greenpathP, greendataP, bg_time_P, nd_time_P, deltat, num_sta_P, ndepth);
		read_green_data_SH(greenpathSH, greendataSH, bg_time_SH, nd_time_SH, deltat, num_sta_SH, ndepth);
		read_green_data_Rayl(greenpathrayl, stainforayl, greendatarayl, bg_time_rayl, nd_time_rayl, deltat, ndistloc, ndepth);
		*fnptsp=nptsP;
                *btimep=bg_time_P;
                *etimep=nd_time_P;
		*fnptsrayl=nptsrayl;
                *btimerayl=bg_time_rayl;
                *etimerayl=nd_time_rayl;
                *stfbtime=stf_btime;
                *stfetime=stf_etime;
                *delta=deltat;
                *nstap=num_sta_P;
		*nstarayl=num_sta_rayl;
                *nsubev=num_subevent;
                *nsectp=snsectp;
		*nsectrayl=snsectrayl;
                *alpha=alphatemp;
		*evlo=fevlo;
		*evla=fevla;
		*evvp=fevvp;
                *evvs=fevvs;
                *cmtscaling=cmttemp;
                *dcconstrain=dctemp;
                *nunknown=num_subevent*5;
                *begin_dep=bg_dep;
                *end_dep=nd_dep;
                *interval_dep=dep_interval;
		*n_dep=ndepth;
		*begin_localdist=bg_distloc;
		*end_localdist=nd_distloc;
		*interval_localdist=distloc_interval;
		*n_distloc=ndistloc;
		*fnptssh=nptsSH;
                *btimesh=bg_time_SH;
                *etimesh=nd_time_SH;
                *nstash=num_sta_SH;
                *nsectsh=snsectsh;
		for (i=0;i<maxflen;i++) { 
                        fltsnp[i]=sfltsnp[i];
                        fltsdp[i]=sfltsdp[i];
			fltsnrayl[i]=sfltsnrayl[i];
                        fltsdrayl[i]=sfltsdrayl[i];
			fltsnsh[i]=sfltsnsh[i];
                        fltsdsh[i]=sfltsdsh[i];
                }
                for (i=0;i<num_subevent*7;i++) {
                        initmodel[i]=mdllist[i];
                }
                char stastnamep[maxstanum][64], seisstnamep[maxstanum][64], stastnamesh[maxstanum][64], seisstnamesh[maxstanum][64];
		char stastnamerayl[maxlocnum][64], seisstnamerayl[maxlocnum][64];
                float sstdistp[maxstanum], sstazp[maxstanum], sstcpp[maxstanum], sstcsp[maxstanum], sbtimep[maxstanum], st1p[maxstanum], sstdistsh[maxstanum], sstazsh[maxstanum], sstcpsh[maxstanum], sstcssh[maxstanum], sbtimesh[maxstanum], st1sh[maxstanum], sraylwt[maxlocnum][3];
		float sstlorayl[maxlocnum], sstlarayl[maxlocnum], sbtimerayl[maxlocnum], st1rayl[maxlocnum];
                convert_staformat_arrays(stainfoP, stastnamep, sstdistp, sstazp, sstcpp, sstcsp);
                convert_seisformat_arrays(seisdataP, seisstnamep, sbtimep, st1p, sxdatap);
		convert_staformat_arrays(stainfoSH, stastnamesh, sstdistsh, sstazsh, sstcpsh, sstcssh);
                convert_seisformat_arrays(seisdataSH, seisstnamesh, sbtimesh, st1sh, sxdatash);
		convert_stalocformat_arrays(stainforayl, sstlorayl, sstlarayl, sraylwt);
                convert_seislocformat_arrays(seisdatarayl, seisstnamerayl, sbtimerayl, st1rayl, sxdatarayl);
                for (i=0;i<maxstanum;i++) {
                        for (j=0;j<64;j++) {
                                stainfop_stname[i][j]=stastnamep[i][j];
                                seisdatap_stname[i][j]=seisstnamep[i][j];
				stainfosh_stname[i][j]=stastnamesh[i][j];
                                seisdatash_stname[i][j]=seisstnamesh[i][j];
                        }
                        stainfop_stdist[i]=sstdistp[i];
                        stainfop_staz[i]=sstazp[i];
                        stainfop_stcp[i]=sstcpp[i];
                        stainfop_stcs[i]=sstcsp[i];
                        seisdatap_btime[i]=sbtimep[i];
                        seisdatap_t1[i]=st1p[i];
			stainfosh_stdist[i]=sstdistsh[i];
                        stainfosh_staz[i]=sstazsh[i];
                        stainfosh_stcp[i]=sstcpsh[i];
                        stainfosh_stcs[i]=sstcssh[i];
                        seisdatash_btime[i]=sbtimesh[i];
                        seisdatash_t1[i]=st1sh[i];
                        for (j=0;j<maxnpts;j++) {
                                seisdatap_xdata[i][j]=sxdatap[i][j];
				seisdatash_xdata[i][j]=sxdatash[i][j];
                        }
                }
		for (i=0;i<maxlocnum;i++) {
			stainforayl_wt[i][0]=sraylwt[i][0];
			stainforayl_wt[i][1]=sraylwt[i][1];
			stainforayl_wt[i][2]=sraylwt[i][2];
			for (j=0;j<64;j++) {
				seisdatarayl_stname[i][j]=seisstnamerayl[i][j];
			}
			stainforayl_stlo[i]=sstlorayl[i];
                        stainforayl_stla[i]=sstlarayl[i];
                        seisdatarayl_btime[i]=sbtimerayl[i];
                        seisdatarayl_t1[i]=st1rayl[i];
			for (j=0;j<maxnpts;j++) {
                                seisdatarayl_xdata[3*i][j]=sxdatarayl[3*i][j];
                                seisdatarayl_xdata[3*i+1][j]=sxdatarayl[3*i+1][j];
                                seisdatarayl_xdata[3*i+2][j]=sxdatarayl[3*i+2][j];
                        }
		}
                free2d(sxdatap);
                free(seisdataP);
                free(greenpathP);
		free2d(sxdatash);
                free(seisdataSH);
                free(greenpathSH);
		free2d(sxdatarayl);
                free(seisdatarayl);
                free(greenpathrayl);
                sxdatap=NULL;
                seisdataP=NULL;
                greenpathP=NULL;
                sxdatash=NULL;
                seisdataSH=NULL;
                greenpathSH=NULL;
		sxdatarayl=NULL;
                seisdatarayl=NULL;
                greenpathrayl=NULL;
	}
}

