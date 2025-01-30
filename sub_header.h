#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sacsubc.h"
#define rad2deg 57.29578
#define maxsubeve 7
#define maxstanum 200
#define maxdistnum 250
#define maxlocnum 70
#define maxnum3d 1
#define maxnpts 500
#define maxflen 30
#define maxndepth 41
#define min(a,b) ( ((a)>(b)) ? (b):(a) )
#define max(a,b) ( ((a)>(b)) ? (a):(b) )

typedef struct
{
        char stname[64];
        float stdist;
        float staz;
        float stcp;
        float stcs;
	char greendist[64];
} staformat;

typedef struct
{
        char stename[64];
	char stnname[64];
	char stzname[64];
        float stlo;
        float stla;
	char stind[64];
	float wte;
	float wtn;
	float wtz;
} stalocformat;

typedef struct
{
        float stlo;
        float stla;
        char stind[64];
} green3dformat;

typedef struct
{
        float ct; // Centroid time of sub-event
        float cx;
        float cy;
        float L; // Rupture length of sub-event
        float Vr; // Rupture velocity of sub-event (~directivity)
        float theta_rupture; // Direction to which rupture of current sub-event goes (~directivity)
	float dep;
} modelformat;

typedef struct
{
        char stname[64];
        float btime;
        float t1;
        float xdata[maxnpts];
} seisdataformat;

typedef struct
{
        char stname[64];
        float btime;
        float t1;
        float edata[maxnpts];
	float ndata[maxnpts];
	float zdata[maxnpts];
} seisdatalocformat;

typedef struct
{
        char ddpzname[512];
        char ddszname[512];
        char dspzname[512];
        char dsszname[512];
        char sspzname[512];
        char ssszname[512];
} greenpathformatP;

typedef struct
{
        char dsstname[512];
        char ssstname[512];
} greenpathformatSH;

typedef struct
{
	char ddprname[512];
        char ddsrname[512];
        char dsprname[512];
        char dssrname[512];
        char ssprname[512];
        char sssrname[512];
	char ddpzname[512];
        char ddszname[512];
        char dspzname[512];
        char dsszname[512];
        char sspzname[512];
        char ssszname[512];
	char dsstname[512];
        char ssstname[512];
} greenpathformat1d;
void read_parameter(float *bg_time_P, float *nd_time_P, float *bg_time_SH, float *nd_time_SH, float *bg_time_Rayl, float *nd_time_Rayl, float *stf_btime, float *stf_etime, float *deltat, int *nptsP, int *nptsSH, int *nptsRayl, int *num_sta_P, int *num_sta_SH,  int *num_sta_Rayl, char stainfofileP[64], char stainfofileSH[64], char stainfofileRayl[64], float *weightP, float *weightSH, float *weightRayl, char inputmodelfile[64], char greenchar1[256], char greenchar2[64], char greenchar3[64], float *evlo, float *evla,float *evvp, float *evvs, float *bg_distloc, float *nd_distloc, float *distloc_interval, float *bg_dep, float *nd_dep, float *dep_interval, int *num_subevent, float *flp, float *fhp, float *fls, float *fhs, float *flr, float *fhr, float *alpha, float *cmtscaling, float *dcconstrain);
void read_stainfo(staformat stainfo[maxstanum], char stainfofile[], int num_sta);
void read_stalocinfo(stalocformat stainfo[maxlocnum], char stainfofile[], int num_sta);
void read_greeninfo3d(green3dformat greeninfo3d[maxnum3d], char greenfile3d[]) ;
void read_initialmodel(modelformat modelvec[maxstanum], float *modellist, char inputmodelfile[], int num_subevent);
void read_greenpath_P(greenpathformatP greenpath[maxndepth*maxstanum], staformat stainfo[maxstanum], char greenchar1[], char greenchar2[], char greenchar3[], float bg_dep, float nd_dep, float dep_interval, int *ndepth, int num_sta);
void read_greenpath_SH(greenpathformatSH greenpath[maxndepth*maxstanum], staformat stainfo[maxstanum], char greenchar1[], char greenchar2[], char greenchar3[], float bg_dep, float nd_dep, float dep_interval, int *ndepth, int num_sta);
void read_1dgreenpath(greenpathformat1d greenpath[maxndepth*maxdistnum], char greenchar1[], char greenchar2[], char greenchar3[],float bg_distloc, float nd_distloc, float distloc_interval, int *ndistloc, float bg_dep, float nd_dep, float dep_interval, int *ndepth);
void read_CMT(char CMTfile[], float *CMT);
void rsac(char filename[], float bg_time, float nd_time, float deltat,  float *btime, float *t1, float xdata[maxnpts], int npts);
void read_seismic_data(staformat stainfo[maxstanum], seisdataformat seisdata[maxstanum], float bg_time, float nd_time, float deltat, int gnpts, float gfl, float gfh, int num_sta, float *flt_sn, float *flt_sd, int *nsects, int velmul) ;
void read_seismic_dataloc(stalocformat stainfo[maxlocnum], seisdatalocformat seisdata[maxlocnum], float bg_time, float nd_time, float deltat, int gnpts, float gfl, float gfh, int num_sta, float *flt_sn, float *flt_sd, int *nsects);
void read_green_data_P(greenpathformatP greenpath[maxndepth*maxstanum], float greendata[6*maxndepth*maxstanum][maxnpts], float greenbtime, float greenetime, float delta, int num_sta, int ndepth);
void read_green_data_SH(greenpathformatSH greenpath[maxndepth*maxstanum], float greendata[2*maxndepth*maxstanum][maxnpts], float greenbtime, float greenetime, float delta, int num_sta, int ndepth);
void read_green_data_Rayl(greenpathformat1d greenpath[maxndepth*maxdistnum], stalocformat stainfo[maxlocnum], float greendata[14*maxndepth*maxdistnum][maxnpts], float greenbtime, float greenetime, float delta, int ndistloc, int ndepth);
void sub_init_(float *btimep, float *etimep, float *btimesh, float *etimesh, float *btimerayl, float *etimerayl, float *stfbtime, float *stfetime, int *fnptsp, int *fnptssh, int *fnptsrayl, float *delta, int *nstap, int *nstash, int *nstarayl, int *nsubev, float *initmodel, char stainfop_stname[maxstanum][64], float *stainfop_stdist, float *stainfop_staz, float *stainfop_stcp, float *stainfop_stcs, char stainfosh_stname[maxstanum][64], float *stainfosh_stdist, float *stainfosh_staz, float *stainfosh_stcp, float *stainfosh_stcs, float *stainforayl_stlo, float *stainforayl_stla, float stainforayl_wt[maxlocnum][3], char seisdatap_stname[maxstanum][64], float *seisdatap_btime, float *seisdatap_t1, float seisdatap_xdata[maxstanum][maxnpts], char seisdatash_stname[maxstanum][64], float *seisdatash_btime, float *seisdatash_t1, float seisdatash_xdata[maxstanum][maxnpts], char seisdatarayl_stname[maxlocnum][64], float *seisdatarayl_btime, float *seisdatarayl_t1, float seisdatarayl_xdata[maxlocnum*3][maxnpts], float greendataP[6*maxndepth*maxstanum][maxnpts], float greendataSH[2*maxndepth*maxstanum][maxnpts], float greendatarayl[14*maxndepth*maxdistnum][maxnpts], float *fltsnp, float *fltsdp, float *fltsnsh, float *fltsdsh, float *fltsnrayl, float *fltsdrayl, int *nsectp, int *nsectsh, int *nsectrayl, float *wtp, float *wtsh, float *wtrayl, float *alpha, float *cmtscaling, float *dcconstrain, float *evlo, float *evla, float *evvp, float *evvs, float *begin_localdist, float *end_localdist, float *interval_localdist, int *n_distloc, float *begin_dep, float *end_dep, float *interval_dep, int *n_dep, int *nunknown);
void sub_forward_(float *btimep, float *etimep, float *btimesh, float *etimesh, float *btimerayl, float *etimerayl, float *stfbtime, float *stfetime, int *ptsp, int *ptssh, int *ptsrayl, float *delta, int *nstap, int *nstash, int *nstarayl, int *nsubev, float *initmodel, char stainfop_stname[maxstanum][64], float *stainfop_stdist, float *stainfop_staz, float *stainfop_stcp, float *stainfop_stcs, char stainfosh_stname[maxstanum][64], float *stainfosh_stdist, float *stainfosh_staz, float *stainfosh_stcp, float *stainfosh_stcs, float *stainforayl_stlo, float *stainforayl_stla, float stainforayl_wt[maxlocnum][3], char seisdatap_stname[maxstanum][64], float *seisdatap_btime, float *seisdatap_t1, float seisdatap_xdata[maxstanum][maxnpts], char seisdatash_stname[maxstanum][64], float *seisdatash_btime, float *seisdatash_t1, float seisdatash_xdata[maxstanum][maxnpts], char seisdatarayl_stname[maxlocnum][64], float *seisdatarayl_btime, float *seisdatarayl_t1, float seisdatarayl_xdata[maxlocnum*3][maxnpts], float greendatap[maxndepth*maxstanum*6][maxnpts], float greendatash[maxndepth*maxstanum*2][maxnpts], float greendatarayl[14*maxndepth*maxdistnum][maxnpts], float *residualall, float *fltsnp, float *fltsdp, float *fltsnsh, float *fltsdsh, float *fltsnrayl, float *fltsdrayl, int *nsectp, int *nsectsh, int *nsectrayl, float *wtp, float *wtsh, float *wtrayl, float *CMTsolution, float *alphat, float *cmtscaling, float *dc, float *evlo, float *evla, float *evvp, float *evvs, float *begin_localdist, float *end_localdist, float *interval_localdist, int *n_distloc, float *begin_dep, float *end_dep, float *interval_dep, int *n_dep);
void read_modelvec(float *modellist, int num_subevent, modelformat modelvec[maxstanum]);
void convolve(float* kernel, float* kernels, float* in1, float* in2, float* out1, float* out2, int kernel_length, int length);
void simpleconvolve(float* kernel, float* in1, float* out1, int kernel_length, int length);
void genGaussP(int gnpts, float gdeltat, float gsig, float gmu, float gsigs, float gmus, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float *grnddpz, float *grnddsz, float *grndspz, float *grndssz, float *grnsspz, float *grnsssz, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, int nn, float dura);
void genGaussP_interp(int gnpts, float gdeltat, float gsig, float gmu, float gsigs, float gmus, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float dz, float *grnddpz0, float *grnddpz1, float *grnddsz0, float *grnddsz1, float *grndspz0, float *grndspz1, float *grndssz0, float *grndssz1, float *grnsspz0, float *grnsspz1, float *grnsssz0, float *grnsssz1, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, int nn, float dura);
void genGaussSH(int gnpts, float gdeltat, float gsig, float gmu, float gsigs, float gmus, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float *grndsst, float *grnssst, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, int nn, float dura);
void genGaussSH_interp(int gnpts, float gdeltat, float gsig, float gmu, float gsigs, float gmus, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5,  float dz, float *grndsst0, float *grndsst1, float *grnssst0, float *grnssst1, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, int nn, float dura) ;
void genGaussRayle_interp(int gnpts, float gdeltat, float gsig, float gmu, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float dx, float dz, float *grnddpr000, float *grnddpr010, float *grnddpr001, float *grnddpr011, float *grnddsr000, float *grnddsr010, float *grnddsr001, float *grnddsr011, float *grndspr000, float *grndspr010, float *grndspr001, float *grndspr011, float *grndssr000, float *grndssr010, float *grndssr001, float *grndssr011, float *grnsspr000, float *grnsspr010, float *grnsspr001, float *grnsspr011, float *grnsssr000, float *grnsssr010, float *grnsssr001, float *grnsssr011, float *grndsst000, float *grndsst010, float *grndsst001, float *grndsst011, float *grnssst000, float *grnssst010, float *grnssst001, float *grnssst011, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, float wt, float dura, int nn);
void genGaussRayln_interp(int gnpts, float gdeltat, float gsig, float gmu, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float dx, float dz, float *grnddpr000, float *grnddpr010, float *grnddpr001, float *grnddpr011, float *grnddsr000, float *grnddsr010, float *grnddsr001, float *grnddsr011, float *grndspr000, float *grndspr010, float *grndspr001, float *grndspr011, float *grndssr000, float *grndssr010, float *grndssr001, float *grndssr011, float *grnsspr000, float *grnsspr010, float *grnsspr001, float *grnsspr011, float *grnsssr000, float *grnsssr010, float *grnsssr001, float *grnsssr011, float *grndsst000, float *grndsst010, float *grndsst001, float *grndsst011, float *grnssst000, float *grnssst010, float *grnssst001, float *grnssst011, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, float wt, float dura, int nn);
void genGaussRaylz_interp(int gnpts, float gdeltat, float gsig, float gmu, int stfnpts, float stfbg, float stfnd, float databg, float datand, float *tdata1, float az, float *g1, float *g2, float *g3, float *g4, float *g5, float dx, float dz, float *grnddpz000, float *grnddpz010, float *grnddpz001, float *grnddpz011, float *grnddsz000, float *grnddsz010, float *grnddsz001, float *grnddsz011, float *grndspz000, float *grndspz010, float *grndspz001, float *grndspz011, float *grndssz000, float *grndssz010, float *grndssz001, float *grndssz011, float *grnsspz000, float *grnsspz010, float *grnsspz001, float *grnsspz011, float *grnsssz000, float *grnsssz010, float *grnsssz001, float *grnsssz011, float flt_sn[maxflen], float flt_sd[maxflen], int nsects, float wt, float dura, int nn);
void convert_staformat_arrays(staformat sstainfo[maxstanum], char sstname[maxstanum][64], float sstdist[maxstanum], float sstaz[maxstanum], float sstcp[maxstanum], float sstcs[maxstanum]);
void convert_stalocformat_arrays(stalocformat stainfo[maxlocnum], float stlo[maxlocnum], float stla[maxlocnum], float wt[maxlocnum][3]) ;
void convert_green3dformat_arrays(green3dformat greeninfo3d[maxnum3d], float lo3d[maxnum3d], float la3d[maxnum3d], char ind3d[maxnum3d][64]);
void convert_seisformat_arrays(seisdataformat sseisdata[maxstanum], char sstname[maxstanum][64], float sbtime[maxstanum], float st1[maxstanum], float **sxdata) ;
void convert_seislocformat_arrays(seisdatalocformat sseisdata[maxlocnum], char sstname[maxlocnum][64], float sbtime[maxlocnum], float st1[maxlocnum], float **sxdata);
void convert_arrays_green3dformat(float lo3d[maxnum3d], float la3d[maxnum3d], char ind3d[maxnum3d][64], green3dformat greeninfo3d[maxnum3d]);
void convert_arrays_staformat(char sstname[maxstanum][64], float sstdist[maxstanum], float sstaz[maxstanum], float sstcp[maxstanum], float sstcs[maxstanum], staformat sstainfo[maxstanum]) ;
void convert_arrays_stalocformat(float sstlo[maxlocnum], float sstla[maxlocnum], float wt[maxlocnum][3], stalocformat sstainfo[maxlocnum]);
void convert_arrays_seisformat(char sstname[maxstanum][64], float sbtime[maxstanum], float st1[maxstanum], float sxdata[][maxnpts], seisdataformat sseisdata[maxstanum]);
void convert_arrays_seislocformat(char sstname[maxlocnum][64], float sbtime[maxlocnum], float st1[maxlocnum], float sxdata[][maxnpts], seisdatalocformat sseisdata[maxlocnum]);
void outputsac(int npts, float *arr, float dt, char *filename);
void genstf(float gdeltat, float gsig, float gmu, int stfnpts, float *tdata1, float dura, int nn);
void output_forward(int num_subevent, int num_sta, modelformat modelvec[maxstanum], int npts, int stfnpts, float deltat, float *tdata, double *tempgxm, double *rclvd, double *GG, double *tempdx, float *Cross_Corre_Coefficient, float *l2res, float residual, float normresidual, float finalresidual, float *T_Cross);
void output_forward_joint(int num_subevent, int num_stap, int num_stash, int num_starayl, modelformat modelvec[maxstanum], int nptsp, int nptssh, int nptsrayl, int stfnpts, float deltat, float *tdata, double *tempgxm, double *rclvd, double *GG, double *tempdx, float *pCCC, float *shCCC, float *raylCCC, float *pT_Cross, float *shT_Cross, float *raylT_Cross, float *l2resp, float *l2ressh, float *l2resrayl, float residualp, float residualsh, float residualrayl, float normresidualp, float normresidualsh, float normresidualrayl, float finalresidual);
float **alloc2d(int m, int n) ;
double **alloc2ddouble(int m, int n);
void free2d(float **p);
void free2ddouble(double **p);
void Search_Cross_Correlation(float *U0,float *U1,int N0,float Delta0,int N1,float T_shift,float dT_Search,float * T_Cross,float *Cross_Corre_Coefficient,float *Scale_Factor,float *Rms_Factor);
void        design        (int        iord,
                           char      *type,
                           char      *aproto,
                           double     a,
                           double     trbndw,
                           double     fl,
                           double     fh,
                           double     ts,
                           float     *sn,
                           float     *sd,
                           int       *nsects);

void        apply         (float     *data,
                           int        nsamps,
                           int        zp,
                           float     *sn,
                           float     *sd,
                           int        nsects);
