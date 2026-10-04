#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sacsubc.h"
#include "sub_header.h"

void convert_staformat_arrays(staformat sstainfo[maxstanum], char sstname[maxstanum][64], float sstdist[maxstanum], float sstaz[maxstanum], float sstcp[maxstanum], float sstcs[maxstanum]) {
	int subi,subj;
	for (subi=0;subi<maxstanum;subi++) {
		for (subj=0;subj<64;subj++){
			sstname[subi][subj]=sstainfo[subi].stname[subj];
		}
		sstdist[subi]=sstainfo[subi].stdist;
		sstaz[subi]=sstainfo[subi].staz;
		sstcp[subi]=sstainfo[subi].stcp;
		sstcs[subi]=sstainfo[subi].stcs;
	}
}

void convert_stalocformat_arrays(stalocformat stainfo[maxlocnum], float stlo[maxlocnum], float stla[maxlocnum], float wt[maxlocnum][3]) {
        int subi,subj;
        for (subi=0;subi<maxlocnum;subi++) {
		wt[subi][0]=stainfo[subi].wte;
		wt[subi][1]=stainfo[subi].wtn;
		wt[subi][2]=stainfo[subi].wtz;
                stlo[subi]=stainfo[subi].stlo;
                stla[subi]=stainfo[subi].stla;
        }
}

void convert_green3dformat_arrays(green3dformat greeninfo3d[maxnum3d], float lo3d[maxnum3d], float la3d[maxnum3d], char ind3d[maxnum3d][64]){
	int subi, subj;
	for (subi=0;subi<maxnum3d;subi++) {
		for (subj=0;subj<64;subj++){
			ind3d[subi][subj]=greeninfo3d[subi].stind[subj];
		}
		lo3d[subi]=greeninfo3d[subi].stlo;
		la3d[subi]=greeninfo3d[subi].stla;
	}
}

void convert_seisformat_arrays(seisdataformat sseisdata[maxstanum], char sstname[maxstanum][64], float sbtime[maxstanum], float st1[maxstanum], float **sxdata) {
	int subi, subj;
	for (subi=0;subi<maxstanum;subi++) {
		for (subj=0;subj<64;subj++) {
			sstname[subi][subj]=sseisdata[subi].stname[subj];
		}
		sbtime[subi]=sseisdata[subi].btime;
		st1[subi]=sseisdata[subi].t1;
		for (subj=0;subj<maxnpts;subj++) {
			sxdata[subi][subj]=sseisdata[subi].xdata[subj];
		}
	}
}

void convert_seislocformat_arrays(seisdatalocformat sseisdata[maxlocnum], char sstname[maxlocnum][64], float sbtime[maxlocnum], float st1[maxlocnum], float **sxdata) {
        int subi, subj;
        for (subi=0;subi<maxlocnum;subi++) {
                for (subj=0;subj<64;subj++) {
                        sstname[subi][subj]=sseisdata[subi].stname[subj];
                }
                sbtime[subi]=sseisdata[subi].btime;
                st1[subi]=sseisdata[subi].t1;
                for (subj=0;subj<maxnpts;subj++) {
			sxdata[subi*3][subj]=sseisdata[subi].edata[subj];
                        sxdata[subi*3+1][subj]=sseisdata[subi].ndata[subj];
                        sxdata[subi*3+2][subj]=sseisdata[subi].zdata[subj];
                }
        }
}

void convert_arrays_staformat(char sstname[maxstanum][64], float sstdist[maxstanum], float sstaz[maxstanum], float sstcp[maxstanum], float sstcs[maxstanum], staformat sstainfo[maxstanum]) {
	int subi,subj;
        for (subi=0;subi<maxstanum;subi++) {
                for (subj=0;subj<64;subj++){
                        sstainfo[subi].stname[subj]=sstname[subi][subj];
                }
                sstainfo[subi].stdist=sstdist[subi];
                sstainfo[subi].staz=sstaz[subi];
                sstainfo[subi].stcp=sstcp[subi];
		sstainfo[subi].stcs=sstcs[subi];
        }
}

void convert_arrays_stalocformat(float sstlo[maxlocnum], float sstla[maxlocnum], float wt[maxlocnum][3], stalocformat sstainfo[maxlocnum]) {
        int subi,subj;
        for (subi=0;subi<maxlocnum;subi++) {
		sstainfo[subi].wte=wt[subi][0];
		sstainfo[subi].wtn=wt[subi][1];
		sstainfo[subi].wtz=wt[subi][2];
                sstainfo[subi].stlo=sstlo[subi];
                sstainfo[subi].stla=sstla[subi];
        }
}

void convert_arrays_green3dformat(float lo3d[maxnum3d], float la3d[maxnum3d], char ind3d[maxnum3d][64], green3dformat greeninfo3d[maxnum3d]) {
	int subi,subj;
	for (subi=0;subi<maxnum3d;subi++) {
		for (subj=0;subj<64;subj++){
			greeninfo3d[subi].stind[subj]=ind3d[subi][subj];
		}
		greeninfo3d[subi].stlo=lo3d[subi];
		greeninfo3d[subi].stla=la3d[subi];
	}
}

void convert_arrays_seisformat(char sstname[maxstanum][64], float sbtime[maxstanum], float st1[maxstanum], float sxdata[][maxnpts], seisdataformat sseisdata[maxstanum]) {
        int subi, subj;
        for (subi=0;subi<maxstanum;subi++) {
                for (subj=0;subj<64;subj++) {
                        sseisdata[subi].stname[subj]=sstname[subi][subj];
                }
                sseisdata[subi].btime=sbtime[subi];
                sseisdata[subi].t1=st1[subi];
                for (subj=0;subj<maxnpts;subj++) {
                        sseisdata[subi].xdata[subj]=sxdata[subi][subj];
                }
        }
}

void convert_arrays_seislocformat(char sstname[maxlocnum][64], float sbtime[maxlocnum], float st1[maxlocnum], float sxdata[][maxnpts], seisdatalocformat sseisdata[maxlocnum]) {
	int subi, subj;
	for (subi=0;subi<maxlocnum;subi++) {
                for (subj=0;subj<64;subj++) {
                        sseisdata[subi].stname[subj]=sstname[subi][subj];
                }
                sseisdata[subi].btime=sbtime[subi];
                sseisdata[subi].t1=st1[subi];
                for (subj=0;subj<maxnpts;subj++) {
                        sseisdata[subi].edata[subj]=sxdata[subi*3][subj];
			sseisdata[subi].ndata[subj]=sxdata[subi*3+1][subj];
			sseisdata[subi].zdata[subj]=sxdata[subi*3+2][subj];
                }
        }
}

