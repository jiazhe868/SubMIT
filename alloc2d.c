#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sub_header.h"

float **alloc2d(int m, int n) {
        int i;
        float **p;
        p=(float **)malloc(sizeof(float *)*m);
        p[0]=(float *)malloc(sizeof(float)*m*n);
        for (i=1;i<m;i++) {
                p[i]=p[i-1]+n;
        }
        return p;
}

double **alloc2ddouble(int m, int n) {
        int i;
        double **p;
        p=(double **)malloc(sizeof(double *)*m);
        p[0]=(double *)malloc(sizeof(double)*m*n);
        for (i=1;i<m;i++) {
                p[i]=p[i-1]+n;
        }
        return p;
}

void free2d(float **p) {
        free(p[0]);
        free(p);
}

void free2ddouble(double **p) {
        free(p[0]);
        free(p);
}
