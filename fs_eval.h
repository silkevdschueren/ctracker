#ifndef FS_EVAL_H
#define FS_EVAL_H

#include <stddef.h>  /* for size_t if needed */

/* full struct definition */
typedef struct {
    int ny;
    int ncoef;
    int na, nb, deg;
    int mmin, mmax, moff, nm;
    int qemin, nq;
    double h;
    double *c;
    double *V;
    double *D;
    double *Q;
} FSField;

/* output struct */
typedef struct {
    double phi;
    double Bx, By, Bs;
    double Ax, Ay, As;
} FSValue;

/* API */
FSField *fs_build(double h, int ny, int na, int nb, int deg, 
                  const double *bs,
                  const double *a, 
                  const double *b);

int fs_eval(FSField *f, double x, double y, double s, FSValue *out);

void fs_free(FSField *f);

#endif