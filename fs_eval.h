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
    double *D1;
    double *D2;
    double *Q;
} FSField;

/* output struct */
typedef struct {
    double phi;
    double Bx, By, Bs;
    double Ax, Ay, As;
    double dAx_dx, dAx_dy, dAx_ds;
    double dAs_dx, dAs_dy, dAs_ds;
} FSValue;

/* API */
FSField *fs_build(double h, int ny, int na, int nb, int deg, 
                  const double *bs,
                  const double *a, 
                  const double *b);

int fs_eval(FSField *f, double x, double y, double s, FSValue *out);

void fs_free(FSField *f);


typedef struct {
    double beta0;
    int delta_mode;
} FSHamiltonianParams;

typedef struct {
    double H;
    double delta;
    double one_plus_delta;
    double radicand;
    double root;
    double grad[6];  /* dH/d{x,px,y,py,tau,ptau} */
    double rhs[6];   /* canonical flow dz/ds */
    double dH_ds;    /* explicit derivative at fixed canonical variables */
    FSValue pot;
} FSHamiltonianFlow;

void fs_hamiltonian_params_default(FSHamiltonianParams *p, double beta0);

int fs_hamiltonian_flow(FSField *f, const FSHamiltonianParams *params_in,
                        double s, const double z[6], FSHamiltonianFlow *flow);

#endif