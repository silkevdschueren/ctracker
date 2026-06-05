#ifndef fs_eval_H
#define fs_eval_H

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
} Expansion;

/* output struct */
typedef struct {
    double phi;
    double Bx, By, Bs;
    double Ax, Ay, As;
    double dAx_dx, dAx_dy, dAx_ds;
    double dAs_dx, dAs_dy, dAs_ds;
} FieldValue;

/* API */
Expansion *build_expansion(double h, int ny, int na, int nb, int deg, 
                  const double *bs,
                  const double *a, 
                  const double *b);

int evaluate_expansion(Expansion *f, double x, double y, double s, FieldValue *out);

void fs_free(Expansion *f);


typedef struct {
    double beta0;
    int delta_mode;
    int newton_max_iter;
    double newton_tol;
    double newton_fd_eps;
} HamiltonianParams;

typedef struct {
    double H;
    double delta;
    double one_plus_delta;
    double radicand;
    double root;
    double grad[6];  /* dH/d{x,px,y,py,tau,ptau} */
    double rhs[6];   /* canonical flow dz/ds */
    double dH_ds;    /* explicit derivative at fixed canonical variables */
    FieldValue pot;
} HamiltonianFlow;

void hamiltonian_params_default(HamiltonianParams *p, double beta0);

int hamiltonian_flow(Expansion *f, const HamiltonianParams *params_in,
                        double s, const double z[6], HamiltonianFlow *flow);

void track_euler(Expansion *f,
                    const HamiltonianParams *params,
                    double z[6],
                    double s0,
                    double ds,
                    int nstep);

void track_rk4(Expansion *f,
                  const HamiltonianParams *params,
                  double z[6],
                  double s0,
                  double ds,
                  int nstep);
                  
void step_gauss_legendre4(Expansion *f,
                            const HamiltonianParams *params,
                            double s0,
                            double ds,
                            const double z0[6],
                            double z1[6]);


void integrate_gauss_legendre4_array(Expansion *f,
                                       const HamiltonianParams *params,
                                       double s0,
                                       double ds,
                                       int nstep,
                                       int ntraj,
                                       double *z);


#endif