#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
  Frenet-Serret scalar/vector potential evaluator for constant h != 0.

  Convention used here:
      B = -grad(phi) = -(d_x phi, d_y phi, (1/q) d_s phi),  q = 1 + h x.

  Input polynomial layout:
      bs[k] = coefficient of s^k in b_s(s) = a_0'(s), k=0..deg.
      a[(n-1)*(deg+1)+k] = coefficient of s^k in a_n(s), n=1..na.
      b[(n-1)*(deg+1)+k] = coefficient of s^k in b_n(s), n=1..nb.

  The builder stores enough y-coefficients to evaluate Phi, B and A through
  order ny in y, with B_y also accurate through order ny.
*/

#define FS_DELTA_RELATIVISTIC 0
#define FS_DELTA_EQUALS_PTAU 1

/* -- STRUCTS -- */

/* Keeps recursion relation coefficients c[i,m,k] and their relevant quantities */
typedef struct {
    int ny;      /* requested output order in y */
    int ncoef;   /* stored phi_i coefficients: 0..ny+1 */
    int na, nb, deg;
    int mmin, mmax, moff, nm;
    int qemin, nq;
    double h;
    double *c;   /* c[i,m,k], polynomial coeff of s^k in q^m term */
    double *V;   /* scratch: c[i,m](s)   */
    double *D1;   /* scratch: d_s c[i,m]  */
    double *D2;   /* scratch: d2_s c[i,m]  */
    double *Q;   /* scratch: q^e, e=qemin.. */
} FSField;

/* Keeps value of the scalar potential, field, vector potential at given position */
typedef struct {
    double phi;
    double Bx, By, Bs;
    double Ax, Ay, As;
    double dAx_dx, dAx_dy, dAx_ds;
    double dAs_dx, dAs_dy, dAs_ds;
} FSValue;

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


/* -- HELPER FUNCTIONS -- */

static void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n, sz);
    if (!p) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    return p;
}

static inline size_t coff(const FSField *f, int i, int j, int k) {
    return (((size_t)i * (size_t)f->nm + (size_t)j) * (size_t)(f->deg + 1) + (size_t)k);
}

/* pointer to coefficient */
static inline double *cptr(FSField *f, int i, int j) {
    return f->c + (((size_t)i * (size_t)f->nm + (size_t)j) * (size_t)(f->deg + 1));
}

static inline const double *ccptr(const FSField *f, int i, int j) {
    return f->c + (((size_t)i * (size_t)f->nm + (size_t)j) * (size_t)(f->deg + 1));
}

/* Horner evaluation of p(s) and p'(s) for ascending coeffs p[k] s^k. 
:const double *p: a pointer to polynomial coefficients, p[0] = constant term, p[1] = coeff of s, ..., p[deg] = coeff of s^deg
:int deg: degree of the polynomial
:double s: the point at which to evaluate the polynomial and its derivative
:double *v: pointer where the value of the polynomial at s will be stored
:double *dv: pointer where the value of the derivative of the polynomial at s will be stored */
static inline void poly_eval_d2(const double *p, int deg, double s, double *v, double *d1, double *d2) {
    double a = p[deg], b = 0.0, c=0.0;
    for (int k = deg - 1; k >= 0; --k) {
        c = c * s + 2.0 * b;
        b = b * s + a;
        a = a * s + p[k];
    }
    *v = a;
    *d1 = b;
    *d2 = c;
}


/* -- FIELD EXPANSION EVALUATION -- */

/* Will return a pointer to an FSField struct 
:double h: curvature
:int ny: max power in y to keep
:int na, nb: number of a, b coefficients
:int deg: degree of the polynomials describing a, b
:const double *a, *b: pointers to the first element of a, b polynomial coefficients */
FSField *fs_build(double h, int ny, int na, int nb, int deg,
                  const double *bs,
                  const double *a, 
                  const double *b) {
    if (h == 0.0) {
        fprintf(stderr, "fs_build: h=0 needs the Cartesian limit; this evaluator assumes h != 0.\n");
        return NULL;
    }

    FSField *f = (FSField *)xcalloc(1, sizeof(*f));
    f->h = h;
    f->ny = ny;
    f->ncoef = ny + 2;      /* store phi_0..phi_{ny+1} so By is also order ny */
    f->na = na;
    f->nb = nb;
    f->deg = deg;
    f->mmax = (na > nb - 1) ? na : (nb - 1);  /* max power in q^m, so max power in x => either an x^n or bn x^n-1*/
    f->mmin = -2 * ((f->ncoef - 1) / 2);  /* minimal m, starting from c[0,2] and using recursion relation until c[ncoeff,mmin] */
    f->moff = -f->mmin;
    f->nm = f->mmax - f->mmin + 1;  /* number of m values */
    f->qemin = f->mmin - 1;
    f->nq = (f->mmax + 2) - f->qemin + 1;
    f->c = (double *)xcalloc((size_t)f->ncoef * (size_t)f->nm * (size_t)(deg + 1), sizeof(double));
    f->V = (double *)xcalloc((size_t)f->ncoef * (size_t)f->nm, sizeof(double));
    f->D1 = (double *)xcalloc((size_t)f->ncoef * (size_t)f->nm, sizeof(double));
    f->D2 = (double *)xcalloc((size_t)f->ncoef * (size_t)f->nm, sizeof(double));
    f->Q = (double *)xcalloc((size_t)f->nq, sizeof(double));

    int nmax = (na > nb) ? na : nb;
    double *invfact = (double *)xcalloc((size_t)nmax + 1, sizeof(double));  /* 1/n! */
    double *invhpow = (double *)xcalloc((size_t)nmax + 1, sizeof(double));  /* 1/h^n*/
    invfact[0] = 1.0;
    invhpow[0] = 1.0;
    for (int n = 1; n <= nmax; ++n) {
        invfact[n] = invfact[n - 1] / (double)n;
        invhpow[n] = invhpow[n - 1] / h;
    }

    /* phi_0(s) = sum_m c[0,m](s) q^m
    c[0,m] = - sum_(n>=m) (-1)^(n-m) / (h^n m! (n-m)!) a_n(s) 
    */
    for (int m = 0; m <= na; ++m) {
        double *dst = cptr(f, 0, m + f->moff);
        /* a_0(s)=int_0^s b_s(u)du contributes -a_0 to phi_0. 
        CAREFUL: this will neglect the highest order in the polynomial, 
        only up to given degree in a0 is kept */
        if (m==0) {
        for (int k = 0; k < deg; ++k)
            dst[k + 1] = bs[k] / (double)(k + 1);
        }
        for (int n = (m > 1 ? m : 1); n <= na; ++n) {
            double sgn = ((n - m) & 1) ? -1.0 : 1.0;
            double fac = -sgn * invhpow[n] * invfact[m] * invfact[n - m];
            const double *an = a + (size_t)(n - 1) * (size_t)(deg + 1);
            for (int k = 0; k <= deg; ++k) dst[k] += fac * an[k];
        }
    }
    
    /* phi_1(q,s) = sum_m c[1,m](s) q^m
    c[1,m] = - sum_(n>=m+1) (-1)^(n-1-m) / (h^(n-1) m! (n-1-m)!) b_n(s) */
    if (f->ncoef > 1) {
        for (int m = 0; m <= nb - 1; ++m) {
            double *dst = cptr(f, 1, m + f->moff);
            for (int n = m + 1; n <= nb; ++n) {
                double sgn = ((n - 1 - m) & 1) ? -1.0 : 1.0;
                double fac = -sgn * invhpow[n - 1] * invfact[m] * invfact[n - 1 - m];
                const double *bn = b + (size_t)(n - 1) * (size_t)(deg + 1);
                for (int k = 0; k <= deg; ++k) dst[k] += fac * bn[k];
            }
        }
    }

    /* Recursion: c[i+2,m] = -(d_s^2 + h^2 (m+2)^2) c[i,m+2] 
    implemented for polynomial expansion of c[i,m] in powers of s 
    C[i+2,m,k] = -(C[i,m+2,k+2]*(k+2)*(k+1) + C[i,m+2,k]*h^2*(m+2)^2)
    */
    for (int i = 0; i + 2 < f->ncoef; ++i) {
        for (int m = f->mmin; m <= f->mmax - 2; ++m) {
            const double *src = ccptr(f, i, (m + 2) + f->moff);  /* memory location of C[i,m+2,0], taking into account that the minimal value is not zero by moff */
            double *dst = cptr(f, i + 2, m + f->moff);  /* memory location of C[i+2,m,0] */
            double lam = h * h * (double)(m + 2) * (double)(m + 2);
            for (int k = 0; k <= deg; ++k) {
                double v = lam * src[k];
                if (k + 2 <= deg) v += (double)(k + 2) * (double)(k + 1) * src[k + 2];
                dst[k] = -v;
            }
        }
    }

    free(invfact);
    free(invhpow);
    return f;
}

void fs_free(FSField *f) {
    if (!f) return;
    free(f->c);
    free(f->V);
    free(f->D1);
    free(f->D2);
    free(f->Q);
    free(f);
}

/* Evaluate all c[i,m](s) and d_s c[i,m](s)
and keep them in V, D */
static void fs_prepare_s(FSField *f, double s) {
    for (int i = 0; i < f->ncoef; ++i) {
        for (int j = 0; j < f->nm; ++j) {
            poly_eval_d2(ccptr(f, i, j), f->deg, s,
                         &f->V[i * f->nm + j],
                         &f->D1[i * f->nm + j],
                         &f->D2[i * f->nm + j]);
        }
    }
}

/* Evaluate fields
:FSField *f: pointer to field to evaluate
:double x, y, s: positions at which to evaluate field
:FSValue *out: pointer to object where to keep evaluation, 
keeping scalar potential, vector potential, fields */
int fs_eval(FSField *f, double x, double y, double s, FSValue *out) {
    const double q = 1.0 + f->h * x;
    if (q == 0.0) return -1; /* singular chart */

    memset(out, 0, sizeof(*out));
    out->Ay = 0.0;

    const int nv = f->ncoef * f->nm;
    double *V = f->V;
    double *D1 = f->D1;
    double *D2 = f->D2;
    double *Q = f->Q;

    fs_prepare_s(f, s);

    /* q powers from e = mmin-1 .. mmax+2 */
    Q[0] = pow(q, (double)f->qemin);
    for (int t = 1; t < f->nq; ++t) Q[t] = Q[t - 1] * q;
#define QPOW(E) Q[(E) - f->qemin]

    /* As(x,0,s) 
    = 1/(1+hx) int_0^x dx' *(1+hx) By(x',0,s) 
    = 1/qh int_1^q dq' q' phi_1(q',s)
    = 1/qh sum_m c[1,m] q^(m+2)/(m+2) - 1/q sum_m c[1,m] 1/(m+2) */
    if (f->ncoef > 1) {
        for (int m = 0; m <= f->mmax; ++m) {
            int j = m + f->moff;
            const double c1 = V[1 * f->nm + j];  /* c[1,m] */
            const double d1 = f->D1[1 * f->nm + j];
            const double g = QPOW(m + 1) - QPOW(-1);
            const double den = f->h * (double)(m + 2);
            if (c1 != 0.0) {
                out->As += c1 * g / den;
                out->dAs_dx += c1 * (((double)(m + 1)) * QPOW(m) + QPOW(-2)) / (double)(m + 2);
            }
            if (d1 != 0.0) out->dAs_ds += d1 * g / den;            
        }
    }

    double t = 1.0; /* y^i / i! */
    for (int i = 0; i <= f->ny; ++i) {
        double sphi = 0.0, gx = 0.0, gs = 0.0, gy = 0.0;
        double dgx_dx = 0.0, dgx_ds = 0.0;
        double dgs_dx = 0.0, dgs_ds = 0.0;

        for (int m = f->mmin; m <= f->mmax; ++m) {
            const int j = m + f->moff;
            const double v = V[i * f->nm + j];      /* c[i,m]' */
            const double d1 = D1[i * f->nm + j];    /* c[i,m] */
            const double d2 = D2[i * f->nm + j];    /* c[i,m] */
            const double qm = QPOW(m);              /* q^m */
            const double qm1 = QPOW(m - 1);         /* q^(m-1) */
            const double qm2 = QPOW(m - 2);         /* q^(m-2) */

            sphi += v * qm;                         /* c[i,m] q^m */
            gx   += f->h * (double)m * v * qm1;     /* h m c[i,m] q^(m-1) */
            gy   += V[(i + 1) * f->nm + j] * qm;    /* c[i+1,m] q^m */
            gs   += d1 * qm1;                       /* c[i,m]' q^(m-1) */

            dgx_dx += f->h * f->h * (double)m * (double)(m-1) * v * qm2;
            dgx_ds += f->h * (double)m * d1 * qm1;
            dgs_dx += f->h * (double)(m-1) * d1 * qm2;
            dgs_ds += d2 * qm1;
        }

        out->phi += sphi * t;  /* c[i,m] q^m y^i/i!*/
        out->Bx  -= gx   * t;  /* -h m c[i,m] q^(m-1) y^i/i! */
        out->By  -= gy   * t;  /* -c[i+1,m] q^m y^i/i! */
        out->Bs  -= gs   * t;  /* -c[i,m]' q^(m-1) y^i/i! */

        /* A_x, A_s through order ny in y: need i = 0..ny-1 */
        if (i < f->ny) {
            double u = t * y / (double)(i + 1);  /* y^(i+1)/(i+1)! */
            out->Ax += gs * u;  /* -c[i,m]' q^(m-1) y^(i+1)/(i+1)! */
            out->As -= gx * u;  /* -h m c[i,m] q^(m-1) y^i/i! *//* -h m c[i,m] q^(m-1) y^i/i! */
            out->dAx_dx += dgs_dx  * u;
            out->dAx_ds += dgs_ds  * u;
            out->dAs_dx += -dgx_dx * u;
            out->dAs_ds += -dgx_ds * u;
        }

        t *= y / (double)(i + 1);
    }

    out->dAx_dy = -out->Bs;
    out->dAs_dy =  out->Bx;

    (void)nv;
    return 0;
#undef QPOW
}


/* -- HAMILTONIAN -- */

void fs_hamiltonian_params_default(FSHamiltonianParams *p, double beta0) {
    if (!p) return;
    p->beta0 = beta0;
    p->delta_mode = FS_DELTA_RELATIVISTIC;
}

static FSHamiltonianParams fs_params_value(const FSHamiltonianParams *p) {
    FSHamiltonianParams q;
    fs_hamiltonian_params_default(&q, 1.0);
    if (p) q = *p;
    if (!(q.beta0 > 0.0)) q.beta0 = 1.0;
    return q;
}

static int fs_delta_from_ptau(const FSHamiltonianParams *params, double ptau,
                              double *delta, double *delta1, double *ddelta1) {
    if (params->delta_mode == FS_DELTA_EQUALS_PTAU) {
        *delta = ptau;
        *delta1 = 1.0 + ptau;  // delta+1
        *ddelta1 = 1.0;  // d(delta+1)/dptau
    }{
        const double r = 1.0 + 2.0 * ptau / params->beta0 + ptau * ptau;
        *delta1 = sqrt(r);
        *delta = *delta1 - 1.0;
        *ddelta1 = (1.0 / params->beta0 + ptau) / (*delta1);
    }
}


int fs_hamiltonian_flow(FSField *f, const FSHamiltonianParams *params_in,
                        double s, const double z[6], FSHamiltonianFlow *flow) {
    FSHamiltonianParams params = fs_params_value(params_in);
    double delta1, delta, ddelta1;
    double q, pix, piy, rad, root;

    memset(flow, 0, sizeof(*flow));

    fs_eval(f, z[0], z[2], s, &flow->pot);
    fs_delta_from_ptau(&params, z[5], &delta, &delta1, &ddelta1);

    q = 1.0 + f->h * z[0];
    pix = z[1] - flow->pot.Ax;
    piy = z[3] - flow->pot.Ay;  /* A_y is zero in this gauge. */
    rad = delta1 * delta1 - pix * pix - piy * piy;
    root = sqrt(rad);

    flow->delta = delta;
    flow->one_plus_delta = delta1;
    flow->radicand = rad;
    flow->root = root;
    flow->H = z[5] / params.beta0 - q * (root + flow->pot.As);

    flow->rhs[0] = q * pix / root;  // dH/dpx
    flow->rhs[2] = q * piy / root;  // dH/dpy
    flow->rhs[4] = 1.0 / params.beta0 - q * delta1 * ddelta1 / root;  // dH/dptau

    flow->rhs[1] = f->h * (root + flow->pot.As)
        + q * (pix * flow->pot.dAx_dx / root + flow->pot.dAs_dx);  // -dH/dx
    flow->rhs[3] = q * (pix * flow->pot.dAx_dy / root + flow->pot.dAs_dy);  // -dH/dy
    flow->rhs[5] = 0.0;  /* -dH/dtau, H has no tau-dependence for these static fields. */

    flow->grad[0] = -flow->rhs[1];
    flow->grad[1] =  flow->rhs[0];
    flow->grad[2] = -flow->rhs[3];
    flow->grad[3] =  flow->rhs[2];
    flow->grad[4] = -flow->rhs[5];
    flow->grad[5] =  flow->rhs[4];
    flow->dH_ds = -q * (pix * flow->pot.dAx_ds / root + flow->pot.dAs_ds);
}