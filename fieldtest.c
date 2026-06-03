#include <stdio.h>
#include <string.h>

#include "fs_eval.h"

int main(void) {
    /* Example: a_0(s)=0, a_1(s)=1, others 0; b_n=0. */
    const int deg = 1, na = 3, nb = 2, ny = 6;
    double bs[deg + 1];
    double a[na * (deg + 1)];
    double b[nb * (deg + 1)];
    memset(bs, 0, sizeof(bs));
    memset(a, 0, sizeof(a));
    memset(b, 0, sizeof(b));
    a[0 * (deg + 1) + 0] = 1.0;  /* a_1 s^0 */
    a[0 * (deg + 1) + 1] = 0.1;  /* a_1 s^1 */
    a[1 * (deg + 1) + 0] = 0.2;  /* a_2 s^0 */
    a[2 * (deg + 1) + 0] = 0.3;  /* a_3 s^0 */
    b[0 * (deg + 1) + 0] = 1.0;  /* b_1 s^0 */
    b[1 * (deg + 1) + 0] = 0.5;  /* b_2 s^0 */
    bs[0] = 0.1;  /* b_s s^0 BE CAREFUL, CONTRIBUTION ORDER deg NEGLECTED SINCE ORDER IS FOR a0 */
    const double h = 0.2;
    const double x = 0.01, y = 0.02, s = 0.001;

    FSField *f = fs_build(h, ny, na, nb, deg, bs, a, b);
    if (!f) return 1;

    FSValue v;
    if (fs_eval(f, x, y, s, &v) != 0) {
        fprintf(stderr, "q=1+h x hit the chart singularity\n");
        fs_free(f);
        return 1;
    }

    printf("phi = %.15g\n", v.phi);
    printf("B   = (%.15g, %.15g, %.15g)\n", v.Bx, v.By, v.Bs);
    printf("A   = (%.15g, %.15g, %.15g)\n", v.Ax, v.Ay, v.As);


    const double px = 0.005, py = 0.007, tau = 0.008, ptau = 0.004;
    const double z[6] = {x, px, y, py, tau, ptau};

    FSHamiltonianParams hp;
    fs_hamiltonian_params_default(&hp, 1.0);

    FSHamiltonianFlow flow;
    fs_hamiltonian_flow(f, &hp, s, z, &flow);

    printf("Hamiltonian values:\n");
    printf("dH/dx = %.15g\n", flow.grad[0]);
    printf("dH/dpx = %.15g\n", flow.grad[1]);
    printf("dH/dy = %.15g\n", flow.grad[2]);
    printf("dH/dpy = %.15g\n", flow.grad[3]);
    printf("dH/dtau = %.15g\n", flow.grad[4]);
    printf("dH/dptau = %.15g\n", flow.grad[5]);
    

    fs_free(f);
    return 0;
}
