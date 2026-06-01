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
    b[1 * (deg + 1) + 0] = 0.1;  /* b_2 s^0 */
    const double h = 0.2;
    const double x = 0.01, y = 0.1, s = 0;

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

    fs_free(f);
    return 0;
}
