import bpmeth
import xtrack as xt
import numpy as np
import sympy as sp

s = 0
x, y, tau = 0.01, 0.02, 0.001
px, py, ptau = 0.0005, 0.0002, 0.0003
beta0=1.

h = 0.2
comp = np.zeros((9, 4), dtype=np.float64)
# comp[2, 0] = 1.0  # a_1 s^0 
# comp[2, 1] = 0.1  # a_1 s^1 
# comp[4, 0] = 0.2  # a_2 s^0 
# comp[6, 0] = 0.3  # a_3 s^0 
comp[1, 0] = 1.0  # b_1 s^0 
# comp[1, 1] = 0.1  # b_1 s^1
# comp[3, 0] = 0.5  # b_2 s^0 
# comp[0, 0] = 0.1  # bs s^0


def comp_to_sympy(comp):
    s = sp.symbols("s")
    out = []
    for row in comp:
        pol = row[0]
        for n, v in enumerate(row[1:]):
            pol += v * s ** (n + 1)
        out.append(pol)
    bs = out[0]
    b = out[1::2]
    a = out[2::2]
    return bs, b, a

bs, b, a = comp_to_sympy(comp)

fexp = bpmeth.FieldExpansion(a=a, b=b, h=h, bs=bs)
# phi = fexp.get_phi()
# print(f"phi = {phi.subs({fexp.x: x, fexp.y: y, fexp.s: s})}")

# B = fexp.get_Bfield()
# print(f"B = {B[0](x, y, s), B[1](x, y, s), B[2](x, y, s)}")

# A = fexp.get_A(lambdify=True)
# print(f"A = {A[0](x, y, s), A[1](x, y, s), A[2](x, y, s)}")


ham = bpmeth.Hamiltonian(0.001, h, fexp, s_start=s)
# vf = ham.get_vectorfield(compile=False)
# vf_eval = vf(s, (x, y, tau, px, py, ptau), beta0)  # xdot, ydot, taudot, pxdot, pydot, ptaudot

# print("Hamiltonian values:")
# print(f"dH/dx = {-vf_eval[3]}")    # -pxdot
# print(f"dH/dpx = {vf_eval[0]}")    #  xdot
# print(f"dH/dy = {-vf_eval[4]}")    # -pydot
# print(f"dH/dpy = {vf_eval[1]}")    #  ydot
# print(f"dH/dtau = {-vf_eval[5]}")  # -ptaudot
# print(f"dH/dptau = {vf_eval[2]}")  #  taudot

pp = xt.Particles(x=x, px=px, y=y, py=py, tau=tau, ptau=ptau, beta0=1)
ham.track(pp)
print("After tracking:")
print(f"x = {pp.x}, px = {pp.px}, y = {pp.y}, py = {pp.py}, zeta = {pp.zeta}, ptau = {pp.ptau}")

pp = xt.Particles(x=x, px=px, y=y, py=py, tau=tau, ptau=ptau, beta0=1)
length = 0.001
    
PolyMagnet = bpmeth.PolySegment(length=length, h=h, comp=comp, s_start=s)
PolyMagnet.track(pp, ds=0.001)
print("After tracking through PolySegment:")
print(f"x = {pp.x}, px = {pp.px}, y = {pp.y}, py = {pp.py}, zeta = {pp.zeta}, ptau = {pp.ptau}")
