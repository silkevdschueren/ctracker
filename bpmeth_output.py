import bpmeth


x, y, s = 0.01, 0.1, 0.
h = 0.2
a = (0,)#("1 + 0.1*s", 0.2, 0.3)
b = (1,)#0.5)
bs = 0#.1

fexp = bpmeth.FieldExpansion(a=a, b=b, h=h, bs=bs)
phi = fexp.get_phi()
print(f"phi = {phi.subs({fexp.x: x, fexp.y: y, fexp.s: s})}")

B = fexp.get_Bfield()
print(f"B = {B[0](x, y, s), B[1](x, y, s), B[2](x, y, s)}")

A = fexp.get_A(lambdify=True)
print(f"A = {A[0](x, y, s), A[1](x, y, s), A[2](x, y, s)}")


px = 0.01
py = 0.02
tau = 0.
ptau = 0.03
beta0=1.


ham = bpmeth.Hamiltonian(1, h, fexp)
vf = ham.get_vectorfield()
vf_eval = vf(s, (x, y, tau, px, py, ptau), beta0)  # xdot, ydot, taudot, pxdot, pydot, ptaudot

print("Hamiltonian values:")
print(f"dH/dx = {-vf_eval[3]}")    # -pxdot
print(f"dH/dpx = {vf_eval[0]}")    #  xdot
print(f"dH/dy = {-vf_eval[4]}")    # -pydot
print(f"dH/dpy = {vf_eval[1]}")    #  ydot
print(f"dH/dtau = {-vf_eval[5]}")  # -ptaudot
print(f"dH/dptau = {vf_eval[2]}")  #  taudot


