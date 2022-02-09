import matplotlib.pyplot as plt
from math import exp,log,sinh,sqrt
from scipy import special
import sympy as sp
import numpy as np
import cmath


def as4(Vgs,Vds):
    q = 1.6e-19
    Na = 2e18
    Eo = 8.85e-12
    Eox = 3.9 * Eo
    Esi = 11.8 * Eo
    tox = 2e-9
    K = 1.38e-23
    x_s = 4.17
    T = 300
    Nd = 1e26
    E_g_eff = 1.12
    phi_m = 4.35
    BGN1 = 0.5
    BGN2 = 1e23
    BGN3 = 19e-3
    ni = 1.45e16
    Ld = 50e-9
    Rd = Ld / 2
    Lch = 32e-9
    Rs = Rd + Lch

    Eg1 = BGN3 * (log(Na / BGN2) + pow((pow(log(Na / BGN2), 2) + BGN1), 0.5))
    E_g = E_g_eff - Eg1
    Ni = ni * pow(exp((Eg1 * q) / (K * T)), 0.5)

    Vbi = (K * T / q) * (log(Na * Nd / Ni ** 2))
    phi_c = (K * T / q) * (log(Na / Ni))
    Vfb = phi_m - (x_s + (E_g / 2) + phi_c)

    def phi_gs1(Vgs):
        a = Vgs - Vfb
        return a

    Wd = 17e-09
    Cox = Eox / (Wd * log(1 + tox / Wd))

    def K1():
        a1, b1, v1 = (Cox * tox) / Esi, (1 + (tox / Wd)), tox ** 2
        a2, b2 = Cox * v1, 2 * Esi * Wd
        return (a1 * (b1) + (a2 / b2))

    K2 = (Cox / Esi) * (1 + (tox / Wd))
    K3 = Cox / (2 * Esi * Wd)
    var = tox ** 2
    K4 = 1 - K1() + (K2 * tox) - K3 * var
    A = (-2 * K3) / K4

    def B(Vgs):
        v1, v2 = (q * Na) / Esi, 2 * K3 * phi_gs1(Vgs)
        a = (1 / K4) * (v1 - v2)
        return a

    # print(Cox,K4,K3)

    r01, r02 = Ld / 2, (Ld / 2) + Lch
    r11 = cmath.sqrt(A)

    v01, v02 = r01 * r11, r02 * r11

    def j0():
        alz = v01
        J0 = special.jv(0, alz)
        # J0 = special.jv(0, alz)
        return J0

    Jo_r1, Jo_r2 = special.jv(0, v01), special.jv(0, v02)
    Yo_r1, Yo_r2 = special.yv(0, v01), special.yv(0, v02)

    # print(A , r01 , r11 , v01 ,   j0())
    # print(Jo_r1,Jo_r2,Yo_r1,Yo_r2)

    def P1():
        a = ((Jo_r1 * Yo_r2) - (Jo_r2 * Yo_r1))
        return a

    def C1(Vgs, Vds):
        up1, up2, up3, up4 = (Vbi + Vds) * Yo_r2, Vbi * Yo_r1, B(Vgs) / A, (Yo_r2 - Yo_r1)
        up = (up1 - up2 - (up3 * up4))
        a = up / P1()
        return a

    P2 = special.yv(0, r02 * r11)

    def C2(Vgs, Vds):
        v1, v2 = Jo_r2 * C1(Vgs, Vds), B(Vgs) / A
        up = ((Vbi - v1) - v2) / (P2)
        return up

    # print(C1(0,0),C2(0,0))

    def phi_s(r, Vgs, Vds):
        v1, v2, v3 = special.jv(0, r * (r11)), special.yv(0, r * (r11)), B(Vgs) / A
        a = ((C1(Vgs, Vds) * v1) + (C2(Vgs, Vds) * v2) + (v3))
        return a

    def A3(r, Vgs, Vds):
        v1, v2 = phi_gs1(Vgs) - phi_s(r, Vgs, Vds), Cox / (Esi * Wd)
        a = v1 * v2
        return a

    def A2(r, Vgs, Vds):
        v1, v2 = phi_gs1(Vgs) - phi_s(r, Vgs, Vds), -Cox / (Esi)
        a = v1 * v2 * (1 + (tox / Wd))
        return a

    def A1(r, Vgs, Vds):
        v01 = tox ** 2
        v1, v2 = phi_s(r, Vgs, Vds) - (A2(r, Vgs, Vds) * tox), A3(r, Vgs, Vds) * v01
        return (v1 - v2)

    def phi(r, z, Vgs, Vds):
        v01 = z ** 2
        v1, v2 = A1(r, Vgs, Vds) + (A2(r, Vgs, Vds) * z), A3(r, Vgs, Vds) * v01
        return (v1 - v2 - phi_c)

    def phi_1(r, Vgs, Vds):
        return (phi_s(r, Vgs, Vds) - phi_c)

    # def EF(r,Vgs,Vds)://////////////////////////////////

    def Wa():
        v1, v2 = r01 ** 2, r02 ** 2
        a = 3.14 * (v2 - v1)
        return a

    xx1 = 2 * 3.14 * (r02)

    ap = np.arange(25e-9, 58e-9, 1e-9)

    I = []

    for r in ap:
        I.append(phi(r, 2e-9, Vgs, Vds))

    ax = plt.plot(I, marker='', color='black', linewidth=3, linestyle='dashed', label="D_Gate")
    plt.legend()

    return  ax
#plt.xlim([0,40e-9])









#I1 = I1[:: -1]
#plt.xlabel('Lc,(m) --> ')
#plt.ylabel(' Surface potential --> ')
#plt.title('Double Gate MOSFET device')
#plt.show()'''