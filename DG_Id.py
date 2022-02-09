import matplotlib.pyplot as plt
from math import exp,log,sinh,sqrt
from scipy import special
import sympy as sp
import numpy as np
import cmath
from scipy import special
from scipy.integrate import trapz

def DG_I():
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
    # print(C2(1,0))

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

    # def EF(r,Vgs,Vds)://////////////////////////////////

    def Wa():
        v1, v2 = r01 ** 2, r02 ** 2
        a = 3.14 * (v2 - v1)
        return a

    xx1 = 2 * 3.14 * (r02)

    Ra = Wa() / xx1

    Vth1 = 0.203
    Vth12 = 0.0954
    a1, a2 = 0.2, 0.001
    theta = 0.5
    Vn = 1.03e5
    vs = Vn
    mu_n = 180e-4
    Voff = 0.1
    Keff = 0.5

    Ec = Vn / mu_n
    theta_short = -0.051
    beta = 10
    Vsb = 0

    def Vdsat(Vgs):
        v1 = Vgs - Vth1
        v2 = v1 / (Lch * Ec)
        return v1 / (1 + v2)

    def alpha_n(Vgs, Vds):
        v1, v2 = 0.5 * Keff, (Vsb) + phi_s(Lch / 2, Vgs, Vds)
        a01, a02 = (0.5 * Keff) / (cmath.sqrt(v2)), a1 + a2 * (v2)
        return 1 + a01 * (1 - (1 / a02))

    def mu_s1(r, Vgs, Vds):
        v1, v2, v3 = theta * (Vgs - Vth1), 2 * Keff * theta, cmath.sqrt(phi_s(r, Vgs, Vds) + Vsb)
        return (mu_n / (1 + v1 + v2 * v3))

    def b():
        v1, v2 = Esi * Wd * tox, 2 * Eox
        a = sqrt(v1 / v2)
        return a

    lamda = 60e-11
    k = 0.55
    ks = 0.55

    def Lsat(Vds, Vgs):
        v1 = (Vds - Vdsat(Vgs)) / (b() * Ec)
        v2 = cmath.sqrt((v1 ** 2) + 1)
        a = 0.6 * b() * (np.log(v1 + v2))
        return a

    def Idlin(Vgs, Vds):
        v1, v2, v3 = Ra * Cox, Lch - Lsat(Vgs, Vds), (Vgs - Vth1)
        a1 = (mu_n * v1) / (v2 + (Vds / Ec))
        a2 = (lamda * v1) / (v2 ** 2)
        a3 = ((k * v3) - 0.5 * Vds) * Vds
        return (a1 + a2) * (a3)

    # print(Idlin(1,0.05))
    # print(Wa() , xx1)

    # print(Idlin(1,0.05))

    def Idsat(Vgs, Vds):
        v1, v2, v3 = Ra * Cox, Lch - Lsat(Vgs, Vds), (Vgs - Vth1)
        a1 = (mu_n * v1) / (v2 + (Vdsat(Vgs) / Ec))
        a2 = (lamda * v1) / (v2 ** 2)
        a3 = ((k * v3) - 0.5 * Vdsat(Vgs)) * Vdsat(Vgs)
        return (a1 + a2) * (a3)

    # print(Idsat(0,0))

    Vgs_range = np.arange(0, 0.3, 0.02)

    r_range = np.arange(25e-9, 58e-9, 1e-9)
    z = np.arange(0, Wd, 16e-9)

    vgs = sp.symbols('vgs')

    def Id_sub(vgs):
        phi_at_different_z = []
        int1_at_diffrent_z = []
        term1 = q / (K * T)
        int1_at_diffrent_z = []
        for R in r_range:
            phi_at_different_z = []

            for z_val in z:
                phi_at_different_z.append(phi(R, z_val, vgs, 0.05))

            final_term = np.exp(term1 * np.array(phi_at_different_z))

            int1 = trapz(final_term, z)

            int1_at_diffrent_z.append(int1)

        term1_for_int2 = 1 / (np.array(int1_at_diffrent_z))

        int2 = trapz(term1_for_int2, r_range)

        def temp1(Vds):
            c = (-q * Vds) / (K * T)
            return c

        temp2 = 1 - exp(temp1(0.05))

        Id_sub1 = (mu_n * Ra * K * T * Ni * temp2) / np.array(int2)  # mu_n*Ra*K*T*Ni(1-(np.exp((-q*Vds)/K*T)))
        return Id_sub1

    Vgs1 = np.arange(0, 1, 0.02)
    Id = []
    for Vgs in Vgs1:
        if Vgs < Vth1 + 0.02:
            Id.append(Id_sub(Vgs))
        elif Vgs >= Vth1 + 0.02:
            if 0.05 <= Vdsat(Vgs):
                Id.append(Idlin(Vgs, 0.05))
        elif Vgs >= Vth1 + 0.02:
            if 0.05 > Vdsat(Vgs):
                Id.append(Idsat(Vgs, 0.05))

    Vgs1 = Vgs1[:-2]
    # Vgs1 = Vgs1[:-1]
    ax = plt.plot(Vgs1, Id , marker='', color='green', linewidth=3, linestyle='dashed', label="D_Gate")
    plt.legend()
    return ax

















#I1 = I1[:: -1]
#plt.xlabel('Lc,(m) --> ')
#plt.ylabel(' Surface potential --> ')
#plt.title('Double Gate MOSFET device')
#plt.show()'''