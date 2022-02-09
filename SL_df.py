import matplotlib.pyplot as plt               #////////////// import modules
import numpy as np
from sympy import *
from math import log, exp, sinh
from scipy import special
from scipy.integrate import trapz

           # /////////////////////////////////////////////////    surface _potential vs. position along channel {RingFET}
def L_diff(Vds,Vgs):
    r_range = np.arange(25e-9, 58e-9, 1e-9)  # //////////// constant data values
    # Z_range = np.arange(3e-9, 20e-9, 1e-9)
    Wd_range = np.arange(3e-9, 18e-9, 1e-9)
    phi_m = 4.35
    T = 300
    q = 1.6e-19
    K = 1.38e-23
    ni = 1.45e16
    epsilon_o = 8.85e-12
    epsilon_ox = 3.9 * epsilon_o
    epsilon_si = 11.8 * epsilon_o
    tox = 2e-9
    Ld = 50e-9
    Rd = Ld / 2
    Nd = 1e26
    x_s = 4.17
    E_g_eff = 1.12
    Na = 2e24
    BGN1 = 0.5
    BGN2 = 1e23
    BGN3 = 19e-3
    Eg1 = BGN3 * (log(Na / BGN2) + pow((pow(log(Na / BGN2), 2) + BGN1), 0.5))
    E_g = E_g_eff - Eg1
    # print(Eg1, K, T, q)
    # print((Eg1*q)/(K*T))
    Ni = ni * pow(exp((Eg1 * q) / (K * T)), 0.5)
    # print(Ni)
    Vbi = (K * T / q) * (log(Na * Nd / Ni ** 2))
    phi_c = (K * T / q) * (log(Na / Ni))
    Vfb = phi_m - (x_s + (E_g / 2) + phi_c)

    L = np.arange(30e-9, 90e-9, 15e-9)

    for Lch in L:
        Rs = Rd + Lch

        def phi_gs(Vgs):
            phi_gs1 = Vgs - Vfb
            return phi_gs1

        # print (phi_gs(Vgs))
        phi_f = (K * T / q) * log(Na / Ni)
        # print(phi_f)
        Vsub = 0
        Wd = 18e-09
        Cox = epsilon_ox / (Wd * log(1 + tox / Wd))

        # print(Cox, tox, epsilon_si, Wd)
        def C2(Vgs):
            CC2 = phi_gs(Vgs) - ((q * Na * (Wd ** 2)) / (2 * epsilon_si)) - ((q * Na * Wd) / (Cox))
            return CC2

        # print(C2(Vgs))

        def M(z, Vgs):
            M = (((q * Na) / (2 * epsilon_si)) * z ** 2) + C2(Vgs)
            return M

        # print(M(Z_range,0))

        alpha = 1e8

        def j0(z):
            alz = alpha * z
            J0 = special.jv(0, alz)
            # J0 = special.jv(0, alz)
            return J0

        def j1(z):
            alz = alpha * z
            J1 = special.jv(1, alz)
            return J1

        def Drain(Vgs, Vds, z):
            D = j0(z) * z * (Vbi + Vds - M(z, Vgs))
            return D

        def Source(Vgs, z):
            S = j0(z) * z * (Vbi - M(z, Vgs))
            return S

        # print("Drain",Drain(Vgs,Vds,Z_range))
        # print("Source",Source(Vgs,Z_range))

        # print(Drain(Vgs,Vds,Wd_range),Wd_range)

        W1 = trapz(Drain(Vgs, Vds, Wd_range), Wd_range)

        W11 = (2 / ((Wd ** 2) * ((j1(Wd)) ** 2))) * W1
        # print("w1 is : ", W1,  "\n", W11)

        W2 = trapz(Source(Vgs, Wd_range), Wd_range)

        W22 = (2 / ((Wd ** 2) * ((j1(Wd)) ** 2))) * W2

        # print("HEY", W2)

        def B(Vgs, Vds):
            B1 = (W22 * np.exp(alpha * Rd) - W11 * np.exp(alpha * Rs)) / (2 * sinh(alpha * (Rd - Rs)))
            return B1

        # print("hello", B(Vgs,Vds))

        def A(Vgs, Vds):
            A1 = (W22 - B(Vgs, Vds) * np.exp(-alpha * Rs)) / (np.exp(alpha * Rs))
            return A1

        # print("hi", A(Vgs,Vds))

        def U(r, z, Vgs, Vds):
            alr = alpha * r
            alrn = -1 * alpha * r
            U1 = j0(z) * (A(Vgs, Vds) * np.exp(alr) + B(Vgs, Vds) * np.exp(alrn))
            return U1

        # print(U(r_range,3e-9,Vgs,Vds), r_range)

        def Phi(r, z, Vgs, Vds):
            Phi1 = M(z, Vgs) + U(r, z, Vgs, Vds) - phi_c
            return Phi1

        # print(Phi(r_range,3*10**-9,Vgs,Vds))

        Lch1 = 58e-9 + (Lch - 32e-9 )
        p = Lch * 1e9
        c = '{0}{1}{2}'.format('Lch= ', np.round(p,2), 'nm')
        r_range1 = np.arange(25e-9, Lch1, 1e-9)
        data = []
        for i in r_range1:
            data.append(Phi(i, 17 * 10 ** -9, Vgs, Vds))

        ax = plt.plot(data, marker='', linewidth=3, linestyle='dashed', label=  c )

    plt.legend()

    return ax