import numpy as n
import matplotlib.pyplot as plt
from math import *
import os,time


'''def A2():
    v1 , v2 = 0.5*q*Na*tsi*tox , 0.5*q*Na*(tsi - y)*y
    return phi_gs - v1/Eox - v2/Eox



def lmda2(x):
    a1 ,b1  = Esi*tox*tsi , 2*Eox
    a2  , a3 = 0.5*x*tsi ,x**2

def shi(x,y):
    v1 , v2 ,v3 = sinh(L/lmda2(x)) , Vbi + Vds - A2() , (L - x)/lmda2(x)

    return  (1/v3)*(v2)*sinh(x/lmda2(x))*(Vbi - A2())*sinh(v3)*v1'''


def as3(Vgs,Vds):
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
    ni = 1.45e-9

    Eg1 = BGN3 * (log(Na / BGN2) + pow((pow(log(Na / BGN2), 2) + BGN1), 0.5))
    E_g = E_g_eff - Eg1
    Ni = ni * pow(exp((Eg1 * q) / (K * T)), 0.5)

    Vbi = (K * T / q) * (log(Na * Nd / Ni ** 2))
    phi_c = (K * T / q) * (log(Na / Ni))
    Vfb = phi_m - (x_s + (E_g / 2) + phi_c)

    def phi_gs(Vgs):
        phi_gs1 = Vgs - Vfb
        return phi_gs1

    tsi = 5e-9

    gma = Esi / Eox
    lmda = sqrt(gma * tox * tsi)
    x = n.arange(32e-9, 0e-9, -1e-9)

    L = 33e-9
    I = []
    for x1 in x:
        a1 = sinh(L / lmda)
        a2 = q * Na / Esi
        a3 = lmda ** 2
        a4 = (L - x1) / lmda
        a5 = Vbi - phi_gs(Vgs)
        a6 = sinh(x1 / lmda)
        U = 1 / a1 * ((a5 + a2 * a3) * sinh(a4) + (Vds + a5 + a2 * a3) * a6)
        I.append(U)

    #I = I[:: -1]
    ax = plt.plot(I,  marker='o', markerfacecolor='blue', markersize=2, color='blue', linewidth=3, label = 'B_Mosfet')
    plt.xlabel('Lc,(m) --> ')
    plt.ylabel(' Surface potential --> ')
    plt.title('Bulk MOSFET device')
    plt.legend()
    return ax
