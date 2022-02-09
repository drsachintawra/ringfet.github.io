import matplotlib.pyplot as plt
from math import log,sqrt
from numpy import *
from sympy import *
import numpy as np
import time,os


def as0():
    # //////////////values/////////////////////////
    ni = 1.45e16
    N1, N2, N3, N4 = 1e26, ni, ni, 1e26
    ts, tr, phi_m, zi = 10e-9, 3e-9, 4.5, 4.17
    q, T, k, Eo, Eg = 1.6e-19, 300, 1.38e-23, 8.854e-12, 1.08
    Lh, Lch = 25e-9, 25e-9
    Es, Er = 11.8 * Eo, 3.9 * Eo

    shi_o, phi_f3 = (-k * T / q) * log(N1 / ni), (k * T / q) * log(N3 / ni)

    def shi_4(Vds0):
        return (k * T / q) * (log(N4 / ni)) + Vds0

    phi_s3 = zi + Eg / 2 + phi_f3
    Vfb = phi_m - phi_s3

    def shi_G3(Vgs0):
        return (Vgs0 - Vfb)

    Cr1, Cr2, Cr3, Cs = (2 / math.pi) * (Er / tr), (Er / tr), (Er / tr), (Es / ts)
    Cr4 = (2 / math.pi) * (Cr1)

    eta_1, eta_2, eta_3, eta_4 = Cr1 / Cs, Cr2 / Cs, Cr3 / Cs, Cr4 / Cs

    def shi_d1(Vgs1):
        v1, v2, v3 = ts ** 2, 2 * eta_1, (q * N1) / Es
        a = (v2 / v1) * shi_G3(Vgs1)
        return (-v1 / v2) * (v3 - a)

    def shi_d2(Vgs2):
        v1, v2, v3 = ts ** 2, 2 * eta_2, (q * N2) / Es
        a = (v2 / v1) * shi_G3(Vgs2)
        return (-v1 / v2) * (v3 - a)

    def shi_d3(Vgs3):
        v1, v2, v3 = ts ** 2, 2 * eta_3, (q * N3) / Es
        a = (v2 / v1) * shi_G3(Vgs3)
        return (-v1 / v2) * (v3 - a)

    def shi_d4(Vgs4):
        v1, v2, v3 = ts ** 2, 2 * eta_4, (q * N4) / Es
        a = (v2 / v1) * shi_G3(Vgs4)
        return (-v1 / v2) * (-v3 - a)

    # //////////////////////////////////
    v1, v2, v3 = ((k * T / q) * log(N1 * N2 / ni ** 2)), 2 * Es * N2, q * N1 * (N1 + N2)
    L1 = sqrt((v1 * v2) / v3)

    def L3(Vgs1, Vds1):
        v1, v2, v3 = (shi_4(Vds1) - shi_d3(Vgs1)), 2 * Es, q * N4
        a = sqrt((v1 * v2) / v3)
        return a

    '''Vd = np.arange(0,1.2,0.2)

    vds = sp.Symbol('vds')
    I = []

    for x in Vd:
        I.append(L3(1,vds))'''

    # plt.plot(Vds , I)
    # plt.show()

    b1 = ts ** 2
    a1, a2, a3, a4 = 2 * eta_1, 2 * eta_2, 2 * eta_3, 2 * eta_4
    lmda_1 = sqrt(a1 / b1)
    lmda_2 = sqrt(a2 / b1)
    lmda_3 = sqrt(a3 / b1)
    lmda_4 = sqrt(a4 / b1)

    def R(Vgs_R):
        return lmda_2 * shi_d2(Vgs_R)

    a_Q, b_Q, c_Q = lmda_1 * coth(lmda_1 * L1), (lmda_2 * exp(lmda_2 * Lh)), sinh(lmda_2 * Lh)

    Q = a_Q - lmda_2 + b_Q / c_Q

    def P(Vgs_P):
        return ((lmda_1 * exp(lmda_1 * L1) * (shi_o - shi_d1(Vgs_P))) - (
                    lmda_1 * coth(lmda_1 * L1) * shi_o * exp(lmda_1 * L1)) -
                lmda_1 * coth(lmda_1 * L1) * shi_d1(Vgs_P) * (1 - exp(lmda_1 * L1)) +
                ((lmda_2 * shi_d2(Vgs_P) * (1 - exp(lmda_2 * Lh))) / (sinh(lmda_2 * Lh))))

    def S(Vgs_S):
        return (-lmda_2 * exp(lmda_2 * Lh) * (R(Vgs_S) / Q + P(Vgs_S) / Q + shi_d2(Vgs_S))) + \
               (lmda_2 * coth(lmda_2 * Lh) * ((R(Vgs_S) * exp(lmda_2 * Lh)) / Q + P(Vgs_S) * exp(lmda_2 * Lh) / Q -
                                              shi_d2(Vgs_S) * (1 - exp(lmda_2 * Lh)))) + \
               (lmda_3 * shi_d3(Vgs_S) * (1 - exp(lmda_3 * Lch))) / (sinh(lmda_3 * Lch)) + lmda_3 * shi_d3(Vgs_S)

    def M():
        v1, v2, v3 = lmda_2 ** 2, lmda_2 * Lh, lmda_3 * Lch  # //////////variable declare////////
        a1, b1 = v1 * exp(v2), Q * sinh(v2)
        a2, b2, a3 = a1 * coth(v2), b1, (lmda_2 * coth(v2) - lmda_3)
        a4, b4 = lmda_3 * exp(v3), sinh(v3)
        return ((a1 / b1) - (a2 / b2) + a3 + (a4 / b4))

    def N(Vgs_N, Vds_N):
        v1, v2, v3 = lmda_3 ** 2, lmda_3 * Lch, lmda_4 * L3(Vgs_N, Vds_N)
        a1, b1 = v1 * exp(v2), M() * sinh(v2)
        a2, b2 = v1 * coth(v2) * exp(v2), M() * sinh(v2)
        a3 = lmda_3 * coth(v2) - lmda_4
        a4, b4 = lmda_4 * exp(v3), sinh(v3)
        return a1 / b1 - a2 / b2 + a3 + a4 / b4

    def L(Vgs_L, Vds_L):
        v1, v2, v3 = lmda_3 * Lch, lmda_4 * shi_d4(Vgs_L), lmda_4 * L3(Vgs_L, Vds_L)
        a011 = S(Vgs_L) / M()
        a01, a02, a04 = a011 + shi_d3(Vgs_L), 1 - exp(v1), 1 - exp(v3)
        a1, a2, a3 = -lmda_3 * exp(v1) * a01, lmda_3 * coth(v1) * shi_d3(Vgs_L) * a02, v2
        a4, b4 = v2 * (1 - exp(v3)), sinh(v3)
        a5, b5 = lmda_4 * shi_4(Vds_L), sinh(v3)
        return a1 - a2 + a3 + a4 / b4 - a5 / b5

    def shi_3(Vgs_shi_3, Vds_shi_3):
        v1 = lmda_3 * Lch
        a1, b1 = -L(Vgs_shi_3, Vds_shi_3), N(Vgs_shi_3, Vds_shi_3)
        a2, b2 = lmda_3 * coth(v1) * exp(v1) * S(Vgs_shi_3), N(Vgs_shi_3, Vds_shi_3) * M()
        return a1 / b1 - a2 / b2

    def shi_2(Vgs_shi_2, Vds_shi_2):
        a1, b1 = lmda_3 * shi_3(Vgs_shi_2, Vds_shi_2), M() * sinh(lmda_3 * Lch)
        a2, b2 = S(Vgs_shi_2), M()
        return a1 / b1 - a2 / b2

    def shi_1(Vgs_shi_1, Vds_shi_1):
        a1, b1 = lmda_2 * shi_2(Vgs_shi_1, Vds_shi_1), Q * sinh(lmda_2 * Lh)
        a2, b2 = R(Vgs_shi_1), Q
        a3, b3 = P(Vgs_shi_1), Q
        return a1 / b1 - a2 / b2 - a3 / b3

    def B4(Vgs_B4, Vds_B4):
        v1, v2 = lmda_4 * L3(Vgs_B4, Vds_B4), lmda_4 * (L1 + Lh + Lch)
        a01, a02 = -(shi_4(Vds_B4) - shi_d4(Vgs_B4)), (shi_3(Vgs_B4, Vds_B4) - shi_d4(Vgs_B4)) * exp(v1)
        a1, b1 = (a01 + a02) * exp(v2), 2 * sinh(v1)
        return a1 / b1

    def A4(Vgs_A4, Vds_A4):
        v1 = lmda_4 * (L1 + Lh + Lch)
        a01, a02 = shi_3(Vgs_A4, Vds_A4) - shi_d4(Vgs_A4), B4(Vgs_A4, Vds_A4) * exp(-v1)
        a1 = (a01 - a02) * exp(-v1)
        return a1

    def B3(Vgs_B3, Vds_B3):
        v1, v2 = lmda_3 * Lch, lmda_3 * (L1 + Lh)
        a011 = shi_2(Vgs_B3, Vds_B3) - shi_d3(Vgs_B3)
        a01, a02 = a011 * exp(v1), shi_3(Vgs_B3, Vds_B3) - shi_d3(Vgs_B3)
        a1, b1 = (a01 - a02) * exp(v2), 2 * sinh(v1)
        return a1 / b1

    def A3(Vgs_A3, Vds_A3):
        v1 = lmda_3 * (L1 + Lh)
        a01, a02 = shi_2(Vgs_A3, Vds_A3) - shi_d2(Vgs_A3), B3(Vgs_A3, Vds_A3) * exp(-v1)
        a1 = (a01 - a02) * exp(-v1)
        return a1

    def B2(Vgs_B2, Vds_B2):
        v1, v2, v3 = lmda_2 * Lh, lmda_2 * L1, lmda_2 * Lh
        a011 = shi_1(Vgs_B2, Vds_B2) - shi_d2(Vgs_B2)
        a01, a02 = a011 * exp(v1), shi_2(Vgs_B2, Vds_B2) - shi_d2(Vgs_B2)
        a1, b1 = (a01 - a02) * exp(v2), 2 * sinh(v3)
        return a1 / b1

    def A2(Vgs_A2, Vds_A2):
        v1 = lmda_2 * L1
        a01, a02 = shi_1(Vgs_A2, Vds_A2) - shi_d2(Vgs_A2), B2(Vgs_A2, Vds_A2) * exp(-v1)
        a1 = (a01 - a02) * exp(-v1)
        return a1

    def B1(Vgs_B1, Vds_B1):
        v1 = lmda_1 * L1
        a011 = shi_o - shi_d1(Vgs_B1)
        a01, a02 = a011 * exp(v1), (shi_1(Vgs_B1, Vds_B1) - shi_d1(Vgs_B1))
        a1, b1 = a01 - a02, 2 * sinh(v1)
        return a1 / b1

    def A1(Vgs_A1, Vds_A1):
        return (shi_d1(Vgs_A1) - shi_d1(Vgs_A1)) - B1(Vgs_A1, Vds_A1)

    def shi_s1(y, Vgs_shi_s1, Vds_shi_s1):
        v1 = lmda_1 * y
        return A1(Vgs_shi_s1, Vds_shi_s1) * exp(v1) + B1(Vgs_shi_s1, Vds_shi_s1) * exp(-v1) + shi_d1(Vgs_shi_s1)

    def shi_s2(y, Vgs_shi_s2, Vds_shi_s2):
        v1 = lmda_2 * y
        return A2(Vgs_shi_s2, Vds_shi_s2) * exp(v1) + B2(Vgs_shi_s2, Vds_shi_s2) * exp(-v1) + shi_d2(Vgs_shi_s2)

    def shi_s3(y, Vgs_shi_s3, Vds_shi_s3):
        v1 = lmda_3 * y
        return A3(Vgs_shi_s3, Vds_shi_s3) * exp(v1) + B3(Vgs_shi_s3, Vds_shi_s3) * exp(-v1) + shi_d3(Vgs_shi_s3)

    def shi_s4(y, Vgs_shi_s4, Vds_shi_s4):
        v1 = lmda_4 * y
        return A4(Vgs_shi_s4, Vds_shi_s4) * exp(v1) + B4(Vgs_shi_s4, Vds_shi_s4) * exp(-v1) + shi_d4(Vgs_shi_s4)

    vds = np.arange(0, 1, 0.1)

    Ln1 = L1 + Lh
    Ln2 = L1 + Lh + Lch



    def shi_t1():
        Vds = 1
        vgs = np.arange(0, 1.4, 0.2)
        y1 = np.arange(-10e-9, 60e-9, 5e-9)
        # print(vgs , y1)
        y2 = y1

        for Vgs in vgs:
            I = []
            for y in y1:
                Vds = 1
                if (y <= 0):
                    I.append(shi_o)
                elif (0 < y <= L1):
                    I.append(shi_s1(y, Vgs, Vds))
                elif (L1 < y <= (Ln1)):
                    I.append(shi_s2(y, Vgs, Vds))
                elif (Ln1 < y <= Ln2):
                    I.append(shi_s3(y, Vgs, Vds))
                elif (Ln2 < y <= (Ln2 + L3(Vgs, Vds))):
                    I.append(shi_s4(y, Vgs, Vds))
                elif (y > (Ln2 + L3(Vgs, Vds))):
                    I.append(shi_4(Vds))

            plt.plot(y2, I)
        plt.xlabel('Lc , (m) --> ')
        plt.ylabel(' Ψs--1(y,Vgs,Vds) --> ')
        plt.title('TFET device {for different Vgs(V)}')
        # plt.suptitle('for different Vgs(v) ')
        plotfile1 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile1)
        res1 = plotfile1
        return res1

    def shi_t2():
        Vgs = 0
        vds1 = np.arange(0, 1.6, 0.2)
        y1 = np.arange(-20e-9, 70e-9, 5e-9)
        # print(vgs , y1)
        y2 = y1

        for Vds in vds1:
            I = []
            for y in y1:
                if (y <= 0):
                    I.append(shi_o)
                elif (0 < y <= L1):
                    I.append(shi_s1(y, Vgs, Vds))
                elif (L1 < y <= (Ln1)):
                    I.append(shi_s2(y, Vgs, Vds))
                elif (Ln1 < y <= Ln2):
                    I.append(shi_s3(y, Vgs, Vds))
                elif (Ln2 < y <= (Ln2 + (L3(Vgs, Vds)))):
                    I.append(shi_s4(y, Vgs, Vds))
                elif (y > (Ln2 + (L3(Vgs, Vds)))):
                    I.append(shi_4(Vds))

            plt.plot(y2, (I))

        plt.xlabel('y --> ')
        plt.ylabel(' Ψs--2(y,Vgs,Vds) --> ')
        plotfile2 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile2)
        res2 = plotfile2
        return res2

    return shi_t1(),shi_t2()


#////////////////////

'''def a11(y,Vgs,Vds):
    v1 = eta_1/ts
    a1 = shi_s1(y,Vgs,Vds)  - shi_G3(Vgs)
    return v1*a1


def a12(y,Vgs,Vds):
    v1 = eta_2/ts
    a1 = shi_s2(y,Vgs,Vds)  - shi_G3(Vgs)
    return v1*a1

def a13(y,Vgs,Vds):
    v1 = eta_3/ts
    a1 = shi_s3(y,Vgs,Vds)  - shi_G3(Vgs)
    return v1*a1

def a14(y,Vgs,Vds):
    v1 = eta_4/ts
    a1 = shi_s4(y,Vgs,Vds)  - shi_G3(Vgs)
    return v1*a1

#///////
def a21(y,Vgs,Vds):

    v0,v01 = -eta_1 , ts**2
    v1 = v0/v01
    a1 = shi_s1(y,Vgs,Vds)  - shi_G3(Vgs)
    return v1*a1


def a22(y,Vgs,Vds):
    v0, v01 = -eta_2, ts ** 2
    v1 = v0 / v01
    a1 = shi_s2(y,Vgs,Vds)  - shi_G3(Vgs)
    return v1*a1

def a23(y,Vgs,Vds):
    v0, v01 = -eta_3, ts ** 2
    v1 = v0 / v01
    a1 = shi_s3(y,Vgs,Vds)  - shi_G3(Vgs)
    return v1*a1

def a24(y,Vgs,Vds):
    v0, v01 = -eta_4, ts ** 2
    v1 = v0 / v01
    a1 = shi_s4(y,Vgs,Vds)  - shi_G3(Vgs)
    return v1*a1


def Ex1(y,Vgs,Vds,x):
    return a11(y,Vgs,Vds) + (2*a21(y,Vgs,Vds)*x)

def Ex2(y,Vgs,Vds,x):
    return a12(y,Vgs,Vds) + (2*a22(y,Vgs,Vds)*x)

def Ex3(y,Vgs,Vds,x):
    return a13(y,Vgs,Vds) + (2*a23(y,Vgs,Vds)*x)
def Ex4(y,Vgs,Vds,x):
    return a14(y,Vgs,Vds) + (2*a24(y,Vgs,Vds)*x)


Vds = 1
vgs = [0,1.2]
x1 = np.arange(0,1.4,0.2)
y1 = np.arange(-10e-9,60e-9,5e-9)
#print(vgs , y1)
y2 = y1


for Vgs in vgs:
    for x in x1:
        Ex = []
        for y in y1:
            Vds = 1
            if (y <= 0):
                Ex.append(0)
            elif (0 < y <= L1):
                Ex.append(Ex1(y, Vgs, Vds, x))
            elif (L1 < y <= (Ln1)):
                Ex.append(Ex2(y, Vgs, Vds, x))
            elif (Ln1 < y <= Ln2):
                Ex.append(Ex3(y, Vgs, Vds, x))
            elif (Ln2 < y <= (Ln2 + L3(Vgs, Vds))):
                Ex.append(Ex4(y, Vgs, Vds, x))
            elif (y > (Ln2 + L3(Vgs, Vds))):
                Ex.append(0)

        plt.plot(y2, Ex)
plt.xlabel('y --> ')
plt.ylabel(' Ex(y,Vgs,Vds) --> ')
plt.show()

#//////////////////////////

def Ey1(y,Vgs,Vds):
    v1 , v2 = lmda_1*exp(lmda_1*y)  ,  lmda_1*exp(-lmda_1*y)
    a1 , a2 = A1(Vgs,Vds)*v1  , B1(Vgs,Vds)*v2
    return (a1 - a2)

def Ey2(y,Vgs,Vds):
    v1 , v2 = lmda_2*exp(lmda_2*y)  ,  lmda_2*exp(-lmda_2*y)
    a1 , a2 = A2(Vgs,Vds)*v1  , B2(Vgs,Vds)*v2
    return (a1 - a2)

def Ey3(y,Vgs,Vds):
    v1 , v2 = lmda_3*exp(lmda_3*y)  ,  lmda_3*exp(-lmda_3*y)
    a1 , a2 = A3(Vgs,Vds)*v1  , B3(Vgs,Vds)*v2
    return (a1 - a2)

def Ey4(y,Vgs,Vds):
    v1 , v2 = lmda_4*exp(lmda_4*y)  ,  lmda_4*exp(-lmda_4*y)
    a1 , a2 = A4(Vgs,Vds)*v1  , B4(Vgs,Vds)*v2
    return (a1 - a2)

Vds = 1
vgs = np.arange(0,1.2,0.2)
y1 = np.arange(-10e-9,60e-9,5e-9)
#print(vgs , y1)
y2 = y1

for Vgs in vgs:
    Ey = []
    for y in y1:
        Vds = 1
        if (y <= 0):
            Ey.append(0)
        elif (0 < y <= L1):
            Ey.append(Ey1(y,Vgs,Vds))
        elif (L1 < y <= (Ln1)):
            Ey.append(Ey2(y,Vgs,Vds))
        elif (Ln1 < y <= Ln2):
            Ey.append(Ey3(y,Vgs,Vds))
        elif (Ln2 < y <= (Ln2 + L3(Vgs, Vds))):
            Ey.append(Ey4(y,Vgs,Vds))
        elif (y > (Ln2 + L3(Vgs, Vds))):
            Ey.append(0)


    plt.plot(y2,Ey)
plt.xlabel('y --> ')
plt.ylabel(' Ey(y,Vgs,Vds) --> ')
plt.show()


//////////////////////////////////////////////////////////////////2D potential


////////////////////a01 not declared in math cad file

def shi_t_1(x,y,Vgs,Vds):
    v1 , v2 = x**2  , x
    return a01(y,Vgs,Vds) +a11(y,Vgs,Vds)*v2 + a21(y,Vgs,Vds)*v1
def shi_t_2(x,y,Vgs,Vds):
    v1 , v2 = x**2  , x
    return a02(y,Vgs,Vds) +a11(y,Vgs,Vds)*v2 + a22(y,Vgs,Vds)*v1
def shi_t_3(x,y,Vgs,Vds):
    v1, v2 = x ** 2, x
    return a03(y,Vgs,Vds) +a11(y,Vgs,Vds)*v2 + a23(y,Vgs,Vds)*v1
def shi_t_4(x,y,Vgs,Vds):
    v1, v2 = x ** 2, x
    return a04(y,Vgs,Vds) +a11(y,Vgs,Vds)*v2 + a24(y,Vgs,Vds)*v1

Vgs = 0.5
Vds = 1
x1 = np.arange(-1e-9,7e-9,1e-9)
y1 = np.arange(-20e-9,100e-9,5e-9)
#print(vgs , y1)
y2 = y1

for x in x1:
    I = []
    for y in y1:
        Vds = 1
        if (y <= 0):
            I.append(0)
        elif (0 < y <= L1):
            I.append(shi_t_1(x,y,Vgs,Vds))
        elif (L1 < y <= (Ln1)):
            I.append(shi_t_2(x,y,Vgs,Vds))
        elif (Ln1 < y <= Ln2):
            I.append(shi_t_3(x,y,Vgs,Vds))
        elif (Ln2 < y <= (Ln2 + L3(Vgs, Vds))):
            I.append(shi_t_4(x,y,Vgs,Vds))
        elif (y > (Ln2 + L3(Vgs, Vds))):
            I.append(shi_4(Vds))


    plt.plot(y2,I)
plt.xlabel('y --> ')
plt.ylabel(' Ψ(y,Vgs,Vds) --> ')
plt.show()

#//////////////////////  EFVB ///////////////

Vgs = 0.2
Vds = 1
x1 = np.arange(-1e-9,7e-9,1e-9)
y1 = np.arange(-20e-9,100e-9,5e-9)
#print(vgs , y1)
y2 = y1

J = Eg/2

for x in x1:
    I = []
    for y in y1:
        Vds = 1
        if (y <= 0):
            I.append(-shi_o - J)
        elif (0 < y <= L1):
            I.append(-shi_t_1(x,y,Vgs,Vds) - J )
        elif (L1 < y <= (Ln1)):
            I.append(-shi_t_2(x,y,Vgs,Vds) - J)
        elif (Ln1 < y <= Ln2):
            I.append(-shi_t_3(x,y,Vgs,Vds) - J)
        elif (Ln2 < y <= (Ln2 + L3(Vgs, Vds))):
            I.append(-shi_t_4(x,y,Vgs,Vds) - J)
        elif (y > (Ln2 + L3(Vgs, Vds))):
            I.append(-shi_4(Vds) - J)


    plt.plot(y2,I)
plt.xlabel('y --> ')
plt.ylabel(' EFVB(x, y,Vgs,Vds) --> ')
plt.show()


#/////////////////////////// EFCB ////////////

Vgs = 0.2
Vds = 1
x1 = np.arange(-1e-9,7e-9,1e-9)
y1 = np.arange(-20e-9,100e-9,5e-9)
#print(vgs , y1)
y2 = y1

J = Eg/2

for x in x1:
    I = []
    for y in y1:
        Vds = 1
        if (y <= 0):
            I.append(-shi_o + J)
        elif (0 < y <= L1):
            I.append(-shi_t_1(x,y,Vgs,Vds) + J )
        elif (L1 < y <= (Ln1)):
            I.append(-shi_t_2(x,y,Vgs,Vds) + J)
        elif (Ln1 < y <= Ln2):
            I.append(-shi_t_3(x,y,Vgs,Vds) + J)
        elif (Ln2 < y <= (Ln2 + L3(Vgs, Vds))):
            I.append(-shi_t_4(x,y,Vgs,Vds) + J)
        elif (y > (Ln2 + L3(Vgs, Vds))):
            I.append(-shi_4(Vds) + J)


    plt.plot(y2,I)
plt.xlabel('y --> ')
plt.ylabel(' EFCB(x,y,Vgs,Vds) --> ')
plt.show()


def E_total():
    a = np.array(Ex) ** 2
    b = np.array(Ey) ** 2
    return np.sqrt(a + b)'''






