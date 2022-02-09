import matplotlib.pyplot as plt
from numpy import *
from sympy import *
import numpy as np
import os,time


#//////////////values/////////////////////////
def as2(Vds):
    ni = 1.45e16
    tox = 2e-9
    N1, N2, N3, N4 = 1e26, ni, ni, 1e26
    # tox = 30e-9
    q, T, k, Eo, Eg = 1.6e-19, 300, 1.38e-23, 8.854e-12, 1.12
    W, L = 18e-9, 32e-9
    Es, Eox = 11.8 * Eo, 3.9 * Eo
    mu_n, Cox = 80 * 10 ** -4, Eo * Eox / tox
    kpa = 60e-5  # mu_n*Cox
    beta = kpa * (W / L)

    #Vgs, Vt, Vds = symbols('Vgs Vt Vds')

    def Idlin(Vgs, Vds, Vt):
        return beta * (Vgs - Vt - 0.5 * Vds) * Vds

    def Idsat(Vgs, Vt, Vds):
        return 0.5 * beta * ((Vgs - Vt) ** 2) * (1 + 0.01 * Vds)

    x = np.arange(0.5, 0.6, 1)
    vgs = np.arange(0, 1, 0.1)

    for x1 in x:
        I = []
        for Vgs in vgs:
            Vt = 0.2
            if Vgs <= Vt:
                I.append(0)
            elif Vgs > Vt and Vds <= (Vgs - Vt):
                I.append(Idlin(Vgs, Vds, Vt))
            elif Vgs > Vt and Vds > (Vgs - Vt):
                I.append(Idsat(Vgs, Vt, Vds))
        plt.plot(vgs, I ,  marker='o', markerfacecolor='skyblue', markersize=3, color='orange', linewidth=1, label = 'B_Mosfet')



    plt.xlabel('Gate Volyage(V) --> ')
    plt.ylabel(' Drain Current,Ids --> ')
    plt.legend()
    plotfile6 = os.path.join('static', str(time.time()) + '.png')
    plt.savefig(plotfile6)
    a6 = plotfile6
    return a6


'''n = 1.45e10
Nb = 1e15
Vsb = 0

Vto = 1
tox1 = 60e-9
Cox1 = Eox/tox1
Vfb = -1
VT = k*T/q
phi_f = 2*VT*log(Nb/n)
gamma = (sqrt(2*Es*q*Nb))/Cox1
#p = phi_si + Vsb
p1 = np.arange(0.5,4,0.5)
def Vth_V(p):
    return Vto + gamma*(sqrt(phi_f + Vsb) - sqrt(phi_f)'''



#plt.show()