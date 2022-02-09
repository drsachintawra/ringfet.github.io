import matplotlib.pyplot as plt               #////////////// import modules
import numpy as np
from sympy import *
from math import log, exp, sinh
from scipy import special
from scipy.integrate import trapz
import os,time
from mosfet import as2
from demo2 import as3
#from tcad import as0


def ash(amp1, amp2):                        #//////////////// input by user
    Vgs = amp1
    Vds = amp2

    r_range = np.arange(25e-9, 58e-9, 1e-9)        #//////////// constant data values
    #Z_range = np.arange(3e-9, 20e-9, 1e-9)
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
    Lch = 32e-9
    Rs = Rd + Lch
    Nd = 1e26
    x_s = 4.17
    E_g_eff = 1.12
    Na = 2e18
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
    data = []
    for i in r_range:
        data.append(Phi(i, 17 * 10 ** -9, Vgs, Vds))

        #/////////////////////////////////////////////////    surface _potential vs. position along channel {RingFET}
    def result1():
        plt.figure()
        plt.plot(data)
        plt.xlabel('Position along the Channel (nm)')
        plt.ylabel('Surface Potential (V)')
        plt.title('RING FET device')
        plotfile1 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile1)
        result = plotfile1
        return result

    #/////////////////////////////////////////////////////////////////    Drain current vs. Vgs {RingFET}
    wd = 13 * 10 ** -9
    z = np.arange(0, wd, 10 ** -9)

    phi_at_different_z = []
    int1_at_diffrent_z = []
    term1 = q / (K * T)
    Vgs = 1
    Vds = 0.05
    final_integration = []
    Vgs_range = np.arange(0, 1, 0.1)
    # Vds_range = np.arange(0,Vds,10**-1)
    # print(phi(25*10**-9,.5,.5, 0.05))
    data = []
    for i in r_range:
        data.append(Phi(i, 3 * 10 ** -9, 0, 0))

    # plt.plot(data)
    # plt.show()
    for vgs in Vgs_range:
        int1_at_diffrent_z = []
        for R in r_range:
            phi_at_different_z = []

            for z_val in z:
                phi_at_different_z.append(Phi(R, z_val, vgs, 0.05))

            final_term = np.exp(term1 * np.array(phi_at_different_z))

            int1 = trapz(final_term, z)

            int1_at_diffrent_z.append(int1)

        term1_for_int2 = 1 / (np.array(int1_at_diffrent_z))

        int2 = trapz(term1_for_int2, r_range)

        final_integration.append(int2)

    mu_n = 330 * 10 ** -4
    Ra = 2.302 * 10 ** -8
    Id_sub = 3.626 * 10 ** -13 / np.array(final_integration)  # mu_n*Ra*K*T*Ni(1-(np.exp((-q*Vds)/K*T)))

    # plt.plot(Vgs_range,Id_sub)
    # plt.show()

    ######################
    Ra1 = 5.341 * 10 ** -7
    mu_n1 = 80 * 10 ** -4
    E_c = (2 * 10 ** 4) / mu_n1
    Vt = 0.2

    def Idlin(Vgs, Vds):
        const1 = (Ra1 * Cox * mu_n1 * E_c) / (E_c * Lch + Vds)
        const2 = Vds * (Vgs - Vt) ** 1.45 - 0.4755 * Vds ** 2
        Lin = const1 * const2
        return Lin

    idlin = []
    for G_bias in Vgs_range:
        idlin.append(Idlin(G_bias, Vds))

    # plt.plot(Vgs_range,idlin)
    # plt.show()

    ############

    I_total1 = Id_sub[Vgs_range <= 0.2]
    I_total2 = np.array(idlin)[Vgs_range > 0.2]
    final_I = np.append(I_total1, I_total2)



    ooo = Vgs_range


    # print(final_I)
    def result2():
        plt.figure(figsize=(8,8))
        f1 = plt.plot(Vgs_range, final_I)
        plt.xlabel('Gate Bias (Vgs)')
        plt.ylabel('Drain Current Id (A)')
        plt.title('RING FET device')
        plotfile2 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile2)
        result = plotfile2
        return result


    #////////////////////////////////////////////////////////////////   Threshold voltage vs gate length
    def result3():
        Len = [32, 45, 60, 90]
        Vth1 = [0.2032, 0.3407, 0.4105, 0.4651]
        Vth2 = [0.1199, 0.2918, 0.3836, 0.451]
        plt.figure()
        plt.plot(Len, Vth1, 'r--', Len, Vth2, 'bs')
        plt.xlabel('Gate Length (Lch)')
        plt.ylabel('Threshold Voltage (Vth)')
        plt.title('RING FET device')
        plotfile3 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile3)
        result = plotfile3
        return result



    #///////////////////////////////////////////////////    transconductance vs. gate bias
    def result4():
        Vgs_range = np.arange(0.00, 0.92, 0.02)
        Vgs_range1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        gm1 = [2.18E-05, 3.56E-05, 5.81E-05, 9.50E-05, 1.55E-04, 2.53E-04, 4.12E-04, 6.69E-04, 1.08E-03, 1.74E-03,
               2.76E-03,
               4.34E-03, 6.71E-03, 1.01E-02, 1.49E-02, 2.11E-02, 2.86E-02, 3.70E-02, 4.58E-02, 5.41E-02, 6.15E-02,
               6.76E-02,
               7.23E-02, 7.59E-02, 7.84E-02, 8.00E-02, 8.11E-02, 8.17E-02, 8.19E-02, 8.19E-02, 8.17E-02, 8.13E-02,
               8.08E-02,
               8.03E-02, 7.96E-02, 7.90E-02, 7.83E-02, 7.75E-02, 7.68E-02, 7.60E-02, 7.52E-02, 7.44E-02, 7.36E-02,
               7.28E-02,
               7.19E-02, 7.11E-02]
        gm2 = [7.22E-05, 9.14E-04, 7.21E-03, 3.28E-02, 6.71E-02, 8.27E-02, 8.28E-02, 7.98E-02, 7.61E-02, 7.22E-02]
        plt.figure()
        plt.plot(Vgs_range, gm1, 'r--', Vgs_range1, gm2, 'bs')
        # plt.plot(Vgs_range, data)
        plt.xlabel('Gate Bias (Vgs)')
        plt.ylabel('Transconductance (gm)')
        plt.title('RING FET device')
        plotfile4 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile4)
        result = plotfile4
        return result


     #///////////////////////////////////////////////////  SS vs. gate length
    def result5():
        Len = [32, 45, 60, 90]
        SS1 = [93, 83.43, 79.09, 77.6]
        SS2 = [97.47, 83.34, 78.48, 77]
        plt.figure()
        plt.plot(Len, SS1, 'r--', Len, SS2, 'bs')
        plt.xlabel('Gate Length (Lch)')
        plt.ylabel('Sub-threshold Slope (SS)')
        plt.title('RING FET device')
        plotfile5 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile5)
        result5 = plotfile5
        return result5

    from DG_Id import DG_I            #import DG_I function from DG_Id(python file)
    def as10():
        ax = plt.plot(ooo, final_I , marker='', color='pink', linewidth=3 , label = 'R_FET')
        return ax

     #////////////////////////////////////////  drain current vs. Vgs (for comparison with )

    def result6():
        plt.figure(figsize= (8,8))
        #plt.subplot(221)
        as10()
        as2(amp2)
        # plt.subplot(222)
        DG_I()
        plt.xlabel('gate votlage, Vgs --> ')
        plt.ylabel(' Drain Current(Ids) --> ')
        plt.title('FET device')
        #plt.subplots_adjust(wspace=0.3, hspace=1, left=0.2)
        plotfile6 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile6)
        result = plotfile6
        return  result

    #////////////////////////////////////// ///////// Ids vs. Vgs (ringfet and bulk mosfet)
    def result61():
        plt.figure(figsize=(8, 8))
        # plt.subplot(221)
        as10()
        as2(amp2)
        plt.xlabel('gate votlage, Vgs --> ')
        plt.ylabel(' Drain Current(Ids) --> ')
        plt.title('FET device')
        # plt.subplots_adjust(wspace=0.3, hspace=1, left=0.2)
        plotfile61 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile61)
        result = plotfile61
        return result

        # ////////////////////////////////////////////// Ids vs. Vgs( for comparison with Ring and  DG_ MOSFET)

    def result62():
        plt.figure(figsize=(8, 8))
        # plt.subplot(221)
        as10()
        # plt.subplot(222)
        DG_I()
        plt.xlabel('gate votlage, Vgs --> ')
        plt.ylabel(' Drain Current(Ids) --> ')
        plt.title('FET device')
        #plt.subplots_adjust(wspace=0.3, hspace=1, left=0.2)
        plotfile62 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile62)
        result = plotfile62
        return result

    from DG_mosfet_SP import as4         #import DG_I function from DG_Id(python file)

    def as1():
        ax = plt.plot(data, marker='', color='red', linewidth=3 , label = 'R_FET')
        return ax


     #////////////////////////////////////  surface potential vs. channel length (for comp. with other devices)
    def result7():
        #plt.subplot(221)
        plt.figure(figsize=(8,8))
        as1()
        as3(amp1, amp2)
        plt.xlabel('Position along the Channel (nm)')
        plt.ylabel('Surface Potential (V)')
        plt.title('FET device')
        # plt.xlim([0,58e-9])
        # plt.subplot(223)
        # as0()

        # bottom=0.1, right=0.9, top=1.0 , left = 0.125 ,
        # plt.subplots_adjust( wspace = 1 , hspace =1)
        plotfile7 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile7)
        result = plotfile7
        return result

    #///////////////////////////////RingFET comp. with DG_FET /////////////////////////
    def result71():
        plt.figure(figsize=(8,8))
        as1()
        as4(amp1, amp2)
        plt.xlabel('Position along the Channel (nm)')
        plt.ylabel('Surface Potential (V)')
        plt.title('FET device')
        plotfile71 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile71)
        result = plotfile71
        return result

    # ///////////////////////////////RingFET comp. with DG_FET/Bulk mosfet /////////////////////////

    def result72():
        # plt.subplot(221)
        plt.figure(figsize=(8, 8))
        as1()
        as4(amp1,amp2)
        as3(amp1, amp2)
        plt.xlabel('Position along the Channel (nm)')
        plt.ylabel('Surface Potential (V)')
        plt.title('FET device')
        # plt.xlim([0,58e-9])
        # plt.subplot(223)
        # as0()

        # bottom=0.1, right=0.9, top=1.0 , left = 0.125 ,
        # plt.subplots_adjust( wspace = 1 , hspace =1)
        plotfile72 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile72)
        result = plotfile72
        return result

    #//////////////////////////////////// RingFET with different Lch /////////////////
    from SL_df import L_diff
    def result73():
        plt.figure(figsize=(8,8))
        L_diff(amp1,amp2)
        plt.xlabel('Position along the Channel (nm)')
        plt.ylabel('Surface Potential (V)')
        plt.title('FET device')
        plotfile8 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile8)
        result = plotfile8
        return result

    def result8():
        plt.figure(figsize=(8, 8))
        as4(amp1, amp2)
        plt.xlabel('Position along the Channel (nm)')
        plt.ylabel('Surface Potential (V)')
        plt.title('FET device')
        plotfile71 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile71)
        result = plotfile71
        return result


    #////////////////////////////////////////// comp. between bulk and DG_mosfet
    def result81():
        plt.figure(figsize=(8, 8))
        as3(amp1,amp2)
        as4(amp1, amp2)
        plt.xlabel('Position along the Channel (nm)')
        plt.ylabel('Surface Potential (V)')
        plt.title('FET device')
        plotfile71 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile71)
        result = plotfile71
        return result



    #////////////////////////////////////////// DG_MOSFET Ids vs. Vgs
    def result9():
        plt.figure(figsize=(8, 8))
        # plt.subplot(221)
        DG_I()
        plt.xlabel('gate votlage, Vgs --> ')
        plt.ylabel(' Drain Current(Ids) --> ')
        plt.title('FET device')
        # plt.subplots_adjust(wspace=0.3, hspace=1, left=0.2)
        plotfile62 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile62)
        result = plotfile62
        return result


    #///////////////////////////////////////// comp.between DG_ and RingFET  ......Ids vs. Vgs
    def result91():
        plt.figure(figsize=(8, 8))
        # plt.subplot(221)
        DG_I()
        as2(amp2)
        plt.xlabel('gate votlage, Vgs --> ')
        plt.ylabel(' Drain Current(Ids) --> ')
        plt.title('FET device')
        # plt.subplots_adjust(wspace=0.3, hspace=1, left=0.2)
        plotfile62 = os.path.join('static', str(time.time()) + '.png')
        plt.savefig(plotfile62)
        result = plotfile62
        return result




    return result1(),result2(),result3(),result4(),result5(),result6(),result7(),result73()\
        ,result71(),result72(),result61(),result62(),result8(),result81(),result9(),result91()







