# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 11:50:20 2022

@author: cqf52775
"""

import numpy as np
from scipy.optimize import fsolve
import math as ma

#input measured velocities (km/s) and density (10^3 kg/m^3)
v_1=5.130
v_2=1.900
v_3=2.520
v_4=2.460
v_5=4.460
v_6=1.940
v_7=2.510
v_8=2.490
v_9=4.050
v_10=2.270
v_11=2.280
v_12=4.530
v_13=4.560
v_14=2.460
v_15=2.460
v_16=4.550
v_17=1.840
v_18=2.210

rho=1.680

'''v_1=4.01
v_2=2.627
v_3=2.12
v_4=3.011
v_5=6.038
v_6=2.601
v_7=2.162
v_8=2.981
v_9=3.802
v_10=1.909
v_11=2.694
v_12=4.102
v_13=5.085
v_14=2.602
v_15=2.69
v_16=5.427
v_17=2.424
v_18=2.262
rho=1.722'''

#input errors for velocities (km/s) and density (10^3 kg/m^3)
e_1=0.05
e_2=0.05
e_3=0.05
e_4=0.05
e_5=0.05
e_6=0.05
e_7=0.05
e_8=0.05
e_9=0.05
e_10=0.05
e_11=0.05
e_12=0.05
e_13=0.05
e_14=0.05
e_15=0.05
e_16=0.05
e_17=0.05
e_18=0.05

e_rho=0.05

#first, do internal consistency checks on data (Dunk, Saunders 1984)
LHS1 = rho*(v_4**2+v_6**2)
RHS1 = rho*(v_2**2+v_8**2)
e_LHS1 = ma.sqrt((v_4**2+(v_6**2*e_rho))**2+(2*rho*v_4*e_4)**2+(2*rho*v_6*e_6)**2)
e_RHS1 = ma.sqrt((v_2**2+(v_8**2*e_rho))**2+(2*rho*v_2*e_2)**2+(2*rho*v_8*e_8)**2)

LHS2 = rho*(v_2**2 + 0.5*(v_1**2+v_3**2+v_5**2+v_8**2))
RHS2 = rho*(v_13**2+v_14**2+v_15**2)

LHS3 = rho*(v_8**2+0.5*(v_2**2+v_5**2+v_7**2+v_9**2))
RHS3 = rho*(v_16**2+v_17**2+v_18**2)

LHS4 = rho**2*v_2**2*v_8**2 - (rho*v_11**2 - 0.5*rho*(v_2**2+v_8**2))**2
RHS4 = rho**2*v_4**2*v_6**2

print('Internal consistency checks for inputted data:')
print('LHS1 =', LHS1, 'RHS1 =', RHS1)
print('LHS2 =', LHS2, 'RHS2 =', RHS2)
print('LHS3 =', LHS3, 'RHS3 =', RHS3)
print('LHS4 =', LHS4, 'RHS4 =', RHS4)

#calculate first 4 elastic constants directly from velocities
c_66=rho*(v_2**2)
c_22=rho*(v_5**2)
c_44=rho*(v_8**2)
c_46=rho*(v_11**2)-(0.5*(c_66+c_44))

#calculate errors for first 4 elastic constants
ec66=c_66*ma.sqrt((2*e_2/v_2)**2+(e_rho/rho)**2)
ec22=c_22*ma.sqrt((2*e_5/v_5)**2+(e_rho/rho)**2)
ec44=c_44*ma.sqrt((2*e_8/v_8)**2+(e_rho/rho)**2)
ec46=ma.sqrt((v_11**2*e_rho)**2+(2*rho*v_11*e_11)**2+(-0.5*ec66)**2+(-0.5*ec44)**2)

print('c_66 =',c_66,'+/-',ec66)
print('c_22 =',c_22, '+/-',ec22)
print('c_44 =',c_44, '+/-', ec44)
print('c_46 =',c_46, '+/-',ec46)

#use optimising/minimising method to find best guesses for next 5 elastic constants using initial estimates
def F(p):
    c_11 = p[0]
    c_15 = p[1]
    c_33 = p[2]
    c_55 = p[3]
    c_35 = p[4]

    f0=c_11+c_55-rho*((v_1**2)+(v_3**2))
    f1=c_11*c_55-c_15**2-(rho*v_1*v_3)**2
    f2=c_55+c_33-rho*((v_7**2)+(v_9**2))
    f3=c_33*c_55-c_35**2-(rho*v_7*v_9)**2
    f4=0.5*(c_11+2*c_55+c_33+2*(c_15+c_35))-rho*((v_10**2)+(v_12**2))

    return np.array([f0, f1, f2, f3, f4])

#input initial estimates
p0 = np.array([40,-1,23,14,-1])

p = fsolve(F, p0)
c_11 = p[0]
c_15 = p[1]
c_33 = p[2]
c_55 = p[3]
c_35 = p[4]

#loop to find max/min elastic constants for diff combinations of max/min input values
v1 = [v_1-e_1,v_1+e_1]
v3 = [v_3-e_3,v_3+e_3]
v7 = [v_7-e_7,v_7+e_7]
v9 = [v_9-e_9,v_9+e_9]
v10 = [v_10-e_10,v_10+e_10]
v12 = [v_12-e_12,v_12+e_12]
rho1 = [rho-e_rho,rho+e_rho]

ls_c11=[]
ls_c15=[]
ls_c33=[]
ls_c55=[]
ls_c35=[]
for i in range(2):
    v1a = v1[i]
    for j in range(2):
        v3a = v3[j]
        for k in range(2):
            v7a = v7[k]
            for l in range(2):
                v9a = v9[l]
                for m in range(2):
                    v10a = v10[m]
                    for n in range(2):
                        v12a = v12[n]
                        for o in range(2):
                            rhoa = rho1[o]
                            
                            #use optimising/minimising method to find best guesses for next 5 elastic constants using initial estimates
                            def f(p):
                                c11 = p[0]
                                c15 = p[1]
                                c33 = p[2]
                                c55 = p[3]
                                c35 = p[4]

                                f0=c11+c55-rhoa*((v1a**2)+(v3a**2))
                                f1=c11*c55-c15**2-(rhoa*v1a*v3a)**2
                                f2=c55+c33-rhoa*((v7a**2)+(v9a**2))
                                f3=c33*c55-c35**2-(rhoa*v7a*v9a)**2
                                f4=0.5*(c11+2*c55+c33+2*(c15+c35))-rhoa*((v10a**2)+(v12a**2))

                                return np.array([f0, f1, f2, f3, f4])

                            p0 = np.array([c_11,c_15,c_33,c_55,c_35])

                            p = fsolve(f, p0)
                            c11 = p[0]
                            c15 = p[1]
                            c33 = p[2]
                            c55 = p[3]
                            c35 = p[4]

                            ls_c11.append(c11)
                            ls_c15.append(c15)
                            ls_c33.append(c33)
                            ls_c55.append(c55)
                            ls_c35.append(c35)
                            

max_c11=max(ls_c11)
min_c11=min(ls_c11)
max_c15=max(ls_c15)
min_c15=min(ls_c15)
max_c33=max(ls_c33)
min_c33=min(ls_c33)
max_c55=max(ls_c55)
min_c55=min(ls_c55)
max_c35=max(ls_c35)
min_c35=min(ls_c35)

#half difference between max and min elastic constant values to find errors
ec11=(max_c11-min_c11)*0.5
ec15=(max_c15-min_c15)*0.5
ec33=(max_c33-min_c33)*0.5
ec55=(max_c55-min_c55)*0.5
ec35=(max_c35-min_c35)*0.5

print('c_11 =',c_11, '+/-',ec11)
print('c_15 =',c_15, '+/-', ec15)
print('c_33 =',c_33,'+/-',ec33)
print('c_55 =',c_55,'+/-',ec55)
print('c_35 =',c_35,'+/-',ec35)

#now use new elastic constants to find c_13 
def G(q):
    c_13 = q[0]
    f0 = (0.5*c_11 + 0.5*c_55 + c_15)*(0.5*c_55 + 0.5*c_33 + c_35) - (rho*v_10*v_12)**2 - (0.5*c_15 + 0.5*c_35 + 0.5*(c_13+c_55))**2
    return np.array([f0])

#input initial estimate for c_13
q0 = np.array([1])
q = fsolve(G, q0)
c_13 = q[0]

#loop to find max/min c_13 for diff combinations of max/min input values
c11c = [min_c11,max_c11]
c15c = [min_c15,max_c15]
c33c = [min_c33,max_c33]
c35c = [min_c35,max_c35]
c55c = [min_c55,max_c55]

ls_c13=[]
for i in range(2):
    c11a = c11c[i]
    for j in range(2):
        c15a = c15c[j]
        for k in range(2):
            c35a = c35c[k]
            for l in range(2):
                c55a = c55c[l]
                for m in range(2):
                    v10a = v10[m]
                    for n in range(2):
                        v12a = v12[n]
                        for o in range(2):
                            rhoa = rho1[o]
                            for p in range(2):
                                c33a = c33c[p]
                                def g(q):
                                    c13 = q[0]
                                    f0 = (0.5*c11a + 0.5*c55a + c15a)*(0.5*c55a + 0.5*c33a + c35a) - (rhoa*v10a*v12a)**2 - (0.5*c15a + 0.5*c35a + 0.5*(c13+c55a))**2
                                    return np.array([f0])

                                q0 = np.array([c_13])
                                q = fsolve(g, q0)
                                c13 = q[0]
                                
                                ls_c13.append(c13)

max_c13=max(ls_c13)
min_c13=min(ls_c13)
ec13=(max_c13-min_c13)*0.5

print('c_13 =',c_13, '+/-',ec13)

#now finding c_12 and c_25; start by checking that RHS and LHS of following eqn within 1% agreement
RHS = 2*rho*(v_13**2+v_14**2+v_15**2)
LHS = c_11 + c_22 + 2*c_66 + c_55 + c_44
pcnt_diff = 100*(RHS-LHS)/RHS
print('% difference for c_12,25 calculation = ', pcnt_diff)

#if <1% difference, continuing to find c_12 and c_25 with sim eqns
A_13 = (c_11+c_66-2*rho*v_13**2)*(c_22+c_66-2*rho*v_13**2)*(c_55+c_44-2*rho*v_13**2)
B_13 = c_11+c_66-2*rho*v_13**2
D_13 = c_55+c_44-2*rho*v_13**2
E_13 = c_15+c_46
F_13 = c_22+c_66-2*rho*v_13**2

A_14 = (c_11+c_66-2*rho*v_14**2)*(c_22+c_66-2*rho*v_14**2)*(c_55+c_44-2*rho*v_14**2)
B_14 = c_11+c_66-2*rho*v_14**2
D_14 = c_55+c_44-2*rho*v_14**2
E_14 = c_15+c_46
F_14 = c_22+c_66-2*rho*v_14**2

A_15 = (c_11+c_66-2*rho*v_15**2)*(c_22+c_66-2*rho*v_15**2)*(c_55+c_44-2*rho*v_15**2)
B_15 = c_11+c_66-2*rho*v_15**2
D_15 = c_55+c_44-2*rho*v_15**2
E_15 = c_15+c_46
F_15 = c_22+c_66-2*rho*v_15**2

def h1(r1):
    c_12_1=r1[0]
    c_25_1=r1[1]

    #f0 = A_13 - B_13*(c_25_1+c_46)**2 - D_13*(c_12_1+c_66)**2 + 2*E_13*(c_12_1+c_66)*(c_25_1+c_46)-(E_13**2)*F_13
    f1 = A_14 - B_14*(c_25_1+c_46)**2 - D_14*(c_12_1+c_66)**2 + 2*E_14*(c_12_1+c_66)*(c_25_1+c_46)-(E_14**2)*F_14
    f2 = A_15 - B_15*(c_25_1+c_46)**2 - D_15*(c_12_1+c_66)**2 + 2*E_15*(c_12_1+c_66)*(c_25_1+c_46)-(E_15**2)*F_15

    return np.array([f1, f2])

#input initial estimates
r0 = np.array([19,-2])

r1 = fsolve(h1, r0)
c_12_1=r1[0]
c_25_1=r1[1]

def h2(r2):
    c_12_2=r2[0]
    c_25_2=r2[1]

    f0 = A_13 - B_13*(c_25_2+c_46)**2 - D_13*(c_12_2+c_66)**2 + 2*E_13*(c_12_2+c_66)*(c_25_2+c_46)-(E_13**2)*F_13
    #f1 = A_14 - B_14*(c_25_2+c_46)**2 - D_14*(c_12_2+c_66)**2 + 2*E_14*(c_12_2+c_66)*(c_25_2+c_46)-(E_14**2)*F_14
    f2 = A_15 - B_15*(c_25_2+c_46)**2 - D_15*(c_12_2+c_66)**2 + 2*E_15*(c_12_2+c_66)*(c_25_2+c_46)-(E_15**2)*F_15

    return np.array([f0, f2])

r2 = fsolve(h2, r0)
c_12_2=r2[0]
c_25_2=r2[1]

def h3(r3):
    c_12_3=r3[0]
    c_25_3=r3[1]

    f0 = A_13 - B_13*(c_25_3+c_46)**2 - D_13*(c_12_3+c_66)**2 + 2*E_13*(c_12_3+c_66)*(c_25_3+c_46)-(E_13**2)*F_13
    f1 = A_14 - B_14*(c_25_3+c_46)**2 - D_14*(c_12_3+c_66)**2 + 2*E_14*(c_12_3+c_66)*(c_25_3+c_46)-(E_14**2)*F_14
    #f2 = A_15 - B_15*(c_25_3+c_46)**2 - D_15*(c_12_3+c_66)**2 + 2*E_15*(c_12_3+c_66)*(c_25_3+c_46)-(E_15**2)*F_15

    return np.array([f0, f1])


r3 = fsolve(h3, r0)
c_12_3=r3[0]
c_25_3=r3[1]

c_12=(c_12_1+c_12_2+c_12_3)/3
c_25=(c_25_1+c_25_2+c_25_3)/3

#define max and min arrays for velocities / elastic constants
v13 = [v_13-e_13,v_13+e_13]
v14 = [v_14-e_14,v_14+e_14]
v15 = [v_15-e_15,v_15+e_15]
c66m = [c_66-ec66,c_66+ec66]
c44m = [c_44-ec44,c_44+ec44]
c22m = [c_22-ec22,c_22+ec22]
c46m = [c_46-ec46,c_46+ec46]

ls_c12=[]
ls_c25=[]
for i in range(2):
    v13a = v13[i]
    for j in range(2):
        v14a = v14[j]
        for k in range(2):
            v15a = v15[k]
            for l in range(2):
                c55a = c55c[l]
                for m in range(2):
                    c66a = c66m[m]
                    for n in range(2):
                        c44a = c44m[n]
                        for o in range(2):
                            c22a = c22m[o]
                            for p in range(2):
                                c46a = c46m[p]
                                for q in range(2):
                                    c15a = c15c[q]
                                    for r in range(2):
                                        rhoa = rho1[r]
                                        for s in range(2):
                                            c11a = c11c[s]
                                            A_13a = (c11a+c66a-2*rhoa*v13a**2)*(c22a+c66a-2*rhoa*v13a**2)*(c55a+c44a-2*rhoa*v13a**2)
                                            B_13a = c11a+c66a-2*rhoa*v13a**2
                                            D_13a = c55a+c44a-2*rhoa*v13a**2
                                            E_13a = c15a+c46a
                                            F_13a = c22a+c66a-2*rhoa*v13a**2
                                            
                                            A_14a = (c11a+c66a-2*rhoa*v14a**2)*(c22a+c66a-2*rhoa*v14a**2)*(c55a+c44a-2*rhoa*v14a**2)
                                            B_14a = c11a+c66a-2*rhoa*v14a**2
                                            D_14a = c55a+c44a-2*rhoa*v14a**2
                                            E_14a = c15a+c46a
                                            F_14a = c22a+c66a-2*rhoa*v14a**2
                                            
                                            A_15a = (c11a+c66a-2*rhoa*v15a**2)*(c22a+c66a-2*rhoa*v15a**2)*(c55a+c44a-2*rhoa*v15a**2)
                                            B_15a = c11a+c66a-2*rhoa*v15a**2
                                            D_15a = c55a+c44a-2*rhoa*v15a**2
                                            E_15a = c15a+c46a
                                            F_15a = c22a+c66a-2*rhoa*v15a**2

                                            def h1(r1):
                                                c_12_1=r1[0]
                                                c_25_1=r1[1]

                                                #f0 = A_13a - B_13a*(c_25_1+c46a)**2 - D_13a*(c_12_1+c66a)**2 + 2*E_13a*(c_12_1+c66a)*(c_25_1+c46a)-(E_13a**2)*F_13a
                                                f1 = A_14a - B_14a*(c_25_1+c46a)**2 - D_14a*(c_12_1+c66a)**2 + 2*E_14a*(c_12_1+c66a)*(c_25_1+c46a)-(E_14a**2)*F_14a
                                                f2 = A_15a - B_15a*(c_25_1+c46a)**2 - D_15a*(c_12_1+c66a)**2 + 2*E_15a*(c_12_1+c66a)*(c_25_1+c46a)-(E_15a**2)*F_15a

                                                return np.array([f1, f2])

                                            #input initial estimates
                                            r0 = np.array([c_12,c_25])

                                            r1 = fsolve(h1, r0)
                                            c_12_1=r1[0]
                                            c_25_1=r1[1]

                                            def h2(r2):
                                                c_12_2=r2[0]
                                                c_25_2=r2[1]

                                                f0 = A_13a - B_13a*(c_25_2+c46a)**2 - D_13a*(c_12_2+c66a)**2 + 2*E_13a*(c_12_2+c66a)*(c_25_2+c46a)-(E_13a**2)*F_13a
                                                #f1 = A_14a - B_14a*(c_25_2+c46a)**2 - D_14a*(c_12_2+c66a)**2 + 2*E_14a*(c_12_2+c66a)*(c_25_2+c46a)-(E_14a**2)*F_14a
                                                f2 = A_15a - B_15a*(c_25_2+c46a)**2 - D_15a*(c_12_2+c66a)**2 + 2*E_15a*(c_12_2+c66a)*(c_25_2+c46a)-(E_15a**2)*F_15a

                                                return np.array([f0, f2])

                                            r2 = fsolve(h2, r0)
                                            c_12_2=r2[0]
                                            c_25_2=r2[1]

                                            def h3(r3):
                                                c_12_3=r3[0]
                                                c_25_3=r3[1]

                                                f0 = A_13a - B_13a*(c_25_3+c46a)**2 - D_13a*(c_12_3+c66a)**2 + 2*E_13a*(c_12_3+c66a)*(c_25_3+c46a)-(E_13a**2)*F_13a
                                                f1 = A_14a - B_14a*(c_25_3+c46a)**2 - D_14a*(c_12_3+c66a)**2 + 2*E_14a*(c_12_3+c66a)*(c_25_3+c46a)-(E_14a**2)*F_14a
                                                #f2 = A_15a - B_15a*(c_25_3+c46a)**2 - D_15a*(c_12_3+c66a)**2 + 2*E_15a*(c_12_3+c66a)*(c_25_3+c46a)-(E_15a**2)*F_15a
                                                
                                                return np.array([f0, f1])

                                            r3 = fsolve(h3, r0)
                                            c_12_3=r3[0]
                                            c_25_3=r3[1]

                                            c12=(c_12_1+c_12_2+c_12_3)/3
                                            c25=(c_25_1+c_25_2+c_25_3)/3
                                            
                                            ls_c12.append(c12)
                                            ls_c25.append(c25)

max_c12=max(ls_c12)
min_c12=min(ls_c12)
ec12=(max_c12-min_c12)*0.5

max_c25=max(ls_c25)
min_c25=min(ls_c25)
ec25=(max_c25-min_c25)*0.5 

print('c_12 =',c_12, '+/-',ec12)
print('c_25 =',c_25, '+/-',ec25)

#finally finding c_23, first check velocities
RHS_23 = 2*rho*(v_16**2+v_17**2+v_18**2)
LHS_23 = c_33 + c_22 + 2*c_44 + c_55 + c_66
pcnt_diff_23 = 100*(RHS_23-LHS_23)/RHS_23
print('% difference for c_23 calculation = ',pcnt_diff_23)

v=v_17
A_v = (c_33+c_44-2*rho*v**2)*(c_22+c_44-2*rho*v**2)*(c_55+c_66-2*rho*v**2)
B_v = (c_33+c_44-2*rho*v**2)*(c_25+c_46)**2
D_v = c_55+c_66-2*rho*v**2
E_v = (c_35+c_46)*(c_25+c_46)
F_v = (c_22+c_44-2*rho*v**2)*(c_35+c_46)**2

def j(s):
    c_23=s[0]

    f0 = A_v - B_v - D_v*(c_23+c_44)**2 + 2*E_v*(c_23+c_44) - F_v

    return np.array([f0])

#input initial estimates
s0 = np.array([20])

s = fsolve(j, s0)
c_23=s[0]

v16 = [v_16-e_16,v_16+e_16]
v17 = [v_17-e_17,v_17+e_17]
v18 = [v_18-e_18,v_18+e_18]
c25c = [c_25-ec25,c_25+ec25]

va=v17

ls_c23=[]
for i in range(2):
    v16a = v16[i]
    for j in range(2):
        v17a = v17[j]
        for k in range(2):
            v18a = v18[k]
            for l in range(2):
                c55a = c55c[l]
                for m in range(2):
                    c66a = c66m[m]
                    for n in range(2):
                        c44a = c44m[n]
                        for o in range(2):
                            c22a = c22m[o]
                            for p in range(2):
                                c46a = c46m[p]
                                for q in range(2):
                                    c35a = c35c[q]
                                    for r in range(2):
                                        rhoa = rho1[r]
                                        for s in range(2):
                                            c33a = c33c[s]
                                            for t in range(2):
                                                c25a = c25c[t]
                                                
                                                A_va = (c33a+c44a-2*rhoa*v16a**2)*(c22a+c44a-2*rhoa*v16a**2)*(c55a+c66a-2*rhoa*v16a**2)
                                                B_va = (c33a+c44a-2*rhoa*v16a**2)*(c25a+c46a)**2
                                                D_va = c55a+c66a-2*rhoa*v16a**2
                                                Ea = (c35a+c46a)*(c25a+c46a)
                                                F_va = (c22a+c44a-2*rhoa*v16a**2)*(c35a+c46a)**2
                                                
                                                s0 = np.array([c_23])
                                                def j1(s1):
                                                    c23=s1[0]
                                                    f0 = A_va - B_va - D_va*(c23+c44a)**2 + 2*Ea*(c23+c44a) - F_va
                                                    return np.array([f0])
                                                s1 = fsolve(j1, s0)
                                                c23=s1[0]

                                                ls_c23.append(c23)
                          

max_c23=max(ls_c23)
min_c23=min(ls_c23)
ec23=(max_c23-min_c23)*0.5

print('c_23 =',c_23, '+/-',ec23)

#Checking results using crystal stability conditions (found in Dunk 1987)
if c_11>0 and c_22>0 and c33>0 and c_44>0 and c_55>0 and c_66>0:
    print('Crystal stability check #1 successful (c_ii>0)')
else:
    print('Crystal stability check #1 unsuccessful (c_ii>0)')

if c_11*c_22 > c_12**2 and c_11*c_33 > c_13**2 and c_11*c_55 > c_15**2 and c_22*c_33>c_23**2 and c_22*c_55>c25**2 and c_33*c_55>c_35**2 and c_44*c_66>c_46**2:
    print('Crystal stability check #2 successful (c_ii*cjj>cij**2)')
else:
    print('Crystal stability check #2 unsuccessful (c_ii*cjj>cij**2)')

if c_11*c_22*c_33+2*c_12*c_13*c_23>c_11*c_23**2+c_33*c_12**2+c_22*c_12**2 and c_11*c_22*c_55+2*c_12*c_25*c_12>c_11*c_25**2+c_55*c_12**2+c_22*c_15**2 and c_22*c_33*c_55+2*c_23*c_35*c_25>c_22*c_35**2+c_55*c_23**2+c_33*c_35**2:
    print('Crystal stability check #3 successful (c_ii*c_jj*c_kk + 2*c_ij*c_ik*c_jk > c_ii*c_jk**2 + c_kk*c_ij**2 + c_jj*c_ik**2)')
else:
    print('Crystal stability check #3 unsuccessful (c_ii*c_jj*c_kk + 2*c_ij*c_ik*c_jk > c_ii*c_jk**2 + c_kk*c_ij**2 + c_jj*c_ik**2)')