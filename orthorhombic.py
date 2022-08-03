# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:18:52 2022

@author: cqf52775
"""

print('ORTHORHOMBIC ELASTIC CONSTANT CALCULATIONS')
print('- Using a simplified version of the method in K.S. Aleksandrov (1958) and velocity numeration from Dunk, Saunders (1984)')
print('- Requires a crystal cut for 6 orientations ([100],[010],[001],[101],[110],[011]) each with 3 polarizations.')
print('- NOTE: These calculations have not been checked against literature, so could be incorrect. Hard to check since there is no literature listing these exact velocities for orthorhombic - may be more useful to try a different method for orthorhombic with fewer velocities')

import numpy as np
from scipy.optimize import fsolve
import math as ma

#input measured velocities (km/s) and density (10^3 kg/m^3)
#velocity numeration is according to that in Dunk, Saunders (1984)
#tested with velocities from [] vlst=2.52,1,1.56,0.81,2.63,1.02,1.57,0.81,3.07,1.14,1.14,2.96,2.4,1.39,1.39,2.38,1.70,1.70; rho=2.98
print(' ')
print('Start by entering measured velocities (km/s) and their errors for each orientation/polarization according to the following numeration (Dunk, Saunders (1984)):')
print('-----------------------------------------------------------------')
print('Velocity #\t\tPropagation direction\t\tPolarization direction')
print('-----------------------------------------------------------------')
print('v_1' + ':\t\t\t', '[100]', '\t\t\t\t\t\t', '[100]')
print('v_2' + ':\t\t\t', '[100]', '\t\t\t\t\t\t', '[010]')
print('v_3' + ':\t\t\t', '[100]', '\t\t\t\t\t\t', '[001]')
print(' ')
print('v_4' + ':\t\t\t', '[010]', '\t\t\t\t\t\t', '[100]')
print('v_5' + ':\t\t\t', '[010]', '\t\t\t\t\t\t', '[010]')
print('v_6' + ':\t\t\t', '[010]', '\t\t\t\t\t\t', '[001]')
print(' ')
print('v_7' + ':\t\t\t', '[001]', '\t\t\t\t\t\t', '[100]')
print('v_8' + ':\t\t\t', '[001]', '\t\t\t\t\t\t', '[010]')
print('v_9' + ':\t\t\t', '[001]', '\t\t\t\t\t\t', '[001]')
print(' ')
print('v_10' + ':\t\t\t', '[101]', '\t\t\t\t\t\t', '[-101]')
print('v_11' + ':\t\t\t', '[101]', '\t\t\t\t\t\t', '[010]')
print('v_12' + ':\t\t\t', '[101]', '\t\t\t\t\t\t', '[101]')
print(' ')
print('v_13' + ':\t\t\t', '[110]', '\t\t\t\t\t\t', '[110]')
print('v_14' + ':\t\t\t', '[110]', '\t\t\t\t\t\t', '[001]')
print('v_15' + ':\t\t\t', '[110]', '\t\t\t\t\t\t', '[1-10]')
print(' ')
print('v_16' + ':\t\t\t', '[011]', '\t\t\t\t\t\t', '[011]')
print('v_17' + ':\t\t\t', '[011]', '\t\t\t\t\t\t', '[0-11]')
print('v_18' + ':\t\t\t', '[011]', '\t\t\t\t\t\t', '[100]')
print('-----------------------------------------------------------------')


print('NOTE: In the next part, the program will ask you to input these 18 velocities. For orthorhombic symmetry, v_10=v_11,v_14=v_15 and v_17=v_18 --> if one of each pair is not measured, just enter its equivalent velocity again.')

vlst=input('1. Input the 18 velocities in the order shown above (km/s), separated by commas : ').split(',')
elst=input('2. Input errors for the 18 velocities separated by commas: ').split(',')
rho=float(input('3. Enter density of crystal in 10^3 kg/m^3 or g/cm^3 = '))
e_rho=float(input('4. Enter error in density of crystal = '))

v_1=float(vlst[0])
v_2=float(vlst[1])
v_3=float(vlst[2])
v_4=float(vlst[3])
v_5=float(vlst[4])
v_6=float(vlst[5])
v_7=float(vlst[6])
v_8=float(vlst[7])
v_9=float(vlst[8])
v_10=float(vlst[9])
v_11=float(vlst[10])
v_12=float(vlst[11])
v_13=float(vlst[12])
v_14=float(vlst[13])
v_15=float(vlst[14])
v_16=float(vlst[15])
v_17=float(vlst[16])
v_18=float(vlst[17])

e_1=float(elst[0])
e_2=float(elst[1])
e_3=float(elst[2])
e_4=float(elst[3])
e_5=float(elst[4])
e_6=float(elst[5])
e_7=float(elst[6])
e_8=float(elst[7])
e_9=float(elst[8])
e_10=float(elst[9])
e_11=float(elst[10])
e_12=float(elst[11])
e_13=float(elst[12])
e_14=float(elst[13])
e_15=float(elst[14])
e_16=float(elst[15])
e_17=float(elst[16])
e_18=float(elst[17])

#calculate first 6 elastic constants directly from velocities
c_11=rho*(v_1**2)
c_66=(rho*(v_2**2)+rho*(v_6**2))/2
c_55=(rho*(v_3**2)+rho*(v_7**2))/2
c_44=(rho*(v_4**2)+rho*(v_8**2))/2
c_22=rho*(v_5**2)
c_33=rho*(v_9**2)

#calculate errors for first 6 elastic constants
ec11=c_11*ma.sqrt((2*e_1/v_1)**2+(e_rho/rho)**2)
ec66=c_66*ma.sqrt((2*e_2/v_2)**2+(e_rho/rho)**2)
ec55=c_55*ma.sqrt((2*e_3/v_3)**2+(e_rho/rho)**2)
ec44=c_44*ma.sqrt((2*e_4/v_4)**2+(e_rho/rho)**2)
ec22=c_22*ma.sqrt((2*e_5/v_5)**2+(e_rho/rho)**2)
ec33=c_33*ma.sqrt((2*e_9/v_9)**2+(e_rho/rho)**2)

print('c_11 =',c_11,'+/-',ec11)
print('c_22 =',c_22, '+/-',ec22)
print('c_33 =',c_33, '+/-', ec33)
print('c_44 =',c_44,'+/-',ec44)
print('c_55 =',c_55, '+/-',ec55)
print('c_66 =',c_66, '+/-', ec66)

min_c11=c_11-ec11
min_c22=c_22-ec22
min_c33=c_33-ec33
min_c44=c_44-ec44
min_c55=c_55-ec55
min_c66=c_66-ec66
max_c11=c_11-ec11
max_c22=c_22-ec22
max_c33=c_33-ec33
max_c44=c_44-ec44
max_c55=c_55-ec55
max_c66=c_66-ec66

c_131=ma.sqrt((c_11+c_55-rho*v_12**2)*(c_55+c_33-rho*v_12**2))-c_55

def G(q):
    c_13 = q[0]
    f0 = (0.5*c_11 + 0.5*c_55 )*(0.5*c_55 + 0.5*c_33) - (rho*v_10*v_12)**2 - (0.5*(c_13+c_55))**2
    return np.array([f0])

#input initial estimate for c_13
q0 = np.array([1])
q = fsolve(G, q0)
c_13 = q[0]

#loop to find max/min c_13 for diff combinations of max/min input values
c11c = [min_c11,max_c11]
c33c = [min_c33,max_c33]
c55c = [min_c55,max_c55]
v10 = [v_10-e_10,v_10+e_10]
v12 = [v_12-e_12,v_12+e_12]
rho1 = [rho-e_rho,rho+e_rho]

ls_c13=[]
for i in range(2):
    c11a = c11c[i]
    for j in range(2):
        c55a = c55c[j]
        for k in range(2):
            v10a = v10[k]
            for l in range(2):
                v12a = v12[l]
                for m in range(2):
                    rhoa = rho1[m]
                    for n in range(2):
                        c33a = c33c[n]
                        def g(q):
                            c13 = q[0]
                            f0 = (0.5*c11a + 0.5*c55a)*(0.5*c55a + 0.5*c33a) - (rhoa*v10a*v12a)**2 - (0.5*(c13+c55a))**2 
                            return np.array([f0])
                        q0 = np.array([c_13])
                        q = fsolve(g, q0)
                        c13 = q[0]
                        ls_c13.append(c13)

max_c13=max(ls_c13)
min_c13=min(ls_c13)
ec13=(max_c13-min_c13)*0.5

print('c_13 =',c_13, '+/-',ec13)

#now finding c_12, start by checking that RHS and LHS of following eqn within 1% agreement
RHS = 2*rho*(v_13**2+v_14**2+v_15**2)
LHS = c_11 + c_22 + 2*c_66 + c_55 + c_44
pcnt_diff = 100*(RHS-LHS)/RHS
print('% difference for c_12,25 calculation = ', pcnt_diff)

#if <1% difference, continuing to find c_12:
A_13 = (c_11+c_66-2*rho*v_13**2)*(c_22+c_66-2*rho*v_13**2)*(c_55+c_44-2*rho*v_13**2)
D_13 = c_55+c_44-2*rho*v_13**2

A_14 = (c_11+c_66-2*rho*v_14**2)*(c_22+c_66-2*rho*v_14**2)*(c_55+c_44-2*rho*v_14**2)
D_14 = c_55+c_44-2*rho*v_14**2

A_15 = (c_11+c_66-2*rho*v_15**2)*(c_22+c_66-2*rho*v_15**2)*(c_55+c_44-2*rho*v_15**2)
D_15 = c_55+c_44-2*rho*v_15**2

c12_13=ma.sqrt(A_13/D_13)-c_66
c12_14=ma.sqrt(A_14/D_14)-c_66
'''c12_15=ma.sqrt(A_15/D_15)-c_66'''
c_12=(c12_13+c12_14)/2

#define max and min arrays for velocities / elastic constants
v13 = [v_13-e_13,v_13+e_13]
v14 = [v_14-e_14,v_14+e_14]
v15 = [v_15-e_15,v_15+e_15]
c66m = [c_66-ec66,c_66+ec66]
c44m = [c_44-ec44,c_44+ec44]
c22m = [c_22-ec22,c_22+ec22]

ls_c12=[]
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
                            for r in range(2):
                                rhoa = rho1[r]
                                for s in range(2):
                                    c11a = c11c[s]
                                    
                                    A_13a = (c11a+c_66-2*rho*v_13**2)*(c_22+c_66-2*rho*v_13**2)*(c_55+c_44-2*rho*v_13**2)
                                    D_13a = c_55+c_44-2*rho*v_13**2

                                    A_14a = (c_11+c_66-2*rho*v_14**2)*(c_22+c_66-2*rho*v_14**2)*(c_55+c_44-2*rho*v_14**2)
                                    D_14a = c_55+c_44-2*rho*v_14**2

                                    A_15a = (c_11+c_66-2*rho*v_15**2)*(c_22+c_66-2*rho*v_15**2)*(c_55+c_44-2*rho*v_15**2)
                                    D_15a = c_55+c_44-2*rho*v_15**2

                                    c12_13a=ma.sqrt(A_13a/D_13a)-c66a
                                    c12_14a=ma.sqrt(A_14a/D_14a)-c66a
                                    '''c12_15a=ma.sqrt(A_15a/D_15a)-c66a'''
                                    c12=(c12_13a+c12_14a)/2
                                    
                                    ls_c12.append(c12)
                                              

max_c12=max(ls_c12)
min_c12=min(ls_c12)
ec12=(max_c12-min_c12)*0.5

print('c_12 =',c_12, '+/-',ec12)

#finally finding c_23, first check velocities
RHS_23 = 2*rho*(v_16**2+v_17**2+v_18**2)
LHS_23 = c_33 + c_22 + 2*c_44 + c_55 + c_66
pcnt_diff_23 = 100*(RHS_23-LHS_23)/RHS_23
print('% difference for c_23 calculation = ',pcnt_diff_23)

A_16 = (c_33+c_44-2*rho*v_16**2)*(c_22+c_44-2*rho*v_16**2)*(c_55+c_66-2*rho*v_16**2)
D_16 = c_55+c_66-2*rho*v_16**2

A_17 = (c_33+c_44-2*rho*v_17**2)*(c_22+c_44-2*rho*v_17**2)*(c_55+c_66-2*rho*v_17**2)
D_17 = c_55+c_66-2*rho*v_17**2

A_18 = (c_33+c_44-2*rho*v_18**2)*(c_22+c_44-2*rho*v_18**2)*(c_55+c_66-2*rho*v_18**2)
D_18 = c_55+c_66-2*rho*v_18**2

c23_16=ma.sqrt(A_16/D_16)-c_44
c23_17=ma.sqrt(A_17/D_17)-c_44
'''c23_18=ma.sqrt(A_18/D_18)-c_44'''
c_23=(c23_16+c23_17)/2

v16 = [v_16-e_16,v_16+e_16]
v17 = [v_17-e_17,v_17+e_17]
v18 = [v_18-e_18,v_18+e_18]

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
                            for r in range(2):
                                rhoa = rho1[r] 
                                for s in range(2):
                                    c33a = c33c[s]
                                    
                                    A_16a = (c33a+c44a-2*rhoa*v16a**2)*(c22a+c44a-2*rho*v16a**2)*(c55a+c66a-2*rhoa*v16a**2)
                                    D_16a = c55a+c66a-2*rhoa*v16a**2

                                    A_17a = (c33a+c44a-2*rhoa*v17a**2)*(c22a+c44a-2*rho*v17a**2)*(c55a+c66a-2*rhoa*v17a**2)
                                    D_17a = c55a+c66a-2*rhoa*v17a**2
                                    
                                    A_18a = (c33a+c44a-2*rhoa*v18a**2)*(c22a+c44a-2*rho*v18a**2)*(c55a+c66a-2*rhoa*v18a**2)
                                    D_18a = c55a+c66a-2*rhoa*v18a**2

                                    c23_16a=ma.sqrt(A_16a/D_16a)-c44a
                                    c23_17a=ma.sqrt(A_17a/D_17a)-c44a
                                    ''''c23_18a=ma.sqrt(A_18a/D_18a)-c44a'''
                                    c23=(c23_16a+c23_17a)/2
                                    
                                    ls_c23.append(c23)
                          

max_c23=max(ls_c23)
min_c23=min(ls_c23)
ec23=(max_c23-min_c23)*0.5

print('c_23 =',c_23, '+/-',ec23)

print('NOTE: These calculations have not been checked against literature, so could be incorrect. Hard to check since there is no literature listing these exact velocities for orthorhombic - may be more useful to try a different method for orthorhombic with fewer velocities')