# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:25:44 2022

@author: cqf52775
"""

#HEXAGONAL CALCS

print('HEXAGONAL ELASTIC CONSTANT CALCULATIONS')
print('- Using relations from Koizumi et al (2006)')
print('- Requires a crystal cut for 3 orientations ([001], [110], in (010) plane, at angle theta to [001] axis) each with up to 3 polarizations.')
print('- Tested with velocities from Koizumi (2006) (3.097,1.518,1.448,3.149,1.505,3.197,1.481); density=1.290; theta=22.8')

import numpy as np
import math as ma

#input measured velocities (km/s) and density (10^3 kg/m^3)
#velocity numeration is according to that in Koizumi, et al (2006)
#tested with velocities from Koizumi: 3.097,1.518,1.448,3.149,1.505,3.197,1.481; rho=1.290; angle=22.8
print(' ')
print('Start by entering measured velocities (km/s) and their errors for each orientation/polarization according to the following numeration:')
print('----------------------------------------------------------------------------------------------')
print('Velocity #\t\tPropagation direction\t\t\t\t\t\t\t\tPolarization direction')
print('----------------------------------------------------------------------------------------------')
print('v_1' + ':\t\t\t', '[110]', '\t\t\t\t\t\t\t\t\t\t\t\t', 'Longitudinal, [110]')
print('v_2' + ':\t\t\t', '[110]', '\t\t\t\t\t\t\t\t\t\t\t\t', 'Transverse, [001]')
print('v_3' + ':\t\t\t', '[110]', '\t\t\t\t\t\t\t\t\t\t\t\t', 'Transverse, [1-10]')
print('v_4' + ':\t\t\t', '[001]', '\t\t\t\t\t\t\t\t\t\t\t\t', 'Longitudinal, [001]')
print('v_5' + ':\t\t\t', 'In (010) plane, at angle theta to [001] axis', '\t\t', 'Transverse, [010]')
print('v_6' + ':\t\t\t', 'In (010) plane, at angle theta to [001] axis', '\t\t', 'Quasilongitudinal')
print('----------------------------------------------------------------------------------------------')

vlst=input('1. Input the 6 velocities in the order shown above (km/s), separated by commas. One of v_3 and v_5 not necessary for calculation of elastic constants, if one was not measured then enter 0: ').split(',')
elst=input('2. Input errors for the 6 velocities separated by commas: ').split(',')
rho=float(input('3. Enter density of crystal in 10^3 kg/m^3 or g/cm^3 = '))
e_rho=float(input('4. Enter error in density of crystal = '))
ang=float(input('5. Enter angle theta between propagation direction for v_5-v_7 and [001] axis = '))

th=(ang*np.pi)/180
v_1=float(vlst[0])
v_2=float(vlst[1])
v_3=float(vlst[2])
v_4=float(vlst[3])
v_5=float(vlst[4])
v_6=float(vlst[5])
e_1=float(elst[0])
e_2=float(elst[1])
e_3=float(elst[2])
e_4=float(elst[3])
e_5=float(elst[4])
e_6=float(elst[5])

c_44_1=rho*(v_2**2)
c_44_2=rho*(v_5**2)
c_44_a=[c_44_1,c_44_2]
n44=np.count_nonzero(c_44_a)
c_44_av=sum(c_44_a)/n44
c_44=c_44_av

c_33=rho*(v_4**2)
c_11=rho*(v_1**2-v_3**2)
c_12=rho*(v_1**2-3*v_3**2)

c_13=ma.sqrt(((2*rho*v_6**2-(c_11*(ma.sin(th))**2+c_33*(ma.cos(th))**2+c_44))**2-(c_11*(ma.sin(th))**2-c_33*(ma.cos(th))**2+c_44*(ma.cos(2*th)))**2)/(ma.sin(2*th))**2)-c_44

#loop to find max/min elastic constants for diff combinations of max/min input values
v1 = [v_1-e_1,v_1+e_1]
v2 = [v_2-e_2,v_2+e_2]
v3 = [v_3-e_3,v_3+e_3]
v4 = [v_4-e_4,v_4+e_4]
v5 = [v_5-e_5,v_5+e_5]
v6 = [v_6-e_6,v_6+e_6]
rho1 = [rho-e_rho,rho+e_rho]

ls_c11=[]
ls_c12=[]
ls_c44=[]
ls_c33=[]
ls_c13=[]
for i in range(2):
    v1a = v1[i]
    for j in range(2):
        v2a = v2[j]
        for k in range(2):
            v3a = v3[k]
            for l in range(2):
                v4a = v4[l]
                for m in range(2):
                    v5a = v5[m]
                    for n in range(2):
                        v6a = v6[n]
                        for p in range(2):
                            rhoa = rho1[p]
                            
                            c_44_1a=rhoa*(v2a**2)
                            c_44_2a=rhoa*(v5a**2)
                            c_44_aa=[c_44_1a,c_44_2a]
                            n44a=np.count_nonzero(c_44_aa)
                            c_44_ava=sum(c_44_aa)/n44a
                            c_44a=c_44_ava

                            c_33a=rhoa*(v4a**2)
                            c_11a=rhoa*(v1a**2-v3a**2)
                            c_12a=rhoa*(v1a**2-3*v3a**2)
                            
                            c_13a=ma.sqrt(((2*rhoa*v6a**2-(c_11a*(ma.sin(th))**2+c_33a*(ma.cos(th))**2+c_44a))**2-(c_11a*(ma.sin(th))**2-c_33a*(ma.cos(th))**2+c_44a*(ma.cos(2*th)))**2)/(ma.sin(2*th))**2)-c_44a
                            
                            ls_c11.append(c_11a)
                            ls_c12.append(c_12a)
                            ls_c44.append(c_44a)
                            ls_c33.append(c_33a)
                            ls_c13.append(c_13a)
                            

max_c11=max(ls_c11)
min_c11=min(ls_c11)
max_c12=max(ls_c12)
min_c12=min(ls_c12)
max_c44=max(ls_c44)
min_c44=min(ls_c44)
max_c33=max(ls_c33)
min_c33=min(ls_c33)
max_c13=max(ls_c13)
min_c13=min(ls_c13)

#half difference between max and min elastic constant values to find errors
ec11=(max_c11-min_c11)*0.5
ec12=(max_c12-min_c12)*0.5
ec44=(max_c44-min_c44)*0.5
ec33=(max_c33-min_c33)*0.5
ec13=(max_c13-min_c13)*0.5

print(' ')
print('--------------------------------------------')
print('Calculated elastic constants . . .')
print('--------------------------------------------')
print('c_11 =',c_11,'+/-',ec11)
print('c_33 =',c_33,'+/-',ec33)
print('c_44 =',c_44,'+/-',ec44)
print('c_12 =',c_12,'+/-',ec12)
print('c_13 =',c_13,'+/-',ec13)
print('--------------------------------------------')