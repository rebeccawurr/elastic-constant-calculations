# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 12:00:17 2022

@author: cqf52775
"""
#CUBIC CALCS (http://www.ioffe.ru/SVA/NSM/Semicond/Ge/mechanic.html#Acoustic)

print('CUBIC ELASTIC CONSTANT CALCULATIONS')
print('- Using relations from R.G. Leisure - Ultrasonic Spectroscopy (2017)')
print('- Requires a crystal cut for 2 orientations ([100] and [110]) each with 2/3 polarizations.')
print('- Tested with Ge velocities from http://www.ioffe.ru/SVA/NSM/Semicond/Ge/mechanic.html#Acoustic (4.87,3.57,5.36,3.57,2.77) and density=5.323')

from scipy.optimize import fsolve
import numpy as np

#tested with velocities from ioffe: vlst=4.87,3.57,5.36,3.57,2.77; rho=5.323
print(' ')
print('Start by entering measured velocities (km/s) and their errors for each orientation/polarization according to the following numeration:')
print('------------------------------------------------------------------------')
print('Velocity #\t\tPropagation direction\t\tPolarization direction')
print('------------------------------------------------------------------------')
print('v_1' + ':\t\t\t', '[100]', '\t\t\t\t\t\t', 'Longitudinal, [100]')
print('v_2' + ':\t\t\t', '[100]', '\t\t\t\t\t\t', 'Transverse [010]/[001]')
print('v_3' + ':\t\t\t', '[110]', '\t\t\t\t\t\t', 'Longitudinal, [110]')
print('v_4' + ':\t\t\t', '[110]', '\t\t\t\t\t\t', 'Transverse [001]')
print('v_5' + ':\t\t\t', '[110]', '\t\t\t\t\t\t', 'Transverse [1-10]')
print('------------------------------------------------------------------------')

vlst=input('1. Input the 5 velocities in the order shown above (km/s), separated by commas. v_4 and v_5 not necessary for calculation of elastic constants, if any not measured then enter 0: ').split(',')
elst=input('2. Input errors for the 5 velocities separated by commas = ').split(',')
rho=float(input('3. Enter density of crystal in 10^3 kg/m^3 or g/cm^3 = '))
e_rho=float(input('4. Enter error in density of crystal = '))

v_1=float(vlst[0])
v_2=float(vlst[1])
v_3=float(vlst[2])
v_4=float(vlst[3])
v_5=float(vlst[4])
e_1=float(elst[0])
e_2=float(elst[1])
e_3=float(elst[2])
e_4=float(elst[3])
e_5=float(elst[4])

c_11_g=rho*(v_1**2)

c_44_1=rho*(v_2**2)
c_44_2=rho*(v_4**2)
c_44_a=[c_44_1,c_44_2]
n44=np.count_nonzero(c_44_a)
c_44_av=sum(c_44_a)/n44

c_12_1=2*rho*v_3**2-c_11_g-2*c_44_av
c_12_2=c_11_g-2*rho*v_5**2
c_12_a=[c_12_1,c_12_2]
n12=np.count_nonzero(c_12_a)
c_12_av=sum(c_12_a)/n12

#use optimising/minimising method to find best guesses for elastic constants using initial estimates
def F(p):
    c_11 = p[0]
    c_12 = p[1]
    c_44 = p[2]

    f0=c_11-rho*(v_1**2)+c_11-c_12-2*rho*v_5**2
    f1=c_44-rho*(v_2**2)+c_44-rho*(v_4**2)
    f2=c_11+c_12+2*c_44-2*rho*(v_3**2)

    return np.array([f0, f1, f2])

#input initial estimates
p0 = np.array([c_11_g,c_12_av,c_44_av])

p = fsolve(F, p0)
c_11 = p[0]
c_12 = p[1]
c_44 = p[2]

#loop to find max/min elastic constants for diff combinations of max/min input values
v1 = [v_1-e_1,v_1+e_1]
v2 = [v_2-e_2,v_2+e_2]
v3 = [v_3-e_3,v_3+e_3]
v4 = [v_4-e_4,v_4+e_4]
v5 = [v_5-e_5,v_5+e_5]
rho1 = [rho-e_rho,rho+e_rho]

ls_c11=[]
ls_c12=[]
ls_c44=[]
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
                    for o in range(2):
                        rhoa = rho1[o]
                        
                        #use optimising/minimising method to find best guesses for elastic constants using initial estimates
                        def f(p):
                            c11 = p[0]
                            c12 = p[1]
                            c44 = p[2]

                            f0=c11-rhoa*(v1a**2)+c11-c12-2*rhoa*v5a**2
                            f1=c44-rhoa*(v2a**2)+c44-rhoa*(v4a**2)
                            f2=c11+c12+2*c44-2*rhoa*(v3a**2)
                            
                            return np.array([f0, f1, f2])

                        #input initial estimates
                        p0 = np.array([c_11,c_12,c_44])

                        p = fsolve(f, p0)
                        c11 = p[0]
                        c12 = p[1]
                        c44 = p[2]
                        
                        ls_c11.append(c11)
                        ls_c12.append(c12)
                        ls_c44.append(c44)
                            

max_c11=max(ls_c11)
min_c11=min(ls_c11)
max_c12=max(ls_c12)
min_c12=min(ls_c12)
max_c44=max(ls_c44)
min_c44=min(ls_c44)

#half difference between max and min elastic constant values to find errors
ec11=(max_c11-min_c11)*0.5
ec12=(max_c12-min_c12)*0.5
ec44=(max_c44-min_c44)*0.5

print(' ')
print('--------------------------------------------')
print('Calculated elastic constants . . .')
print('--------------------------------------------')
print('c_11 =',c_11,'+/-',ec11)
print('c_12 =',c_12,'+/-',ec12)
print('c_44 =',c_44,'+/-',ec44)
print('--------------------------------------------')
