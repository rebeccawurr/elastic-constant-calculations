# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 17:21:47 2022

@author: cqf52775
"""

print('TETRAGONAL (4, -4, 4/m) ELASTIC CONSTANT CALCULATIONS')
print('- Using relations from Chung, Li (1971)')
print('- Requires a crystal cut for 6 orientations ([001], [100], [1/root2 1/root2 0], [0 1/root2 1/root2], [root3/2 1/2 0], [1/2 root3/2 0]) each with 2/3 polarizations.')
print('- Tested with velocities from Chung, Li (revised) (4.77,2.77,5.10,2.77,2.96,5.55,2.77,2.41,5.00,2.61,3.08,2.78,5.12,3.01,2.77,5.52,2.17); density=4.54')
print('- Note that c_13 value calculation seems to have some error and needs to be checked. All other elastic constants should be accurate.')

import numpy as np
import math as ma
print(' ')
print('Start by entering measured velocities (km/s) and their errors for each orientation/polarization according to the following numeration (Chung, Li (1971)):')
print('-----------------------------------------------------------------')
print('Velocity #\t\tPropagation direction\t\tPolarization direction')
print('-----------------------------------------------------------------')
print('v_1' + ':\t\t\t', '[001]', '\t\t\t\t\t\t', '[001]')
print('v_2' + ':\t\t\t', '[001]', '\t\t\t\t\t\t', '[100]')
print(' ')
print('v_3' + ':\t\t\t', '[100]', '\t\t\t\t\t\t', '[100]')
print('v_4' + ':\t\t\t', '[100]', '\t\t\t\t\t\t', '[001]')
print('v_5' + ':\t\t\t', '[100]', '\t\t\t\t\t\t', '[010]')
print(' ')
print('v_6' + ':\t\t\t', '[1/root2 1/root2 0]', '\t\t', 'quasilongitudinal')
print('v_7' + ':\t\t\t', '[1/root2 1/root2 0]', '\t\t', 'shear [001]')
print('v_8' + ':\t\t\t', '[1/root2 1/root2 0]', '\t\t', 'quasishear')
print(' ')
print('v_9' + ':\t\t\t', '[0 1/root2 1/root2]', '\t\t', 'quasilongitudinal')
print('v_10' + ':\t\t\t', '[0 1/root2 1/root2]', '\t\t', 'quasishear')
print('v_11' + ':\t\t\t', '[0 1/root2 1/root2]', '\t\t', 'quasishear')
print(' ')
print('v_12' + ':\t\t\t', '[root3/2 1/2 0]', '\t\t\t', 'shear [001]')
print('v_13' + ':\t\t\t', '[root3/2 1/2 0]', '\t\t\t', 'quasilongitudinal')
print('v_14' + ':\t\t\t', '[root3/2 1/2 0]', '\t\t\t', 'quasishear')
print(' ')
print('v_15' + ':\t\t\t', '[1/2 root3/2 0]', '\t\t\t', 'shear [001]')
print('v_16' + ':\t\t\t', '[1/2 root3/2 0]', '\t\t\t', 'quasilongitudinal')
print('v_17' + ':\t\t\t', '[1/2 root3/2 0]', '\t\t\t', 'quasishear')
print('-----------------------------------------------------------------')

vlst=input('1. Input the 17 velocities in the order shown above (km/s), separated by commas. v_11,v_14,v_16 and v_17 not necessary for calculation of elastic constants, if one was not measured then enter 0: ').split(',')
elst=input('2. Input errors for the 17 velocities separated by commas: ').split(',')
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

#finding c33 and 3 values of c44 from 3 velocities
c_33=rho*v_1**2
c_44_1=rho*v_2**2
c_44_2=rho*v_4**2
c_44_3=rho*v_7**2
c_44_a=[c_44_1,c_44_2,c_44_3]
n44=np.count_nonzero(c_44_a)
c_44_av=sum(c_44_a)/n44
c_44=c_44_av

#calculating c11, c66, c16, c12
d=rho*((v_3**2)+(v_5**2))
g=(rho*v_3**2-rho*v_5**2)**2
h=(rho*v_6**2-rho*v_8**2)**2
f=(rho*v_13**2-rho*v_14**2)**2-(rho*v_16**2-rho*v_17**2)**2

t1=(h+g)+2*ma.sqrt(g*h-((f**2)/3))
t2=(h+g)-2*ma.sqrt(g*h-((f**2)/3))

ww=3**(-1/2)
w1a=ww*ma.sqrt(3*g-((f**2)/t1))
w1b=-ww*ma.sqrt(3*g-((f**2)/t1))
w2a=ww*ma.sqrt(3*g-((f**2)/t2))
w2b=-ww*ma.sqrt(3*g-((f**2)/t2))

c_11_1a=0.5*(d-w1a)
c_11_1b=0.5*(d-w1b)
c_11_2a=0.5*(d-w2a)
c_11_2b=0.5*(d-w2b)

c_66_1a=0.5*(d+w1a)
c_66_1b=0.5*(d+w1b)
c_66_2a=0.5*(d+w2a)
c_66_2b=0.5*(d+w2b)

c_16_1a1=0.5*ma.sqrt(g-w1a**2)
c_16_1a2=-0.5*ma.sqrt(g-w1a**2)
c_16_2a1=0.5*ma.sqrt(g-w2a**2)
c_16_2a2=-0.5*ma.sqrt(g-w2a**2)

c12v=2*ma.sqrt(3)

c11l=[c_11_1a,c_11_1b,c_11_2a,c_11_2b]
c66l=[c_66_1a,c_66_1b,c_66_2a,c_66_2b]
c16l=[c_16_1a1,c_16_1a2,c_16_2a1,c_16_2a2]
c12l=[(-i+(f/(c12v*j))) for i,j in zip(c11l,c16l)]

G=[((c_33+5*(c_44)**2+2*i**2+j**2+l**2+2*c_44*(l+c_33+j)-(4*((rho*v_9**2)**2+(rho*v_10**2)**2+(rho*v_11**2)**2)))/2) for i,j,l in zip(c16l,c66l,c11l)]
c13l1=[(-c_44+ma.sqrt(c_44**2-m)) for m in G]
c13l2=[(-c_44-ma.sqrt(c_44**2-m)) for m in G]

#loop to find max/min elastic constants for diff combinations of max/min input values
v1 = [v_1-e_1,v_1+e_1]
v2 = [v_2-e_2,v_2+e_2]
v3 = [v_3-e_3,v_3+e_3]
v4 = [v_4-e_4,v_4+e_4]
v5 = [v_5-e_5,v_5+e_5]
v6 = [v_6-e_6,v_6+e_6]
v7 = [v_7-e_7,v_7+e_7]
v8 = [v_8-e_8,v_8+e_8]
v9 = [v_9-e_9,v_9+e_9]
v10 = [v_10-e_10,v_10+e_10]
v11 = [v_11-e_11,v_11+e_11]
v12 = [v_12-e_12,v_12+e_12]
v13 = [v_13-e_13,v_13+e_13]
v14 = [v_14-e_14,v_14+e_14]
v15 = [v_15-e_15,v_15+e_15]
v16 = [v_16-e_16,v_16+e_16]
v17 = [v_17-e_17,v_17+e_17]
rho1 = [rho-e_rho,rho+e_rho]

ls_c44=[]
ls_c33=[]
ls_c111=[]
ls_c121=[]
ls_c1311=[]
ls_c1312=[]
ls_c161=[]
ls_c661=[]
ls_c112=[]
ls_c122=[]
ls_c1321=[]
ls_c1322=[]
ls_c162=[]
ls_c662=[]
ls_c113=[]
ls_c123=[]
ls_c1331=[]
ls_c1332=[]
ls_c163=[]
ls_c663=[]
ls_c114=[]
ls_c124=[]
ls_c1341=[]
ls_c1342=[]
ls_c164=[]
ls_c664=[]

for a1 in range(2):
    v1a = v1[a1]
    for a2 in range(2):
        v2a = v2[a2]
        for a3 in range(2):
            v3a = v3[a3]
            for a4 in range(2):
                v4a = v4[a4]
                for a5 in range(2):
                    v5a = v5[a5]
                    for a6 in range(2):
                        v6a = v6[a6]
                        for a7 in range(2):
                            v7a = v7[a7]
                            for a8 in range(2):
                                v8a = v8[a8]
                                for a9 in range(2):
                                    v9a = v9[a9]
                                    for a10 in range(2):
                                        v10a = v10[a10]
                                        for a11 in range(2):
                                            v11a = v11[a11]
                                            for a12 in range(2):
                                                v12a = v12[a12]
                                                for a13 in range(2):
                                                    v13a = v13[a13]
                                                    for a14 in range(2):
                                                        v14a = v14[a14]
                                                        for a15 in range(2):
                                                            v15a = v15[a15]
                                                            for a16 in range(2):
                                                                v16a = v16[a16]
                                                                for a17 in range(2):
                                                                    v17a = v17[a17]
                                                                    for a18 in range(2):
                                                                        rhoa = rho1[a18]
                                                                        
                                                                        #finding c33 and 3 values of c44 from 3 velocities
                                                                        c_33a=rhoa*v1a**2
                                                                        c_44_1a=rhoa*v2a**2
                                                                        c_44_2a=rhoa*v4a**2
                                                                        c_44_3a=rhoa*v7a**2
                                                                        c_44_aa=[c_44_1a,c_44_2a,c_44_3a]
                                                                        n44a=np.count_nonzero(c_44_aa)
                                                                        c_44_ava=sum(c_44_aa)/n44a
                                                                        c_44a=c_44_ava
                                                                        
                                                                        ls_c44.append(c_44a)
                                                                        ls_c33.append(c_33a)

                                                                        #calculating c11, c66, c16, c12
                                                                        da=rhoa*((v3a**2)+(v5a**2))
                                                                        ga=(rhoa*v3a**2-rhoa*v5a**2)**2
                                                                        ha=(rhoa*v6a**2-rhoa*v8a**2)**2
                                                                        fa=(rhoa*v13a**2-rhoa*v14a**2)**2-(rhoa*v16a**2-rhoa*v17a**2)**2
                                                                        
                                                                        t1a=(ha+ga)+2*ma.sqrt(ga*ha-((fa**2)/3))
                                                                        t2a=(ha+ga)-2*ma.sqrt(ga*ha-((fa**2)/3))

                                                                        w1aa=ww*ma.sqrt(3*ga-((fa**2)/t1a))
                                                                        w1ba=-ww*ma.sqrt(3*ga-((fa**2)/t1a))
                                                                        w2aa=ww*ma.sqrt(3*ga-((fa**2)/t2a))
                                                                        w2ba=-ww*ma.sqrt(3*ga-((fa**2)/t2a))
                                                                        
                                                                        c_111a=0.5*(da-w1aa)
                                                                        c_661a=0.5*(da+w1aa)
                                                                        c_161a=0.5*ma.sqrt(ga-w1aa**2)
                                                                        c_121a=-c_111a+(fa/(c12v*c_161a))
                                                                        Ga1=(c_33a+5*(c_44a)**2+2*c_161a**2+c_661a**2+c_111a**2+2*c_44a*(c_111a+c_33a+c_661a)-(4*((rhoa*v9a**2)**2+(rhoa*v10a**2)**2+(rhoa*v11a**2)**2)))/2
                                                                        
                                                                        c_1311a=-c_44a+ma.sqrt(c_44a**2-Ga1)
                                                                        c_1312a=-c_44a-ma.sqrt(c_44a**2-Ga1)
                                                                        
                                                                        c_112a=0.5*(da-w1ba)
                                                                        c_662a=0.5*(da+w1ba)
                                                                        c_162a=-0.5*ma.sqrt(ga-w1aa**2)
                                                                        c_122a=-c_112a+(fa/(c12v*c_162a))
                                                                        Ga2=(c_33a+5*(c_44a)**2+2*c_162a**2+c_662a**2+c_112a**2+2*c_44a*(c_112a+c_33a+c_662a)-(4*((rhoa*v9a**2)**2+(rhoa*v10a**2)**2+(rhoa*v11a**2)**2)))/2
                                                            
                                                                        c_1321a=-c_44a+ma.sqrt(c_44a**2-Ga2)
                                                                        c_1322a=-c_44a-ma.sqrt(c_44a**2-Ga2)
                                                                        
                                                                        c_113a=0.5*(da-w2aa)
                                                                        c_663a=0.5*(da+w2aa)
                                                                        c_163a=0.5*ma.sqrt(ga-w2aa**2)
                                                                        c_123a=-c_113a+(fa/(c12v*c_163a))
                                                                        Ga3=(c_33a+5*(c_44a)**2+2*c_163a**2+c_663a**2+c_113a**2+2*c_44a*(c_113a+c_33a+c_663a)-(4*((rhoa*v9a**2)**2+(rhoa*v10a**2)**2+(rhoa*v11a**2)**2)))/2
                                                                        
                                                                        c_1331a=-c_44a+ma.sqrt(c_44a**2-Ga3)
                                                                        c_1332a=-c_44a-ma.sqrt(c_44a**2-Ga3)
                                                                        
                                                                        c_114a=0.5*(da-w2ba)
                                                                        c_664a=0.5*(da+w2ba)
                                                                        c_164a=-0.5*ma.sqrt(ga-w2aa**2)
                                                                        c_124a=-c_114a+(fa/(c12v*c_164a))
                                                                        
                                                                        Ga4=(c_33a+5*(c_44a)**2+2*c_164a**2+c_664a**2+c_114a**2+2*c_44a*(c_114a+c_33a+c_664a)-(4*((rhoa*v9a**2)**2+(rhoa*v10a**2)**2+(rhoa*v11a**2)**2)))/2
                                                                        
                                                                        c_1341a=-c_44a+ma.sqrt(c_44a**2-Ga4)
                                                                        c_1342a=-c_44a-ma.sqrt(c_44a**2-Ga4)
                                                                        
                                                                        ls_c111.append(c_111a)
                                                                        ls_c121.append(c_121a)
                                                                        ls_c1311.append(c_1311a)
                                                                        ls_c1312.append(c_1312a)
                                                                        ls_c161.append(c_161a)
                                                                        ls_c661.append(c_661a)
                                                                        
                                                                        ls_c112.append(c_112a)
                                                                        ls_c122.append(c_122a)
                                                                        ls_c1321.append(c_1321a)
                                                                        ls_c1322.append(c_1322a)
                                                                        ls_c162.append(c_162a)
                                                                        ls_c662.append(c_662a)
                                                                        
                                                                        ls_c113.append(c_113a)
                                                                        ls_c123.append(c_123a)
                                                                        ls_c1331.append(c_1331a)
                                                                        ls_c1332.append(c_1332a)
                                                                        ls_c163.append(c_163a)
                                                                        ls_c663.append(c_663a)
                                                                        
                                                                        ls_c114.append(c_114a)
                                                                        ls_c124.append(c_124a)
                                                                        ls_c1341.append(c_1341a)
                                                                        ls_c1342.append(c_1342a)
                                                                        ls_c164.append(c_164a)
                                                                        ls_c664.append(c_664a)
                                                                      
max_c44=max(ls_c44)
min_c44=min(ls_c44)
max_c33=max(ls_c33)
min_c33=min(ls_c33)
max_c11=[max(ls_c111),max(ls_c112),max(ls_c113),max(ls_c114)]
min_c11=[min(ls_c111),min(ls_c112),min(ls_c113),min(ls_c114)]
max_c66=[max(ls_c661),max(ls_c662),max(ls_c663),max(ls_c664)]
min_c66=[min(ls_c661),min(ls_c662),min(ls_c663),min(ls_c664)]
max_c16=[max(ls_c161),max(ls_c162),max(ls_c163),max(ls_c164)]
min_c16=[min(ls_c161),min(ls_c162),min(ls_c163),min(ls_c164)]
max_c12=[max(ls_c121),max(ls_c122),max(ls_c123),max(ls_c124)]
min_c12=[min(ls_c121),min(ls_c122),min(ls_c123),min(ls_c124)]
max_c131=[max(ls_c1311),max(ls_c1321),max(ls_c1331),max(ls_c1341)]
min_c131=[min(ls_c1311),min(ls_c1321),min(ls_c1331),min(ls_c1341)]
max_c132=[max(ls_c1312),max(ls_c1322),max(ls_c1332),max(ls_c1342)]
min_c132=[min(ls_c1312),min(ls_c1322),min(ls_c1332),min(ls_c1342)]

print('. . .')

for i in range(4):
    c11a=c11l[i]
    c66a=c66l[i]
    c16a=c16l[i]
    c12a=c12l[i]
    c13a1=c13l1[i]
    c13a2=c13l2[i]
    c13a12=[c13a1,c13a2]
    print(' ')
    print('Performing stability checks on results set #',i+1,'...')
    print('Set #',i+1,' : c_11 = ',c11a,', c_66 = ',c66a,', c_16 = ',c16a,', c_12 = ',c12a,'c_13 value #1 =', c13a1,'c_13 value #2 =', c13a2,'.')
    if c11a>0 and c66a>0 and c16a**2<c11a*c66a and c11a>abs(c12a):
        print('Checks #1-4 successful')
    else: 
        print('Some checks unsuccessful:')
        if c11a>0: print('Check #1 (c_11>0) successful')
        else: print('Check #1 (c_11>0) unsuccessful')
        if c66a>0: print('Check #2 (c_66>0) successful')
        else: print('Check #2 (c_66>0) unsuccessful')
        if c16a**2 < c11a*c66a: print('Check #3 (c_16**2 < c_11*c_66) successful')
        else: print('Check #3 (c_16**2 < c_11*c_66) unsuccessful')
        if c11a > abs(c12a): print('Check #4 (c_11 > |c_12|) successful')
        else: print('Check #4 (c_11 > |c_12|) unsuccesful')
    for l in range(2):
        c13a=c13a12[l]
        if c13a**2<c11a*c_33:
            print('Check #5 successful for c_13 value #',l+1)
        else: print('Check #5 unsuccessful for c_13 value #',l+1)
            
print('. . .')

j=int(input('Type the number of the set you want to use (1/2/3/4):'))-1
print('. . .')
print('For set #',j+1,'c_13 value #1=',c13l1[j],', c_13 value #2=',c13l2[j])
n=int(input('Type the number of the c_13 value you want to use (1/2):'))-1
print('. . .')


#half difference between max and min elastic constant values to find errors
ec11=(max_c11[j]-min_c11[j])*0.5
ec12=(max_c12[j]-min_c12[j])*0.5
ec44=(max_c44-min_c44)*0.5
ec33=(max_c33-min_c33)*0.5
ec131=(max_c131[j]-min_c131[j])*0.5
ec132=(max_c132[j]-min_c132[j])*0.5
ec66=(max_c66[j]-min_c66[j])*0.5
ec16=(max_c16[j]-min_c16[j])*0.5

c13a1=c13l1[j]
c13a2=c13l2[j]
c13a12=[c13a1,c13a2]
ec13l=[ec131,ec132]

print(' ')
print('----------------------------------------------------')
print('Calculated elastic constants . . .')
print('----------------------------------------------------')
print('c_11 =',c11l[j],'+/-',ec11)
print('c_12 =',c12l[j],'+/-',ec12)
print('c_13* =',c13a12[n],'+/-',ec13l[n])
print('c_16 =',c16l[j],'+/-',ec16)
print('c_33 =',c_33,'+/-',ec33)
print('c_44 =',c_44,'+/-',ec44)
print('c_66 =',c66l[j],'+/-',ec66)
print('----------------------------------------------------')

print('* Note that c_13 value calculation seems to have some error and needs to be checked. All other elastic constants should be accurate.')