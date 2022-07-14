# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 13:24:18 2022

@author: cqf52775
"""
import matplotlib.pyplot as plt 
import numpy as np
from scipy.optimize import fsolve, minimize
import math as ma
import uncertainties as uc

#input measured velocities (km/s) and density (10^3 kg/m^3)
v_1=4.01
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
rho=1.722

'''v_1=5.130
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
rho=1.680'''

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

pk1=1/(e_1)*(1/2)
pk2=1/(e_2)*(1/2)
pk3=1/(e_3)*(1/2)
pk4=1/(e_4)*(1/2)
pk5=1/(e_5)*(1/2)
pk6=1/(e_6)*(1/2)
pk7=1/(e_7)*(1/2)
pk8=1/(e_8)*(1/2)
pk9=1/(e_9)*(1/2)
pk10=1/(e_10)*(1/2)
pk11=1/(e_11)*(1/2)
pk12=1/(e_12)*(1/2)
pk13=1/(e_13)*(1/2)
pk14=1/(e_14)*(1/2)
pk15=1/(e_15)*(1/2)
pk16=1/(e_16)*(1/2)
pk17=1/(e_17)*(1/2)
pk18=1/(e_18)*(1/2)


ve_1=uc.ufloat(v_1,e_1)
ve_2=uc.ufloat(v_2,e_2)
ve_3=uc.ufloat(v_3,e_3)
ve_4=uc.ufloat(v_4,e_4)
ve_5=uc.ufloat(v_5,e_5)
ve_6=uc.ufloat(v_6,e_6)
ve_7=uc.ufloat(v_7,e_7)
ve_8=uc.ufloat(v_8,e_8)
ve_9=uc.ufloat(v_9,e_9)
ve_10=uc.ufloat(v_10,e_10)
ve_11=uc.ufloat(v_11,e_11)
ve_12=uc.ufloat(v_12,e_12)
ve_13=uc.ufloat(v_13,e_13)
ve_14=uc.ufloat(v_14,e_14)
ve_15=uc.ufloat(v_15,e_15)
ve_16=uc.ufloat(v_16,e_16)
ve_17=uc.ufloat(v_17,e_17)
ve_18=uc.ufloat(v_18,e_18)
rhoe=uc.ufloat(rho,e_rho)


#first, do internal consistency checks on data (Dunk, Saunders 1984)
LHS1 = rhoe*(ve_4**2+ve_6**2)
RHS1 = rhoe*(ve_2**2+ve_8**2)

LHS2 = rhoe*(ve_2**2 + 0.5*(ve_1**2+ve_3**2+ve_5**2+ve_8**2))
RHS2 = rhoe*(ve_13**2+ve_14**2+ve_15**2)

LHS3 = rhoe*(ve_8**2+0.5*(ve_2**2+ve_5**2+ve_7**2+ve_9**2))
RHS3 = rho*(ve_16**2+ve_17**2+ve_18**2)

LHS4 = rhoe**2*ve_2**2*ve_8**2 - (rhoe*ve_11**2 - 0.5*rhoe*(ve_2**2+ve_8**2))**2
RHS4 = rhoe**2*ve_4**2*ve_6**2

print('Internal consistency checks for inputted data:')
print('LHS1 =', LHS1, 'RHS1 =', RHS1)
print('LHS2 =', LHS2, 'RHS2 =', RHS2)
print('LHS3 =', LHS3, 'RHS3 =', RHS3)
print('LHS4 =', LHS4, 'RHS4 =', RHS4)

cx=[1,2,3,4,5,6,7,8]
cy=[LHS1.n,RHS1.n,LHS2.n,RHS2.n,LHS3.n,RHS3.n,LHS4.n,RHS4.n]
cye=[LHS1.s,RHS1.s,LHS2.s,RHS2.s,LHS3.s,RHS3.s,LHS4.s,RHS4.s]

plt.errorbar(cx, cy, yerr = cye,fmt='o',ecolor = 'red',color='black')
plt.title("Internal Consistency Checks")
plt.show()
print("Ensure pairs are within error of each other before continuing")
input("Press enter to continue...")


#finding initial estimates for elastic constants
#calculate first 4 elastic constants directly from velocities
c660=rho*(v_2**2)
c220=rho*(v_5**2)
c440=rho*(v_8**2)
c460=rho*(v_11**2)-(0.5*(c660+c440))

#calculate errors for first 4 elastic constants
ec660=c660*ma.sqrt((2*e_2/v_2)**2+(e_rho/rho)**2)
ec220=c220*ma.sqrt((2*e_5/v_5)**2+(e_rho/rho)**2)
ec440=c440*ma.sqrt((2*e_8/v_8)**2+(e_rho/rho)**2)
ec460=ma.sqrt((v_11**2*e_rho)**2+(2*rho*v_11*e_11)**2+(-0.5*ec660)**2+(-0.5*ec440)**2)


f0=rho*((v_1**2)+(v_3**2))
f1=(rho*v_1*v_3)**2
f2=rho*((v_7**2)+(v_9**2))
f3=(rho*v_7*v_9)**2
f4=rho*((v_10**2)+(v_12**2)-0.5*(v_1**2+v_3**2+v_7**2+v_9**2))

def fxn1(x):
    return(-(-x**2+f0*x-f1)**(1/2)+(-x**2+x*f2-f3)**(1/2)-f4)
def fxn2(x):
    return(-(-x**2+f0*x-f1)**(1/2)-(-x**2+x*f2-f3)**(1/2)-f4)
def fxn3(x):
    return((-x**2+f0*x-f1)**(1/2)-(-x**2+x*f2-f3)**(1/2)-f4)
def fxn4(x):
    return((-x**2+f0*x-f1)**(1/2)+(-x**2+x*f2-f3)**(1/2)-f4)

xlist=np.linspace(8,15,num=1000)

plt.figure(0)
plt.plot(xlist,fxn1(xlist),label="Fxn1")
plt.plot(xlist,fxn2(xlist),label="Fxn2")
plt.plot(xlist,fxn3(xlist),label="Fxn3")
plt.plot(xlist,fxn4(xlist),label="Fxn4")
plt.plot(xlist,xlist*0,"--",label="y=0")
plt.ylim([-1,1])
plt.legend()
plt.title("Roots of function to find c55")
plt.show()

print("Enter initial x values close to roots for each function (enter 0 if function doesn't have root):")
guess_1a=int(input("Fxn1 root guess 1 = "))
guess_1b=int(input("Fxn1 root guess 2 = "))
guess_2a=int(input("Fxn2 root guess 1 = "))
guess_2b=int(input("Fxn2 root guess 2 = "))
guess_3a=int(input("Fxn3 root guess 1 = "))
guess_3b=int(input("Fxn3 root guess 2 = "))
guess_4a=int(input("Fxn4 root guess 1 = "))
guess_4b=int(input("Fxn4 root guess 2 = "))

root1a = fsolve(fxn1,guess_1a)
root1b = fsolve(fxn1,guess_1b)
root2a = fsolve(fxn2,guess_2a)
root2b = fsolve(fxn2,guess_2b)
root3a = fsolve(fxn3,guess_3a)
root3b = fsolve(fxn3,guess_3b)
root4a = fsolve(fxn4,guess_4a)
root4b = fsolve(fxn4,guess_4b)
c55=0.5*(v_7**2+v_3**2)*rho

print("root 1a = ",root1a)
print("root 1b = ",root1b)
print("root 2a = ",root2a)
print("root 2b = ",root2b)
print("root 3a = ",root3a)
print("root 3b = ",root3b)
print("root 4a = ",root4a)
print("root 4b = ",root4b)

rootlist = [root1a,root1b,root2a,root2b,root3a,root3b,root4a,root4b]
rootlist = [i for i in rootlist if i != 0]
c550 = min(rootlist,key=lambda y:abs(y-c55))[0]
print("Root closest to", c55,":",c550)
input("Press enter to continue finding rest of elastic constants")
#calculating the other constants from c550
c110 = f0 - c550
c330 = f2 - c550
c150a = (c110*c550 - f1)**(1/2)
c350a = (c330*c550 - f3)**(1/2)
f4_1 = c150a + c350a
f4_2 = -c150a + c350a
f4_3 = c150a - c350a
f4_4 = -c150a - c350a

f4_list = [f4_1,f4_2,f4_3,f4_4]
f4_a = min(f4_list,key=lambda z:abs(z-f4))
if f4_a==f4_1:
    c150=c150a
    c350=c350a
elif f4_a==f4_2:
    c150=-c150a
    c350=c350a
elif f4_a==f4_3:
    c150=c150a
    c350=-c350a
elif f4_a==f4_4:
    c150=-c150a
    c350=-c350a

c130 = (-(c150+c350+c550)+((c110+c550+2*c150)*(c550+c330+2*c350)-(v_10*v_11)**2)**(1/2))/rho


a=c550+c440;b=c220+c660;d=c150+c460;l=c110+c660;u=c660+c550;t=c220+c440;h=c350+c460;g=c440+c330
Q=b*(a*l-d)**2-8*rho**3*(v_13*v_14*v_15)**2
G=t*(u*g-h**2)-8*rho**3*(v_16*v_17*v_18)**2
M=b*(a+l)+a*l-d**2-4*rho**2*((v_13*v_14)**2+(v_13*v_15)**2+(v_15*v_14)**2)
N=t*(u+g)+u*g-h**2-4*rho**2*((v_16*v_17)**2+(v_16*v_18)**2+(v_17+v_18)**2)

x=np.linspace(-50,50,num=1000)
Sprime1=2*(-2*h*(N-(x+c460)**2)**(1/2)+((2*h*(x+c460)**2)/(N-(x+c460)**2)**(1/2))+(2*g-2*u)*(x+c460))*(-2*h*(x+c460)*(N-(x+c460)**2)**(1/2)+g*(x+c460)**2+u*(N-(x+c460)**2)-G)+2*(-2*d*(M-(x+c460)**2)**(1/2)+((2*d*(x+c460)**2)/(M-(x+c460)**2)**(1/2))+(x+c460)*(2*l-2*a))*(-2*d*(x+c460)*(M-(x+c460)**2)**(1/2)+l*(x+c460)**2+a*(M-(x+c460)**2)-Q)


fig=plt.figure(1)
ax=fig.add_subplot(1,1,1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(x,Sprime1,label="FxnS")
plt.ylim([-1,1])
plt.legend()
plt.title("Roots of function to find c25")
plt.show()

print("Enter initial x value close to one of the roots:")
guess_a=int(input("Initial root guess = "))

def fxnS(x):
    return(2*(-2*h*(N-(x+c460)**2)**(1/2)+((2*h*(x+c460)**2)/(N-(x+c460)**2)**(1/2))+(2*g-2*u)*(x+c460))*(-2*h*(x+c460)*(N-(x+c460)**2)**(1/2)+g*(x+c460)**2+u*(N-(x+c460)**2)-G)+2*(-2*d*(M-(x+c460)**2)**(1/2)+((2*d*(x+c460)**2)/(M-(x+c460)**2)**(1/2))+(x+c460)*(2*l-2*a))*(-2*d*(x+c460)*(M-(x+c460)**2)**(1/2)+l*(x+c460)**2+a*(M-(x+c460)**2)-Q))

c250 = fsolve(fxnS,guess_a)[0]
c120=(M-(c250+c460)**2)**(1/2)-c660
c230=(N-(c250+c460)**2)**(1/2)-c440

input("Press enter to continue...")

print("c_11 initial estimate = ",c110)
print("c_22 initial estimate = ",c220)
print("c_33 initial estimate = ",c330)
print("c_44 initial estimate = ",c440)
print("c_55 initial estimate = ",c550)
print("c_66 initial estimate = ",c660)
print("c_12 initial estimate = ",c120)
print("c_13 initial estimate = ",c130)
print("c_23 initial estimate = ",c230)
print("c_15 initial estimate = ",c150)
print("c_25 initial estimate = ",c250)
print("c_35 initial estimate = ",c350)
print("c_46 initial estimate = ",c460)

c_11=c110
c_22=c220
c_33=c330
c_44=c440
c_55=c550
c_66=c660
c_12=c120
c_13=c130
c_23=c230
c_15=c150
c_25=c250
c_35=c350
c_46=c460

'''
#finding velocities from elastic constants - cubic roots
p1=(c660+c220)*(c550+c440)+(c110+c660)*(c660+c220+c550+c440)-(c150+c460)**2-(c120+c660)**2-(c250+c460)**2
p2=(1/8)*((c110+c660)*(c660+c220)*(c550+c440)+2*(c150+c460)*(c250+c460)*(c120+c660)-(c110+c660)*(c250+c460)**2-(c550+c440)*(c120+c660)**2-(c150+c460)**2*(c660+c220))
p3=c660+(1/2)*(c110+c220+c550+c440)
p4=(c220+2*c440+c330)*(c660+c550)+(c220+c440)*(c440+c330)-(c460+c350)**2-(c460+c250)**2-(c440+c230)**2
p5=(1/8)*((c660+c550)*(c220+c440)*(c440+c330)+2*(c460+c350)*(c440+c230)*(c460+c250)-(c660+c550)*(c440+c230)**2-(c460+c250)**2*(c440+c330)-(c220+c440)*(c460+c350)**2)
p6=c440+(1/2)*(c660+c550+c220+c330)

coeff=[1,-p3,(p1/4),-p2]
roots=np.roots(coeff)
X1=(roots[np.isreal(roots)]).real
X10=X1[0];X11=X1[1];X12=X1[2]
V10=(X10/rho)**(1/2);V11=(X11/rho)**(1/2);V12=(X12/rho)**(1/2)
coeff2=[1,-p6,(p4/4),-p5]
roots2=np.roots(coeff2)
X2=(roots2[np.isreal(roots2)]).real
X20=X2[0];X21=X2[1];X22=X2[2]
V20=(X20/rho)**(1/2);V21=(X21/rho)**(1/2);V22=(X22/rho)**(1/2)


#working backwards to find associated velocities
v1a=((c110+c550+((c110-c550)**2+4*c150**2)**(1/2))/(2*rho))**(1/2)
v2a=(c660/rho)**(1/2)
v3a=((c110+c550-((c110-c550)**2+4*c150**2)**(1/2))/(2*rho))**(1/2)
v4a=((c440+c660+((c440-c660)**2+4*c460**2)**(1/2))/(2*rho))**(1/2)
v5a=(c220/rho)**(1/2)
v6a=((c440+c660-((c440-c660)**2+4*c460**2)**(1/2))/(2*rho))**(1/2)
v7a=((c330+c550-((c330-c550)**2+4*c350**2)**(1/2))/(2*rho))**(1/2)
v8a=(c440/rho)**(1/2)
v9a=((c330+c550+((c330-c550)**2+4*c350**2)**(1/2))/(2*rho))**(1/2)
v10a=((c110+c330+2*(c150+c350+c550)-(4*(2*c350+c330-c110-2*c150)**2+(c150+c350+c130+c550)**2)**(1/2))/(4*rho))**(1/2)
v11a=((c660+c440+2*c460)/(2*rho))*(1/2)
v12a=((c110+c330+2*(c150+c350+c550)+(4*(2*c350+c330-c110-2*c150)**2+(c150+c350+c130+c550)**2)**(1/2))/(4*rho))**(1/2)
v13a=V10
v14a=V12
v15a=V11
v16a=V20
v17a=V22
v18a=V21

print(v1a,v2a,v3a,v4a,v5a,v6a,v7a,v9a,v10a,v11a,v12a,v13a,v14a,v15a,v16a,v17a,v18a)
print(v_1,v_2,v_3,v_4,v_5,v_6,v_7,v_8,v_9,v_10,v_11,v_12,v_13,v_14,v_15,v_16,v_17,v_18)
'''
print("After optimization process, elastic constants:")

#define c_ij = cij0 for elastic constants you want fixed during optimization process:

def F(params):
    c_11,c_22,c_33,c_44,c_55,c_66,c_12,c_13,c_23,c_15,c_25,c_35,c_46=params
    Fv1=((c_11+c_55+((c_11-c_55)**2+4*c_15**2)**(1/2))/(2*rho))**(1/2)-v_1
    Fv2=(c_66/rho)**(1/2)-v_2
    Fv3=((c_11+c_55-((c_11-c_55)**2+4*c_15**2)**(1/2))/(2*rho))**(1/2)-v_3
    Fv4=((c_44+c_66+((c_44-c_66)**2+4*c_46**2)**(1/2))/(2*rho))**(1/2)-v_4
    Fv5=(c_22/rho)**(1/2)-v_5
    Fv6=((c_44+c_66-((c_44-c_66)**2+4*c_46**2)**(1/2))/(2*rho))**(1/2)-v_6
    Fv7=((c_33+c_55-((c_33-c_55)**2+4*c_35**2)**(1/2))/(2*rho))**(1/2)-v_7
    Fv8=(c_44/rho)**(1/2)-v_8
    Fv9=((c_33+c_55+((c_33-c_55)**2+4*c_35**2)**(1/2))/(2*rho))**(1/2)-v_9
    Fv10=((c_11+c_33+2*(c_15+c_35+c_55)-(4*(2*c_35+c_33-c_11-2*c_15)**2+(c_15+c_35+c_13+c_55)**2)**(1/2))/(4*rho))**(1/2)-v_10
    Fv11=((c_66+c_44+2*c_46)/(2*rho))*(1/2)-v_11
    Fv12=((c_11+c_33+2*(c_15+c_35+c_55)+(4*(2*c_35+c_33-c_11-2*c_15)**2+(c_15+c_35+c_13+c_55)**2)**(1/2))/(4*rho))**(1/2)-v_12
    
    p1=(c_66+c_22)*(c_55+c_44)+(c_11+c_66)*(c_66+c_22+c_55+c_44)-(c_15+c_46)**2-(c_12+c_66)**2-(c_25+c_46)**2
    p2=(1/8)*((c_11+c_66)*(c_66+c_22)*(c_55+c_44)+2*(c_15+c_46)*(c_25+c_46)*(c_12+c_66)-(c_11+c_66)*(c_25+c_46)**2-(c_55+c_44)*(c_12+c_66)**2-(c_15+c_46)**2*(c_66+c_22))
    p3=c_66+(1/2)*(c_11+c_22+c_55+c_44)
    p4=(c_22+2*c_44+c_33)*(c_66+c_55)+(c_22+c_44)*(c_44+c_33)-(c_46+c_35)**2-(c_46+c_25)**2-(c_44+c_23)**2
    p5=(1/8)*((c_66+c_55)*(c_22+c_44)*(c_44+c_33)+2*(c_46+c_35)*(c_44+c_23)*(c_46+c_25)-(c_66+c_55)*(c_44+c_23)**2-(c_46+c_25)**2*(c_44+c_33)-(c_22+c_44)*(c_46+c_35)**2)
    p6=c_44+(1/2)*(c_66+c_55+c_22+c_33)

    coeff=[1,-p3,(p1/4),-p2]
    roots=np.roots(coeff)
    X1=(roots[np.isreal(roots)]).real
    X10=X1[0];X11=X1[1];X12=X1[2]
    Fv13=(X10/rho)**(1/2)-v_13;Fv15=(X11/rho)**(1/2)-v_15;Fv14=(X12/rho)**(1/2)-v_14

    coeff2=[1,-p6,(p4/4),-p5]
    roots2=np.roots(coeff2)
    X2=(roots2[np.isreal(roots2)]).real
    X20=X2[0];X21=X2[1];X22=X2[2]
    Fv16=(X20/rho)**(1/2)-v_16;Fv18=(X21/rho)**(1/2)-v_18;Fv17=(X22/rho)**(1/2)-v_17
    
    return((pk1*(Fv1)**2)+(pk2*(Fv2)**2)+(pk3*(Fv3)**2)+(pk4*(Fv4)**2)+(pk5*(Fv5)**2)+(pk6*(Fv6)**2)+(pk7*(Fv7)**2)+(pk8*(Fv8)**2)+(pk9*(Fv9)**2)+(pk10*(Fv10)**2)+(pk11*(Fv11)**2)+(pk12*(Fv12)**2)+(pk13*(Fv13)**2)+(pk14*(Fv14)**2)+(pk15*(Fv15)**2)+(pk16*(Fv16)**2)+(pk17*(Fv17)**2)+(pk18*(Fv18)**2))

initial_guess=[c110,c220,c330,c440,c550,c660,c120,c130,c230,c150,c250,c350,c460]
result=minimize(F,initial_guess)
if result.success:
    fitted_params = result.x
    print(fitted_params)
else:
    print("Poor result, fitted parameters:",result.x)


#Checking results using crystal stability conditions (found in Dunk 1987)
if c_11>0 and c_22>0 and c_33>0 and c_44>0 and c_55>0 and c_66>0:
    print('Crystal stability check #1 successful (c_ii>0)')
else:
    print('Crystal stability check #1 unsuccessful (c_ii>0)')

if c_11*c_22 > c_12**2 and c_11*c_33 > c_13**2 and c_11*c_55 > c_15**2 and c_22*c_33>c_23**2 and c_22*c_55>c_25**2 and c_33*c_55>c_35**2 and c_44*c_66>c_46**2:
    print('Crystal stability check #2 successful (c_ii*cjj>cij**2)')
else:
    print('Crystal stability check #2 unsuccessful (c_ii*cjj>cij**2)')

if c_11*c_22*c_33+2*c_12*c_13*c_23>c_11*c_23**2+c_33*c_12**2+c_22*c_12**2 and c_11*c_22*c_55+2*c_12*c_25*c_12>c_11*c_25**2+c_55*c_12**2+c_22*c_15**2 and c_22*c_33*c_55+2*c_23*c_35*c_25>c_22*c_35**2+c_55*c_23**2+c_33*c_35**2:
    print('Crystal stability check #3 successful (c_ii*c_jj*c_kk + 2*c_ij*c_ik*c_jk > c_ii*c_jk**2 + c_kk*c_ij**2 + c_jj*c_ik**2)')
else:
    print('Crystal stability check #3 unsuccessful (c_ii*c_jj*c_kk + 2*c_ij*c_ik*c_jk > c_ii*c_jk**2 + c_kk*c_ij**2 + c_jj*c_ik**2)')