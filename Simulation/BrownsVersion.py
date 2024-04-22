#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 14:30:59 2024

Browns version of CV

@author: elliothowell
"""

import numpy as np
import matplotlib.pyplot as plt
import math
c0ox= 6.1E-08 #mol/cm^3
D= 1e-5 #[cm^2/s]
k0=0.1 #[cm/s]
n=1
alpha= 0.5
Eeq= 0.0 #[V vs NHE]
F= 96485 #[C/mol]
R= 8.314 #J/(mol*K)
T= 298 #K
f=n*F/(R*T) # [1/V]
l=500
v=0.1 #V/s
Ei=1 #V
Ef=-1
Efwd= np.linspace(Ei, Ef, int(l/2)) #potential range
Erev= np.linspace(Ef, Ei, int(l/2))
E= np.concatenate((Efwd,Erev))
#E=np.linspace(Ei, Ef, l)
eta= E-Eeq
dt= np.abs(Ef-Ei)/v/l #s
dx= 10e-4 #[cm]
Dm= D*dt/dx**2
iteration= np.arange(0,51,1)
cR= np.zeros((l+1,100))
cO= np.zeros((l+1,100))
cO[:,-1]= c0ox
cO[0,:]= c0ox
Jox=[]
Jred=[]
J=[]
#
testList1 = []
testList2 = []

for k in range(1,l+1):
    kf= k0*math.exp(-alpha*f*eta[k-1])
    kb= k0*math.exp((1-alpha)*f*eta[k-1])
    jox=-(kf*cO[k-1,1]-kb*cR[k-1,1])/(1+(kf*dx/D)+(kb*dx/D)) #different from bard
    testList1.append(kf*cO[k-1,1]-kb*cR[k-1,1])
    testList2.append((1+(kf*dx/D)+(kb*dx/D)))
    jred=-jox
    Jox.append(jox)
    Jred.append(jred)
    cO0= cO[k-1,1]+jox*dx/D #one extra space away from boundry
    cR0= cR[k-1,1]+jred*dx/D
    # print(f"kf: {np.round(kf*cO[k-1,1],3)}, kb: {np.round(kb*cR[k-1,1],3)}")
    if cO0<0:
        cO0=0
    if cR0<0:
        cR0=0
    cO[k,0]= cO0
    cR[k,0]= cR0
    J.append(-n*F*jox)
    for j in range(1,99):
        cO[k,j]= cO[k-1,j] + Dm*(cO[k-1, j+1]-2*cO[k-1,j]+cO[k-1,j-1])
        # cO[k,j]=cox
        cR[k,j]= cR[k-1,j] + Dm*(cR[k-1, j+1]-2*cR[k-1,j]+cR[k-1,j-1])
        # cR[k,j]=cred
plt.plot(E,testList1)
# plt.plot(E, cO[1:,0])
# plt.plot(E, cR[1:,0])