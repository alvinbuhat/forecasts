# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 18:49:09 2020

@author: MONICA, ALVIN, BEN PAUL, EDD FRANCIS, ELEANOR, JONATHAN 
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame

R0=1.5811
sigma=0.0130
gamma=0.0069
tau=14
d=0.0021

def StateFunc(t,x,R0,sigma,gamma,d,tau,dt):
    xprime=np.zeros(4)
    S0=5000000
    E0=10000
    I0=1
    beta=(R0/tau)*(S0+E0+I0)/S0
    xprime[0]=-beta*x[0]*x[2]/(x[0]+x[1]+x[2]+x[3])
    xprime[1]=beta*x[0]*x[2]/(x[0]+x[1]+x[2]+x[3])-sigma*x[1]
    xprime[2]=sigma*x[1]-gamma*x[2]-d*x[2]
    xprime[3]=gamma*x[2]
    return xprime


def RK4ForwardState(t,x,R0,sigma,gamma,d,tau,dt,N):
    for i in range(1,N):
        k1=StateFunc(t[i-1],x[i-1],R0,sigma,gamma,d,tau,dt)
        k2=StateFunc(t[i-1]+(dt/2),x[i-1]+(dt/2)*k1,R0,sigma,gamma,d,tau,dt)
        k3=StateFunc(t[i-1]+(dt/2),x[i-1]+(dt/2)*k2,R0,sigma,gamma,d,tau,dt)
        k4=StateFunc(t[i-1]+dt,x[i-1]+dt*k3,R0,sigma,gamma,d,tau,dt)
        x[i]=x[i-1]+(dt/6)*(k1+2*k2+2*k3+k4)
    return x
    
def cov(T):
    test=-1
    tol=0.001
    N=8200
    t=np.linspace(0,T,N)
    dt=T/N
    
    x=np.zeros((N,4))
    S0=5000000
    E0=10000
    I0=1
    x[0,0]=S0
    x[0,1]=E0
    x[0,2]=I0
    x[0,3]=0
#==============================================================================
    while (test<0):
        oldx=x
        
        x=RK4ForwardState(t,x,R0,sigma,gamma,d,tau,dt,N)
        

        temp2=tol*np.sum(np.abs(x))-np.sum(np.abs(oldx-x))       
        test=temp2
 
#==============================================================================

    cases=[2,3,7,21,30,46,49,61,108,137,139,184,199,214,227,304,377,458,548,632,702,804,1071,1414,1541,2079,2306,2628,3013, \
    3088,3240,3653,3757,3863,4069,4188,4421,4641,4925,5215,5445,5652,5870,6078,6249,6446,6582,6692,6963,7174,7276, \
    7561,7758,7939,8192,8466,8747,8902,9197,9459,9658,9978,10315,10434,10581,10764,11052,11316,11583,11840,12054, \
    12268,12475,12680,12903,13179,13390,13553,13732,13989,14271,14621,14997,15534,16576,17162,18018,18566,18924,19670, \
    20300,20541,21246,21796,22369,22880,23608,24050,24658,25261,25796,26279,26642,27098,27658,28312,29249,29898, \
    30526,31669,32138,32909,33908,34644,35288,36257,37308,38159,38414,39905,41287,43626,45626,47118,49599,50957, \
    52141,53483,55567,56396,57016,58317,60749,62453,64719,66935,68450,70382,71908,74088,76163,78157,80232,81880,83508,85376,89288, \
    93351,98200,103127,106282,112269,115692,119188,122512,126686,129727,136620,139540,143903,147842,153998,158305]
    
 
    t1=np.linspace(0,162,163)
##==============================================================================
    plt.subplots_adjust(top=0.8,hspace=0.5, wspace=0.3)

    plt.plot(t+1,x[:,2],linewidth=3,label="Model Output")
    plt.scatter(t1+1,cases,color="red",label="Data",marker=10)
    plt.ylabel('Cumulative Cases',fontsize=12)   
    plt.legend(loc="best",bbox_to_anchor=(0.25,0.95))
    plt.xlabel('Days',fontsize=12)

 
    plt.show()
    
    df=DataFrame({'Cumu':x[:,2],'Recov':x[:,3]}) 
    df.to_excel('OUTPUT.xlsx')  


  