# -*- coding: utf-8 -*-
"""
Created on Thu May 13 23:42:23 2021

@author: shivd
"""
#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import math
import sys
from scipy.optimize import curve_fit
from functions import ExactMCF
from functions import PhiCircle
from functions import PhiEq_1D
from functions import IntegrationScheme
from functions import DBC_2ghostpoints_2D
from mpltools import annotation


def default_plotting():
    plt.figure(1,figsize=[10,8])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

def timevec(dt,T):
    return np.linspace(0,T,int(T/dt))

def mesh(L,ds):
    x=np.arange(-(L), (L), ds)
    y=x
    return np.meshgrid(x, y)

def update_mesh_function(r,L,eps,ds):
    xm, ym = mesh(L,ds)
    phi = PhiCircle(xm,ym,r,eps)
    return xm,ym,phi


def AllenCahn_FD(phi,eps,ds,dt,L,Tmax,label,print_steps=False):
    #
    #label select the information to extract, namely
    #'n_plot': plot some plots during the time iteration (aux_variable is the number of plots)
    #'max_over_time': max(u) as function of t
    #'int_u': integral of u as function of t
    #
    t=0
    j=0
    ntsteps=int(Tmax/dt)
    vecout=np.zeros(ntsteps)
    nplot=100
    #i
    xfit=np.linspace(0,L,int(xm.shape[0]/2))
    for i in range(0,ntsteps):

        #optional printout of the progress
        if(print_steps==True):
            print(i,'/',ntsteps)

        #iteration
        if(i==0):   
            u=phi  #initialization of phi
        else:       
            u=IntegrationScheme(u,eps,ds,dt)
   
        if(label=='r_t'):

            #I exploit the function PhiEq_1D(x,r,eps) to fit the data and extract the best r associated to the current phi distribution
            #As input i take (with no loss of generality the positive x axis and the associated phi values
            #plt.plot(xfit,yfit)
            #plt.plot(xfit,PhiEq_1D(xfit,popt[0],popt[1]))
            #plt.show()
            #print(popt)

            yfit=u[int(xm.shape[0]/2)-1,int(xm.shape[0]/2):int(xm.shape[0])]
            if(i==0):
                popt0, pcov0 = curve_fit(PhiEq_1D, xfit, yfit)
                popt=popt0
            else:
                popt, pcov = curve_fit(PhiEq_1D, xfit, yfit)
            vecout[i]=popt[0]
        #set output type
        if(label=='n_plot'):
            if((i) % nplot == 0):
                plt.figure(1,figsize=[8,6])
                plt.rc('text', usetex=True)
                plt.rc('font', family='serif',size=16)
                plt.title('Progress: '+"{:.2f}".format(t)+'/'+"{:.2f}".format(Tmax))
                val=1
                levelvec=np.linspace(0,val,256,endpoint=False)
                plt.contourf(xm,ym,u[:,:],levels = np.append(levelvec,[val]),cmap='seismic')
                cbar=plt.colorbar()
                cbar.set_ticks(np.append(np.linspace(0,val,10,endpoint=False),val))
                plt.ylabel('$y$',size=20)
                plt.xlabel('$x$',size=20)
                plt.savefig('AllenCahn_MCF'+str(j)+'.png')
                plt.show()
                plt.close()
                j=j+1

        t=t+dt
    #
    return vecout

#Radius of the Circular domain Omega
r=4.0                              
L=2*r #Size L_OmegaDDA

#Parameters for the Initial Condition
eps=1.0





ds=0.2
dt=0.001
eps=1.0
T=((r**2)/2.0)/1.0
xm,ym,phi=update_mesh_function(r,L,eps,ds)
AllenCahn_FD(phi,eps,ds,dt,L,T,label='n_plot',print_steps=True)


