import math
import numpy as np
import pandas as pd

def ExactMCF(r0,t,alpha):
    alpha=1
    return np.sqrt(r0**2-2.0*alpha*t)

def PhiEq_1D(x,r,eps):
    return 0.5*(1-np.tanh(3.0*((abs(x)-r)/eps)))


def PhiEq_x0(x,x0,eps):
    return 0.5*(1-np.tanh(3.0*(x-x0)/eps))

def PhiCircle(xm,ym,r,eps):
    return 0.5*(1-np.tanh(3.0*((xm[:,:]**2+ym[:,:]**2)**0.5-r)/eps))

 
def DBC_2ghostpoints_1D(u,lv,rv):
    u[0]=lv
    u[1]=lv
    u[u.shape[0]-2]=rv
    u[u.shape[0]-1]=rv
    return u

#Option for DBC as at point (ii) in 2D
def DBC_2ghostpoints_2D(u,lv,rv,label='x'):
    if(label=='x'):
        u[0,:]=lv
        u[1,:]=lv
        u[u.shape[0]-2,:]=lv
        u[u.shape[0]-1,:]=lv
    if(label=='y'):
        u[:,0]=lv
        u[:,1]=lv
        u[:,u.shape[1]-2]=lv
        u[:,u.shape[1]-1]=lv
    return u

def DBC(i,n):
        if (i<2):
                return 0
        elif (i>n-3):
                return n-1
        else:
                return int(i)
def PBC(i,n):
        if (i<0):
                return int(i+n-1)
        elif (i>n-2):
                return int(i-n+1)
        else:
                return int(i)

# Numerical approximation of  1st order
def d1f_5pst(u,ds,label='x'):
    n=u.shape[0]
    d1f=np.zeros(u.shape)
    for i in range(0,n):
        for j in range(0,n):
            if(label=='x'):
                d1f[i,j]=(-u[PBC(i+2,n),j]+8*u[PBC(i+1,n),j]-8*u[PBC(i-1,n),j]+u[PBC(i-2,n),j])/(12*ds)
            if(label=='y'):
                d1f[i,j]=(-u[i,PBC(j+2,n)]+8*u[i,PBC(j+1,n)]-8*u[i,PBC(j-1,n)]+u[i,PBC(j-2,n)])/(12*ds)
    return d1f



# Numerical approximation of  2nd order

def d2f_5pst(u,dx,label='x'):
    n=u.shape[0]
    d2f=np.zeros(u.shape)
    for i in range(0,n):
        for j in range(0,n):
            if(label=='x1D'):
                d2f[i]=(-u[DBC(i+2,n)]+16*u[DBC(i+1,n)]-30*u[i]+16*u[DBC(i-1,n)]-u[DBC(i-2,n)])/(12*dx*dx)
            if(label=='x'):
                d2f[i,j]=(-u[PBC(i+2,n),j]+16*u[PBC(i+1,n),j]-30*u[i,j]+16*u[PBC(i-1,n),j]-u[PBC(i-2,n),j])/(12*dx*dx)
            if(label=='y'):
                d2f[i,j]=(-u[i,PBC(j+2,n)]+16*u[i,PBC(j+1,n)]-30*u[i,j]+16*u[i,PBC(j-1,n)]-u[i,PBC(j-2,n)])/(12*dx*dx)
    return d2f


def Bphi(phi):
    return 18*(phi**2)*(1-phi)**2

def d1Bphi(phi):
    return 18*(2*phi-6*phi**2+4*phi**3)

def d2Bphi(phi):
    return 18*(2-12*phi+12*phi**2)





# Definition and initialization of psi

def psi(phi,ds):
    #print(phi.shape)
    angle =  np.arctan(d1f_5pst(phi,ds,label = 'y')/d1f_5pst(phi,ds,label = 'x'))
    angle = np.nan_to_num(angle,nan= 0.000001)
    return angle
    #psi = arctan(gradyphi/gradxphi)
    


# Regularisation parameter for curve smoothning  

def Rnz(z):
    neta = 0.5 
    return neta/5*(np.log(2) + np.log(1+np.cosh(5*z/neta)))
    # Rn = neta/5*(ln(2) + ln(1+cos(5*z/neta)))

# gamma nhat

def gamma(phi,ds):
    gamma2 = 1
    gamma1 = 2
    angle = psi(phi,ds)
    var1 = np.cos(angle)
    var2 = np.sin(angle)
    #print(z1, z2)
    return gamma2*Rnz(var1) + gamma1* Rnz(var2)
    # gamma2 * Rn(cos(psi)) + gamma1 * Rn(sin(psi))

#def something(phi):
#    gamma = Interface_Energy_Density(phi)
#    return gamma

# Mhat  

def Mhat(phi,ds):
    M1 = 1 
    M2 =  1
    angle = psi(phi,ds)
    return M1 * np.cos(angle) **2 + M2* np.sin(angle)**2
    #print(m4.shape)
    # Mhat = M1* cos**2(psi) + M2* sin**2(psi)


# Function that returns 1st order FD numerical approximation for the parameter passed in x and y direction. 
def d1x_d1y_phi(phi,ds,label = 'x'):
    if (label == 'x'):
        d1xphi = d1f_5pst(phi,ds,label = 'x')
        #print("This value is d1phix.shape",d1xphi.shape)
        #print("This value is d1phix[0,:]",d1xphi[0,:])
        return d1xphi
    if (label == 'y'):
        d1yphi = d1f_5pst(phi,ds,label = 'y')
        #print("This value is d2phiy.shape",d1yphi.shape)
        #print("This value is d2phiy[0,:]",d1yphi[0,:])
        return d1yphi
    
    
# Function that returns 2nd order FD numerical approximation for the parameter passed in x and y direction.    
def d2xphi_d2yphi(phi,dx,label = 'x'):
    if (label == 'x'):
        d2phix = d2f_5pst(phi, dx, label = 'x')
        return d2phix
    
    if (label == 'y'):
        d2phiy = d2f_5pst(phi, dx, label = 'y')
        return d2phiy

    
# Function returns divergance of gamma and gradient of phi    
def  Divgamma_grad_psi(phi,ds):
    #div_grad_phi = np.zeros(phi.shape)
    d1xphi = d1f_5pst(phi,ds, label = 'x')
    #print(d1xphi ,"This is d1xphi")
    d1yphi = d1f_5pst(phi,ds, label = 'y')
    #print(d1xphi ,"This is d1yphi")
    return d1f_5pst((gamma(phi,ds)*d1xphi),ds,'x')+d1f_5pst((gamma(phi,ds)*d1yphi),ds,'y')  
    # d( gamma * gradxphi(x,y) ) / dx + d( gamma * grady(x,y)) / dy


# Function returns differentiation of the gamma function defined earlier. This function is used to calculate
# the A and B terms.      
def d1fGamma_psi(phi,ds):
    gamma2 = 1
    gamma1 = 2
    #print('###')
    #print(phi)
    #print(phi.shape)
    angle = psi(phi,ds)
    c1 = -np.sin(angle)
    c2 = np.cos(angle)
    return gamma2*Rnz(c1) + gamma1* Rnz(c2)
    # d1fGamma_psi = dgamma/dpsi
    
    
    
# Function for A term in the formula 
def psi_phi_xcomp(phi,ds):
    xcomp = d1fGamma_psi(phi,ds)*(-d1x_d1y_phi(phi,ds,label = 'y')/((d1x_d1y_phi(phi, ds,label = 'x'))**2 + (d1x_d1y_phi(phi, ds,label = 'y'))**2))
    x_comp = np.nan_to_num(xcomp,nan= 0.000001)
    return x_comp
    # A = dgamma/dpsi* -gradphiy/gradphix**2 + gradphiy**2
    
# Function for B term in the formula
def psi_phi_ycomp(phi,ds):
    ycomp = d1fGamma_psi(phi,ds)*(d1x_d1y_phi(phi,ds,label = 'x')/((d1x_d1y_phi(phi, ds, label = 'x'))**2 + (d1x_d1y_phi(phi, ds, label = 'y'))**2))
    y_comp = np.nan_to_num(ycomp,nan= 0.000001)
    return y_comp
    # B = dgamma/dpsi* gradphix/gradphix**2 +  gradphiy**2



# Function returns divergence of 2nd order num approx with A and B terms
def Div_2_AB(phi,ds,eps):
    return eps*(d1f_5pst(((d1x_d1y_phi(phi,ds,label = 'x') + d1x_d1y_phi(phi,ds,label = 'y')) * psi_phi_xcomp(phi,ds)),ds,'x') + d1f_5pst(((d1x_d1y_phi(phi,ds,label = 'x') + d1x_d1y_phi(phi,ds,label = 'y')) * psi_phi_ycomp(phi,ds)),ds,'y'))
    #Div_2_AB = eps*(dx(|gradxphi|*A) + dy(|gradyphi|*B))
    
    
# Function returns the complete intergartion scheme for anisotropy    
def  Anisotropy_integration_scheme(phi,eps,dt,ds):
    return (1/(1+eps**-2* Mhat(phi,ds)*gamma(phi,ds)*dt*d2Bphi(phi)))*(phi+eps**-1*Mhat(phi,ds)*dt*((-eps**-1*gamma(phi,ds)*d1Bphi(phi))+(eps**-1*gamma(phi,ds)*d2Bphi(phi)*phi) + eps*Divgamma_grad_psi(phi,ds)) + Div_2_AB(phi,ds,eps))
    
    #print('this is result' ,result[0,:])
    