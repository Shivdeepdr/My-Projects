def ExactMCF(r0,t,alpha):
    import numpy as np
    alpha=1
    return np.sqrt(r0**2-2.0*alpha*t)

def PhiEq_1D(x,r,eps):
    import numpy as np
    return 0.5*(1-np.tanh(3.0*((abs(x)-r)/eps)))


def PhiEq_x0(x,x0,eps):
    import numpy as np
    return 0.5*(1-np.tanh(3.0*(x-x0)/eps))

def PhiCircle(xm,ym,r,eps):
    import numpy as np
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

def d2f_5pst(u,dx,label='x1D'):
    import numpy as np
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

def IntegrationScheme(phi,eps,ds,dt):
    factor=1.0/(1.0-(dt/eps)*d2Bphi(phi))
    return factor*(phi+dt*(eps*(d2f_5pst(phi,ds,label='x')+d2f_5pst(phi,ds,label='y'))-(1.0/eps)*d1Bphi(phi)-(phi/eps)*d2Bphi(phi)))

def IntegrationScheme1D(phi,eps,ds,dt):
    factor=1.0/(1.0-(dt/eps)*d2Bphi(phi))
    return factor*(phi+dt*(eps*(d2f_5pst(phi,ds,label='x1D'))-(1.0/eps)*d1Bphi(phi)-(phi/eps)*d2Bphi(phi)))
        





