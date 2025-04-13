# In[ ]:
# import necessary modules
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d, PchipInterpolator
import math
from scipy.integrate import solve_ivp
import pandas as pd
from numpy import log, exp
Gamma=math.gamma

def table(f,t):
    tablef=np.zeros(t.shape[0])
    for i in np.arange(t.shape[0]):
        tablef[i]=f(t[i])
    return tablef

from background import rhodm, rhob, rhor, Ht, alpha, photon_mfp_full, cb2_full, a_rec, a_eq, a_background, rho_tot
omegam=0.12038+0.022032
omegab=0.022032/omegam
omegadm=0.12038/omegam

Hsim=Ht#lambda a: Ht(0.01)*np.sqrt(rhob(a)+rhodm(a)+rhor(a))/np.sqrt(rhob(0.01)+rhodm(0.01)+rhor(0.01))
rho_tot2=lambda a:rhob(a)+rhodm(a)+rhor(a)
rhom=lambda a:rhob(a)+rhodm(a)
#In[]
#scale invariant simplified perturbation eqn post recombination
def ddel_baryon(a,delta): #DM and baryon perturbation equations
        dldm,thdm,dlb,thb=delta[0],delta[1],delta[2],delta[3]
        k2phi=-3./2*(a*Hsim(a))**2*(rhob(a)*dlb+rhodm(a)*dldm)/rho_tot2(a)
        del_dldm=-(thdm/a**2/Hsim(a))
        del_thdm=-thdm/a+k2phi/a**2/Hsim(a)
        del_thb=-thb/a+k2phi/a**2/Hsim(a)+Hsim(a)*rhom(a)/rho_tot(a)*1 #the S0 term captures the impact of PMF on delta_b. S0 is taken to be a function of kd, which in turn is a function of a.
        del_dlb=-(thb/a**2/Hsim(a))

        return [del_dldm,del_thdm,del_dlb,del_thb]

def ddel_baryon_grav(a,delta,k): #DM and baryon perturbation equations
        phi,dldm,thdm,dlb,thb=delta[0],delta[1],delta[2],delta[3],delta[4]
        del_phi=-phi/a+(-k**2*phi-3./2*(a*Ht(a))**2*(rhob(a)*dlb+rhodm(a)*dldm)/rho_tot(a))/(3*(a*Ht(a))**2)/a
        del_dldm=-(thdm/a**2/Ht(a) - 3*del_phi)
        del_thdm=-thdm/a+k**2/a**2/Ht(a)*phi
        del_thb=-(1+alpha(a)/Ht(a))*thb/a+k**2/a**2/Ht(a)*phi+Ht(a)*0 #the S0 term captures the impact of PMF on delta_b. S0 is taken to be a function of kd, which in turn is a function of a.
        del_dlb=-(thb/a**2/Ht(a) - 3*del_phi)

        return [del_phi,del_dldm,del_thdm,del_dlb,del_thb]

delta_b0=0
delta_dm0=0
theta_b0=0#theta_b_pmf[i_100]/(norm)/0.817
theta_dm0=0#-(delta_dm_pmf[i_100]-delta_dm_pmf[i_100-1])/(a_pmf[i_100]-a_pmf[i_100-1])*0.01**2*Ht(0.01)/(norm)/1.62

deltasolve0=[delta_dm0,theta_dm0,delta_b0,theta_b0] #initial condition assumes everythong zero as perturbations are sourced by PMFs. I start my computation from a_mfp at which point silk damping should have erasedd any earlier evolution.
deltasolve=solve_ivp(lambda a,y: ddel_baryon(a,y),[a_rec,1],deltasolve0,method='BDF',dense_output=True,atol=1e-6,rtol=1e-5)

#In[]:
#matter power spectrum from PMFs
def G(nB):
    if nB==-2.5:
        GnB=0.446299
    elif nB==-2:
        GnB=1.178097
    elif nB==-2.9:
        GnB=0.0684258
    elif nB==-2.95:
         GnB=0.0284859
    elif nB==-1.6:
        GnB=5.185
    elif nB==-1.8:
        GnB=1.89
    elif nB==-2.2:
        GnB=0.816046
    else:
        GnB=5
    return GnB
k1mpc=1/0.678

kappa=lambda nB, eta : (G(nB)/1.14/eta)**0.25
kd=lambda B1mpc, nB, eta: (0.1*kappa(nB,eta)*B1mpc)**(-2/(nB+5))*k1mpc
kd_i=lambda B1mpc, nB: (0.01*np.sqrt(nB+5)*B1mpc)**(-2/(nB+5))*k1mpc
#PS0_an=lambda k,B1mpc,nB,eta: (2*np.pi**2/k**3)*10**-4*(k/k1mpc)**(2*nB+10)*B1mpc**4*G(nB)*exp(-2*k**2/kd(B1mpc,nB,eta)**2)
def PS0_an(k,B1mpc,nB,eta):
    if nB<-1.5:
        ans=(2*np.pi**2/k**3)*0.918*10**-4*(k/k1mpc)**(2*nB+10)*B1mpc**4*G(nB)*exp(-2*k**2/kd(B1mpc,nB,eta)**2)
    else:
        ans=(2*np.pi**2/k**3)*(k/k1mpc)**(7)*(0.1*B1mpc)**(14/(5+nB))*exp(-2*k**2/kd(B1mpc,nB,eta)**2)
    
    return ans

Pb_an=lambda a,k,B1mpc,nB,eta: deltasolve.sol(a)[2]**2*PS0_an(k,B1mpc,nB,eta)
Pdm_an=lambda a,k,B1mpc,nB,eta: deltasolve.sol(a)[0]**2*PS0_an(k,B1mpc,nB,eta)
Pm_an=lambda a,k,B1mpc,nB,eta: (omegab*deltasolve.sol(a)[2]+omegadm*deltasolve.sol(a)[0])**2*PS0_an(k,B1mpc,nB,eta)
PB_an=lambda k,B1mpc,nB,eta: 10**-18*4*np.pi**2/Gamma((nB+3)/2)*B1mpc**2/k1mpc**3*(k/k1mpc)**(nB)*exp(-k**2/kd(B1mpc,nB,eta)**2)

xi_b=lambda a: deltasolve.sol(a)[2]
xi_dm=lambda a: deltasolve.sol(a)[0]
xi_m=lambda a: omegab*deltasolve.sol(a)[2]+omegadm*deltasolve.sol(a)[0]
fb_ratio=lambda a: deltasolve.sol(a)[2]/(omegab*deltasolve.sol(a)[2]+omegadm*deltasolve.sol(a)[0])

M_D=lambda B1mpc, nB, eta:5*10**12*(1/kd(B1mpc,nB,eta))**3 #provides the damping mass. Our results breakdown for masses below this.

def findB(nB, eta, Box, Grid):
    knyq=np.pi*Grid/Box #in h/Mpc
    B1mpc_ans=(1.5*knyq)**(-(nB+5)/2)*10/kappa(nB,eta)
    return B1mpc_ans

def findBox(nB, eta, Box, Grid):
    knyq=np.pi*Grid/Box #in h/Mpc
    B1mpc_ans=(1.5*knyq)**(-(nB+5)/2)*10/kappa(nB,eta)
    return B1mpc_ans

#import colibri.cosmology as cc
#C = cc.cosmo()
Pk_m=np.loadtxt('/scratch/pralegan/ngenic_pmf/PM_lcdm.dat')
K_i=Pk_m[:,0]
PM_i=Pk_m[:,1]

def Pm_tot(a,K,B1mpc,nB,eta):
    #kk,pk_class=C.class_Pk(k=K,z=1/a-1)
    pk_extrapolate=np.exp(np.interp(np.log(K),np.log(K_i),np.log(PM_i)))*(a/0.01)**2
    if (B1mpc==0):
        pk_pmf=0*pk_extrapolate
    else:
        pk_pmf=table(lambda x: Pm_an(a,x,B1mpc,nB,eta), K)
    return pk_pmf+pk_extrapolate#pk_class[0]

#HMF
def HMF(zz,B1mpc,nB,eta):
    size=105
    start=-2
    if (B1mpc==0):
        end=2
    else:
        end=log(5*kd(B1mpc,nB,eta))/log(10)
    step=(end-start)/size
    K=10**np.arange(start,end,step)
    logM = np.arange(5.1,15.,0.1)
    
    HMF_pmf = C.halo_mass_function(logM = logM, k = K, pk = 1*Pm_tot(1/(1+zz),K,B1mpc,nB,eta), window = 'th', mass_fun = 'ShethTormen')

    temp=np.stack([10.**logM,10.**logM*HMF_pmf[0]], axis=0)
    
    return temp

#Extrapolate Pk with gravity from z=100
def Pks_z(zz,Pb_i,PDM_i,K):#Pbi_i and PDM_i are power specctra at z=100, zz is the redshifts at which I want the output
    PDM=np.zeros([zz.shape[0],K.shape[0]])
    Pb=np.zeros([zz.shape[0],K.shape[0]])
    Pm=np.zeros([zz.shape[0],K.shape[0]])
    j=0
    for k in K:
        delta_b0=(k**3*Pb_i[j]/(2*np.pi**2))**0.5
        delta_dm0=(k**3*PDM_i[j]/(2*np.pi**2))**0.5
        phi_0=-3/2*(0.01*Ht(0.01))**2/k**2*(rhob(0.01)*delta_b0+rhodm(0.01)*delta_dm0)/rho_tot(0.01)
        theta_b0=-0.01*Ht(0.01)*delta_b0
        theta_dm0=-0.01*Ht(0.01)*delta_dm0

        deltasolve0=[phi_0,delta_dm0,theta_dm0,delta_b0,theta_b0] #initial condition assumes everythong zero as perturbations are sourced by PMFs. I start my computation from a_mfp at which point silk damping should have erasedd any earlier evolution.
        deltasolve=solve_ivp(lambda a,y: ddel_baryon_grav(a,y,k),[0.01,1],deltasolve0,method='BDF',dense_output=True,atol=1e-6,rtol=1e-5)
        
        l=0
        for z in zz:
            if z==0:
                a_z=0.99
            else:
                a_z=1/(1+z)

            i_z=np.where(deltasolve.t>a_z)[0][0]
            PDM[l,j]=deltasolve.y[1,i_z]**2
            Pb[l,j]=deltasolve.y[3,i_z]**2
            Pm[l,j]=((0.022*deltasolve.y[3,i_z]+0.12*deltasolve.y[1,i_z])/0.142)**2
            l=l+1
        j=j+1
    return [Pb,PDM,Pm]
# %%
