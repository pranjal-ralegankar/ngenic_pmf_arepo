
#In[]
import ics_to_hdf5 as ith
import os

################# read ics$###########################
# folder_name='B_1_nB_2_9_box_55_512_eta_0_3_renorm'
folder_name='B_1_nB_2_9_box_7_128_eta_0_3_renorm_test4'
folder_sims='ics'
fullpath_snap = os.path.join(folder_sims, folder_name, 'ics')

gfile = ith.ReadGadgetMultipleFiles(fullpath_snap, fformat=0, longids=False, discard=['u'], opt_fields=['u', 'bfld'],  stop_at='', quiet=False, fill_mass=False, fixVelocity=False, verbose=False, swapbytes=False, max_num_files=0)

# gfile = ith.ReadGadget(fullpath_snap, fformat=0, longids=False, discard=['u'], opt_fields=['u', 'bfld'],  stop_at='', quiet=False, fill_mass=False, fixVelocity=False, verbose=False, swapbytes=False)
#In[]
############### obtain Bfields$$$######################
import numpy as np

Npart=gfile.npart[0]
BoxSize=gfile.BoxSize[0]/1e3 #In Mpc/h
grid=round(Npart**(1/3))

B=gfile.bfld[:]
pos=gfile.pos[:Npart]/1e3 #In Mpc/h
pos_dm=gfile.pos[Npart:]/1e3 #In Mpc/h

Bx = np.zeros((grid,grid,grid), dtype=np.float32)
By = np.zeros((grid,grid,grid), dtype=np.float32)
Bz = np.zeros((grid,grid,grid), dtype=np.float32)

n_b = np.zeros((grid,grid,grid), dtype=np.float32)
n_dm=np.zeros((grid,grid,grid), dtype=np.float32)
#In[]
############Finding density field######################
import MAS_library as MASL
MAS     = 'CIC'
verbose = True 

MASL.MA(pos, n_b, BoxSize, MAS, verbose=verbose)
delta_b=n_b/np.mean(n_b, dtype=np.float64)-1.0

MASL.MA(pos_dm, n_dm, BoxSize, MAS, verbose=verbose)
delta_dm=n_dm/np.mean(n_dm, dtype=np.float64)-1.0

rho_m=n_b*gfile.massarr[0]+n_dm*gfile.massarr[1]
delta_m=rho_m/np.mean(rho_m, dtype=np.float64)-1.0
#In[]
#############Finding weighted Bx,By,Bz fields##########
weight_Bx=B[:,0]
MASL.MA(pos_dm, Bx, BoxSize, MAS, W=weight_Bx, verbose=verbose)

weight_By=B[:,1]
MASL.MA(pos_dm, By, BoxSize, MAS, W=weight_By, verbose=verbose)

weight_Bz=B[:,2]
MASL.MA(pos_dm, Bz, BoxSize, MAS, W=weight_Bz, verbose=verbose)

#In[]
############ Pure B field##################
Bx_field = Bx/n_dm
By_field = By/n_dm
Bz_field = Bz/n_dm

print("Bfield computed")
#In[]
# ########### V_i power spectrum################
import Pk_library as PKL

axis          = 0      #no RSD                                                                                                                                           
cpus          = 1      #number of openmp threads
Pk_Bx = PKL.Pk(Bx_field, BoxSize, axis, MAS, cpus, verbose)
Pk_By = PKL.Pk(By_field, BoxSize, axis, MAS, cpus, verbose)
Pk_Bz = PKL.Pk(Bz_field, BoxSize, axis, MAS, cpus, verbose)

Pk_b_temp = PKL.Pk(delta_b, BoxSize, axis, MAS, cpus, verbose)
Pk_dm_temp = PKL.Pk(delta_dm, BoxSize, axis, MAS, cpus, verbose)
Pk_m_temp = PKL.Pk(delta_m, BoxSize, axis, MAS, cpus, verbose)

#In[]
# ######Obtaining total velocity power spectrum############
Pk_B=Pk_Bx.Pk[:,0]+Pk_By.Pk[:,0]+Pk_Bz.Pk[:,0]
K_B=Pk_Bx.k3D

Pk_b=Pk_b_temp.Pk[:,0]
K_b=Pk_b_temp.k3D

Pk_dm=Pk_dm_temp.Pk[:,0]
K_dm=Pk_dm_temp.k3D

print("power spectra computed")
#In[]
############ Theoretical estimation#################
import sys
from pathlib import Path

current_dir = os.path.abspath(os.getcwd()) #get current directory
# Add base_code directory to sys.path
sys.path.append(os.path.join(current_dir, 'base_code'))

from paper_convention import table, PB_an, Pb_an, Pdm_an, Pm_an, kd, deltasolve, Pm_tot


B1mpc=1
nB=-2.9
eta=0.3
a=0.01

lambda_lcdm_b=1
lambda_lcdm_dm=1
lambda_pmf_b=1
lambda_pmf_dm=1

PB_an_table=table(lambda k: PB_an(k,B1mpc,nB,eta),K_B)

#I approximate baryond and dark matter LCDM power spectrum to be the same as matter power spectrum for simplicity
Pdm_an_table=lambda_lcdm_dm**2*table(lambda k: Pm_tot(a,k,0,nB,eta),K_dm)+lambda_pmf_dm**2*table(lambda k: Pdm_an(a,k,B1mpc,nB,eta),K_dm)

Pb_an_table=lambda_lcdm_b**2*table(lambda k: Pm_tot(a,k,0,nB,eta),K_b)+lambda_pmf_b**2*table(lambda k: Pb_an(a,k,B1mpc,nB,eta),K_b)

#In[]
############# Plotting###########################
print("plotting power spectra")
folder_out='Pk_plot/'+folder_name
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2,2,figsize=(14,12))
plt.subplots_adjust(left=0.06,right=0.96,top=0.96,bottom=0.09,wspace=0.2,hspace=0.05)

ax[0,0].loglog(K_B,Pk_B,'-k')
ax[0,0].loglog(K_B,PB_an_table,'-b')

ax[0,1].loglog(K_B/kd(B1mpc,nB,0.3),Pk_B/PB_an_table,'-k')

ax[1,0].loglog(K_dm,K_dm**2/2/np.pi**2*Pk_dm,'-', linewidth=2)
ax[1,0].loglog(K_dm,K_dm**2/2/np.pi**2*Pdm_an_table,'--', linewidth=2)

ax[1,1].loglog(K_b,K_b**2/2/np.pi**2*Pk_b,'-', linewidth=2)
ax[1,1].loglog(K_b,K_b**2/2/np.pi**2*Pb_an_table,'--', linewidth=2)

fig.savefig(folder_out+".png", format="png", bbox_inches="tight")
# %%
