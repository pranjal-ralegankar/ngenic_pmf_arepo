
#In[]
import ics_to_hdf5 as ith
import os
import numpy as np

################# read snap$###########################
folder_name='B_1_nB_2_9_box_7_128_eta_0_3_renorm_test4'
folder_sims='hdf5/'
fullpath_snap1 = folder_sims+folder_name
# fullpath_snap5 = '/scratch/egaraldi/BSF/pmf_B_test4_128/output/snapdir_013/snap_013.hdf5'#'/scratch/egaraldi/BSF/pmf_B_test/output/snapdir_010/snap_010.hdf5'

gfile = ith.ReadGadget(fullpath_snap1, fformat=3, opt_fields=['bfld'],quiet=True)

UnitLength_in_cm=3.085678e21   # defines length unit of output (in cm/h) 
UnitMass_in_g=1.989e43      # defines mass unit of output (in g/cm)
UnitVelocity_in_cm_per_s=1e5           # defines velocity unit of output (in cm/sec)
Gauss_to_internal = np.sqrt(UnitLength_in_cm**3)/(np.sqrt(UnitMass_in_g) * UnitVelocity_in_cm_per_s)

############### obtain Bfields ######################

z=gfile.redshift
print("\n redshift=",z)
Npart=gfile.npart[0]
BoxSize=gfile.BoxSize/1e3 #In Mpc/h
grid=round(Npart**(1/3))

B=np.array(gfile.bfld[:], dtype=np.float32)/Gauss_to_internal*0.678
pos=np.array(gfile.pos[:Npart], dtype=np.float32)/1e3 #In Mpc/h
pos_dm=np.array(gfile.pos[Npart:], dtype=np.float32)/1e3 #In Mpc/h

Bx = np.zeros((grid,grid,grid), dtype=np.float32)
By = np.zeros((grid,grid,grid), dtype=np.float32)
Bz = np.zeros((grid,grid,grid), dtype=np.float32)
n_b = np.zeros((grid,grid,grid), dtype=np.float32)
rho_b = np.zeros((grid,grid,grid), dtype=np.float32)
rho_dm=np.zeros((grid,grid,grid), dtype=np.float32)

############Finding density field######################
import MAS_library as MASL
MAS     = 'CIC'
verbose = False 

MASL.MA(pos, n_b, BoxSize, MAS, verbose=verbose)
delta_b=n_b/np.mean(n_b, dtype=np.float64)-1.0

MASL.MA(pos_dm, rho_dm, BoxSize, MAS, verbose=verbose)
delta_dm=rho_dm/np.mean(rho_dm, dtype=np.float64)-1.0

#############Finding weighted Bx,By,Bz fields##########
zero_indices = np.where(n_b == 0) #find index values where rho_b is 0
n_b[zero_indices]=0.001 #set rho_b to a small value instead on those locations

weight_Bx=B[:,0]
MASL.MA(pos, Bx, BoxSize, MAS, W=weight_Bx, verbose=verbose)

weight_By=B[:,1]
MASL.MA(pos, By, BoxSize, MAS, W=weight_By, verbose=verbose)

weight_Bz=B[:,2]
MASL.MA(pos, Bz, BoxSize, MAS, W=weight_Bz, verbose=verbose)

############ Pure B field##################
Bx_field = Bx/n_b
By_field = By/n_b
Bz_field = Bz/n_b

print("Bfield computed")
# ########### B_i power spectrum################
import Pk_library as PKL

axis          = 0      #no RSD                                                                                                                                          
cpus          = 1      #number of openmp threads
Pk_Bx = PKL.Pk(Bx_field, BoxSize, axis, MAS, cpus, verbose)
Pk_By = PKL.Pk(By_field, BoxSize, axis, MAS, cpus, verbose)
Pk_Bz = PKL.Pk(Bz_field, BoxSize, axis, MAS, cpus, verbose)

Pk_b_temp = PKL.Pk(delta_b, BoxSize, axis, MAS, cpus, verbose)
Pk_dm_temp = PKL.Pk(delta_dm, BoxSize, axis, MAS, cpus, verbose)

print("power spectra computed")
# ######Obtaining total power spectrum############
Pk_B=Pk_Bx.Pk[:,0]+Pk_By.Pk[:,0]+Pk_Bz.Pk[:,0]
K_B=Pk_Bx.k3D

Pk_b=Pk_b_temp.Pk[:,0]
K_b_temp=Pk_b_temp.k3D
K_b=Pk_b_temp.k3D

Pk_dm=Pk_dm_temp.Pk[:,0]
K_dm_temp=Pk_dm_temp.k3D
K_dm=Pk_dm_temp.k3D

#### power spectrum from ngenic at z=100 ##########
#In[]
############ Theoretical estimation#################
import sys
from pathlib import Path

current_dir = os.path.abspath(os.getcwd()) #get current directory
# Add base_code directory to sys.path
sys.path.append(os.path.join(current_dir, 'base_code'))

from paper_convention import table, PB_an, Pb_an, Pdm_an, Pm_an, kd, deltasolve, Pm_tot, kd_i


B1mpc=1
nB=-2.9
eta=0.3
a=0.01

lambda_lcdm_b=1
lambda_lcdm_dm=1
lambda_pmf_b=1
lambda_pmf_dm=1
PB_kd_set_kdi=0

if PB_kd_set_kdi==1:
    PB_an_table=table(lambda k: PB_an(k,B1mpc,nB,10**7)*np.exp(-k**2/kd_i(B1mpc,nB)**2),K_B)
else:
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

ax[0,0].loglog(K_B,K_B**3/2/np.pi**2*Pk_B,'-k')
ax[0,0].loglog(K_B,K_B**3/2/np.pi**2*PB_an_table,'-b')
ax[0,0].set_ylim([
    (K_B**3/2/np.pi**2*PB_an_table)[0]*10**-5,
    (K_B**3/2/np.pi**2*PB_an_table).max()*100
])
# Create a mask for non-zero values in PB_an_table
mask = PB_an_table/PB_an_table[0] > 10**-30

# Use the mask to filter your data
ax[0,1].loglog(K_B[mask]/kd(B1mpc,nB,0.3),
    Pk_B[mask]/PB_an_table[mask],
    '-k')

ax[1,0].loglog(K_dm,K_dm**3/2/np.pi**2*Pk_dm,'-', linewidth=2)
ax[1,0].loglog(K_dm,K_dm**3/2/np.pi**2*Pdm_an_table,'--', linewidth=2)

ax[1,1].loglog(K_b,K_b**3/2/np.pi**2*Pk_b,'-', linewidth=2)
ax[1,1].loglog(K_b,K_b**3/2/np.pi**2*Pb_an_table,'--', linewidth=2)

fig.savefig(folder_out+".png", format="png", bbox_inches="tight")


# %%
