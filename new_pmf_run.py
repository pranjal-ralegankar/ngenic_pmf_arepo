import numpy as np

#In[]
#Choose PMF parameters
B1mpc=0.8 #in nanoGauss
nB=-2 
eta=0.3 #0 #set eta to zero if want to use kd_i in the initial pmf power spectrum

test=0 #if set to 1 then we minimize lcdm contribution to focus only on PMFs and also turn off lorentz displacements

z_start=99# 1090 #starting redshift

boost_off=0 #set to 1 if lorentz displacement is to be turned off

#Choose box parameters
Box=20 #in Mpc/h
grid=128

PB_kd_set_kdi=1 #set to one if want to compare analytical ks set to kd_I
##############################checks######################################
if grid%64==0:
    fac=int(grid/64)
else:
    raise ValueError("grid size is not multiple of 64. My whole code requires grid size to be a multiple of 64.")

if nB==-2.9:
    boost=2.23
if nB==-2.0:
    boost=1.2
if nB==-2.5:
    boost=1.5
if boost_off==1:
    boost=0
#In[]#############################################################################################################################################################
#I read a basic parameter file and then edit it for non-vanilla changes
with open('NGENIC-PMFs_new/parameterfiles/B_1_nB_2_9_box_55_512_eta_0_3.param', 'r') as file:
    # Read the entire file
    content = file.read()

filename_no_eta=0
if eta==0:
    import sys
    from pathlib import Path

    # Add the subdir directory to sys.path
    sys.path.append(str(Path('../').resolve()))

    from paper_convention import kd_i, G, k1mpc
    eta=G(nB)/1.14*(0.1*B1mpc*(kd_i(B1mpc,nB)/k1mpc)**((nB+5)/2))**4

    filename_no_eta=1

modified_content=content.replace('B1mpc            1','B1mpc            '+str(B1mpc))
modified_content=modified_content.replace('nB    		 -2.9','nB    		 '+str(nB))
modified_content=modified_content.replace('eta    		 0.3','eta    		 '+str(eta))

modified_content=modified_content.replace('Box              55000.0','Box              '+str(Box*1000))

#filename is the name of the parameter file in which the new parameters will be saved
if filename_no_eta==1:
    filename='B_'+str(B1mpc).replace('.','_')+'_nB_'+str(-nB).replace('.','_')+'_box_'+str(Box)+'_'+str(grid)
else:
    filename='B_'+str(B1mpc).replace('.','_')+'_nB_'+str(-nB).replace('.','_')+'_box_'+str(Box)+'_'+str(grid)+'_eta_'+str(eta).replace('.','_')

modified_content=modified_content.replace('boost            2.23','boost            '+str(boost))

if test==1:
    modified_content=modified_content.replace('lambda_lcdm_b    1','lambda_lcdm_b    0.0001')
    modified_content=modified_content.replace('lambda_lcdm_dm   1','lambda_lcdm_dm    0.0001')
    modified_content=modified_content.replace('boost            2.23','boost            0')
    filename='test_'+filename


if grid!=512:
    modified_content=modified_content.replace('Nmesh           512','Nmesh           '+str(grid))
    modified_content=modified_content.replace('Nsample         512','Nsample         '+str(grid))
    modified_content=modified_content.replace('GlassTileFac      8','GlassTileFac      '+str(fac))
    if grid<256 or grid==1024:
        modified_content=modified_content.replace('NumFilesWrittenInParallel 4','NumFilesWrittenInParallel '+str(fac))

if z_start!=99:
    modified_content=modified_content.replace('Redshift         99','Redshift         '+str(z_start))
    filename=filename+'_z_'+str(z_start)

#writing the modified parameters into a file
with open('NGENIC-PMFs_new/parameterfiles/'+filename+'.param', 'w') as file:
    # Write some initial content to the file
    file.write(modified_content)

#In[]########################################################################################################################################
#Now I edit the powerspectra file
with open('hdf5_powerspectra_base.py', 'r') as file:
    # Read the entire file
    content2 = file.read()

modified_content2=content2.replace("folder_name='B_1_nB_2_9_box_7_128_eta_0_3_renorm_test4'","folder_name='"+filename+".hdf5'")
modified_content2=modified_content2.replace("cpus          = 1","cpus          = "+str(fac))

modified_content2=modified_content2.replace("B1mpc=1","B1mpc="+str(B1mpc))
modified_content2=modified_content2.replace("nB=-2.9","nB="+str(nB))
modified_content2=modified_content2.replace("eta=0.3","eta="+str(eta))
modified_content2=modified_content2.replace("a=0.01","a="+str(1/(z_start+1)))

if PB_kd_set_kdi==1:
    modified_content2=modified_content2.replace("PB_kd_set_kdi=0","PB_kd_set_kdi=1")

if test==1:
    modified_content2=modified_content2.replace('lambda_lcdm_b=1','lambda_lcdm_b=0.0001')
    modified_content2=modified_content2.replace('lambda_lcdm_dm=1','lambda_lcdm_dm=0.0001')
    modified_content2=modified_content2.replace('lambda_pmf_b=1','lambda_pmf_b=0')
    modified_content2=modified_content2.replace('lambda_pmf_dm=1','lambda_pmf_dm=0')

#writing the modified parameters into a file
with open('hdf5_powerspectra_run.py', 'w') as file:
    # Write some initial content to the file
    file.write(modified_content2)
#In[]#################################################################################################################################################
#Now I edit the bash script job file

with open('job_bfld_base.sh', 'r') as file:
    # Read the entire file
    content3 = file.read()

modified_content3=content3.replace('#SBATCH --job-name=pranjal_pylians_512','#SBATCH --job-name='+filename)

if grid==128:
    modified_content3=modified_content3.replace('mpirun -np 2 ./N-GenIC parameterfiles/temp.param','mpirun -np 2 ./N-GenIC parameterfiles/'+filename+'.param')
elif grid==256:
    modified_content3=modified_content3.replace('#SBATCH --ntasks-per-node=2','#SBATCH --ntasks-per-node=4')
    modified_content3=modified_content3.replace('mpirun -np 2 ./N-GenIC parameterfiles/temp.param','mpirun -np '+str(fac)+' ./N-GenIC parameterfiles/'+filename+'.param')
elif grid==512:
    modified_content3=modified_content3.replace('#SBATCH --mem=10000mb','#SBATCH --mem=30000mb')
    modified_content3=modified_content3.replace('#SBATCH --ntasks-per-node=2','#SBATCH --ntasks-per-node=8')
    modified_content3=modified_content3.replace('#SBATCH --time=00:5:00','#SBATCH --time=00:10:00')
    modified_content3=modified_content3.replace('mpirun -np 2 ./N-GenIC parameterfiles/temp.param','mpirun -np 8 ./N-GenIC parameterfiles/'+filename+'.param')
elif grid==1024:
    modified_content3=modified_content3.replace('#SBATCH -N 1','#SBATCH -N 2')
    modified_content3=modified_content3.replace('#SBATCH --ntasks-per-node=2','#SBATCH --ntasks-per-node=8')
    modified_content3=modified_content3.replace('#SBATCH --mem=10000mb','#SBATCH --mem=120000mb')
    modified_content3=modified_content3.replace('#SBATCH --time=00:5:00','#SBATCH --time=00:40:00')
    modified_content3=modified_content3.replace('mpirun -np 2 ./N-GenIC parameterfiles/temp.param','mpirun -np '+str(fac)+' ./N-GenIC parameterfiles/'+filename+'.param')
else:
    raise ValueError("For chosen grid size need to add custom bash script parameters")

modified_content3=modified_content3.replace('python3 hdf5_output.py temp','python3 hdf5_output.py '+filename)

with open('job_bfld_run.sh', 'w') as file:
    # Write some initial content to the file
    file.write(modified_content3)

#In[]############################################################################################################################################3
#submit the script to run
import subprocess
subprocess.run(['sbatch', 'job_bfld_run.sh'])
# %%
