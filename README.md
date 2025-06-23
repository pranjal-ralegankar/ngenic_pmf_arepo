# ngenic_pmf_arepo
A modified NGenIC code that generates initial conditions for AREPO with primordial magnetic fields
This code should be used to generate ics at z=1099, i.e. directly after recombination. 

Executing this code requires the installation of FFTW module https://www.fftw.org/, pylians https://pylians3.readthedocs.io/en/master/installation.html, Colibri https://github.com/GabrieleParimbelli/COLIBRI, and CLASS https://github.com/lesgourg/class_public.

The NGENIC-PMFs_new folder contains all the C files that actually generate the ics. 
In the main.c code, lines 131 to 155 and 321-323 inputs the analytical magnetic field power spectrum to Ngenic. So if you want to change the form of the input power spectrum you need to play around that code.

The code outside the NGENIC-PMFs_new folder is intended to automate running Ngenic and is based on python. This code generates a parameter file, then writes a job script to run on computing cluster, and then prints the ngenic output for you to check if everything is running as expected. 

The script is "new_pmf_run.py". In the first 20 lines you need to specify simulation parameters like pmf strength box size etc. To change anything outside of those parameters, you need to edit the B_1_nB_2_9_box_55_512_eta_0_3.param file in \NGENIC-PMFs_new\parameterfiles.
The basic idea is that for different simulations, we only need to change one or two parameters while remaining code remains the same. So the folder "base_code" provide the basic skeleton of codes, which are edited by "new_pmf_run.py" to cater to different simulations.
From lines 44-104,  "new_pmf_run.py" takes a basic .param file in Ngenic folder and edits it according to the simulation specifications provided in the first 20 lines.
From lines 108 to 124, it edits the "powerspecta_run.py" file. This python file is in charge of analyzing the output of Ngenic and comparing it with expectations.
From lines 126 to 160, it edits the job script file "job_bfld_run.sh". It specifies running Ngenic with correct parameter file and specifies machine specification.
Finally it submits the job to the cluster.

Note that the details of job script have to be changed according to your specific machine. To do that you must edit the "job_bfld_base.sh" file in "base_code" folder. The "job_bfld_base.sh" file acts like the basic skeleton of the job script which is later edited by "new_pmf_run.py" for different pmf configs. 

Finally, let me explain the paper_convention.py file in the base_code folder. This file contains all the analytical physics from the paper https://arxiv.org/abs/2410.02676. The functions PB, Pb_an, Pdm_an, etc generates power spectrum at specified redshift and time according to the results of the paper  https://arxiv.org/abs/2410.02676. The functions Pb_an, Pdm_an etc are imported from this code in "powerspectra_run.py" to compare analytical results with Ngenic output. 

To generate initial conditions, first go to new_pmf_run.py and edit the first 20 lines to specify all the desire parameters. Running this code will submit a job to your machine. Once the job_script is executed you would find the initial conditions in the folder "ics". Additionally, in the folder pylians_output you will find the power spectrum of the generated initial condition compared with the analytical expectation.
