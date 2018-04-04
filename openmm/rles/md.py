#!/bin/python2
from __future__ import print_function
#requires the file CHOMM.py, which is simple wrapper function to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
firstrun=1   ;# initial run index
numrun=1    ;# numbef of runs
name='hp36' ; #prefix for output files
#platformName='OpenCL' ; #optional; default is 'CUDA'
#==============================
# parameters required by CHOMM (some have default values)
restart=1 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
restartfile=None ;
restartfile='hp360eq.xml';

psffile='struc/'+name+'_sno3nbx.psf' ;
pdbfile='struc/'+name+'_msn3.pdb' ;
topfile='struc/rlestop.top';
paramfile='struc/rlespar1o3.par' ;

implicitSolvent=0 ;# run OBC2 implicit solvent simulation

xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
xscfile='.xsc' ;         # to obtain cell vectors from last line of xsc file
boxfile='struc/'+name+'.str'; # to obtain cell vectors from str file used in structure solvation

# specify larger box manually
dx=62 ;
dy=62 ;
dz=62 ;

hmass=1;       # amu, can use heavy hydrogens
friction=1     # 1/ps, thermostat coupling
dt=2;        # timestep in fs
pmefreq=1;     # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=9;     # nonbonded cutoff

constraints=0;   # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=pdbfile; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology
conscol=1 ; # 1 beta ; 2 occupancy

shake=1; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=300; # kelvin
andersen=0;    # to use Andersen instead of Langevin ; (Note that I see energy up drifts quite often with andersen)
barostat=1;
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=1; # whether to use PME
pbc=1; # whether periodic boundary conditions are on

# for RLES:
alchfile=pdbfile;
alchcol=1;
rles_nrep=3;

mini=1;          # whether to minimize before dynamics
ministeps=1000;   # number of minimization iterations

nsteps=100000;    # number of simulation steps
outputfreq=1000; # frequency of generating output
dcdfreq=10000;   # frequency of dcd output

flag='eq'
outputName=name+str(firstrun)+flag ;
#
# run MD simulations with different parameters one after the other
#
for i in range(numrun):
 irun=i+firstrun;
 print(" =============================");
 print(" Run ", irun, "(will quit after", firstrun+numrun-1,")");
# set some run-specific options
# constraintscaling = (90-10*irun) ;# turn off gradually by run 10
 if irun > 0:
   flag='nvt'
   constraints=0 ;# to remove equilibration restraints
   nsteps=10000000 ;# increase number of steps
#   barostat=0;
   mini=0;
#   hmass=1.0
#   dt=1.0
#   cutoff=9 ;# to decrease cutoff
#
 outputName=name+str(irun)+flag ;
#
 from os.path import expanduser
 exec(open(expanduser('CHOMMrles.py')).read())
 restart=1
 restartfile=outputName+'.xml';
# check if run was successful, rerun if not
