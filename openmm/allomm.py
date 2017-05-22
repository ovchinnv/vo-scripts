#!/bin/python
#requires the file CHOMM.py, which is simple wrapper function to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
firstrun=1   ;# initial run index
numrun=2
name='traf-lig27'
#platformName='CPU' ; #optional; default os 'CUDA'
#==============================
# parameters required by CHOMM (some have default values)

psffile='../struc/'+name+'_sn.psf' ;
pdbfile='../struc/'+name+'-msn.pdb' ;
topfile='../struc/'+name+'.top';
paramfile='../struc/'+name+'.par' ;

implicitSolvent=0 ;# run OBC2 implicit solvent simulation

xscfile='.xsc' ;         # to obtain cell vectors from last line of xsc file
boxfile='../struc/'+name+'.str'; # to obtain cell vectors from str file used in structure solvation

# specify larger box manually
#dx=112 ;
#dy=112 ;
#dz=112 ;

flag='eq';
hmass=1;       # amu, can use heavy hydrogens
friction=0.1   # 1/ps, thermostat coupling
dt=1;          # timestep in fs
pmefreq=1;     # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=9;      # nonbonded cutoff
switchdist=7.5 ; # (optional) switching distance

constraints=1;   # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=pdbfile; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology
conscol=1 ; # 1 beta ; 2 occupancy

shake=0; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=298; # kelvin
andersen=0;    # to use Andersen instead of Langevin ; (Note that I see energy up drifts quite often with andersen)
barostat=1;
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=1; # whether to use PME
pbc=1; # whether periodic boundary conditions are on

mini=1;          # whether to minimize before dynamics
ministeps=100;   # number of minimization iterations

numeq=           # number of equilibration runs
numeqsteps=100000; # number of equilibration steps
nummdsteps=500000; # number of production steps
outputfreq=10000;  # frequency of generating output
dcdfreq=10000;     # frequency of dcd output

flag='eq'
nsteps=numeqsteps;
#
if (firstrun==0):
 restart=0 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
 restartfile=None ;
else:
 restart=1 ;
 restartfile='scratch/'+name+str(firstrun-1)+flag+'.xml';
#restartfile= ;# to override
 xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
#
# run MD simulations with different parameters one after the other
#
irun=firstrun
while irun < numrun:
 print(" =============================");
 print(" Run ", irun, "(will quit after", numrun-1,")");
# set some run-specific options
# constraintscaling = (90-10*irun) ;# turn off gradually by run 10
 if irun >= numeq:
   flag='nvt'
   constraints=0 ;# to remove equilibration restraints
   nsteps=nummdsteps ;# increase number of steps
   hmass=4.0
   dt=4.0
   mini=0
   shake=1
   barostat=0 ;# turn off barostat
   pmefreq=1 ; # using MTS does not improve speed in my experience
#   cutoff=9 ;# to change cutoff
#
 outputName='scratch/'+name+str(irun)+flag ;
 from os.path import expanduser
 exec(open(expanduser('~/scripts/openmm/CHOMM.py')).read())
# to do : check if run was successful, rerun if not
 i=i+1
 restartfile=outputName+'.xml'
 xmlfile=restartfile
 restart=1
