#!/bin/python
from __future__ import print_function
import fileinput
import sys
from os import path
#=====================================================================================
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
firstrun=5   ;# initial run index
numrun=5    ;# numbef of runs
name='hp35c27' ; #prefix for output files
rles_nrep=4 ; # number of rles replicas
#platformName='OpenCL' ; #optional; default is 'CUDA'
#==============================
psffile='struc/'+name+'_sno'+str(rles_nrep)+'nbx.psf' ;
pdbfile='struc/'+name+'_msn'+str(rles_nrep)+'.pdb' ;
topfile='struc/rlestopc27.top';
paramfile='struc/rlespar1o'+str(rles_nrep)+'c27.par' ;

implicitSolvent=0 ;# run OBC2 implicit solvent simulation

#xscfile='.xsc' ;         # to obtain cell vectors from last line of xsc file
boxfile='struc/'+name+'.str'; # to obtain cell vectors from str file used in structure solvation

# specify larger box manually
#dx=62 ;
#dy=62 ;
#dz=62 ;

hmass=1;       # amu, can use heavy hydrogens
friction=0.1     # 1/ps, thermostat coupling
dt=1;        # timestep in fs
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
barostat=0;
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=1; # whether to use PME
pbc=1; # whether periodic boundary conditions are on

# for RLES:
alchfile=pdbfile;
alchcol=1;

dynamo=1
dynamoTemplate='rles.dyn'
adaptive_restart='NONE'

mini=1;          # whether to minimize before dynamics
ministeps=1000;   # number of minimization iterations

nsteps=10000000;    # number of simulation steps
outputfreq=10000; # frequency of generating output
dcdfreq=10000;   # frequency of dcd output

flag='nvt'
outputName=name+str(firstrun)+flag ;
#
if (firstrun==0):
 restart=0 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
 restartfile=None ;
else:
 restart=1 ;
 if (firstrun>0):
  flag='nvt'
 restartfile='./'+name+str(firstrun-1)+flag+'.xml';
#restartfile= ;# to override
 xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
#
# run MD simulations with different parameters one after the other
#
irun=firstrun
while irun < firstrun + numrun :
 print(" =============================");
 print(" Run ", irun, "(will quit after", firstrun+numrun-1,")");
# set some run-specific options
 if irun > 0:
#   flag='nvt'
#   constraints=0 ;# to remove equilibration restraints
#   nsteps=30000000 ;# increase number of steps
#   barostat=0;
   mini=0;
#   hmass=1.0
#   dt=1.0
#   cutoff=9 ;# to decrease cutoff
#
# dynamo section :
#
 adaptive_output='adaptive.restart'+str(irun)+'.txt'
 if (irun>0):
   adaptive_restart='adaptive.restart'+str(irun-1)+'.txt'
#
 dynamoConfig='rles'+str(irun)+'.in'
# modify config template :
 df=open(dynamoTemplate,'r');
 dd=df.read()
 dd=dd.replace('@{restart_file}',adaptive_restart)
 dd=dd.replace('@{output_file}',adaptive_output)
 dd=dd.replace('@{irun}',str(irun))
 dd=dd.replace('@{iprev}',str(irun-1))
 df=open(dynamoConfig,'w');
 df.write(dd);
 df.close();
 dynamoLog=dynamoConfig+'.log';

 outputName='./'+name+str(irun)+flag ;
 from os.path import expanduser
 exec(open(expanduser('~/scripts/openmm/rles/CHOMMrles.py')).read())
 irun=irun+1
 restartfile=outputName+'.xml'
 xmlfile=restartfile
 restart=1

