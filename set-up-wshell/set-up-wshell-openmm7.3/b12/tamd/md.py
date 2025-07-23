#!/usr/bin/python
from __future__ import print_function
import fileinput
import time
from os import path, mkdir
#requires the file CHOMM.py, which is simple wrapper function to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
firstrun=1   ;# initial run index (initialization at firstrun=1, unlike plain MD)
endrun=2
name='b12'
strucdir='../struc';
#platformName='CPU' ; #optional; default is 'CUDA'
#==============================
if not (path.exists('scratch')):
 mkdir('scratch')
# parameters required by CHOMM (some have default values)

psffile=strucdir+'/'+name+'_sn.psf' ;
pdbfile=strucdir+'/'+name+'_sn.pdb' ;
topfile=strucdir+'/'+name+'.top';
paramfile=strucdir+'/'+name+'.par' ;

implicitSolvent=0 ;# run OBC2 implicit solvent simulation

hmass=4;       # amu, can use heavy hydrogens
friction=1   # 1/ps, thermostat coupling
dt=4;          # timestep in fs
pmefreq=1;     # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=9;      # nonbonded cutoff
switchdist=8.5 ; # (optional) switching distance

constraints=1;   # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=strucdir+'/'+name+'_sn-res.pdb'; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology
conscol=1 ; # 1 beta ; 2 occupancy

shake=1; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=298; # kelvin
andersen=0;    # to use Andersen instead of Langevin ; (Note that I see energy up drifts quite often with andersen)
barostat=0;
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=0; # whether to use PME
pbc=0; # whether periodic boundary conditions are on
removeCOM=0 ; # whether to remove COM drift

dynamo=1
dynamoTemplate='dynamo-init.tmp' ;# initial input ; changed for later runs
watershell_restart='NONE'
tamdtemp=3000; # TAMD temperature (see Maragliano & EVE 2006 publication for discussion)

mini=1;          # whether to minimize before dynamics
ministeps=0;   # number of minimization iterations

numeq=1             # number of equilibration runs
numeqsteps=100000; # number of equilibration steps
nummdsteps=100000; # number of production steps
#nummdsteps=20000
outputfreq=10000;  # frequency of generating output
dcdfreq=10000;     # frequency of dcd output

nsteps=numeqsteps;
#
flag='eq'
if (firstrun==1):
 restart=0 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
 restartfile=None ;
else:
 if (firstrun>2): # run 1 restarts from 0eq, others, from <i>nvt
  flag='nvt';
 restart=1 ;
 restartfile='scratch/'+name+str(firstrun-1)+flag+'.xml';
#restartfile= ;# to override
 xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
#
# run MD simulations with different parameters one after the other
#
irun=firstrun
while irun <= endrun:

 print(" =============================");
 print(" Run ", irun, "(will quit after", endrun,")");
# set some run-specific options
# constraintscaling = (90-10*irun) ;# turn off gradually by run 10
 if irun > numeq:
   flag='nvt'
   constraints=0 ;# to remove equilibration restraints
   nsteps=nummdsteps ;# increase number of steps
   hmass=4.0
   dt=4.0
   mini=0
   shake=1
   barostat=0 ;# turn off barostat
   pmefreq=1 ; # in case PME is on, which it might not be (see above)
   dynamoTemplate='dynamo.tmp'; # main plugin file
#   cutoff=9 ;# to change cutoff
#
# dynamo section :
#
 watershell_output='watershell'+str(irun)+'.restart.txt'
 if (irun>1):
   watershell_restart='watershell'+str(irun-1)+'.restart.txt'
#
 dynamoConfig='dynamo'+str(irun)+'.in'
# modify config template :
 df=open(dynamoTemplate,'r');
 dd=df.read()
 dd=dd.replace('@{restart_file}',watershell_restart)
 dd=dd.replace('@{output_file}',watershell_output)
 dd=dd.replace('@{iprev}',str(irun-1))
 dd=dd.replace('@{irun}',str(irun))
 dd=dd.replace('@{tamdtemp}',str(tamdtemp))
 df=open(dynamoConfig,'w');
 df.write(dd);
 df.close();
 dynamoLog=dynamoConfig+'.log';
 outputName='scratch/'+name+str(irun)+flag ;
 from os.path import expanduser
 exec(open(expanduser('../../CHOMM.py')).read())
# to do : check if run was successful, rerun if not
 irun=irun+1
 restartfile=outputName+'.xml'
 xmlfile=restartfile
 restart=1
 time.sleep(10)
 