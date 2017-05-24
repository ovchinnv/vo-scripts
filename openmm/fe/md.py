#!/bin/python
#requires the file CHOMM.py, which is a simple wrapper to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#platformName='CPU' ; # default is CUDA; CPU is usually too slow ; OpenCL is a bit slower than CUDA
#platformName='OpenCL' ; # default is CUDA; CPU is usually too slow
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
firstrun=64   ;# initial run index
numrun=1    ;# numbef of runs
name='eth-k' ; #prefix for output files

#==============================
# parameters required by CHOMM (some have default values)
restart=1 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
#restartfile=None ;
restartfile='eth-k0nvt.xml';

psffile='struc/eth-k-wat.psf' ;
pdbfile='struc/eth-k-wat.pdb' ;
topfile='param/eth-k.top';
paramfile='param/eth-k.par' ;

xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
#xscfile=None ;         # to obtain cell vectors from last line of xsc file

dx=28.; # specify box size manually
dy=28.;
dz=28.;

hmass=1.       # amu, can use heavy hydrogens
friction=1     # 1/ps, thermostat coupling
dt=2;          # timestep in fs
pmefreq=1;     # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=10;     # nonbonded cutoff

constraints=1;   # harmonic positional restraints for equilibration
constraintscaling=0.25; # to scale hatmonic restraints uniformly
consfile=pdbfile; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology
conscol=2;        # 1--beta ; 2--occupancy

shake=1; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=300; # kelvin
andersen=0;    # to use Andersen instead of Langevin ; (Note that I see energy up drifts quite often with andersen)
barostat=1;
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=1; # whether to use PME
pbc=1; # whether periodic boundary conditions are on

mini=1;          # whether to minimize before dynamics
ministeps=1000;  # number of minimization iterations

#=== alchemical section
alch=0 ;         # alchemical transformation on or off
alchfile=pdbfile;# PDB file with atoms marked
alchcol=1;       # column : 1 -- beta ; 2 -- occupancy
alchout='eth-k'; # base file name for alchemical output
alchfreq=100;    # frequency of alchemical output
alchdecouple=1 ; # whether to decouple or to annihilate

nsteps=10000;    # number of simulation steps
outputfreq=5000;  # frequency of generating output
dcdfreq=1000;    # frequency of dcd output

flag='nvt'
outputName=name+str(firstrun)+flag ;
#
# run MD simulations with different parameters one after the other
#
for i in range(numrun):
 irun=i+firstrun;
 print(" ==========================================================");
 print(" Run ", irun, "(will quit after", firstrun+numrun-1,")");
# set some run-specific options
# constraintscaling = (90-10*irun) ;# turn off gradually by run 10
 switchdist=cutoff-1 ;
 if ( irun > 0 ) :
  mini=0 ;
  nsteps=10000 ;
  alch=1;
  nalchwin=32 ;
  dwin=1/nalchwin;
  if ( irun > 2*nalchwin ) :
   print(" ==========================================================");
   print(" Quitting before run", irun, "because FEP cycle is complete");
   print(" ==========================================================");
   break ; # completed FEP cycle
  if ( irun > nalchwin ) : # backward
   lambda0 = (irun-nalchwin-1) * dwin ;
   lambda1 = lambda0 + dwin ;
  else : # forward
   lambda0 =-(irun-nalchwin-1) * dwin ;
   lambda1 = lambda0 - dwin ;
#
#  print(lambda0)
#  print(lambda1)
#
 outputName=name+str(irun)+flag ;
#
 exec(open('CHOMM.py').read())
 restart=1
 restartfile=outputName+'.xml';
 xmlfile=restartfile;
# to do : check if run was successful, rerun if not


