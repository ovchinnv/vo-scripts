#!/bin/python
#requires the file CHOMM.py, which is simple wrapper function to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
firstrun=1   ;# initial run index
numrun=19    ;# numbef of runs
name='fre' ; #prefix for output files

#==============================
# parameters required by CHOMM (some have default values)
restart=1 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
restartfile=None ;
restartfile='fre0eq.xml';

psffile='struc/fre-fad-x.psf' ;
pdbfile='struc/fre_fad.pdb' ;
paramfile='param/fre.prm' ;

xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
xscfile='fre0.xsc' ;         # to obtain cell vectors from last line of xsc file
boxfile='struc/boxsize.str'; # to obtain cell vectors from str file used in structure solvation

hmass=1;         # amu, can use heavy hydrogens
friction=1     # 1/ps, thermostat coupling
dt=1.0;            # timestep in fs
pmefreq=1;       # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=12;      # nonbonded cutoff

constraints=1;   # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=pdbfile; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology

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

nsteps=20000;    # number of simulation steps
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
 constraintscaling = (90-10*irun) ;# turn off gradually by run 10
 if irun > 0:
  barostat=0
  mini=0
  shake=1
  dt=2
  if irun > 9:
   flag='nvt'
   constraints=0
   nsteps=10000000 ;# increase number of steps
   hmass=4
   dt=3.5
   cutoff=9
#
 outputName=name+str(irun)+flag ;
#
 exec(open('CHOMM.py').read())
 restart=1
 restartfile=outputName+'.xml';
# check if run was successful, rerun if not


