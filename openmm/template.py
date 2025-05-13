#!/bin/python
#requires the file CHOMM.py, which is a simple wrapper to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#platformName='CPU' ; # default is CUDA; CPU is usually too slow
#
#==============================
# parameters required by CHOMM (some have default values)
restart=1 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
restartfile='@{restartfile}' ;

psffile='@{psffile}' ;
pdbfile='@{pdbfile}' ;
paramfile='@{prmfile}';
topfile='@{topfile}' ;

implicitSolvent=1 ;# run OBC2 implicit solvent simulation

xmlfile=restartfile ;        # to obtain cell vectors from xml file produced with OMM (default option if restart file is provided)
xscfile='.xsc' ;         # to obtain cell vectors from last line of xsc file

dx=1; # specify box size manually
dy=1;
dz=1;

hmass=@{hmass};         # amu, can use heavy hydrogens
friction=@{friction}     # 1/ps, thermostat coupling
dt=@{dt};              # timestep in fs
pmefreq=@{pmefreq};       # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=@{cutoff};       # nonbonded cutoff
#switchdist=@{switchdist};   # optional : default is cutoff - 1.5

constraints=@{restraints};   # harmonic positional restraints for equilibration
constraintscaling=@{restraintscale}; # to scale harmonic restraints uniformly
consfile='@{restraintfile}'; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology
conscol=@{restraintcol}; # 1 -- beta , 2 -- occupancy

shake=@{shake}; # whether to constrain bonds involving hydrogens

thermostat=@{thermostat};  # whether to use a thermostat
temperature=@{temperature}; # kelvin
andersen=0;    # to use Andersen instead of Langevin ; (Note that I see energy up drifts quite often with andersen)
barostat=@{barostat};
pressure=1;    # units of atm
membrane_on=0; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=0; # whether to use PME
pbc=0; # whether periodic boundary conditions are on

#struna=1; # whether string plugin is active
#strunaConfig='@{sminput}'; # config name for string plugin
#strunaLog='@{smlog}'; # lof name for string plugin
dynamo=1; # whether string plugin is active
dynamoConfig='@{sminput}'; # config name for string plugin
dynamoLog='@{smlog}'; # lof name for string plugin

mini=0;          # whether to minimize before dynamics
ministeps=@{ministeps};   # number of minimization iterations

nsteps=@{nsteps};    # number of simulation steps
outputfreq=int(@{freq}/10); # frequency of generating output
dcdfreq=@{freq};   # frequency of dcd output

flag=''

outputName='@{outputname}';

exec(open('CHOMM.py').read())
