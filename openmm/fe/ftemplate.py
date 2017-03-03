#!/bin/python
#=====================================================================================
name='traf-lig27'

psffile='./struc/'+name+'_sn.psf' ;
pdbfile=name+'-mini.pdb' ;
topfile='./struc/'+name+'.top';
paramfile='./struc/'+name+'.par' ;

flag='fep';
hmass=4;       # amu, can use heavy hydrogens
friction=0.1   # 1/ps, thermostat coupling
dt=4;          # timestep in fs
pmefreq=1;     # >1 requires multiple timestepping, which slows down the code substantially
cutoff=9;      # nonbonded cutoff
switchdist=7.5 ; # (optional) switching distance

constraints=0;   # harmonic positional restraints for equilibration

shake=1; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=298; # kelvin
barostat=1;
pressure=1;    # units of atm
pme=1; # whether to use PME
pbc=1; # whether periodic boundary conditions are on

dx=72.88126833885093
dy=dx
dz=dx

nsteps=100000;    # number of production steps
outputfreq=10000;  # frequency of generating output
dcdfreq=10000;     # frequency of dcd output

restart=@{restart}
restartfile='@{restartfile}'
outputName='@{outputname}'

#=== alchemical section
alch=1 ;         # alchemical transformation on or off
alchfile=name+'-fep.pdb' ;# PDB file with atoms marked
alchcol=1;       # column : 1 -- beta ; 2 -- occupancy
alchout='@{alchout}'; # file name for alchemical output
alchfreq=100;    # frequency of alchemical output
alchdecouple=0 ; # whether to decouple or to annihilate
lambda0=@{lambda0} ;# current lambda
lambda1=@{lambda1} ;# perturbation lambda

#=== conformational restraints via struna
struna=1
strunaConfig='@{strunaConfig}'
strunaLog='@{strunaLog}'
#
from os.path import expanduser
exec(open(expanduser('~/scripts/openmm/fe/CHOMM.py')).read())
