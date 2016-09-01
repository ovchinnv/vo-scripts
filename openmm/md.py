#!/bin/python
#requires the file CHOMM.py, which is simple wrapper function to run MD using OpenMM using CHARMM parameters
#=====================================================================================
#
# aux parameters (e.g. they help define the required ones, but are not themselves used by CHOMM)
irun=0     ;# initial run index
nrun=1     ;# number of runs
name='emr-apo-dmpc' ;# prefix for output files
#==============================
# parameters required by CHOMM (some have default values)
restart=0 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
restartfile=None ;

psffile='struc/'+name+'.psf' ;
pdbfile='struc/'+name+'.res' ;
paramfile='param/emr-tpp-dmpc.prm' ;
#boxfile='struc/boxsize.str';
dx=67.759903; # can specify manually
dy=dx;
dz=93.37640;

hmass=4;         # amu, can use heavy hydrogens
friction=0.1;    # 1/ps, thermostat coupling
dt=3.5;          # timestep in fs ; cannot go up tp 4 due to crashes as the system approaches 300K; keeping <= 3.5
pmefreq=1;       # >1 requires multiple timestepping (MTS) which actually slows the code down compared wo pmefreq=1 case
cutoff=9;        # nonbonded cutoff
switchdist=7.25;

constraints=0;   # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=pdbfile;# as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology

shake=1 ; # whether to constrain bonds involving hydrogens; shape=2 means AllBonds
#shake=0 ;

thermostat=1;  # whether to use a thermostat
temperature=298; # kelvin
andersen=0;    # to use Andersen instead of Langevin
barostat=0;
pressure=1;    # atm
membrane_on=1; # whether to use a barostat for membrane simulations (z-axis is the membrane normal)
pme=1; # whether to use PME (valid only if pbc=1)
pbc=1; # whether periodic boundary conditions are on

mini=(irun==0);  # whether to minimize before dynamics
ministeps=100;   # number of minimization iterations

nsteps=1000000;  # number of simulation steps
outputfreq=1000; # frequency of generating output
dcdfreq=10000;   #frequency of dcd output
outputName=name+str(irun) ;

# run MD
exec(open('CHOMM.py').read())
