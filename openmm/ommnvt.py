#!/bin/python
#
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, stderr, exit

# input parameters :
irun=0 ;# initial run index
nrun=10 ;
restart=0 ; # 0 -- start from PDB coordinates; 1 -- restart from native xml file
restartfile=None ;

psffile='struc/fre-fad-x.psf' ;
pdbfile='struc/fre_fad.pdb' ;
paramfile='param/fre.prm' ;
boxfile='struc/boxsize.str';

temperature=298*kelvin; # kelvin
hmass=1*amu; # amu, can use heavy hydrogens
friction=1/picosecond; # thermostat coupling
dt=2*femtosecond; # timestep in fs
pmefreq=1; # >1 requires multiple timestepping, which _dramatically_ slows down the code
cutoff=12*angstrom; # nonbonded cutoff

constraints=0; # harmonic positional restraints for equilibration
constraintscaling=1; # to scale hatmonic restraints uniformly
consfile=pdbfile; # as in NAMD/ACEMD, this file must have identical atom ordering to that in the system topology

shake=1; # whether to constrain bonds involving hydrogens
#shake=0 ;

thermostat=1; # whether to use a thermostat
andersen=1; # use Andersen instead of Langevin
barostat=0;
membrane_on=0;# whether to use a barostat for membrane simulations
pme=1; # whether to use PME
pbc=1; # whether periodic boundary conditions are on

mini=(irun==0); # whether to minimize before dynamics
ministeps=100; # number of minimization iterations

nsteps=1000000; #number of simulation steps
freq=1000; # frequency of generating output
dcdfreq=10000; #frequency of dcd output
name='fre' ; #prefix for output files

#=====================================================================#
#============= No modifications are normally needed below this line ==#
#=====================================================================#
outputName='fre'+str(irun) ;
#=====================================================================
# define optional parameters that may not have been defined by user
#
try :
 mini
except NameError:
 mini=0
#
try :
 constraints
except NameError:
 constraints=0
#
try :
 pbc
except NameError:
 pbc=0
#
try :
 pme
except NameError:
 pme=0
#
try :
 thermostat
except NameError:
 thermostat=0
try :
 andersen
except NameError:
 andersen=0
#
try :
 barostat
except NameError:
 barostat=0
try :
 membrane_on
except NameError:
 membrane_on=0
#
try :
 platform
except NameError:
# use CUDA unless variable 'platform' defined
 platform="CUDA"
#==========================
def dprint(*args):
 print(" ===> CHOMMPy : ",end="");
 for arg in args:
  print(arg,end="")
 print(" ...")
#========================== Initialize simulation system
dprint("Reading PSF from file '", psffile, "'");
psf=CharmmPsfFile(psffile);
dprint("Reading parameter file '", paramfile,"'");
params=CharmmParameterSet(paramfile, permissive=True); # need permissive to avoid providing atom types ( a la xplor psf )
#========================================================
if (pbc):
 dprint("Periodic boundary conditions will be used")
 if (not restart):
  try :
   dx; dy; dz;
  except NameError:
   dprint("Setting orthorhombic cell lengths from file '",boxfile,"'")
   with open(boxfile) as f:
    lines=f.readlines()
    i=0;
    for line in lines:
     words = line.split()
     if (words[0].upper() == 'SET') :
      break;
     i=i+1;
    # get  box size
    dx = float(lines[i].split()[2])*angstrom;
    dy = float(lines[i+1].split()[2])*angstrom;
    dz = float(lines[i+2].split()[2])*angstrom;
  else:
   pass
# will need to get dimensions from the xml file !
 dprint("Periodic cell dimensions are (", dx, ")x(", dy, ")x(", dz,")")
 psf.setBox(dx, dy, dz);
 nbondMethod=PME
else :
 nbondMethod=CutoffNonPeriodic
#===================================================== SHAKE
if shake:
 cons=HBonds
 dprint("Will constrain bonds involving hydrogens");
else:
 cons=None
dprint("Initializing simulation system");
dprint("Nonbonded cutoff is ",cutoff)
if (hmass>1*amu):
 dprint("Hydrogen mass is ",hmass)
system=psf.createSystem(params,
                        nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff, switchDistance=cutoff-1.5*angstrom,
                        constraints=cons, removeCMMotion=False, hydrogenMass=hmass,
                        verbose=True);

#================= harmonic restraints from file, a la NAMD/ACEMD
if (constraints) :
 dprint("Adding absolute positional harmonic restraints to atoms marked in beta column of pdf file '"+consfile+"'");
 force=CustomExternalForce("s*0.5*k*periodicdistance(x,y,z,x0,y0,z0)");
# force=CustomExternalForce("s*0.5*k*( (x-x0)^2 + (y-y0)^2 + (z-z0)^2 )");
 force.addPerParticleParameter("k");
 force.addPerParticleParameter("x0");
 force.addPerParticleParameter("y0");
 force.addPerParticleParameter("z0");
 force.addGlobalParameter("s", constraintscaling);
# read per atom restraints :
 res=PDBFile(consfile);
 iatom=0; icons=0;
 for r, o, b  in zip(res.positions, res.occupancy, res.temperature_factor) :
  bnodim=b/angstrom/angstrom; # have to deal with units, which are A^2 for B-factors
  if ( bnodim > 0 ) :
   icons+=1;
   k=bnodim*constraintscaling*kilocalorie/mole/angstrom/angstrom
#   dprint(" Adding restraint on atom ",iatom," with force constant ", k );
   force.addParticle(iatom, [k,r[0],r[1],r[2]]);
  iatom+=1;
 dprint("Added restraints on ", icons, " atoms");
 dprint("Harmonic force constants will be scaled uniformly by x"+str(constraintscaling));
 system.addForce(force)

#================= add integrator :
dprint("Configuring integrator");
if (pme and pmefreq > 1) :
#================= multiple timestepping :
#================= first, add thermostat force
    if (thermostat):
     dprint("Initializing Anderson Thermostat with coupling to bath with friction ",friction," at temperature ",temperature)
     thermostatForce=AndersenThermostat(temperature, friction);
     thermostatForce.setForceGroup(0) ; # make sure to assign a group for MTS
     system.addForce(thermostatForce);
    dprint("Multiple time stepping (MTS) will be used");
# split force objects into a multiple groups
# note that createSystem puts different psf sections into different force groups for ease of energy decomposition;
    for f in system.getForces() :
# put all forces into the same group :
# we should be able to use many groups withe same substep in the RESPA init, but that might slow it down
     f.setForceGroup(0);
# reciprocal forces get a separate group for RESPA
     if isinstance(f,NonbondedForce) :
      f.setReciprocalSpaceForceGroup(31);
    dprint("Initializing RESPA MTS integrator");
    integrator=MTSIntegrator(dt, [(31,1), (0,2)]);
else :
    if (thermostat):
     if (andersen):
      dprint("Initializing Andersen integrator with timestep ",dt," coupled to bath with friction ",friction," at temperature ",temperature);
      integrator=LangevinIntegrator(temperature, friction, dt);
     else:
      dprint("Initializing Langevin integrator with timestep ",dt," coupled to bath with friction ",friction," at temperature ",temperature);
      integrator=LangevinIntegrator(temperature, friction, dt);
    else:
     dprint("Initializing Verlet integrator");
     integrator=VerletIntegrator(dt);
#
#====================================================
#
dprint("Initializing compute platform ",platform);
platform=Platform.getPlatformByName(platform);
properties={'CudaPrecision': 'mixed'};
dprint("Preparing simulation topology");
if (platform=="CUDA") :
 simulation=Simulation(psf.topology, system, integrator, platform, properties);
else :
 simulation=Simulation(psf.topology, system, integrator);

if (restart == 0) :
 dprint("Setting simulation coordinates from file '"+pdbfile+"'");
 pdb=PDBFile(pdbfile);
 simulation.context.setPositions(pdb.positions);
else :
 dprint("Setting simulation coordinates from file '"+restartfile+"'");
 simulation.loadState(restartfile);

#================ Print initial energy
state=simulation.context.getState(getEnergy=True) ;
dprint("Initial energy : ",state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)*kilocalories_per_mole)
#state=simulation.context.getState();
#a,b,c=state.getPeriodicBoxVectors();
#fxsc=open(outputName+'.xsc','w');
#fxsc.write("#CHOMMPy xsc stub\n");
#fxsc.write(str(nsteps)+" "+str(a[0].value_in_unit(angstrom))+" 0 0 0 "+str(b[1].value_in_unit(angstrom))+" 0 0 0 "+str(c[2].value_in_unit(angstrom))+" 0 0 0 0 0 0 0 0 0\n");
#fxsc.close();
#quit()
#================ Energy minimization
if (mini) :
 dprint("Minimizing energy for ",ministeps," steps");
 simulation.minimizeEnergy(maxIterations=ministeps, tolerance = 0*kilocalorie/mole); # optional iterations, tolerance
 state=simulation.context.getState(getEnergy=True) ;
 dprint( "Energy after minimization: ",state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)*kilocalories_per_mole)
#=============== MD simulation
if (nsteps>0):
 simulation.reporters.append(DCDReporter(outputName+'.dcd',dcdfreq));
 simulation.reporters.append(StateDataReporter(stdout, freq, step=True, potentialEnergy=True, kineticEnergy=True, speed=True, temperature=True, separator=' '));
 dprint("Running MD simulation for ",nsteps," steps");
 simulation.step(nsteps);
 dprint("Writing simulation restart files");
 simulation.saveState(outputName+'.xml');
 simulation.saveCheckpoint(outputName+'.chk');
#==== write periodic box vectors
 state=simulation.context.getState();
 a,b,c=state.getPeriodicBoxVectors();
 fxsc=open(outputName+'.xsc','w');
 fxsc.write("#CHOMMPy.xsc stub\n");
 fxsc.write(str(nsteps)+" "+str(a[0].value_in_unit(angstrom))+" 0 0 0 "+str(b[1].value_in_unit(angstrom))+" 0 0 0 "+str(c[2].value_in_unit(angstrom))+" 0 0 0 0 0 0 0 0 0\n");
 fxsc.close();

