#!/bin/python2
from __future__ import print_function
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u
from sys import stdout, stderr, exit, argv
from shutil import move
import random


if len(argv) < 3 :
 raise ValueError("USAGE : ./md.py <structure_name> <param-dir>")

myname=argv[0];
name=argv[1];
prmdir=argv[2];

psffile=name+'.psf' ;
corfile=name+'_m.cor' ;
implicitSolvent=1 ;# run OBC2 implicit solvent simulation

hmass=4;       # amu, can use heavy hydrogens
friction=0.1   # units of 1/ps, Langevin thermostat
dt=4;          # timestep in fs
cutoff=20;      # nonbonded cutoff (should be around 12)
switchdist=cutoff-1.5 ; # switching distance
shake=1; # whether to constrain bonds involving hydrogens

thermostat=1;  # whether to use a thermostat
temperature=298; # kelvin

removeCOM=0; # whether to remove COM motion

mini=1;          # whether to minimize before dynamics
ministeps=0;     # number of minimization iterations

mdsteps=100000000;# number of production steps
outputfreq=10000;  # frequency of generating output
dcdfreq=2500;     # frequency of dcd output

restart=0
#restartfile='./scratch/.xml'; # to restart
outputName=name;

platformName="OpenCL" ;
platformName="CUDA" ;
# ==============================================
# ========================== Aux print functions
def dprint(*args):
  print(" ===> mdvac.py : ",end="");
  for arg in args:
   print(arg,end="")
  print(""); # flush
# ==========================
def printe(simulation):
  forceGroups={'Bond':0, 'Angle':1, 'Dihed':2, 'UB':3, 'IMPR':4, 'CMAP':5, 'NBOND':6};
  ener={};
# evaluate
  for key in forceGroups:
   state=simulation.context.getState(getEnergy=True, groups=1<<forceGroups[key]) ;
   ener[key]=state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole);
# total potential energy
  ener['PE']=simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalories_per_mole);
# print
  for key in ['Bond', 'Angle', 'Dihed', 'UB', 'IMPR', 'CMAP', 'NBOND', 'PE']:
   print(key, end="\t\t\t");
  print(); # newline
  for key in ['Bond', 'Angle', 'Dihed', 'UB', 'IMPR', 'CMAP', 'NBOND', 'PE']:
   print(ener[key], end="\t");
  print();

#========================== Initialize simulation system
dprint("Reading PSF from file '", psffile, "'");
psf=app.CharmmPsfFile(psffile);
toppar=[prmdir+'/'+p for p in ["par_all36m_prot.prm", "top_all36_prot.rtf", "staples-toppar.str", "toppar_water_ions.str" ]];
#toppar=("../../util/par_all36m_prot.prm", "../../util/top_all36_prot.rtf", "../../util/staples-toppar.str" );
dprint( "Reading parameters/topology from ", *toppar)
params=app.CharmmParameterSet(*toppar, permissive=False);

if (cutoff>0):
  nbondMethod=app.CutoffNonPeriodic
else:
  nbondMethod=app.NoCutoff
#===================================================== SHAKE
if (shake==1):
  cons=app.HBonds
  rigidWater=True
  dprint("Will constrain all bonds involving hydrogens");
elif (shake>1):
  cons=app.AllBonds
  rigidWater=True
  dprint("Will constrain all bond lengths");
else:
  cons=None
  rigidWater=False
#
dprint("Initializing simulation system");
dprint("Nonbonded cutoff is ",cutoff*u.angstrom,". Switching is active at ",switchdist*u.angstrom)
#
if (hmass>1):
  dprint("Hydrogen mass is ",hmass*u.amu)
  if (implicitSolvent==1):
   system=psf.createSystem(params,
                         nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff*u.angstrom, switchDistance=switchdist*u.angstrom,
                         constraints=cons, rigidWater=rigidWater, removeCMMotion=removeCOM, hydrogenMass=hmass*u.amu,
                         implicitSolvent=app.OBC2,
                         verbose=True);
  else:
   system=psf.createSystem(params,
                         nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff*u.angstrom, switchDistance=switchdist*u.angstrom,
                         constraints=cons, removeCMMotion=removeCOM, hydrogenMass=hmass*u.amu, rigidWater=rigidWater,
                         verbose=True);
else:
  if (implicitSolvent==1):
   system=psf.createSystem(params,
                         nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff*u.angstrom, switchDistance=switchdist*u.angstrom,
                         constraints=cons, rigidWater=rigidWater, removeCMMotion=removeCOM,
                         implicitSolvent=app.OBC2,
                         verbose=True);
  else:
   system=psf.createSystem(params,
                         nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff*u.angstrom, switchDistance=switchdist*u.angstrom,
                         constraints=cons, removeCMMotion=removeCOM, rigidWater=rigidWater,
                         verbose=True);

# integrator
if (thermostat):
    dprint("Initializing Langevin thermostatted integrator with timestep ",dt*u.femtosecond," coupled to bath with friction ",friction/u.picosecond," at temperature ",temperature*u.kelvin);
    integrator=mm.LangevinIntegrator(temperature*u.kelvin, friction/u.picosecond, dt*u.femtosecond);
else:
   dprint("Initializing Verlet integrator with timestep ",dt*u.femtosecond);
   integrator=mm.VerletIntegrator(dt*u.femtosecond);
#
# compute platform (e.g. CPU v GPU)
#
dprint("Initializing compute platform ",platformName);
platform=mm.Platform.getPlatformByName(platformName);
dprint("Preparing simulation topology");
if (platformName=="CUDA") :
  properties={'CudaPrecision': 'mixed'};
  simulation=app.Simulation(psf.topology, system, integrator, platform, properties);
elif (platformName=="OpenCL") :
  properties={'OpenCLPrecision': 'mixed'};
  simulation=app.Simulation(psf.topology, system, integrator, platform, properties);
else :
  simulation=app.Simulation(psf.topology, system, integrator, platform);

if (restart == 0) :
  if (corfile!=None):
   cor=app.CharmmCrdFile(corfile);
   dprint("Setting simulation coordinates from file '",corfile,"'");
   simulation.context.setPositions(cor.positions);
   try :
    velcorfile
   except NameError:
    velcorfile=None
   if (velcorfile!=None):
    vel=app.CharmmCrdFile(velcorfile);
    dprint("Setting simulation velocities from file '",velcorfile,"'");
    simulation.context.setVelocities(vel.positions);
  else:
   pdb=app.PDBFile(pdbfile);
   dprint("Setting simulation coordinates from file '",pdbfile,"'");
   simulation.context.setPositions(pdb.positions);
   if (velpdbfile!=None):
    vel=app.PDBFile(velpdbfile);
    dprint("Setting simulation velocities from file '",velpdbfile,"'");
    simulation.context.setVelocities(vel.positions);
else :
  dprint("Setting simulation restart data from file '",restartfile,"'");
  with open(restartfile, 'r') as f:
   xml=f.read();
   oldstate=mm.XmlSerializer.deserialize(xml)
   simulation.context.setPositions(oldstate.getPositions());
   simulation.context.setVelocities(oldstate.getVelocities());
   simulation.context.setTime(oldstate.getTime());
#
#================ Print initial energy compoments :
dprint("Initial Potential energy" );
printe(simulation);
#================ Energy minimization
if (mini) :
  dprint("Minimizing energy for ",ministeps," steps");
  simulation.minimizeEnergy(maxIterations=ministeps); # optional maxIterations, tolerance
  dprint("Potential energy after minimization");
  printe(simulation);
#=============== MD simulation
if (mdsteps>0):
  try :
   qrandname
# randomize temporary dcd names to avoid overwrite if running several calculations in the current dir
  except NameError:
   qrandname=1
  if (qrandname):
   outdcd='output_'+str(random.randint(1,10000))+'.dcd';
  else:
   outdcd='output.dcd'
#
  simulation.reporters.append(app.DCDReporter(outdcd,dcdfreq));
  simulation.reporters.append(app.StateDataReporter(stdout, outputfreq, step=True, potentialEnergy=True, kineticEnergy=True, speed=True, temperature=True, 
                                                    volume=0, separator=' \t '));
  dprint("Running MD simulation for ",mdsteps," steps");
  simulation.step(mdsteps);
#==== move dcd file to destination file
  move(outdcd, outputName+'.dcd');

dprint("Writing simulation restart files");
simulation.saveState(outputName+'.xml');
# simulation.saveCheckpoint(outputName+'.chk'); # usually do not need this file
#==== write periodic box vectors
state=simulation.context.getState();
a,b,c=state.getPeriodicBoxVectors();
fxsc=open(outputName+'.xsc','w');
fxsc.write("#CHOMMPy.xsc stub\n");
fxsc.write(str(mdsteps)+" "+str(a[0].value_in_unit(u.angstrom))+" 0 0 0 "+str(b[1].value_in_unit(u.angstrom))+" 0 0 0 "+str(c[2].value_in_unit(u.angstrom))+" 0 0 0 0 0 0 0 0 0\n");
fxsc.close();
# ==== minor cleanup
del switchdist;
del system;
del simulation;
