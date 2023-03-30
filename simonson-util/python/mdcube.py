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

psffile=name+'s.psf' ;
corfile=name+'_ms.cor' ;
implicitSolvent=0 ;# run OBC2 implicit solvent simulation

hmass=4;       # amu, can use heavy hydrogens
friction=0.1   # units of 1/ps, Langevin thermostat
dt=4;          # timestep in fs
cutoff=10;      # nonbonded cutoff (should be around 12 for exp. solvent, >=20 for OBC2)
switchdist=cutoff-1.5 ; # switching distance
shake=1; # whether to constrain bonds involving hydrogens

barostat=1;
pressure=1;    # units of atm
pme=1; # whether to use PME
pbc=1; # whether periodic boundary conditions are on
boxfile=name+'.str'

thermostat=1;  # whether to use a thermostat
temperature=298; # kelvin

removeCOM=0; # whether to remove COM motion

mini=1;          # whether to minimize before dynamics
ministeps=100;     # number of minimization iterations; 0 means use energy criterion

mdsteps=100000000;# number of production steps
outputfreq=1000;  # frequency of generating output
dcdfreq=2500;     # frequency of dcd output

restart=0
#restartfile='./scratch/.xml'; # to restart
outputName=name+'wat';

platformName="OpenCL" ;
platformName="CUDA" ;
# ==============================================
# ========================== Aux functions
def dprint(*args):
  print(" ===>",myname,": ",end="");
  for arg in args:
   print(arg,end="")
  print(""); # flush
#==========================
def derror(*args):
  print(" ===>",myname,"ERROR : ",end="");
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
#========================== box dimensions from .str file produced during system preparation
def get_box_size_str(boxfile):
  with open(boxfile) as f:
     lines=f.readlines()
     i=0;
     for line in lines:
      words = line.split()
      if (words[0].upper() == 'SET') :
       break;
      i=i+1;
    # get box size
     dx = float(lines[i].split()[2])+1.5; # box padding
     dy = float(lines[i+1].split()[2])+1.5;
     dz = float(lines[i+2].split()[2])+1.5;
  return(dx, dy, dz);
#========================== box dimensions from .xml restart file produced with openmm
def get_box_size_xml(xmlfile):
  with open(xmlfile) as f:
     lines=f.readlines()
     i=1;
     for line in lines:
      words = line.split()
      if (words[0].upper() == '<PERIODICBOXVECTORS>') :
       break;
      i=i+1;
    # get  box size
     dx = (float(lines[i].split('"')[1]))*u.nanometers.conversion_factor_to(u.angstrom)  ; # remember that these are in nanoneters
     dy = (float(lines[i+1].split('"')[3]))*u.nanometers.conversion_factor_to(u.angstrom) ;
     dz = (float(lines[i+2].split('"')[5]))*u.nanometers.conversion_factor_to(u.angstrom) ;
  return(dx, dy, dz);
#========================== box dimensions from .xsc file produced with ACEMD/NAMD programs
def get_box_size_xsc(xscfile):
  with open(xscfile) as f:
     lines=f.readlines()
     for line in lines:
      pass
     cell=line.split(' ')
     dx = float(cell[1]);
     dy = float(cell[5]);
     dz = float(cell[9]);
  return(dx, dy, dz);
#
#========================== Initialize simulation system
dprint("Reading PSF from file '", psffile, "'");
psf=app.CharmmPsfFile(psffile);
toppar=[prmdir+'/'+p for p in ["par_all36m_prot.prm", "top_all36_prot.rtf", "staples-toppar.str", "toppar_water_ions.str" ]];
#toppar=("../../util/par_all36m_prot.prm", "../../util/top_all36_prot.rtf", "../../util/staples-toppar.str", \
# "../../util/toppar_water_ions.str" );
dprint( "Reading parameters/topology from ", *toppar)
params=app.CharmmParameterSet(*toppar, permissive=False);
#
#========================================================
if (pbc):
  dprint("Periodic boundary conditions will be used")
  if (not restart or resetcell):
   try :
    dx; dy; dz; # check if dimensions are specified manually
   except NameError:
    try :
     boxfile;
     dprint("Setting orthorhombic cell lengths from file '",boxfile,"'")
     dx, dy, dz=get_box_size_str(boxfile);
    except NameError:
     try:
      xscfile;
      dprint("Setting orthorhombic cell lengths from file '",xscfile,"'")
      dx, dy, dz=get_box_size_xsc(xscfile);
     except NameError:
      try:
       xmlfile;
       dprint("Setting orthorhombic cell lengths from file '",xmlfile,"'")
       dx, dy, dz=get_box_size_xml(xmlfile);
      except NameError:
       derror("Could not set periodic cell size.")
  else:
   try:
    restartfile;
    dprint("Setting orthorhombic cell lengths from file '",restartfile,"'")
    dx, dy, dz=get_box_size_xml(restartfile);
   except NameError:
    try:
     xmlfile;
     dprint("Setting orthorhombic cell lengths from file '",xmlfile,"'")
     dx, dy, dz=get_box_size_xml(xmlfile);
    except NameError:
     derror("Could not set periodic cell size.")
#
  try:
   dprint("Periodic cell dimensions are (", dx*u.angstrom, ")x(", dy*u.angstrom, ")x(", dz*u.angstrom,")")
   psf.setBox(dx*u.angstrom, dy*u.angstrom, dz*u.angstrom);
  except Exception:
   exit(-1)
#===========================================================
  if (pme):
   nbondMethod=app.PME
   dprint("PME is on");
  else:
   nbondMethod=app.CutoffPeriodic
   dprint("PME is off");
else:
  if (cutoff>0):
   nbondMethod=app.CutoffNonPeriodic
#   psf.setBox(1000*u.angstrom, 1000*u.angstrom, 1000*u.angstrom) # set to a very large box to eliminate wrapping
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
    if (barostat):
     dprint("Initializing Monte-Carlo barostat at pressure ",pressure*u.atmosphere)
     barostatForce=mm.MonteCarloBarostat(pressure*u.atmosphere, temperature*u.kelvin)
     system.addForce(barostatForce) ;
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
                                                    volume=pbc, separator=' \t '));
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
fxsc.write("#xsc stub\n");
fxsc.write(str(mdsteps)+" "+str(a[0].value_in_unit(u.angstrom))+" 0 0 0 "+str(b[1].value_in_unit(u.angstrom))+" 0 0 0 "+str(c[2].value_in_unit(u.angstrom))+" 0 0 0 0 0 0 0 0 0\n");
fxsc.close();
# ==== minor cleanup
del switchdist;
del system;
del simulation;
