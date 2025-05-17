#!/bin/python2
from __future__ import print_function
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u
from sys import stdout, stderr, exit, argv
import mdtraj as md

if len(argv) < 3 :
 raise ValueError("USAGE : ./mdene.py <psf_file> <dcd_file> <param-dir>")

myname=argv[0];
psffile=argv[1];
dcdfile=argv[2];
prmdir=argv[3];
pdbfile=argv[4]; # use PDB topology for MDTraj

implicitSolvent=1 ;# run OBC2 implicit solvent simulation
cutoff=20;      # nonbonded cutoff (should be around 12)
switchdist=cutoff-1.5 ; # switching distance

platformName="OpenCL" ;
platformName="CUDA" ;
# ==============================================
# ========================== Aux print functions
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
def printe(simulation, header=False):
# print
  if (header):
   for key in ['Bond', 'Angle', 'Dihed', 'UB', 'IMPR', 'CMAP', 'NBOND', 'PE']:
    print(key, end="\t\t");
   print(); # newline
  else:
   forceGroups={'Bond':0, 'Angle':1, 'Dihed':2, 'UB':3, 'IMPR':4, 'CMAP':5, 'NBOND':6};
   ener={};
# evaluate
   for key in forceGroups:
    state=simulation.context.getState(getEnergy=True, groups=1<<forceGroups[key]) ;
    ener[key]=state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole);
# total potential energy
   ener['PE']=simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalories_per_mole);
   for key in ['Bond', 'Angle', 'Dihed', 'UB', 'IMPR', 'CMAP', 'NBOND', 'PE']:
    print(ener[key], end="\t");
   print();
#========================== Initialize simulation system
dprint("Reading PSF from file '", psffile, "'");
psf=app.CharmmPsfFile(psffile);
#toppar=[prmdir+'/'+p for p in ["par_all36m_prot.prm", "top_all36_prot.rtf", "staples-toppar.str", "toppar_water_ions.str" ]];
toppar=[prmdir+'/'+p for p in ["c36j18.top", "c36j18.par"]];
#toppar=("../../util/par_all36m_prot.prm", "../../util/top_all36_prot.rtf", "../../util/staples-toppar.str" );
dprint( "Reading parameters/topology from ", *toppar)
params=app.CharmmParameterSet(*toppar, permissive=False);

if (cutoff>0):
  nbondMethod=app.CutoffNonPeriodic
else:
  nbondMethod=app.NoCutoff
#===================================================== SHAKE
dprint("Initializing simulation system");
dprint("Nonbonded cutoff is ",cutoff*u.angstrom,". Switching is active at ",switchdist*u.angstrom)
if (implicitSolvent==1):
   system=psf.createSystem(params,
                         nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff*u.angstrom, switchDistance=switchdist*u.angstrom,
                         implicitSolvent=app.OBC2, verbose=True);
else:
   system=psf.createSystem(params,
                         nonbondedMethod=nbondMethod, nonbondedCutoff=cutoff*u.angstrom, switchDistance=switchdist*u.angstrom,
                         verbose=True);

# Verlet integrator (not actually used but needed to initialize simulation)
integrator=mm.VerletIntegrator(1*u.femtosecond);
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


dprint("Reading simulation coordinates from file '",dcdfile,"'");
dcd=md.load(dcdfile,top=pdbfile);
dprint("Computing potential energy for ", dcd.n_frames," frames:" );
printe(simulation,header=True); # print header only
# loop over all frames and compute energy
for i in range(dcd.n_frames) : # number of frames
# dprint("Reading coirdinate frame: ",i);
 simulation.context.setPositions(dcd.xyz[i]);
 printe(simulation);
