#!/usr/bin/env python3

import os, sys
import math
import mdtraj
import openmm
from openmm import unit
import openmmscripts as omms
import numpy as np

import logging
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


class DistanceCentroidForceReporter(object):
    def __init__(self, file, reportInterval, append=False):
        self._reportInterval = reportInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            if append:
                self._out = open(file, 'a')
            else:
                self._out = open(file, 'w')
        else:
            self._out = file
        print("#%17s %18s %18s %18s %18s %18s" % ('Step', 'd0', 'K', 'E', 'd1', 'd2'), file=self._out)

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, False, True)

    def report(self, simulation, state):
        K, d0, energy, distance = DistanceCentroidForce_getDistance(simulation)
        distance2 = DistanceCentroidZeroForce_getDistance(simulation)
        outstr = '%18d %+18.10E %+18.10E %+18.10E %+18.10E %18.10E' % (
            simulation.currentStep, d0.value_in_unit(unit.angstrom), 
            K.value_in_unit(unit.kilocalories_per_mole/unit.angstrom**2),
            energy.value_in_unit(unit.kilocalories_per_mole),
            distance.value_in_unit(unit.angstrom), distance2.value_in_unit(unit.angstrom))
        print(outstr, file=self._out)
        try:
            self._out.flush()
        except AttributeError:
            pass

    def __del__(self):
        if self._openedFile:
            self._out.close()

def DistanceCentroidForce(group1_index, group2_index, usepbc=True):
    energy_expression = '(DistanceCentroidForce_K/2)*(distance(g1, g2)-DistanceCentroidForce_d0)^2;';
    restrain_force = openmm.CustomCentroidBondForce(2, energy_expression)
    restrain_force.setUsesPeriodicBoundaryConditions(usepbc)
    restrain_force.addGroup(group1_index)
    restrain_force.addGroup(group2_index)
    restrain_force.addBond([0,1])
    restrain_force.addGlobalParameter('DistanceCentroidForce_K', 0.0)
    restrain_force.addGlobalParameter('DistanceCentroidForce_d0', 0.0)
    restrain_force.setName('DistanceCentroidForce')
    return restrain_force

def DistanceCentroidForce_getDistance(simulation):
    possible_forces = list()
    for force in simulation.system.getForces():
        if force.getName()==('DistanceCentroidForce'):
            possible_forces.append(force.getForceGroup())
    if len(possible_forces)!=1:
        #raise ValueError("Cannot find the correct force")
        return 0*unit.kilojoules_per_mole/unit.nanometers**2, 0*unit.nanometer, 0*unit.kilojoules_per_mole, 0*unit.nanometer
    force_group = possible_forces[0]
    energy = simulation.context.getState(getEnergy=True, groups={force_group}).getPotentialEnergy()
    d0 = simulation.context.getParameter('DistanceCentroidForce_d0') * unit.nanometer
    K = simulation.context.getParameter('DistanceCentroidForce_K') * unit.kilojoules_per_mole/unit.nanometers**2
    tmp = 2.0*energy/K
    distance = d0 + math.sqrt(tmp.value_in_unit(unit.nanometer**2))*unit.nanometer
    return K, d0, energy, distance

def DistanceCentroidZeroForce(group1_index, group2_index, usepbc=True):
    energy_expression = 'DistanceCentroidZeroForce_K*distance(g1, g2);';
    restrain_force = openmm.CustomCentroidBondForce(2, energy_expression)
    restrain_force.setUsesPeriodicBoundaryConditions(usepbc)
    restrain_force.addGroup(group1_index)
    restrain_force.addGroup(group2_index)
    restrain_force.addBond([0,1])
    restrain_force.addGlobalParameter('DistanceCentroidZeroForce_K', 0.0)
    restrain_force.setName('DistanceCentroidZeroForce')
    return restrain_force

def DistanceCentroidZeroForce_getDistance(simulation):
    possible_forces = list()
    for force in simulation.system.getForces():
        if force.getName()==('DistanceCentroidZeroForce'):
            possible_forces.append(force.getForceGroup())
    if len(possible_forces)!=1:
        #raise ValueError("Cannot find correct force")
        return 0*unit.nanometer
    force_group = possible_forces[0]
    simulation.context.setParameter('DistanceCentroidZeroForce_K', 1.0)
    force = simulation.context.getState(getEnergy=True, groups={force_group}).getPotentialEnergy()
    distance = force / unit.kilojoules_per_mole * unit.nanometer
    simulation.context.setParameter('DistanceCentroidZeroForce_K', 0.0)
    return distance




def make_system_obc2_nocutoff(psffile, corfile):
    #cutoffmethod = openmm.app.PME
    #cutoff = 12*unit.angstrom
    #switchdist = 10*unit.angstrom
    psf = openmm.app.CharmmPsfFile(psffile)
    pdb = openmm.app.CharmmCrdFile(corfile)
    #boxx, boxy, boxz = boxsize
    #psf.setBox(boxx, boxy, boxz)
    #for i, _ in enumerate(pdb.positions):
    #    pdb.positions[i] += [0.5*boxx]*3*unit.nanometer
    params = openmm.app.CharmmParameterSet(os.environ['HOME']+'/opt/ff/charmmff_jul18.prm')
    system = psf.createSystem(params, nonbondedMethod=openmm.app.NoCutoff, constraints=openmm.app.AllBonds, rigidWater=True, removeCMMotion=False, implicitSolvent=openmm.app.OBC2)
    return psf.topology, pdb.positions, system

def make_system_obc2_cutoff(psffile, corfile, cutoff=12*unit.angstrom, switchdist=10*unit.angstrom):
    psf = openmm.app.CharmmPsfFile(psffile)
    pdb = openmm.app.CharmmCrdFile(corfile)
    params = openmm.app.CharmmParameterSet(os.environ['HOME']+'/opt/ff/charmmff_jul18.prm')
    system = psf.createSystem(params, nonbondedMethod=openmm.app.CutoffNonPeriodic, constraints=openmm.app.AllBonds, 
                rigidWater=True, removeCMMotion=False, implicitSolvent=openmm.app.OBC2, nonbondedCutoff=cutoff, switchDistance=switchdist)
    return psf.topology, pdb.positions, system


def make_simulation(topology, coordinates, system, prevxml=None, temperature=300.0*unit.kelvin, timestep=2.0*unit.femtosecond):
    integrator = openmm.LangevinIntegrator(temperature, 10/unit.picosecond, timestep)
    pltform = openmm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed' , "CudaDeviceIndex" : '0'}
    simulation = openmm.app.Simulation(topology, system, integrator, pltform, properties)
    if prevxml:
        simulation.loadState(prevxml)
    else:
        simulation.context.setPositions(coordinates)
        simulation.context.setVelocitiesToTemperature(temperature)
    return simulation

def add_distance_restrain(system, topology, coordinates, sel1, sel2, usepbc=True):
    mdtraj_topology = mdtraj.Topology().from_openmm(topology)
    ndx1 = mdtraj_topology.select(sel1).tolist()
    ndx2 = mdtraj_topology.select(sel2).tolist()
    center1 = np.average(np.array(coordinates.value_in_unit_system(unit.md_unit_system))[ndx1], axis=0)
    center2 = np.average(np.array(coordinates.value_in_unit_system(unit.md_unit_system))[ndx2], axis=0)
    distance = np.linalg.norm(center1 - center2)
    force = DistanceCentroidForce(ndx1, ndx2, usepbc)
    force.setForceGroup(omms.BIAS_FORCE_GROUP)
    system.addForce(force)
    force2 = DistanceCentroidZeroForce(ndx1, ndx2, usepbc)
    force2.setForceGroup(omms.BIAS_FORCE_GROUP+1)
    system.addForce(force2)
    return system, distance*unit.nanometer

def minimize(simulation, ministep=1000, minitol=False):
    energy = omms.get_energy(simulation.context)
    log.info('Energy Before Minimization is %s' % (str(energy)))
    omms.print_energy_decomposition(simulation.context)
    log.info("Running energy minimization (maxsteps = %s tolerance = %s)" % (ministep, minitol))
    simulation.minimizeEnergy(maxIterations=ministep, tolerance=minitol)
    energy = omms.get_energy(simulation.context)
    log.info('Energy After Minimization is %s' % energy)
    omms.print_energy_decomposition(simulation.context)

def add_reporters(simulation, nameout, simultime=1*unit.nanosecond, logtime=10*unit.picosecond, timestep=2*unit.femtoseconds):
    log.info("Adding reporters...")
    totstep = int(round(simultime/timestep, 0))
    logfreq = int(round(logtime/timestep, 0))
    simulation.reporters.append(openmm.app.StateDataReporter(sys.stdout, logfreq, step=True, time=True, 
        progress=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, 
        totalEnergy=True, temperature=True, speed=True, totalSteps=totstep, separator=' '))
    log.info("Adding OpenMM log file reporter to %s every %s (%d steps)" % (nameout+'.log', logtime, logfreq))
    simulation.reporters.append(openmm.app.StateDataReporter(nameout+'.log', logfreq, step=True, time=True, 
        progress=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, 
        totalEnergy=True, temperature=True, speed=True, totalSteps=totstep, separator=' '))
    log.info("Adding DCD reporter to %s every %s (%d steps)" % (nameout+'.dcd', logtime, logfreq))
    simulation.reporters.append(openmm.app.DCDReporter(nameout+'.dcd', logfreq, append=False))
    log.info("Adding Energy log file to %s every %s (%d steps)" % (nameout+'.ene', logtime, logfreq))
    simulation.reporters.append(omms.EnergyReporter(nameout+'.ene', logfreq, append=False))
    simulation.reporters.append(DistanceCentroidForceReporter(nameout+'.dist', logfreq, append=False))

def smd(simulation, simultime, timestep, startdistance, enddistance, timeincrement):
    numincrements = int(round(simultime/timeincrement, 0))
    stepdistance = enddistance/numincrements
    stepstep = int(round(simultime/timestep/numincrements, 0))
    if stepstep*timestep!=timeincrement:
        raise ValueError(":(")
    log.info("Running for a total of %s (%d steps), pulling from %s to %s, incrementing the distance by %s every %s (%d steps)." % (simultime, stepstep*numincrements, startdistance, startdistance+enddistance, stepdistance, timeincrement, stepstep))
    for i in range(numincrements):
        d0 = startdistance + stepdistance*(i+1)
        K = 1000*unit.kilocalories_per_mole/(unit.angstrom**2)
        simulation.context.setParameter('DistanceCentroidForce_d0', d0.value_in_unit_system(unit.md_unit_system))
        simulation.context.setParameter('DistanceCentroidForce_K', K.value_in_unit_system(unit.md_unit_system))
        simulation.step(stepstep)


def main():
    top, cor, system = make_system_obc2_cutoff('model/model-chm.psf', 'model/model-chm.cor')
    system, startdistance = add_distance_restrain(system, top, cor, 'chainid 0', 'chainid 1', usepbc=False)
    simulation = make_simulation(top, cor, system)
    minimize(simulation)
    add_reporters(simulation, nameout='smd/run1')
    simultime = 1*unit.nanosecond
    timestep = 2*unit.femtosecond
    smd(simulation, simultime, timestep, startdistance, 3*unit.nanometers, 1*unit.picoseconds/50)

    #log.info("Saving status as restart file to %s" % (nameout+'.xml'))
    #simulation.saveState(nameout+'.xml')
    #log.info("Normal termination :)")

main()
