#!/usr/bin/env python3

import os, sys, math
import copy
import numpy as np
from pdbfixer import PDBFixer
import openmm
from openmm import unit
import mdtraj

import logging
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

rng = np.random.Generator(np.random.PCG64(4242))

def polyalanine(pdbfile):
    fixer = PDBFixer(pdbfile)
    for residue in fixer.topology.residues():
        if residue.name!='ALA':
            mutate = '%s-%s-ALA' % (residue.name, residue.id)
            fixer.applyMutations([mutate], residue.chain.id)
    fixer.findMissingResidues()                                         
    fixer.findMissingAtoms()                                            
    fixer.addMissingAtoms()                                             
    fixer.addMissingHydrogens(7.0)
    return fixer

def fixpdb(pdbfile):
    fixer = PDBFixer(pdbfile)
    fixer.findMissingResidues()                                         
    fixer.findMissingAtoms()                                            
    fixer.addMissingAtoms()                                             
    fixer.addMissingHydrogens(7.0)
    return fixer

def sequence(pdb):
    AA_3to1 = {
        'ALA':'A', 'CYS':'C', 'GLU':'E', 'ASP':'D', 'GLY':'G',
        'PHE':'F', 'ILE':'I', 'HIS':'H', 'LYS':'K', 'MET':'M',
        'LEU':'L', 'ASN':'N', 'GLN':'Q', 'PRO':'P', 'SER':'S',
        'ARG':'R', 'THR':'T', 'TRP':'W', 'VAL':'V', 'TYR':'Y',
        }
    return "".join([AA_3to1[residue.name] for residue in pdb.topology.residues()])

def solvate(pdbfile, solvname='solvated.pdb', buffer=12*unit.angstrom, ionicstrength=0.1*unit.molar):
    pdb = openmm.app.PDBFile(pdbfile)
    forcefield = openmm.app.ForceField('charmm36.xml', 'charmm36/water.xml')
    coor = np.array(pdb.positions.value_in_unit(unit.nanometers))
    center = np.average(coor, axis=0) * unit.nanometers
    bsize = np.amax(np.amax(coor, axis=0)-np.amin(coor, axis=0)) + buffer.value_in_unit(unit.nanometers)
    bbox = openmm.Vec3(bsize, bsize, bsize)*unit.nanometers
    for i, cor in enumerate(pdb.positions):
        pdb.positions[i] = cor-center + 0.5*bbox
    modeller = openmm.app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, boxSize=bbox, model='tip3p', ionicStrength=ionicstrength)
    with open(solvname, 'w') as fp:
        openmm.app.PDBFile.writeFile(modeller.topology, modeller.positions, fp)
    return

def make_system_pme(solvname='solvated.pdb', cutoff=12*unit.angstrom, switchdist=10*unit.angstrom):
    pdb = openmm.app.PDBFile(solvname)
    forcefield = openmm.app.ForceField('charmm36.xml', 'charmm36/water.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=openmm.app.PME, constraints=openmm.app.AllBonds, 
                rigidWater=True, removeCMMotion=False, nonbondedCutoff=cutoff, switchDistance=switchdist)
    return pdb, system

def minimize(system, pdb='solvated.pdb', pdboutfile='minimized.pdb', outstate='minimized.xml', ministep=1000, minitol=False):
    integrator = openmm.LangevinIntegrator(0*unit.kelvin, 10/unit.picosecond, 2*unit.femtoseconds)
    pltform = openmm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed' , "CudaDeviceIndex" : '0'}
    simulation = openmm.app.Simulation(pdb.topology, system, integrator, pltform, properties)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=ministep, tolerance=minitol)
    with open(pdboutfile, 'w') as fp:
        positions = simulation.context.getState(getPositions=True).getPositions()
        openmm.app.PDBFile.writeFile(pdb.topology, positions, fp)
    simulation.saveState(outstate)
    return

def plain_md(pdbfile, xmlstate, outname='heated', simultime=10*unit.nanoseconds, temperature=300*unit.kelvin, timestep=2*unit.femtoseconds, pressure=1*unit.bar):
    pdb, system = make_system_pme(pdbfile)
    system.addForce(openmm.MonteCarloBarostat(pressure, temperature))
    integrator = openmm.LangevinIntegrator(temperature, 10/unit.picosecond, timestep)
    pltform = openmm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed' , "CudaDeviceIndex" : '0'}
    simulation = openmm.app.Simulation(pdb.topology, system, integrator, pltform, properties)
    simulation.loadState(xmlstate)
    add_reporters(simulation, outname, simultime)
    simulation.currentStep = 0
    steps = int(simultime/timestep)
    simulation.step(steps)
    with open(outname+'.pdb', 'w') as fp:
        positions = simulation.context.getState(getPositions=True).getPositions()
        openmm.app.PDBFile.writeFile(pdb.topology, positions, fp)
    state = simulation.context.getState(getPositions=True, getVelocities=True, getParameters=False)
    with open(outname+'.xml', 'w') as fp:
        fp.write(openmm.XmlSerializer.serialize(state))
    return

def add_reporters(simulation, nameout, simultime=1*unit.nanosecond, logtime=100*unit.picosecond, timestep=2*unit.femtoseconds):
    totstep = int(round(simultime/timestep, 0))
    logfreq = int(round(logtime/timestep, 0))
    simulation.reporters.append(openmm.app.StateDataReporter(nameout+'.log', logfreq, step=True, time=True,
        progress=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, speed=True, totalSteps=totstep, separator=' '))
    simulation.reporters.append(openmm.app.DCDReporter(nameout+'.dcd', logfreq, append=False))




def run_md(pdb):
    seq = sequence(pdb)
    dirname = os.path.join('trpzip2', seq)
    if os.path.isfile(dirname+'/dyna.xml'):
        return
    os.makedirs(dirname, exist_ok=True)
    openmm.app.PDBFile.writeFile(pdb.topology, pdb.positions, open(dirname+'/initial.pdb', 'w'))
    solvate(dirname+'/initial.pdb', solvname=dirname+'/solvated.pdb', buffer=12*unit.angstrom, ionicstrength=0.1*unit.molar)
    pdb, system = make_system_pme(solvname=dirname+'/solvated.pdb', cutoff=12*unit.angstrom, switchdist=10*unit.angstrom)
    minimize(system, pdb, pdboutfile=dirname+'/minimized.pdb', outstate=dirname+'/minimized.xml', ministep=1000, minitol=False)
    plain_md(dirname+'/minimized.pdb', dirname+'/minimized.xml', outname=dirname+'/dyna', simultime=10*unit.nanoseconds, temperature=300*unit.kelvin, timestep=2*unit.femtoseconds, pressure=1*unit.bar)



def compute_rmsd(pdb):
    seq = sequence(pdb)
    dirname = os.path.join('trpzip2', seq)
    trj = mdtraj.load(dirname+'/dyna.dcd', top=dirname+'/minimized.pdb')
    ref = mdtraj.load(dirname+'/minimized.pdb')
    atoms = ref.topology.select("backbone")
    rmsd = 10.0*mdtraj.rmsd(trj, ref, atom_indices=atoms).flatten()
    N = int(len(rmsd)*0.5)
    avg = np.average(rmsd[N:])
    err = np.std(rmsd[N:])
    return avg, err

def mutate(pdb):
    pos = rng.integers(len(list(pdb.topology.residues())))
    residue = list(pdb.topology.residues())[pos]
    targetpdb = openmm.app.PDBFile('trpzip2.pdb')
    targetres = list(targetpdb.topology.residues())[pos].name
    resto = targetres if residue.name=='ALA' else 'ALA'
    mutate = '%s-%s-%s' % (residue.name, residue.id, resto)
    log.info("Mutating %s" % (mutate))
    mutant = copy.deepcopy(pdb)
    mutant.applyMutations([mutate], residue.chain.id)
    mutant.findMissingResidues()
    mutant.findMissingAtoms()
    mutant.addMissingAtoms()
    modeller = openmm.app.Modeller(mutant.topology, mutant.positions)
    modeller.delete(atom for atom in modeller.topology.atoms() if atom.element is openmm.app.element.hydrogen)
    mutant.topology = modeller.topology
    mutant.positions = modeller.positions
    mutant.addMissingHydrogens(7.0)
    return mutant

def p_accept(val1, val2):
    v1, e1 = val1
    v2, e2 = val2
    delta = v2 - v1
    error = math.sqrt(e1**2 + e2**2)
    return 0.5*(1.0-math.erf(delta/(error*math.sqrt(2.0))))

def main():
    #polyala = polyalanine('helix.pdb')
    #polyala = polyalanine('trpzip2.pdb')
    polyala = fixpdb('trpzip2.pdb')
    openmm.app.PDBFile.writeFile(polyala.topology, polyala.positions, open('polyala.pdb', 'w'))

    protein = PDBFixer('polyala.pdb')
    run_md(protein)
    rmsd_old = compute_rmsd(protein)

    log.info("%s  %6.3f +- %6.3f   %6.3f +- %6.3f   %8.6f <? %8.6f   %s" % (sequence(protein), rmsd_old[0], rmsd_old[1], 0.0, 0.0, 0.0, 0.0, ''))

    while True:
        mutant = mutate(protein)
        run_md(mutant)
        rmsd = compute_rmsd(mutant)
        pac = p_accept(rmsd_old, rmsd)
        rfloat = rng.random()
        if rfloat<pac:
            log.info("%s  %6.3f +- %6.3f   %6.3f +- %6.3f   %8.6f <? %8.6f   %s" % (sequence(mutant), rmsd[0], rmsd[1], rmsd_old[0], rmsd_old[1], rfloat, pac, 'Accepted!'))
            protein = copy.deepcopy(mutant)
            rmsd_old = rmsd
        else:
            log.info("%s  %6.3f +- %6.3f   %6.3f +- %6.3f   %8.6f <? %8.6f   %s" % (sequence(mutant), rmsd[0], rmsd[1], rmsd_old[0], rmsd_old[1], rfloat, pac, 'Rejected!'))


main()

