# Imports


```python
from openmm import MonteCarloMembraneBarostat, unit
from openmm.app import CharmmParameterSet, CharmmPsfFile, PME, HBonds
import os, bz2, openmm
```


```python

```

## Setting Variables


```python
psf_path = './inputs/step5_input.psf'
charmm_param_dir = "./inputs/charmm-params/"
state_path = "./inputs/state.xml.bz2"
```


```python
hydrogen_mass = 4.0 * unit.amu
temperature = 303.15 * unit.kelvin
friction = 1 / unit.picoseconds
time_step = 0.002 * unit.picoseconds
pressure = 1 * unit.bar
surface_tension = 0  # units are complicated
```


```python
nonbonded_method = PME
constraints = HBonds
```


```python

```


```python

```

# Loading Files


```python
psf = CharmmPsfFile(psf_path)
```


```python
param_paths = [os.path.join(charmm_param_dir, path) for path in os.listdir(charmm_param_dir)]
params = CharmmParameterSet(*param_paths)
```


```python
with bz2.open(state_path, 'rb') as infile:
    state = openmm.XmlSerializer.deserialize(infile.read().decode())
```


```python
x, y, z = state.getPeriodicBoxVectors()
psf.setBox(x[0], y[1], z[2])
```

# Build System Basics


```python
system = psf.createSystem(params,
                              nonbondedMethod=nonbonded_method,
                              constraints=constraints,
                              removeCMMotion=False,
                              hydrogenMass=hydrogen_mass)

integrator = openmm.LangevinMiddleIntegrator(temperature,
                                             friction,
                                             time_step)

barostat = openmm.MonteCarloMembraneBarostat(pressure,
                                             surface_tension,
                                             temperature,
                                             MonteCarloMembraneBarostat.XYIsotropic,
                                             MonteCarloMembraneBarostat.ZFree
                                             )
barostat.setFrequency(50)

system.addForce(barostat)
```




    8



# Create RMSD Restraint Force


```python
rmsd_cv = openmm.RMSDForce(state.getPositions())
energy_expression = f"(spring_constant/2)*max(0, RMSD-RMSDmax)^2"
# energy_expression = f"RMSD"
restraint_force = openmm.CustomCVForce(energy_expression)
restraint_force.addCollectiveVariable('RMSD', rmsd_cv)
restraint_force.addGlobalParameter('RMSDmax', 0.4)
restraint_force.addGlobalParameter("spring_constant", 0)
# *openmm.unit.kilojoules_per_mole / openmm.unit.nanometers
```




    1




```python
restraint_force.getForceGroup()
```




    0



## set force group of restraint force


```python
force_group = 20
restraint_force.setForceGroup(force_group)
restraint_force.getForceGroup()
```




    20



## Add restraint to system


```python
force_idx = system.addForce(restraint_force)
print(force_idx)
```

    9



```python
system.getForce(9)
```




    <openmm.openmm.CustomCVForce; proxy of <Swig Object of type 'OpenMM::CustomCVForce *' at 0x186aca030> >



# Build Simulation Context


```python
sim = openmm.app.Simulation(psf.topology,
                                system=system,
                                integrator=integrator,
                                )
```

## set state


```python
sim.context.setState(state)
```


```python

```

## Check values of restraint force


```python
state0 = sim.context.getState(getForces=True,
                                 getEnergy=True,
                                 groups=force_group)
```


```python
state0.getPotentialEnergy()
```




    Quantity(value=119478.8203125, unit=kilojoule/mole)




```python
state0.getForces()[0]
```




    Quantity(value=Vec3(x=-37.856407165527344, y=-57.86224365234375, z=-1.0456074476242065), unit=kilojoule/(nanometer*mole))



## check collective variable


```python
system0 = sim.context.getSystem()
force = system0.getForce(9)
```


```python
print(force.getName())
print(force.getCollectiveVariableValues(sim.context))
```

    CustomCVForce
    (0.000806081945833448,)



```python
force.getEnergyFunction()
```




    '(spring_constant/2)*max(0, RMSD-RMSDmax)^2'




```python

```


```python

```
