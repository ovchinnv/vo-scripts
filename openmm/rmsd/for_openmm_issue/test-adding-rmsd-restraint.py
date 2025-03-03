#!/usr/bin/env python
# coding: utf-8

# # Imports

# In[132]:


from openmm import MonteCarloMembraneBarostat, unit
from openmm.app import CharmmParameterSet, CharmmPsfFile, PME, HBonds
import os, bz2, openmm


# In[ ]:





# ## Setting Variables

# In[140]:


psf_path = './inputs/step5_input.psf'
charmm_param_dir = "./inputs/charmm-params/"
state_path = "./inputs/state.xml.bz2"


# In[141]:


hydrogen_mass = 4.0 * unit.amu
temperature = 303.15 * unit.kelvin
friction = 1 / unit.picoseconds
time_step = 0.002 * unit.picoseconds
pressure = 1 * unit.bar
surface_tension = 0  # units are complicated


# In[142]:


nonbonded_method = PME
constraints = HBonds


# In[ ]:





# In[ ]:





# # Loading Files

# In[143]:


psf = CharmmPsfFile(psf_path)


# In[144]:


param_paths = [os.path.join(charmm_param_dir, path) for path in os.listdir(charmm_param_dir)]
params = CharmmParameterSet(*param_paths)


# In[145]:


with bz2.open(state_path, 'rb') as infile:
    state = openmm.XmlSerializer.deserialize(infile.read().decode())


# In[146]:


x, y, z = state.getPeriodicBoxVectors()
psf.setBox(x[0], y[1], z[2])


# # Build System Basics

# In[162]:


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


# # Create RMSD Restraint Force

# In[163]:


rmsd_cv = openmm.RMSDForce(state.getPositions())
energy_expression = f"(spring_constant/2)*max(0, RMSD-RMSDmax)^2"
# energy_expression = f"RMSD"
restraint_force = openmm.CustomCVForce(energy_expression)
restraint_force.addCollectiveVariable('RMSD', rmsd_cv)
restraint_force.addGlobalParameter('RMSDmax', 0.4)
restraint_force.addGlobalParameter("spring_constant", 0)
# *openmm.unit.kilojoules_per_mole / openmm.unit.nanometers


# In[164]:


restraint_force.getForceGroup()


# ## set force group of restraint force

# In[165]:


force_group = 20
restraint_force.setForceGroup(force_group)
restraint_force.getForceGroup()


# ## Add restraint to system

# In[166]:


force_idx = system.addForce(restraint_force)
print(force_idx)


# In[167]:


system.getForce(9)


# # Build Simulation Context

# In[168]:


sim = openmm.app.Simulation(psf.topology,
                                system=system,
                                integrator=integrator,
                                )


# ## set state

# In[169]:


sim.context.setState(state)


# In[ ]:





# ## Check values of restraint force

# In[170]:


state0 = sim.context.getState(getForces=True,
                                 getEnergy=True,
                                 groups=force_group)


# In[171]:


state0.getPotentialEnergy()


# In[172]:


state0.getForces()[0]


# ## check collective variable

# In[173]:


system0 = sim.context.getSystem()
force = system0.getForce(9)


# In[174]:


print(force.getName())
print(force.getCollectiveVariableValues(sim.context))


# In[175]:


force.getEnergyFunction()


# In[ ]:





# In[ ]:




