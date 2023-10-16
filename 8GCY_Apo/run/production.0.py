import time
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, exit, stderr
from parmed import unit as u
from copy import deepcopy


########### Functions ############

# Function to add backbone position restraints
def add_backbone_posres(system, positions, atoms, restraint_force):
  force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
  force_amount = restraint_force * kilocalories_per_mole/angstroms**2
  force.addGlobalParameter("k", force_amount)
  force.addPerParticleParameter("x0")
  force.addPerParticleParameter("y0")
  force.addPerParticleParameter("z0")
  for i, (atom_crd, atom) in enumerate(zip(positions, atoms)):
    if atom.name in  ('CA', 'C', 'N', 'O'):
      force.addParticle(i, atom_crd.value_in_unit(nanometers))
  posres_sys = deepcopy(system)
  posres_sys.addForce(force)
  return posres_sys

############# Set up #################

### Load structure/psf
psf = CharmmPsfFile('../build/apo.psf')
pdb = PDBFile('../build/apo.pdb')

### Parameter files
params = CharmmParameterSet('../../toppar/top_all36_prot.rtf',
                            '../../toppar/par_all36m_prot.prm',
                            '../../toppar/toppar_water_ions.str')

### Set pbc box
coords = pdb.positions
min_crds = [coords[0][0], coords[0][1], coords[0][2]]
max_crds = [coords[0][0], coords[0][1], coords[0][2]]

for coord in coords:
    min_crds[0] = min(min_crds[0], coord[0])
    min_crds[1] = min(min_crds[1], coord[1])
    min_crds[2] = min(min_crds[2], coord[2])
    max_crds[0] = max(max_crds[0], coord[0])
    max_crds[1] = max(max_crds[1], coord[1])
    max_crds[2] = max(max_crds[2], coord[2])

psf.setBox(max_crds[0]-min_crds[0],
                 max_crds[1]-min_crds[1],
                 max_crds[2]-min_crds[2],
            )

### Create System
system = psf.createSystem(params, nonbondedMethod=app.PME,
                                removeCMMotion=True,
                                nonbondedCutoff=12.0*u.angstroms,
                                constraints=app.HBonds,
                                switchDistance=10.0*u.angstroms,
                                hydrogenMass=4*amu
                        )

### Add barostat with 1 atm, 300 K
barostat = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))


############# Simulation parameters ##############

### Simulation time
mdsteps = 250000000 # Run for 1 us

### Integrator
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
integrator.setConstraintTolerance(0.00001)

### Device information
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': f'0,1,2,3', 'Precision': 'mixed'}


############# Restart simulation ################

### Set simulation and load equil checkpoint
simulation = Simulation(psf.topology, system, integrator, platform, properties)
simulation.loadCheckpoint('eq.chk')

### set positions and velocities from previous sim
eq_state = simulation.context.getState(getVelocities=True, getPositions=True)
positions = eq_state.getPositions()
velocities = eq_state.getVelocities()

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)


############ Add reporters #################

### DCD reporter
simulation.reporters.append(
    DCDReporter('output/prod.0.dcd', 5000))

### Data reporter
simulation.reporters.append(
    StateDataReporter(
         'output/prod.0.csv',
         10000,
         step=True,
         potentialEnergy=True,
         temperature=True,
         progress=True,
         remainingTime=True,
        speed=True,
        volume=True,
        totalSteps=mdsteps,
        separator='\t'
        )
    )

### Checkpointer
simulation.reporters.append(
    CheckpointReporter(
        'output/prod.0.restart.chk',
        100000
        )
    )

############# Run simulation! #############

print('Running Production...')
simulation.step(mdsteps)
simulation.saveState('output/prod.0.state')
simulation.saveCheckpoint('output/prod.0.chk')

