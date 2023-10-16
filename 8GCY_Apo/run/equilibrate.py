import time
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, exit, stderr
from parmed import unit as u
from copy import deepcopy


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

# Set up 
# Load structure/psf
psf = CharmmPsfFile('../build/apo.psf')
pdb = PDBFile('../build/apo.pdb')

print(psf)
print(pdb)

# Parameter files
params = CharmmParameterSet('../../toppar/top_all36_prot.rtf',
                            '../../toppar/par_all36m_prot.prm',
                            '../../toppar/toppar_water_ions.str')


# Set pbc box
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

#System preperation 

system = psf.createSystem(params, nonbondedMethod=app.PME,
                                removeCMMotion=True,
                                nonbondedCutoff=12.0*u.angstroms,
                                constraints=app.HBonds,
                                switchDistance=10.0*u.angstroms,
                                hydrogenMass=4*amu
                        )

posres_sys = add_backbone_posres(system, pdb.positions, psf.topology.atoms(), 42)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
integrator.setConstraintTolerance(0.00001)
#print(Platform.getName())
platform = Platform.getPlatformByName('CUDA')
print(platform)
properties = {'DeviceIndex': f'0,1,2,3', 'Precision': 'mixed'}

simulation = Simulation(psf.topology, posres_sys, integrator, platform, properties)
print(simulation)
#simulation = Simulation(psf.topology, posres_sys, integrator, properties)

simulation.context.setPositions(pdb.positions)
print(simulation)
simulation.reporters.append(
    StateDataReporter(
        'equilibrate.log',
        5000,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True
    )
)


#Minimize
print('Minimizing...')
#simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
simulation.step(5000)

#WarmUp with a NVT run.  Slowly warm up temperature - every 1000 steps raise the temperature by 5 K, ending at 300 K
simulation.context.setVelocitiesToTemperature(5*kelvin)
print('Warming up the system...')
T = 5
mdsteps = 50000
for i in range(60):
  simulation.step(int(mdsteps/60) )
  temperature = (T+(i*T))*kelvin 
  integrator.setTemperature(temperature)


#NPT equilibration, reducing backbone constraints
mdsteps = 500000
barostat = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))
simulation.context.reinitialize(True)
print('Running NPT equilibration...')
start_npt = time.time()
for i in range(100):
  simulation.step(int(mdsteps/100))
  simulation.context.setParameter('k', (float(99.02-(i*0.98))*kilocalories_per_mole/angstroms**2))

simulation.context.setParameter('k', 0)
end_npt = time.time()

print(f"ns/day: {2*3600*24/(end_npt - start_npt)}")

# save the equilibration results to file : state is platform independent but less precise, checkpoint file
simulation.saveState('output/eq.state')
simulation.saveCheckpoint('output/eq.chk')
