"""
openMM code to calculate osmotic pressure and osmotic coefficients from flat-bottom potentials

"""

# Imports
import logging
logging.basicConfig(
    level=logging.INFO
)

import polymerist as ps
from polymerist.genutils.importutils import module_hierarchy


import warnings
warnings.filterwarnings('ignore')
import os
from pathlib import Path
from polymerist.genutils.fileutils.pathutils import assemble_path


# Parameterizing PDB system

## Initializing paths and working directory
cwd=os.getcwd()
working_dir = Path(cwd+'/SPW_1m_r1') # can change this to wherever your files are
mol_name = 'SPW_1m_r1'

pdb_path = assemble_path(working_dir, mol_name, extension='pdb') 
sdf_path = assemble_path(working_dir, mol_name, extension='sdf') 
inc_path = assemble_path(working_dir, mol_name, extension='pkl') 

## Load topology from PDB
from openff.toolkit import Molecule, Topology, ForceField
from polymerist.mdtools.openfftools import topology

if sdf_path.exists():
    offtop = topology.topology_from_sdf(sdf_path)
else:
    assert(pdb_path.exists())
    offtop = Topology.from_pdb(pdb_path)

## Parameterize system with OpenFF
import pickle
from openff.interchange import Interchange
from polymerist.unitutils.interop import openff_to_openmm
from polymerist.mdtools.openfftools import boxvectors
from polymerist.mdtools.openfftools import TKREGS, FF_PATH_REGISTRY

### load interchange from file if already extant, otherwise make a new one and save it
if inc_path.exists():
    with inc_path.open('rb') as file:
        inc  = pickle.load(file)
else:
    ff = ForceField('openff-2.1.0.offxml', 'tip3p.offxml') # load generic Sage + TIP3P force fields
    inc = ff.create_interchange(topology=offtop, toolkit_registry=TKREGS['OpenEye Toolkit']) # convert to interchange for export
    inc.box = boxvectors.get_topology_bbox(offtop)

    if not sdf_path.exists():
        topology.topology_to_sdf(sdf_path, inc.topology)

    with inc_path.open('wb') as file:
        pickle.dump(inc, file)

### extract OpenMM-specific objects to initialize simulations later
ommtop = inc.to_openmm_topology()
ommsys = inc.to_openmm_system(combine_nonbonded_forces=False, add_constrained_forces=True)
ommpos = openff_to_openmm(inc.positions)


# Applying restraints - Flat-bottom Potentials

from openmm.unit import kilojoule_per_mole, nanometer
from openmm import CustomExternalForce

k = 4184*(kilojoule_per_mole/nanometer**2)
delta_Z = 2.4 * nanometer
z_center = 7.2 * nanometer

fb_force = CustomExternalForce('0.5*k*(max(0, abs(z-z0)-rbf)^2)')
fb_force.addGlobalParameter('k', k)
fb_force.addGlobalParameter('rbf', delta_Z)
fb_force.addGlobalParameter('z0', z_center)
ommsys.addForce(fb_force)

## Applying costum external force to only ions
from openmm.app.topology import Topology

for atom in ommtop.atoms():
    res=atom.residue.name
    if (res == 'NA') or (res == 'CL'):
        fb_force.addParticle(
            atom.index,
            []# if (res == 'NA') or (res == 'CL') else []
        )


# Creating OpenMM Simulation

## Defining simulation parameters
from polymerist.mdtools.openmmtools.parameters import SimulationParameters, ThermoParameters, IntegratorParameters, ReporterParameters
from polymerist.mdtools.openmmtools.thermo import EnsembleFactory

from openmm.unit import kelvin, atmosphere, nanosecond, picosecond, femtoseconds


all_omm_sims : dict[str, SimulationParameters] = {
    'equil_sim' : SimulationParameters(
        integ_params=IntegratorParameters(
            time_step=2*femtoseconds,
            # total_time=100*femtoseconds,
            total_time=2*nanosecond,
            num_samples=50,
        ),
        thermo_params=ThermoParameters(
            ensemble='NVT',
            temperature=300 * kelvin,
        ),
        reporter_params=ReporterParameters(
            traj_ext='dcd',
        ),
    ),
    'prod_sim' : SimulationParameters(
        integ_params=IntegratorParameters(
            time_step=2*femtoseconds,
            total_time=20*nanosecond,
            num_samples=250,
        ),
        thermo_params=ThermoParameters(
            ensemble='NVT',
            temperature=300 * kelvin,
        ),
        reporter_params=ReporterParameters(),
    ),
}    

## Run defined simulation in series

from polymerist.mdtools.openmmtools.execution import run_simulation_schedule

omm_sim_dir = Path('SPW_1m_r1/OpenMM_FBPs_r1')
omm_sim_dir.mkdir(exist_ok=True)

history = run_simulation_schedule(
    working_dir=omm_sim_dir,
    schedule=all_omm_sims,
    init_top=ommtop,
    init_sys=ommsys,
    init_pos=ommpos,
    return_history=True
)