"""
openMM code to calculate osmotic pressure and osmotic coefficients from harmonic potentials

"""

# Imports
import logging
logging.basicConfig(
    level=logging.INFO
)

import warnings
warnings.filterwarnings('ignore')

from typing import Union, Iterable, Optional

import argparse
from argparse import Namespace

import shutil
import pickle
from pathlib import Path

from openmm import CustomExternalForce, System
from openmm.app import Topology as OMMTopology

import openmm.unit
from openmm.unit import kilojoule_per_mole, nanometer
from openmm.unit import Unit, Quantity

from openff.toolkit import Molecule, Topology, ForceField
from openff.interchange import Interchange

from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.genutils.decorators.functional import allow_string_paths

from polymerist.unitutils.interop import openff_to_openmm

from polymerist.mdtools.openfftools import topology
from polymerist.mdtools.openfftools import boxvectors
from polymerist.mdtools.openfftools import TKREGS, FF_PATH_REGISTRY

from polymerist.mdtools.openmmtools.execution import run_simulation_schedule
from polymerist.mdtools.openmmtools.parameters import SimulationParameters


## Initializing paths and working directory
@allow_string_paths
def produce_interchange(pdb_path : Optional[str], working_dir : Path, file_prefix : str, ff_names : Optional[Iterable[str]]=None) -> Interchange:
    '''Load a valid OpenFF Interchange for the given PDB system, creating a new one and writing relevant files if none is found'''
    # boilerplate to initial default values if none are given
    if ff_names is None:
        ff_names = [
            'openff-2.1.0.offxml',
            'tip3p.offxml',
        ]

    # copy PDB into working directory if not already extant
    if pdb_path.parent != working_dir:
        shutil.copy(pdb_path, working_dir)

    sdf_path = assemble_path(working_dir, file_prefix, extension='sdf') 
    inc_path = assemble_path(working_dir, file_prefix, extension='pkl') 

    ## Load topology from PDB
    if sdf_path.exists():
        offtop = topology.topology_from_sdf(sdf_path)
    else:
        assert(pdb_path.exists())
        offtop = Topology.from_pdb(pdb_path)

    ## Parameterize system with OpenFF  
    ### load interchange from file if already extant, otherwise make a new one and save it
    if inc_path.exists():
        with inc_path.open('rb') as file:
            inc  = pickle.load(file)
    else:
        ff = ForceField(*ff_names) # load generic Sage + TIP3P force fields
        inc = ff.create_interchange(topology=offtop, toolkit_registry=TKREGS['OpenEye Toolkit']) # convert to interchange for export
        inc.box = boxvectors.get_topology_bbox(offtop)

        if not sdf_path.exists():
            topology.topology_to_sdf(sdf_path, inc.topology)

        with inc_path.open('wb') as file:
            pickle.dump(inc, file)

    return inc

# Applying restraints - Harmonic Potential
def insert_harmonic_potential_OpenMM(
        openmm_system : System,
        openmm_topology : OMMTopology,
        k : Quantity=0.68095403*(kilojoule_per_mole/nanometer**2),
        delta_z : Optional[Quantity]= None,
        z_center : Quantity=7.2 * nanometer,
        target_residues : Optional[list[str]]=None,
    ) -> None:
    '''Adds a custom harmonic force to an OpenMM System and registers selected particles to that force'''
    if target_residues is None:
        target_residues = ['NA', 'CL']

    fb_force = CustomExternalForce('0.5*k*((z-z0)^2)')
    fb_force.addGlobalParameter('k', k)
    fb_force.addGlobalParameter('z0', z_center)
    openmm_system.addForce(fb_force)

    ## Applying custom external force to only the specified residues
    for atom in openmm_topology.atoms():
        res=atom.residue.name
        if res in target_residues:
            fb_force.addParticle(atom.index, [])

def insert_flat_bottom_potential_OpenMM(
        openmm_system : System,
        openmm_topology : OMMTopology,
        k : Quantity=4184*(kilojoule_per_mole/nanometer**2),
        z_center : Quantity=7.2 * nanometer,
        delta_z : Quantity= 2.4 * nanometer,
        target_residues : Optional[list[str]]=None,
    ) -> None:
    '''Adds a custom flat-bottomed potential force to an OpenMM System and registers selected particles to that force'''
    fb_force = CustomExternalForce('0.5*k*(max(0, abs(z-z0)-rbf)^2)')
    fb_force.addGlobalParameter('k', k)
    fb_force.addGlobalParameter('rbf', delta_z)
    fb_force.addGlobalParameter('z0', z_center)
    openmm_system.addForce(fb_force)

    ## Applying custom external force to only the specified residues
    for atom in openmm_topology.atoms():
        res=atom.residue.name
        if res in target_residues:
            fb_force.addParticle(atom.index, [])

RESTRAINT_TYPES = {
    'FBP' : insert_flat_bottom_potential_OpenMM,
    'HP' : insert_harmonic_potential_OpenMM,
}
RESTRAINT_TYPE_ALIASES = {
    'FBP' : 'flat-bottomed potential',
    'HP' : 'harmonic potential',
}

# read arguments from CLI to pass into OpenMM runner

def parse_args() -> Namespace:
    '''Read user inputs from the command line and preprocess basic parts of them'''
    parser = argparse.ArgumentParser(description='Run OpenMM simulation schedules according to presets of simulation parameters')

    parser.add_argument('-n', '--name', help='Tag to use to refer to working directory and internal files', required=True)
    parser.add_argument('-p', '--postfix', help='Optional extra modifier to "name"; will be appended to name with an underscore separator', default='')
    parser.add_argument('-c', '--cwd', help='The current working directory that all created files/directories should be made relative to', type=Path, default=Path.cwd())
    parser.add_argument('-pdb', '--pdb_path', help='A path to a PDB file containing the initial structure of the system being simulated', type=Path, required=True)
    parser.add_argument('-omm', '--openmm_dir', help='The directory path into which OpenMM setup, checkpoint, and output files should be written', type=Path)
    parser.add_argument('-spp', '--sim_param_paths',
        help='One or more serialized simulation presets; \
            \nsimulations built from these presets will be run in the order they are provided here',
        type=Path, nargs='+', required=True
    )
    parser.add_argument('-r', '--restraint_type', help='The kind of potential restraint to apply to ions in the systems', choices=['HP', 'FBP'], required=True)
    parser.add_argument('-k', help='The spring constant to apply in the harmonic potential restraint') # this should also default to a float, but want to retain the NoneType value to check if this is unset
    parser.add_argument('--z_center', help='', type=float, default=7.2)
    parser.add_argument('--delta_z', help='', type=float, default=2.4)
    parser.add_argument('-du', '--distance_unit', help='The unit convention to use for units of length', choices=['angstrom', 'angstroms', 'nanometer', 'nanometers'], default='nanometer')

    args = parser.parse_args()

    # define directories
    if args.postfix != '':
        args.postfix = '_' + args.postfix
    args.file_prefix = f'{args.name}{args.postfix}'
    args.working_dir = args.cwd / args.file_prefix # Can change this to wherever your files are

    if args.openmm_dir is None:
        args.openmm_dir = args.working_dir / f'{args.file_prefix}_OpenMM'

    # define Quanitity parameters from restraint
    args.restraint_fn = RESTRAINT_TYPES[args.restraint_type]
    logging.info(f'Using restraint type {args.restraint_type} ({RESTRAINT_TYPE_ALIASES[args.restraint_type]})')

    args.distance_unit = getattr(openmm.unit, args.distance_unit) # convert string name of unit into OpenMM object

    args.delta_z = args.delta_z * args.distance_unit
    args.z_center = args.z_center * args.distance_unit

    # slightly more involved handling of force constants, as this default does differ between restraint implementations
    FORCE_CONST_UNIT = (kilojoule_per_mole/nanometer**2)
    if args.k is not None:
        args.k = float(args.k) * FORCE_CONST_UNIT # enforce Quantity with valid units
    else:
        if args.restraint_type == 'HP':
            args.k = 0.68095403 * FORCE_CONST_UNIT
        elif args.restraint_type == 'FBP':
            args.k = 4184 * FORCE_CONST_UNIT
        
    return args


def main():
    '''Run an OpenMM simulation schedule in a given working directory, with arbitrary PDB starting structure and output directory'''
    '''Ensure that a working directory is set up given a name handle that will be used to refer to it
    Returns the resulting path'''
    args = parse_args()
    args.working_dir.mkdir(exist_ok=True)
    args.openmm_dir.mkdir(exist_ok=True)

    # Build simulation schedule based on serialized simulation parameter presets
    schedule = {
        path.stem : SimulationParameters.from_file(path)
            for path in args.sim_param_paths
    }

    # make Interchange and export to OpenMM objects
    inc = produce_interchange(args.pdb_path, working_dir=args.working_dir, file_prefix=args.file_prefix)

    ommtop = inc.to_openmm_topology()
    ommsys = inc.to_openmm_system(combine_nonbonded_forces=False, add_constrained_forces=True)
    ommpos = openff_to_openmm(inc.positions)

    # apply appropriate restraint
    args.restraint_fn( 
        openmm_system=ommsys, 
        openmm_topology=ommtop, 
        k=args.k, 
        delta_z=args.delta_z, 
        z_center=args.z_center
    )

    history = run_simulation_schedule(
        working_dir=args.openmm_dir,
        schedule=schedule,
        init_top=ommtop,
        init_sys=ommsys,
        init_pos=ommpos,
        return_history=True
    )

if __name__ == '__main__':
    main()