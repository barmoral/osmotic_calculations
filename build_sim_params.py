'''Procedurally build and format simulation parameter presets for polymerist-based OpenMM hooks'''

from pathlib import Path
from openmm.unit import kelvin, atmosphere, nanosecond, picosecond, femtoseconds

from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.mdtools.openmmtools.parameters import SimulationParameters, ThermoParameters, IntegratorParameters, ReporterParameters
from polymerist.mdtools.openmmtools.thermo import EnsembleFactory


all_omm_sims : dict[str, SimulationParameters] = {
    'equil_sim' : SimulationParameters(
        integ_params=IntegratorParameters(
            time_step=2*femtoseconds,
            total_time=2*nanosecond,
            num_samples=100,
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
            num_samples=1001,
        ),
        thermo_params=ThermoParameters(
            ensemble='NVT',
            temperature=300 * kelvin,
        ),
        reporter_params=ReporterParameters(),
    ),
    'test_sim' : SimulationParameters(
        integ_params=IntegratorParameters(
            time_step=2*femtoseconds,
            total_time=2*picosecond,
            num_samples=10,
        ),
        thermo_params=ThermoParameters(
            ensemble='NVT',
            temperature=600 * kelvin,
        ),
        reporter_params=ReporterParameters(),
    ),
}  

if __name__ == '__main__':
    sim_dir = Path('./sim_param_sets')
    sim_dir.mkdir(exist_ok=True)

    for step_name, sim_params in all_omm_sims.items():
        sp_path = assemble_path(sim_dir, step_name, extension='json')
        sim_params.to_file(sp_path)