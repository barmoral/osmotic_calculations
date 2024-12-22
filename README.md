Calculation and Analysis of Osmotic Values
================

Calculations can be performed with two different types of restraints: [Harmonic Potential (HP) and Flat-Bottom Potentials (FBPs)](osmotic_theory_summary.md).

### Simulation Details
* Python codes were written for simulation boxes with dimension of 4.8nm x 4.8nm x 14.4 nm.
* The number of ions in each simulation changes depending on the desired concentration in the box.
* The number of water molecules stays fixed at 11100.
* The sample simulations performed in this repository were equilibrated for 2 ns and their production run lasted 20 ns.

* Input coordinate files are generated using packmol. The packmol inputs and outputs can be found on the directory called ['structures'](structures/).

* The force constant (k, with units of kJ/mol/nm^2) used for HP and FBP restraints is different. While the FBP method is not very sensitive to the force constant value (in here I use the same value as the method developers, [Luo & Roux](https://pubs.acs.org/doi/10.1021/jz900079w)), HP is highly sensitive, and this value has to be estimated beforehand using the python notebook called ['force_constant_hps.ipynb'](https://github.com/barmoral/osmotic_calculations/blob/main/force_constant_hps.ipynb).

* I included example script files to use in a cluster called [run_HP_openmm.sh](https://github.com/barmoral/osmotic_calculations/blob/main/run_HP_openmm.sh) and [run_FBP_openmm.sh](https://github.com/barmoral/osmotic_calculations/blob/main/run_FBP_openmm.sh), which can also give an idea of how to run the calculations from any CLI.


### Steps to run Harmonic Potential calculations and analysis:
1. Define maximum concentration desired and corresponding number of ions needed in the given volume. (I used 3.5 molal as maximum concentration, which requires 217 NaCl molecules (217 of each ion)). 
2. Estimate the force constant using notebook [force_constant_hps.ipynb](https://github.com/barmoral/osmotic_calculations/blob/main/force_constant_hps.ipynb). (0.68095403 kJ/mol/nm^2 was the result I obtained)
3. Run in cluster using script file [run_HP_openmm.sh](https://github.com/barmoral/osmotic_calculations/blob/main/run_HP_openmm.sh). (I indicated the variables that need to be modified. Make sure the variable 'RESTRAINT' is set to HP)
4. The script file will generate the [sim_param_sets](https://github.com/barmoral/osmotic_calculations/tree/main/sim_param_sets) directory, which contains the json files with simulation scheduling instructions. If you want to modify the equilibration or production run setup, it must be done in the code [build_sim_params.py](https://github.com/barmoral/osmotic_calculations/blob/main/build_sim_params.py) directly before running the script file.
5. The script file will then run the calculation code, [osmotic_sim_dispatch.py](https://github.com/barmoral/osmotic_calculations/blob/main/osmotic_sim_dispatch.py), which will generate two directories: 1) holds a copy of the input .pdb file, and generated .sdf and .pkl files. 2) holds all simulation outputs (labeled with postfix '_sims' at the end).
6. Get the .dcd trajectory files from the 'prod_sim' subfolder of each calculation, and use the function mdconvert in the CLI from the mdtraj package to convert these to .xtc trajectory files.
7. Get the .pdb coordinate files from the 'prod_sim' subfolder of each calculation, and copy these along with their corresponding .xtc trajectory files into the [osmotic_calculations_trajs/HP](https://github.com/barmoral/osmotic_calculations/tree/main/osmotic_calculations_trajs/HP) directory, or directly into the HarmonicPotential_analysis directory. (Just make sure to specify the desired working directory)
8. Run the analysis notebook called [harmonic_maxlikelihood_analysis.ipynb](https://github.com/barmoral/osmotic_calculations/blob/main/HarmonicPotential_analysis/harmonic_maxlikelihood_analysis.ipynb), also located in [HarmonicPotential_analysis directory](https://github.com/barmoral/osmotic_calculations/tree/main/HarmonicPotential_analysis). (Note: In this code, a bootstrapping calculation will be a bottleneck, taking ~11 minutes to run if you use 500 bootstraps).

### Steps to run Flat-Bottom Potentials calculations and analysis:
1. Define which concentrations are wanted and their corresponding number of ions for the given volume. (I used 1 molal, 2 molal, and 3 molal, which respectively required 65, 128, and 188 NaCl molecules)
2. Define desired force constant. (I use 4184 kJ/mol/nm^2, which is the one [Luo & Roux](https://pubs.acs.org/doi/10.1021/jz900079w) recommended in their 2010 paper I was attempting to replicate)
3. Run in cluster using script file [run_FBP_openmm.sh](https://github.com/barmoral/osmotic_calculations/blob/main/run_FBP_openmm.sh). (I indicated the variables that need to be modified. Make sure the variable 'RESTRAINT' is set to FBP)
4. The script file will generate the [sim_param_sets](https://github.com/barmoral/osmotic_calculations/tree/main/sim_param_sets) directory, which contains the json files with simulation scheduling instructions. If you want to modify the equilibration or production run setup, it must be done in the code [build_sim_params.py](https://github.com/barmoral/osmotic_calculations/blob/main/build_sim_params.py) directly before running the script file.
5. The script file will then run the calculation code, [osmotic_sim_dispatch.py](https://github.com/barmoral/osmotic_calculations/blob/main/osmotic_sim_dispatch.py), which will generate two directories: 1) holds a copy of the input .pdb file, and generated .sdf and .pkl files. 2) holds all simulation outputs (labeled with postfix '_sims' at the end).
6. Get the .dcd trajectory files from the 'prod_sim' subfolder of each calculation, and use the function mdconvert in the CLI from the mdtraj package to convert these to .xtc trajectory files.
7. Get the .pdb coordinate files from the 'prod_sim' subfolder of each calculation, and copy these along with their corresponding .xtc trajectory files into the [osmotic_calculations_trajs/FBPs](https://github.com/barmoral/osmotic_calculations/tree/main/osmotic_calculations_trajs/FBPs) directory, or directly into the [FlatBottomPotentials_analysis](https://github.com/barmoral/osmotic_calculations/tree/main/FlatBottomPotentials_analysis) directory. (Just make sure to specify the desired working directory)
8. Run the analysis notebook called [FBP_analysis.ipynb](https://github.com/barmoral/osmotic_calculations/blob/main/FlatBottomPotentials_analysis/FBP_analysis.ipynb), also located in FlatBottomPotentials_analysis directory.


### Code to compare FBPs and HP results with experiments:
The directory called FBPs_vs_HP_vs_Experiments merges the results for both calculation types and plots them along with experiments to be able to compare all of them.