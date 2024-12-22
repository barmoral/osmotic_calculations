Theory for the Calculation of Osmotic Values
================

Early molecular dynamics methods for osmotic pressure, introduced by Murad and Powles [1], utilized semipermeable membranes within the simulation to separate solvent and solution regions. Luo and Roux [2] refined this approach by employing flat-bottom potentials as virtual semi-permeable walls to simplify the process and avoid complications from direct membrane modeling. These techniques require multiple simulations to explore the concentration dependence of the osmotic equation-of-state. To streamline this, Milner et al. [3] developed the osmotic force balance method, applying an external harmonic potential to solutes, thereby enabling easier assessment of osmotic pressure gradients. They successfully applied this to study substances like NaCl and mixtures like benzene and pyridine. Hosseini et al. [4] expanded on this in 2023, calculating chemical potentials for 15 alkali halide salts, thereby providing a more efficient method to determine osmotic pressures and propose Lennard-Jones parameters for these salts. I have replicated the calculations of osmotic pressure with respect to concentration using both techniques.

### <ins>Osmotic Pressure calculation using flat-bottom potentials:</ins>

Simulating virtual semi-permeable walls using flat-bottomed potentials (figure 1) allows for the calculation of osmotic pressure for a single concentration. The mean force of the system's walls exerted on the ions is directly related to the osmotic pressure as shown in *equation 1*, where Fwall is the instantaneous force exerted by one wall onto the ions, k is the force constant, N is the number of steps, zi is the position in the z axis of the ions located outside of both the walls, and zwall is the position of the wall. The mean force is averaged between the two half-harmonic walls. The osmotic pressure was calculated as shown in *equation 2*, where A is the cross-sectional area of the simulation box, and with it, the osmotic coefficients are calculated as a ratio between the observed and ideal results.

$$\langle F_{wall} \rangle = \frac{k}{N} \sum_{N} \sum_{i} (|z_{i} - z_{wall}|)$$  &nbsp;&nbsp;*(1)* 	
	
$$\Pi_{observed}=\langle F_{wall} \rangle / A$$  &nbsp;&nbsp;*(2)*		

![Figure 1: 1m NaCl simulation with outlined flat-bottom potential (red line).](https://github.com/barmoral/osmotic_calculations/blob/main/FBP_model.png)
*Figure 1: 1m NaCl simulation with outlined flat-bottom potential (red line).*

### <ins>Osmotic pressure calculation using a harmonic potential:</ins>

A harmonic potential allows for the calculation of osmotic pressure for a nonuniform equilibrium concentration profile for ions (figure 2). The potential curve used in this method uses the form of   , where k is the spring constant and z is the z coordinate centered in the box. In equilibrium, slices of a concentration profile are regarded as stationary under the sum of three forces: the gradient of external potential (FU, *eq. 3*), and the osmotic pressure (FΠ, *eq. 4*) on the two sides of the slice, which will differ slightly because of the nonuniform concentration. C(z) is the total concentration at height z, and Δz is the length of the slice.

$$F_{U} = C(z) \frac{dU(z)}{dz} \Delta z $$  &nbsp;&nbsp;*(3)*		

$$F_{\Pi} = \frac{d\Pi(z)}{dz} \Delta z$$   &nbsp;&nbsp;*(4)*

When the system is in equilibrium these forces sum to 0. The change in osmotic pressure across a slice is proportional to the change in the harmonic potential multiplied by the concentration of the slice, as seen in *equation 5*. One can integrate the force balance equation after plugging in it the harmonic potential, as seen in *equation 6*, and with it solve for the osmotic pressure at the coordinates of both sides of a slice by averaging over their corresponding concentration.

$$\frac{d\Pi(z)}{dz}= -\nu C(z) \frac{dU(z)}{dz}$$  &nbsp;&nbsp;*(5)*
	
$$\Pi(z)= - \nu \int_{\infty}^{z} C(z) \frac{dU(z)}{dz} dz$$   &nbsp;&nbsp;*(6)*


![Figure 2: 3m NaCl simulation with outlines harmonic potential (red line) and concentration gradient (orange line).](https://github.com/barmoral/osmotic_calculations/blob/main/HP_model.png)
*Figure 2: 3m NaCl simulation with outlines harmonic potential (red line) and concentration gradient (orange line).*


### References:
1. S. Murad and J. G. Powles, “Computer simulation of osmosis and reverse osmosis in solutions,” Chemical Physics Letters, vol. 225, no. 4, pp. 437–440, 1994, doi: 10.1016/0009-2614(94)87108-6.
2. Y. Luo and B. Roux, “Simulation of Osmotic Pressure in Concentrated Aqueous Salt Solutions,” J. Phys. Chem. Lett., vol. 1, no. 1, pp. 183–189, 2010, doi: 10.1021/jz900079w.
3. C. Gillespie and S. T. Milner, “Using osmotic pressure simulations to test potentials for ions,” Soft Matter, vol. 16, no. 42, pp. 9816–9821, 2020, doi: 10.1039/D0SM00957A.
4. A. Hosseini and H. S. Ashbaugh, “Osmotic Force Balance Evaluation of Aqueous Electrolyte Osmotic Pressures and Chemical Potentials,” J. Chem. Theory Comput., vol. 19, no. 23, pp. 8826–8838, 2023, doi: 10.1021/acs.jctc.3c00982.
