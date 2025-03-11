# symmetry-FEM-for-optical-nonsymmorphic-problem
nonsymmorphic symmetry FEM for light band structure problem
![image](Figure/PCG/png)
![image](Figure/PCG/TQTPHC)
![image](Figure/PCG/NLS)
## Overview
The Open-source MATLAB package nonsymmorphic symmetry adapted FEM (NSA-FEM) is an optical finite element program used for band structure problems. Based on the symmetry in group theory, it can divide the original problem into several sub-tasks and truncate the computational domain. This program significantly enhances the computational efficiency for band structure analysis problems. We provided examples of photonic crystal gratings with plane group pg, twisted quadrupole topological photonic crystals with plane group pg and photonic nodal line space group, respectively.
## Usage
### function for matrix assembly
1. `Assembel.m`:Used for standard FEM matrix assembly with Irreps of the little group at k. Output Ai,Aj,Av,Bi,Bj,Bv for the assembly of sparse matrices A and B.
2. `Assembel2.m`:Used for matrix assembly with full-group Irreps induced at k∗. Output A2i,A2j,A2v,B2i,B2j,B2v for the assembly of sparse matrices A2 and B2.
### Main function of FEM
1. `Bandgap.m`: Main function of the NSA-FEM to compute the band structure of the model in nonsymmorphic symmetry configuration.
### Run examples
Run the `Assembel.m`(`Assembel2.m`) to calculate sparse matrices A and B. Run the `Bandgap.m`to calculate the band structure.
