# Boundary Conditions in Supersonic Airfoil Simulations
## Reproducibility Repository

In this repo, we have the full code required to run all the simulations
mentioned in the text _Boundary Conditions in Supersonic Airfoil Simulations_.

Below are the mappings of code files to their references in the text:

- `airfoil_1_mach10_crash.jl`: Corresponds to Equation (3) and Listing 2 in Section 3.2.1
- `airfoil_2_mach9_steady.jl`: Mentioned in Section 3.2.1
- `airfoil_3_rampup.jl`: Discussed in Listing 3, used in Figures 4 and 5 in Section 3.2.3
- `airfoil_4_dirichlet.jl`: Listing 4, Figure 6 in Section 3.2.4
- `airfoil_5_dirichlet_high_resolution.jl`: Listing 4, Figure 7 in Section 3.2.4
- `NACA6412airfoil.geo`: A geometry file used by all code samples. Also used to create Figure 1

This code was adapted from the code _FVM Example: Euler-Eqns with Trixi.jl_ by Manuel Torrilhon (`torrilhon@rwth-aachen.de`) for the class 'Numerical Methods for Partial Differential Equations' at RWTH Aachen. 