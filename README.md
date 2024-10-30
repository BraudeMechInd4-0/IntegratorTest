# Orbital Integrator Evaluation Code

This repository contains the MATLAB code used to evaluate different numerical integrators for space situational awareness applications. The code compares various integration methods including RK4, RK8, ODE45, ODE78, ODE113, and Modified Picard Chebyshev Iterations (MPCI) under both Keplerian and perturbed orbital conditions.

## Installation

1. Clone or download this repository to your local machine
2. Add all subfolders to your MATLAB path:


## Main Test Files

The primary test files are:

- `testODEJ2.m`: Tests integrator performance with J2 and atmospheric drag perturbations
- `testODEKep.m`: Tests integrator performance under Keplerian motion (two-body problem)


## Dependencies

- MATLAB (tested on R2022b)
- All submodules given here (matlab2tikz, table2latex)


## Output

The test scripts will generate runtime and accuracy comparisons between different integrators. Runtime results are saved as latex tables, the figures are saved as pdf and tikz files.

## Citation

If you use this code in your research, please cite:
```
AbedAlkream, S., AbedEllatif, M., Peretz, G., & Denenberg, E. (2024). 
Evaluating Integrators for Autonomous Space Situational Awareness Module. 
Acta Astronautica.
```