# Orbital Integrator Evaluation Code

This repository contains the MATLAB code used to evaluate different numerical integrators for space situational awareness applications. The code compares various integration methods including RK4, RK8, ODE45, ODE78, ODE113, and Modified Picard Chebyshev Iterations (MPCI) under both Keplerian and perturbed orbital conditions.

## MATLAB

### MATLAB Installation

1. Clone or download this repository to your local machine
2. Add all subfolders to your MATLAB path:


### Main Test Files

The primary test files are:

- `main_testODEJ2.m`: Tests integrator performance with J2 and atmospheric drag perturbations
- `main_testODEKep.m`: Tests integrator performance under Keplerian motion (two-body problem)
- `main_GenerateGraphs_cpp`: Use the saved csv from the c++ implementation to check the algorithms' accuracy


### Dependencies

- MATLAB (tested on R2022b)
- All submodules given here (matlab2tikz, table2latex)


### Output

The test scripts will generate runtime and accuracy comparisons between different integrators. Runtime results are saved as latex tables, the figures are saved as pdf and tikz files.

## C++

### C++ Installation
You can use `\.install.sh` given in the folder or use cmake to build.

The primary file is `Main_to_print_executable_time.cpp`

### Dependencies

The build requires a compiler (such as g++ or gcc), cmake and the eigen library. If on linux or mac systems one can install all prerequisits with
`sudo apt install build-essential cmake libeigen-dev`

Running the executable from `./build/Main_to_print_executable_time` would create the csv files in the `results` folder. In addition the results folder contains a handy python script to incorporate different csv files into tables that have later been used in the paper.


## Citation

If you use this code in your research, please cite:
```
AbedAlkream, S., AbedEllatif, M., Peretz, G., & Denenberg, E. 
_Evaluating Integrators for Autonomous Space Situational Awareness Module._ 
Astrodynamics. 2025
```