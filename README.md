## ChemistrySolver
### Nine Species Chemistry Simulator

This is a Nine Species Chemistry simulator for ASTR498 - Computational Astrophysics.
<br><br>
#### Solver.py  
The main simulator. This will integrate until the specified time and produce an output "simulation.png" showing how the values change for each species (and Temperature).

Add the "-v" flag to specify initial values, otherwise these default values will be used:
* n_total = 5
* h_ionized_frac = -6
* he_ionized_frac = -5
* h_mol_ionized_frac = -2
* T = 150000
* final_t = 1000000
* safety_factor = 100000
<br><br>
#### CalcReactionRates.py
Helper functions for Solver.py to compute each of the reaction rates for a given temperature.

#### RHSGenerator.py
A right hand side generator for the nine species of H and He, based on reactions 1-19, 57, and 58 (numbering from with Abel et al. (1997)).
