## ChemistrySolver
### Nine Species Chemistry Simulator

This is a Nine Species Chemistry simulator for ASTR498 - Computational Astrophysics.

RHSGenerator.py - a right hand side generator for the nine species of H and He, based on reactions 1-19, 57, and 58 (numbering from with Abel et al. (1997)).

Solver.py - a simulator that takes H, He, and H2 ionization fractions, a time to evolve to, and density. This will integrate until the specified time and produce an output "simulation.png" showing how the values change for each species (and Temperature).
