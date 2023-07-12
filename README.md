# Github for Appendix E and F in Master's thesis
This repository includes code used to conduct my master's thesis research at SFSU. 


## Appendix E. CSV file for ML algorithm

The csv file `FINAL_metal_halide_full_moments.csv` includes the cation, anion, structure, temperature, metal radius/ halide radius, geometric mean of the metal and halide radii, all of the mathematical moments of polarization and bond angle, number of real jumps, number of jump sites, diffusion coefficient (cm^2/s), and activation energy (eV). This CSV file was used in the machine learning algorithm where the features (1st 20 columns, excluding `Rm/Rx` and `geo mean`) were used to predict the targets (either `Diffusion coefficient` or `Ea(eV)`). The number real jumps and number of jump sites were not included in the ML algorithm for they were seen as too insignificant in the algorithm. The number of real jumps was also used in the scatter plot in Figure 11. 


## Appendix F 
