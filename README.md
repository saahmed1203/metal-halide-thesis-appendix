# Github for Appendix E and F in Master's thesis
This repository includes code used to conduct my master's thesis research at SFSU. 


## Appendix E. CSV file for ML algorithm

The csv file `FINAL_metal_halide_full_moments.csv` includes the cation, anion, structure, temperature, metal radius/ halide radius, geometric mean of the metal and halide radii, all of the mathematical moments of polarization and bond angle, number of real jumps, number of jump sites, diffusion coefficient (cm^2/s), and activation energy (eV). This CSV file was used in the machine learning algorithm where the features (1st 20 columns, excluding `Rm/Rx` and `geo mean`) were used to predict the targets (either `Diffusion coefficient` or `Ea(eV)`). The number real jumps and number of jump sites were not included in the ML algorithm for they were seen as too insignificant in the algorithm. The number of real jumps was also used in the scatter plot in Figure 11. 


## Appendix F . Code for bond analyses

### Diffusivity calculation
The mean squared displacement of the cations was calculated using `msd.py` in order to calculate cation diffusivity. Licensed by MIT. 

### Jump analysis (NEWcreateSamosFile.py and findjumps.py)
The program SITATOR was used to calcuate cation jumps. `NEWcreateSamosFile.py`found all the jump sites while `findjumps.py` filtered out the real jumps. 

### PDF analysis
The pair distribution functions were created using the executable `RDFpbc_gen_v4.x`

### Bond angle analysis
The bond angles were calculated using `click_AngleAnalysis.py`.

### Polarization analysis
The halide polarizations were calculated using `halide.pol.newNL.py.V2`

### ML algorithm
The Jupyter Notebook `FINAL_Algorithm_for_Ea_and_diffusion.ipynb` includes all of the ML algorithm and bar plot (Figure 15) illustrating which features were most important. The CSV file `FINAL_metal_halide_full_moments.csv` was used for this notebook.

### Density plots
The isosurface density diagrams (Figure 18) were created using `xyz2dens_pbc.x` and `xyz_transform.alpha.x`. `xyz_transform.alpha.x` collapses the supercell into a single unit cell, and `xyz2den_pbc.x` conducts symmetry operations to produce the final isosurface density diagrams. 
