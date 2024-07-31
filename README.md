# wellmixediobl-jpo2024
last update: 31/07/2024

This repository contains the python scripts that were used to generate the figures in the paper:
 L.-A. Couston, Turbulent ice-ocean boundary layers in the well-mixed regime: insights from direct numerical simulations. Journal of Physical Oceanography, in press (2024)

The simulation data is available on Zenodo at 


- python3 jpo_3d_view.py xxx.npz
with xxx=last_200_1_1_10_1, last_400_1_1_10_1, last_800_1_1_10_1.npz generates plots that can be rearranged to produce Fig. (2) of the paper.
The data sets are numpy arrays of model variables originally saved as .hdf5 checkpoint files by the Dedalus code (final time of the simulation).

- python3 jpo_figures_march2024.py
generates Fig. (3) to (8) of the paper. It uses different data sets:
  - (spatially-)scalar variables stored as numpy arrays in files of the form data_scalars_REV_1_PRV_LEV_1.npz with REV, PRV, and LEV, the values for the Reynolds, Prandtl, and Lewis numbers, respectively (there are as many files as simulations/model configurations).
  - vertically-varying variables (plane averages) stored as numpy arrays in files of the form data_profiles_REV_1_PRV_LEV_1.npz (same notation code as above).
  - variables related to the calculation of the TKE dissipation rates (including plane averages) and stored as numpy arrays in files of the form data_dissipationrate_REV_1_PRV_LEV_1.npz (same notation code as above).

The data_scalars_...npz, data_profiles_...npz, and data_dissipationrate_...npz files can be found in the compressed folders scalars-140624.tar.xz, profiles-140624.tar.xz, and dissipationrate-140624.tar.xz, respectively. 

Some filenames have prefixes that denote simulation results obtained with model configurations slightly different from the default one. Most often, these are results obtained from simulations with different resolutions, described in the Supplementary Material (SM). Specifically:
- SP prefix means that the results were obtained from a full T-S simulation initialized from the final checkpoint of a thermal driving simulation
- data-wide prefix denotes results obtained for the wide domain (last column of Table 3 in the paper)
- data-wide prefix denotes results obtained for the wide domain (see SM)
- LRSP prefix denotes results obtained for an overall lower resolution than default (see SM)
- LYRSP prefix denotes results obtained for a lower spanwise resolution than default (see SM)
- HXRSP prefix denotes results obtained for a higher streamwise resolution than default (see SM)
- HRSP_Z prefix denotes results obtained for a higher wall-normal resolution than default (see SM)
- HRSP prefix denotes results obtained for an overall higher resolution than default (see SM)
