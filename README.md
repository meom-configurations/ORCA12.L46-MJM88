# ORCA12.L46-MJM88
ORCA12.L46-MJM88 is the interannual counter part of the climatological run ORCA12.L46-GJM02 performed in the frame of GENCI Grands Défis on IDRIS ADA super computer.

This repository hold the code, and configuration files usefull for running this configuration.

## REFERENCE CODE:
### NEMO 
   DCM revision 1278, which corresponds to NEMO version 3.4. As the exact revision of NEMO on the IPSL forge is not sure, The full reference NEMO code used in this config is reproduce in the repository under NEMO directory.

### XIOS
  No XIOS used for this run (not ready when running this experiment). For the output, the key_dimgout was used 

## BRIEF DESCRIPTION:
### Overview
  This global configuration is the interannual simulation corresponding to ORCA12.L46-GJM02 climatological simulation. The code is exactly the same than ORCA12.L46-GJM02 except for the atmospheric forcing. In this simulation we use the standard NEMO sbcblk_core module, instead of the modified version used in GJM02.  The run covers the period 1958-2012 and has been performed on jade CINES super computer. It has been performed during summer 2013, from June to September.

### Parameterizations:
 1. use filtered free surface (key_dynspg_flt)
 2. use vector form advection scheme for dynamics, with EEN vorticity component.
 3. use TVD advection scheme for tracers.
 4. use biharmonic viscosity scaled with the cube of the mesh size.
 5. use laplacian isopycnal diffusivity for tracers.
 6. use TKE vertical mixing parameterization with enhanced vertical diffusion for deep convection. use tidal mixing parameterization.
 7. use LIM2 ice model, VP rheology.
 8. use BBL (bottom boundary layer) parameterization.
 9. use free slip lateral condition except in the Idonesian through-flow area, Mediterannean Sea, West coast of Greenland ( near cape desolation) and in the northern part of Nares Strait.
10. Use TOP with CFC passive tracers.

### Forcing:
 1. This run uses DFS4.4 forcing, with __absolute__ winds and CORE bulk formulae.
 2. SSS restoring toward Levitus 98, with a piston velocity of 167 mm/day ( 60 days/10 meters).
 3. Run-off from Dai-Trenberth including climatological iceberg contribution from Da Silva

### TOP settings
  TOP was turned on after 2 years of spin-up, and uses initial conditions for CFC11 concentration and inventory from previous ORCA025-G70. A TOP/CFC11 restart file for TRC was created from ORCA025-G70 output on 1960m01d05.

### Output data set: 5d output only [556Gb/year ]
  The data set has been converted to netcdf4/hdf5 with deflation level 1, and compacted into tar files for archiving. Available file types are :
 * **gridT** files : votemper, vosaline, sossheig, somxl010, somixhgt,sohefldo, etc...
 * **gridU** files : vozocrtx, vozotaux
 * **gridV** files : vomecrty, vometauy
 * **gridW** files : vovecrtz, voavttke
 * **flxT** files : 15 atmosperic fluxes and atmospheric variables.
 * **icemod** files : 19 ice model variables
 * **ptrcT** files : CFC11
 * **diadT** files : QTRCCFX, QINTCFC, INVCFC

### Run time files
   The run time files not  indicated in the namelist are:

 * Bathymetry : ```bathymetry_ORCA12_V3.3.nc```
 * Coordinates : ```coordinates_ORCA_R12_lbclnk_no_z.nc```
 * Ice Initialisation : ```ORCA12.L46-MAL95_y1998-2007m01_icemod_initMAL101.nc```
 * Bottom Friction enhancement : ```orca12_bfr_coef_MAL101.nc```
 * AABW damping mask : ```ORCA12.L46_dmp_mask.nc```
 * CFC atmospherique : ```cfc1112_updated4.atm```



### Bibliography:
Molines J.M., B. Barnier, T. Penduff, A.M. Treguier, J. Le Sommer. "ORCA12.L46 climatological and inter annual simulations forced with DFS4.4: GJM02 and MJM88.",  Drakkar Group Experiment report GDRI-DRAKKAR-2014-03-19 (2014) [Technical Report](https://www.drakkar-ocean.eu/publications/reports/orca12_reference_experiments_2014)
