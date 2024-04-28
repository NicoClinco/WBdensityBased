# WBdensityBased

This directory contains the code for a finite-volume density-based solver
with the well-balanced property for the gravity source term.

The directory 'DCfoam' contains the application for the Density-current, while
The directory 'RTBfoam' contains the application for the Smooth-rising bubble.
The directory 'HydroFoam' contains the application for the hydrostatic atmosphere
over a flat terrain

The application DCFoam must be compiled after having compiled the code in
'AtmosDBasedEnergyDC' that contains the numeric flux formulation for the density
current.

The application RTBFoam must be compiled after having compiled the code
'AtmosDBasedEnergyRTB'.

The application HydroFoam must be compiled after having compiled the code
'AtmosDBasedHydro'.
