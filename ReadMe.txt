========================================================================
    			Project Overview
========================================================================
This is multigroup, 3D neutron diffusion equation solver using FDM scheme
multiple sample input files are available to see the input
program reads the inp.txt

Program File "MAIN.f90" calls the other subroutines for calulation i.e.,
module.f90
SubReadData.f90
SubCoeffs.f90
SubPowerMethod.f90
SubSOR.f90

Results are critaclity value and Flux distribution