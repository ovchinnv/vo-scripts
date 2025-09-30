
## Description

 This repository contains a Linux bash shell script to illustrate MD
 simulations of two proteins using the open source software OpenMM and a
 solvent boundary potential plugin.  Tested under ArchLinux, the script
 "./make" will download and compile OpenMM and the dynamo plugin, and
 run the test cases in folders 'dhfr' and 'b12'. The output of the test
 cases is included in the repository as a .txz archive.

## Required software:
 Name:	|	Purpose:	|	Source:
 ---------------|----------------|----------------------------------
 Python	|	To run test cases|	Linux distribution
 GNU c, c++, gfortran | to compile code | Linux distribution
 Nvidia CUDA environment |to compile & run code | Linux distribution

## Required hardware:
 A CUDA-capable Nvidia GPU

## How to Run

The procedure has been tested on an Arch Linux workstation, where

`./make`

was used to download and compile the codes and run the simulations.

The complete output from make upon successful execution is provided for
reference in the tarball set-up-wshell-openmm7.3-output.txz

