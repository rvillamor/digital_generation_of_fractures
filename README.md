# digital_generation_of_fractures
A series of MATLAB scripts to generate self-affine fractures and export them as STL bodies

This repository contains a series of scripts that I developed* to generate STL files of fracture geometries. These could be use for digital fabrication (e.g., 3D printing) or numerical simulations. Please, have a look to the following files:

% driver_ExecuteMe.m: This script shows, step-by-step, how to generate the STL files of a rough fractures using the RSG_brown1995() and fields2bodies() functions.

% RSG_brown1995(): This function implements the model described in Brown (1995) --> 
ABSTRACT: "A simple mathematical model of rough-walled fractures in rock is described which requires the specification of only three main parameters: the fractal dimension, the rms roughness at a reference length scale, anda length scale describing the degree of mismatch between the two fracture surfaces"

% fields2bodies(): A simple function to convert the 2D topograpgy fields generated with RSG_brown1995() into STL files. Specifically one can choose to use such 2D fields to generate "prismatic" and "cylindrical" specimens that can be used for 3D-printing + experimens (which are much more fun than numerical simulations).

******************
* Important note *
******************
The functions stlwrite() and surf2solid() included here are the sole property of Sven Holcombe (please see their respective
lincenses attached). I HAVE NOT developed those, and they are only included here for convenience. These two functions available from:
(1) https://www.mathworks.com/matlabcentral/fileexchange/42876-surf2solid-make-a-solid-volume-from-a-surface-for-3d-printing and (2) https://www.mathworks.com/matlabcentral/fileexchange/20922-stlwrite-write-ascii-or-binary-stl-files

