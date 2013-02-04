Metamaterial Superlens: Finite Size Effects.
===============================================
Code should be run in fortran90/95 which generates data files which are then plotted by the MATLAB scripts.
Controlled also by shell scripts

Notes
--------
* The gaussian codes are newer still and use a gaussian waveform instead of single ray
* The parametised code is newer than the unparametised.
* The unparametised code is probably useless for the dielectric but may be useful for the metal.

* Data (fortran output) is stored in /data as .dat files with four columns which are the x,z,abs(field) and realpart(field) values respectively the filenames illustrate the parameter values.

* Plots (MATLAB output) are stored in /plots as .png files with parameters recorded in the titles of the plots and filenames


History
--------

* 28/03/2012 - correctly parametrised the gaussian (i.e. using g=g_unparam/lambda resulting in a division by eta**2), got negative refraction towork correctly by fixing signs of kz2, took plots with various angles of incidence and g values and demonstrated negative refraction 

* 27/03/2012 - set up directories for data and plots, made some basic plots of positive refraction in dielectric (parametised), tried and failed to code for negative refractive index and imaginary epsilon (losses/metals).

* 26/02/2012 - set up github, uploaded fortran code and matlab scripts
