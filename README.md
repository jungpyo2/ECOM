# ECOM
A MHD axisymmetric Equilibrium solver via COnformal Mapping (ECOM)

This directory contains source code for ECOM and several examples.

If you want to compile ECOM, just type

   make

If you run the execution file (./ecom), it will find the namelist input file "ecom.in",
and produce some output files including "ecom.out".
If you need more documatation for ECOM or input namelist, please see the website
http://cims.nyu.edu/~jungpyo/projects/ECOM.html

A example for Alcator-Cmod tokamak is in the subdirectory "Cmod_example".
Another simple example of ecom.in for Solovev equilbria can be

&GSparameter
ibytpe=0
iptype=0
eps=1e-14
nt2=128
nsub=8
kcehb=16
maxiter=20
R0=1.0d0
Z0=0.0d0
F0=1.0d0
q0=1.0d0
reps=0.32d0
rkappa=1.7d0
delta=0.33d0
/

Please email  Jungpyo Lee (jungpyo@psfc.mit.edu) and Antoine Cerfon (cerfon@cims.nyu.edu), if you have any question on the code.
