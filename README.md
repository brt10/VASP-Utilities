# VASP-Utilities

These codes use as input the files from VASP = Vienna Ab initio Simulation Package.

To compile these codes on a unix machine use the Fortran compiler available eg ifort.

mk_super.f : input is a CONTCAR unit cell file and output is a POSCAR file which is a supercell of the unit cell

nn.f : input is a CONTCAR file and returns nearest neighbor lists with distances
