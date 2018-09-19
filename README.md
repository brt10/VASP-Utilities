# VASP-Utilities

These codes use as input the files from VASP = Vienna Ab initio Simulation Package.

To compile these codes on a unix machine use the Fortran compiler available eg ifort.

mk_neb.f : input two POSCAR files, output : 16 POSCAR files for a NEB VASP simulation

CONTCAR is the input file for VASP which includes lattice vectors and basis atom vectors

mk_super.f : input is a CONTCAR unit cell file, output is a POSCAR file which is a supercell of the unit cell

nn.f : input is a CONTCAR file, outputs nearest neighbor lists with distances
