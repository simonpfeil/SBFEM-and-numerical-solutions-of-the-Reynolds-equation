# SBFEM-and-numerical-solutions-of-the-Reynolds-equation
# AUTHOR, AFFILIATION, DATE: Simon Pfeil, OvGU Magdeburg (Germany), 30.06.2025

# This repository contains MATLAB and Fortran algorithms for solving the Reynolds equation for hydrodynamic bearings.
# Concretely, two solution methods are provided:
# - a semi-analytical solution based on the scaled boundary finite element method (SBFEM) and
# - a numerical solution that can be motivated by the finite element method (FEM), the finite volume method (FVM), or the finite difference method (FDM).
# The numerical solution will be called FVM in the following (omitting the FEM or FDM interpretations).

# Further clarifications are given by the comments in the respective source files, in combination with the illustrations in "sketches.pdf".
# The computational models are discussed thoroughly in the thesis "Simulating Hydrodynamic Bearings with the Scaled Boundary Finite Element Method" by Simon Pfeil, but that document is not publically available yet.

# The Fortran variants of the SBFEM and FVM algorithms use BLAS and LAPACK routines, which are available, for example, through the Intel oneAPI Math Kernel Library (oneMKL) included in the Intel oneAPI Base Toolkit.
# The Fortran variant of the FVM algorithm relies on further subprograms by other authors, including an implementation of a CG solver by Steffen Nitzschke (OvGU Magdeburg), an ILU factorization and BiCGStab solver by Yousef Saad, and an IC factorization by Mark T. Jones and Paul E. Plassmann. The CG algorithm is included in - and may be redistributed along with the other content of - this repository. The other mentioned algorithms (by Yousef Saad as well as by Mark T. Jones and Paul E. Plassmann) are NOT included in this repository; please follow the instructions given in FORTRAN/DEMONSTRATION_FVM_ELROD/MODULE_SPARSKIT2_BCGSTAB.f90 and FORTRAN/DEMONSTRATION_FVM_ELROD/MODULE_JPICC.f90 to obtain and include them.

# Disclaimer:
# This selection of algorithms is offered freely in the spirit of being helpful, but it comes with no assurances whatsoever. There is no guarantee that it will function as intended or meet any specific needs. No liability is accepted, and there are no promises — explicit or implied — regarding its quality, reliability, or suitability for any purpose.
