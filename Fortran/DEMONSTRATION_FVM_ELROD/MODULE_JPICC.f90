! SOURCE OF THE PROGRAMS CONTAINED IN THIS MODULE IS THE FOLLOWING PUBLICATION:
! authors: Mark T. Jones, Paul E. Plassmann
! title: Algorithm 740: Fortran subroutines to compute improved incomplete Cholesky factorizations
! published in: ACM Transactions on Mathematical Software, vol. 21, no. 1, March 1995, p. 18-19
! publisher: Association for Computing Machinery
! address: New York, NY, USA
! issn: 0098-3500
! url: https://doi.org/10.1145/200979.200986
! doi: 10.1145/200979.200986

! COPYRIGHT, AS STATED IN THE WORK REFERENCED ABOVE:
! "Permission to copy without fee all or part of this material is granted provided that the copies 
! are not made or distributed for direct commercial advantage, the ACM copyright notice and the 
! title of the publication and its date appear, and notice is given that copying is by permission 
! of the Association for Computing Machinery. To copy otherwise, or to republish, requires a fee 
! and/or specific permission. Copyright 1995 ACM 0098-3500/95/0300-0018 $03,50"

MODULE MODULE_JPICC
  
  CONTAINS
  
  INCLUDE '.\jones_plassmann_jpicc\jpicc.f'
  INCLUDE '.\jones_plassmann_jpicc\dbsort.f'
  INCLUDE '.\jones_plassmann_jpicc\dhsort.f'
  INCLUDE '.\jones_plassmann_jpicc\ibsort.f'
  INCLUDE '.\jones_plassmann_jpicc\ihsort.f'
  
END MODULE
