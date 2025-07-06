! This module relies on an IC factorization available under
!
!  authors: Mark T. Jones, Paul E. Plassmann
!  title: Algorithm 740: Fortran subroutines to compute improved incomplete Cholesky factorizations
!  published in: ACM Transactions on Mathematical Software, vol. 21, no. 1, March 1995, p. 18-19
!  publisher: Association for Computing Machinery
!  address: New York, NY, USA
!  issn: 0098-3500
!  url: https://doi.org/10.1145/200979.200986
!  doi: 10.1145/200979.200986
!
! Please download the subprograms listed under "CONTAINS" and place them in the folder
! "./jones_plassmann_jpicc".

MODULE MODULE_JPICC
  
  CONTAINS
  
  INCLUDE '.\jones_plassmann_jpicc\jpicc.f'
  INCLUDE '.\jones_plassmann_jpicc\dbsort.f'
  INCLUDE '.\jones_plassmann_jpicc\dhsort.f'
  INCLUDE '.\jones_plassmann_jpicc\ibsort.f'
  INCLUDE '.\jones_plassmann_jpicc\ihsort.f'
  
END MODULE
