! This module relies on a selection of algorithms for solving equation systems by means of
! BCGSTAB. The files listed under "CONTAINS" correspond to particular subprograms available
! under
!
!  Yousef Saad: SPARSKIT2. University of Minnesota, Department of Computer Science and Engineering, 
!  200 Union Street S.E., Minneapolis, MN 55455 USA, [saad -at- umn -dot- edu],
!  https://www-users.cse.umn.edu/~saad/software/SPARSKIT/
!
! Please download these subprograms from there and place them in the folder "./sparskit2_bcgstab" 
! using the corresponding file names.

MODULE MODULE_SPARSKIT2_BCGSTAB

IMPLICIT NONE

CONTAINS
  
  include '.\sparskit2_bcgstab\amux.f90'
  include '.\sparskit2_bcgstab\atmux.f90'
  include '.\sparskit2_bcgstab\bcgstab.f90'
  include '.\sparskit2_bcgstab\bisinit.f90'
  include '.\sparskit2_bcgstab\blas1.f90'
  include '.\sparskit2_bcgstab\brkdn.f90'
  include '.\sparskit2_bcgstab\cg.f90'
  include '.\sparskit2_bcgstab\coocsr.f90'
  include '.\sparskit2_bcgstab\distdot.f90'
  include '.\sparskit2_bcgstab\ilu0.f90'
  include '.\sparskit2_bcgstab\ilut.f90'
  include '.\sparskit2_bcgstab\lusol.f90'
  include '.\sparskit2_bcgstab\lutsol.f90'
  include '.\sparskit2_bcgstab\qsplit.f90'
  include '.\sparskit2_bcgstab\runrc.f90'
  include '.\sparskit2_bcgstab\stopbis.f90'
  include '.\sparskit2_bcgstab\tidycg.f90'
  
END MODULE
