!-----------------------------------------------------------------------
      subroutine lusol(n, y, x, alu, jlu, ju) 
      integer:: n, jlu(*), ju(*)   
      real*8 x(n), y(n), alu(*) 
        
        
!-----------------------------------------------------------------------
!                                                                       
! This routine solves the system (LU) x = y,                            
! given an LU decomposition of a matrix stored in (alu, jlu, ju)        
! modified sparse row format                                            
!                                                                       
!-----------------------------------------------------------------------
! on entry:                                                             
! n   = dimension of system                                             
! y   = the right-hand-side vector                                      
! alu, jlu, ju                                                          
!     = the LU matrix as provided from the ILU routines.                
!                                                                       
! on return                                                             
! x   = solution of LU x = y.                                           
!-----------------------------------------------------------------------
!                                                                       
! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)         
!       will solve the system with rhs x and overwrite the result on x .
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
      integer i,k 
!                                                                       
! forward solve                                                         
!                                                                       
      do i = 1, n 
        x(i) = y(i) 
        do k=jlu(i),ju(i)-1 
          x(i) = x(i) - alu(k)* x(jlu(k)) 
        enddo 
      enddo 
!                                                                       
!     backward solve.                                                   
!                                                                       
      do i = n, 1, -1 
        do k=ju(i),jlu(i+1)-1 
              x(i) = x(i) - alu(k)*x(jlu(k)) 
        enddo 
        x(i) = alu(i)*x(i) 
      enddo 
!                                                                       
      return 
      END          subroutine                                 
