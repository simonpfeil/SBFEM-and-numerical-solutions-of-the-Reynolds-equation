!-----------------------------------------------------------------------
      subroutine lutsol(n, y, x, alu, jlu, ju) 
      
        integer n, jlu(*), ju(*)  
        real*8 x(n), y(n), alu(*) 
       
!-----------------------------------------------------------------------
!                                                                       
! This routine solves the system  Transp(LU) x = y,                     
! given an LU decomposition of a matrix stored in (alu, jlu, ju)        
! modified sparse row format. Transp(M) is the transpose of M.          
!-----------------------------------------------------------------------
! on entry:                                                             
! n   = dimension of system                                             
! y   = the right-hand-side vector                                      
! alu, jlu, ju                                                          
!     = the LU matrix as provided from the ILU routines.                
!                                                                       
! on return                                                             
! x   = solution of transp(LU) x = y.                                   
!-----------------------------------------------------------------------
!                                                                       
! Note: routine is in place: call lutsol (n, x, x, alu, jlu, ju)        
!       will solve the system with rhs x and overwrite the result on x .
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
        integer i,k 
!                                                                       
        do 10 i = 1, n 
           x(i) = y(i) 
   10   continue 
!                                                                       
! forward solve (with U^T)                                              
!                                                                       
        do i = 1, n 
           x(i) = x(i) * alu(i) 
           do k=ju(i),jlu(i+1)-1 
              x(jlu(k)) = x(jlu(k)) - alu(k)* x(i) 
           enddo 
        enddo 
!                                                                       
!     backward solve (with L^T)                                         
!                                                                       
      do i = n, 1, -1 
        do k=jlu(i),ju(i)-1 
              x(jlu(k)) = x(jlu(k)) - alu(k)*x(i) 
        enddo 
      enddo 
!                                                                       
      return 
      END         subroutine                                  
