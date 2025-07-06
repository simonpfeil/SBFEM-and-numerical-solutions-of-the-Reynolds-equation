!----------------------------------------------------------------------c
      subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*) 
!-----------------------------------------------------------------------
!         A times a vector                                              
!-----------------------------------------------------------------------
! multiplies a matrix by a vector using the dot product form            
! Matrix A is stored in compressed sparse row storage.                  
!                                                                       
! on entry:                                                             
!----------                                                             
! n     = row dimension of A                                            
! x     = real array of length equal to the column dimension of         
!         the A matrix.                                                 
! a, ja,                                                                
!    ia = input matrix in compressed sparse row format.                 
!                                                                       
! on return:                                                            
!-----------                                                            
! y     = real array of length n, containing the product y=Ax           
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
      real*8 t 
      integer i, k 
!-----------------------------------------------------------------------
      do 100 i = 1,n 
!                                                                       
!     compute the inner product of row i with vector x                  
!                                                                       
         t = 0.0d0 
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k)) 
   99    continue 
!                                                                       
!     store result in y(i)                                              
!                                                                       
         y(i) = t 
  100 continue 
!                                                                       
      return 
!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------
      END           subroutine                                
