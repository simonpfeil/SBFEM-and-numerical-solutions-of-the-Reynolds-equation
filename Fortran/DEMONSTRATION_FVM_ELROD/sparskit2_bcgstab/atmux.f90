!-----------------------------------------------------------------------
      subroutine atmux (n, x, y, a, ja, ia) 
      real*8 x(*), y(*), a(*) 
      integer n, ia(*), ja(*) 
!-----------------------------------------------------------------------
!         transp( A ) times a vector                                    
!-----------------------------------------------------------------------
! multiplies the transpose of a matrix by a vector when the original    
! matrix is stored in compressed sparse row storage. Can also be        
! viewed as the product of a matrix by a vector when the original       
! matrix is stored in the compressed sparse column format.              
!-----------------------------------------------------------------------
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
! y     = real array of length n, containing the product y=transp(A)*x  
!                                                                       
!-----------------------------------------------------------------------
!     local variables                                                   
!                                                                       
      integer i, k 
!-----------------------------------------------------------------------
!                                                                       
!     zero out output vector                                            
!                                                                       
      do 1 i=1,n 
         y(i) = 0.0 
    1 continue 
!                                                                       
! loop over the rows                                                    
!                                                                       
      do 100 i = 1,n 
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k) 
   99    continue 
  100 continue 
!                                                                       
      return 
!-------------end-of-atmux----------------------------------------------
!-----------------------------------------------------------------------
      END            subroutine                               
