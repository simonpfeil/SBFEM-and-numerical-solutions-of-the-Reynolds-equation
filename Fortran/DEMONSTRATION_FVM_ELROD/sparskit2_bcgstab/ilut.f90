!----------------------------------------------------------------------c
      subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr) 
!-----------------------------------------------------------------------
      implicit none 
      integer n 
      real*8 a(*),alu(*),w(n+1),droptol 
      integer ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),lfil,iwk,ierr 
!----------------------------------------------------------------------*
!                      *** ILUT preconditioner ***                     *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
!     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
!----------------------------------------------------------------------*
! PARAMETERS                                                            
!-----------                                                            
!                                                                       
! on entry:                                                             
!==========                                                             
! n       = integer. The row dimension of the matrix A. The matrix      
!                                                                       
! a,ja,ia = matrix stored in Compressed Sparse Row format.              
!                                                                       
! lfil    = integer. The fill-in parameter. Each row of L and each row  
!           of U will have a maximum of lfil elements (excluding the    
!           diagonal element). lfil must be .ge. 0.                     
!           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO 
!           EARLIER VERSIONS.                                           
!                                                                       
! droptol = real*8. Sets the threshold for dropping small terms in the  
!           factorization. See below for details on dropping strategy.  
!                                                                       
!                                                                       
! iwk     = integer. The lengths of arrays alu and jlu. If the arrays   
!           are not big enough to store the ILU factorizations, ilut    
!           will stop with an error message.                            
!                                                                       
! On return:                                                            
!===========                                                            
!                                                                       
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in       
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
!           contains the i-th row of L (excluding the diagonal entry=1) 
!           followed by the i-th row of U.                              
!                                                                       
! ju      = integer array of length n containing the pointers to        
!           the beginning of each row of U in the matrix alu,jlu.       
!                                                                       
! ierr    = integer. Error message with the following meaning.          
!           ierr  = 0    --> successful return.                         
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.          
!                            (The elimination process has generated a   
!                            row in L or U whose length is .gt.  n.)    
!           ierr  = -2   --> The matrix L overflows the array al.       
!           ierr  = -3   --> The matrix U overflows the array alu.      
!           ierr  = -4   --> Illegal value for lfil.                    
!           ierr  = -5   --> zero row encountered.                      
!                                                                       
! work arrays:                                                          
!=============                                                          
! jw      = integer work array of length 2*n.                           
! w       = real work array of length n+1.                              
!                                                                       
!---------------------------------------------------------------------- 
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]       
! jw(n+1:2n)  stores nonzero indicators                                 
!                                                                       
! Notes:                                                                
! ------                                                                
! The diagonal elements of the input matrix must be  nonzero (at least  
! 'structurally').                                                      
!                                                                       
!----------------------------------------------------------------------*
!---- Dual drop strategy works as follows.                             *
!                                                                      *
!     1) Theresholding in L and U as set by droptol. Any element whose *
!        magnitude is less than some tolerance (relative to the abs    *
!        value of diagonal element in u) is dropped.                   *
!                                                                      *
!     2) Keeping only the largest lfil elements in the i-th row of L   *
!        and the largest lfil elements in the i-th row of U (excluding *
!        diagonal elements).                                           *
!                                                                      *
! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
! keeping  the largest  elements in  each row  of L  and U.   Taking   *
! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
! (however, fill-in is then mpredictible).                             *
!----------------------------------------------------------------------*
!     locals                                                            
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len 
      real*8 tnorm, t, abs, s, fact 
      
      ! ju0
      ! k
      ! j1,j2
      ! j
      ! ii = loop index for rows
      ! i
      ! lenl = length of lower diagonal part of current row ii
      ! lenu = length of upper diagonal part of current row ii including the diagonal entry
      ! jj
      ! jrow jrow = row index < ii to eliminate nonzero in lower part
      ! jpos
      ! len 
      
      ! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]       
      ! jw(n+1:2n)  stores nonzero indicators 
      
      ! tnorm = norm_1 of current row divided by number of symbolic nonzeros i this row
      ! t
      ! abs
      ! s
      ! fact 
!       
      if (lfil .lt. 0) goto 998  ! exit with ierr = -4
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)    
!     and pointer array.                                                
!-----------------------------------------------------------------------
      ju0 = n+2 
      jlu(1) = ju0 
!                                                                       
!     initialize nonzero indicator array.                               
!                                                                       
      do 1 j=1,n 
         jw(n+j)  = 0 
    1 continue 
!-----------------------------------------------------------------------
!     beginning of main loop.                                           
!-----------------------------------------------------------------------
      do 500 ii = 1, n 
         ! compute norm_1 of ii-th row: (JK) 
         j1 = ia(ii) 
         j2 = ia(ii+1) - 1 
         tnorm = 0.0d0 
         do 501 k=j1,j2 
            tnorm = tnorm+abs(a(k)) 
  501    continue 
         if (tnorm .eq. 0.0) goto 999 ! zero-row encountered, exit with ierr = -5 (JK)
         tnorm = tnorm/real(j2-j1+1)  ! divide norm_1 by nr. of nonzero elements (JK)
!                                                                       
!     unpack L-part and U-part of row of A in arrays w:
!     - compute lenghts of lower and upper part of row ii, lenl and lenu, and 
!       store lower part in w(1:lenl), jw(1:lenl) 
!       and upper part in w(ii:ii+lenu-1), jw(ii:ii+lenu-1)
!       and in jw(n+1:2n) indices such that jw(jw(n+k)) = k for k /= ii
!       
         lenu = 1 
         lenl = 0 
         jw(ii) = ii 
         w(ii) = 0.0 
         jw(n+ii) = ii 
!                                                                       
         do 170  j = j1, j2 
            k = ja(j) 
            t = a(j) 
            if (k .lt. ii) then  ! element of lower part of A (JK)
               lenl = lenl+1 
               jw(lenl) = k 
               w(lenl) = t 
               jw(n+k) = lenl 
            else if (k .eq. ii) then  ! diagonal element of A (JK)
               w(ii) = t 
            else  ! ! element of upper part of A (JK)
               lenu = lenu+1 
               jpos = ii+lenu-1 
               jw(jpos) = k 
               w(jpos) = t 
               jw(n+k) = jpos 
            endif 
170      continue 
         
         jj = 0 
         len = 0 
!                                                                       
!     eliminate previous rows                                           
!      
! bubble sort lower part of row ii in increasing j-index oder (loop 150):
  150    jj = jj+1  !
         if (jj .gt. lenl) goto 160  ! exit row elimintion loop (and continue withreset double-pointer to zero (U-part)) (JK)
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.         
!-----------------------------------------------------------------------
         jrow = jw(jj)  ! jrow = to be determined row index < ii to eliminate nonzero in lower part
         k = jj  ! to determine index of smallest jw(k), k=jj+1, ..., lenl
!                                                                       
!     determine smallest column index of lower part of row ii                                  
!                                                                       
         do 151 j=jj+1,lenl 
            if (jw(j) .lt. jrow) then 
               jrow = jw(j) 
               k = j 
            endif 
  151    continue 
!                                                                       
         if (k .ne. jj) then ! exchange elements
!     exchange in jw                                                    
            j = jw(jj) 
            jw(jj) = jw(k) 
            jw(k) = j 
!     exchange in jr                                                    
            jw(n+jrow) = jj 
            jw(n+j) = k 
!     exchange in w                                                     
            s = w(jj) 
            w(jj) = w(k) 
            w(k) = s 
         endif 
!                                                                       
!     zero out element in row by setting jw(n+jrow) to zero.            
!                                                                       
         jw(n+jrow) = 0 
!                                                                       
!     get the multiplier for row to be eliminated (jrow).               
!                                                                       
         fact = w(jj)*alu(jrow) 
         if (abs(fact) .le. droptol) goto 150  ! cycle (JK)
!                                                                       
!     combine current row and row jrow                                  
!                                                                       
         do 203 k = ju(jrow), jlu(jrow+1)-1 
            s = fact*alu(k) 
            j = jlu(k) 
            jpos = jw(n+j) 
            if (j .ge. ii) then 
!                                                                       
!     dealing with upper part.                                          
!                                                                       
               if (jpos .eq. 0) then 
!                                                                       
!     this is a fill-in element                                         
!                                                                       
                  lenu = lenu+1 
                  if (lenu .gt. n) goto 995  ! return with ierr = -1: incomprehensible error. Matrix must be wrong.
                  i = ii+lenu-1 
                  jw(i) = j 
                  jw(n+j) = i 
                  w(i) = - s 
               else 
!                                                                       
!     this is not a fill-in element                                     
!                                                                       
                  w(jpos) = w(jpos) - s 
                                                                        
               endif 
            else 
!                                                                       
!     dealing  with lower part.                                         
!                                                                       
               if (jpos .eq. 0) then 
!                                                                       
!     this is a fill-in element                                         
!                                                                       
                  lenl = lenl+1 
                  if (lenl .gt. n) goto 995 
                  jw(lenl) = j 
                  jw(n+j) = lenl 
                  w(lenl) = - s 
               else 
!                                                                       
!     this is not a fill-in element                                     
!                                                                       
                  w(jpos) = w(jpos) - s 
               endif 
            endif 
  203    continue 
!                                                                       
!     store this pivot element -- (from left to right -- no danger of   
!     overlap with the working elements in L (pivots).                  
!                                                                       
         len = len+1 
         w(len) = fact 
         jw(len)  = jrow 
         goto 150  ! row elimination loop (JK)
         
160      continue 
         
!                                                                       
!     reset double-pointer to zero (U-part)                             
!                                                                       
         do 308 k=1, lenu 
            jw(n+jw(ii+k-1)) = 0 
  308    continue 
!                                                                       
!     update L-matrix                                                   
!                                                                       
         lenl = len 
         len = min0(lenl,lfil) 
!                                                                       
!     sort by quick-split                                               
!                                                                       
         call qsplit (w,jw,lenl,len) 
!                                                                       
!     store L-part                                                      
!                                                                       
         do 204 k=1, len 
            if (ju0 .gt. iwk) goto 996 
            alu(ju0) =  w(k) 
            jlu(ju0) =  jw(k) 
            ju0 = ju0+1 
  204    continue 
!                                                                       
!     save pointer to beginning of row ii of U                          
!                                                                       
         ju(ii) = ju0 
!                                                                       
!     update U-matrix -- first apply dropping strategy                  
!                                                                       
         len = 0 
         do k=1, lenu-1 
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               len = len+1 
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif 
         enddo 
         lenu = len+1 
         len = min0(lenu,lfil) 
!                                                                       
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len) 
!                                                                       
!     copy                                                              
!                                                                       
         t = abs(w(ii)) 
         if (len + ju0 .gt. iwk) goto 997 
         do 302 k=ii+1,ii+len-1 
            jlu(ju0) = jw(k) 
            alu(ju0) = w(k) 
            t = t + abs(w(k) ) 
            ju0 = ju0+1 
  302    continue 
!                                                                       
!     store inverse of diagonal element of u                            
!                                                                       
         if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm 
!                                                                       
         alu(ii) = 1.0d0/ w(ii) 
!                                                                       
!     update pointer to beginning of next row of U.                     
!                                                                       
         jlu(ii+1) = ju0 
!-----------------------------------------------------------------------
!     end main loop                                                     
!-----------------------------------------------------------------------
  500 continue 
      ierr = 0 
      return 
!                                                                       
!     incomprehensible error. Matrix must be wrong.                     
!                                                                       
  995 ierr = -1 
      return 
!                                                                       
!     insufficient storage in L.                                        
!                                                                       
  996 ierr = -2 
      return 
!                                                                       
!     insufficient storage in U.                                        
!                                                                       
  997 ierr = -3 
      return 
!                                                                       
!     illegal lfil entered.                                             
!                                                                       
  998 ierr = -4 
      return 
!                                                                       
!     zero row encountered                                              
!                                                                       
  999 ierr = -5 
      return 
!----------------end-of-ilut--------------------------------------------
!-----------------------------------------------------------------------
      END             subroutine                              
