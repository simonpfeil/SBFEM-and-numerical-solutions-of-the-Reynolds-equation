!-----------------------------------------------------------------------
        subroutine qsplit(a,ind,n,ncut) 
        integer n, ind(n),  ncut 
        real*8 a(n) 
        
!-----------------------------------------------------------------------
!     does a quick-sort split of a real array.                          
!     on input a(1:n). is a real array                                  
!     on output a(1:n) is permuted such that its elements satisfy:      
!                                                                       
!     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and                   
!     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut                       
!                                                                       
!     ind(1:n) is an integer array which permuted in the same way as a(*
!-----------------------------------------------------------------------
        real*8 tmp, abskey 
        integer itmp, first, last ,mid,j
!-----                                                                  
        first = 1 
        last = n 
        if (ncut .lt. first .or. ncut .gt. last) return 
!                                                                       
!     outer loop -- while mid .ne. ncut do                              
!                                                                       
    1   mid = first 
        abskey = abs(a(mid)) 
        do 2 j=first+1, last 
           if (abs(a(j)) .gt. abskey) then 
              mid = mid+1 
!     interchange                                                       
              tmp = a(mid) 
              itmp = ind(mid) 
              a(mid) = a(j) 
              ind(mid) = ind(j) 
              a(j)  = tmp 
              ind(j) = itmp 
           endif 
    2   continue 
!                                                                       
!     interchange                                                       
!                                                                       
        tmp = a(mid) 
        a(mid) = a(first) 
        a(first)  = tmp 
!                                                                       
        itmp = ind(mid) 
        ind(mid) = ind(first) 
        ind(first) = itmp 
!                                                                       
!     test for while loop                                               
!                                                                       
        if (mid .eq. ncut) return 
        if (mid .gt. ncut) then 
           last = mid-1 
        else 
           first = mid+1 
        endif 
        goto 1 
!----------------end-of-qsplit------------------------------------------
!-----------------------------------------------------------------------
      END             subroutine                              
