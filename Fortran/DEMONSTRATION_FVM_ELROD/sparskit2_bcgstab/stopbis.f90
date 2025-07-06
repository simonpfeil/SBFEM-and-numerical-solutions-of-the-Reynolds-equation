!-----------------------------------------------------------------------
      logical function stopbis(n,ipar,mvpi,fpar,r,delx,sx) 
      implicit none 
      integer n,mvpi,ipar(16) 
      real*8 fpar(16), r(n), delx(n), sx!, distdot 
      !external distdot 
!-----------------------------------------------------------------------
!     function for determining the stopping criteria. return value of   
!     true if the stopbis criteria is satisfied.                        
!-----------------------------------------------------------------------
      if (ipar(11) .eq. 1) then 
         stopbis = .true. 
      else 
         stopbis = .false. 
      endif 
      if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then 
         ipar(1) = -1 
         stopbis = .true. 
      endif 
      if (stopbis) return 
!                                                                       
!     computes errors                                                   
!                                                                       
      fpar(5) = sqrt(distdot(n,r,1,r,1)) 
      fpar(11) = fpar(11) + 2 * n 
      if (ipar(3).lt.0) then 
!                                                                       
!     compute the change in the solution vector                         
!                                                                       
         fpar(6) = sx * sqrt(distdot(n,delx,1,delx,1)) 
         fpar(11) = fpar(11) + 2 * n 
         if (ipar(7).lt.mvpi+mvpi+1) then 
!                                                                       
!     if this is the end of the first iteration, set fpar(3:4)          
!                                                                       
            fpar(3) = fpar(6) 
            if (ipar(3).eq.-1) then 
               fpar(4) = fpar(1) * fpar(3) + fpar(2) 
            endif 
         endif 
      else 
         fpar(6) = fpar(5) 
      endif 
!                                                                       
!     .. the test is struct this way so that when the value in fpar(6)  
!       is not a valid number, STOPBIS is set to .true.                 
!                                                                       
      if (fpar(6).gt.fpar(4)) then 
         stopbis = .false. 
         ipar(11) = 0 
      else 
         stopbis = .true. 
         ipar(11) = 1 
      endif 
!                                                                       
      return 
      END        function                                   
