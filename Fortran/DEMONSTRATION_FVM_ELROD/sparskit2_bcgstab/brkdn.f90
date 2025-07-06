!-----------------------------------------------------------------------
      logical function brkdn(alpha, ipar) 
      implicit none 
      integer ipar(16) 
      real*8 alpha, beta, zero, one 
      parameter (zero=0.0D0, one=1.0D0) 
!-----------------------------------------------------------------------
!     test whether alpha is zero or an abnormal number, if yes,         
!     this routine will return .true.                                   
!                                                                       
!     If alpha == 0, ipar(1) = -3,                                      
!     if alpha is an abnormal number, ipar(1) = -9.                     
!-----------------------------------------------------------------------
      brkdn = .false. 
      if (alpha.gt.zero) then 
         beta = one / alpha 
         if (.not. beta.gt.zero) then 
            brkdn = .true. 
            ipar(1) = -9 
         endif 
      else if (alpha.lt.zero) then 
         beta = one / alpha 
         if (.not. beta.lt.zero) then 
            brkdn = .true. 
            ipar(1) = -9 
         endif 
      else if (alpha.eq.zero) then 
         brkdn = .true. 
         ipar(1) = -3 
      else 
         brkdn = .true. 
         ipar(1) = -9 
      endif 
      return 
      END        function                                   
