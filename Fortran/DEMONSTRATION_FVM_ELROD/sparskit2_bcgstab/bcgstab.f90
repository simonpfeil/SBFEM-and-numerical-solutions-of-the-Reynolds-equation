!-----------------------------------------------------------------------
      subroutine bcgstab(n, rhs, sol, ipar, fpar, w) 
      implicit none 
      integer n, ipar(16) 
      real*8 rhs(n), sol(n), fpar(16), w(n,8) 
!-----------------------------------------------------------------------
!     BCGSTAB --- Bi Conjugate Gradient stabilized (BCGSTAB)            
!     This is an improved BCG routine. (1) no matrix transpose is       
!     involved. (2) the convergence is smoother.                        
!                                                                       
!                                                                       
!     Algorithm:                                                        
!     Initialization - r = b - A x, r0 = r, p = r, rho = (r0, r),       
!     Iterate -                                                         
!     (1) v = A p                                                       
!     (2) alpha = rho / (r0, v)                                         
!     (3) s = r - alpha v                                               
!     (4) t = A s                                                       
!     (5) omega = (t, s) / (t, t)                                       
!     (6) x = x + alpha * p + omega * s                                 
!     (7) r = s - omega * t                                             
!     convergence test goes here                                        
!     (8) beta = rho, rho = (r0, r), beta = rho * alpha / (beta * omega)
!         p = r + beta * (p - omega * v)                                
!                                                                       
!     in this routine, before successful return, the fpar's are         
!     fpar(3) == initial (preconditionied-)residual norm                
!     fpar(4) == target (preconditionied-)residual norm                 
!     fpar(5) == current (preconditionied-)residual norm                
!     fpar(6) == current residual norm or error                         
!     fpar(7) == current rho (rhok = <r, r0>)                           
!     fpar(8) == alpha                                                  
!     fpar(9) == omega                                                  
!                                                                       
!     Usage of the work space W                                         
!     w(:, 1) = r0, the initial residual vector                         
!     w(:, 2) = r, current residual vector                              
!     w(:, 3) = s                                                       
!     w(:, 4) = t                                                       
!     w(:, 5) = v                                                       
!     w(:, 6) = p                                                       
!     w(:, 7) = tmp, used in preconditioning, etc.                      
!     w(:, 8) = delta x, the correction to the answer is accumulated    
!               here, so that the right-preconditioning may be applied  
!               at the end                                              
!-----------------------------------------------------------------------
!     external routines used                                            
!                                                                       
     ! real*8 distdot 
     ! logical stopbis, brkdn 
     ! external distdot, stopbis, brkdn 
!                                                                       
      real*8 one 
      parameter(one=1.0D0) 
!                                                                       
!     local variables                                                   
!                                                                       
      integer i 
      real*8 alpha,beta,rho,omega 
      logical lp, rp 
      save lp, rp 
!                                                                       
!     where to go                                                       
!                                                                       
      if (ipar(1).gt.0) then 
         goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110) ipar(10) 
      else if (ipar(1).lt.0) then 
         goto 900 
      endif 
!                                                                       
!     call the initialization routine                                   
!                                                                       
      call bisinit(ipar,fpar,8*n,1,lp,rp,w) 
      if (ipar(1).lt.0) return 
!                                                                       
!     perform a matvec to compute the initial residual                  
!                                                                       
      ipar(1) = 1 
      ipar(8) = 1 
      ipar(9) = 1 + n 
      do i = 1, n 
         w(i,1) = sol(i) 
      enddo 
      ipar(10) = 1 
      return 
   10 ipar(7) = ipar(7) + 1 
      ipar(13) = ipar(13) + 1 
      do i = 1, n 
         w(i,1) = rhs(i) - w(i,2) 
      enddo 
      fpar(11) = fpar(11) + n 
      if (lp) then 
         ipar(1) = 3 
         ipar(10) = 2 
         return 
      endif 
!                                                                       
   20 if (lp) then 
         do i = 1, n 
            w(i,1) = w(i,2) 
            w(i,6) = w(i,2) 
         enddo 
      else 
         do i = 1, n 
            w(i,2) = w(i,1) 
            w(i,6) = w(i,1) 
         enddo 
      endif 
!                                                                       
      fpar(7) = distdot(n,w,1,w,1) 
      fpar(11) = fpar(11) + 2 * n 
      fpar(5) = sqrt(fpar(7)) 
      fpar(3) = fpar(5) 
      if (abs(ipar(3)).eq.2) then 
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2) 
         fpar(11) = fpar(11) + 2 * n 
      else if (ipar(3).ne.999) then 
         fpar(4) = fpar(1) * fpar(3) + fpar(2) 
      endif 
      if (ipar(3).ge.0) fpar(6) = fpar(5) 
      if (ipar(3).ge.0 .and. fpar(5).le.fpar(4) .and.                   &
     &     ipar(3).ne.999) then                                         
         goto 900 
      endif 
!                                                                       
!     beginning of the iterations                                       
!                                                                       
!     Step (1), v = A p                                                 
   30 if (rp) then 
         ipar(1) = 5 
         ipar(8) = 5*n+1 
         if (lp) then 
            ipar(9) = 4*n + 1 
         else 
            ipar(9) = 6*n + 1 
         endif 
         ipar(10) = 3 
         return 
      endif 
!                                                                       
   40 ipar(1) = 1 
      if (rp) then 
         ipar(8) = ipar(9) 
      else 
         ipar(8) = 5*n+1 
      endif 
      if (lp) then 
         ipar(9) = 6*n + 1 
      else 
         ipar(9) = 4*n + 1 
      endif 
      ipar(10) = 4 
      return 
   50 if (lp) then 
         ipar(1) = 3 
         ipar(8) = ipar(9) 
         ipar(9) = 4*n + 1 
         ipar(10) = 5 
         return 
      endif 
!                                                                       
   60 ipar(7) = ipar(7) + 1 
!                                                                       
!     step (2)                                                          
      alpha = distdot(n,w(1,1),1,w(1,5),1) 
      fpar(11) = fpar(11) + 2 * n 
      if (brkdn(alpha, ipar)) goto 900 
      alpha = fpar(7) / alpha 
      fpar(8) = alpha 
!                                                                       
!     step (3)                                                          
      do i = 1, n 
         w(i,3) = w(i,2) - alpha * w(i,5) 
      enddo 
      fpar(11) = fpar(11) + 2 * n 
!                                                                       
!     Step (4): the second matvec -- t = A s                            
!                                                                       
      if (rp) then 
         ipar(1) = 5 
         ipar(8) = n+n+1 
         if (lp) then 
            ipar(9) = ipar(8)+n 
         else 
            ipar(9) = 6*n + 1 
         endif 
         ipar(10) = 6 
         return 
      endif 
!                                                                       
   70 ipar(1) = 1 
      if (rp) then 
         ipar(8) = ipar(9) 
      else 
         ipar(8) = n+n+1 
      endif 
      if (lp) then 
         ipar(9) = 6*n + 1 
      else 
         ipar(9) = 3*n + 1 
      endif 
      ipar(10) = 7 
      return 
   80 if (lp) then 
         ipar(1) = 3 
         ipar(8) = ipar(9) 
         ipar(9) = 3*n + 1 
         ipar(10) = 8 
         return 
      endif 
   90 ipar(7) = ipar(7) + 1 
!                                                                       
!     step (5)                                                          
      omega = distdot(n,w(1,4),1,w(1,4),1) 
      fpar(11) = fpar(11) + n + n 
      if (brkdn(omega,ipar)) goto 900 
      omega = distdot(n,w(1,4),1,w(1,3),1) / omega 
      fpar(11) = fpar(11) + n + n 
      if (brkdn(omega,ipar)) goto 900 
      fpar(9) = omega 
      alpha = fpar(8) 
!                                                                       
!     step (6) and (7)                                                  
      do i = 1, n 
         w(i,7) = alpha * w(i,6) + omega * w(i,3) 
         w(i,8) = w(i,8) + w(i,7) 
         w(i,2) = w(i,3) - omega * w(i,4) 
      enddo 
      fpar(11) = fpar(11) + 6 * n + 1 
!                                                                       
!     convergence test                                                  
      if (ipar(3).eq.999) then 
         ipar(1) = 10 
         ipar(8) = 7*n + 1 
         ipar(9) = 6*n + 1 
         ipar(10) = 9 
         return 
      endif 
      if (stopbis(n,ipar,2,fpar,w(1,2),w(1,7),one))  goto 900 
  100 if (ipar(3).eq.999.and.ipar(11).eq.1) goto 900 
!                                                                       
!     step (8): computing new p and rho                                 
      rho = fpar(7) 
      fpar(7) = distdot(n,w(1,2),1,w(1,1),1) 
      omega = fpar(9) 
      beta = fpar(7) * fpar(8) / (fpar(9) * rho) 
      do i = 1, n 
         w(i,6) = w(i,2) + beta * (w(i,6) - omega * w(i,5)) 
      enddo 
      fpar(11) = fpar(11) + 6 * n + 3 
      if (brkdn(fpar(7),ipar)) goto 900 
!                                                                       
!     end of an iteration                                               
!                                                                       
      goto 30 
!                                                                       
!     some clean up job to do                                           
!                                                                       
  900 if (rp) then 
         if (ipar(1).lt.0) ipar(12) = ipar(1) 
         ipar(1) = 5 
         ipar(8) = 7*n + 1 
         ipar(9) = ipar(8) - n 
         ipar(10) = 10 
         return 
      endif 
  110 if (rp) then 
         call tidycg(n,ipar,fpar,sol,w(1,7)) 
      else 
         call tidycg(n,ipar,fpar,sol,w(1,8)) 
      endif 
!                                                                       
      return 
!-----end-of-bcgstab                                                    
      END           subroutine                                
