      subroutine runrc(n,rhs,sol,ipar,fpar,wk,guess,a,ja,ia,            &
     &     au,jau,ju,solver)                                            
      implicit none 
      integer n,ipar(16),ia(n+1),ja(*),ju(*),jau(*) 
      real*8 fpar(16),rhs(n),sol(n),guess(n),wk(*),a(*),au(*) 
      external solver 
!-----------------------------------------------------------------------
!     the actual tester. It starts the iterative linear system solvers  
!     with a initial guess suppied by the user.                         
!                                                                       
!     The structure {au, jau, ju} is assumed to have the output from    
!     the ILU* routines in ilut.f.                                      
!                                                                       
!-----------------------------------------------------------------------
!     local variables                                                   
!                                                                       
      integer i, iou, its 
      real*8 res!, dnrm2 
!     real dtime, dt(2), time                                           
!     external dtime                                                    
                     ! print information ? (JK)                         
      logical prinfo 
      !external dnrm2 
      save its,res 
                                                                        
!---- print info:                                                       
     prinfo = .false.                                                 
!      prinfo = .true. 

!                                                                       
!     ipar(2) can be 0, 1, 2, please don't use 3                        
!                                                                       
      if (ipar(2).gt.2) then 
         print *, 'I can not do both left and right preconditioning.' 
         return 
      endif 
!                                                                       
!     normal execution                                                  
!                                                                       
      its = 0 
      res = 0.0D0 
!                                                                       
      do i = 1, n 
         sol(i) = guess(i) 
      enddo 
!                                                                       
      
      ipar(1) = 0                                            
      ! Start Iteration
   10 call solver(n,rhs,sol,ipar,fpar,wk) 
!                                                                                                                                             
      if (ipar(1).eq.1) then 
         call amux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia) 
         goto 10 
      else if (ipar(1).eq.2) then 
         call atmux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia) 
         goto 10 
      else if (ipar(1).eq.3 .or. ipar(1).eq.5) then 
         call lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju) 
         goto 10 
      else if (ipar(1).eq.4 .or. ipar(1).eq.6) then 
         call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju) 
         goto 10 
      else if (ipar(1).le.0) then 
         if (ipar(1).eq.0) then
             if (prinfo) then
               print *, 'Iterative sovler has satisfied convergence test.' 
            endif  ! prinfo 
         else if (ipar(1).eq.-1) then 
            print *, 'Iterative solver has iterated too many times.' 
         else if (ipar(1).eq.-2) then 
            print *, 'Iterative solver was not given enough work space.' 
            print *, 'The work space should at least have ', ipar(4),   &
     &           ' elements.'                                           
         else if (ipar(1).eq.-3) then 
            print *, 'Iterative sovler is facing a break-down.' 
         else 
            print *, 'Iterative solver terminated. code =', ipar(1) 
         endif 
      endif 
      ! Ende Iteration
      
!      if (prinfo) then
!        write (iou, *) ipar(7), real(fpar(6)) 
!        write (iou, *) '# retrun code =', ipar(1),                        &
!     &                 '	convergence rate =', fpar(7)                               
!      endif  ! prinfo         
!!                                                                       
!!     check the error                                                   
!!                                                                       
!      call amux(n,sol,wk,a,ja,ia) 
!      do i = 1, n 
!         wk(n+i) = sol(i) -1.0D0 
!         wk(i) = wk(i) - rhs(i) 
!      enddo 
!      
!      if (prinfo) then         
!        write (iou, *) '# the actual residual norm is', dnrm2(n,wk,1) 
!        write (iou, *) '# the error norm is', dnrm2(n,wk(1+n),1) 
!     endif
!!                                                                       
!      if (iou.ne.6) close(iou) 
                                                                        
      return 
      END         subroutine                                  
