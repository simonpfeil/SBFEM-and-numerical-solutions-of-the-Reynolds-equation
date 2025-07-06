      subroutine  dcopy(n,dx,incx,dy,incy) 
!                                                                       
!     copies a vector, x, to a vector, y.                               
!     uses unrolled loops for increments equal to one.                  
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision dx(1),dy(1) 
      integer i,incx,incy,ix,iy,m,mp1,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!        code for unequal increments or equal increments                
!          not equal to 1                                               
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        dy(iy) = dx(ix) 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!        code for both increments equal to 1                            
!                                                                       
!                                                                       
!        clean-up loop                                                  
!                                                                       
   20 m = mod(n,7) 
      if( m .eq. 0 ) go to 40 
      do 30 i = 1,m 
        dy(i) = dx(i) 
   30 continue 
      if( n .lt. 7 ) return 
   40 mp1 = m + 1 
      do 50 i = mp1,n,7 
        dy(i) = dx(i) 
        dy(i + 1) = dx(i + 1) 
        dy(i + 2) = dx(i + 2) 
        dy(i + 3) = dx(i + 3) 
        dy(i + 4) = dx(i + 4) 
        dy(i + 5) = dx(i + 5) 
        dy(i + 6) = dx(i + 6) 
   50 continue 
      return 
      END   subroutine                                        
                                                                        
      double precision function ddot(n,dx,incx,dy,incy) 
!                                                                       
!     forms the dot product of two vectors.                             
!     uses unrolled loops for increments equal to one.                  
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision dx(1),dy(1),dtemp 
      integer i,incx,incy,ix,iy,m,mp1,n 
!                                                                       
      ddot = 0.0d0 
      dtemp = 0.0d0 
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!        code for unequal increments or equal increments                
!          not equal to 1                                               
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        dtemp = dtemp + dx(ix)*dy(iy) 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      ddot = dtemp 
      return 
!                                                                       
!        code for both increments equal to 1                            
!                                                                       
!                                                                       
!        clean-up loop                                                  
!                                                                       
   20 m = mod(n,5) 
      if( m .eq. 0 ) go to 40 
      do 30 i = 1,m 
        dtemp = dtemp + dx(i)*dy(i) 
   30 continue 
      if( n .lt. 5 ) go to 60 
   40 mp1 = m + 1 
      do 50 i = mp1,n,5 
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +             &
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue 
   60 ddot = dtemp 
      return 
      END        function                                   
!                                                                       
      double precision function dasum(n,dx,incx) 
!                                                                       
!     takes the sum of the absolute values.                             
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision dx(1),dtemp 
      integer i,incx,m,mp1,n,nincx 
!                                                                       
      dasum = 0.0d0 
      dtemp = 0.0d0 
      if(n.le.0)return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      nincx = n*incx 
      do 10 i = 1,nincx,incx 
        dtemp = dtemp + dabs(dx(i)) 
   10 continue 
      dasum = dtemp 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
!                                                                       
!        clean-up loop                                                  
!                                                                       
   20 m = mod(n,6) 
      if( m .eq. 0 ) go to 40 
      do 30 i = 1,m 
        dtemp = dtemp + dabs(dx(i)) 
   30 continue 
      if( n .lt. 6 ) go to 60 
   40 mp1 = m + 1 
      do 50 i = mp1,n,6 
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2)) &
     &  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))           
   50 continue 
   60 dasum = dtemp 
      return 
      END               function                            
                                                                        
      subroutine daxpy(n,da,dx,incx,dy,incy) 
!                                                                       
!     constant times a vector plus a vector.                            
!     uses unrolled loops for increments equal to one.                  
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision dx(1),dy(1),da 
      integer i,incx,incy,ix,iy,m,mp1,n 
!                                                                       
      if(n.le.0)return 
      if (da .eq. 0.0d0) return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!        code for unequal increments or equal increments                
!          not equal to 1                                               
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        dy(iy) = dy(iy) + da*dx(ix) 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!        code for both increments equal to 1                            
!                                                                       
!                                                                       
!        clean-up loop                                                  
!                                                                       
   20 m = mod(n,4) 
      if( m .eq. 0 ) go to 40 
      do 30 i = 1,m 
        dy(i) = dy(i) + da*dx(i) 
   30 continue 
      if( n .lt. 4 ) return 
   40 mp1 = m + 1 
      do 50 i = mp1,n,4 
        dy(i) = dy(i) + da*dx(i) 
        dy(i + 1) = dy(i + 1) + da*dx(i + 1) 
        dy(i + 2) = dy(i + 2) + da*dx(i + 2) 
        dy(i + 3) = dy(i + 3) + da*dx(i + 3) 
   50 continue 
      return 
      END        subroutine                                   
      double precision function dnrm2 ( n, dx, incx) 
      integer          next ,n, incx, nn, i, j
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one , float
      data   zero, one /0.0d0, 1.0d0/ 
!                                                                       
!     euclidean norm of the n-vector stored in dx() with storage        
!     increment incx .                                                  
!     if    n .le. 0 return with result = 0.                            
!     if n .ge. 1 then incx must be .ge. 1                              
!                                                                       
!           c.l.lawson, 1978 jan 08                                     
!                                                                       
!     four phase method     using two built-in constants that are       
!     hopefully applicable to all machines.                             
!         cutlo = maximum of  dsqrt(u/eps)  over all known machines.    
!         cuthi = minimum of  dsqrt(v)      over all known machines.    
!     where                                                             
!         eps = smallest no. such that eps + 1. .gt. 1.                 
!         u   = smallest positive no.   (underflow limit)               
!         v   = largest  no.            (overflow  limit)               
!                                                                       
!     brief outline of algorithm..                                      
!                                                                       
!     phase 1    scans zero components.                                 
!     move to phase 2 when a component is nonzero and .le. cutlo        
!     move to phase 3 when a component is .gt. cutlo                    
!     move to phase 4 when a component is .ge. cuthi/m                  
!     where m = n for x() real and m = 2*n for complex.                 
!                                                                       
!     values for cutlo and cuthi..                                      
!     from the environmental parameters listed in the imsl converter    
!     document the limiting values are as follows..                     
!     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
!                   univac and dec at 2**(-103)                         
!                   thus cutlo = 2**(-51) = 4.44089e-16                 
!     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.          
!                   thus cuthi = 2**(63.5) = 1.30438e19                 
!     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.             
!                   thus cutlo = 2**(-33.5) = 8.23181d-11               
!     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19                    
!     data cutlo, cuthi / 8.232d-11,  1.304d19 /                        
!     data cutlo, cuthi / 4.441e-16,  1.304e19 /                        
      data cutlo, cuthi / 8.232d-11,  1.304d19 / 
!                                                                       
      if(n .gt. 0) go to 10 
         dnrm2  = zero 
         go to 300 
!                                                                       
   10 assign 30 to next 
      sum = zero 
      nn = n * incx 
!                                                 begin main loop       
      i = 1 
   20    go to next,(30, 50, 70, 110) 
   30 if( dabs(dx(i)) .gt. cutlo) go to 85 
      assign 50 to next 
      xmax = zero 
!                                                                       
!                        phase 1.  sum is zero                          
!                                                                       
   50 if( dx(i) .eq. zero) go to 200 
      if( dabs(dx(i)) .gt. cutlo) go to 85 
!                                                                       
!                                prepare for phase 2.                   
      assign 70 to next 
      go to 105 
!                                                                       
!                                prepare for phase 4.                   
!                                                                       
  100 i = j 
      assign 110 to next 
      sum = (sum / dx(i)) / dx(i) 
  105 xmax = dabs(dx(i)) 
      go to 115 
!                                                                       
!                   phase 2.  sum is small.                             
!                             scale to avoid destructive underflow.     
!                                                                       
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75 
!                                                                       
!                     common code for phases 2 and 4.                   
!                     in phase 4 sum is large.  scale to avoid overflow.
!                                                                       
  110 if( dabs(dx(i)) .le. xmax ) go to 115 
         sum = one + sum * (xmax / dx(i))**2 
         xmax = dabs(dx(i)) 
         go to 200 
!                                                                       
  115 sum = sum + (dx(i)/xmax)**2 
      go to 200 
!                                                                       
!                                                                       
!                  prepare for phase 3.                                 
!                                                                       
   75 sum = (sum * xmax) * xmax 
!                                                                       
!                                                                       
!     for real or d.p. set hitest = cuthi/n                             
!     for complex      set hitest = cuthi/(2*n)                         
!                                                                       
   85 hitest = cuthi/float( n ) 
!                                                                       
!                   phase 3.  sum is mid-range.  no scaling.            
!                                                                       
      do 95 j =i,nn,incx 
      if(dabs(dx(j)) .ge. hitest) go to 100 
   95    sum = sum + dx(j)**2 
      dnrm2 = dsqrt( sum ) 
      go to 300 
!                                                                       
  200 continue 
      i = i + incx 
      if ( i .le. nn ) go to 20 
!                                                                       
!              end of main loop.                                        
!                                                                       
!              compute square root and adjust for scaling.              
!                                                                       
      dnrm2 = xmax * dsqrt(sum) 
  300 continue 
      return 
      END             function                              
                                                                        
      subroutine  dscal(n,da,dx,incx) 
!     scales a vector by a constant.                                    
!     uses unrolled loops for increment equal to one.                   
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision da,dx(1) 
      integer i,incx,m,mp1,n,nincx 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      nincx = n*incx 
      do 10 i = 1,nincx,incx 
        dx(i) = da*dx(i) 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
!                                                                       
!        clean-up loop                                                  
!                                                                       
   20 m = mod(n,5) 
      if( m .eq. 0 ) go to 40 
      do 30 i = 1,m 
        dx(i) = da*dx(i) 
   30 continue 
      if( n .lt. 5 ) return 
   40 mp1 = m + 1 
      do 50 i = mp1,n,5 
        dx(i) = da*dx(i) 
        dx(i + 1) = da*dx(i + 1) 
        dx(i + 2) = da*dx(i + 2) 
        dx(i + 3) = da*dx(i + 3) 
        dx(i + 4) = da*dx(i + 4) 
   50 continue 
      return 
      END                subroutine                           
                                                                        
      subroutine  dswap (n,dx,incx,dy,incy) 
!                                                                       
!     interchanges two vectors.                                         
!     uses unrolled loops for increments equal one.                     
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision dx(1),dy(1),dtemp 
      integer i,incx,incy,ix,iy,m,mp1,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!       code for unequal increments or equal increments not equal       
!         to 1                                                          
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        dtemp = dx(ix) 
        dx(ix) = dy(iy) 
        dy(iy) = dtemp 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!       code for both increments equal to 1                             
!                                                                       
!                                                                       
!       clean-up loop                                                   
!                                                                       
   20 m = mod(n,3) 
      if( m .eq. 0 ) go to 40 
      do 30 i = 1,m 
        dtemp = dx(i) 
        dx(i) = dy(i) 
        dy(i) = dtemp 
   30 continue 
      if( n .lt. 3 ) return 
   40 mp1 = m + 1 
      do 50 i = mp1,n,3 
        dtemp = dx(i) 
        dx(i) = dy(i) 
        dy(i) = dtemp 
        dtemp = dx(i + 1) 
        dx(i + 1) = dy(i + 1) 
        dy(i + 1) = dtemp 
        dtemp = dx(i + 2) 
        dx(i + 2) = dy(i + 2) 
        dy(i + 2) = dtemp 
   50 continue 
      return 
      END          subroutine                                   
                                                                        
      integer function idamax(n,dx,incx) 
!                                                                       
!     finds the index of element having max. absolute value.            
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision dx(1),dmax 
      integer i,incx,ix,n 
!                                                                       
      idamax = 0 
      if( n .lt. 1 ) return 
      idamax = 1 
      if(n.eq.1)return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      ix = 1 
      dmax = dabs(dx(1)) 
      ix = ix + incx 
      do 10 i = 2,n 
         if(dabs(dx(ix)).le.dmax) go to 5 
         idamax = i 
         dmax = dabs(dx(ix)) 
    5    ix = ix + incx 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 dmax = dabs(dx(1)) 
      do 30 i = 2,n 
         if(dabs(dx(i)).le.dmax) go to 30 
         idamax = i 
         dmax = dabs(dx(i)) 
   30 continue 
      return 
      END            function                               
!                                                                       
      subroutine  drot (n,dx,incx,dy,incy,c,s) 
!                                                                       
!     applies a plane rotation.                                         
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision dx(1),dy(1),dtemp,c,s 
      integer i,incx,incy,ix,iy,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!       code for unequal increments or equal increments not equal       
!         to 1                                                          
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        dtemp = c*dx(ix) + s*dy(iy) 
        dy(iy) = c*dy(iy) - s*dx(ix) 
        dx(ix) = dtemp 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!       code for both increments equal to 1                             
!                                                                       
   20 do 30 i = 1,n 
        dtemp = c*dx(i) + s*dy(i) 
        dy(i) = c*dy(i) - s*dx(i) 
        dx(i) = dtemp 
   30 continue 
      return 
      END         subroutine                                    
!                                                                       
      subroutine drotg(da,db,c,s) 
!                                                                       
!     construct givens plane rotation.                                  
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision da,db,c,s,roe,scale,r,z 
!                                                                       
      roe = db 
      if( dabs(da) .gt. dabs(db) ) roe = da 
      scale = dabs(da) + dabs(db) 
      if( scale .ne. 0.0d0 ) go to 10 
         c = 1.0d0 
         s = 0.0d0 
         r = 0.0d0 
         go to 20 
   10 r = scale*dsqrt((da/scale)**2 + (db/scale)**2) 
      r = dsign(1.0d0,roe)*r 
      c = da/r 
      s = db/r 
   20 z = 1.0d0 
      if( dabs(da) .gt. dabs(db) ) z = s 
      if( dabs(db) .ge. dabs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c 
      da = r 
      db = z 
      return 
      END     subroutine                                        
!                                                                       
      subroutine  ccopy(n,cx,incx,cy,incy) 
!                                                                       
!     copies a vector, x, to a vector, y.                               
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      complex cx(1),cy(1) 
      integer i,incx,incy,ix,iy,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!        code for unequal increments or equal increments                
!          not equal to 1                                               
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        cy(iy) = cx(ix) 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!        code for both increments equal to 1                            
!                                                                       
   20 do 30 i = 1,n 
        cy(i) = cx(i) 
   30 continue 
      return 
      END          subroutine                                   
      subroutine  cscal(n,ca,cx,incx) 
!                                                                       
!     scales a vector by a constant.                                    
!     jack dongarra, linpack,  3/11/78.                                 
!                                                                       
      complex ca,cx(1) 
      integer i,incx,n,nincx 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      nincx = n*incx 
      do 10 i = 1,nincx,incx 
        cx(i) = ca*cx(i) 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 do 30 i = 1,n 
        cx(i) = ca*cx(i) 
   30 continue 
      return 
      END      subroutine                                       
!                                                                       
      subroutine  csrot (n,cx,incx,cy,incy,c,s) 
!                                                                       
!     applies a plane rotation, where the cos and sin (c and s) are real
!     and the vectors cx and cy are complex.                            
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      complex cx(1),cy(1),ctemp 
      real c,s 
      integer i,incx,incy,ix,iy,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!       code for unequal increments or equal increments not equal       
!         to 1                                                          
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        ctemp = c*cx(ix) + s*cy(iy) 
        cy(iy) = c*cy(iy) - s*cx(ix) 
        cx(ix) = ctemp 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!       code for both increments equal to 1                             
!                                                                       
   20 do 30 i = 1,n 
        ctemp = c*cx(i) + s*cy(i) 
        cy(i) = c*cy(i) - s*cx(i) 
        cx(i) = ctemp 
   30 continue 
      return 
      END      subroutine                                       
      subroutine  cswap (n,cx,incx,cy,incy) 
!                                                                       
!     interchanges two vectors.                                         
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      complex cx(1),cy(1),ctemp 
      integer i,incx,incy,ix,iy,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!       code for unequal increments or equal increments not equal       
!         to 1                                                          
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        ctemp = cx(ix) 
        cx(ix) = cy(iy) 
        cy(iy) = ctemp 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!       code for both increments equal to 1                             
   20 do 30 i = 1,n 
        ctemp = cx(i) 
        cx(i) = cy(i) 
        cy(i) = ctemp 
   30 continue 
      return 
      END           subroutine                                  
      subroutine  csscal(n,sa,cx,incx) 
!                                                                       
!     scales a complex vector by a real constant.                       
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      complex cx(1) 
      real sa 
      integer i,incx,n,nincx 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      nincx = n*incx 
      do 10 i = 1,nincx,incx 
        cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i))) 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 do 30 i = 1,n 
        cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i))) 
   30 continue 
      return 
      END       subroutine                                      
