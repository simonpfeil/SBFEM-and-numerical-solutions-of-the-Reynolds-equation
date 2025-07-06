      function distdot(n,x,ix,y,iy) 
      integer n, ix, iy 
      !real*8 distdot, x(*), y(*), ddot 
	  real*8 distdot, x(*), y(*)
      !external ddot 
      distdot = ddot(n,x,ix,y,iy) 
      return 
      END             function                              
