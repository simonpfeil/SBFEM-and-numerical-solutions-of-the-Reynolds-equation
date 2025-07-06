	subroutine ilu0(n, a, ja, ia, alu, jlu, ju, iw, ierr)
	implicit real*8 (a-h,o-z)
	real*8 a(*), alu(*)
        integer ja(*), ia(*), ju(*), jlu(*), iw(*)
        integer :: n,ierr, ju0,i,ii,js,j,jcol,jf,jm,jrow,tl,jj,jw
!------------------ right preconditioner ------------------------------*
!                    ***   ilu(0) preconditioner.   ***                *
!----------------------------------------------------------------------*
! Note that this has been coded in such a way that it can be used
! with pgmres. Normally, since the data structure of the L+U matrix is
! the same as that the A matrix, savings can be made. In fact with
! some definitions (not correct for general sparse matrices) all we
! need in addition to a, ja, ia is an additional diagonal.
! ILU0 is not recommended for serious problems. It is only provided
! here for comparison purposes.
!-----------------------------------------------------------------------
!
! on entry:
!---------
! n       = dimension of matrix
! a, ja,
! ia      = original matrix in compressed sparse row storage.
!
! on return:
!-----------
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju	  = pointer to the diagonal elements in alu, jlu.
!
! ierr	  = integer indicating error code on return
!	     ierr = 0 --> normal return
!	     ierr = k --> code encountered a zero pivot at step k.
! work arrays:
!-------------
! iw	    = integer work array of length n.
!------------
! IMPORTANT
!-----------
! it is assumed that the the elements in the input matrix are stored
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L sorted by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling ilu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!
!-----------------------------------------------------------------------
        ju0 = n+2
        jlu(1) = ju0
!
! initialize work vector to zero's
!
	do 31 i=1, n
           iw(i) = 0
 31     continue
!
! main loop
!
	do 500 ii = 1, n
           js = ju0
!
! generating row number ii of L and U.
!
           do 100 j=ia(ii),ia(ii+1)-1
!
!     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
!
              jcol = ja(j)
              if (jcol .eq. ii) then
                 alu(ii) = a(j)
                 iw(jcol) = ii
                 ju(ii)  = ju0
              else
                 alu(ju0) = a(j)
                 jlu(ju0) = ja(j)
                 iw(jcol) = ju0
                 ju0 = ju0+1
              endif
 100       continue
           jlu(ii+1) = ju0
           jf = ju0-1
           jm = ju(ii)-1
!
!     exit if diagonal element is reached.
!
           do 150 j=js, jm
              jrow = jlu(j)
              tl = alu(j)*alu(jrow)
              alu(j) = tl
!
!     perform  linear combination
!
              do 140 jj = ju(jrow), jlu(jrow+1)-1
                 jw = iw(jlu(jj))
                 if (jw .ne. 0) alu(jw) = alu(jw) - tl*alu(jj)
 140          continue
 150       continue
!
!     invert  and store diagonal element.
!
           if (alu(ii) .eq. 0.0d0) goto 600
           alu(ii) = 1.0d0/alu(ii)
!
!     reset pointer iw to zero
!
           iw(ii) = 0
           do 201 i = js, jf
 201          iw(jlu(i)) = 0
 500       continue
           ierr = 0
           return
!
!     zero pivot :
!
 600       ierr = ii
!
           return
!------- end-of-ilu0 ---------------------------------------------------
!-----------------------------------------------------------------------
           end subroutine