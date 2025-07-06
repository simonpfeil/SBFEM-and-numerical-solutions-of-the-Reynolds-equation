!****************************************************************
!      SORTS A DOUBLE PRECISION VECTOR INDIRECTLY ADDRESSED BY
!      AN INTEGER VECTOR(USES HEAP SORT)
!      THE INDVEC IS REARRANGED SUCH THAT INDVEC(1) ADDRESSES
!      THE LARGEST ELEMENT IN AKEYS, INDVEC(2) ADDRESSES THE
!      NEXT LARGEST ....
!****************************************************************
       SUBROUTINE IHSORT(LEN,INDVEC,AKEYS)
!      THE LENGTH OF THE INTEGER ARRAY
       INTEGER LEN
!      THE INTEGER ARRAY THAT INDIRECTLY ADDRESSES THE D.P. ARRAY
       INTEGER INDVEC(*)
!      THE ARRAY TO BE SORTED
       DOUBLE PRECISION AKEYS(*)

!      THE REST ARE INTERNAL VARIABLES
       INTEGER K, M, LHEAP, RHEAP, MID
       INTEGER X

!****************************************
!      START OF EXECUTABLE STATEMENTS
!****************************************
       IF (LEN.LE.1) RETURN

!      BUILD THE HEAP
       MID = LEN/2
       DO 300 K = MID, 1, -1
          X = INDVEC(K)
          LHEAP = K
          RHEAP = LEN
          M = LHEAP*2
100       CONTINUE
             IF (M.GT.RHEAP) THEN
                INDVEC(LHEAP) = X
                GOTO 200
             ENDIF
             IF (M.LT.RHEAP) THEN
                IF (ABS(AKEYS(INDVEC(M))).GT.ABS(AKEYS(INDVEC(M+1)))) &
                    M = M + 1
             ENDIF
             IF (ABS(AKEYS(X)).LE.ABS(AKEYS(INDVEC(M)))) THEN
                M = RHEAP + 1
             ELSE
                INDVEC(LHEAP) = INDVEC(M)
                LHEAP = M
                M = 2*LHEAP
             ENDIF
             GOTO 100
200       CONTINUE
300    CONTINUE

!      SORT THE HEAP
       DO 600 K = LEN, 2, -1
          X = INDVEC(K)
          INDVEC(K) = INDVEC(1)
          LHEAP = 1
          RHEAP = K-1
          M = 2
400       CONTINUE
             IF (M.GT.RHEAP) THEN
                INDVEC(LHEAP) = X
                GOTO 500
             ENDIF
             IF (M.LT.RHEAP) THEN
                IF (ABS(AKEYS(INDVEC(M))).GT.ABS(AKEYS(INDVEC(M+1))))  &
                    M = M + 1
             ENDIF
             IF (ABS(AKEYS(X)).LE.ABS(AKEYS(INDVEC(M)))) THEN
                M = RHEAP + 1
             ELSE
                INDVEC(LHEAP) = INDVEC(M)
                LHEAP = M
                M = 2*LHEAP
             ENDIF
             GOTO 400
500       CONTINUE
600    CONTINUE

       RETURN
!****************************************************************
!      END OF IHSORT
!****************************************************************
       END SUBROUTINE