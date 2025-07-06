!****************************************************************
!      SORTS AN INTEGER VECTOR (USES HEAP SORT)
!      ASCENDING ORDER
!****************************************************************
       SUBROUTINE DHSORT(LEN,KEYS)
!      THE LENGTH OF THE ARRAY
       INTEGER LEN
!      THE ARRAY TO BE SORTED
       INTEGER KEYS(*)

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
          X = KEYS(K)
          LHEAP = K
          RHEAP = LEN
          M = LHEAP*2
100       CONTINUE
             IF (M.GT.RHEAP) THEN
                KEYS(LHEAP) = X
                GOTO 200
             ENDIF
             IF (M.LT.RHEAP) THEN
                IF (KEYS(M) .LT. KEYS(M+1)) M = M+1
             ENDIF
             IF (X.GE.KEYS(M)) THEN
                M = RHEAP + 1
             ELSE
                KEYS(LHEAP) = KEYS(M)
                LHEAP = M
                M = 2*LHEAP
             ENDIF
             GOTO 100
200       CONTINUE
300    CONTINUE

!      SORT THE HEAP
       DO 600 K = LEN, 2, -1
          X = KEYS(K)
          KEYS(K) = KEYS(1)
          LHEAP = 1
          RHEAP = K-1
          M = 2
400       CONTINUE
             IF (M.GT.RHEAP) THEN
                KEYS(LHEAP) = X
                GOTO 500
             ENDIF
             IF (M.LT.RHEAP) THEN
                IF (KEYS(M) .LT. KEYS(M+1)) M = M+1
             ENDIF
             IF (X.GE.KEYS(M)) THEN
                M = RHEAP + 1
             ELSE
                KEYS(LHEAP) = KEYS(M)
                LHEAP = M
                M = 2*LHEAP
             ENDIF
             GOTO 400
500       CONTINUE
600    CONTINUE

       RETURN
!****************************************************************
!      END OF DHSORT
!****************************************************************
       END SUBROUTINE