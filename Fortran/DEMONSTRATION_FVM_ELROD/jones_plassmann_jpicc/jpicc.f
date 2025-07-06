!****************************************************************
!     PERFORM JONES/PLASSMANN INCOMPLETE CHOLESKI:COLUMN ORIENTED
!****************************************************************
      INTEGER FUNCTION JPICC(N,DIAG,A,IA,IE,JA,TA,ITCOL,IFIRST,LIST)
      !INTEGER FUNCTION JPICC(N,DIAG,A,IA,JA,TA,ITCOL,IFIRST,LIST)
      
!     IF THE FACTORIZATION WAS P.D. THEN 0 IS RETURNED
!     OTHERWISE A NEGATIVE VALUE IS RETURNED THAT INDICATES
!     THE COLUMN NUMBER WHERE A NEGATIVE DIAGONAL WAS ENCOUNTERED

!     THE ORDER OF THE MATRIX
!     INPUT ONLY
      INTEGER N
!     THE DIAGONALS OF A
!     INPUT/OUTPUT
      DOUBLE PRECISION DIAG(*)
!     THE OFF-DIAGONALS OF A
!     INPUT/OUTPUT
      DOUBLE PRECISION A(*)
!*     POINTERS TO THE COLUMNS OF A
!*     IA(K) IS THE INDEX IN A() AND JA() WHERE COLUMN K STARTS
!*     ONLY THE STRICTLY LOWER TRIANGLE OF A IS STORED
!*     IA IS LENGTH N+1 (POSITION N+1 INDICATES WHERE COLUMN N+1
!*     WOULD START IF IT EXISTED)
!*     INPUT
!
!     POINTERS TO THE COLUMNS OF A
!     IA(K) IS THE INDEX IN A() AND JA() WHERE COLUMN K STARTS
!     ONLY THE STRICTLY LOWER TRIANGLE OF A IS STORED
!     IA IS LENGTH N      
!     INPUT      
      INTEGER IA(*)
!     Index Ende der Zeilen, auch length(IE)=N
!     INPUT
      INTEGER IE(*)
!     THE ROW NUMBERS OF THE OFF-DIAGONALS OF A
!     INPUT/OUTPUT      
      INTEGER JA(*)
!     A TEMPORARY WORK VECTOR OF LENGTH N TO KEEP THE CURRENT COLUMN
!     CONTENTS DESTROYED
      DOUBLE PRECISION TA(*)
!     A TEMPORARY WORK VECTOR OF LENGTH N TO KEEP TRACK OF THE ROW
!     VALUES IN THE CURRENT COLUMN
!     CONTENTS DESTROYED
      INTEGER ITCOL(*)
!     IFIRST(J) POINTS TO THE NEXT VALUE IN COLUMN J TO USE (LENGTH N)
!     IFIRST ALSO HAS A DUAL USE.  AT STEP K, ONLY THE FIRST K-1 
!     ELEMENTS ARE USED FOR THE ABOVE PURPOSE.  FOR THE LAST N-K 
!     ELEMENTS, IFIRST(J) INDICATES IF IF A NONZERO VALUE EXISTS IN 
!     POSITION J OF COLUMN K.
!     CONTENTS DESTROYED
      INTEGER IFIRST(*)
!     LIST(J) POINTS TO A LINKED LIST OF COLUMNS THAT WILL UPDATE
!     COLUMN J (LENGTH N)
!     CONTENTS DESTROYED
      INTEGER LIST(*)

!     SUBROUTINES USED 
!      EXTERNAL IBSORT, DBSORT, DHSORT

!     VARIABLES USED
      INTEGER ISJ, IEJ, ISK, IEK
      INTEGER I, J, K
      INTEGER TALEN
      INTEGER ROW, COUNT
      DOUBLE PRECISION LVAL
      INTEGER IPTR
 
!****************************************
!      START OF EXECUTABLE STATEMENTS
!****************************************

      DO 100 J = 1, N
        IFIRST(J) = 0
        LIST(J) = 0
100    CONTINUE

!     LOOP OVER ALL COLUMNS
      DO 900 K = 1,N
!        LOAD COLUMN K INTO TA
         TALEN = 0
         ISK = IA(K)
         !IEK = IA(K+1)-1
         IEK= IE(K)
         DO 200 J = ISK, IEK
           ROW = JA(J)
           TA(ROW) = A(J)
           TALEN = TALEN + 1
           ITCOL(TALEN) = ROW
           IFIRST(ROW) = 1
200       CONTINUE

!        MAKE SURE THE DIAGONAL OF K IS OKAY AND THEN TAKE THE SQRT
         IF (DIAG(K).LE.0.0D0) THEN
           DIAG(K) = -1.0D0
           GOTO 1000
         ENDIF
         DIAG(K) = SQRT(DIAG(K))

!        UPDATE COLUMN K USING THE PREVIOUS COLUMNS
         J = LIST(K)
300      CONTINUE
         IF (J.EQ.0) GOTO 500
           ISJ = IFIRST(J)
           !IEJ = IA(J+1)-1
           IEJ = IE(j)
           LVAL = A(ISJ)
           ISJ = ISJ + 1
           IF (ISJ.LT.IEJ) THEN
             IFIRST(J) = ISJ
             IPTR = J
             J = LIST(J)
             LIST(IPTR) = LIST(JA(ISJ))
             LIST(JA(ISJ)) = IPTR
           ELSE
             J = LIST(J)
           ENDIF
           DO 400 I = ISJ, IEJ
             ROW = JA(I)
             IF (IFIRST(ROW).NE.0) THEN
               TA(ROW) = TA(ROW) - LVAL*A(I)
             ELSE
               IFIRST(ROW) = 1
               TALEN = TALEN + 1
               ITCOL(TALEN) = ROW
               TA(ROW) = - LVAL*A(I)
             ENDIF
400        CONTINUE
         GOTO 300
500      CONTINUE

!        UPDATE REMAINING DIAGONALS USING COLUMN K
         DO 600 J = 1, TALEN
           ROW = ITCOL(J)
           TA(ROW) = TA(ROW)/DIAG(K)
           DIAG(ROW) = DIAG(ROW) - TA(ROW)*TA(ROW)
600      CONTINUE

!        FIND THE LARGEST ELEMENTS IN COLUMN K NOW
         COUNT = MIN(IEK-ISK+1,TALEN)
         CALL IBSORT(TALEN,COUNT,TA,ITCOL)
         IF (COUNT.LT.20) THEN
           CALL DBSORT(COUNT,ITCOL)
         ELSE
           CALL DHSORT(COUNT,ITCOL)
         ENDIF

!        PUT THE LARGEST ELEMENTS BACK INTO THE SPARSE DATA STRUCTURE
         COUNT = 1
         DO 700 J = ISK, IEK
           A(J) = TA(ITCOL(COUNT))
           JA(J) = ITCOL(COUNT)
           COUNT = COUNT + 1
700      CONTINUE

!        IFIRST AND LIST KEEP TRACK OF WHERE IN COLUMN K WE ARE
         IF (ISK.LT.IEK) THEN
           IPTR = JA(ISK)
           LIST(K) = LIST(IPTR)
           LIST(IPTR) = K
           IFIRST(K) = ISK
         ENDIF

         DO 800 J = 1, TALEN
           IFIRST(ITCOL(J)) = 0
800      CONTINUE

900   CONTINUE

      JPICC = 0
      RETURN

!     IF AN ERROR OCCURED, RETURN A NEGATIVE VALUE
1000  CONTINUE
      JPICC = -K
      RETURN
!****************************************************************
!     END OF JPICC
!****************************************************************
      END FUNCTION