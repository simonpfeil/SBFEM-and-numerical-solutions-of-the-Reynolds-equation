!****************************************************************
!      GET THE K LARGEST NONZEROES IN AKEYS INDIRECTLY ADDRESSED 
!      BY INDVEC; UPON EXIT THE FIRST K ELEMENTS IN INDVEC WILL
!      CONTAIN THE INDICES OF THE K LARGEST ELEMENTS IN AKEYS
!****************************************************************
       SUBROUTINE IBSORT(N,K,AKEYS,INDVEC)

!      THE LENGTH OF THE INTEGER VECTOR
       INTEGER N
!      THE NUMBER WANTED
       INTEGER K
!      THE DOUBLE PRECISION VECTOR TO BE SORTED
       DOUBLE PRECISION AKEYS(*)
!      THE INTEGER VECTOR ASSOCIATED WITH AKEYS
!      INDVEC(I) GIVES THE POSITION IN AKEYS OF THE ITH ELEMENT
       INTEGER INDVEC(*)
 
!      THE REST ARE INTERNAL VARIABLES
       INTEGER I,J
       INTEGER ITEMP, CURPTR, RIGHT, LEFT
       DOUBLE PRECISION CURMIN, NEWVAL, CURVAL, LVAL

       !EXTERNAL IHSORT
 
!****************************************
!      START OF EXECUTABLE STATEMENTS
!****************************************

!      IF THE LIST IS SMALL OR THE NUMBER REQUIRED IS 0 THEN 
!      RETURN
       IF ((N.LE.1).OR.(K.LE.0)) RETURN

!      HEAP SORT THE FIRST K ELEMENTS OF THE VECTOR
       CALL IHSORT(K,INDVEC,AKEYS)

!      LOOP THROUGH THE REST OF THE VECTOR AND FIND ANY ELEMENTS
!      THAT ARE LARGER THAN ANY OF THE FIRST K ELEMENTS
       CURMIN = ABS(AKEYS(INDVEC(K)))
       DO 400 I = K+1, N
          ITEMP = INDVEC(I)
          NEWVAL = ABS(AKEYS(ITEMP))
          IF (NEWVAL.GT.CURMIN) THEN
!           FIND POSITION FOR NEW VALUE
            LEFT = 1
            LVAL = ABS(AKEYS(INDVEC(1)))
            IF (NEWVAL.GT.LVAL) THEN
              CURPTR = 1
              GOTO 200
            ENDIF
            RIGHT = K
            CURPTR = (K+1)/2
100         CONTINUE            
            IF (RIGHT.GT.LEFT+1) THEN
              CURVAL = ABS(AKEYS(INDVEC(CURPTR)))
              IF (CURVAL.LT.NEWVAL) THEN
                RIGHT = CURPTR
              ELSE
                LEFT = CURPTR
                LVAL = CURVAL
              ENDIF
              CURPTR = (RIGHT+LEFT)/2
              GOTO 100
            ENDIF
            CURPTR = RIGHT

!           SHIFT SORTED VALUES AND INSERT NEW VALUE
200         CONTINUE
            INDVEC(I) = INDVEC(K)
            DO 300 J = K, CURPTR+1, -1
              INDVEC(J) = INDVEC(J-1)
300         CONTINUE
            INDVEC(CURPTR) = ITEMP
            CURMIN = ABS(AKEYS(INDVEC(K)))
          ENDIF
400    CONTINUE
 
       RETURN
!****************************************************************
!      END OF IBSORT
!****************************************************************
      END SUBROUTINE