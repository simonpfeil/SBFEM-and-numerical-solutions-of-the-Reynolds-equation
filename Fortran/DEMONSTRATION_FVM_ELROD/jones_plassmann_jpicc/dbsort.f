!****************************************************************
!      SORTS AN INTEGER VECTOR (USES BUBBLE SORT)
!      ASCENDING ORDER
!****************************************************************
       SUBROUTINE DBSORT(N,KEYVEC)

!      THE LENGTH OF THE VECTOR
       INTEGER N
!      THE INTEGER VECTOR TO BE SORTED
       INTEGER KEYVEC(*)
 
!      THE REST ARE INTERNAL VARIABLES
       INTEGER I,J
       INTEGER TEMP
 
!****************************************
!      START OF EXECUTABLE STATEMENTS
!****************************************
       DO 200 I = 1, N-1
          DO 100 J = I+1, N
             IF (KEYVEC(I).GT.KEYVEC(J)) THEN
                TEMP = KEYVEC(I)
                KEYVEC(I) = KEYVEC(J)
                KEYVEC(J) = TEMP
             ENDIF
100       CONTINUE
200    CONTINUE
 
       RETURN
!****************************************************************
!      END OF DBSORT
!****************************************************************
       END SUBROUTINE