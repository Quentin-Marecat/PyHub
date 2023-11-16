SUBROUTINE TEMPERATURE(NVAL,VAL,T,CONV,P,N)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NVAL 
    REAL*8, INTENT(IN) :: VAL(NVAL), CONV , T
    INTEGER, INTENT(OUT) :: N
    REAL*8, INTENT(OUT) :: P(NVAL)
    REAL*8 :: ZERO, V,V2
    INTEGER :: I
    ZERO=1.E-14
    P=0.
    IF (T .GE. ZERO) THEN
        V=0.
        DO I = 1,NVAL
            V2=EXP((VAL(1)-VAL(I))/T)
            IF (V2 .GE. CONV ) THEN
                V=V+V2
            ELSE 
                GOTO 10
            ENDIF
        ENDDO
10       N=I
        IF (N>NVAL) N = NVAL
        DO I = 1,N 
            P(I)=EXP((VAL(1)-VAL(I))/T)/SQRT(V)
        ENDDO
    ELSE 
        DO I = 1,NVAL
            V2=ABS(VAL(1)-VAL(I))
            IF (V2 .GE. CONV ) GOTO 20
        ENDDO
20       N=I-1
        IF (N>NVAL) N = NVAL
        DO I = 1,N 
            P(I)=1/SQRT(DBLE(N))
        ENDDO
    ENDIF
    RETURN
END SUBROUTINE
