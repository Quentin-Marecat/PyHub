MODULE FUNCMOD
    ! ---
    ! --- ALL FUNCTIONS USED
    ! ---
    CONTAINS
    
    recursive function BINOM(n, k) result(res)
    integer*4, intent(in) :: n, k
    integer*4 :: res
    
    if (k == 0 .or. k == n) then
    res = 1
    else
    res = BINOM(n-1, k-1) + BINOM(n-1, k)
    end if
    end function BINOM

    INTEGER*4 FUNCTION KINOP(STATEI,POSI,POSF)
        INTEGER*4,INTENT(IN)::STATEI,POSI,POSF
        INTEGER*4,SAVE::A,B
        KINOP=0
        A=ISHFT(1,POSI-1)
        B=ISHFT(1,POSF-1)
        IF(IAND(STATEI,A)==A .AND. IAND(STATEI,B)==0) KINOP = STATEI - A + B
        IF (IAND(STATEI,A)==A .AND. POSI==POSF) KINOP=STATEI
        RETURN
    END FUNCTION

    REAL*8 FUNCTION ANTICOM(NO,POSI,STATEI,POSF,STATEF)
        INTEGER*4,INTENT(IN)::STATEF,STATEI,NO,POSF,POSI
        INTEGER*4,SAVE:: A2,I2
        ANTICOM=1.
        IF (POSI/=POSF) THEN
            DO I2 =NO,1,-1
                A2=ISHFT(1,I2-1)
                IF (I2==POSI) THEN
                    EXIT
                ELSE IF (IAND(STATEI,A2)==A2) THEN
                    ANTICOM=-ANTICOM
                ENDIF
            ENDDO
            DO I2 =NO,1,-1
                A2=ISHFT(1,I2-1)
                IF (I2==POSF) THEN
                    EXIT
                ELSE IF (IAND(STATEF,A2)==A2) THEN
                    ANTICOM=-ANTICOM
                ENDIF
            ENDDO
        ENDIF
        RETURN
    END FUNCTION


    REAL*8 FUNCTION N_J(STATE,POS)
        INTEGER, INTENT(IN) :: STATE,POS
        INTEGER :: A
        N_J = 0.
        A = ISHFT(1,POS-1)
        IF(IAND(STATE,A) .EQ. A) N_J = 1.
        RETURN 
    END function N_J


    REAL*8 FUNCTION ANTICOM2(NO,STATEI,POSF)
        INTEGER*4,INTENT(IN)::STATEI,NO,POSF
        INTEGER*4 :: A,I
        ANTICOM2=1.
        DO I =NO,1,-1
            A=ISHFT(1,I-1)
            IF (I==POSF) THEN
                EXIT
            ELSE IF (IAND(STATEI,A)==A) THEN
                ANTICOM2=-ANTICOM2
            ENDIF
        ENDDO
        RETURN
    END FUNCTION

END MODULE FUNCMOD