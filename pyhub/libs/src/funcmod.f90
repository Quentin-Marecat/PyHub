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

    INTEGER*4 FUNCTION KINOP(NORB,STATEI,POSI,POSF)
        INTEGER*4,INTENT(IN)::STATEI,POSI,POSF,NORB
        INTEGER*4,SAVE::A,B
        KINOP=0
        A=ISHFT(1,NORB+1-POSI-1)
        B=ISHFT(1,NORB+1-POSF-1)
        IF(IAND(STATEI,A)==A .AND. IAND(STATEI,B)==0) KINOP = STATEI - A + B
        IF (IAND(STATEI,A)==A .AND. POSI==POSF) KINOP=STATEI
        RETURN
    END FUNCTION

    REAL*8 FUNCTION ANTICOM(NORB,POS,STATE)
        INTEGER*4,INTENT(IN):: NORB,STATE,POS
        INTEGER*4,SAVE:: A2,I2
        ANTICOM=1.
        IF (POS>1) THEN
!            DO I2 =NORB,POS+1,-1
            DO I2 =1,POS-1
                A2=ISHFT(1,NORB+1-I2-1)
                IF (IAND(STATE,A2)==A2) ANTICOM=-ANTICOM
            ENDDO
        ENDIF
        RETURN
    END FUNCTION

    ! REAL*8 FUNCTION ANTICOM(NO,POSI,STATEI,POSF,STATEF)
    !     INTEGER*4,INTENT(IN)::STATEF,STATEI,NO,POSF,POSI
    !     INTEGER*4,SAVE:: A2,I2
    !     ANTICOM=1.
    !     IF (POSI/=POSF) THEN
    !         IF (POSI>1) THEN
    !             DO I2 =1,POSI-1
    !                 A2=ISHFT(1,I2-1)
    !                 IF (IAND(STATEI,A2)==A2) ANTICOM=-ANTICOM
    !             ENDDO
    !         ENDIF
    !         IF (POSF>1) THEN
    !             DO I2 =1,POSF-1
    !                 A2=ISHFT(1,I2-1)
    !                 IF (IAND(STATEF,A2)==A2) ANTICOM=-ANTICOM
    !             ENDDO
    !         ENDIF
    !     ENDIF
    !     RETURN
    ! END FUNCTION

    INTEGER*4 FUNCTION CREATION(NORB,POS,STATE)
        INTEGER, INTENT(IN) :: NORB,STATE,POS
        CREATION = STATE + ISHFT(1,NORB+1-POS-1)
        RETURN 
    END function CREATION

    INTEGER*4 FUNCTION ANNIHILATION(NORB,POS,STATE)
        INTEGER, INTENT(IN) :: NORB,STATE,POS
        ANNIHILATION = STATE - ISHFT(1,NORB+1-POS-1)
        RETURN 
    END function ANNIHILATION

    REAL*8 FUNCTION SZ(NORB,POS,STATE_UP,STATE_DOWN)
        INTEGER, INTENT(IN) :: NORB,POS,STATE_UP,STATE_DOWN
        INTEGER :: A
        SZ=0
        A = ISHFT(1,NORB+1-POS-1)
        IF(IAND(STATE_UP,A) .EQ. A) SZ = SZ+0.5
        IF(IAND(STATE_DOWN,A) .EQ. A) SZ = SZ-0.5
        RETURN 
    END function SZ


    LOGICAL FUNCTION IS_PART(NORB,POS,STATE)
        INTEGER, INTENT(IN) :: NORB,STATE,POS
        INTEGER :: A
        IS_PART = .FALSE.
        A = ISHFT(1,NORB+1-POS-1)
        IF(IAND(STATE,A) .EQ. A) IS_PART = .TRUE.
        RETURN 
END function IS_PART


END MODULE FUNCMOD