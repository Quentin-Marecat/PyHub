SUBROUTINE OPERATOR_BASIS()
    USE BASISMOD
    USE OPEMOD
    USE FUNCMOD
    USE HDF5
    IMPLICIT NONE 
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER :: I,J,L1,NEW,POS,A,A2,I2,B2, NEW2,POS2,L2, SPINK
    REAL*8 :: VAL2, VAL3, VAL2_

    DO I = 1,NB_OPE 
        BUP_ = BU_
        POSU_ = POSU
        BDOWN_ = BD_
        POSD_ = POSD
        VAL2 = 1.
        IF (ABS(COEFF(I)) .GT. 1.E-14) THEN
            DO J = 1,MAX_LEN_OPE
                IF ((OPE(I,J) .NE. 0) .AND. (ABS(VAL2)>1.E-14)) THEN 
                    IF (OPE(I,J) .EQ. 1) THEN
                        CALL C_DAGGER_C_BASIS(SITEI(I,J),SITEJ(I,J),SPINI(I,J),POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 2) THEN
                        CALL NI_BASIS(SITEI(I,J),SPINI(I,J),POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 3) THEN
                        CONTINUE
                    ELSE IF (OPE(I,J) .EQ. 4) THEN
                        CALL NIM1_BASIS(SITEI(I,J),SPINI(I,J),POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 5) THEN
                        CALL SZ_BASIS(SITEI(I,J),BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 6) THEN
                        SPINK=2
                        CALL C_DAGGER_C_UNRES_BASIS(SITEI(I,J),SITEI(I,J),SPINK,POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 7) THEN
                        SPINK=1
                        CALL C_DAGGER_C_UNRES_BASIS(SITEI(I,J),SITEI(I,J),SPINK,POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 8) THEN
                        VAL2_=1.
                        CALL SZ_BASIS(SITEI(I,J),BUP_,BDOWN_,VAL2_)
                        CALL SZ_BASIS(SITEJ(I,J),BUP_,BDOWN_,VAL2_)
                        PSI_OUT((POSD_-1)*NSUP+POSU_) = PSI_OUT((POSD_-1)*NSUP+POSU_) + COEFF(I)*VAL2_
                        CALL SISJ_BASIS(SITEI(I,J),SITEJ(I,J),POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 9) THEN
                        CALL C_DAGGER_C_UNRES_BASIS(SITEI(I,J),SITEJ(I,J),SPINJ(I,J),POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 10) THEN
                        CALL C_DAGGER_BASIS(SITEI(I,J),SPINI(I,J),POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    ELSE IF (OPE(I,J) .EQ. 11) THEN
                        CALL C_BASIS(SITEI(I,J),SPINI(I,J),POSU_,POSD_,BUP_,BDOWN_,VAL2)
                    else
                        CONTINUE
                    ENDIF
                ENDIF
            ENDDO
            PSI_OUT((POSD_-1)*NSUP+POSU_) = PSI_OUT((POSD_-1)*NSUP+POSU_) + COEFF(I)*VAL2
        ENDIF
    ENDDO
END SUBROUTINE


SUBROUTINE C_DAGGER_C_BASIS(SI,SJ,SPI,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SJ,SPI
    REAL*8 :: VAL
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER :: L1,NEW
    IF (SPI.EQ.1) THEN
        NEW=KINOP(BUP_,SJ,SI)
        IF (NEW .NE. 0) THEN
            DO L1 = 1,NSUP 
                IF (NEW==BUP(L1)) THEN
                    POSU_=L1
                    EXIT 
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,SJ,BUP_,SI,NEW)
            BUP_ = NEW
        ELSE 
            VAL = 0.
        ENDIF
    ELSE
        NEW=KINOP(BDOWN_,SJ,SI)
        IF (NEW.NE. 0) THEN
            DO L1 = 1,NSDOWN 
                IF (NEW==BDOWN(L1)) THEN
                    POSD_=L1
                    EXIT 
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,SJ,BDOWN_,SI,NEW)
            BDOWN_ = NEW
        ELSE 
            VAL = 0.
        ENDIF
    ENDIF
END subroutine


SUBROUTINE NI_BASIS(SI,SPI,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SPI
    REAL*8 :: VAL
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER*4 :: A
    A=ISHFT(1,SI-1)
    IF (SPI .EQ. 1) THEN
        IF (IAND(BUP_,A) .NE. A) VAL = 0.
    ELSE
        IF (IAND(BDOWN_,A) .NE. A) VAL = 0.
    ENDIF
END subroutine


SUBROUTINE NIM1_BASIS(SI,SPI,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SPI
    REAL*8 :: VAL
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER*4 :: A
    A=ISHFT(1,SI-1)
    IF (SPI .EQ. 1) THEN
        IF (IAND(BUP_,A) .EQ. A) VAL = 0.
    ELSE
        IF (IAND(BDOWN_,A) .EQ. A) VAL = 0.
    ENDIF
END subroutine


SUBROUTINE SZ_BASIS(SI,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI
    REAL*8 :: VAL2_,VAL
    INTEGER*4 :: BUP_,BDOWN_
    INTEGER*4 :: A
    VAL2_=0.
    A=ISHFT(1,SI-1)
    IF (IAND(BUP_,A) .EQ. A) VAL2_=VAL2_+0.5
    IF (IAND(BDOWN_,A) .EQ. A) VAL2_=VAL2_-0.5
    VAL = VAL*VAL2_
END subroutine



SUBROUTINE SISJ_BASIS(SI,SJ,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SJ
    REAL*8 :: VAL,VAL_
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER*4 :: A, SPIN
    INTEGER :: L2,L3,NEW, NEW2
    IF (SI==SJ) THEN
        A=ISHFT(1,SJ-1)
        ! --- (n_i\uparrow(1 - ni_\downarrow) + n_i\downarrow(1 - ni_\uparrow))/2
        IF (IAND(BUP_,A) .EQ. A .AND. IAND(BDOWN_,A) .NE. A) THEN
            VAL=VAL/2
        ELSE IF (IAND(BUP_,A) .NE. A .AND. IAND(BDOWN_,A) .EQ. A) THEN
            VAL=VAL/2
        ELSE 
            VAL = 0.
        ENDIF
    ELSE IF (SI .NE. SJ) THEN
        VAL_ = VAL
        ! --- (-(c^\dagger_k\downarrow c_j\downarrow) * (c^\dagger_j\uparrow  c_k\uparrow)  - (c^\dagger_k\uparrow  c_j\uparrow)  * (c^\dagger_j\downarrow c_k\downarrow))/2
        SPIN = 1
        CALL C_DAGGER_C_BASIS(SI,SJ,SPIN,POSU_,POSD_,BUP_,BDOWN_,VAL_)
        SPIN = 2
        CALL C_DAGGER_C_BASIS(SJ,SI,SPIN,POSU_,POSD_,BUP_,BDOWN_,VAL_)
        IF (ABS(VAL_)<1.E-14) THEN
            VAL_ = VAL
            BUP_ = BU_ 
            BDOWN_ = BD_
            SPIN = 2
            CALL C_DAGGER_C_BASIS(SI,SJ,SPIN,POSU_,POSD_,BUP_,BDOWN_,VAL_)
            SPIN = 1
            CALL C_DAGGER_C_BASIS(SJ,SI,SPIN,POSU_,POSD_,BUP_,BDOWN_,VAL_)
        ENDIF
        VAL = -VAL_/2
    ENDIF
    RETURN

END subroutine



SUBROUTINE C_DAGGER_C_UNRES_BASIS(SI,SJ,SPJ,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SJ,SPJ
    REAL*8 :: VAL
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER*4 :: A, A2, I2, B2
    INTEGER*4 :: L1,L2,L3,NEW, NEW2
    A=ISHFT(1,SI-1)
    A2=ISHFT(1,SJ-1)
    IF (SPJ.EQ.1) THEN
        DO I2 =NORB,1,-1 ! --- ANTICOMMUTATION FOR ELECTRON UP
            B2=ISHFT(1,I2-1)
            IF (IAND(BDOWN_,B2)==B2) VAL=-VAL
        ENDDO
        ! --- REMOVE UP
        IF (IAND(BUP_,A2) .EQ. A2) THEN
            NEW = BUP_ - A2
            DO L1 = 1,NSUP 
                IF (NEW==BUP(L1)) THEN
                    POSU_=L1
                    EXIT
                ENDIF
            ENDDO
            IF (SJ<NORB) THEN
                DO I2 =NORB,SJ+1,-1
                    B2=ISHFT(1,I2-1)
                    IF (IAND(BUP_,B2)==B2) VAL=-VAL
                ENDDO
            ENDIF
            BUP_ = NEW
            IF (IAND(BDOWN_,A) .EQ. 0) THEN
                NEW2=BDOWN_ + A
                DO L2 = 1,NSDOWN 
                    IF (NEW2==BDOWN(L2)) THEN
                        POSD_=L2
                        EXIT
                    ENDIF 
                ENDDO
                IF (SI<NORB) THEN
                    DO I2 =NORB,SI+1,-1
                        B2=ISHFT(1,I2-1)
                        IF (IAND(BDOWN_,B2)==B2) VAL=-VAL
                    ENDDO
                ENDIF
                BDOWN_=NEW2
            ELSE
                VAL=0.
            ENDIF
        else
            VAL=0.
        ENDIF
    ELSE
        ! --- REMOVE DOWN
        IF (IAND(BDOWN_,A2) .EQ. A2) THEN
            NEW = BDOWN_-A2
            DO I2 =NORB,1,-1 ! --- ANTICOMMUTATION FOR ELECTRON UP
                B2=ISHFT(1,I2-1)
                IF (IAND(NEW,B2)==B2) VAL=-VAL
            ENDDO
            DO L1 = 1,NSDOWN
                IF (NEW==BDOWN(L1)) THEN
                    POSD_=L1
                    EXIT 
                ENDIF
            ENDDO
            IF (SJ<NORB) THEN
                DO I2 =NORB,SJ+1,-1
                    B2=ISHFT(1,I2-1)
                    IF (IAND(BDOWN_,B2)==B2) VAL=-VAL
                ENDDO
            ENDIF
            BDOWN_ = NEW
            IF (IAND(BUP_,A) .EQ. 0) THEN
                NEW2 = BUP_+A
                DO L2 = 1,NSUP
                    IF (NEW2==BUP(L2)) THEN
                        POSU_=L2
                        EXIT 
                    ENDIF
                ENDDO
                IF (SI<NORB) THEN
                    DO I2 =NORB,SI+1,-1
                        B2=ISHFT(1,I2-1)
                        IF (IAND(BUP_,B2)==B2) VAL=-VAL
                    ENDDO
                ENDIF
                BUP_=NEW2
            ELSE 
                VAL = 0.
            ENDIF
        ELSE 
            VAL = 0.
        ENDIF
    ENDIF
    RETURN
end subroutine



SUBROUTINE C_DAGGER_BASIS(SI,SPI,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SPI
    REAL*8 :: VAL
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER*4 :: A, A2, I2, B2
    INTEGER*4 :: L1,L2,L3,NEW, NEW2
    A=ISHFT(1,SI-1)
    IF (SPI.EQ.1) THEN
        DO I2 =NORB,1,-1 ! --- ANTICOMMUTATION FOR ELECTRON UP
            B2=ISHFT(1,I2-1)
            IF (IAND(BDOWN_,B2)==B2) VAL=-VAL
        ENDDO
        ! --- ADD UP
        IF (IAND(BUP_,A) .EQ. 0) THEN
            NEW = BUP_ + A
            DO L1 = 1,NSUP 
                IF (NEW==BUP(L1)) THEN
                    POSU_=L1
                    EXIT
                ENDIF
            ENDDO
            IF (SI<NORB) THEN
                DO I2 =NORB,SI+1,-1
                    B2=ISHFT(1,I2-1)
                    IF (IAND(BUP_,B2)==B2) VAL=-VAL
                ENDDO
            ENDIF
            BUP_ = NEW
        else
            VAL=0.
        ENDIF
    ELSE
        ! --- ADD DOWN
        IF (IAND(BDOWN_,A) .EQ. 0) THEN
            NEW = BDOWN_+A
            DO L1 = 1,NSDOWN
                IF (NEW==BDOWN(L1)) THEN
                    POSD_=L1
                    EXIT 
                ENDIF
            ENDDO
            IF (SI<NORB) THEN
                DO I2 =NORB,SI+1,-1
                    B2=ISHFT(1,I2-1)
                    IF (IAND(BDOWN_,B2)==B2) VAL=-VAL
                ENDDO
            ENDIF
            BDOWN_ = NEW
        ELSE 
            VAL = 0.
        ENDIF
    ENDIF
    RETURN
end subroutine


SUBROUTINE C_BASIS(SI,SPI,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SPI
    REAL*8 :: VAL
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER*4 :: A, A2, I2, B2
    INTEGER*4 :: L1,L2,L3,NEW, NEW2
    A=ISHFT(1,SI-1)
    IF (SPI.EQ.1) THEN
        DO I2 =NORB,1,-1 ! --- ANTICOMMUTATION FOR ELECTRON UP
            B2=ISHFT(1,I2-1)
            IF (IAND(BDOWN_,B2)==B2) VAL=-VAL
        ENDDO
        ! --- REMOVE UP
        IF (IAND(BUP_,A) .EQ. A) THEN
            NEW = BUP_ - A
            DO L1 = 1,NSUP 
                IF (NEW==BUP(L1)) THEN
                    POSU_=L1
                    EXIT
                ENDIF
            ENDDO
            IF (SI<NORB) THEN
                DO I2 =NORB,SI+1,-1
                    B2=ISHFT(1,I2-1)
                    IF (IAND(BUP_,B2)==B2) VAL=-VAL
                ENDDO
            ENDIF
            BUP_ = NEW
        else
            VAL=0.
        ENDIF
    ELSE
        ! --- REMOVE DOWN
        IF (IAND(BDOWN_,A) .EQ. A) THEN
            NEW = BDOWN_-A
            DO L1 = 1,NSDOWN
                IF (NEW==BDOWN(L1)) THEN
                    POSD_=L1
                    EXIT 
                ENDIF
            ENDDO
            IF (SI<NORB) THEN
                DO I2 =NORB,SI+1,-1
                    B2=ISHFT(1,I2-1)
                    IF (IAND(BDOWN_,B2)==B2) VAL=-VAL
                ENDDO
            ENDIF
            BDOWN_ = NEW
        ELSE 
            VAL = 0.
        ENDIF
    ENDIF
    RETURN
end subroutine