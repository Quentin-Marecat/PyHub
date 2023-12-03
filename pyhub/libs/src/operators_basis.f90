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
        NEW=KINOP(NORB,BUP_,SJ,SI)
        IF (NEW .NE. 0) THEN
            DO L1 = 1,NSUP 
                IF (NEW==BUP(L1)) THEN
                    POSU_=L1
                    EXIT 
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,SJ,BUP_)
            VAL=VAL*ANTICOM(NORB,SI,NEW)
            BUP_ = NEW
        ELSE 
            VAL = 0.
        ENDIF
    ELSE
        NEW=KINOP(NORB,BDOWN_,SJ,SI)
        IF (NEW.NE. 0) THEN
            DO L1 = 1,NSDOWN 
                IF (NEW==BDOWN(L1)) THEN
                    POSD_=L1
                    EXIT 
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,SJ,BDOWN_)
            VAL=VAL*ANTICOM(NORB,SI,NEW)
            BDOWN_ = NEW
        ELSE 
            VAL = 0.
        ENDIF
    ENDIF
END subroutine


SUBROUTINE NI_BASIS(SI,SPI,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE FUNCMOD
    USE OPEMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SPI
    REAL*8 :: VAL
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    IF (SPI .EQ. 1) THEN
        IF (.NOT. IS_PART(NORB,SI,BUP_)) VAL = 0.
    ELSE
        IF (.NOT. IS_PART(NORB,SI,BDOWN_)) VAL = 0.
    ENDIF
END subroutine


SUBROUTINE NIM1_BASIS(SI,SPI,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE FUNCMOD
    USE OPEMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SPI
    REAL*8 :: VAL
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    IF (SPI .EQ. 1) THEN
        IF (IS_PART(NORB,SI,BUP_)) VAL = 0.
    ELSE
        IF (IS_PART(NORB,SI,BDOWN_)) VAL = 0.
    ENDIF
END subroutine


SUBROUTINE SZ_BASIS(SI,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,BUP_,BDOWN_
    REAL*8 :: VAL
    VAL = VAL*SZ(NORB,SI,BUP_,BDOWN_)
    RETURN
END subroutine



SUBROUTINE SISJ_BASIS(SI,SJ,POSU_,POSD_,BUP_,BDOWN_,VAL)
    USE BASISMOD 
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SI,SJ
    REAL*8 :: VAL,VAL_
    INTEGER*4 :: POSU_,POSD_,BUP_,BDOWN_
    INTEGER*4 :: SPIN
    INTEGER :: L2,L3,NEW, NEW2
    IF (SI==SJ) THEN
        ! --- (n_i\uparrow(1 - ni_\downarrow) + n_i\downarrow(1 - ni_\uparrow))/2
        IF (IS_PART(NORB,SJ,BUP_) .AND. .NOT. IS_PART(NORB,SJ,BDOWN_)) THEN
            VAL=VAL/2
        ELSE IF (.NOT. IS_PART(NORB,SJ,BUP_) .AND. IS_PART(NORB,SJ,BDOWN_)) THEN
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
    INTEGER*4 :: I2, B2
    INTEGER*4 :: L1,L2,L3,NEW, NEW2,A
    IF (SPJ.EQ.1) THEN
        ! --- REMOVE UP
        IF (IS_PART(NORB,SJ,BUP_)) THEN
            NEW = ANNIHILATION(NORB,SJ,BUP_)
            DO L1 = 1,NSUP 
                IF (NEW==BUP(L1)) THEN
                    POSU_=L1
                    EXIT
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,SJ,BUP_)
            BUP_ = NEW
            IF (.NOT. IS_PART(NORB,SI,BDOWN_)) THEN
                NEW2=CREATION(NORB,SI,BDOWN_)
                DO L2 = 1,NSDOWN 
                    IF (NEW2==BDOWN(L2)) THEN
                        POSD_=L2
                        EXIT
                    ENDIF 
                ENDDO
                VAL=VAL*ANTICOM(NORB,0,BDOWN_)
                VAL=VAL*ANTICOM(NORB,SI,BDOWN_)
                BDOWN_=NEW2
            ELSE
                VAL=0.
            ENDIF
        else
            VAL=0.
        ENDIF
    ELSE
        ! --- REMOVE DOWN
        IF (IS_PART(NORB,SJ,BDOWN_)) THEN
            NEW = ANNIHILATION(NORB,SJ,BDOWN_)
            DO L1 = 1,NSDOWN
                IF (NEW==BDOWN(L1)) THEN
                    POSD_=L1
                    EXIT 
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,SJ,BDOWN_)
            BDOWN_ = NEW
            IF (.NOT. IS_PART(NORB,SI,BUP_)) THEN
                VAL=VAL*ANTICOM(NORB,0,NEW)
                NEW2 = CREATION(NORB,SI,BUP_)
                DO L2 = 1,NSUP
                    IF (NEW2==BUP(L2)) THEN
                        POSU_=L2
                        EXIT 
                    ENDIF
                ENDDO
                VAL=VAL*ANTICOM(NORB,SI,BUP_)
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
    INTEGER*4 :: I2, B2
    INTEGER*4 :: L1,L2,L3,NEW, NEW2
    IF (SPI.EQ.1) THEN
        ! --- ADD UP
        IF (.NOT. IS_PART(NORB,SI,BUP_)) THEN
            NEW = CREATION(NORB,SI,BUP_)
            DO L1 = 1,NSUP 
                IF (NEW==BUP(L1)) THEN
                    POSU_=L1
                    EXIT
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,0,BDOWN_)
            VAL=VAL*ANTICOM(NORB,SI,BUP_)
            BUP_ = NEW
        else
            VAL=0.
        ENDIF
    ELSE
        ! --- ADD DOWN
        IF (.NOT. IS_PART(NORB,SI,BDOWN_)) THEN
            NEW = CREATION(NORB,SI,BDOWN_)
            DO L1 = 1,NSDOWN
                IF (NEW==BDOWN(L1)) THEN
                    POSD_=L1
                    EXIT 
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,SI,BDOWN_)
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
    INTEGER*4 :: A2, I2, B2
    INTEGER*4 :: L1,L2,L3,NEW, NEW2
    IF (SPI.EQ.1) THEN
        ! --- REMOVE UP
        IF (IS_PART(NORB,SI,BUP_)) THEN
            NEW = ANNIHILATION(NORB,SI,BUP_)
            DO L1 = 1,NSUP 
                IF (NEW==BUP(L1)) THEN
                    POSU_=L1
                    EXIT
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,0,BDOWN_)
            VAL=VAL*ANTICOM(NORB,SI,BUP_)
            BUP_ = NEW
        else
            VAL=0.
        ENDIF
    ELSE
        ! --- REMOVE DOWN
        IF (IS_PART(NORB,SI,BDOWN_)) THEN
            NEW = ANNIHILATION(NORB,SI,BDOWN_)
            DO L1 = 1,NSDOWN
                IF (NEW==BDOWN(L1)) THEN
                    POSD_=L1
                    EXIT 
                ENDIF
            ENDDO
            VAL=VAL*ANTICOM(NORB,SI,BDOWN_)
            BDOWN_ = NEW
        ELSE 
            VAL = 0.
        ENDIF
    ENDIF
    RETURN
end subroutine