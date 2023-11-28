MODULE HUBMOD
    ! ---------------------------------- !
    ! --- GLOBAL DATA FOR HUB SOLVER --- !
    ! ---------------------------------- !
    IMPLICIT NONE

    ! --- INTEGER*4 --- !
    INTEGER*4 :: MAXLCZ, NDEG, NB_COMP_STATES, EXC_STATE,NW, N2, NBLCZ, IS_ULOC
    INTEGER*4 :: NORB,NUP,NDOWN,NSUP,NSDOWN,NSTATES,NELEC
    INTEGER*4,ALLOCATABLE :: BUP(:),BDOWN(:)
    INTEGER*4 :: NORB2COMP, NB_POLES(2,2)

    ! --- FLOAT*8 --- !
    REAL*8, ALLOCATABLE ::  UL(:,:,:,:),T(:,:),J_MATRIX(:,:), UPARTIAL(:,:)
    REAL*8, ALLOCATABLE ::  H(:,:), ULDIAG(:), PSI(:)
    REAL*8 :: ACC_LCZ,TP,Z,E0, U_FLOAT
    REAL*8, ALLOCATABLE :: W(:),VEC(:,:),VAL(:)

    ! --- LOGICAL --- !
    LOGICAL ::  RENORM,DO_BASIS,DO_SOLVE,DO_RQ, DO_SPGF, ISLANCZOS, DO_RQ_TB
    logical :: STORE_H

    ! --- variable protected  --- !
    INTEGER*4 ::  NSUP_,NSDOWN_, NELEC_, NSTATES_, NUP_, NDOWN_
    INTEGER*4, ALLOCATABLE :: BUP_(:), BDOWN_(:)
    REAL*8 :: E0_,SZ_
    REAL*8, ALLOCATABLE :: VAL_(:),VEC_(:,:), LCZVEC_(:,:), PSI_(:,:)

    CONTAINS


    SUBROUTINE ORDER(NORB,NSTATES,Q,L)
        IMPLICIT NONE 
        INTEGER*4, INTENT(IN) :: NORB
        INTEGER*4 :: NSTATES
        REAL*8,allocatable :: Q(:,:),L(:)
        INTEGER*4 :: I,J, TO_DEL(NSTATES), NB_DEL
        REAL*8, ALLOCATABLE :: Q_O(:,:), L_O(:)
        LOGICAL, ALLOCATABLE :: MK(:)
        NB_DEL=0
        TO_DEL(1)=0
        DO I = 1,NSTATES
            IF (NORM2(Q(I,:))<1.E-14) THEN 
                NB_DEL=NB_DEL+1
                TO_DEL(NB_DEL)=I
            ENDIF 
        ENDDO 
        ALLOCATE(Q_O(NSTATES-NB_DEL,NORB),L_O(NSTATES-NB_DEL))
        NB_DEL=1
        J=1
        DO I = 1,NSTATES
            IF (I==TO_DEL(NB_DEL)) THEN 
                NB_DEL=NB_DEL+1
                CONTINUE
            ELSE 
                Q_O(J,:)=Q(I,:)
                L_O(J)=L(I)
                J=J+1
            ENDIF
        ENDDO
        NSTATES=NSTATES-(NB_DEL-1)
        DEALLOCATE(Q,L)
        ALLOCATE(Q(NSTATES,NORB),L(NSTATES),MK(NSTATES))
        MK=.TRUE.
        DO I = 1,NSTATES
            J=MINLOC(L_O,NSTATES,MK)
            L(I)=L_O(J)
            Q(I,:)=Q_O(J,:)
            MK(J)=.FALSE.
        ENDDO
        DEALLOCATE(Q_O,L_O)
        RETURN 
    END SUBROUTINE

    SUBROUTINE SET_POLE(N1,POLE1,N2,POLE2,NB)
        implicit NONE 
        INTEGER, INTENT(IN) :: N1,N2
        REAL*8, INTENT(IN) :: POLE1(N1,2)
        REAL*8 :: POLE2(N2,2)
        INTEGER, INTENT(OUT) :: NB
        INTEGER :: I,J,K
        NB = 0
        DO I = 1,N1
            IF (I < N1) THEN
                IF (ABS(POLE1(I,1)) .LT. 1.E-14 .AND. ABS(POLE1(I+1,1)) .LT. 1.E-14) EXIT
            ENDIF
            if (ABS(POLE1(I,2)) .GT. 1.E-12) THEN
                DO J = 1,N2-1
                    IF (ABS(POLE1(I,1)-POLE2(J,1)) .LT. 1.E-12) THEN
                        POLE2(J,2) = POLE2(J,2) + POLE1(I,2)
                        EXIT
                    ELSE IF ((POLE1(I,1) .GT. POLE2(J,1) .AND. POLE1(I,1) .LT. POLE2(J+1,1))) THEN
                        DO K = N2-1,J+1,-1
                            POLE2(K+1,1) = POLE2(K,1)
                            POLE2(K+1,2) = POLE2(K,2)
                        ENDDO
                        POLE2(J+1,1) = POLE1(I,1)
                        POLE2(J+1,2) = POLE1(I,2)
                        EXIT 
                    ELSE IF (ABS(POLE2(J,1)) .LT. 1.E-14 .AND. ABS(POLE2(J+1,1)) .LT. 1.E-14) THEN
                        POLE2(J,1) = POLE1(I,1)
                        POLE2(J,2) = POLE1(I,2)
                        EXIT
                    ELSE IF (POLE1(I,1) .LT. POLE2(J,1) .AND. J .EQ. 1) THEN
                        DO K = N2-1,J,-1
                            POLE2(K+1,1) = POLE2(K,1)
                            POLE2(K+1,2) = POLE2(K,2)
                        ENDDO
                        POLE2(J,1) = POLE1(I,1)
                        POLE2(J,2) = POLE1(I,2)
                        EXIT
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        DO I = 1,N2
            IF (ABS(POLE2(I,1)) .LT. 1.E-14 .AND. ABS(POLE2(I+1,1)) .LT. 1.E-14) THEN
                EXIT
            ELSE IF (ABS(POLE2(I,2)) .LT. 1.E-12) THEN
                DO K = I,N2-1
                    POLE2(K,1) = POLE2(K+1,1)
                    POLE2(K,2) = POLE2(K+1,2)
                ENDDO
            ELSE
                NB = NB + 1
            ENDIF
        ENDDO
        RETURN 
    END subroutine SET_POLE

END MODULE HUBMOD
