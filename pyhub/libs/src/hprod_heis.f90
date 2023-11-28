SUBROUTINE HPROD(V,V1)
    USE HUBMOD 
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: V(NSTATES)
    REAL*8,INTENT(OUT) :: V1(NSTATES)
    INTEGER :: I,J,K,L2,NEW,NEW2,POS,A,BDOWNELEM
    REAL*8 :: AC,AC2
    V1 = 0.
    DO J = 1,NORB
        A=ISHFT(1,J-1)
        IF (ABS(J_MATRIX(J,J))>1.E-14) THEN
        ! --- (n_i\uparrow(1 - ni_\downarrow) + n_i\downarrow(1 - ni_\uparrow))/2
            DO I = 1,NSUP
                BDOWNELEM = 2**NORB-1 - BUP(I)
                IF (IAND(BUP(I),A) .EQ. A .AND. IAND(BDOWNELEM,A) .EQ. 0) THEN
                    V1(I)=V1(I)+V(I)*0.5*J_MATRIX(J,K)
                ELSE IF (IAND(BUP(I),A) .EQ. 0 .AND. IAND(BDOWNELEM,A) .EQ. A) THEN
                    V1(I)=V1(I)+V(I)*0.5*J_MATRIX(J,K)
                ENDIF
            ENDDO
        ENDIF
        DO K = 1,NORB
            IF ((ABS(J_MATRIX(J,K))>1.E-14) .AND. (J.NE.K)) THEN
                IF (J .NE. K) THEN
                ! --- (-(c^\dagger_k\downarrow c_j\downarrow) * (c^\dagger_j\uparrow  c_k\uparrow)  - (c^\dagger_k\uparrow  c_j\uparrow)  * (c^\dagger_j\downarrow c_k\downarrow))/2
                    DO I = 1,NSUP
                    ! --- KINETIC UP
                        NEW=KINOP(BUP(I),J,K)
                        IF (NEW>0) THEN
                            POS=0
                            DO L2 = 1,NSUP 
                                IF (NEW==BUP(L2)) THEN
                                POS=L2
                                    EXIT 
                                ENDIF
                            ENDDO
                            IF (POS>0) THEN
                                BDOWNELEM = 2**NORB-1 - BUP(I)
                                AC=ANTICOM(NORB,J,BUP(I),K,NEW)
                                NEW2=KINOP(BDOWNELEM,K,J)
                                IF (NEW2 .eq. (2**NORB-1-BUP(POS))) THEN
                                    AC2=ANTICOM(NORB,K,BDOWNELEM,J,NEW2)
                                    V1(POS)=V1(POS)-AC*AC2*V(I)*J_MATRIX(J,K)
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
                DO I = 1,NSUP 
                    BDOWNELEM = 2**NORB-1 - BUP(I)
                    V1(I)=V1(I)+V(I)*&
                    (N_J(BUP(I),J)-N_J(BDOWNELEM,J))*(N_J(BUP(I),K)-N_J(BDOWNELEM,K))*0.25*J_MATRIX(J,K)
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    END subroutine
