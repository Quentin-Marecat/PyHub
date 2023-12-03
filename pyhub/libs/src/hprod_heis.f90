SUBROUTINE HPROD(V,V1)
    USE HUBMOD 
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: V(NSTATES)
    REAL*8,INTENT(OUT) :: V1(NSTATES)
    INTEGER :: I,J,K,L2,NEW,NEW2,POS,BDOWNELEM
    REAL*8 :: AC
    V1 = 0.
    DO J = 1,NORB
        IF (ABS(J_MATRIX(J,J))>1.E-14) THEN
        ! --- (n_i\uparrow(1 - ni_\downarrow) + n_i\downarrow(1 - ni_\uparrow))/2
            DO I = 1,NSUP
                BDOWNELEM = 2**NORB-1 - BUP(I)
                IF (IS_PART(NORB,J,BUP(I)) .AND. .NOT. IS_PART(NORB,J,BDOWNELEM)) THEN
                    V1(I)=V1(I)+V(I)*0.5*J_MATRIX(J,J)
                ELSE IF (.NOT. IS_PART(NORB,J,BUP(I)) .AND. IS_PART(NORB,J,BDOWNELEM)) THEN
                    V1(I)=V1(I)+V(I)*0.5*J_MATRIX(J,J)
                ENDIF
            ENDDO
        ENDIF
        DO K = 1,NORB
            IF ((ABS(J_MATRIX(J,K))>1.E-14) .AND. (J.NE.K)) THEN
                IF (J .NE. K) THEN
                ! --- (-(c^\dagger_k\downarrow c_j\downarrow) * (c^\dagger_j\uparrow  c_k\uparrow)  - (c^\dagger_k\uparrow  c_j\uparrow)  * (c^\dagger_j\downarrow c_k\downarrow))/2
                    DO I = 1,NSUP
                    ! --- KINETIC UP
                        NEW=KINOP(NORB,BUP(I),J,K)
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
                                AC = ANTICOM(NORB,J,BUP(I))
                                AC = AC*ANTICOM(NORB,K,NEW)
!                                AC=ANTICOM(NORB,J,BUP(I),K,NEW)
                                NEW2=KINOP(NORB,BDOWNELEM,K,J)
                                IF (NEW2 .eq. (2**NORB-1-BUP(POS))) THEN
                                    AC = AC*ANTICOM(NORB,J,BDOWNELEM)
                                    AC = AC*ANTICOM(NORB,K,NEW2)
!                                    AC2=ANTICOM(NORB,K,BDOWNELEM,J,NEW2)
                                    V1(POS)=V1(POS)-AC*V(I)*J_MATRIX(J,K)
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
            endif
            DO I = 1,NSUP 
                BDOWNELEM = 2**NORB-1 - BUP(I)
                V1(I)=V1(I)+V(I)*&
                SZ(NORB,J,BUP(I),BDOWNELEM)*SZ(NORB,K,BUP(I),BDOWNELEM)*J_MATRIX(J,K)
            ENDDO
        ENDDO
    ENDDO
    END subroutine
