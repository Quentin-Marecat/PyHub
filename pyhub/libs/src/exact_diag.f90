
SUBROUTINE EXACT_DIAG()
    USE BASISMOD
    USE HUBMOD
    IMPLICIT NONE
    INTEGER*4 :: I
    REAL*8,ALLOCATABLE :: V(:),V2(:)
    IF (ALLOCATED(VAL)) DEALLOCATE(VAL)
    IF (ALLOCATED(VEC)) DEALLOCATE(VEC)
    IF (ALLOCATED(H)) DEALLOCATE(H)
    ALLOCATE(H(NSTATES,NSTATES))
    H=0.
    ALLOCATE(V(NSTATES),V2(NSTATES))
    V=0.
    DO I = 1,NSTATES
        V(I)=1.
        CALL HPROD(V,V2)
        H(:,I)=V2
        V(I)=0.
    ENDDO
    DEALLOCATE(V,V2)
    ALLOCATE(VAL(NSTATES),VEC(NSTATES,NSTATES))
    CALL DIAGMAT(NSTATES,H,VAL,VEC)
!    DEALLOCATE(H)
    NDEG=1 
    DO I = 2,NSTATES 
        IF (ABS(VAL(I)-VAL(1))<1.E-14) THEN 
            NDEG=NDEG+1
        ELSE 
            EXIT
        ENDIF 
    ENDDO
    RETURN
END SUBROUTINE


SUBROUTINE EXACT_DIAG_()
    USE BASISMOD
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4 :: I,J,K,L1,NEW,POS
    INTEGER*4 :: A,A2,M,N,POS2,NEW2,L2,U,D
    REAL*8 :: AC,AC2
    IF (ALLOCATED(VAL)) DEALLOCATE(VAL)
    IF (ALLOCATED(VEC)) DEALLOCATE(VEC)
    IF (ALLOCATED(H)) DEALLOCATE(H)
    ALLOCATE(H(NSTATES,NSTATES))
    H=0.
    DO J = 1,NORB
        A=ISHFT(1,J-1)
        IF (IS_ULOC) THEN
            DO U = 1,NSUP 
                IF (IAND(BUP(U),A) .EQ. A) THEN
                    DO D = 1,NSDOWN
                        IF (IAND(BDOWN(D),A) .EQ. A) &
                        H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+ULDIAG(J)
                    ENDDO 
                ENDIF 
            ENDDO 
        ENDIF
        DO K=1,NORB
            A2=ISHFT(1,K-1)
            ! --- ELECTRON UP
            IF (ABS(T(J,K)) .GE. 1.E-14) THEN
                DO I = 1,NSUP
                    ! --- KINETIC UP
                    NEW=KINOP(BUP(I),J,K)
                    IF (NEW>0) THEN
                        POS=0
                        DO L1 = 1,NSUP 
                            IF (NEW==BUP(L1)) THEN
                                POS=L1
                                EXIT 
                            ENDIF
                        ENDDO
                        IF (POS>0) THEN
                            AC=ANTICOM(NORB,J,BUP(I),K,NEW)
                            DO L1 = 1,NSDOWN 
                                H((L1-1)*NSUP+POS,(L1-1)*NSUP+I) = H((L1-1)*NSUP+POS,(L1-1)*NSUP+I)+AC*T(J,K)
                                IF (ABS(SZ) .LT. 1.E-14) H((POS-1)*NSUP+L1,(I-1)*NSUP+L1) = &
                                H((POS-1)*NSUP+L1,(I-1)*NSUP+L1)+AC*T(J,K)
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
            IF (ABS(SZ) .GE. 1.E-14) THEN
                ! --- ELECTRON DOWN 
                IF (ABS(T(J,K)) .GE. 1.E-14) THEN
                    DO L1 = 1,NSDOWN
                    ! --- KINETIC DOWN
                        NEW=KINOP(BDOWN(L1),J,K)
                        IF (NEW>0) THEN
                            POS=0
                            DO I = 1,NSDOWN 
                                IF (NEW==BDOWN(I)) THEN
                                    POS=I
                                    EXIT 
                                ENDIF
                            ENDDO
                            IF (POS>0) THEN
                                AC=ANTICOM(NORB,J,BDOWN(L1),K,NEW)
                                DO I = 1,NSUP 
                                    H((POS-1)*NSUP+I,(L1-1)*NSUP+I)=H((POS-1)*NSUP+I,(L1-1)*NSUP+I)+AC*T(J,K)
                                ENDDO
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
            IF (.NOT. IS_ULOC) THEN
                DO M = 1,NORB
                    DO N=1,NORB
                        IF (ABS(UL(J,K,M,N)) .GE. 1.E-14) THEN
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
                                        AC=ANTICOM(NORB,J,BUP(I),K,NEW)
                                        DO L1 = 1,NSDOWN
                                            ! --- KINETIC DOWN
                                            NEW2=KINOP(BDOWN(L1),M,N)
                                            IF (NEW2>0) THEN
                                                POS2=0
                                                DO L2 = 1,NSDOWN
                                                    IF (NEW2==BDOWN(L2)) THEN
                                                        POS2=L2
                                                        EXIT 
                                                    ENDIF
                                                ENDDO
                                                IF (POS2>0) THEN
                                                    AC2=ANTICOM(NORB,M,BDOWN(L1),N,NEW2)
                                                    H((POS2-1)*NSUP+POS,(L1-1)*NSUP+I) = H((POS2-1)*NSUP+POS,(L1-1)*NSUP+I)+&
                                                    AC*AC2*UL(J,K,M,N)
                                                ENDIF
                                            ENDIF
                                        ENDDO
                                    ENDIF
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    ALLOCATE(VAL(NSTATES),VEC(NSTATES,NSTATES))
    CALL DIAGMAT(NSTATES,H,VAL,VEC)
!    DEALLOCATE(H)
    NDEG=1 
    DO I = 2,NSTATES 
        IF (ABS(VAL(I)-VAL(1))<1.E-14) THEN 
            NDEG=NDEG+1
        ELSE 
            EXIT
        ENDIF 
    ENDDO
    RETURN
END SUBROUTINE