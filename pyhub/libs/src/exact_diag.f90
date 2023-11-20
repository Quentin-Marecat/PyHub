SUBROUTINE EXACT_DIAG_()
    USE BASISMOD
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER :: I,J,K,L,U,D,TMP,NEW,NEW2,POS,POS2,A
    REAL*8 :: AC, AC2
    IF (ALLOCATED(VAL)) DEALLOCATE(VAL)
    IF (ALLOCATED(VEC)) DEALLOCATE(VEC)
    IF (ALLOCATED(H)) DEALLOCATE(H)
    ALLOCATE(H(NSTATES,NSTATES))
    H=0.
    IF (NUP .NE. NDOWN) THEN
        DO U = 1,NSUP
            DO I = 1,NORB
                A = ISHFT(1,I-1)
                DO J=1,NORB
                    DO D = 1,NSDOWN
                        IF (I .EQ. J .AND. ABS(J_MATRIX(I,I)) .GE. 1.E-14) THEN
                            IF (IAND(BUP(U),A) .EQ. A .AND. IAND(BDOWN(D),A) .EQ. 0) &
                            H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+0.5*J_MATRIX(I,I)
                            IF (IAND(BUP(U),A) .EQ. 0 .AND. IAND(BDOWN(D),A) .EQ. A) &
                            H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+0.5*J_MATRIX(I,I)
                        ENDIF
                        H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+J_MATRIX(I,J)*&
                        (SZ_J(BUP(U),I)-SZ_J(BDOWN(D),I))*(SZ_J(BUP(U),J)-SZ_J(BDOWN(D),J))
                    ENDDO
                    NEW=KINOP(BUP(U),I,J)
                    IF (NEW>0) THEN
                        POS=0
                        DO TMP = 1,NSUP 
                            IF (NEW==BUP(TMP)) THEN
                                POS=TMP
                                EXIT 
                            ENDIF
                        ENDDO
                        IF (POS>0) THEN
                            AC=ANTICOM(NORB,I,BUP(U),J,NEW)
                            DO D = 1,NSDOWN 
                                H((D-1)*NSUP+POS,(D-1)*NSUP+U)=H((D-1)*NSUP+POS,(D-1)*NSUP+U)+AC*T(I,J)
                                IF (I .NE. J .AND. ABS(J_MATRIX(I,J)) .GE. 1.E-14) THEN
                                    NEW2=KINOP(BDOWN(D),J,I)
                                    IF (NEW2>0) THEN
                                        POS2=0
                                        DO TMP = 1,NSDOWN 
                                            IF (NEW2==BDOWN(TMP)) THEN
                                                POS2=TMP
                                                EXIT 
                                            ENDIF
                                        ENDDO
                                        IF (POS2>0) THEN
                                            AC2=ANTICOM(NORB,J,BDOWN(D),I,NEW2)
                                            H((POS2-1)*NSUP+POS,(D-1)*NSUP+U)=H((POS2-1)*NSUP+POS,(D-1)*NSUP+U)-&
                                            AC*AC2*J_MATRIX(I,J)
                                        ENDIF
                                    ENDIF
                                ENDIF
                                IF (.NOT. IS_ULOC) THEN
                                    DO K = 1,NORB
                                        DO L = 1,NORB
                                            NEW2=KINOP(BDOWN(D),K,L)
                                            IF (NEW2>0) THEN
                                                POS2=0
                                                DO TMP = 1,NSDOWN
                                                    IF (NEW2==BDOWN(TMP)) THEN
                                                        POS2=TMP
                                                        EXIT 
                                                    ENDIF
                                                ENDDO
                                                IF (POS2>0) THEN
                                                    AC2=ANTICOM(NORB,K,BDOWN(D),L,NEW2)
                                                    H((POS2-1)*NSUP+POS,(D-1)*NSUP+U)=&
                                                    H((POS2-1)*NSUP+POS,(D-1)*NSUP+U)+AC2*AC*UL(I,J,K,L)
                                                ENDIF
                                            ENDIF
                                        ENDDO
                                    ENDDO
                                ELSE IF (I .EQ. J .AND. IAND(BDOWN(D),A).EQ.A) THEN
                                    H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+ULDIAG(I)
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        DO D = 1,NSDOWN
            DO I = 1,NORB
                DO J=1,NORB
                    NEW=KINOP(BDOWN(D),I,J)
                    IF (NEW>0) THEN
                        POS=0
                        DO TMP = 1,NSDOWN
                            IF (NEW==BDOWN(TMP)) THEN
                                POS=TMP
                                EXIT 
                            ENDIF
                        ENDDO
                        IF (POS>0) THEN
                            AC=ANTICOM(NORB,I,BDOWN(D),J,NEW)
                            DO U = 1,NSUP
                                H((POS-1)*NSUP+U,(D-1)*NSUP+U)=H((POS-1)*NSUP+U,(D-1)*NSUP+U)+AC*T(I,J)
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ELSE
        DO U = 1,NSUP
            DO I = 1,NORB
                A = ISHFT(1,I-1)
                DO J=1,NORB
                    DO D = 1,NSDOWN
                        IF (I .EQ. J .AND. ABS(J_MATRIX(I,I)) .GE. 1.E-14) THEN
                            IF (IAND(BUP(U),A) .EQ. A .AND. IAND(BDOWN(D),A) .EQ. 0) &
                            H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+0.5*J_MATRIX(I,I)
                            IF (IAND(BUP(U),A) .EQ. 0 .AND. IAND(BDOWN(D),A) .EQ. A) &
                            H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+0.5*J_MATRIX(I,I)
                        ENDIF
                        H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+J_MATRIX(I,J)*&
                        (SZ_J(BUP(U),I)-SZ_J(BDOWN(D),I))*(SZ_J(BUP(U),J)-SZ_J(BDOWN(D),J))
                    ENDDO
                    NEW=KINOP(BUP(U),I,J)
                    IF (NEW>0) THEN
                        POS=0
                        DO TMP = 1,NSUP 
                            IF (NEW==BUP(TMP)) THEN
                                POS=TMP
                                EXIT 
                            ENDIF
                        ENDDO
                        IF (POS>0) THEN
                            AC=ANTICOM(NORB,I,BUP(U),J,NEW)
                            DO D = 1,NSDOWN 
                                H((D-1)*NSUP+POS,(D-1)*NSUP+U)=H((D-1)*NSUP+POS,(D-1)*NSUP+U)+AC*T(I,J)
                                H((POS-1)*NSUP+D,(U-1)*NSUP+D)=H((POS-1)*NSUP+D,(U-1)*NSUP+D)+AC*T(I,J)
                                IF (I .NE. J .AND. ABS(J_MATRIX(I,J)) .GE. 1.E-14) THEN
                                    NEW2=KINOP(BDOWN(D),J,I)
                                    IF (NEW2>0) THEN
                                        POS2=0
                                        DO TMP = 1,NSDOWN 
                                            IF (NEW2==BDOWN(TMP)) THEN
                                                POS2=TMP
                                                EXIT 
                                            ENDIF
                                        ENDDO
                                        IF (POS2>0) THEN
                                            AC2=ANTICOM(NORB,J,BDOWN(D),I,NEW2)
                                            H((POS2-1)*NSUP+POS,(D-1)*NSUP+U)=&
                                            H((POS2-1)*NSUP+POS,(D-1)*NSUP+U)-AC*AC2*J_MATRIX(I,J)
                                        ENDIF
                                    ENDIF
                                ENDIF
                                IF (.NOT. IS_ULOC) THEN
                                    DO K = 1,NORB
                                        DO L = 1,NORB
                                            NEW2=KINOP(BDOWN(D),K,L)
                                            IF (NEW2>0) THEN
                                                POS2=0
                                                DO TMP = 1,NSDOWN
                                                    IF (NEW2==BDOWN(TMP)) THEN
                                                        POS2=TMP
                                                        EXIT 
                                                    ENDIF
                                                ENDDO
                                                IF (POS2>0) THEN
                                                    AC2=ANTICOM(NORB,K,BDOWN(D),L,NEW2)
                                                    H((POS2-1)*NSUP+POS,(D-1)*NSUP+U)=&
                                                    H((POS2-1)*NSUP+POS,(D-1)*NSUP+U)+AC2*AC*UL(I,J,K,L)
                                                ENDIF
                                            ENDIF
                                        ENDDO
                                    ENDDO
                                ELSE IF (I .EQ. J .AND. IAND(BDOWN(D),A).EQ.A) THEN
                                    H((D-1)*NSUP+U,(D-1)*NSUP+U)=H((D-1)*NSUP+U,(D-1)*NSUP+U)+ULDIAG(I)
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDIF 
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