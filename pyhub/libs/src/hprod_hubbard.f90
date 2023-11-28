SUBROUTINE HPROD_(V,V2)
    USE HUBMOD
    ! ------------------------------------------- !
    ! --- COMPUTE THE APPLICATION OF H OVER V --- !
    ! ------------------------------------------- !
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: V(NSTATES)
    REAL*8,INTENT(OUT) :: V2(NSTATES)
    V2 = 0.
    CALL TPROD(V,V2)
    CALL UPROD(V,V2)
!    CALL JPROD(V,V2)
    RETURN
END SUBROUTINE


SUBROUTINE HPROD(V,V1)
    USE HUBMOD 
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: V(NSTATES)
    REAL*8,INTENT(OUT) :: V1(NSTATES)
    INTEGER :: I,J,K,M,N,L1,L2,NEW,NEW2,POS,POS2,A,U,D
    REAL*8 :: AC,AC2
    V1 = 0.
    DO J = 1,NORB 
        A = ISHFT(1,J-1)
        IF (IS_ULOC .EQ. 2) THEN
            DO U = 1,NSUP 
                IF (IAND(BUP(U),A) .EQ. A) THEN
                    DO D = 1,NSDOWN
                        IF (IAND(BDOWN(D),A) .EQ. A) V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*ULDIAG(J)
                    ENDDO
                ENDIF
            ENDDO 
        ENDIF
        DO K = 1,NORB
            IF ((ABS(T(J,K)) .GE. 1.E-14) .OR. (UPARTIAL(J,K) .GE. 1.E-14)) THEN
                DO U = 1,NSUP 
                    NEW=KINOP(BUP(U),J,K)
                    IF (NEW>0) THEN
                        POS=0
                        DO L1 = 1,NSUP 
                            IF (NEW==BUP(L1)) THEN
                                POS=L1
                                EXIT 
                            ENDIF
                        ENDDO
                        IF (POS>0) THEN
                            AC=ANTICOM(NORB,J,BUP(U),K,NEW)
                            IF (NUP .EQ. NDOWN) THEN
                                DO D = 1,NSDOWN
                                    V1((D-1)*NSUP+POS)=V1((D-1)*NSUP+POS)+AC*V((D-1)*NSUP+U)*T(J,K)
                                    V1((POS-1)*NSUP+D)=V1((POS-1)*NSUP+D)+AC*V((U-1)*NSUP+D)*T(J,K)
                                ENDDO
                            ELSE 
                                DO D = 1,NSDOWN
                                    V1((D-1)*NSUP+POS)=V1((D-1)*NSUP+POS)+AC*V((D-1)*NSUP+U)*T(J,K)
                                ENDDO
                            ENDIF
                            IF (UPARTIAL(J,K) .GE. 1.E-14) THEN
                                DO M = 1,NORB 
                                    DO N = 1,NORB
                                        IF (ABS(UL(J,K,M,N)) .GE. 1.E-14) THEN
                                            DO D = 1,NSDOWN
                                                ! --- KINETIC DOWN
                                                NEW2=KINOP(BDOWN(D),M,N)
                                                IF (NEW2>0) THEN
                                                    POS2=0
                                                    DO L2 = 1,NSDOWN
                                                        IF (NEW2==BDOWN(L2)) THEN
                                                            POS2=L2
                                                            EXIT 
                                                        ENDIF
                                                    ENDDO
                                                    IF (POS2>0) THEN
                                                        AC2=ANTICOM(NORB,M,BDOWN(D),N,NEW2)
                                                        V1((POS2-1)*NSUP+POS) = V1((POS2-1)*NSUP+POS)+&
                                                        AC*AC2*V((D-1)*NSUP+U)*UL(J,K,M,N)
                                                    ENDIF
                                                ENDIF
                                            ENDDO
                                        ENDIF
                                    ENDDO 
                                ENDDO
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
                IF (NUP .NE. NDOWN) THEN
                    DO D = 1,NSDOWN
                        NEW=KINOP(BDOWN(D),J,K)
                        IF (NEW>0) THEN
                            POS=0
                            DO L1 = 1,NSDOWN 
                                IF (NEW==BDOWN(L1)) THEN
                                    POS=L1
                                    EXIT 
                                ENDIF
                            ENDDO
                            IF (POS>0) THEN
                                AC=ANTICOM(NORB,J,BDOWN(D),K,NEW)
                                DO U = 1,NSUP
                                    V1((POS-1)*NSUP+U)=V1((POS-1)*NSUP+U)+AC*V((D-1)*NSUP+U)*T(J,K)
                                ENDDO
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    ENDDO 
    END subroutine

! SUBROUTINE HPROD_(V,V1)
!     USE HUBMOD
!     USE FUNCMOD
!     ! --- MUCH SLOWER THAN HPROD
!     ! ------------------------------------------- !
!     ! --- COMPUTE THE APPLICATION OF H OVER V --- !
!     ! ------------------------------------------- !
!     IMPLICIT NONE
!     REAL*8, INTENT(IN) :: V(NSTATES)
!     REAL*8, INTENT(OUT) :: V1(NSTATES)
!     INTEGER :: I,J,K,L,U,D,TMP,NEW,NEW2,POS,POS2,A
!     REAL*8 :: AC, AC2
!     V1 = 0.
!     IF (NUP .NE. NDOWN) THEN
!         DO U = 1,NSUP
!             DO I = 1,NORB
!                 A = ISHFT(1,I-1)
!                 DO J=1,NORB
!                     DO D = 1,NSDOWN
!                         IF (I .EQ. J .AND. ABS(J_MATRIX(I,I)) .GE. 1.E-14) THEN
!                             IF (IAND(BUP(U),A) .EQ. A .AND. IAND(BDOWN(D),A) .EQ. 0) &
!                             V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*0.5*J_MATRIX(I,I)
!                             IF (IAND(BUP(U),A) .EQ. 0 .AND. IAND(BDOWN(D),A) .EQ. A) &
!                             V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*0.5*J_MATRIX(I,I)
!                         ENDIF
!                         ! V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*J_MATRIX(I,J)*&
!                         ! (SZ_J(BUP(U),I)-SZ_J(BDOWN(D),I))*(SZ_J(BUP(U),J)-SZ_J(BDOWN(D),J))
!                     ENDDO
!                     NEW=KINOP(BUP(U),I,J)
!                     IF (NEW>0) THEN
!                         POS=0
!                         DO TMP = 1,NSUP 
!                             IF (NEW==BUP(TMP)) THEN
!                                 POS=TMP
!                                 EXIT 
!                             ENDIF
!                         ENDDO
!                         IF (POS>0) THEN
!                             AC=ANTICOM(NORB,I,BUP(U),J,NEW)
!                             DO D = 1,NSDOWN 
!                                 V1((D-1)*NSUP+POS)=V1((D-1)*NSUP+POS)+AC*V((D-1)*NSUP+U)*T(I,J)
!                                 IF (I .NE. J .AND. ABS(J_MATRIX(I,J)) .GE. 1.E-14) THEN
!                                     NEW2=KINOP(BDOWN(D),J,I)
!                                     IF (NEW2>0) THEN
!                                         POS2=0
!                                         DO TMP = 1,NSDOWN 
!                                             IF (NEW2==BDOWN(TMP)) THEN
!                                                 POS2=TMP
!                                                 EXIT 
!                                             ENDIF
!                                         ENDDO
!                                         IF (POS2>0) THEN
!                                             AC2=ANTICOM(NORB,J,BDOWN(D),I,NEW2)
!                                             V1((POS2-1)*NSUP+POS)=V1((POS2-1)*NSUP+POS)-&
!                                             AC*AC2*V((D-1)*NSUP+U)*J_MATRIX(I,J)
!                                         ENDIF
!                                     ENDIF
!                                 ENDIF
!                                 IF (.NOT. IS_ULOC) THEN
!                                     DO K = 1,NORB
!                                         DO L = 1,NORB
!                                             IF (ABS(UL(I,J,K,L)) .GT. 1.E-14) THEN
!                                                 NEW2=KINOP(BDOWN(D),K,L)
!                                                 IF (NEW2>0) THEN
!                                                     POS2=0
!                                                     DO TMP = 1,NSDOWN
!                                                         IF (NEW2==BDOWN(TMP)) THEN
!                                                             POS2=TMP
!                                                             EXIT 
!                                                         ENDIF
!                                                     ENDDO
!                                                     IF (POS2>0) THEN
!                                                         AC2=ANTICOM(NORB,K,BDOWN(D),L,NEW2)
!                                                         V1((POS2-1)*NSUP+POS)=V1((POS2-1)*NSUP+POS)+&
!                                                         AC2*AC*V((D-1)*NSUP+U)*UL(I,J,K,L)
!                                                     ENDIF
!                                                 ENDIF
!                                             ENDIF
!                                         ENDDO
!                                     ENDDO
!                                 ELSE IF (I .EQ. J .AND. IAND(BDOWN(D),A).EQ.A) THEN
!                                     V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*ULDIAG(I)
!                                 ENDIF
!                             ENDDO
!                         ENDIF
!                     ENDIF
!                 ENDDO
!             ENDDO
!         ENDDO
!         DO D = 1,NSDOWN
!             DO I = 1,NORB
!                 DO J=1,NORB
!                     NEW=KINOP(BDOWN(D),I,J)
!                     IF (NEW>0) THEN
!                         POS=0
!                         DO TMP = 1,NSDOWN
!                             IF (NEW==BDOWN(TMP)) THEN
!                                 POS=TMP
!                                 EXIT 
!                             ENDIF
!                         ENDDO
!                         IF (POS>0) THEN
!                             AC=ANTICOM(NORB,I,BDOWN(D),J,NEW)
!                             DO U = 1,NSUP
!                                 V1((POS-1)*NSUP+U)=V1((POS-1)*NSUP+U)+AC*V((D-1)*NSUP+U)*T(I,J)
!                             ENDDO
!                         ENDIF
!                     ENDIF
!                 ENDDO
!             ENDDO
!         ENDDO
!     ELSE
!         DO U = 1,NSUP
!             DO I = 1,NORB
!                 A = ISHFT(1,I-1)
!                 DO J=1,NORB
!                     DO D = 1,NSDOWN
!                         IF (I .EQ. J .AND. ABS(J_MATRIX(I,I)) .GE. 1.E-14) THEN
!                             IF (IAND(BUP(U),A) .EQ. A .AND. IAND(BDOWN(D),A) .EQ. 0) &
!                             V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*0.5*J_MATRIX(I,I)
!                             IF (IAND(BUP(U),A) .EQ. 0 .AND. IAND(BDOWN(D),A) .EQ. A) &
!                             V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*0.5*J_MATRIX(I,I)
!                         ENDIF
!                         ! V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*J_MATRIX(I,J)*&
!                         ! (SZ_J(BUP(U),I)-SZ_J(BDOWN(D),I))*(SZ_J(BUP(U),J)-SZ_J(BDOWN(D),J))
!                     ENDDO
!                     NEW=KINOP(BUP(U),I,J)
!                     IF (NEW>0) THEN
!                         POS=0
!                         DO TMP = 1,NSUP 
!                             IF (NEW==BUP(TMP)) THEN
!                                 POS=TMP
!                                 EXIT 
!                             ENDIF
!                         ENDDO
!                         IF (POS>0) THEN
!                             AC=ANTICOM(NORB,I,BUP(U),J,NEW)
!                             DO D = 1,NSDOWN 
!                                 V1((D-1)*NSUP+POS)=V1((D-1)*NSUP+POS)+AC*V((D-1)*NSUP+U)*T(I,J)
!                                 V1((POS-1)*NSUP+D)=V1((POS-1)*NSUP+D)+AC*V((U-1)*NSUP+D)*T(I,J)
!                                 IF (I .NE. J .AND. ABS(J_MATRIX(I,J)) .GE. 1.E-14) THEN
!                                     NEW2=KINOP(BDOWN(D),J,I)
!                                     IF (NEW2>0) THEN
!                                         POS2=0
!                                         DO TMP = 1,NSDOWN 
!                                             IF (NEW2==BDOWN(TMP)) THEN
!                                                 POS2=TMP
!                                                 EXIT 
!                                             ENDIF
!                                         ENDDO
!                                         IF (POS2>0) THEN
!                                             AC2=ANTICOM(NORB,J,BDOWN(D),I,NEW2)
!                                             V1((POS2-1)*NSUP+POS)=V1((POS2-1)*NSUP+POS)-&
!                                             AC*AC2*V((D-1)*NSUP+U)*J_MATRIX(I,J)
!                                         ENDIF
!                                     ENDIF
!                                 ENDIF
!                                 IF (.NOT. IS_ULOC) THEN
!                                     DO K = 1,NORB
!                                         DO L = 1,NORB
!                                             IF (ABS(UL(I,J,K,L)) .GT. 1.E-14) THEN
!                                                 NEW2=KINOP(BDOWN(D),K,L)
!                                                 IF (NEW2>0) THEN
!                                                     POS2=0
!                                                     DO TMP = 1,NSDOWN
!                                                         IF (NEW2==BDOWN(TMP)) THEN
!                                                             POS2=TMP
!                                                             EXIT 
!                                                         ENDIF
!                                                     ENDDO
!                                                     IF (POS2>0) THEN
!                                                         AC2=ANTICOM(NORB,K,BDOWN(D),L,NEW2)
!                                                         V1((POS2-1)*NSUP+POS)=V1((POS2-1)*NSUP+POS)+&
!                                                         AC2*AC*V((D-1)*NSUP+U)*UL(I,J,K,L)
!                                                     ENDIF
!                                                 ENDIF
!                                             ENDIF
!                                         ENDDO
!                                     ENDDO
!                                 ELSE IF (I .EQ. J .AND. IAND(BDOWN(D),A).EQ.A) THEN
!                                     V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*ULDIAG(I)
!                                 ENDIF
!                             ENDDO
!                         ENDIF
!                     ENDIF
!                 ENDDO
!             ENDDO
!         ENDDO
!     ENDIF 
!     RETURN
! END SUBROUTINE



SUBROUTINE TPROD(V,V1)
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: V(NSTATES)
    REAL*8 :: V1(NSTATES)
    INTEGER :: I,J,K,L1,NEW,POS,A,A2
    REAL*8 :: AC
    DO J = 1,NORB
        A=ISHFT(1,J-1)
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
                                V1((L1-1)*NSUP+POS)=V1((L1-1)*NSUP+POS)+AC*V((L1-1)*NSUP+I)*T(J,K)
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
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
                                V1((POS-1)*NSUP+I)=V1((POS-1)*NSUP+I)+AC*V((L1-1)*NSUP+I)*T(J,K)
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    RETURN
END subroutine TPROD


SUBROUTINE UPROD(V,V1)
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: V(NSTATES)
    REAL*8 :: V1(NSTATES)
    INTEGER :: I,J,K,M,N,L1,L2,NEW,NEW2,POS,POS2,A,U,D, NDBL
    REAL*8 :: AC,AC2
    DO J = 1,NORB
        IF (IS_ULOC .EQ. 3) THEN
            DO K=1,NORB
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
                                                    V1((POS2-1)*NSUP+POS) = V1((POS2-1)*NSUP+POS)+&
                                                    AC*AC2*V((L1-1)*NSUP+I)*UL(J,K,M,N)
                                                ENDIF
                                            ENDIF
                                        ENDDO
                                    ENDIF
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        else IF (IS_ULOC .EQ. 2) THEN
            A = ISHFT(1,J-1)
            DO U = 1,NSUP 
                IF (IAND(BUP(U),A) .EQ. A) THEN
                    DO D = 1,NSDOWN
                        IF (IAND(BDOWN(D),A) .EQ. A) &
                        V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*ULDIAG(J)
                    ENDDO 
                ENDIF 
            ENDDO 
        ENDIF
    ENDDO
    IF (IS_ULOC .EQ. 1) THEN
        DO U = 1,NSUP 
            DO D = 1,NSDOWN
                NDBL = 0
                DO J = 1,NORB
                    A=ISHFT(1,J-1)
                    IF (((IAND(BUP(U),A) .EQ. A) .AND. (IAND(BDOWN(D),A) .EQ. A))) NDBL=NDBL+1
                ENDDO
                V1((D-1)*NSUP+U)=V1((D-1)*NSUP+U)+V((D-1)*NSUP+U)*NDBL*U_FLOAT
            ENDDO
        ENDDO
    ENDIF
    RETURN
END subroutine UPROD


SUBROUTINE JPROD(V,V1)
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: V(NSTATES)
    REAL*8 :: V1(NSTATES)
    INTEGER :: I,J,K,L1,L2,NEW,NEW2,POS,POS2,A
    REAL*8 :: AC, AC2
    DO J = 1,NORB
        A=ISHFT(1,J-1)
        DO K=1,NORB
            IF (J .EQ. K .AND. ABS(J_MATRIX(J,J)) .GE. 1.E-14) THEN
                ! --- J*0.5* n_i\uparrow(1 - ni_\downarrow) + J*0.5* n_i\downarrow(1 - ni_\uparrow)
                DO I = 1,NSUP
                    DO L1 = 1,NSDOWN 
                        IF (IAND(BUP(I),A) .EQ. A .AND. IAND(BDOWN(L1),A) .EQ. 0) &
                        V1((L1-1)*NSUP+I)=V1((L1-1)*NSUP+I)+V((L1-1)*NSUP+I)*0.5*J_MATRIX(J,J)
                        IF (IAND(BUP(I),A) .EQ. 0 .AND. IAND(BDOWN(L1),A) .EQ. A) &
                        V1((L1-1)*NSUP+I)=V1((L1-1)*NSUP+I)+V((L1-1)*NSUP+I)*0.5*J_MATRIX(J,J)
                    ENDDO
                ENDDO
            ELSE IF (J .NE. K .AND. ABS(J_MATRIX(J,K)) .GE. 1.E-14) THEN
                ! --- SPIN 0.5*S^+KS^-J--- !
                ! --- = -J*0.5*c^\dagger_j\downarrow c_k\downarrow c^\dagger_k\uparrow  c_j\uparrow
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
                                NEW2=KINOP(BDOWN(L1),K,J)
                                IF (NEW2>0) THEN
                                    POS2=0
                                    DO L2 = 1,NSDOWN 
                                        IF (NEW2==BDOWN(L2)) THEN
                                            POS2=L2
                                            EXIT 
                                        ENDIF
                                    ENDDO
                                    IF (POS2>0) THEN
                                        AC2=ANTICOM(NORB,K,BDOWN(L1),J,NEW2)
                                        V1((POS2-1)*NSUP+POS)=V1((POS2-1)*NSUP+POS)-&
                                        AC*AC2*V((L1-1)*NSUP+I)*J_MATRIX(J,K)
                                    ENDIF
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
            ! IF (ABS(J_MATRIX(J,K)) .GE. 1.E-14) THEN
            !     DO I = 1,NSUP 
            !         DO L1 = 1,NSDOWN
            !             V1((L1-1)*NSUP+I)=V1((L1-1)*NSUP+I)+V((L1-1)*NSUP+I)*J_MATRIX(J,K)*&
            !             (SZ_J(BUP(I),J)-SZ_J(BDOWN(L1),J))*(SZ_J(BUP(I),K)-SZ_J(BDOWN(L1),K))
            !         ENDDO
            !     ENDDO
            ! ENDIF
        ENDDO
    ENDDO
    RETURN
END subroutine JPROD