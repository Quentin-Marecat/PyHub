SUBROUTINE C_DAGGER_C(J,K,SPIN)
    USE BASISMOD
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8 :: PSI_WORK(NSTATES)
    INTEGER, INTENT(IN) :: J,K,SPIN
    INTEGER :: I,L1,NEW,POS
    REAL*8 :: AC
    PSI_WORK = 0.
    IF (SPIN.EQ.1) THEN
        DO I = 1,NSUP
            ! --- KINETIC UP
            NEW=KINOP(BUP(I),K,J)
            IF (NEW>0) THEN
                POS=0
                DO L1 = 1,NSUP 
                    IF (NEW==BUP(L1)) THEN
                        POS=L1
                        EXIT 
                    ENDIF
                ENDDO
                IF (POS>0) THEN
                    AC=ANTICOM(NORB,K,BUP(I),J,NEW)
                    DO L1 = 1,NSDOWN 
                        PSI_WORK((L1-1)*NSUP+POS)=PSI_WORK((L1-1)*NSUP+POS)+AC*PSI_IN((L1-1)*NSUP+I)
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    ELSE
        ! --- ELECTRON DOWN 
        DO L1 = 1,NSDOWN
        ! --- KINETIC DOWN
            NEW=KINOP(BDOWN(L1),K,J)
            IF (NEW>0) THEN
                POS=0
                DO I = 1,NSDOWN 
                    IF (NEW==BDOWN(I)) THEN
                        POS=I
                        EXIT 
                    ENDIF
                ENDDO
                IF (POS>0) THEN
                    AC=ANTICOM(NORB,K,BDOWN(L1),J,NEW)
                    DO I = 1,NSUP 
                        PSI_WORK((POS-1)*NSUP+I)=PSI_WORK((POS-1)*NSUP+I)+AC*PSI_IN((L1-1)*NSUP+I)
                     ENDDO
                ENDIF
            ENDIF
        ENDDO
    ENDIF
    PSI_IN = PSI_WORK
    RETURN
END subroutine C_DAGGER_C


SUBROUTINE C_DAGGER_C_UNRES(J,K,SPINK)
    ! --- c^\dag_j\bar{spink} c_k\spink
    USE BASISMOD
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8 :: PSI_WORK(NSTATES)
    INTEGER, INTENT(IN) :: J,K,SPINK
    INTEGER :: I,L1,NEW,POS,A,A2,I2,B2, NEW2,POS2,L2
    REAL*8 :: AC,AC2
    PSI_WORK = 0.
    A=ISHFT(1,J-1)
    A2=ISHFT(1,K-1)
    IF (SPINK.EQ.1) THEN
        DO I = 1,NSUP
            ! --- REMOVE UP
            NEW=IAND(BUP(I),A2)
            IF (NEW .NE. 0) THEN
                POS=0
                DO L1 = 1,NSUP 
                    IF (NEW==BUP(L1)) THEN
                        POS=L1
                        AC=1
                        DO I2 =NORB,1,-1
                            B2=ISHFT(1,I2-1)
                            IF (I2==K) THEN
                                EXIT
                            ELSE IF (IAND(BUP(I),B2)==B2) THEN
                                AC=-AC
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO
                IF (POS .NE. 0) THEN
                    AC2=1
                    B2=ISHFT(1,I2-1)
                    DO I2 =NORB,1,-1
                        B2=ISHFT(1,I2-1)
                        IF (IAND(BUP(POS),B2)==B2) AC2 = -1*AC2
                    ENDDO
                    DO L1 = 1,NSDOWN 
                        NEW2=IAND(BDOWN(L1),A)
                        IF (NEW2 .EQ. 0) THEN
                            POS2=0
                            DO L2 = 1,NSDOWN 
                                IF (NEW2==BDOWN(L2)) THEN
                                    POS2=L2
                                    DO I2 =NORB,1,-1
                                        B2=ISHFT(1,I2-1)
                                        IF (I2==J) THEN
                                            EXIT
                                        ELSE IF (IAND(BDOWN(L1),B2)==B2) THEN
                                            AC2=-AC2
                                        ENDIF
                                    ENDDO
                                ENDIF
                            ENDDO
                            IF (POS2 .NE. 0) THEN
                                PSI_WORK((POS2-1)*NSUP+POS)=PSI_WORK((POS2-1)*NSUP+POS)+AC*AC2*PSI_IN((L1-1)*NSUP+I)
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    ELSE
        DO I = 1,NSDOWN
            ! --- REMOVE DOWN
            NEW=IAND(BDOWN(I),A2)
            IF (NEW .NE. 0) THEN
                POS=0
                DO L1 = 1,NSDOWN
                    IF (NEW==BDOWN(L1)) THEN
                        POS=L1
                        AC=1
                        DO I2 =NORB,1,-1
                            B2=ISHFT(1,I2-1)
                            IF (I2==K) THEN
                                EXIT
                            ELSE IF (IAND(BDOWN(I),B2)==B2) THEN
                                AC=-AC
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO
                IF (POS .NE. 0) THEN
                    AC2=1
                    B2=ISHFT(1,I2-1)
                    DO I2 =NORB,1,-1
                        B2=ISHFT(1,I2-1)
                        IF (IAND(BDOWN(POS),B2)==B2) AC2 = -1*AC2
                    ENDDO
                    DO L1 = 1,NSUP
                        NEW2=IAND(BUP(L1),A)
                        IF (NEW2 .EQ. 0) THEN
                            POS2=0
                            DO L2 = 1,NSUP
                                IF (NEW2==BUP(L2)) THEN
                                    POS2=L2
                                    DO I2 =NORB,1,-1
                                        B2=ISHFT(1,I2-1)
                                        IF (I2==J) THEN
                                            EXIT
                                        ELSE IF (IAND(BUP(L1),B2)==B2) THEN
                                            AC2=-AC2
                                        ENDIF
                                    ENDDO
                                ENDIF
                            ENDDO
                            IF (POS2 .NE. 0) THEN
                                PSI_WORK((POS-1)*NSUP+POS2)=PSI_WORK((POS-1)*NSUP+POS2)+AC*AC2*PSI_IN((I-1)*NSUP+L1)
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    ENDIF
    PSI_IN = PSI_WORK
    RETURN
END subroutine C_DAGGER_C_UNRES





SUBROUTINE NI(J,SPIN)
    USE BASISMOD
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: J,SPIN
    INTEGER :: A,U,D
    A=ISHFT(1,J-1)
    IF (SPIN .EQ. 1) THEN
        DO U = 1,NSUP 
            IF (IAND(BUP(U),A) .NE. A) THEN
                DO D = 1,NSDOWN
                    PSI_IN((D-1)*NSUP+U)=0.
                ENDDO 
            ENDIF 
        ENDDO 
    ELSE
        DO D = 1,NSDOWN 
            IF (IAND(BDOWN(D),A) .NE. A) THEN
                DO U = 1,NSUP
                    PSI_IN((D-1)*NSUP+U)=0.
                ENDDO 
            ENDIF 
        ENDDO 
    ENDIF
    RETURN
END subroutine NI




SUBROUTINE NIM1(J,SPIN)
    USE BASISMOD
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: J,SPIN
    INTEGER :: A,U,D
    A=ISHFT(1,J-1)
    IF (SPIN.EQ.1) THEN
        DO U = 1,NSUP 
            IF (IAND(BUP(U),A) .EQ. A) THEN
                DO D = 1,NSDOWN
                    PSI_IN((D-1)*NSUP+U)=0.
                ENDDO 
            ENDIF 
        ENDDO 
    ELSE
        DO D = 1,NSDOWN 
            IF (IAND(BDOWN(D),A) .EQ. A) THEN
                DO U = 1,NSUP
                    PSI_IN((D-1)*NSUP+U)=0.
                ENDDO 
            ENDIF 
        ENDDO 
    ENDIF
    RETURN
END subroutine NIM1



SUBROUTINE SZ_OPE(J)
    USE BASISMOD
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8 :: PSI_WORK(NSTATES)
    INTEGER,INTENT(IN) :: J
    INTEGER :: A,U,D
    PSI_WORK = 0.
    A=ISHFT(1,J-1)
    DO U = 1,NSUP 
        IF (IAND(BUP(U),A) .EQ. A) THEN
            DO D = 1,NSDOWN
                PSI_WORK((D-1)*NSUP+U)=PSI_WORK((D-1)*NSUP+U)+PSI_IN((D-1)*NSUP+U)/2.
            ENDDO 
        ENDIF 
    ENDDO 
    DO D = 1,NSDOWN 
        IF (IAND(BDOWN(D),A) .EQ. A) THEN
            DO U = 1,NSUP
                PSI_WORK((D-1)*NSUP+U)=PSI_WORK((D-1)*NSUP+U)-PSI_IN((D-1)*NSUP+U)/2.
            ENDDO 
        ENDIF 
    ENDDO 
    PSI_IN = PSI_WORK
    RETURN
END subroutine SZ_OPE

! SUBROUTINE SPLUS(J)
!     ! --- SPIN S^+K--- !
!     USE BASISMOD
!     USE OPEMOD
!     USE FUNCMOD
!     IMPLICIT NONE
!     REAL*8 :: PSI_WORK(NSTATES)
!     INTEGER,INTENT(IN) :: J
!     INTEGER :: I,L1,L2,NEW,NEW2,POS,POS2,A
!     REAL*8 :: AC, AC2
!     PSI_WORK = 0.
!     A=ISHFT(1,J-1)
!     ! --- c^\dagger_j\uparrow  c_j\downarrow) 
!     DO I = 1,NSUP
!         ! --- KINETIC UP
!         NEW=KINOP(BUP(I),J,K)
!         IF (NEW>0) THEN
!             POS=0
!             DO L2 = 1,NSUP 
!                 IF (NEW==BUP(L2)) THEN
!                     POS=L2
!                     EXIT 
!                 ENDIF
!             ENDDO
!             IF (POS>0) THEN
!                 AC=ANTICOM(NORB,J,BUP(I),K,NEW)
!                 DO L1 = 1,NSDOWN 
!                     NEW2=KINOP(BDOWN(L1),K,J)
!                     IF (NEW2>0) THEN
!                         POS2=0
!                         DO L2 = 1,NSDOWN 
!                             IF (NEW2==BDOWN(L2)) THEN
!                                 POS2=L2
!                                 EXIT 
!                             ENDIF
!                         ENDDO
!                         IF (POS2>0) THEN
!                             AC2=ANTICOM(NORB,K,BDOWN(L1),J,NEW2)
!                             PSI_WORK((POS2-1)*NSUP+POS)=PSI_WORK((POS2-1)*NSUP+POS)+&
!                             AC*AC2*PSI_IN((L1-1)*NSUP+I)*0.5
!                         ENDIF
!                     ENDIF
!                 ENDDO
!             ENDIF
!         ENDIF
!     ENDDO
!     PSI_IN = PSI_WORK
!     RETURN
! END subroutine SPLUS

SUBROUTINE SISJ(K,J)
    ! --- SPIN (S^+KS^-J + S^-KS^+J)/2 + S^Z_KS^Z_J--- !
    USE BASISMOD
    USE OPEMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8 :: PSI_WORK(NSTATES)
    INTEGER,INTENT(IN) :: J,K
    INTEGER :: I,L1,L2,NEW,NEW2,POS,POS2,A,L3
    REAL*8 :: AC, AC2
    PSI_WORK = 0.
    A=ISHFT(1,J-1)
    IF (K==J) THEN
        ! --- (n_i\uparrow(1 - ni_\downarrow) + n_i\downarrow(1 - ni_\uparrow))/2
        DO I = 1,NSUP
            DO L1 = 1,NSDOWN 
                IF (IAND(BUP(I),A) .EQ. A .AND. IAND(BDOWN(L1),A) .EQ. 0) &
                PSI_WORK((L1-1)*NSUP+I)=PSI_WORK((L1-1)*NSUP+I)+PSI_IN((L1-1)*NSUP+I)*0.5
                IF (IAND(BUP(I),A) .EQ. 0 .AND. IAND(BDOWN(L1),A) .EQ. A) &
                PSI_WORK((L1-1)*NSUP+I)=PSI_WORK((L1-1)*NSUP+I)+PSI_IN((L1-1)*NSUP+I)*0.5
            ENDDO
        ENDDO
    ELSE IF (J .NE. K) THEN
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
                    AC=ANTICOM(NORB,J,BUP(I),K,NEW)
                    DO L1 = 1,NSDOWN 
                        NEW2=KINOP(BDOWN(L1),K,J)
                        IF (NEW2>0) THEN
                            POS2=0
                            DO L3 = 1,NSDOWN 
                                IF (NEW2==BDOWN(L3)) THEN
                                    POS2=L3
                                    EXIT 
                                ENDIF
                            ENDDO
                            IF (POS2>0) THEN
                                AC2=ANTICOM(NORB,K,BDOWN(L1),J,NEW2)
                                PSI_WORK((POS2-1)*NSUP+POS)=PSI_WORK((POS2-1)*NSUP+POS)-&
                                AC*AC2*PSI_IN((L1-1)*NSUP+I)*0.5
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        DO I = 1,NSDOWN
            ! --- KINETIC DOWN
            NEW=KINOP(BDOWN(I),J,K)
            IF (NEW>0) THEN
                POS2=0
                DO L3 = 1,NSDOWN
                    IF (NEW==BDOWN(L3)) THEN
                        POS2=L3
                        EXIT 
                    ENDIF
                ENDDO
                IF (POS2>0) THEN
                    AC=ANTICOM(NORB,J,BDOWN(I),K,NEW)
                    DO L1 = 1,NSUP
                        NEW2=KINOP(BUP(L1),K,J)
                        IF (NEW2>0) THEN
                            POS=0
                            DO L2 = 1,NSUP 
                                IF (NEW2==BUP(L2)) THEN
                                    POS=L2
                                    EXIT 
                                ENDIF
                            ENDDO
                            IF (POS>0) THEN
                                AC2=ANTICOM(NORB,K,BUP(L1),J,NEW2)
                                PSI_WORK((POS2-1)*NSUP+POS)=PSI_WORK((POS2-1)*NSUP+POS)-&
                                AC*AC2*PSI_IN((I-1)*NSUP+L1)*0.5
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    ENDIF
    DO I = 1,NSUP 
        DO L1 = 1,NSDOWN
            PSI_WORK((L1-1)*NSUP+I)=PSI_WORK((L1-1)*NSUP+I)+PSI_IN((L1-1)*NSUP+I)*&
            (N_J(BUP(I),J)-N_J(BDOWN(L1),J))*(N_J(BUP(I),K)-N_J(BDOWN(L1),K))*0.25
        ENDDO
    ENDDO
    PSI_IN = PSI_WORK
    RETURN
END subroutine SISJ