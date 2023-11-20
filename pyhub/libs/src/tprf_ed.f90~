SUBROUTINE TPRF_ED()
    ! --- COMPUTE SINGLE PARTICLE GREEN FUCNTION
    USE COMMOD 
    USE RQMOD
    USE FUNCMOD
    USE HDF5
    IMPLICIT NONE
    INTEGER :: ERROR,J,I,K,M
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID, DSET_ID, SSPACE_ID, SGRP_ID,GRP_ID
    INTEGER(HSIZE_T), DIMENSION(6) :: D6, START_6, COUNT_6
    INTEGER(HSIZE_T), DIMENSION(2) :: D2

    IF (.NOT. ALLOCATED(SIZ)) ALLOCATE(SIZ(4,NW))
    NELEC_ = NELEC
    NUP_ = NUP 
    NDOWN_ = NDOWN
    SZ_ = SZ
    NSTATES_ = NSTATES
    NSUP_ = NSUP
    NSDOWN_ = NSDOWN
    ALLOCATE(BUP_(NSUP_),BDOWN_(NSDOWN_))
    BUP_=BUP 
    BDOWN_ = BDOWN
    E0_ = E0

    D6 = (/N2,2,4,NORB2COMP,NORB2COMP,NW/)
    ! --- QL MATRIX --- !
    DO J = 1,4
        ! --- Q L MATRIX ASSIGNATION --- !
        IF (J==1) THEN
            IF (allocated(Q)) DEALLOCATE(Q,WK)
            CALL QL_E_UP()
        ELSE IF (J==2) THEN
            IF (NUP_ .NE. NDOWN_) THEN 
                DEALLOCATE(Q,WK)
                CALL QL_E_DOWN()
            ENDIF
        ELSE IF (J==3) THEN
            DEALLOCATE(Q,WK)
            CALL QL_H_UP()
        ELSE IF (J==4) THEN
            IF (NUP_ .NE. NDOWN_) THEN 
                DEALLOCATE(Q,WK)
                CALL QL_H_DOWN()
            ENDIF
        ENDIF
        SIZ(J,EXC_STATE) = NSTATES
        IF (ALLOCATED(POLE)) DEALLOCATE(POLE)
        ALLOCATE(POLE(N2,2))
        POLE=0.
        
        COUNT_6=(/N2,2,1,1,1,1/)
        D2=(/N2,2/)
        DO I = 1,NORB2COMP
            DO K = 1,I
                POLE=0.
                DO M = 1,NSTATES
                    POLE(M,1) = WK(M)
                    POLE(M,2) = Q(M,I)*Q(M,K)*W(EXC_STATE)/Z
                ENDDO
                START_6=(/0,0,J-1,I-1,K-1,EXC_STATE-1/)
                CALL h5open_f(ERROR)
                CALL h5fopen_f('hubbard.h5', H5F_ACC_RDWR_F, FILE_ID, ERROR)
                CALL h5gopen_f(FILE_ID,'output',GRP_ID, ERROR)
                CALL h5gopen_f(GRP_ID,'spgf',SGRP_ID, ERROR)
                DO M = 1,2
                    CALL h5screate_simple_f(6,D6,SPACE_ID,ERROR)
                    CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
                    CALL h5dopen_f(SGRP_ID, 'pole_full', DSET_ID, ERROR)
                    call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START_6,COUNT_6,ERROR)
                    call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, POLE, D2, ERROR,SSPACE_ID,SPACE_ID)
                    call h5dclose_f(DSET_ID,ERROR)
                    call h5sclose_f(SSPACE_ID, ERROR)
                    call h5sclose_f(SPACE_ID, ERROR)
                    IF (I==K) EXIT
                    ! --- transpose --- !
                    START_6(5)=I-1
                    START_6(4)=K-1
                ENDDO
                call h5gclose_f(SGRP_ID, ERROR)
                call h5gclose_f(GRP_ID, ERROR)
                call h5fclose_f(FILE_ID, ERROR)
                call h5close_f(ERROR)
                IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
            ENDDO
        ENDDO
    ENDDO
    NELEC = NELEC_
    NUP = NUP_ 
    NDOWN = NDOWN_
    SZ = SZ_
    NSTATES = NSTATES_
    NSUP = NSUP_
    NSDOWN = NSDOWN_
    IF (ALLOCATED(BUP)) DEALLOCATE(BUP)
    IF (ALLOCATED(BDOWN)) DEALLOCATE(BDOWN)
    ALLOCATE(BUP(NSUP),BDOWN(NSDOWN))
    BUP=BUP_ 
    BDOWN = BDOWN_
    DEALLOCATE(BUP_,BDOWN_)
    E0 = E0_
    RETURN
END subroutine TPRF_ED



SUBROUTINE QL_E_UP()
    USE COMMOD
    USE HUBMOD
    USE RQMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4 :: I
    REAL*8,ALLOCATABLE :: Q_L(:)
    ! --- ELEC ADD UP --- !
    NELEC=NELEC_+1
    SZ=SZ_+0.5
    CALL BASIS()
    CALL EXACT_DIAG()
    ALLOCATE(Q(NSTATES,NORB2COMP),WK(NSTATES),Q_L(NSTATES))
    WK(:)=-E0+VAL(:)
    DO I = 1,NORB2COMP
        CALL GF_E_UP(NSTATES,I,VEC,Q_L)
        Q(:,I)=Q_L
    ENDDO
    DEALLOCATE(Q_L)
    CALL ORDER(NORB2COMP,NSTATES,Q,WK)
    RETURN
END SUBROUTINE

SUBROUTINE GF_E_UP(NS,SITEI,V,Q_L)
    USE COMMOD
    USE RQMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SITEI
    INTEGER*4,INTENT(IN) :: NS 
    REAL*8,INTENT(IN) :: V(NS,NS)
    REAL*8,INTENT(OUT) :: Q_L(NS)
    REAL*8,ALLOCATABLE :: PSII(:)
    REAL*8::AC
    INTEGER*4 :: A,I,J,K
    ALLOCATE(PSII(NS))
    PSII=0.
    DO I = 1,NSUP_
        A=ISHFT(1,SITEI-1)
        IF (IAND(BUP_(I),A)==0) THEN
            DO J=1,NSUP
                IF (BUP(J)==BUP_(I)+A) THEN
                    AC=ANTICOM2(NORB,BUP_(I),SITEI)
                    DO K = 1,NSDOWN
                        PSII((K-1)*NSUP+J)=AC*PSI((K-1)*NSUP_+I)
                    ENDDO
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    DO I =1,NS
        Q_L(I)=DOT_PRODUCT(V(:,I),PSII)
    ENDDO
    DEALLOCATE(PSII)
    RETURN 
END SUBROUTINE



SUBROUTINE QL_E_DOWN()
    USE COMMOD
    USE HUBMOD
    USE RQMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4 :: I
    REAL*8,ALLOCATABLE :: Q_L(:)
    ! --- ELEC ADD DOWN --- !
    NELEC=NELEC_+1
    SZ=SZ_-0.5
    CALL BASIS()
    CALL EXACT_DIAG()
    ALLOCATE(Q(NSTATES,NORB2COMP),WK(NSTATES),Q_L(NSTATES))
    WK(:)=-E0+VAL(:)
    DO I = 1,NORB2COMP
        CALL GF_E_DOWN(NSTATES,I,VEC,Q_L)
        Q(:,I)=Q_L
    ENDDO
    DEALLOCATE(Q_L)
    CALL ORDER(NORB2COMP,NSTATES,Q,WK)
    RETURN
END SUBROUTINE
    

SUBROUTINE GF_E_DOWN(NS,SITEI,V,Q_L)
    USE COMMOD
    USE RQMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SITEI
    INTEGER*4, INTENT(IN) :: NS
    REAL*8,INTENT(IN) :: V(NS,NS)
    REAL*8,INTENT(OUT) :: Q_L(NS)
    REAL*8,ALLOCATABLE :: PSII(:)
    REAL*8::AC
    INTEGER*4 :: A,I,J,K
    ALLOCATE(PSII(NS))
    PSII=0.
    DO I = 1,NSDOWN_
        A=ISHFT(1,SITEI-1)
        IF (IAND(BDOWN_(I),A)==0) THEN
            DO J=1,NSDOWN
                IF (BDOWN_(I)+A==BDOWN(J)) THEN
                    AC=ANTICOM2(NORB,BDOWN_(I),SITEI)
                    DO K = 1,NSUP
                        PSII((J-1)*NSUP+K)=AC*PSI((I-1)*NSUP_+K)
                    ENDDO
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    DO I =1,NS
        Q_L(I) = DOT_PRODUCT(V(:,I),PSII)
    ENDDO
    DEALLOCATE(PSII)
    RETURN 
END SUBROUTINE
    

SUBROUTINE QL_H_UP()
    USE COMMOD
    USE HUBMOD
    USE RQMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4 :: I
    REAL*8,ALLOCATABLE :: Q_L(:)
    ! --- ELEC REM UP --- !
    NELEC=NELEC_-1
    SZ=SZ_-0.5
    CALL BASIS()
    CALL EXACT_DIAG()
    ALLOCATE(Q(NSTATES,NORB2COMP),WK(NSTATES),Q_L(NSTATES))
    WK(:)=E0-VAL(:)
    DO I = 1,NORB2COMP
        CALL GF_H_UP(NSTATES,I,VEC,Q_L)
        Q(:,I)=Q_L
    ENDDO
    DEALLOCATE(Q_L)
    CALL ORDER(NORB2COMP,NSTATES,Q,WK)
    RETURN
END SUBROUTINE

SUBROUTINE GF_H_UP(NS,SITEI,V,Q_L)
    USE COMMOD
    USE RQMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SITEI
    INTEGER*4, INTENT(IN) :: NS
    REAL*8,INTENT(IN) :: V(NS,NS)
    REAL*8,INTENT(OUT) :: Q_L(NS)
    REAL*8,ALLOCATABLE :: PSII(:)
    REAL*8::AC
    INTEGER*4 :: A,I,J,K
    ALLOCATE(PSII(NS))
    PSII=0.
    DO I = 1,NSUP_
        A=ISHFT(1,SITEI-1)
        IF (IAND(BUP_(I),A)==A) THEN
            DO J=1,NSUP
                IF (BUP(J)==BUP_(I)-A) THEN
                    AC=ANTICOM2(NORB,BUP_(I),SITEI)
                    DO K = 1,NSDOWN
                        PSII((K-1)*NSUP+J)=AC*PSI((K-1)*NSUP_+I)
                    ENDDO
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    DO I = 1,NS
        Q_L(I) = DOT_PRODUCT(V(:,I),PSII)
    ENDDO
    DEALLOCATE(PSII)
    RETURN 
END SUBROUTINE
    


SUBROUTINE QL_H_DOWN()
    USE COMMOD
    USE HUBMOD
    USE RQMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4 :: I
    REAL*8,ALLOCATABLE :: Q_L(:)
    IF (ALLOCATED(Q)) DEALLOCATE(Q)
    IF (ALLOCATED(WK)) DEALLOCATE(WK)
    ! --- ELEC REM DOWN --- !
    NELEC=NELEC_-1
    SZ=SZ_+0.5
    CALL BASIS()
    CALL EXACT_DIAG()
    ALLOCATE(Q(NSTATES,NORB2COMP),WK(NSTATES),Q_L(NSTATES))
    WK(:)=E0-VAL(:)
    DO I = 1,NORB2COMP
        CALL GF_H_DOWN(NSTATES,I,VEC,Q_L)
        Q(:,I)=Q_L
    ENDDO
    DEALLOCATE(Q_L)
    CALL ORDER(NORB2COMP,NSTATES,Q,WK)
    RETURN
END SUBROUTINE

SUBROUTINE GF_H_DOWN(NS,SITEI,V,Q_L)
    USE COMMOD
    USE RQMOD
    USE FUNCMOD
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: SITEI
    INTEGER*4, INTENT(IN) :: NS
    REAL*8,INTENT(IN) :: V(NS,NS)
    REAL*8,INTENT(OUT) :: Q_L(NS)
    REAL*8,ALLOCATABLE :: PSII(:)
    REAL*8::AC
    INTEGER*4 :: A,I,J,K
    ALLOCATE(PSII(NS))
    PSII=0.
    DO I = 1,NSDOWN_
        A=ISHFT(1,SITEI-1)
        IF (IAND(BDOWN_(I),A)==A) THEN
            DO J=1,NSDOWN
                IF (BDOWN_(I)-A==BDOWN(J)) THEN
                    AC=ANTICOM2(NORB,BDOWN_(I),SITEI)
                    DO K = 1,NSUP
                        PSII((J-1)*NSUP+K)=AC*PSI((I-1)*NSUP_+K)
                    ENDDO
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    DO I = 1,NS 
        Q_L(I) = DOT_PRODUCT(V(:,I),PSII)
    ENDDO
    DEALLOCATE(PSII)
    RETURN 
END SUBROUTINE
