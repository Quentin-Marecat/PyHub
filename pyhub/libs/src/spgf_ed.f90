SUBROUTINE SPGF_ED()
    USE COMMOD
    USE FUNCMOD
    USE HDF5
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: Q(:,:)
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID,GRP_ID,DSET_ID,SPACE_ID,SSPACE_ID
    INTEGER(HSIZE_T), DIMENSION(3) :: D3,COUNT3,START3
    INTEGER(HSIZE_T), DIMENSION(2) :: D2
    INTEGER(HSIZE_T), DIMENSION(1) :: D1

    ! --- create dataset
!    IF (NSTATES<MAXLCZ) THEN
    NB_POLES(1,1) = INT(BINOM(NORB,NUP-1))*INT(BINOM(NORB,NDOWN))
    NB_POLES(2,1) = INT(BINOM(NORB,NUP+1))*INT(BINOM(NORB,NDOWN))
    NB_POLES(1,2) = INT(BINOM(NORB,NDOWN-1))*INT(BINOM(NORB,NUP))
    NB_POLES(2,2) = INT(BINOM(NORB,NDOWN+1))  *INT(BINOM(NORB,NUP)) 

    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gcreate_f(FILE_ID, 'spgf',GRP_ID, ERROR )


    D3=(/NW, NB_POLES(1,1),NORB2COMP/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'Q_lesser_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D3=(/NW, NB_POLES(2,1),NORB2COMP/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'Q_greater_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D3=(/NW, NB_POLES(1,2),NORB2COMP/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'Q_lesser_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D3=(/NW, NB_POLES(2,2),NORB2COMP/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'Q_greater_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"

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
    ! ---
    ! --- LESSER UP
    NELEC=NELEC_-1
    SZ=SZ_-0.5
    CALL BASIS()
    CALL EXACT_DIAG()

    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
    D1=(/NB_POLES(1,1)/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'positions_lesser_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    IF (NUP_ .EQ. NDOWN_) THEN
        CALL h5dcreate_f(GRP_ID, 'positions_lesser_down',H5T_NATIVE_DOUBLE, SPACE_ID,DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
        CALL h5dclose_f(DSET_ID, ERROR)
    ENDIF
    DEALLOCATE(VAL)
    CALL h5sclose_f(SPACE_ID, ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    ALLOCATE(Q(NSTATES,NORB2COMP))
    DO EXC_STATE = 1,NW
        CALL READ_PSI(EXC_STATE)
        CALL GF_H_UP(Q)

        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )

        D3=(/NW, NB_POLES(1,1),NORB2COMP/)
        D2=(/ NB_POLES(1,1),NORB2COMP /)
        CALL h5screate_simple_f(3,D3,SPACE_ID,ERROR)
        CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
        START3=(/EXC_STATE-1,0,0/)
        COUNT3=(/1,NB_POLES(1,1),NORB2COMP/)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
        CALL h5dopen_f(GRP_ID, 'Q_lesser_up', DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D3, ERROR, SSPACE_ID, SPACE_ID)
        call h5dclose_f(DSET_ID,ERROR)
        IF (NUP_ .EQ. NDOWN_) THEN
            CALL h5dopen_f(GRP_ID, 'Q_lesser_down', DSET_ID, ERROR)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D1, ERROR,SSPACE_ID,SPACE_ID)
            call h5dclose_f(DSET_ID,ERROR)
        endif
        call h5sclose_f(SSPACE_ID, ERROR)
        call h5sclose_f(SPACE_ID, ERROR)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    ENDDO
    DEALLOCATE(Q)
    ! --- LESSER DOWN
    IF (NUP_.NE.NDOWN_) THEN
        IF (allocated(Q)) DEALLOCATE(Q)
        NELEC=NELEC_-1
        SZ=SZ_+0.5
        CALL BASIS()
        CALL EXACT_DIAG()

        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
        D1=(/NB_POLES(1,2)/)
        CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
        CALL h5dcreate_f(GRP_ID, 'positions_lesser_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        DEALLOCATE(VAL)
        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
        ALLOCATE(Q(NSTATES,NORB2COMP))
        DO EXC_STATE = 1,NW
            CALL READ_PSI(EXC_STATE)
            CALL GF_H_DOWN(Q)
            ! --- LESSER DOWN
            CALL h5open_f(ERROR)
            CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
            call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )

            D3=(/NW, NB_POLES(1,2),NORB2COMP/)
            D2=(/ NB_POLES(1,2),NORB2COMP /)
            CALL h5screate_simple_f(3,D3,SPACE_ID,ERROR)
            CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
            START3=(/EXC_STATE-1,0,0/)
            COUNT3=(/1,NB_POLES(1,2),NORB2COMP/)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            CALL h5dopen_f(GRP_ID, 'Q_lesser_down', DSET_ID, ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D3, ERROR, SSPACE_ID, SPACE_ID)
            call h5dclose_f(DSET_ID,ERROR)
            call h5sclose_f(SSPACE_ID, ERROR)
            call h5sclose_f(SPACE_ID, ERROR)

            call h5gclose_f(GRP_ID,ERROR)
            call h5fclose_f(FILE_ID,ERROR)
            CALL h5close_f(ERROR)
            IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
        ENDDO 
        DEALLOCATE(Q)
    ENDIF
    ! --- GREATER UP
    NELEC=NELEC_+1
    SZ=SZ_+0.5
    CALL BASIS()
    CALL EXACT_DIAG()
    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
    D1=(/NB_POLES(2,1)/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'positions_greater_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    IF (NUP_ .EQ. NDOWN_) THEN
        CALL h5dcreate_f(GRP_ID, 'positions_greater_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
        CALL h5dclose_f(DSET_ID, ERROR)
    ENDIF
    DEALLOCATE(VAL)
    CALL h5sclose_f(SPACE_ID, ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    ALLOCATE(Q(NSTATES,NORB2COMP))
    DO EXC_STATE = 1,NW
        CALL READ_PSI(EXC_STATE)
        CALL GF_E_UP(Q)

        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )

        D3=(/NW, NB_POLES(2,1),NORB2COMP/)
        D2=(/ NB_POLES(2,1),NORB2COMP /)
        CALL h5screate_simple_f(3,D3,SPACE_ID,ERROR)
        CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
        START3=(/EXC_STATE-1,0,0/)
        COUNT3=(/1,NB_POLES(2,1),NORB2COMP/)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
        CALL h5dopen_f(GRP_ID, 'Q_greater_up', DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D3, ERROR, SSPACE_ID, SPACE_ID)
        call h5dclose_f(DSET_ID,ERROR)
        IF (NUP_ .EQ. NDOWN_) THEN
            CALL h5dopen_f(GRP_ID, 'Q_greater_down', DSET_ID, ERROR)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D1, ERROR,SSPACE_ID,SPACE_ID)
            call h5dclose_f(DSET_ID,ERROR)
        endif
        call h5sclose_f(SSPACE_ID, ERROR)
        call h5sclose_f(SPACE_ID, ERROR)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    ENDDO
    DEALLOCATE(Q)
    ! --- GREATER DOWN
    IF (NUP_.NE.NDOWN_) THEN
        NELEC=NELEC_+1
        SZ=SZ_-0.5
        CALL BASIS()
        CALL EXACT_DIAG()
        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
        D1=(/NB_POLES(2,2)/)
        CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
        CALL h5dcreate_f(GRP_ID, 'positions_greater_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        DEALLOCATE(VAL)
        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
        ALLOCATE(Q(NSTATES,NORB2COMP))
        DO EXC_STATE = 1,NW
            CALL READ_PSI(EXC_STATE)
            CALL GF_E_DOWN(Q)

            CALL h5open_f(ERROR)
            CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
            call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )

            D3=(/NW, NB_POLES(2,2),NORB2COMP/)
            D2=(/ NB_POLES(2,2),NORB2COMP /)
            CALL h5screate_simple_f(3,D3,SPACE_ID,ERROR)
            CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
            START3=(/EXC_STATE-1,0,0/)
            COUNT3=(/1,NB_POLES(2,2),NORB2COMP/)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            CALL h5dopen_f(GRP_ID, 'Q_greater_down', DSET_ID, ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D3, ERROR, SSPACE_ID, SPACE_ID)
            call h5dclose_f(DSET_ID,ERROR)
            call h5sclose_f(SSPACE_ID, ERROR)
            call h5sclose_f(SPACE_ID, ERROR)
        ENDDO
        DEALLOCATE(Q)
        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    ENDIF


    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
    D2=(/2,2/)
    CALL h5screate_simple_f(2,D2, SPACE_ID, ERROR)
    call h5acreate_f(GRP_ID, 'nb_poles', H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    call h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NB_POLES, D2, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN
END SUBROUTINE SPGF_ED





SUBROUTINE GF_E_UP(Q)
    USE COMMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8,INTENT(OUT) :: Q(NSTATES,NORB2COMP)
    REAL*8,ALLOCATABLE :: PSII(:)
    REAL*8::AC
    INTEGER*4 :: A,I,J,K,SITEI
    ALLOCATE(PSII(NSTATES))
    DO SITEI = 1,NORB2COMP
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
        Q(:,SITEI) = MATMUL(PSII,VEC)
    ENDDO
    DEALLOCATE(PSII)
    RETURN 
END SUBROUTINE



SUBROUTINE GF_E_DOWN(Q)
    USE COMMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8,INTENT(OUT) :: Q(NSTATES,NORB2COMP)
    REAL*8,ALLOCATABLE :: PSII(:)
    REAL*8::AC
    INTEGER*4 :: A,I,J,K,SITEI
    ALLOCATE(PSII(NSTATES))
    DO SITEI = 1,NORB2COMP
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
        Q(:,SITEI) = MATMUL(PSII,VEC)
    ENDDO
    DEALLOCATE(PSII)
    RETURN 
END SUBROUTINE



SUBROUTINE GF_H_UP(Q)
    USE COMMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8,INTENT(OUT) :: Q(NSTATES,NORB2COMP)
    REAL*8,ALLOCATABLE :: PSII(:)
    REAL*8::AC
    INTEGER*4 :: A,I,J,K,SITEI
    ALLOCATE(PSII(NSTATES))
    DO SITEI = 1,NORB2COMP
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
        Q(:,SITEI) = MATMUL(PSII,VEC)
    ENDDO
    DEALLOCATE(PSII)
    RETURN 
END SUBROUTINE


SUBROUTINE GF_H_DOWN(Q)
    USE COMMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8,INTENT(OUT) :: Q(NSTATES,NORB2COMP)
    REAL*8,ALLOCATABLE :: PSII(:)
    REAL*8::AC
    INTEGER*4 :: A,I,J,K,SITEI
    ALLOCATE(PSII(NSTATES))
    DO SITEI = 1,NORB2COMP
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
        Q(:,SITEI) = MATMUL(PSII,VEC)
    ENDDO
    DEALLOCATE(PSII)
    RETURN 
END SUBROUTINE
