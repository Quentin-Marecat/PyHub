SUBROUTINE STATIC_RQ()
    USE COMMOD
    USE HDF5
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: DENS_MAT(:,:,:), DBL_OCC(:,:,:,:)
    real*8,allocatable :: NI(:), SMATRIX(:,:)
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID, DSET_ID, GRP_ID, MEMSPACE_ID
    INTEGER(HSIZE_T), DIMENSION(5) :: D5,START5,COUNT5
    INTEGER(HSIZE_T), DIMENSION(4) :: D4
    INTEGER(HSIZE_T), DIMENSION(3) :: D3,START3,COUNT3
    INTEGER(HSIZE_T), DIMENSION(2) :: D2,START2,COUNT2
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gcreate_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )
    D3=(/NW,NORB,NORB/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'density_matrix_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'density_matrix_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D5=(/NW,NORB,NORB,NORB,NORB/)
    CALL h5screate_simple_f(5,D5, SPACE_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'dbl_occ',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D2=(/NW,NORB/)
    CALL h5screate_simple_f(2,D2, SPACE_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'ni',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D3=(/NW,NORB,NORB/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'sisj',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"

    DO EXC_STATE = 1,NW
        CALL READ_PSI(EXC_STATE)

        ALLOCATE(DENS_MAT(NORB,NORB,2))
        CALL ONE_RDM(PSI,PSI,DENS_MAT)

        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )

        D3=(/NW,NORB,NORB/)
        D2=(/NORB,NORB/)
        START3=(/EXC_STATE-1,0,0/)
        COUNT3=(/1,NORB,NORB/)
        CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
        call h5screate_simple_f(2, D2, MEMSPACE_ID,  ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3, COUNT3, ERROR)
        CALL h5dopen_f(GRP_ID, 'density_matrix_up',DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, DENS_MAT(:,:,1), D2, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5dopen_f(GRP_ID, 'density_matrix_down',DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, DENS_MAT(:,:,2), D2, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        CALL h5sclose_f(MEMSPACE_ID, ERROR)
        DEALLOCATE(DENS_MAT)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"

        ALLOCATE(DBL_OCC(NORB,NORB,NORB,NORB),NI(NORB))
        CALL TWO_RDM(PSI,PSI,DBL_OCC,NI)

        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )

        D5=(/NW,NORB,NORB,NORB,NORB/)
        D4=(/NORB,NORB,NORB,NORB/)
        START5=(/EXC_STATE-1,0,0,0,0/)
        COUNT5=(/1,NORB,NORB,NORB,NORB/)
        CALL h5screate_simple_f(5,D5, SPACE_ID, ERROR)
        call h5screate_simple_f(4, D4, MEMSPACE_ID,  ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START5, COUNT5, ERROR)
        CALL h5dopen_f(GRP_ID, 'dbl_occ',DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, DBL_OCC, D4, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        CALL h5sclose_f(MEMSPACE_ID, ERROR)
        DEALLOCATE(DBL_OCC)
        D2=(/NW,NORB/)
        D1=(/NORB/)
        START2=(/EXC_STATE-1,0/)
        COUNT2=(/1,NORB/)
        CALL h5screate_simple_f(2,D2, SPACE_ID, ERROR)
        call h5screate_simple_f(1, D1, MEMSPACE_ID,  ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START2, COUNT2, ERROR)
        CALL h5dopen_f(GRP_ID, 'ni', DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, NI, D1, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        CALL h5sclose_f(MEMSPACE_ID, ERROR)
        DEALLOCATE(NI)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)

        ALLOCATE(SMATRIX(NORB,NORB))
        CALL J_SPIN(PSI,PSI,SMATRIX)

        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )

        D3=(/NW,NORB,NORB/)
        D2=(/NORB,NORB/)
        START3=(/EXC_STATE-1,0,0/)
        COUNT3=(/1,NORB,NORB/)
        CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
        call h5screate_simple_f(2, D2, MEMSPACE_ID,  ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3, COUNT3, ERROR)
        CALL h5dopen_f(GRP_ID, 'sisj', DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, SMATRIX, D2, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        CALL h5sclose_f(MEMSPACE_ID, ERROR)
        DEALLOCATE(SMATRIX)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    ENDDO
    RETURN
END SUBROUTINE



SUBROUTINE ONE_RDM(BRA,KET,DENS_MAT)
    USE COMMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: BRA(NSTATES),KET(NSTATES)
    REAL*8 :: DENS_MAT(NORB,NORB,2)
    INTEGER :: I,J,K,L1,NEW,POS
    REAL*8 :: AC
    IF (NUP .NE. NDOWN) THEN
        DO J = 1,NORB
            DO K=1,J
                DENS_MAT(J,K,1) = 0.
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
                                DENS_MAT(J,K,1)=DENS_MAT(J,K,1)+AC*BRA((L1-1)*NSUP+I)*KET((L1-1)*NSUP+POS)
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
                DENS_MAT(K,J,1) = DENS_MAT(J,K,1)
                DENS_MAT(J,K,2) = 0.
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
                                DENS_MAT(J,K,2)=DENS_MAT(J,K,2)+AC*BRA((L1-1)*NSUP+I)*KET((POS-1)*NSUP+I)
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
                DENS_MAT(K,J,2) = DENS_MAT(J,K,2)
            ENDDO
        ENDDO
    ELSE 
        DO J = 1,NORB
            DO K=1,J
                DENS_MAT(J,K,1) = 0.
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
                                DENS_MAT(J,K,1)=DENS_MAT(J,K,1)+AC*BRA((L1-1)*NSUP+I)*KET((L1-1)*NSUP+POS)
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
                DENS_MAT(K,J,1) = DENS_MAT(J,K,1)
                DENS_MAT(J,K,2) = DENS_MAT(J,K,1)
                DENS_MAT(K,J,2) = DENS_MAT(J,K,1)
            ENDDO
        ENDDO
    ENDIF
    RETURN
END subroutine ONE_RDM


SUBROUTINE TWO_RDM(BRA,KET,DBL_OCC,NI)
    USE COMMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: BRA(NSTATES),KET(NSTATES)
    REAL*8 :: DBL_OCC(NORB,NORB,NORB,NORB), NI(NORB)
    INTEGER :: I,J,K,M,N,L1,L2,NEW,NEW2,POS,POS2,A,U,D
    REAL*8 :: AC,AC2
    DO J = 1,NORB
        IF (DO_RQ_TB) THEN
            DO K=1,J
                DO M = 1,NORB
                    DO N=1,M
                        DBL_OCC(J,K,M,N) = 0.
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
                                                DBL_OCC(J,K,M,N) = DBL_OCC(J,K,M,N)+&
                                                AC*AC2*BRA((L1-1)*NSUP+I)*KET((POS2-1)*NSUP+POS)
                                            ENDIF
                                        ENDIF
                                    ENDDO
                                ENDIF
                            ENDIF
                        ENDDO
                        DBL_OCC(J,K,N,M)=DBL_OCC(J,K,M,N)
                        DBL_OCC(K,J,M,N)=DBL_OCC(J,K,M,N)
                        DBL_OCC(K,J,N,M)=DBL_OCC(J,K,M,N)
                        IF (J.EQ.K .AND. J.EQ.M .AND. J.EQ.N) &
                        NI(J) = DBL_OCC(J,J,J,J)
                    ENDDO
                ENDDO
            ENDDO
        else
            NI(J) = 0.
            A = ISHFT(1,J-1)
            DO U = 1,NSUP 
                IF (IAND(BUP(U),A) .EQ. A) THEN
                    DO D = 1,NSDOWN
                        IF (IAND(BDOWN(D),A) .EQ. A) &
                        NI(J)=NI(J)+BRA((D-1)*NSUP+U)*KET((D-1)*NSUP+U)
                    ENDDO 
                ENDIF 
            ENDDO 
        ENDIF
    ENDDO
    RETURN
END subroutine TWO_RDM


SUBROUTINE J_SPIN(BRA,KET,SMATRIX)
    USE COMMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: BRA(NSTATES),KET(NSTATES)
    REAL*8 :: SMATRIX(NORB,NORB)
    INTEGER :: I,J,K,L1,L2,NEW,NEW2,POS,POS2,A
    REAL*8 :: AC, AC2
    DO J = 1,NORB
        A=ISHFT(1,J-1)
        DO K=1,J
            SMATRIX(J,K) = 0.
            IF (J .EQ. K) THEN
                ! --- 0.5* n_i\uparrow(1 - ni_\downarrow) + 0.5* n_i\downarrow(1 - ni_\uparrow)
                DO I = 1,NSUP
                    DO L1 = 1,NSDOWN 
                        IF (IAND(BUP(I),A) .EQ. A .AND. IAND(BDOWN(L1),A) .EQ. 0) &
                        SMATRIX(J,J)=SMATRIX(J,J)+0.5*BRA((L1-1)*NSUP+I)*KET((L1-1)*NSUP+I)
                        IF (IAND(BUP(I),A) .EQ. 0 .AND. IAND(BDOWN(L1),A) .EQ. A) &
                        SMATRIX(J,J)=SMATRIX(J,J)+0.5*BRA((L1-1)*NSUP+I)*KET((L1-1)*NSUP+I)
                    ENDDO
                ENDDO
            ELSE
                ! --- SPIN 0.5*\sum_JK (S^+KS^-J+S^-KS^+J) = \sum_JK S^+KS^-J--- !
                ! --- = \sum_JK  -c^\dagger_j\downarrow c_k\downarrow c^\dagger_k\uparrow  c_j\uparrow
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
                                        SMATRIX(J,K)=SMATRIX(J,K)-&
                                        AC*AC2*BRA((L1-1)*NSUP+I)*KET((POS2-1)*NSUP+POS)
                                    ENDIF
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
            DO I = 1,NSUP 
                DO L1 = 1,NSDOWN
                    SMATRIX(J,K)=SMATRIX(J,K)+BRA((L1-1)*NSUP+I)*KET((L1-1)*NSUP+I)*&
                    (SZ_J(BUP(I),J)-SZ_J(BDOWN(L1),J))*(SZ_J(BUP(I),K)-SZ_J(BDOWN(L1),K))
                ENDDO
            ENDDO
            SMATRIX(K,J)=SMATRIX(J,K)
        ENDDO
    ENDDO
    RETURN
END subroutine J_SPIN





SUBROUTINE STATIC_RQ_()
    USE COMMOD
    USE FUNCMOD
    USE HDF5
    IMPLICIT NONE
    REAL*8 :: DENS_MAT(NORB,NORB,2), DBL_OCC(NORB,NORB,NORB,NORB)
    real*8 :: NI(NORB), SMATRIX(NORB,NORB)
    INTEGER :: I,J,K,L,U,D,TMP,NEW,NEW2,POS,POS2,A
    REAL*8 :: AC, AC2
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID, DSET_ID, GRP_ID, MEMSPACE_ID
    INTEGER(HSIZE_T), DIMENSION(5) :: D5,START5,COUNT5
    INTEGER(HSIZE_T), DIMENSION(4) :: D4
    INTEGER(HSIZE_T), DIMENSION(3) :: D3,START3,COUNT3
    INTEGER(HSIZE_T), DIMENSION(2) :: D2,START2,COUNT2
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gcreate_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )
    D3=(/NORB,NORB,NW/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'density_matrix_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'density_matrix_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D5=(/NORB,NORB,NORB,NORB,NW/)
    CALL h5screate_simple_f(5,D5, SPACE_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'dbl_occ',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D2=(/NORB,NW/)
    CALL h5screate_simple_f(2,D2, SPACE_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'ni',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D3=(/NORB,NORB,NW/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    call h5dcreate_f(GRP_ID, 'sisj',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    DO EXC_STATE = 1,NW
        CALL READ_PSI(EXC_STATE)
        DBL_OCC = 0.
        DENS_MAT = 0.
        NI = 0.
        SMATRIX = 0.
        IF (NUP .NE. NDOWN) THEN
            DO U = 1,NSUP
                DO I = 1,NORB
                    A = ISHFT(1,I-1)
                    DO J=1,I
                        DO D = 1,NSDOWN
                            IF (I .EQ. J) THEN
                                IF (IAND(BUP(U),A) .EQ. A .AND. IAND(BDOWN(D),A) .EQ. 0) &
                                SMATRIX(I,I)=SMATRIX(I,I)+PSI((D-1)*NSUP+U)**2*0.5
                                IF (IAND(BUP(U),A) .EQ. 0 .AND. IAND(BDOWN(D),A) .EQ. A) &
                                SMATRIX(I,I)=SMATRIX(I,I)+PSI((D-1)*NSUP+U)**2*0.5
                            ENDIF
                            SMATRIX(I,J)=SMATRIX(I,J)+PSI((D-1)*NSUP+U)**2*&
                            (SZ_J(BUP(U),I)-SZ_J(BDOWN(D),I))*(SZ_J(BUP(U),J)-SZ_J(BDOWN(D),J))
                        ENDDO
                        SMATRIX(J,I) = SMATRIX(I,J)
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
                                    DENS_MAT(I,J,1)=DENS_MAT(I,J,1)+AC*PSI((D-1)*NSUP+U)*PSI((D-1)*NSUP+POS)
                                    IF (I .NE. J) THEN
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
                                                SMATRIX(I,J)=SMATRIX(I,J)-&
                                                AC*AC2*PSI((D-1)*NSUP+U)*PSI((POS2-1)*NSUP+POS)
                                            ENDIF
                                        ENDIF
                                    ENDIF
                                    DO K = 1,NORB
                                        DO L = 1,K
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
                                                    DBL_OCC(I,J,K,L)=DBL_OCC(I,J,K,L)+&
                                                    AC2*AC*PSI((D-1)*NSUP+U)*PSI((POS2-1)*NSUP+POS)
                                                ENDIF
                                            ENDIF
                                            IF (I.EQ.J .AND. I.EQ.K .AND. I.EQ.L) &
                                            NI(I) = DBL_OCC(I,I,I,I)
                                        ENDDO
                                    ENDDO
                                ENDDO
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
            DO D = 1,NSDOWN
                DO I = 1,NORB
                    DO J=1,I
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
                                    DENS_MAT(I,J,2)=DENS_MAT(I,J,2)+AC*PSI((D-1)*NSUP+U)*PSI((POS-1)*NSUP+U)
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
                    DO J=1,I
                        DO D = 1,NSDOWN
                            IF (I .EQ. J) THEN
                                IF (IAND(BUP(U),A) .EQ. A .AND. IAND(BDOWN(D),A) .EQ. 0) &
                                SMATRIX(I,I)=SMATRIX(I,I)+PSI((D-1)*NSUP+U)**2*0.5
                                IF (IAND(BUP(U),A) .EQ. 0 .AND. IAND(BDOWN(D),A) .EQ. A) &
                                SMATRIX(I,I)=SMATRIX(I,I)+PSI((D-1)*NSUP+U)**2*0.5
                            ENDIF
                            SMATRIX(I,J)=SMATRIX(I,J)+PSI((D-1)*NSUP+U)**2*&
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
                                    DENS_MAT(I,J,1)=DENS_MAT(I,J,1)+AC*PSI((D-1)*NSUP+U)*PSI((D-1)*NSUP+POS)
                                    DENS_MAT(I,J,2) = DENS_MAT(I,J,1)
                                    IF (I .NE. J) THEN
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
                                                SMATRIX(I,J)=SMATRIX(I,J)-&
                                                AC*AC2*PSI((D-1)*NSUP+U)*PSI((POS2-1)*NSUP+POS)
                                            ENDIF
                                        ENDIF
                                    ENDIF
                                    DO K = 1,NORB
                                        DO L = 1,K
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
                                                    DBL_OCC(I,J,K,L)=DBL_OCC(I,J,K,L)+&
                                                    AC2*AC*PSI((D-1)*NSUP+U)*PSI((POS2-1)*NSUP+POS)
                                                ENDIF
                                            ENDIF
                                            IF (I.EQ.J .AND. I.EQ.K .AND. I.EQ.L) &
                                            NI(I) = DBL_OCC(I,I,I,I)
                                        ENDDO
                                    ENDDO
                                ENDDO
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ENDIF 
        DO I = 1,NORB 
            DO J = 1,I 
                DENS_MAT(J,I,1) = DENS_MAT(I,J,1)
                DENS_MAT(J,I,2) = DENS_MAT(I,J,2)
                SMATRIX(J,I) = SMATRIX(I,J)
                DO K = 1,NORB 
                    DO L = 1,K
                        DBL_OCC(I,J,L,K)=DBL_OCC(I,J,K,L)
                        DBL_OCC(J,I,K,L)=DBL_OCC(I,J,K,L)
                        DBL_OCC(J,I,L,K)=DBL_OCC(I,J,K,L)
                    ENDDO
                ENDDO 
            ENDDO 
        ENDDO
        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )

        D3=(/NORB,NORB,NW/)
        D2=(/NORB,NORB/)
        START3=(/0,0,EXC_STATE-1/)
        COUNT3=(/NORB,NORB,1/)
        CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
        call h5screate_simple_f(2, D2, MEMSPACE_ID,  ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3, COUNT3, ERROR)
        CALL h5dopen_f(GRP_ID, 'density_matrix_up',DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, DENS_MAT(:,:,1), D2, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5dopen_f(GRP_ID, 'density_matrix_down',DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, DENS_MAT(:,:,2), D2, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        CALL h5sclose_f(MEMSPACE_ID, ERROR)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"

        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )

        D5=(/NORB,NORB,NORB,NORB,NW/)
        D4=(/NORB,NORB,NORB,NORB/)
        START5=(/0,0,0,0,EXC_STATE-1/)
        COUNT5=(/NORB,NORB,NORB,NORB,1/)
        CALL h5screate_simple_f(5,D5, SPACE_ID, ERROR)
        call h5screate_simple_f(4, D4, MEMSPACE_ID,  ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START5, COUNT5, ERROR)
        CALL h5dopen_f(GRP_ID, 'dbl_occ',DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, DBL_OCC, D4, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        CALL h5sclose_f(MEMSPACE_ID, ERROR)
        D2=(/NORB,NW/)
        D1=(/NORB/)
        START2=(/0,EXC_STATE-1/)
        COUNT2=(/NORB,1/)
        CALL h5screate_simple_f(2,D2, SPACE_ID, ERROR)
        call h5screate_simple_f(1, D1, MEMSPACE_ID,  ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START2, COUNT2, ERROR)
        CALL h5dopen_f(GRP_ID, 'ni', DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, NI, D1, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        CALL h5sclose_f(MEMSPACE_ID, ERROR)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)

        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )

        D3=(/NORB,NORB,NW/)
        D2=(/NORB,NORB/)
        START3=(/0,0,EXC_STATE-1/)
        COUNT3=(/NORB,NORB,1/)
        CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
        call h5screate_simple_f(2, D2, MEMSPACE_ID,  ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3, COUNT3, ERROR)
        CALL h5dopen_f(GRP_ID, 'sisj', DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, SMATRIX, D2, ERROR, MEMSPACE_ID, SPACE_ID)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
        CALL h5sclose_f(MEMSPACE_ID, ERROR)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    ENDDO
    RETURN
END SUBROUTINE


