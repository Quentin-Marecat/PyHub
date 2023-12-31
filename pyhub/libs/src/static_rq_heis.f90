SUBROUTINE STATIC_RQ()
    USE HUBMOD
    USE HDF5
    IMPLICIT NONE
    real*8,allocatable :: SMATRIX(:,:)
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID, DSET_ID, GRP_ID, MEMSPACE_ID
    INTEGER(HSIZE_T), DIMENSION(3) :: D3,START3,COUNT3
    INTEGER(HSIZE_T), DIMENSION(2) :: D2
    CALL h5open_f(ERROR)
    CALL h5fopen_f('solver.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gcreate_f(FILE_ID, 'reduced_quantities',GRP_ID, ERROR )

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
        CALL READ_PSI_HUBBARD(EXC_STATE)

        ALLOCATE(SMATRIX(NORB,NORB))
        CALL J_SPIN(PSI,PSI,SMATRIX)

        CALL h5open_f(ERROR)
        CALL h5fopen_f('solver.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
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


SUBROUTINE J_SPIN(BRA,KET,SMATRIX)
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: BRA(NSTATES),KET(NSTATES)
    REAL*8 :: SMATRIX(NORB,NORB)
    INTEGER :: I,J,K,L2,NEW,NEW2,POS,BDOWNELEM
    REAL*8 :: AC
    DO J = 1,NORB
        DO K = 1,NORB
            SMATRIX(J,K) = 0.
            IF (J.EQ.K) THEN
                ! --- (n_i\uparrow(1 - ni_\downarrow) + n_i\downarrow(1 - ni_\uparrow))/2
                DO I = 1,NSUP
                    BDOWNELEM = 2**NORB-1 - BUP(I)
                    IF (IS_PART(NORB,J,BUP(I)) .AND. .NOT. IS_PART(NORB,J,BDOWNELEM)) THEN
                        SMATRIX(J,J)=SMATRIX(J,J)+0.5*BRA(I)*KET(I)
                    ELSE IF (.NOT. IS_PART(NORB,J,BUP(I)) .AND. IS_PART(NORB,J,BDOWNELEM)) THEN
                        SMATRIX(J,J)=SMATRIX(J,J)+0.5*BRA(I)*KET(I)
                    ENDIF
                ENDDO
            ELSE 
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
                            AC=ANTICOM(NORB,J,BUP(I))
                            AC=AC*ANTICOM(NORB,K,NEW)
!                            AC=ANTICOM(NORB,J,BUP(I),K,NEW)
                            NEW2=KINOP(NORB,BDOWNELEM,K,J)
                            IF (NEW2 .eq. (2**NORB-1-BUP(POS))) THEN
                                AC=AC*ANTICOM(NORB,K,BDOWNELEM)
                                AC=AC*ANTICOM(NORB,J,NEW2)
!                                AC2=ANTICOM(NORB,K,BDOWNELEM,J,NEW2)
                                SMATRIX(J,K)=SMATRIX(J,K)-AC*BRA(POS)*KET(I)
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
            DO I = 1,NSUP 
                BDOWNELEM = 2**NORB-1 - BUP(I)
                SMATRIX(J,K)=SMATRIX(J,K) + BRA(I)*KET(I)*&
                SZ(NORB,J,BUP(I),BDOWNELEM)*SZ(NORB,K,BUP(I),BDOWNELEM)
            ENDDO
        ENDDO
    ENDDO
    RETURN
END SUBROUTINE