SUBROUTINE SOLVE()
    USE COMMOD 
    USE HDF5
    IMPLICIT NONE
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID,GRP_ID
    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gcreate_f(FILE_ID, 'solve',GRP_ID, ERROR )
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"

    IF (NSTATES<MAXLCZ) THEN
        CALL EXACT_DIAG()
    ELSE 
        CALL LANCZOS()
    ENDIF
    CALL BOLTZMANN()
    CALL WRITE_solve()
    RETURN 
END SUBROUTINE

SUBROUTINE WRITE_solve()
    USE COMMOD 
    USE HUBMOD 
    USE HDF5
    ! ---
    ! --- WRITE BASIS SET FOR HUBBARD CALCULATION
    ! ---
    IMPLICIT NONE
    INTEGER :: ERROR, SIZE, LOG,NDIM_
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID, DSET_ID, GRP_ID
    INTEGER(HSIZE_T), DIMENSION(2) :: D2
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    REAL*8 :: CONV
    ! --- hdf5 solutions --- !
    SIZE=NSTATES 
    CONV=1.E-14
    IF (ISLANCZOS) SIZE=NBLCZ
    IF (ISLANCZOS) CONV=ACC_LCZ
    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'solve',GRP_ID, ERROR )
    NDIM_ = NSTATES
    IF (ISLANCZOS) NDIM_ = NBLCZ
    IF (STORE_H) THEN
        D2=(/NDIM_,NDIM_/)
        CALL h5screate_simple_f(2,D2, SPACE_ID, ERROR)
        CALL h5dcreate_f(GRP_ID, 'H',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
        CALL h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, H,D2, ERROR)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
    ENDIF
    DEALLOCATE(H)
    D2=(/NDIM_,NDIM_/)
    CALL h5screate_simple_f(2,D2, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'eigenvectors',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VEC,D2, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    DEALLOCATE(VEC)
    D1=(/NDIM_/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'eigenvalues',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL,D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    DEALLOCATE(VAL)
    ! --- number Boltzmann weight
    D1=(/1/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5acreate_f(GRP_ID, 'nw',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NW,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    ! --- temperature Boltzmann weight
    D1=(/NW/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'boltzmann',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, W,D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    ! --- nb written states
!    D1=(/1/)
!    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
!    CALL h5acreate_f(GRP_ID, 'nb_comp_states',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
!    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NB_COMP_STATES,D1, ERROR)
!    CALL h5aclose_f(DSET_ID, ERROR)
!    CALL h5sclose_f(SPACE_ID, ERROR)
    ! --- nb degenerated gs
    D1=(/1/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5acreate_f(GRP_ID, 'deg',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NDEG,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    ! --- is lanczos
    LOG = 0 
    IF (ISLANCZOS) LOG=1
    D1=(/1/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5acreate_f(GRP_ID, 'is_lanczos',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, LOG,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN 
END SUBROUTINE



SUBROUTINE READ_PSI(EXC_STATE_)
    USE COMMOD
    USE HDF5
    ! ---
    ! --- READ PSI FROM HUBBARD CALCULATION
    ! ---
    IMPLICIT NONE
    INTEGER*4 :: ERROR,LOG
    INTEGER, INTENT(IN) :: EXC_STATE_
    REAL*8,ALLOCATABLE :: V(:),V2(:,:)
    INTEGER(HID_T)  :: FILE_ID, DSET_ID, GRP_ID,SGRP_ID,SPACE_ID, SSPACE_ID
    INTEGER(HSIZE_T), DIMENSION(2) :: D2, START_2, COUNT_2
    INTEGER(HSIZE_T), DIMENSION(1) :: D1, START_1, COUNT_1
    IF (ALLOCATED(PSI)) DEALLOCATE(PSI)

    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)

    ! --- nstates
    D1=(/1/)
    call h5gopen_f(FILE_ID, 'basis',GRP_ID, ERROR )
    call h5aopen_f(GRP_ID, 'nstates', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSTATES_, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    call h5gclose_f(GRP_ID, ERROR )

    call h5gopen_f(FILE_ID, 'solve',GRP_ID, ERROR )
    ISLANCZOS = .FALSE.
    D1=(/1/)
    CALL h5aopen_f(GRP_ID, 'is_lanczos',DSET_ID, ERROR)
    CALL h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, LOG,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    IF (LOG==1) ISLANCZOS= .TRUE.

    ALLOCATE(PSI(NSTATES_))
    IF (.NOT. ISLANCZOS) THEN
        COUNT_2=(/NSTATES_,1/)
        D2=(/NSTATES_,NSTATES_/)
        CALL h5screate_simple_f(2,D2,SPACE_ID,ERROR)
        CALL h5screate_simple_f(2,COUNT_2,SSPACE_ID,ERROR)
        call h5dopen_f(GRP_ID, 'eigenvectors', DSET_ID, ERROR)
        START_2=(/0,EXC_STATE_-1/)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START_2, COUNT_2, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE,PSI, COUNT_2, ERROR,SSPACE_ID,SPACE_ID)
        call h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SSPACE_ID,ERROR)
        call h5sclose_f(SPACE_ID, ERROR)

        COUNT_1=(/1/)
        D1=(/NSTATES_/)
        CALL h5screate_simple_f(1,D1,SPACE_ID,ERROR)
        CALL h5screate_simple_f(1,COUNT_1,SSPACE_ID,ERROR)
        call h5dopen_f(GRP_ID, 'eigenvalues', DSET_ID, ERROR)
        START_1=(/EXC_STATE_-1/)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START_1, COUNT_1, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE,E0, COUNT_1, ERROR,SSPACE_ID,SPACE_ID)
        call h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SSPACE_ID,ERROR)
        call h5sclose_f(SPACE_ID, ERROR)
    ELSE
        D1=(/1/)
        call h5gopen_f(GRP_ID, 'lanczos',SGRP_ID, ERROR )
        call h5aopen_f(SGRP_ID, 'nb_lcz', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NBLCZ, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        call h5gclose_f(SGRP_ID, ERROR )
        ALLOCATE(V(NBLCZ))
        COUNT_2=(/NBLCZ,1/)
        D2=(/NBLCZ,NBLCZ/)
        CALL h5screate_simple_f(2,D2,SPACE_ID,ERROR)
        CALL h5screate_simple_f(2,COUNT_2,SSPACE_ID,ERROR)
        call h5dopen_f(GRP_ID, 'eigenvectors', DSET_ID, ERROR)
        START_2=(/0,EXC_STATE_-1/)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START_2, COUNT_2, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE,V, COUNT_2, ERROR,SSPACE_ID,SPACE_ID)
        call h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SSPACE_ID,ERROR)
        call h5sclose_f(SPACE_ID, ERROR)
        COUNT_1=(/1/)
        D1=(/NBLCZ/)
        CALL h5screate_simple_f(1,D1,SPACE_ID,ERROR)
        CALL h5screate_simple_f(1,COUNT_1,SSPACE_ID,ERROR)
        call h5dopen_f(GRP_ID, 'eigenvalues', DSET_ID, ERROR)
        START_1=(/EXC_STATE_-1/)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START_1, COUNT_1, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE,E0, COUNT_1, ERROR,SSPACE_ID,SPACE_ID)
        call h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SSPACE_ID,ERROR)
        call h5sclose_f(SPACE_ID, ERROR)
        ALLOCATE(V2(NSTATES_,NBLCZ))
        call h5gopen_f(GRP_ID, 'lanczos',SGRP_ID, ERROR )
        D2=(/NSTATES_,MAXLCZ/)
        START_2=(/0,0/)
        COUNT_2=(/NSTATES_,NBLCZ/)
        CALL h5screate_simple_f(2,D2,SPACE_ID,ERROR)
        CALL h5screate_simple_f(2,COUNT_2,SSPACE_ID,ERROR)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START_2, COUNT_2, ERROR)
        call h5dopen_f(SGRP_ID, 'lanczos_vectors', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE,V2, D2, ERROR,SSPACE_ID,SPACE_ID)
        call h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SSPACE_ID,ERROR)
        call h5sclose_f(SPACE_ID, ERROR)
        call h5gclose_f(SGRP_ID,ERROR)
        PSI = MATMUL(V2,V)       
        DEALLOCATE(V2,V) 
!        PSI = PSI/norm2(PSI)
    ENDIF
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    RETURN
END SUBROUTINE READ_PSI