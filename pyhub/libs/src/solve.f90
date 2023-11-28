SUBROUTINE SOLVE()
    USE HUBMOD
    USE HDF5
    IMPLICIT NONE
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID,GRP_ID
    CALL h5open_f(ERROR)
    CALL h5fopen_f('solver.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
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
    CALL WRITE_SOLVE_HUBBARD()
    RETURN 
END SUBROUTINE

SUBROUTINE WRITE_SOLVE_HUBBARD()
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
    CALL h5fopen_f('solver.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'solve',GRP_ID, ERROR )
    ! --- nstates
    D1=(/1/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5acreate_f(GRP_ID, 'nstates',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NSTATES,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
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