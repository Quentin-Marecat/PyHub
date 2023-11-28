MODULE BASISMOD 

    INTEGER*4 :: LENELEM2COMP
    INTEGER*4,ALLOCATABLE :: ELEM2COMP(:)

    CONTAINS 

    SUBROUTINE INIT_BASIS(INDEX)
        USE HDF5
        ! ---
        IMPLICIT NONE
        INTEGER :: INDEX
        INTEGER*4 :: ERROR
        INTEGER(HID_T)  :: FILE_ID, DSET_ID
        INTEGER(HSIZE_T), DIMENSION(1) :: D1
        
        CALL h5open_f(ERROR)
        CALL h5fopen_f('basis.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)
        D1=(/1/)
        call h5aopen_f(FILE_ID, 'index', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, INDEX, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        D1 = (/1/)
        call h5aopen_f(FILE_ID, 'lenelem2comp', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, LENELEM2COMP, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        D1 = (/LENELEM2COMP/)
        ALLOCATE(ELEM2COMP(LENELEM2COMP))
        call h5aopen_f(FILE_ID, 'elem2comp', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, ELEM2COMP, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
        RETURN
    END subroutine INIT_BASIS

    SUBROUTINE READ_BASIS(NORB,NELEC,NSTATES,BASIS)
        USE HDF5
        ! ---
        ! --- READ BASIS FROM HUBBARD CALCULATION
        ! --- ASSUME THAT THE BASIS MATCH WITH INPUT
        ! ---
        IMPLICIT NONE
        INTEGER*4, INTENT(IN) :: NORB,NELEC
        INTEGER*4,INTENT(OUT) :: NSTATES
        INTEGER*4, ALLOCATABLE,INTENT(OUT) :: BASIS(:)
        CHARACTER(2) :: BASISINDEX, BASISINDEX2
        INTEGER*4 :: ERROR
        INTEGER(HID_T)  :: FILE_ID, DSET_ID, GRP_ID, SGRP_ID
        INTEGER(HSIZE_T), DIMENSION(1) :: D1
        
        CALL h5open_f(ERROR)
        CALL h5fopen_f('basis.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)
        write(BASISINDEX, '(I0.2)') NORB
        call h5gopen_f(FILE_ID, BASISINDEX,GRP_ID, ERROR )
        write(BASISINDEX2, '(I0.2)') NELEC
        call h5gopen_f(GRP_ID, BASISINDEX2,SGRP_ID, ERROR )
        call h5aopen_f(SGRP_ID, 'nstates', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSTATES, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        D1 = (/NSTATES/)
        IF (ALLOCATED(BASIS)) DEALLOCATE(BASIS)
        allocate(BASIS(NSTATES))
        call h5dopen_f(SGRP_ID, 'basis', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, BASIS, D1, ERROR)
        call h5dclose_f(DSET_ID, ERROR)
    
        call h5gclose_f(SGRP_ID, ERROR )
        call h5gclose_f(GRP_ID, ERROR )
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
        RETURN
    END SUBROUTINE READ_BASIS

END MODULE BASISMOD