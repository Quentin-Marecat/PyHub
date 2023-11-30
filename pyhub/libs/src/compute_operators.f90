PROGRAM OPERATOR
    USE BASISMOD
    USE OPEMOD
    IMPLICIT NONE 
    logical :: BOOL
    INTEGER :: I,J,SPINK,K, NSTATES_
    INTEGER,ALLOCATABLE :: BASIS(:)

    CALL READ_OPERATOR()
    CALL INIT_BASIS(NORB)
    ALLOCATE(BUP(NSUP),BDOWN(NSDOWN))
    K=1
    DO I = 1,NELEMUP
        CALL READ_BASIS(NORB,ELEM2COMP(I),NSTATES_,BASIS)
        DO J = 1,NSTATES_ 
            BUP(K) = BASIS(J)
            K=K+1
        ENDDO
    ENDDO
    K=1
    DO I = 1,NELEMDOWN
        CALL READ_BASIS(NORB,ELEM2COMP(NELEMUP+I),NSTATES_,BASIS)
        DO J = 1,NSTATES_ 
            BDOWN(K) = BASIS(J)
            K=K+1
        ENDDO
    ENDDO
    ALLOCATE(PSI_OUT(NSTATES))
    PSI_OUT = 0.
    IF (.NOT. BELEM) THEN
        ALLOCATE(PSI_IN(NSTATES))
        PSI_OUT = 0.
        DO I = 1,NB_OPE 
            IF (ABS(COEFF(I)) .GT. 1.E-14) THEN
                CALL READ_PSI_IN()
                DO J = 1,MAX_LEN_OPE
                    IF (OPE(I,J) .NE. 0) THEN 
                        IF (OPE(I,J) .EQ. 1) THEN
                            CALL C_DAGGER_C(SITEI(I,J),SITEJ(I,J),SPINI(I,J))
                        ELSE IF (OPE(I,J) .EQ. 2) THEN
                            CALL NI(SITEI(I,J),SPINI(I,J))
                        ELSE IF (OPE(I,J) .EQ. 3) THEN
                            CONTINUE
                        ELSE IF (OPE(I,J) .EQ. 4) THEN
                            CALL NIM1(SITEI(I,J),SPINI(I,J))
                        ELSE IF (OPE(I,J) .EQ. 5) THEN
                            CALL SZ_OPE(SITEI(I,J))
                        ELSE IF (OPE(I,J) .EQ. 6) THEN
                            SPINK=2
                            CALL C_DAGGER_C_UNRES(SITEI(I,J),SITEI(I,J),SPINK)
                        ELSE IF (OPE(I,J) .EQ. 7) THEN
                            SPINK=1
                            CALL C_DAGGER_C_UNRES(SITEI(I,J),SITEI(I,J),SPINK)
                        ELSE IF (OPE(I,J) .EQ. 8) THEN
                            CALL SISJ(SITEI(I,J),SITEJ(I,J))                     
                        ELSE IF (OPE(I,J) .EQ. 9) THEN
                            CALL C_DAGGER_C_UNRES(SITEI(I,J),SITEJ(I,J),SPINJ(I,J))
                        ELSE IF (OPE(I,J) .EQ. 10) THEN
                            CALL C_DAGGER(SITEI(I,J),SPINI(I,J))
                        ELSE IF (OPE(I,J) .EQ. 11) THEN
                            CALL C_(SITEI(I,J),SPINI(I,J))
                        else
                            CONTINUE
                        ENDIF
                    ENDIF
                ENDDO
                PSI_OUT = PSI_OUT + COEFF(I)*PSI_IN
            ENDIF
        ENDDO
    ELSE
        CALL OPERATOR_BASIS()
    ENDIF
    CALL WRITE_SOLUTION()
END PROGRAM





SUBROUTINE WRITE_SOLUTION()
    USE BASISMOD
    USE OPEMOD
    USE HDF5
    IMPLICIT NONE
    INTEGER*4 :: ERROR, OPEINDEX_INT
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID,DSET_ID, GRP_ID, SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    CHARACTER(3) :: OPEINDEX

    CALL h5open_f(ERROR)
    CALL h5fopen_f('operators.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5aopen_f(FILE_ID, 'index', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, OPEINDEX_INT, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    write(OPEINDEX, '(I0.3)') OPEINDEX_INT
    call h5gopen_f(FILE_ID, OPEINDEX,GRP_ID, ERROR )

    call h5gopen_f(GRP_ID, 'psi',SGRP_ID, ERROR )
    IF (.NOT. AVG) THEN
        D1=(/NSTATES/)
        CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
        CALL h5dcreate_f(SGRP_ID,'output',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
        CALL h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, PSI_OUT,D1, ERROR)
        CALL h5dclose_f(DSET_ID, ERROR)
        CALL h5sclose_f(SPACE_ID, ERROR)
    ELSE
        IF (.NOT. BELEM) THEN
            D1=(/NSTATES/)
            call h5dopen_f(SGRP_ID, 'input', DSET_ID, ERROR)
            call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, PSI_IN, D1, ERROR)
            call h5dclose_f(DSET_ID, ERROR)
            AVERAGE = dot_product(PSI_IN,PSI_OUT)
            DEALLOCATE(PSI_IN)
        ELSE
            AVERAGE = PSI_OUT((POSD-1)*NSUP+POSU)
        ENDIF
        D1=(/1/)
        call h5aopen_f(SGRP_ID, 'avg_value', DSET_ID, ERROR)
        call h5awrite_f(DSET_ID, H5T_NATIVE_DOUBLE, AVERAGE, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
    ENDIF
    IF (ALLOCATED(PSI_OUT)) DEALLOCATE(PSI_OUT)
    call h5gclose_f(GRP_ID,ERROR)
    call h5gclose_f(SGRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN
END SUBROUTINE


SUBROUTINE READ_PSI_IN()
    USE BASISMOD
    USE OPEMOD
    USE HDF5
    INTEGER*4 :: ERROR, OPEINDEX_INT
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID,DSET_ID, GRP_ID, SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    CHARACTER(3) :: OPEINDEX

    CALL h5open_f(ERROR)
    CALL h5fopen_f('operators.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)

    call h5aopen_f(FILE_ID, 'index', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, OPEINDEX_INT, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    write(OPEINDEX, '(I0.3)') OPEINDEX_INT
    call h5gopen_f(FILE_ID, OPEINDEX,GRP_ID, ERROR )

    call h5gopen_f(GRP_ID, 'psi',SGRP_ID, ERROR )
    D1=(/NSTATES/)
    call h5dopen_f(SGRP_ID, 'input', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, PSI_IN, D1, ERROR)
    call h5dclose_f(DSET_ID, ERROR)

    call h5gclose_f(SGRP_ID,ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
END SUBROUTINE READ_PSI_IN



SUBROUTINE READ_OPERATOR()
    USE OPEMOD
    USE HDF5
    ! ---
    ! ---
    IMPLICIT NONE
    INTEGER*4 :: ERROR, MV, OPEINDEX_INT
    INTEGER(HID_T)  :: FILE_ID, DSET_ID, GRP_ID, SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    INTEGER(HSIZE_T), DIMENSION(2) :: D2
    CHARACTER(3) :: OPEINDEX
    
    CALL h5open_f(ERROR)
    CALL h5fopen_f('operators.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)
    call h5aopen_f(FILE_ID, 'index', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, OPEINDEX_INT, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    write(OPEINDEX, '(I0.3)') OPEINDEX_INT
    call h5gopen_f(FILE_ID, OPEINDEX,GRP_ID, ERROR )

    CALL h5aopen_f(GRP_ID, 'basis_elem', DSET_ID, ERROR)
    CALL h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    BELEM = .FALSE. 
    IF (MV .EQ. 1) BELEM = .TRUE.

    call h5gopen_f(GRP_ID, 'psi',SGRP_ID, ERROR )
    ! --- kind of psi
    D1=(/1/)
    call h5aopen_f(SGRP_ID, 'nstates', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSTATES, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- kind of psi
    D1=(/1/)
    call h5aopen_f(SGRP_ID, 'nsup', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSUP, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- kind of psi
    D1=(/1/)
    call h5aopen_f(SGRP_ID, 'nsdown', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSDOWN, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- kind of psi
    D1=(/1/)
    call h5aopen_f(SGRP_ID, 'nelemup', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NELEMUP, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- kind of psi
    D1=(/1/)
    call h5aopen_f(SGRP_ID, 'nelemdown', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NELEMDOWN, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- kind of psi
    D1=(/1/)
    call h5aopen_f(SGRP_ID, 'avg', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    AVG = .FALSE. 
    IF (MV .EQ. 1) AVG = .TRUE.
    IF (BELEM) THEN
        call h5aopen_f(SGRP_ID, 'posu', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, POSU, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        call h5aopen_f(SGRP_ID, 'posd', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, POSD, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        call h5aopen_f(SGRP_ID, 'up', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, BU_, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        call h5aopen_f(SGRP_ID, 'down', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, BD_, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
    ENDIF

    call h5gclose_f(SGRP_ID,ERROR)


    call h5gopen_f(GRP_ID, 'operators',SGRP_ID, ERROR )
    ! -- nb operations
    D1=(/1/)
    call h5aopen_f(SGRP_ID, 'nb_ope', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NB_OPE, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! -- nb operations
    D1=(/1/)
    call h5aopen_f(SGRP_ID, 'max_nb_ope', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MAX_LEN_OPE, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)

    ALLOCATE(COEFF(NB_OPE),SITEI(NB_OPE,MAX_LEN_OPE), SITEJ(NB_OPE,MAX_LEN_OPE), SPINI(NB_OPE,MAX_LEN_OPE))
    ALLOCATE(SPINJ(NB_OPE,MAX_LEN_OPE), OPE(NB_OPE,MAX_LEN_OPE))

    D1=(/NB_OPE/)
    call h5dopen_f(SGRP_ID, 'coeff', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, COEFF, D1, ERROR)
    call h5dclose_f(DSET_ID, ERROR)

    D2=(/NB_OPE,MAX_LEN_OPE/)
    call h5dopen_f(SGRP_ID, 'site1_index', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, SITEI, D1, ERROR)
    call h5dclose_f(DSET_ID, ERROR)
    call h5dopen_f(SGRP_ID, 'spin1_index', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, SPINI, D1, ERROR)
    call h5dclose_f(DSET_ID, ERROR)
    call h5dopen_f(SGRP_ID, 'site2_index', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, SITEJ, D1, ERROR)
    call h5dclose_f(DSET_ID, ERROR)
    call h5dopen_f(SGRP_ID, 'spin2_index', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, SPINJ, D1, ERROR)
    call h5dclose_f(DSET_ID, ERROR)
    call h5dopen_f(SGRP_ID, 'str_index', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, OPE, D1, ERROR)
    call h5dclose_f(DSET_ID, ERROR)

    call h5gclose_f(SGRP_ID,ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN
END SUBROUTINE READ_OPERATOR