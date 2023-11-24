PROGRAM HUBBARD
    USE BASISMOD
    USE HUBMOD
    USE FUNCMOD
    USE HDF5
    ! ------------------------------------------------ !
    ! --- MAIN PROGRAM FOR SOLVING HUBBARD CLUSTER --- !
    ! ------------------------------------------------ !
    IMPLICIT NONE
    LOGICAL :: BOOL

    ! --- read input
    BOOL = .FALSE.
    CALL READ_BASIS(BOOL)
    CALL READ_INPUT()
    IF (DO_SOLVE) CALL SOLVE_HUBBARD()
    CALL READ_BOLTZMANN()
    IF (DO_RQ) CALL STATIC_RQ()

    IF (DO_SPGF) THEN
        IF (NSTATES>MAXLCZ .AND. TP .GT. 1.E-14) THEN
            WRITE(*,*) 'SET TEMPERATURE TO 0 FOR LANCZOS GF'
            STOP
        ENDIF
        IF (NSTATES<MAXLCZ) THEN
            CALL SPGF_ED()
        ELSE
            CALL SPGF_BANDLANCZOS()
        ENDIF
    ENDIF
END PROGRAM


SUBROUTINE READ_INPUT()
    USE BASISMOD
    USE HUBMOD
    USE HDF5
    ! ---
    ! --- READ INPUT PARAMETERS FOR HUBBARD CALCULATION
    ! ---
    IMPLICIT NONE
    INTEGER*4 :: ERROR, INTEG, MV,I
    INTEGER(HID_T)  :: FILE_ID, DSET_ID, GRP_ID, SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(4) :: D4
    INTEGER(HSIZE_T), DIMENSION(2) :: D2
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    CALL h5open_f(ERROR)

    CALL h5fopen_f('basis.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)
    D1=(/1/)
    call h5aopen_f(FILE_ID, 'index', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, BASISINDEX_INT, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    write(BASISINDEX, '(I0.3)') BASISINDEX_INT
    call h5gopen_f(FILE_ID, BASISINDEX,GRP_ID, ERROR )
    call h5gopen_f(GRP_ID, 'input',SGRP_ID, ERROR )
    D1=(/1/)
    ! --- nb_sites
    call h5aopen_f(SGRP_ID, 'nb_sites', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NORB, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ALLOCATE(UL(NORB,NORB,NORB,NORB),T(NORB,NORB),J_MATRIX(NORB,NORB))
    ! --- nb_elec
    call h5aopen_f(SGRP_ID, 'nb_elec', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NELEC, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- sz
    call h5aopen_f(SGRP_ID, 'sz', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_DOUBLE, SZ, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)

    ! ---
    call h5gclose_f(SGRP_ID,ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)


    CALL h5fopen_f('hubbard.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'input',GRP_ID, ERROR )

    ! --- SOLUTION CALC
    call h5aopen_f(GRP_ID, 'do_solution', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    DO_SOLVE=.FALSE.
    IF (MV/=0) THEN
        do_SOLVE=.TRUE.
    ENDIF
    ! --- STORE HAMILTONIAN
    call h5aopen_f(GRP_ID, 'store_H', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    STORE_H = .FALSE.
    IF (MV/=0) THEN
        STORE_H=.TRUE.
    ENDIF
    ! --- 1/2 BODY CALC
    call h5aopen_f(GRP_ID, 'do_rq', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    DO_RQ=.FALSE.
    IF (MV/=0) THEN
        DO_RQ=.TRUE.
    ENDIF
    ! --- 2 BODY CALC
    call h5aopen_f(GRP_ID, 'do_rq_two_body', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    DO_RQ_TB=.FALSE.
    IF (MV/=0) THEN
        DO_RQ_TB=.TRUE.
    ENDIF
    ! --- SPGF CALC
    call h5aopen_f(GRP_ID, 'do_spgf', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    DO_SPGF=.FALSE.
    IF (MV/=0) THEN
        DO_SPGF=.TRUE.
    ENDIF
    ! --- nb_sites_comp
    call h5aopen_f(GRP_ID, 'nb_sites_comp', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NORB2COMP, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- Temperature
    call h5aopen_f(GRP_ID, 'T', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_DOUBLE, TP, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- nb_comp_states
    call h5aopen_f(GRP_ID, 'nb_comp_states', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NB_COMP_STATES, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- excited_state
    call h5aopen_f(GRP_ID, 'excited_state', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, EXC_STATE, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ! --- 
    ! --- is U local
    call h5aopen_f(GRP_ID, 'is_Uloc', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    IS_ULOC=.FALSE.
    IF (MV/=0) THEN
        IS_ULOC=.TRUE.
        ALLOCATE(ULDIAG(NORB))
        D1=(/NORB/)
        call h5dopen_f(GRP_ID, 'U', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, ULDIAG, D1, ERROR)
        call h5dclose_f(DSET_ID, ERROR)
    ELSE
        D4=(/NORB,NORB,NORB,NORB/)
        call h5dopen_f(GRP_ID, 'U', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, UL, D4, ERROR)
        call h5dclose_f(DSET_ID, ERROR)
    ENDIF
    ! --- 
    D2=(/NORB,NORB/)
    ! --- t_matrix
    call h5dopen_f(GRP_ID, 't_matrix', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, T, D2, ERROR)
    call h5dclose_f(DSET_ID, ERROR)
    ! ! --- J matrix
    ! call h5dopen_f(GRP_ID, 'J', DSET_ID, ERROR)
    ! call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, J_MATRIX, D2, ERROR)
    ! call h5dclose_f(DSET_ID, ERROR)
    ! --- 
    call h5gopen_f(GRP_ID, 'lanczos',SGRP_ID, ERROR )
    D1=(/1/)
     ! --- max_lcz
    call h5aopen_f(SGRP_ID, 'max_lcz', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MAXLCZ, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)   
     ! --- acc_lcz
    call h5aopen_f(SGRP_ID, 'acc_lcz', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_DOUBLE, ACC_LCZ, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)   
    ! --- Gram-Schmidt lanczos renormalization
    call h5aopen_f(SGRP_ID, 'renorm', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, INTEG, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    RENORM=.FALSE.
    IF (INTEG/=0) THEN
        RENORM=.TRUE.
    ENDIF
    ! ---
    call h5gclose_f(SGRP_ID,ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    call h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN
END SUBROUTINE READ_INPUT

