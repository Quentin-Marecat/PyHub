PROGRAM HUBBARD
    USE BASISMOD
    USE HUBMOD
    USE HDF5
    ! ------------------------------------------------ !
    ! --- MAIN PROGRAM FOR SOLVING HUBBARD CLUSTER --- !
    ! ------------------------------------------------ !
    IMPLICIT NONE

    ! --- read input
    CALL INIT_BASIS(NORB)
    NUP = ELEM2COMP(1)
    NDOWN = ELEM2COMP(2)
    DEALLOCATE(ELEM2COMP)
    CALL READ_BASIS(NORB,NUP,NSUP,BUP)
    CALL READ_BASIS(NORB,NDOWN,NSDOWN,BDOWN)
    NSTATES = NSUP*NSDOWN
    NELEC = NUP + NDOWN
    CALL READ_INPUT()
    IF (DO_SOLVE) CALL SOLVE()
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
    USE HUBMOD
    USE HDF5
    ! ---
    ! --- READ INPUT PARAMETERS FOR HUBBARD CALCULATION
    ! ---
    IMPLICIT NONE
    INTEGER*4 :: ERROR, INTEG, MV
    INTEGER(HID_T)  :: FILE_ID, DSET_ID, GRP_ID, SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(4) :: D4
    INTEGER(HSIZE_T), DIMENSION(2) :: D2
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    CALL h5open_f(ERROR)

    CALL h5fopen_f('solver.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)
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
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, IS_ULOC, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    ALLOCATE(UPARTIAL(NORB,NORB))
    UPARTIAL = 0.
    IF (IS_ULOC .EQ. 1) THEN
        call h5aopen_f(GRP_ID, 'U', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_DOUBLE, U_FLOAT, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
    ELSE IF (IS_ULOC .EQ. 2) THEN
        ALLOCATE(ULDIAG(NORB))
        D1=(/NORB/)
        call h5dopen_f(GRP_ID, 'U', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, ULDIAG, D1, ERROR)
        call h5dclose_f(DSET_ID, ERROR)
    ELSE IF (IS_ULOC .EQ. 3) THEN
        D4=(/NORB,NORB,NORB,NORB/)
        ALLOCATE(UL(NORB,NORB,NORB,NORB))
        call h5dopen_f(GRP_ID, 'U', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, UL, D4, ERROR)
        call h5dclose_f(DSET_ID, ERROR)
        D2=(/NORB,NORB/)
        call h5dopen_f(GRP_ID, 'Upartial', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, UPARTIAL, D2, ERROR)
        call h5dclose_f(DSET_ID, ERROR)
    ELSE 
        CONTINUE
    ENDIF
    ! --- 
    D2=(/NORB,NORB/)
    ALLOCATE(T(NORB,NORB))
    ! --- t_matrix
    call h5dopen_f(GRP_ID, 't_matrix', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, T, D2, ERROR)
    call h5dclose_f(DSET_ID, ERROR)
    ! --- J matrix
    ALLOCATE(J_MATRIX(NORB,NORB))
    call h5dopen_f(GRP_ID, 'J', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, J_MATRIX, D2, ERROR)
    call h5dclose_f(DSET_ID, ERROR)
    ! --- 
    call h5gopen_f(GRP_ID, 'lanczos',SGRP_ID, ERROR )
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



SUBROUTINE READ_PSI_HUBBARD(EXC_STATE_)
    USE HUBMOD
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
    CALL h5fopen_f('solver.h5',H5F_ACC_RDONLY_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'solve',GRP_ID, ERROR )
    D1=(/1/)
    CALL h5aopen_f(GRP_ID, 'nstates',DSET_ID, ERROR)
    CALL h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSTATES_,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    ISLANCZOS = .FALSE.
    D1=(/1/)
    CALL h5aopen_f(GRP_ID, 'is_lanczos',DSET_ID, ERROR)
    CALL h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, LOG,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    IF (LOG==1) ISLANCZOS= .TRUE.

    IF (ALLOCATED(PSI)) DEALLOCATE(PSI)
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
        PSI = PSI/norm2(PSI)
    ENDIF
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    RETURN
END SUBROUTINE READ_PSI_HUBBARD


! SUBROUTINE READ_BASIS(NELEC,NSTATES,BASIS)
!     USE HUBMOD
!     USE HDF5
!     ! ---
!     ! --- READ BASIS FROM HUBBARD CALCULATION
!     ! --- ASSUME THAT THE BASIS MATCH WITH INPUT
!     ! ---
!     IMPLICIT NONE
!     CHARACTER(2) :: BASISINDEX, BASISINDEX2
!     INTEGER*4 :: ERROR
!     INTEGER(HID_T)  :: FILE_ID, DSET_ID, GRP_ID, SGRP_ID
!     INTEGER(HSIZE_T), DIMENSION(1) :: D1
    
!     CALL h5open_f(ERROR)
!     CALL h5fopen_f('basis.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)
!     D1=(/1/)
!     write(BASISINDEX, '(I0.2)') NORB
!     call h5gopen_f(FILE_ID, BASISINDEX,GRP_ID, ERROR )
!     write(BASISINDEX2, '(I0.2)') NELEC
!     call h5gopen_f(GRP_ID, BASISINDEX2,SGRP_ID, ERROR )
!     D1 = (/1/)
!     call h5aopen_f(FILE_ID, 'nstates', DSET_ID, ERROR)
!     call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSTATES, D1, ERROR)
!     call h5aclose_f(DSET_ID, ERROR)
!     D1 = (/NSTATES/)
!     IF (ALLOCATED(BASIS)) DEALLOCATE(BASIS)
!     allocate(BASIS(NSTATES))
!     call h5dopen_f(SGRP_ID, 'basis', DSET_ID, ERROR)
!     call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, BASIS, D1, ERROR)
!     call h5dclose_f(DSET_ID, ERROR)

!     call h5gclose_f(SGRP_ID, ERROR )
!     call h5gclose_f(GRP_ID, ERROR )
!     call h5fclose_f(FILE_ID,ERROR)
!     CALL h5close_f(ERROR)
!     IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
!     RETURN
! END SUBROUTINE READ_BASIS