SUBROUTINE LANCZOS()
    USE BASISMOD
    USE HUBMOD
    USE HDF5
    ! -------------------------------------- !
    ! --- LANCZOS DIAG FOR HUBBARD MODEL --- !
    ! -------------------------------------- !
    IMPLICIT NONE
    INTEGER*4, PARAMETER :: NADD=5
    INTEGER*4 :: I
    REAL*8 :: EV,NR
    REAL*8, ALLOCATABLE :: ALPHA(:),BETA(:),V0(:), V(:), V2(:)
    LOGICAL :: CV
    ! --- hdf5 --- !
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID, MEMSPACE, DSET_ID, GRP_ID, SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(1:2) :: D2, START, COUNT
    INTEGER(HSIZE_T), DIMENSION(1) :: D1

    ISLANCZOS=.TRUE.
    ALLOCATE(ALPHA(MAXLCZ),BETA(MAXLCZ))
    ALLOCATE(V0(NSTATES),V(NSTATES),V2(NSTATES))
    DO I = 1,NSTATES
        V(I)=RAND()-0.5
    !    V(I)=1.
    ENDDO
    NR = NORM2(V)
    V=(1/NR)*V
    BETA(1)=NR
    V0=0.
    CV=.FALSE.
    NBLCZ=0
    EV=0.
    NDEG=0
    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'solve',GRP_ID, ERROR )
    call h5gcreate_f(GRP_ID, 'lanczos',SGRP_ID, ERROR )
    D2=(/NSTATES,MAXLCZ/)
    COUNT=(/INT8(NSTATES),INT8(1)/)
    D1=(/INT8(NSTATES)/)
    CALL h5screate_simple_f(2,D2,SPACE_ID,ERROR)
    CALL h5dcreate_f(SGRP_ID,'lanczos_vectors',H5T_NATIVE_DOUBLE,SPACE_ID,DSET_ID,ERROR)
    call h5screate_simple_f(1, D1, MEMSPACE,  ERROR)
    DO WHILE(.NOT. CV .AND. NBLCZ < MAXLCZ-10)
        DO I=1,NADD
            ! --- hdf5 write vector
            START=(/INT8(0),INT8(NBLCZ)/)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START, COUNT, ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, V, D1, ERROR, MEMSPACE, SPACE_ID)
            NBLCZ=NBLCZ+1
            ! --- lcz step
            CALL HPROD(V,V2)
            ALPHA(NBLCZ)=DOT_PRODUCT(V,V2)
            V2 = V2 - ALPHA(NBLCZ)*V - BETA(NBLCZ)*V0
            NR = NORM2(V2)
            BETA(NBLCZ+1) = NR
            V0=V
            V=(1/NR)*V2
        ENDDO
        DEALLOCATE(V2)
        ALLOCATE(H(NBLCZ,NBLCZ),VAL(NBLCZ),VEC(NBLCZ,NBLCZ))
        H=0.
        DO I = 1,NBLCZ-1
            H(I,I)=ALPHA(I)
            H(I,I+1)=BETA(I+1)
            H(I+1,I)=BETA(I+1)
        ENDDO
!        write(*,'(6F10.6)') (alpha(i),i=nblcz-5+1,nblcz) 
!        write(*,'(6F10.6)') (beta(i),i=nblcz-5+2,nblcz+1) 
        H(NBLCZ,NBLCZ)=ALPHA(NBLCZ)
        CALL DIAGMAT(NBLCZ,H,VAL,VEC)
!        write(*,'(6F10.6)') (VAL(i),i=1,5) 
        IF(ABS(VAL(NDEG+1)-EV)<ACC_LCZ) THEN
            IF (ABS(VAL(NDEG+1)-VAL(1))>ACC_LCZ .AND. (NDEG+1)>1) THEN 
                CV=.TRUE.
            ELSE
                DEALLOCATE(VAL,VEC,H)
                NDEG = NDEG+1
            ENDIF
        ELSE
            EV=VAL(NDEG+1)
            DEALLOCATE(VAL,VEC,H)
        ENDIF
        ALLOCATE(V2(NSTATES))
    ENDDO
    IF (NBLCZ>MAXLCZ-10) THEN
        WRITE(*,*) 'FORTRAN LANCZOS ERROR: MAX LCZ STEP REACH'
    ENDIF
! --- hdf5 lanczos --- !
    call h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(MEMSPACE,ERROR)
    CALL h5sclose_f(SPACE_ID,ERROR)
    D1=(/NBLCZ/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    ! --- alpha
    CALL h5dcreate_f(SGRP_ID,'alpha',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, ALPHA(:NBLCZ),D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    DEALLOCATE(ALPHA)
    ! --- beta
    CALL h5dcreate_f(SGRP_ID,'beta',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, (/BETA(:NBLCZ)/),D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    DEALLOCATE(BETA)
    ! --- 
    CALL h5sclose_f(SPACE_ID, ERROR)
    ! --- conv
    D1=(/1/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5acreate_f(SGRP_ID,'lcz_conv',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL(NDEG+1)-EV,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    ! --- nb lcz steps
    CALL h5acreate_f(SGRP_ID,'nb_lcz',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NBLCZ,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    ! --- max lcz steps
    CALL h5acreate_f(SGRP_ID,'max_lcz',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, MAXLCZ,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    CALL h5gclose_f(SGRP_ID, ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    CALL h5fclose_f(FILE_ID,ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in LANCZOS hdf5 files"
    RETURN 
END SUBROUTINE


