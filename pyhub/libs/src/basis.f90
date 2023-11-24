SUBROUTINE BASIS()
    USE BASISMOD
    USE FUNCMOD
    IMPLICIT NONE 
    integer :: I

    IF (ALLOCATED(BUP)) DEALLOCATE(BUP)
    IF (ALLOCATED(BDOWN)) DEALLOCATE(BDOWN)

    IF (.NOT. FOCK) THEN
        NUP=INT((((2*SZ)+NELEC))/2)
        NDOWN=NELEC-NUP
        NSUP=INT(BINOM(NORB,NUP))
        NSDOWN=INT(BINOM(NORB,NDOWN))
        NSTATES=NSUP*NSDOWN
        ALLOCATE(BUP(NSUP),BDOWN(NSDOWN))
        CALL HILBERT(NORB,NUP,NSUP,BUP)
        IF (INT(2*SZ)/=0) THEN 
            CALL HILBERT(NORB,NDOWN,NSDOWN,BDOWN)
        ELSE 
            BDOWN=BUP
        ENDIF
    ELSE
        NSUP = 2**NORB
        NSDOWN = 2**NORB
        ALLOCATE(BUP(NSUP),BDOWN(NSDOWN))
        NSTATES=NSUP*NSDOWN
        DO I = 0,2**NORB-1
            BUP(I+1) = I 
            BDOWN(I+1) = I
        ENDDO 
    ENDIF
    RETURN 
END SUBROUTINE


SUBROUTINE WRITE_BASIS()
    USE BASISMOD
    USE HDF5
    ! ---
    ! --- WRITE BASIS SET 
    ! ---
    IMPLICIT NONE
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID, DSET_ID, GRP_ID,  SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    CALL h5open_f(ERROR)
    CALL h5fopen_f('basis.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    D1=(/1/)
    call h5aopen_f(FILE_ID, 'index', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, BASISINDEX_INT, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    write(BASISINDEX, '(I0.3)') BASISINDEX_INT
    call h5gopen_f(FILE_ID, BASISINDEX,GRP_ID, ERROR )
        D1=(/1/)
    call h5aopen_f(FILE_ID, 'index', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, BASISINDEX_INT, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    write(BASISINDEX, '(I0.3)') BASISINDEX_INT
    call h5gopen_f(FILE_ID, BASISINDEX,GRP_ID, ERROR )
    call h5gcreate_f(GRP_ID, 'basis',SGRP_ID, ERROR )
    D1=(/1/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    ! --- nstates
    CALL h5acreate_f(SGRP_ID, 'nstates',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NSTATES,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    ! --- nup
    CALL h5acreate_f(SGRP_ID, 'nsup',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NSUP,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    ! --- ndown
    CALL h5acreate_f(SGRP_ID, 'nsdown',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NSDOWN,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    ! --- basis up
    D1=(/NSUP/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5dcreate_f(SGRP_ID,'basis_up',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5dwrite_f(DSET_ID, H5T_NATIVE_INTEGER, BUP,D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    ! --- basis down
    D1=(/NSDOWN/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5dcreate_f(SGRP_ID, 'basis_down',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5dwrite_f(DSET_ID, H5T_NATIVE_INTEGER, BDOWN,D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    CALL h5gclose_f(SGRP_ID, ERROR)
    CALL h5gclose_f(GRP_ID, ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN 
END SUBROUTINE


SUBROUTINE READ_BASIS(INIT)
    USE BASISMOD
    USE HDF5
    ! ---
    ! --- READ BASIS FROM HUBBARD CALCULATION
    ! --- ASSUME THAT THE BASIS MATCH WITH INPUT
    ! ---
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: INIT
    INTEGER*4 :: ERROR, MV
    INTEGER(HID_T)  :: FILE_ID, DSET_ID, GRP_ID, SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    
    CALL h5open_f(ERROR)
    CALL h5fopen_f('basis.h5', H5F_ACC_RDONLY_F, FILE_ID, ERROR)
    D1=(/1/)
    call h5aopen_f(FILE_ID, 'index', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, BASISINDEX_INT, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    write(BASISINDEX, '(I0.3)') BASISINDEX_INT
    call h5gopen_f(FILE_ID, BASISINDEX,GRP_ID, ERROR )
    IF (INIT) THEN
        call h5gopen_f(GRP_ID, 'input',SGRP_ID, ERROR )
        ! --- nb_sites
        D1=(/1/)
        call h5aopen_f(SGRP_ID, 'nb_sites', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NORB, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        ! --- nb_elec
        D1=(/1/)
        call h5aopen_f(SGRP_ID, 'nb_elec', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NELEC, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        ! --- sz
        D1=(/1/)
        call h5aopen_f(SGRP_ID, 'sz', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_DOUBLE, SZ, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        ! --- fock
        D1=(/1/)
        call h5aopen_f(SGRP_ID, 'fock', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, MV, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        FOCK=.FALSE.
        IF (MV/=0) FOCK=.TRUE.
        IF (.NOT. FOCK) THEN
            NUP=INT((((2*SZ)+NELEC))/2)
            NDOWN=NELEC-NUP
        ENDIF
        call h5gclose_f(SGRP_ID, ERROR )
    ELSE
        call h5gopen_f(GRP_ID, 'input',SGRP_ID, ERROR )
        ! --- nb_sites
        D1=(/1/)
        call h5aopen_f(SGRP_ID, 'nb_sites', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NORB, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        call h5gclose_f(SGRP_ID, ERROR )

        IF (ALLOCATED(BUP)) DEALLOCATE(BUP)
        IF (ALLOCATED(BDOWN)) DEALLOCATE(BDOWN)
        IF (ALLOCATED(BFULL)) DEALLOCATE(BFULL)

        call h5gopen_f(GRP_ID, 'basis',SGRP_ID, ERROR )
    !    ISLANCZOS = .FALSE.
    !    IF (MAXLCZ<NSTATES) ISLANCZOS = .TRUE.

        ! --- nsup
        D1=(/1/)
        call h5aopen_f(SGRP_ID, 'nsup', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSUP, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        ! --- basis_up
        ALLOCATE(BUP(NSUP))
        D1=(/NSUP/)
        call h5dopen_f(SGRP_ID, 'basis_up', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, BUP, D1, ERROR)
        call h5dclose_f(DSET_ID, ERROR)
        ! --- nsdown
        call h5aopen_f(SGRP_ID, 'nsdown', DSET_ID, ERROR)
        call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NSDOWN, D1, ERROR)
        call h5aclose_f(DSET_ID, ERROR)
        ! --- basis_down
        ALLOCATE(BDOWN(NSDOWN))
        D1=(/NSDOWN/)
        call h5dopen_f(SGRP_ID, 'basis_down', DSET_ID, ERROR)
        call h5dread_f(DSET_ID, H5T_NATIVE_INTEGER, BDOWN, D1, ERROR)
        call h5dclose_f(DSET_ID, ERROR)
        call h5gclose_f(sGRP_ID, ERROR )
        NSTATES = NSUP * NSDOWN
    ENDIF
    call h5gclose_f(GRP_ID, ERROR )
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN
END SUBROUTINE READ_BASIS


SUBROUTINE HILBERT(NO,NE,NS,BASIS)
    ! ---------------------------------------------------------- !
    ! --- RETURN HILBERT BASIS FOR NE ELECTRONS OF SAME SPIN --- !
    ! ---------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: NO,NE,NS
    INTEGER*4, INTENT(OUT) :: BASIS(NS)
    INTEGER*4 :: CPT,I,MOVE,POS,POSF
    BASIS=0
    DO I =1,NE
        BASIS(1)=BASIS(1)+ISHFT(1,I-1)
    ENDDO
    CPT=1
    POS=NE 
    DO WHILE(CPT<NS)
        DO WHILE(POS<NO)
            BASIS(CPT+1)=MOVE(BASIS(CPT),POS)
            CPT=CPT+1
            POS=POS+1
        ENDDO
        IF (CPT<NS) THEN
            CALL BMOVE(NO,BASIS(CPT),POSF,POS)
            BASIS(CPT+1)=POSF
            CPT=CPT+1
        ENDIF
    ENDDO
    RETURN
! --- RETURN NUMBER OF ELEMENTS IN THE BASIS AND A LIST OF INTEGER REPRESENTING THE BASIS --- !
END SUBROUTINE


INTEGER*4 FUNCTION MOVE(STATEF,POSF)
    INTEGER*4, INTENT(IN) :: STATEF,POSF
    INTEGER*4 :: A,B
    A=ISHFT(1,POSF-1)
    B=IEOR(STATEF,A)
    A=ISHFT(1,POSF)
    MOVE=IEOR(B,A)
    RETURN
END FUNCTION 



SUBROUTINE BMOVE(NO,STATEF,STATEM,POSM)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: NO,STATEF
    INTEGER*4, INTENT(OUT) :: STATEM,POSM
    INTEGER*4 :: A,STATEF_,I,NSTATEM,POSLE
    NSTATEM=0
    STATEF_=STATEF 
    A=ISHFT(1,NO-NSTATEM-1)
    DO WHILE(IAND(STATEF,A) == A)
        NSTATEM=NSTATEM+1
        A=ISHFT(1,NO-NSTATEM-1)
    ENDDO
    DO I=1,NSTATEM
        STATEF_=STATEF_-ISHFT(1,NO-I)
    ENDDO
    POSLE=NO
    A=ISHFT(1,POSLE-1)
    DO WHILE(IAND(STATEF_,A)==0)
        POSLE=POSLE-1
        A=ISHFT(1,POSLE-1)
    ENDDO
    STATEF_=STATEF_-ISHFT(1,POSLE-1)
    STATEM=STATEF_ 
    DO I =1,NSTATEM+1
        STATEM=STATEM+ISHFT(1,POSLE-1+I)
    ENDDO
    POSM=POSLE+1+NSTATEM 
    RETURN
END SUBROUTINE
