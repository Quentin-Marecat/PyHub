PROGRAM MAIN
    use BASISMOD
    ! -------------------------------------------------- !
    ! --- MAIN PROGRAM FOR COMPUTING MANY-BODY BASIS --- !
    ! -------------------------------------------------- !
    IMPLICIT NONE
    INTEGER*4 :: I,NORB,NSTATES
    INTEGER*4,ALLOCATABLE ::  BASIS(:)
    INTEGER*4 :: BINOM

    ! --- read input
    CALL INIT_BASIS(NORB)
    DO I = 1,LENELEM2COMP
        NSTATES = BINOM(NORB,ELEM2COMP(I))
        ALLOCATE(BASIS(NSTATES))
        CALL HILBERT(NORB,ELEM2COMP(I),NSTATES,BASIS)
        CALL WRITE_BASIS(NORB,ELEM2COMP(I),NSTATES,BASIS)
        DEALLOCATE(BASIS)
    ENDDO

    CONTAINS

END PROGRAM



SUBROUTINE WRITE_BASIS(NORB,NELEC,NSTATES,BASIS)
    USE HDF5
    ! ---
    ! --- WRITE BASIS SET 
    ! ---
    IMPLICIT NONE
    INTEGER*4,INTENT(IN) :: NORB,NELEC,NSTATES 
    INTEGER*4,INTENT(IN) :: BASIS(NSTATES)
    INTEGER :: ERROR
    CHARACTER(2) :: BASISINDEX, BASISINDEX2
    INTEGER(HID_T)  :: FILE_ID, SPACE_ID, DSET_ID, GRP_ID,  SGRP_ID
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    CALL h5open_f(ERROR)
    CALL h5fopen_f('basis.h5', H5F_ACC_RDWR_F, FILE_ID, ERROR)
    write(BASISINDEX, '(I0.2)') NORB
    call h5gopen_f(FILE_ID, BASISINDEX,GRP_ID, ERROR )
    write(BASISINDEX2, '(I0.2)') NELEC
    D1 = (/1/)
    call h5gcreate_f(GRP_ID, BASISINDEX2,SGRP_ID, ERROR )
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5acreate_f(SGRP_ID, 'nstates',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NSTATES,D1, ERROR)
    CALL h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    D1 = (/NSTATES/)
    CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
    CALL h5dcreate_f(SGRP_ID, 'basis',H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    CALL h5dwrite_f(DSET_ID, H5T_NATIVE_INTEGER, BASIS,D1, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    CALL h5gclose_f(SGRP_ID, ERROR)
    CALL h5gclose_f(GRP_ID, ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN 
END SUBROUTINE


SUBROUTINE HILBERT(NO,NE,NS,BASIS)
    ! ---------------------------------------------------------- !
    ! --- RETURN HILBERT BASIS FOR NE ELECTRONS OF SAME SPIN --- !
    ! ---------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: NO,NE
    INTEGER*4,INTENT(IN) :: NS
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



! INTEGER*4 FUNCTION FACTO(K)
! INTEGER*4,INTENT(IN)::K
! INTEGER*4::I
! FACTO=1 
! DO I = 2,K
!     FACTO=FACTO*I 
! ENDDO
! END FUNCTION FACTO

! INTEGER*8 FUNCTION BINOM(N,K)
!     INTEGER*4,INTENT(IN) :: N,K
!     BINOM=FACTO(N)/(FACTO(K)*FACTO(N-K))
!     IF (BINOM==0) BINOM=1
! END FUNCTION BINOM

recursive function BINOM(n, k) result(res)
integer*4, intent(in) :: n, k
integer*4 :: res

if (k == 0 .or. k == n) then
res = 1
else
res = BINOM(n-1, k-1) + BINOM(n-1, k)
end if
end function BINOM


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