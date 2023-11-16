SUBROUTINE BOLTZMANN()
    USE COMMOD
    USE HUBMOD
    IMPLICIT NONE
    INTEGER :: NVAL
    REAL*8 :: ZERO
    INTEGER :: I
    NVAL = SIZE(VAL)
    ZERO=1.E-8
    DO NW = 2,NVAL
        IF (EXP((VAL(1)-VAL(NW))/TP) .LE. ZERO ) GOTO 10
    ENDDO 
10  NW=NW-1
    IF (NW .LT. NB_COMP_STATES) NW = NB_COMP_STATES
    IF (ALLOCATED(W)) DEALLOCATE(W)
    ALLOCATE(W(NW))
    W = 0.

    IF (TP .GE. ZERO) THEN
        W = 0.
        DO I = 1,NW
            W(I)=EXP((VAL(1)-VAL(I))/TP)
        ENDDO
        Z = SUM(W)
    ELSE 
        W(1) = 1.
        Z = 1.
    ENDIF
    RETURN
END SUBROUTINE


SUBROUTINE READ_BOLTZMANN()
    USE COMMOD
    USE HDF5
    ! ---
    ! --- READ BOLTZMANN FROM HUBBARD CALCULATION
    ! ---
    IMPLICIT NONE
    INTEGER*4 :: ERROR
    INTEGER(HID_T)  :: FILE_ID, DSET_ID, GRP_ID
    INTEGER(HSIZE_T), DIMENSION(1) :: D1
    IF (ALLOCATED(W)) DEALLOCATE(W)

    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)

    call h5gopen_f(FILE_ID, 'solve',GRP_ID, ERROR )
    D1=(/1/)
    call h5aopen_f(GRP_ID, 'nw', DSET_ID, ERROR)
    call h5aread_f(DSET_ID, H5T_NATIVE_INTEGER, NW, D1, ERROR)
    call h5aclose_f(DSET_ID, ERROR)

    ALLOCATE(W(NW))
    D1=(/NW/)
    call h5dopen_f(GRP_ID, 'boltzmann', DSET_ID, ERROR)
    call h5dread_f(DSET_ID, H5T_NATIVE_DOUBLE, W, D1, ERROR)
    call h5dclose_f(DSET_ID, ERROR)
    Z = SUM(W)

    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    
    IF (EXC_STATE .NE. 0) THEN
        W = 0.
        W(EXC_STATE) = 1.
        Z = 1.
    ENDIF
    RETURN
END SUBROUTINE READ_BOLTZMANN
