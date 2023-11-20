PROGRAM MAIN
    USE BASISMOD 
    USE FUNCMOD
    use HDF5
    ! -------------------------------------------------- !
    ! --- MAIN PROGRAM FOR COMPUTING MANY-BODY BASIS --- !
    ! -------------------------------------------------- !
    IMPLICIT NONE
    logical :: BOOL

    ! --- read input
    BOOL = .TRUE.
    CALL READ_BASIS(BOOL)
    CALL BASIS()
    CALL WRITE_BASIS()
END PROGRAM
