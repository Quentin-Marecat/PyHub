MODULE BASISMOD
    ! ---------------------------------- !
    ! --- GLOBAL DATA FOR HUB SOLVER --- !
    ! ---------------------------------- !
    IMPLICIT NONE

    ! --- INTEGER*4 --- !
    INTEGER*4 :: NELEC, NORB, BASISINDEX_INT
    INTEGER*4 :: NDOWN, NUP, NSTATES, NSDOWN, NSUP, ORDER
    INTEGER*4, ALLOCATABLE :: BUP(:), BDOWN(:), BFULL(:)

    ! --- FLOAT*8 --- !
    REAL*8 :: SZ
    REAL*8, ALLOCATABLE :: PSI(:)

    ! --- LOCGICAL --- !
    LOGICAL :: FOCK, SPIN_RES

    ! --- character --- !   
    CHARACTER(len=3) :: BASISINDEX

END MODULE BASISMOD