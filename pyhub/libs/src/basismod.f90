MODULE BASISMOD
    ! ---------------------------------- !
    ! --- GLOBAL DATA FOR HUB SOLVER --- !
    ! ---------------------------------- !
    IMPLICIT NONE

    ! --- INTEGER*4 --- !
    INTEGER*4 :: NELEC, NORB
    INTEGER*4 :: NDOWN, NUP, NSTATES, NSDOWN, NSUP, ORDER
    INTEGER*4, ALLOCATABLE :: BUP(:), BDOWN(:)

    ! --- FLOAT*8 --- !
    REAL*8 :: SZ
    REAL*8, ALLOCATABLE :: PSI(:)

    ! --- LOCGICAL --- !
    LOGICAL :: FOCK

END MODULE BASISMOD