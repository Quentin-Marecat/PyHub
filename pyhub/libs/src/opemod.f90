MODULE OPEMOD
    ! --------------------------------- !
    ! --- GLOBAL DATA FOR OPERATORS --- !
    ! --------------------------------- !
    IMPLICIT NONE

    ! --- INTEGER*4 --- !
    INTEGER*4, ALLOCATABLE :: SITEI(:,:), SITEJ(:,:), SPINI(:,:), SPINJ(:,:), OPE(:,:)
    INTEGER*4, ALLOCATABLE :: bup(:),bdown(:)
    INTEGER*4 :: NB_OPE, MAX_LEN_OPE, PSI_KIND, POSU,POSD,BU_,BD_
    INTEGER*4 :: NUP,NSUP,NDOWN,NSDOWN,NELEC,NSTATES,NORB, NELEMUP,NELEMDOWN

    ! --- FLOAT*8 --- !
    REAL*8, ALLOCATABLE :: PSI_OUT(:), PSI_IN(:), COEFF(:)
    REAL*8 :: AVERAGE

    LOGICAL :: AVG, BELEM

END MODULE OPEMOD