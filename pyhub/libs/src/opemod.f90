MODULE OPEMOD
    ! --------------------------------- !
    ! --- GLOBAL DATA FOR OPERATORS --- !
    ! --------------------------------- !
    IMPLICIT NONE

    ! --- INTEGER*4 --- !
    INTEGER*4, ALLOCATABLE :: SITEI(:,:), SITEJ(:,:), SPINI(:,:), SPINJ(:,:), OPE(:,:)
    INTEGER*4, ALLOCATABLE :: bup(:),bdown(:)
    INTEGER*4 :: NB_OPE, MAX_LEN_OPE, PSI_KIND,  OPEINDEX_INT
    INTEGER*4 :: NUP,NSUP,NDOWN,NSDOWN,NELEC,NSTATES,NORB

    ! --- FLOAT*8 --- !
    REAL*8, ALLOCATABLE :: PSI_OUT(:), PSI_IN(:), COEFF(:)
!    REAL*8 :: 

    LOGICAL :: AVG

    CHARACTER(3) :: OPEINDEX

END MODULE OPEMOD