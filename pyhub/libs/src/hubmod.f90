MODULE HUBMOD
    ! ---------------------------------- !
    ! --- GLOBAL DATA FOR HUB SOLVER --- !
    ! ---------------------------------- !
    IMPLICIT NONE

    ! --- INTEGER*4 --- !
    INTEGER*4 :: MAXLCZ, NDEG, NB_COMP_STATES, EXC_STATE,NW, N2, NBLCZ
    INTEGER*4 :: NORB2COMP, NB_POLES(2,2)

    ! --- FLOAT*8 --- !
    REAL*8, ALLOCATABLE ::  UL(:,:,:,:),T(:,:),J_MATRIX(:,:)
    REAL*8, ALLOCATABLE ::  H(:,:), ULDIAG(:)
    REAL*8 :: ACC_LCZ,TP,Z,E0
    REAL*8, ALLOCATABLE :: W(:),VEC(:,:),VAL(:)

    ! --- LOGICAL --- !
    LOGICAL ::  RENORM,DO_BASIS,DO_SOLVE,DO_RQ, DO_SPGF, ISLANCZOS, IS_ULOC, DO_RQ_TB
    logical :: STORE_H

    ! --- variable protected  --- !
    INTEGER*4 ::  NSUP_,NSDOWN_, NELEC_, NSTATES_, NUP_, NDOWN_
    INTEGER*4, ALLOCATABLE :: BUP_(:), BDOWN_(:)
    REAL*8 :: E0_,SZ_
    REAL*8, ALLOCATABLE :: VAL_(:),VEC_(:,:), LCZVEC_(:,:), PSI_(:,:)

END MODULE HUBMOD
