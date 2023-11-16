SUBROUTINE DIAGMAT(natom,H0,Val,Vect)
    ! diagonalization matrix H0
    !  Val = transpose(Vect)@H0@Vect
    IMPLICIT NONE
    !
    ! I/O 
    INTEGER*4, intent(in) :: natom
    REAL (KIND=8), intent(in) :: H0(natom,natom)
    REAL (KIND=8), intent(out) :: Val(natom), Vect(natom,natom)
    !
    ! Work arrays
    INTEGER :: info, il, iu, m
    REAL*8,ALLOCATABLE :: XH0(:,:), XORK(:)
    INTEGER, ALLOCATABLE :: IWORK(:), IFAIL(:)
    REAL (KIND=8) :: abstol, vl, vu
    !
    ! ==== code section ====
    ALLOCATE(XH0(NATOM,NATOM),XORK(8*NATOM),IWORK(5*NATOM),IFAIL(NATOM))
    xork = 0.0d0
    iwork = 0
    ifail = 0
    info = 0
    abstol = 1.0d-15
    vl = 0.0d0
    vu = 0.0d0
    il = 1
    iu = natom
    m = 0
    !
    xH0 = H0
    Val = 0.d0
    Vect = 0.d0
    !
    CALL DSYEVX('V','A','U',natom, xH0, natom, vl,vu, il, iu,abstol, m, &
         & Val, Vect, natom,xork, 8*natom, iwork, ifail, info)
    !
    DEALLOCATE(XH0,XORK,IWORK,IFAIL)
!    IF (INFO /= 0) THEN
!      PRINT*, 'Error from LAPAC lib '
!      STOP
!   ENDIF
   !
   RETURN
    !
  END SUBROUTINE DIAGMAT