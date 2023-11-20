SUBROUTINE SPGF_BANDLANCZOS()
    USE BASISMOD
    USE HUBMOD
    USE FUNCMOD
    USE HDF5
    IMPLICIT NONE
    INTEGER :: I,J
    REAL*8,ALLOCATABLE :: Q(:,:)
    REAL*8,ALLOCATABLE :: PSI_LCZ(:,:)
    real*8 :: COEFF_INIT(NORB2COMP,NORB2COMP)
    INTEGER :: ERROR
    INTEGER(HID_T)  :: FILE_ID,GRP_ID,DSET_ID,SPACE_ID,SSPACE_ID
    INTEGER(HSIZE_T), DIMENSION(3) :: D3,COUNT3,START3
    INTEGER(HSIZE_T), DIMENSION(2) :: D2
    INTEGER(HSIZE_T), DIMENSION(1) :: D1

    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gcreate_f(FILE_ID, 'spgf',GRP_ID, ERROR )

    D3=(/NW, MAXLCZ,NORB2COMP/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'Q_lesser_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D3=(/NW, MAXLCZ,NORB2COMP/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'Q_greater_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D3=(/NW, MAXLCZ,NORB2COMP/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'Q_lesser_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    D3=(/NW, MAXLCZ,NORB2COMP/)
    CALL h5screate_simple_f(3,D3, SPACE_ID, ERROR)
    CALL h5dcreate_f(GRP_ID, 'Q_greater_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
    CALL h5dclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)

    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"

    NELEC_ = NELEC
    NUP_ = NUP 
    NDOWN_ = NDOWN
    SZ_ = SZ
    NSTATES_ = NSTATES
    NSUP_ = NSUP
    NSDOWN_ = NSDOWN
    ALLOCATE(BUP_(NSUP_),BDOWN_(NSDOWN_))
    BUP_=BUP 
    BDOWN_ = BDOWN
    E0_ = E0
    ! ---
    IF(ALLOCATED(PSI)) DEALLOCATE(PSI)
    ALLOCATE(PSI(NSTATES_))
    DO EXC_STATE = 1,1
        CALL READ_PSI_HUBBARD(EXC_STATE)
        ! --- LESSER UP
        NELEC=NELEC_-1
        SZ=SZ_-0.5
        CALL BASIS()
        ALLOCATE(PSI_LCZ(NSTATES,NORB2COMP))
        CALL LCZ_H_UP(PSI_LCZ)
        CALL GF_BANDLANCZOS(PSI_LCZ,COEFF_INIT)
        DEALLOCATE(PSI_LCZ)
        ALLOCATE(Q(NBLCZ,NORB2COMP))
        Q(:,1)=COEFF_INIT(1,1)*VEC(1,:)
        DO I = 1,NORB2COMP
            Q(:,I)=COEFF_INIT(I,I)*VEC(I,:)
            DO J = 1,I-1
                Q(:,I)=Q(:,I)+COEFF_INIT(I,J)*VEC(J,:)
            ENDDO
        ENDDO
        DEALLOCATE(VEC)
        NB_POLES(1,1) = NBLCZ
        ! --- 
        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
        D1=(/ NBLCZ/)
        CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
        CALL h5dcreate_f(GRP_ID, 'positions_lesser_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
        CALL h5dclose_f(DSET_ID, ERROR)
        IF (NUP_ .EQ. NDOWN_) THEN
            NB_POLES(1,2) = NBLCZ
            CALL h5dcreate_f(GRP_ID, 'positions_lesser_down',H5T_NATIVE_DOUBLE, SPACE_ID,DSET_ID, ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
            CALL h5dclose_f(DSET_ID, ERROR)
        ENDIF
        DEALLOCATE(VAL)

        D3=(/NW,MAXLCZ,NORB2COMP/)
        D2=(/NBLCZ,NORB2COMP /)
        CALL h5screate_simple_f(3,D3,SPACE_ID,ERROR)
        CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
        START3=(/EXC_STATE-1,0,0/)
        COUNT3=(/1,NBLCZ,NORB2COMP/)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
        CALL h5dopen_f(GRP_ID, 'Q_lesser_up', DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D3, ERROR, SSPACE_ID, SPACE_ID)
        call h5dclose_f(DSET_ID,ERROR)
        IF (NUP_ .EQ. NDOWN_) THEN
            CALL h5dopen_f(GRP_ID, 'Q_lesser_down', DSET_ID, ERROR)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D1, ERROR,SSPACE_ID,SPACE_ID)
            call h5dclose_f(DSET_ID,ERROR)
        endif
        call h5sclose_f(SSPACE_ID, ERROR)
        call h5sclose_f(SPACE_ID, ERROR)
        DEALLOCATE(Q)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"

        ! --- LESSER DOWN
        IF (NUP_.NE.NDOWN_) THEN
            NELEC=NELEC_-1
            SZ=SZ_+0.5
            CALL BASIS()
            ALLOCATE(PSI_LCZ(NSTATES,NORB2COMP))
            CALL LCZ_H_DOWN(PSI_LCZ)
            CALL GF_BANDLANCZOS(PSI_LCZ,COEFF_INIT)
            DEALLOCATE(PSI_LCZ)
            ALLOCATE(Q(NBLCZ,NORB2COMP))
            Q(:,1)=COEFF_INIT(1,1)*VEC(1,:)
            DO I = 2,NORB2COMP
                Q(:,I)=COEFF_INIT(I,I)*VEC(I,:)
                DO J = 1,I-1
                    Q(:,I)=Q(:,I)+COEFF_INIT(I,J)*VEC(J,:)
                ENDDO
            ENDDO
            DEALLOCATE(VEC)
            NB_POLES(1,2) = NBLCZ

                   ! --- 
       ! --- 
            CALL h5open_f(ERROR)
            CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
            call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
            D1=(/ NBLCZ/)
            CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
            CALL h5dcreate_f(GRP_ID, 'positions_lesser_down',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
            CALL h5dclose_f(DSET_ID, ERROR)
            DEALLOCATE(VAL)
    
            D3=(/NW,MAXLCZ,NORB2COMP/)
            D2=(/NBLCZ,NORB2COMP /)
            CALL h5screate_simple_f(3,D3,SPACE_ID,ERROR)
            CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
            START3=(/EXC_STATE-1,0,0/)
            COUNT3=(/1,NBLCZ,NORB2COMP/)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            CALL h5dopen_f(GRP_ID, 'Q_lesser_down', DSET_ID, ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D3, ERROR, SSPACE_ID, SPACE_ID)
            call h5dclose_f(DSET_ID,ERROR)
            call h5sclose_f(SSPACE_ID, ERROR)
            call h5sclose_f(SPACE_ID, ERROR)
            DEALLOCATE(Q)
    
            call h5gclose_f(GRP_ID,ERROR)
            call h5fclose_f(FILE_ID,ERROR)
            CALL h5close_f(ERROR)
            IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
        ENDIF
        ! --- GREATER UP
        NELEC=NELEC_+1
        SZ=SZ_+0.5
        CALL BASIS()
        ALLOCATE(PSI_LCZ(NSTATES,NORB2COMP))
        CALL LCZ_E_UP(PSI_LCZ)
        CALL GF_BANDLANCZOS(PSI_LCZ,COEFF_INIT)
        DEALLOCATE(PSI_LCZ)
        ALLOCATE(Q(NBLCZ,NORB2COMP))
        Q(:,1)=COEFF_INIT(1,1)*VEC(1,:)
        DO I = 2,NORB2COMP
            Q(:,I)=COEFF_INIT(I,I)*VEC(I,:)
            DO J = 1,I-1
                Q(:,I)=Q(:,I)+COEFF_INIT(I,J)*VEC(J,:)
            ENDDO
        ENDDO
        DEALLOCATE(VEC)
        NB_POLES(2,1) = NBLCZ

       ! --- 
        CALL h5open_f(ERROR)
        CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
        call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
        D1=(/ NBLCZ/)
        CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
        CALL h5dcreate_f(GRP_ID, 'positions_greater_up',H5T_NATIVE_DOUBLE, SPACE_ID, DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
        CALL h5dclose_f(DSET_ID, ERROR)
        IF (NUP_ .EQ. NDOWN_) THEN
            NB_POLES(2,2) = NBLCZ
            CALL h5dcreate_f(GRP_ID, 'positions_greater_down',H5T_NATIVE_DOUBLE, SPACE_ID,DSET_ID, ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
            CALL h5dclose_f(DSET_ID, ERROR)
        ENDIF
        DEALLOCATE(VAL)

        D3=(/NW,MAXLCZ,NORB2COMP/)
        D2=(/NBLCZ,NORB2COMP /)
        CALL h5screate_simple_f(3,D3,SPACE_ID,ERROR)
        CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
        START3=(/EXC_STATE-1,0,0/)
        COUNT3=(/1,NBLCZ,NORB2COMP/)
        call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
        CALL h5dopen_f(GRP_ID, 'Q_greater_up', DSET_ID, ERROR)
        call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D3, ERROR, SSPACE_ID, SPACE_ID)
        call h5dclose_f(DSET_ID,ERROR)
        IF (NUP_ .EQ. NDOWN_) THEN
            CALL h5dopen_f(GRP_ID, 'Q_greater_down', DSET_ID, ERROR)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D1, ERROR,SSPACE_ID,SPACE_ID)
            call h5dclose_f(DSET_ID,ERROR)
        endif
        call h5sclose_f(SSPACE_ID, ERROR)
        call h5sclose_f(SPACE_ID, ERROR)
        DEALLOCATE(Q)

        call h5gclose_f(GRP_ID,ERROR)
        call h5fclose_f(FILE_ID,ERROR)
        CALL h5close_f(ERROR)
        IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"

        ! --- GREATER DOWN
        IF (NUP_.NE.NDOWN_) THEN
            NELEC=NELEC_+1
            SZ=SZ_-0.5
            CALL BASIS()
            ALLOCATE(PSI_LCZ(NSTATES,NORB2COMP))
            CALL LCZ_E_DOWN(PSI_LCZ)
            CALL GF_BANDLANCZOS(PSI_LCZ,COEFF_INIT)
            DEALLOCATE(PSI_LCZ)
            ALLOCATE(Q(NBLCZ,NORB2COMP))
            Q(:,1)=COEFF_INIT(1,1)*VEC(1,:)
            DO I = 2,NORB2COMP
                Q(:,I)=COEFF_INIT(I,I)*VEC(I,:)
                DO J = 1,I-1
                    Q(:,I)=Q(:,I)+COEFF_INIT(I,J)*VEC(J,:)
                ENDDO
            ENDDO
            DEALLOCATE(VEC)
            NB_POLES(2,2) = NBLCZ

        ! --- 
            CALL h5open_f(ERROR)
            CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
            call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
            D1=(/ NBLCZ/)
            CALL h5screate_simple_f(1,D1, SPACE_ID, ERROR)
            CALL h5dcreate_f(GRP_ID, 'positions_greater_down',H5T_NATIVE_DOUBLE, SPACE_ID,DSET_ID, ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, VAL, D1, ERROR)
            CALL h5dclose_f(DSET_ID, ERROR)
            DEALLOCATE(VAL)

            D3=(/NW,MAXLCZ,NORB2COMP/)
            D2=(/NBLCZ,NORB2COMP /)
            CALL h5screate_simple_f(3,D3,SPACE_ID,ERROR)
            CALL h5screate_simple_f(2,D2,SSPACE_ID,ERROR)
            START3=(/EXC_STATE-1,0,0/)
            COUNT3=(/1,NBLCZ,NORB2COMP/)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            CALL h5dopen_f(GRP_ID, 'Q_greater_down', DSET_ID, ERROR)
            call h5sselect_hyperslab_f(SPACE_ID, H5S_SELECT_SET_F, START3,COUNT3,ERROR)
            call h5dwrite_f(DSET_ID, H5T_NATIVE_DOUBLE, Q, D1, ERROR,SSPACE_ID,SPACE_ID)
            call h5dclose_f(DSET_ID,ERROR)
            call h5sclose_f(SSPACE_ID, ERROR)
            call h5sclose_f(SPACE_ID, ERROR)
            DEALLOCATE(Q)

            call h5gclose_f(GRP_ID,ERROR)
            call h5fclose_f(FILE_ID,ERROR)
            CALL h5close_f(ERROR)
            IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
        ENDIF
    ENDDO
    CALL h5open_f(ERROR)
    CALL h5fopen_f('hubbard.h5',H5F_ACC_RDWR_F, FILE_ID, ERROR)
    call h5gopen_f(FILE_ID, 'spgf',GRP_ID, ERROR )
    D2=(/2,2/)
    CALL h5screate_simple_f(2,D2, SPACE_ID, ERROR)
    call h5acreate_f(GRP_ID, 'nb_poles', H5T_NATIVE_INTEGER, SPACE_ID, DSET_ID, ERROR)
    call h5awrite_f(DSET_ID, H5T_NATIVE_INTEGER, NB_POLES, D2, ERROR)
    call h5aclose_f(DSET_ID, ERROR)
    CALL h5sclose_f(SPACE_ID, ERROR)
    call h5gclose_f(GRP_ID,ERROR)
    call h5fclose_f(FILE_ID,ERROR)
    CALL h5close_f(ERROR)
    IF (ERROR/=0) WRITE(6,*)" *** Error in solutions hdf5 files"
    RETURN
END SUBROUTINE SPGF_BANDLANCZOS





SUBROUTINE GF_BANDLANCZOS(V,COEFF_INIT)
    USE BASISMOD 
    USE HUBMOD
    IMPLICIT NONE
    REAL*8, PARAMETER :: DF=1.E-10
    REAL*8:: V(NSTATES,NORB2COMP)
    REAL*8,INTENT(OUT) :: COEFF_INIT(NORB2COMP,NORB2COMP)
    REAL*8,ALLOCATABLE :: HT(:,:),H_DEF(:,:), COEFF(:,:)
    REAL*8, ALLOCATABLE :: BETA1(:,:), BETA2(:,:), ALPHA(:)
    INTEGER*4 :: I,J
    INTEGER*4 ::PC_J,CPT,K,NB_DEF
    INTEGER*4,ALLOCATABLE :: DEF(:)
    REAL*8 ::EV, NR
    REAL*8, ALLOCATABLE :: V0(:,:), V2(:)
    LOGICAL :: CV,TEST
    IF (ALLOCATED(VAL)) DEALLOCATE(VAL)
    IF (ALLOCATED(VEC)) DEALLOCATE(VEC)
    ALLOCATE(COEFF(MAXLCZ,NORB2COMP))
    ALLOCATE(V0(NSTATES,NORB2COMP),V2(NSTATES))
    ALLOCATE(BETA1(MAXLCZ,NORB2COMP),BETA2(MAXLCZ,NORB2COMP),ALPHA(MAXLCZ))
    ALLOCATE(DEF(NORB2COMP))
    DEF=0
    COEFF_INIT=0.
    ! --- initialisation of the orthonormal basis set --- !
    NR=NORM2(V(:,1))
    COEFF_INIT(1,1)=NR
    V(:,1)=(1/NR)*V(:,1)
    DO I = 2,NORB2COMP 
        V2=V(:,I)
        DO J = 1,I-1
            COEFF_INIT(I,J)=DOT_PRODUCT(V(:,J),V2)
        ENDDO
        DO J = 1,I-1
            V2=V2-COEFF_INIT(I,J)*V(:,J)
        ENDDO
        NR=NORM2(V2)
        COEFF_INIT(I,I)=NR
        V(:,I)=(1/NR)*V2
    ENDDO
    V0=0.
    CV=.FALSE.
    NBLCZ=0
    EV=0.
    NB_DEF=0
    CPT=0
    DO WHILE(.NOT. CV .AND. NBLCZ < (MAXLCZ-10))
        DO I = 1,NORB2COMP
            NBLCZ=NBLCZ+1
            CALL HPROD(V(:,1),V2)
            ALPHA(NBLCZ)=dot_product(V(:,1),V2)
            DO J = 1,NORB2COMP
                NR=dot_product(V0(:,J),V2)
                V2=V2-NR*V0(:,J)
            ENDDO
            DO J = 1,NORB2COMP-1
                BETA2(NBLCZ,J)=dot_product(V(:,J+1),V2)
                V2=V2-BETA2(NBLCZ,J)*V(:,J+1)
            ENDDO
            V2=V2-ALPHA(NBLCZ)*V(:,1)
            BETA2(NBLCZ,NORB2COMP)=NORM2(V2)
            ! --- deflation test --- !
            IF (BETA2(NBLCZ,NORB2COMP) .LE. DF) THEN
!                WRITE(*,*)NBLCZ+NORB2COMP, BETA2(NBLCZ,NORB2COMP)
                IF (DEF(I)==0) THEN
                    DEF(I)=NBLCZ+NORB2COMP
                    NB_DEF=NB_DEF+1
                ENDIF
            ENDIF
            DO J = 1,NORB2COMP-1
                V0(:,J)=V0(:,J+1)
            ENDDO
            V0(:,NORB2COMP)=V(:,1)
            DO J = 1,NORB2COMP-1
                V(:,J)=V(:,J+1)
            ENDDO
            V(:,NORB2COMP)=(1/BETA2(NBLCZ,NORB2COMP))*V2
        ENDDO
        ALLOCATE(HT(NBLCZ,NBLCZ))
        HT=0.
        DO I = 1,NBLCZ
            HT(I,I)=ALPHA(I)
        ENDDO
        DO I = 1,NBLCZ-NORB2COMP
            DO J = 1,NORB2COMP 
                HT(I+J,I)=BETA2(I,J)
                HT(I,I+J)=HT(I+J,I)  
            ENDDO
        ENDDO 
        DO I = NBLCZ-NORB2COMP+1,NBLCZ-1
            DO J = 1,NBLCZ-I
                HT(I+J,I)=BETA2(I,J)
                HT(I,I+J)=HT(I+J,I)
            ENDDO
        ENDDO   
        CPT=0
        IF (NB_DEF.GE.1) THEN
            ! --- H def --- !
            ALLOCATE(H_DEF(NBLCZ,NBLCZ))
            DO I = 1,NBLCZ 
                TEST=.TRUE.
                DO K = 1,NORB2COMP
                    IF (MOD(I-DEF(K),NORB2COMP)==0 .AND. I.GE.DEF(K) .AND. DEF(K)/=0) THEN
                        CPT=CPT+1
                        TEST=.FALSE.
                    ENDIF
                ENDDO
                IF (TEST) THEN 
                    PC_J=0
                    DO J = 1,NBLCZ
                        PC_J=PC_J+1
                        DO K = 1,NORB2COMP
                            IF (MOD(J-DEF(K),NORB2COMP)==0 .AND. J.GE.DEF(K) .AND. DEF(K)/=0) PC_J=PC_J-1
                        ENDDO
                        H_DEF(PC_J,I-CPT)=HT(J,I)
                    ENDDO
                ENDIF
            ENDDO
            DEALLOCATE(HT)
            ALLOCATE(HT(NBLCZ-CPT,NBLCZ-CPT))
            HT=H_DEF(:NBLCZ-CPT,:NBLCZ-CPT) 
            DEALLOCATE(H_DEF)   
        ENDIF            
        IF (ALLOCATED(VAL)) DEALLOCATE(VAL)
        IF (ALLOCATED(VEC)) DEALLOCATE(VEC)
        ALLOCATE(VAL(NBLCZ-CPT),VEC(NBLCZ-CPT,NBLCZ-CPT))
        CALL DIAGMAT(NBLCZ-CPT,HT,VAL,VEC)
!        WRITE(*,*)NBLCZ,ABS(VAL(1)-EV),BETA2(NBLCZ,NORB2COMP), ALPHA(NBLCZ)
        IF(ABS(VAL(1)-EV)<ACC_LCZ) THEN
            CV=.TRUE.
        ELSE
            EV=VAL(1)
        ENDIF
        DEALLOCATE(HT)
    ENDDO
    IF (NBLCZ>MAXLCZ-10) THEN
        WRITE(*,*) 'FORTRAN BANDLANCZOS ERROR: MAX LCZ STEP REACH'
    ENDIF
    DEALLOCATE(V0,V2,COEFF,ALPHA,BETA1,BETA2,DEF)
    NBLCZ=NBLCZ-CPT
    RETURN 
END SUBROUTINE




SUBROUTINE LCZ_E_UP(PSI_LCZ)
    USE BASISMOD
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(OUT) :: PSI_LCZ(NSTATES,NORB2COMP)
    REAL*8 :: AC
    INTEGER*4 :: A,I,J,K,POS
    PSI_LCZ = 0.
    DO POS = 1,NORB2COMP
        DO I = 1,NSUP_
            A=ISHFT(1,POS-1)
            IF (IAND(BUP_(I),A)==0) THEN
                DO J=1,NSUP
                    IF (BUP(J)==BUP_(I)+A) THEN
                        AC=ANTICOM2(NORB,BUP_(I),POS)
                        DO K = 1,NSDOWN
                            PSI_LCZ((K-1)*NSUP+J,POS)=AC*PSI((K-1)*NSUP_+I)
                        ENDDO
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    RETURN 
END SUBROUTINE


SUBROUTINE LCZ_E_DOWN(PSI_LCZ)
    USE BASISMOD
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(OUT) :: PSI_LCZ(NSTATES,NORB2COMP)
    REAL*8 :: AC
    INTEGER*4 :: A,I,J,K,POS
    PSI_LCZ = 0.
    DO POS = 1,NORB2COMP
        DO I = 1,NSDOWN_
            A=ISHFT(1,POS-1)
            IF (IAND(BDOWN_(I),A)==0) THEN
                DO J=1,NSDOWN
                    IF (BDOWN_(I)+A==BDOWN(J)) THEN
                        AC=ANTICOM2(NORB,BDOWN_(I),POS)
                        DO K = 1,NSUP
                            PSI_LCZ((J-1)*NSUP+K,POS)=AC*PSI((I-1)*NSUP_+K)
                        ENDDO
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    RETURN 
END SUBROUTINE


SUBROUTINE LCZ_H_UP(PSI_LCZ)
    USE BASISMOD
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(OUT) :: PSI_LCZ(NSTATES,NORB2COMP)
    REAL*8 :: AC
    INTEGER*4 :: A,I,J,K,POS
    PSI_LCZ = 0.
    DO POS = 1,NORB2COMP
        DO I = 1,NSUP_
            A=ISHFT(1,POS-1)
            IF (IAND(BUP_(I),A)==A) THEN
                DO J=1,NSUP
                    IF (BUP(J)==BUP_(I)-A) THEN
                        AC=ANTICOM2(NORB,BUP_(I),POS)
                        DO K = 1,NSDOWN
                            PSI_LCZ((K-1)*NSUP+J,POS)=AC*PSI((K-1)*NSUP_+I)
                        ENDDO
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    RETURN 
END SUBROUTINE


SUBROUTINE LCZ_H_DOWN(PSI_LCZ)
    USE BASISMOD
    USE HUBMOD
    USE FUNCMOD
    IMPLICIT NONE
    REAL*8, INTENT(OUT) :: PSI_LCZ(NSTATES,NORB2COMP)
    REAL*8 :: AC
    INTEGER*4 :: A,I,J,K,POS
    PSI_LCZ = 0.
    DO POS = 1,NORB2COMP
        DO I = 1,NSDOWN_
            A=ISHFT(1,POS-1)
            IF (IAND(BDOWN_(I),A)==A) THEN
                DO J=1,NSDOWN
                    IF (BDOWN_(I)-A==BDOWN(J)) THEN
                        AC=ANTICOM2(NORB,BDOWN_(I),POS)
                        DO K = 1,NSUP
                            PSI_LCZ((J-1)*NSUP+K,POS)=AC*PSI((I-1)*NSUP_+K)
                        ENDDO
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    RETURN 
END SUBROUTINE