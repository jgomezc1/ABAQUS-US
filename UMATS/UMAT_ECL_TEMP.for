C23456789012345678901234567890123456789012345678901234567890123456789012
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     UMAT_ECL_TEMP.for                                                C
C     Elasticity CLasical theory                                       C
C     Temperature dependent material properties                        C
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C                                                                      C
C     ISOTROPIC ELASTICITY                                             C
C     CANNOT BE USED FOR PLANE STRESS                                  C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C                                                                      C
C     UNIVERSIDAD EAFIT                                                C
C     LABORATORIO DE MECANICA APLICADA                                 C
C     BLOQUE 14-PISO 2                                                 C
C     MEDELLIN, COLOMBIA                                               C
C                                                                      C
C     LAST UPDATED APRIL 16/2015                                       C
C                                                                      C
C     LOCAL ARRAYS                                                     C
C                                                                      C
C     PROPS(1) - E1                                                    C
C     PROPS(2).- E2                                                    C
C     PROPS(3) - G1                                                    C
C     PROPS(4) - G2                                                    C
C     PROPS(5) - CTE                                                   C
C                                                                      C
C     DDSDDE() - MATERIAL JACOBIAN                                     C
C     STRESS() - UPDATED STRESS VECTOR                                 C
C                                                                      C
C     STATE VARIABLES DEFINITION                                       C
C                                                                      C
C      1-4: Elastic strain vector                                      C
C      5-8: Thermal strain vector                                      C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,
     2DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
     3PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0,
     4DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C 
      DIMENSION STRESS(NTENS), STATEV(NSTATV),DDSDDE(NTENS, NTENS),
     1DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),
     2PREDEF(1), DPRED(1), PROPS(NPROPS),COORDS(3),DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3)
C
      DIMENSION DS(NTENS),DSTRTHER(NTENS),DSTRANM(NTENS),EELAS(NTENS),
     1ETHER(NTENS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=30,
     2           FOUR=4.D0)
C
C**********************************************************************
C               C L A S S I C A L    E L A S T I C I T Y
C**********************************************************************
C
      DO K1=1, NTENS
        EELAS(K1)=STATEV(K1)
        ETHER(K1)=STATEV(K1+NTENS)
      END DO
C
C     Temperature at the end of the step
C
      THETA=TEMP
      DTHETA=DTEMP
C
C     Temperature Dependent Elastic properties
C
      EMOD=(PROPS(1)+PROPS(2)*THETA)*1000.0
C                              Shear modulus
      GMOD=(PROPS(3)+PROPS(4)*THETA)*1000.0
C                              Poisson's ratio
      ENU=(EMOD/(TWO*GMOD))-ONE
C                              Coefficient of thermal expansion
      CTE=PROPS(5)
C
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=two*GMOD
      EG=GMOD
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
C
C     Elastic stiffness
C
      DO K1=1, 3
        DO K2=1, 3
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=4, 4
        DDSDDE(K1, K1)=EG
      END DO
C
C     Thermal strain increment
C
      CALL KCLEAR(DSTRTHER,NTENS,1)
      DO K1=1,NDI
        DSTRTHER(K1)=CTE*DTHETA
      END DO
C
C     Mechanical strain increment
C
      DO K1=1,NDI
        DSTRANM(K1)=DSTRAN(K1)-DSTRTHER(K1)
      END DO
      DSTRANM(NDI+1)=DSTRAN(NDI+1)
C
C     Elastic strain
C
      DO K1=1,NDI
        EELAS(K1)=EELAS(K1)+DSTRANM(K1)
        ETHER(K1)=ETHER(K1)+DSTRTHER(K1)
      END DO
C
C
C     Calculate the stress increment and
C     updates the stress vector.
C
      CALL KMAVEC(DDSDDE,NTENS,NTENS,DSTRANM,DS)
      CALL KUPDVEC(STRESS,NTENS,DS)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C             M A T R I X   H A N D L I N G                            C
C-------------U T I L I T I E S   B L O C K--------------              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KCLEAR(A,N,M)                                        C
C      Clear a real matrix                                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KCLEAR(A,N,M)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
      DIMENSION A(N,M)
C
      DO I=1,N
        DO J=1,M
          A(I,J)=ZERO
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)                        C
C      Real matrix product                                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.D0)
      DIMENSION A(NRA,NCA),B(NRB,NCB),C(NRA,NCB)
C
      CALL KCLEAR(C,NRA,NCB)
      DUM=ZERO
      DO I=1,NRA
        DO J=1,NCB
         DO K=1,NCA
           DUM=DUM+A(I,K)*B(K,J)
          END DO
          C(I,J)=DUM
          DUM=ZERO
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KSMULT(A,NR,NC,S)                                    C
C      Matrix times a scalar.                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSMULT(A,NR,NC,S)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NR,NC)
C
      DO I=1,NR
        DO J=1,NC
          DUM=A(I,J)
          A(I,J)=S*DUM
          DUM=0.D0
        END DO  
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KUPDMAT(A,NR,NC,B)                                   C
C      Updates an existing matrix with an incremental matrix.          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KUPDMAT(A,NR,NC,B)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.D0)
C
      DIMENSION A(NR,NC),B(NR,NC)
C
      DO I=1,NR
        DO J=1,NC
          DUM=A(I,J)
          A(I,J)=ZERO
          A(I,J)=DUM+B(I,J)
          DUM=ZERO
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMTRAN(A,NRA,NCA,B)                                  C      
C      Matrix transpose                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KMTRAN(A,NRA,NCA,B)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NRA,NCA),B(NCA,NRA)
C
      CALL KCLEAR(B,NCA,NRA)
      DO I=1,NRA
       DO J=1,NCA
         B(J,I)=A(I,J)
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)                                C
C      Real matrix times vector                                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.D0)
      DIMENSION A(NRA,NCA),B(NCA),C(NRA)
C
      CALL KCLEARV(C,NRA)
C
      DO K1=1,NRA
        DO K2=1,NCA
          C(K1)=C(K1)+A(K1,K2)*B(K2)	    
        END DO
      END DO     
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KCLEARV(A,N)                                         C
C      Clear a real vector                                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KCLEARV(A,N)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
C
      DIMENSION A(N)
C
      DO I=1,N
        A(I)=ZERO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KUPDVEC(A,NR,B)                                      C
C      Updates an existing vector with an incremental vector.          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KUPDVEC(A,NR,B)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.D0)
C
      DIMENSION A(NR),B(NR)
C
      DO I=1,NR
        DUM=A(I)
        A(I)=ZERO
        A(I)=DUM+B(I)
        DUM=ZERO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KVECSUB(A,NRA,B,NRB,C)                               C
C      Substracts one column vector from another column vector         C
C      IFLAG=0 for substraction                                        C
C      IFLAG=1 for addition                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KVECSUB(A,NRA,B,NRB,C,IFLAG)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0, ONENEG=-1.0D0)
C
      DIMENSION A(NRA,1),B(NRB,1),C(NRB,1)
C
      SCALAR=ONENEG
C
      IF (IFLAG.EQ.1) SCALAR=ONE
C
      DO I=1,NRA
        C(I,1)=A(I,1)+B(I,1)*SCALAR
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMATSUB(A,NRA,NCA,B,C,IFLAG)                         C
C      Substracts one rectangular matrix from another rectangular      C
C      matrix                                                          C
C      IFLAG=0 for substraction                                        C
C      IFLAG=1 for addition                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KMATSUB(A,NRA,NCA,B,C,IFLAG)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0, ONENEG=-1.0D0)
C
      DIMENSION A(NRA,NCA),B(NRA,NCA),C(NRA,NCA)
C
      CALL KCLEAR(C,NRA,NCA)
C
      SCALAR=ONENEG
C
      IF (IFLAG.EQ.1) SCALAR=ONE
C
      DO I=1,NRA
        DO J=1,NCA
          C(I,J)=A(I,J)+B(I,J)*SCALAR
        END DO
      END DO
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE IDENTITY                                              C
C     CREATES AN IDENTITY MATRIX OF DIMENSIONS NDIM,NDIM               C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KIDENTITY(DEL,NDIM)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0)
C
      DIMENSION DEL(NDIM,NDIM)
C
      CALL KCLEAR(DEL,NDIM,NDIM)
C
      DO K1=1,NDIM
        DEL(K1,K1)=ONE
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KINVERSE                                                 C
C                                                                      C
C   IVEERSE OF A MATRIX USING LU DECOMPOSITION                         C
C   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
C                                                                      C
C   A   Matrix to be inverted.                                         C
C   Y   Inverse of A                                                   C
C   N   Dimension                                                      C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KINVERSE(A,Y,NP,N)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      DIMENSION A(NP,NP),Y(NP,NP),INDX(NP),AUX(NP,NP)
C
      CALL KCLEAR(AUX,NP,NP)
      CALL KCOPYMAT(A,AUX,N)
C
      DO I=1,N
        DO J=1,N
          Y(I,J)=ZERO
        END DO
        Y(I,I)=ONE
      END DO
      CALL KLUDCMP(AUX,N,NP,INDX,D)
      DO J=1,N
        CALL KLUBKSB(AUX,N,NP,INDX,Y(1,J))
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KLUDCMP                                                 C
C                                                                      C
C   LU MATRIX DECOMPOSITION                                            C
C   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KLUDCMP(A,N,NP,INDX,D)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NMAX=500,TINY=1.0E-20,ZERO=0.D0,ONE=1.D0)
C
      DIMENSION INDX(N),A(NP,NP),VV(NMAX)
C
      D=ONE
      DO I=1,N
        AAMAX=ZERO
        DO J=1,N
          IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        END DO
        IF(AAMAX.EQ.0.) PAUSE 'SINGULAR MATRIX IN LUDCMP'
        VV(I)=ONE/AAMAX
      END DO
C
      DO J=1,N
        DO I=1,J-1
          SUM=A(I,J)
          DO K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
        END DO
        AAMAX=ZERO
        DO I=J,N
          SUM=A(I,J)
          DO K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
          DUM=VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX) THEN
             IMAX=I
             AAMAX=DUM
          END IF
        END DO
        IF(J.NE.IMAX) THEN
           DO K=1,N
             DUM=A(IMAX,K)
             A(IMAX,K)=A(J,K)
             A(J,K)=DUM
           END DO
           D=-D
           VV(IMAX)=-VV(J)
        END IF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.) A(J,J)=TINY
        IF(J.NE.N) THEN
           DUM=ONE/A(J,J)
           DO I=J+1,N
             A(I,J)=A(I,J)*DUM
           END DO
        END IF
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KLUBKSB                                                 C
C                                                                      C
C   FORWARD SUBSTITUTION                                               C
C   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KLUBKSB(A,N,NP,INDX,B)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0)
C
      DIMENSION INDX(N),A(NP,NP),B(NP)
C
      II=0
      DO I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF(II.NE.0) THEN
           DO J=II,I-1
             SUM=SUM-A(I,J)*B(J)
           END DO
        ELSE IF(SUM.NE.ZERO) THEN
           II=I
        END IF
        B(I)=SUM
      END DO
C
      DO I=N,1,-1
        SUM=B(I)
        DO J=I+1,N
          SUM=SUM-A(I,J)*B(J)
        END DO
        B(I)=SUM/A(I,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE KCOPYMAT                                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KCOPYMAT(A,B,N)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION A(N,N),B(N,N)
C
      CALL KCLEAR(B,N,N)
C
      DO K1=1,N
        DO K2=1,N
          B(K1,K2)=A(K1,K2)
        END DO
      END DO
C
      RETURN
C
      END