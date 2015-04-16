CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     UMAT_PCLI.for     Plasticity Classical     Theory                C
C     Requires a UEL.for subroutine                                    C
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C     RATE INDEPENDENT ISOTROPIC HARDENING PLASTICITY                  C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C     UNIVERSIDAD EAFIT                                                C
C     LABORATORIO DE MECANICA APLICADA                                 C
C     BLOQUE 14-PISO 2                                                 C
C     MEDELLIN, COLOMBIA                                               C
C                                                                      C
C     LAST UPDATED APRIL 16/2015                                       C
C                                                                      C
C     LOCAL ARRAYS                                                     C
C                                                                      C
C     EELAS - ELASTIC STRAINS STATEV(1..NTENS)                         C
C     EPLAS - PLASTIC STRAINS STATEV(NTENS+1...2*NTENS)                C
C     FLOW  - PLASTIC FLOW DIRECTION                                   C
C     HARD  - HARDENING MODULUS                                        C
C                                                                      C
C     PROPS(1) - E                                                     C
C     PROPS(2) - NU                                                    C
C     PROPS(5..) - SYIELD AN HARDENING DATA                            C
C     CALLS KUHARD FOR CURVE OF YIELD STRESS VS. PLASTIC STRAINS       C
C     EQPLAS - EQUIVALENT PLASTIC STRAIN STATEV(2*NTENS+1)             C
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
     3DFGRD0(3, 3), DFGRD1(3, 3)
C
      DIMENSION EELAS(NTENS), EPLAS(NTENS),DS(NTENS),DSTRESS(NTENS)
C
C
      DIMENSION AUX1(1,NTENS),AUX2(NTENS,NTENS),AUX3(NTENS,NTENS),
     1AUX4(NTENS,NTENS),AUX5(NTENS,NTENS),AUX6(NTENS,1),AUX7(1,NTENS),
     2AUX8(1,1),AUX9(1,NTENS),AUX10(1,1),AUX11(NTENS,NTENS),
     3AUX12(NTENS,1),AUX13(1,NTENS),AUX14(1,1),AUX15(1,NTENS),
     4AUX16(NTENS,NTENS),AUX17(NTENS,NTENS),STRESST(1,NTENS),
     5P(NTENS,NTENS),SINVAR(1,1),BI(NTENS,NTENS),STRESSUPD(NTENS,1),
     6SDEV(NTENS,1),EM(NTENS,NTENS),B(NTENS,NTENS),Q(NTENS,NTENS),
     7QT(NTENS,NTENS),DP(NTENS,NTENS),DC(NTENS,NTENS),DEL(NTENS,NTENS),
     8BT(NTENS,NTENS),DIAG(NTENS,NTENS),GDIA(NTENS,NTENS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=30,
     2           FOUR=4.D0)
C
C     Inititalize arrays
C
      CALL KCLEAR(AUX1,1,NTENS)
      CALL KCLEAR(AUX2,NTENS,NTENS)
      CALL KCLEAR(AUX3,NTENS,NTENS)
      CALL KCLEAR(AUX4,NTENS,NTENS)
      CALL KCLEAR(AUX5,NTENS,NTENS)
      CALL KCLEAR(AUX6,NTENS,1)
      CALL KCLEAR(AUX7,1,NTENS)
      CALL KCLEAR(AUX8,1,1)
      CALL KCLEAR(AUX9,1,NTENS)
      CALL KCLEAR(AUX10,1,1)
      CALL KCLEAR(AUX11,NTENS,NTENS)
      CALL KCLEAR(AUX12,NTENS,1)
      CALL KCLEAR(AUX13,1,NTENS)
      CALL KCLEAR(AUX14,1,1)
      CALL KCLEAR(AUX15,1,NTENS)
      CALL KCLEAR(AUX16,NTENS,NTENS)
      CALL KCLEAR(AUX17,NTENS,NTENS)
C
      CALL KCLEAR(P,NTENS,NTENS)
      CALL KCLEAR(Q,NTENS,NTENS)
      CALL KCLEAR(DP,NTENS,NTENS)
      CALL KCLEAR(DC,NTENS,NTENS)
      CALL KCLEAR(QT,NTENS,NTENS)
      CALL KCLEAR(DEL,NTENS,NTENS)
C
      CALL KCLEAR(STRESST,1,NTENS)
      CALL KCLEAR(SINVAR,1,1)
C
      CALL KCLEAR(STRESSUPD,NTENS,1)
      CALL KCLEAR(AUX17,NTENS,NTENS)
      CALL KCLEAR(B,NTENS,NTENS)
      CALL KCLEAR(BT,NTENS,NTENS)
C 
C     Recover equivalent plastic strain, elastic strains, and plastic
C     strains. Also initialize user definde data sets.
C
      DO K1=1, NTENS
        EELAS(K1)=STATEV(K1)
        EPLAS(K1)=STATEV(K1+NTENS)
      END DO
      EQPLAS=STATEV(1+2*NTENS)
C
C     Elastic properties
C
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2),ENUMAX)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
      CBETA2=EG2
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
C     Form proyector and diagonal decomposition matrices
C
      CALL KPROYECTOR(P)
C
C     Calculate predictor stress and elastic strains
C
      CALL KMAVEC(DDSDDE,NTENS,NTENS,DSTRAN,DS)
      CALL KUPDVEC(STRESS,NTENS,DS)
      CALL KUPDVEC(EELAS,NTENS,DSTRAN)
C
C     Calculate equivalent mises stress
C
      CALL KMTRAN(STRESS,NTENS,1,STRESST)
      CALL KMMULT(P,NTENS,NTENS,STRESS,NTENS,1,SDEV)
      CALL KMMULT(STRESST,1,NTENS,SDEV,NTENS,1,SINVAR)
      FBAR=DSQRT(SINVAR(1,1))
C
C     Get the yield stress from the specifid hardening function.
C
      CALL KUHARD(SYIEL0,EHARD,EQPLAS,2,PROPS(3))
C
C     Determine if actively yielding
C
      SYIELD=SYIEL0
      IF(FBAR.GT.(ONE+TOLER)*SYIEL0) THEN
C
C       Actively yielding-Perform local Newton iterations
C       to find consistncy parameter and equivalent plastic
C       strain
C
C       Starts iterations
C
        ITER=1
        GAM_PAR=ZERO
        IFLAG=0
        DO
          CALL KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,EMOD,ENU)
          CALL KMTRAN(Q,NTENS,NTENS,QT)
          CALL KMMULT(STRESST,1,NTENS,Q,NTENS,NTENS,AUX1)
          CALL KMMULT(QT,NTENS,NTENS,P,NTENS,NTENS,AUX2)
          CALL KMMULT(AUX2,NTENS,NTENS,Q,NTENS,NTENS,AUX3)
          CALL KMMULT(AUX3,NTENS,NTENS,DIAG,NTENS,NTENS,AUX4)
          CALL KMMULT(AUX4,NTENS,NTENS,QT,NTENS,NTENS,AUX5)
          CALL KMMULT(AUX5,NTENS,NTENS,STRESS,NTENS,1,AUX6)
          CALL KMMULT(AUX1,1,NTENS,DIAG,NTENS,NTENS,AUX7)
          CALL KMMULT(AUX7,1,NTENS,AUX6,NTENS,1,AUX8)
          CALL KMMULT(AUX1,1,NTENS,GDIA,NTENS,NTENS,AUX9)
          CALL KMMULT(AUX9,1,NTENS,AUX6,NTENS,1,AUX10)
          FBAR=DSQRT(AUX8(1,1))
          TETA2=ONE-(TWO/THREE)*EHARD*GAM_PAR
          FJAC=TETA2*AUX10(1,1)/FBAR-(TWO/THREE)*EHARD*FBAR
          FGAM=FBAR-SYIELD
C
C         Updates
C
          GAM_PAR=GAM_PAR-FGAM/FJAC
          EQPLAS1=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR*FBAR
          CALL KUHARD(SYIELD,EHARD,EQPLAS1,2,PROPS(3))
C
          IF(ABS(FGAM/FJAC).LT.TOLER) THEN
            IFLAG=0
            GOTO 801
          ELSE
            IF(ITER.GT.MAXITER) THEN
              IFLAG=1
              GOTO 802
            END IF
          END IF
C
          ITER=ITER+1
        END DO
C
  801   CONTINUE
C
C       Local Newton algorithm converged
C       Update stresses, elastic and plastic strains, equivalent plastic
C       strains
C
        CALL KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,EMOD,ENU)
        CALL KMMULT(Q,NTENS,NTENS,DIAG,NTENS,NTENS,AUX17)
        CALL KMMULT(AUX17,NTENS,NTENS,QT,NTENS,NTENS,B)
        CALL KMMULT(B,NTENS,NTENS,STRESS,NTENS,1,STRESSUPD)
C
        CALL KCLEAR(STRESS,NTENS,1)
        DO K1=1,NTENS
          STRESS(K1)=STRESSUPD(K1,1)
        END DO
C
        CALL KCLEAR(STRESST,1,NTENS)
        CALL KCLEAR(SDEV,NTENS,1)
        CALL KMTRAN(STRESS,NTENS,1,STRESST)
        CALL KMMULT(P,NTENS,NTENS,STRESS,NTENS,1,SDEV)
        CALL KMMULT(STRESST,1,NTENS,SDEV,NTENS,1,SINVAR)
        FBAR=DSQRT(SINVAR(1,1))
C
        DO K1=1,NTENS
          EPLAS(K1)=EPLAS(K1)+GAM_PAR*SDEV(K1,1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
        END DO
C
        EQPLAS=EQPLAS1
C
C       Formulate the consistent material Jacobian (tangent)
C
        CALL KCLEAR(EM,NTENS,NTENS)
        CALL KMMULT(B,NTENS,NTENS,DDSDDE,NTENS,NTENS,EM)
        CALL KMMULT(EM,NTENS,NTENS,P,NTENS,NTENS,AUX11)
        CALL KMMULT(AUX11,NTENS,NTENS,STRESS,NTENS,1,AUX12)
        CALL KMMULT(STRESST,1,NTENS,P,NTENS,NTENS,AUX13)
        CALL KMMULT(AUX13,1,NTENS,AUX12,NTENS,1,AUX14)
        CALL KMTRAN(AUX12,NTENS,1,AUX15)
        CALL KMMULT(AUX12,NTENS,1,AUX15,1,NTENS,AUX16)
        SCALAR1=ONE/AUX14(1,1)
        CALL KSMULT(AUX16,NTENS,NTENS,SCALAR1)
        TETA2=ONE-(TWO/THREE)*EHARD*GAM_PAR
        CBETA=(TWO/THREE/TETA2/AUX14(1,1))*FBAR*FBAR*EHARD
        SCALAR2=ONE/(ONE+CBETA)
        CALL KSMULT(AUX16,NTENS,NTENS,SCALAR2)
        CALL KCLEAR(DDSDDE,NTENS,NTENS)
        CALL KMATSUB(EM,NTENS,NTENS,AUX16,DDSDDE,0)
C
      END IF
C
C     Store elastic strains, (equivalent) plastic strains
C     in state variable array
C
      DO K1=1,NTENS
        STATEV(      K1)=EELAS(K1)
        STATEV(NTENS+K1)=EPLAS(K1)
      END DO
      STATEV(2*NTENS+1)=EQPLAS
      STATEV(2*NTENS+2)=DSQRT(THREE/TWO)*FBAR
C
  802 IF (IFLAG.EQ.1) THEN
         WRITE(*,*)
         WRITE(*,*) 'LOCAL PLASTICITY ALGORITHM DID NOT CONVREGED'
         WRITE(*,*) 'AT GAUSS POINT=',NPT, 'ELEMENT=',NOEL
         WRITE(*,*) 'AFTER=',ITER,' ITERATIONS'
         WRITE(*,*) 'LAST CORRECTION=',fgam/fjac
         CALL XIT
      END IF
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE UHARD                                                 C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KUHARD(SYIELD,EHARD,EQPLAS,NVALUE,TABLE)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION TABLE(2,NVALUE)
C
      PARAMETER(ZERO=0.D0,TWO=2.D0,THREE=3.D0)
C
      SYIEL0=TABLE(1,1)
C
C     Compute hardening modulus
C
      EHARD=(TABLE(1,2)-TABLE(1,1))/TABLE(2,2)
C
C     Compute yield stress corresponding to EQPLAS
C
      SYIELD=DSQRT(TWO/THREE)*(SYIEL0+EHARD*EQPLAS)
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE SPECTRAL                                              C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,E,ENU)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1          SIX=6.D0)
C
      DIMENSION Q(4,4),DP(4,4),DC(4,4),DIAG(4,4),GDIA(4,4)
C
      CALL KCLEAR(Q,4,4)
      CALL KCLEAR(DP,4,4)
      CALL KCLEAR(DC,4,4)
      CALL KCLEAR(DIAG,4,4)
      CALL KCLEAR(GDIA,4,4)
C
      EG2=E/(ONE+ENU)
      EG=EG2/TWO
      ELAM=EG2*ENU/(ONE-TWO*ENU)
      CBETA2=EG2
      CALFA2=EG2
C
      Q(1,1)=ZERO
      Q(1,2)=TWO/DSQRT(SIX)
      Q(1,3)=ONE/DSQRT(THREE)
      Q(2,1)=-DSQRT(TWO)/TWO
      Q(2,2)=-ONE/DSQRT(SIX)
      Q(2,3)=ONE/DSQRT(THREE)
      Q(3,1)=DSQRT(TWO)/TWO
      Q(3,2)=-ONE/DSQRT(SIX)
      Q(3,3)=ONE/DSQRT(THREE)
      Q(4,4)=ONE
C
      DP(1,1)=ONE
      DP(2,2)=ONE
      DP(3,3)=ZERO
      DP(4,4)=TWO
C
      DC(1,1)=EG2
      DC(2,2)=EG2
      DC(3,3)=THREE*ELAM+EG2
      DC(4,4)=EG2
C
      DIAG(1,1)=ONE/(ONE+EG2*GAM_PAR)
      DIAG(2,2)=ONE/(ONE+EG2*GAM_PAR)
      DIAG(3,3)=ONE
      DIAG(4,4)=ONE/(ONE+EG2*GAM_PAR)
      GDIA(1,1)=-(EG2/((ONE+EG2*GAM_PAR)**2))
      GDIA(2,2)=-(EG2/((ONE+EG2*GAM_PAR)**2))
      GDIA(3,3)=ZERO
      GDIA(4,4)=-(EG2/((ONE+EG2*GAM_PAR)**2))
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE PROYECTOR                                             C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KPROYECTOR(P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0)
C
      DIMENSION P(4,4)
C
      CALL KCLEAR(P,4,4)
C
      P(1,1)=TWO/THREE
      P(1,2)=-ONE/THREE
      P(1,3)=-ONE/THREE
      P(2,1)=-ONE/THREE
      P(2,2)=TWO/THREE
      P(2,3)=-ONE/THREE
      P(3,1)=-ONE/THREE
      P(3,2)=-ONE/THREE
      P(3,3)=TWO/THREE
      P(4,4)=TWO
C
      RETURN
C
      END
C
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)                        C
C      Real matrix product                                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KSMULT(A,NR,NC,S)                                    C
C      Matrix times a scalar.                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KUPDMAT(A,NR,NC,B)                                   C
C      Updates an existing matrix with an incremental matrix.          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMTRAN(A,NRA,NCA,B)                                  C      
C      Matrix transpose                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)                                C
C      Real matrix times vector                                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KCLEARV(A,N)                                         C
C      Clear a real vector                                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KUPDVEC(A,NR,B)                                      C
C      Updates an existing vector with an incremental vector.          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
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