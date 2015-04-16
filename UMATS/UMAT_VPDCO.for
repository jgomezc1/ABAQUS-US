CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     UMAT_VPDCO.for    Viscoplasticity Damage Coupled Cosserat Couple C
C                       Stress Theory                                  C
C     Requires a UEL.for subroutine                                    C
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C     DAMAGE MECHANICS COUPLED RATE DEPENDENT HARDENING PLASTICITY     C
C     COMBINED ISOTROPIC/KINEMATIC HARDENING                           C
C     RETURN MAPPING ALGORITHM                                         C
C                                                                      C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C     CREATED BY JUAN GOMEZ                                            C
C     UNIVERSIDAD EAFIT                                                C
C     LABORATORIO DE MECANICA APLICADA                                 C
C     BLOQUE 14-PISO 2                                                 C
C     MEDELLIN, COLOMBIA                                               C
C                                                                      C
C     LAST UPDATED APRIL 16/2015                                       C
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
      DIMENSION EELAS(NTENS), EPLAS(NTENS),DS(NTENS),DSTRESS(NTENS)
C
      DIMENSION AUX1(1,NTENS),AUX2(NTENS,NTENS),AUX3(NTENS,NTENS),
     1AUX4(NTENS,NTENS),AUX5(NTENS,NTENS),AUX6(NTENS,1),AUX7(1,NTENS),
     2AUX8(1,1),AUX9(1,NTENS),AUX10(1,1),AUX11(NTENS,NTENS),
     3AUX12(NTENS,1),AUX13(1,NTENS),AUX14(1,1),AUX15(1,NTENS),
     4AUX16(NTENS,NTENS),AUX17(NTENS,NTENS),STSRELT(1,NTENS),
     5P(NTENS,NTENS),SINVAR(1,1),BI(NTENS,NTENS),STRESSUPD(NTENS,1),
     6SDEV(NTENS,1),EM(NTENS,NTENS),B(NTENS,NTENS),Q(NTENS,NTENS),
     7QT(NTENS,NTENS),DP(NTENS,NTENS),DC(NTENS,NTENS),DEL(NTENS,NTENS),
     8BT(NTENS,NTENS),DIAG(NTENS,NTENS),GDIA(NTENS,NTENS),
     9PM(NTENS,NTENS),XBACK(NTENS),STSREL(NTENS),EPL(NTENS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=16)
C
C**********************************************************************
C               C O S S E R A T  E L A S T I C I T Y
C**********************************************************************
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
      CALL KCLEAR(STSREL,NTENS,1)
      CALL KCLEAR(P,NTENS,NTENS)
      CALL KCLEAR(Q,NTENS,NTENS)
      CALL KCLEAR(DP,NTENS,NTENS)
      CALL KCLEAR(DC,NTENS,NTENS)
      CALL KCLEAR(QT,NTENS,NTENS)
C 
C     Recover state variables
C     Elastic strain,Plastic strain,equivalent plastic strain
C     Back stress,Damage, Helmholtz free energy term.
C
      DO K1=1, NTENS
        EELAS(K1)=STATEV(K1)
        EPLAS(K1)=STATEV(K1+NTENS)
      END DO
      EQPLAS=STATEV(1+2*NTENS)
      JJ=1
      DO II=14,2*NTENS+7
        XBACK(JJ)=STATEV(II)
        JJ=JJ+1
      END DO
      D=STATEV(20)
      ESTE=STATEV(21)
      SEQUIV=STATEV(22)
C
C     Compute viscosity prameter sneta and assimilates material properties
C
      CALL VSPRATE(PROPS,NPROPS,EMOD,ENU,SNETA,SIG0,SIGSAT,
     1             HMOD,HRDRATE)
C
C     Elastic properties
C
      SN=PROPS(20)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
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
      DDSDDE(5,5)=EG2
      DDSDDE(6,6)=EG2
C
C     Form proyector and diagonal decomposition matrices
C
      CALL KPROYECTOR(P)
C
C     Calculate predictor stress(using damage) and elastic strains
C
      CALL KMAVEC(DDSDDE,NTENS,NTENS,DSTRAN,DS)
      CALL KSMULT(DS,NTENS,1,(ONE-D))
      CALL KUPDVEC(STRESS,NTENS,DS)
      CALL KUPDVEC(EELAS,NTENS,DSTRAN)
C
      DO K1=1,NTENS
        STSREL(K1)=STRESS(K1)-XBACK(K1)
      END DO
C
C     Calculate equivalent mises stress
C
      CALL KCLEAR(STSRELT,1,NTENS)
      CALL KCLEAR(SINVAR,1,1)
      CALL KMTRAN(STSREL,NTENS,1,STSRELT)
      CALL KMMULT(P,NTENS,NTENS,STSREL,NTENS,1,SDEV)
      CALL KMMULT(STSRELT,1,NTENS,SDEV,NTENS,1,SINVAR)
      FBAR=DSQRT(SINVAR(1,1))
C
C     Get the yield stress from the specifid hardening function.
C
      CALL KUHARDNLIN(SYIEL0,SYIELDK,HMOD,SIG0,SIGSAT,HRDRATE,
     1                EHARDI,EHARDK,EQPLAS)
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
        GAM_PAR=ZERO+1.0d-25
        PHIINV=0.0D0
        DETDG=0.D0
        IFLAG=0
        DO
          IF(SNETA.NE.0.0) THEN
             PHIINV=((SNETA/DTIME)**(ONE/SN))*(GAM_PAR**(ONE/SN))
             DETDG=(ONE/SN)*((SNETA/DTIME)**(ONE/SN))
     $             *(GAM_PAR**((ONE/SN)-ONE))
          END IF
          CALL KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,EMOD,ENU,EHARDK,D)
          CALL KMTRAN(Q,NTENS,NTENS,QT)
          CALL KMMULT(STSRELT,1,NTENS,Q,NTENS,NTENS,AUX1)
          CALL KMMULT(QT,NTENS,NTENS,P,NTENS,NTENS,AUX2)
          CALL KMMULT(AUX2,NTENS,NTENS,Q,NTENS,NTENS,AUX3)
          CALL KMMULT(AUX3,NTENS,NTENS,DIAG,NTENS,NTENS,AUX4)
          CALL KMMULT(AUX4,NTENS,NTENS,QT,NTENS,NTENS,AUX5)
          CALL KMMULT(AUX5,NTENS,NTENS,STSREL,NTENS,1,AUX6)
          CALL KMMULT(AUX1,1,NTENS,DIAG,NTENS,NTENS,AUX7)
          CALL KMMULT(AUX7,1,NTENS,AUX6,NTENS,1,AUX8)
          CALL KMMULT(AUX1,1,NTENS,GDIA,NTENS,NTENS,AUX9)
          CALL KMMULT(AUX9,1,NTENS,AUX6,NTENS,1,AUX10)
          FBAR=DSQRT(AUX8(1,1))
          TETA2=ONE-(TWO/THREE)*EHARDI*GAM_PAR*(ONE-D)
          FJAC=TETA2*AUX10(1,1)/FBAR-(TWO/THREE)*EHARDI*(ONE-D)
     &    *FBAR-DETDG
          FGAM=FBAR-SYIELD-PHIINV
C
C         Updates
C
          GAM_PAR=GAM_PAR-FGAM/FJAC
          EQPLAS1=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR*(ONE-D)*FBAR
          CALL KUHARDNLIN(SYIELD,SYIELDK,HMOD,SIG0,SIGSAT,HRDRATE,
     1                    EHARDI,EHARDK,EQPLAS1)
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
C       strains, Back stress.
C
        CALL KCLEAR(STRESSUPD,NTENS,1)
        CALL KCLEAR(AUX17,NTENS,NTENS)
        CALL KCLEAR(B,NTENS,NTENS)
        CALL KCLEAR(BT,NTENS,NTENS)
        CALL KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,EMOD,ENU,EHARDK,D)
        CALL KMMULT(Q,NTENS,NTENS,DIAG,NTENS,NTENS,AUX17)
        CALL KMMULT(AUX17,NTENS,NTENS,QT,NTENS,NTENS,B)
        CALL KMTRAN(B,NTENS,NTENS,BT)
        CALL KMMULT(B,NTENS,NTENS,STSREL,NTENS,1,STRESSUPD)
C
        CALL KCLEAR(STSREL,NTENS,1)
        DO K1=1,NTENS
          STSREL(K1)=STRESSUPD(K1,1)
          XBACK(K1 )=XBACK(K1)+(TWO/THREE)*GAM_PAR*EHARDK*(ONE-D)
     &    *STSREL(K1)
          STRESS(K1)=STSREL(K1)+XBACK(K1)
        END DO
C
        CALL KCLEAR(STSRELT,1,NTENS)
        CALL KCLEAR(SDEV,NTENS,1)
        CALL KMTRAN(STSREL,NTENS,1,STSRELT)
        CALL KMMULT(P,NTENS,NTENS,STSREL,NTENS,1,SDEV)
        CALL KMMULT(STSRELT,1,NTENS,SDEV,NTENS,1,SINVAR)
        FBAR=DSQRT(SINVAR(1,1))
C
        DO K1=1,NTENS
          EPL(K1)=GAM_PAR*SDEV(K1,1)/(ONE-D)
          EPLAS(K1)=EPLAS(K1)+(GAM_PAR*SDEV(K1,1)/(ONE-D))
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
        END DO
        EQPLAS=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR*(ONE-D)*FBAR
C
C       Formulate the consistent material Jacobian (tangent)
C
        CALL KDAMACAL(STRESS,EPL,NTENS,PROPS,NPROPS,ESTE,D)
        CALL KCLEAR(EM,NTENS,NTENS)
        TETA1=ONE+(TWO/THREE)*EHARDK*GAM_PAR*(ONE-D)
        TETA2=ONE-(TWO/THREE)*EHARDI*GAM_PAR*(ONE-D)
        CALL KSMULT(B,NTENS,NTENS,TETA1)
        CALL KMMULT(B,NTENS,NTENS,DDSDDE,NTENS,NTENS,EM)
        CALL KMMULT(EM,NTENS,NTENS,P,NTENS,NTENS,AUX11)
        CALL KMMULT(AUX11,NTENS,NTENS,STSREL,NTENS,1,AUX12)
        CALL KMMULT(STSRELT,1,NTENS,P,NTENS,NTENS,AUX13)
        CALL KMMULT(AUX13,1,NTENS,AUX12,NTENS,1,AUX14)
        CALL KMTRAN(AUX12,NTENS,1,AUX15)
        CALL KMMULT(AUX12,NTENS,1,AUX15,1,NTENS,AUX16)
        SCALAR1=ONE/AUX14(1,1)
        CALL KSMULT(AUX16,NTENS,NTENS,SCALAR1)
        VPCONS=((TETA1**2)*FBAR)/(TETA2*AUX14(1,1))
        CBETA=(TWO*TETA1*FBAR*FBAR*(ONE-D)*(EHARDI*TETA1+EHARDK*TETA2))
     &  /(THREE*TETA2*AUX14(1,1))
        CBETA=CBETA+VPCONS
        SCALAR2=ONE/(ONE+CBETA)
        CALL KSMULT(AUX16,NTENS,NTENS,SCALAR2)
        CALL KCLEAR(DDSDDE,NTENS,NTENS)
        CALL KMATSUB(EM,NTENS,NTENS,AUX16,DDSDDE,0)
        CALL KSMULT(DDSDDE,NTENS,NTENS,(ONE-D))
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
      JJ=1
      DO II=14,2*NTENS+7
        STATEV(II)=XBACK(JJ)
        JJ=JJ+1
      END DO
      STATEV(20)=D
      STATEV(21)=ESTE
      SEQUIV=FBAR
      STATEV(22)=SEQUIV
C
  802 IF (IFLAG.EQ.1) THEN
         WRITE(*,*)
         WRITE(*,*) 'LOCAL PLASTICITY ALGORITHM DID NOT CONVERGED'
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
C     SUBROUTINE VSPRATE                                               C
C                                                                      C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE VSPRATE(PROPS,NPROPS,EMOD,ENU,SNETA,SIG0,SIGSAT,HMOD,
     1                   HRDRATE)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.0d0, TWO=2.0d0)
C
      DIMENSION PROPS(NPROPS)
C
C                              Temperature
      THETA=PROPS(1)
C                              Young's modulus
      EMOD=(PROPS(2)+PROPS(3)*THETA)*1000.0
C                              Shear modulus
      GMOD=(PROPS(4)+PROPS(5)*THETA)*1000.0
C                              Poisson's ratio
      ENU=(EMOD/(TWO*GMOD))-ONE
C                              Initial yield stress
      SIG0=PROPS(6)+PROPS(7)*THETA
C                              Coefficient of thermal expansion
      CTE=PROPS(8)
C                              Kinematic hardening parameters
      XINFI=PROPS(9)
      GAMHARD=PROPS(10)
C                              Isotropic hardening parameters
      RINFI=PROPS(11)+PROPS(12)*THETA
      CHARDI=PROPS(13)
C                              Dimensionless strain rate constant
      A=PROPS(14)
C                              Frequency factor
      D0=PROPS(15)
C                              Burger's vector magnitude
      B=PROPS(16)
C                              Boltzman's constant
      BK=PROPS(17)
C                              Average grain size
      DS=PROPS(18)
C                              Grain size eponent
      PS=props(19)
C                              Stress exponent
      SN=PROPS(20)
C                              Creep activation energy
      CQ=PROPS(21)
C                              Universal gas constant
      RGAS=PROPS(22)
C                              Avogadro's number
      AVNUM=PROPS(23)
C                              Atomic weight
      AWEIGTH=PROPS(24)
C                              Relative length scale
      PLS=PROPS(25)
      PENNUM=PROPS(26)
C
      SNETA=0.D0
      POWER=EXP(-CQ/(RGAS*THETA))
      FLUIDITY=((A*D0*EMOD*B)/(BK*THETA)*(EMOD**sn))*((B/DS)**PS)*POWER
      IF(FLUIDITY.NE.0.0) SNETA=ONE/FLUIDITY
      HMOD=GAMHARD*XINFI
      SIGSAT=RINFI
      HRDRATE=CHARDI
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C SUBROUTINE KDAMACAL                                                  C
C                                                                      C
C  Calculates the energy and damage of constitutive model used         C
C  for viscoplastic damage                                             C
C                                                                      C
C       INPUT ARGUMENTS------                                          C
C                                                                      C
C STR(4)       :Stress vector                                          C
C STRAVP(4)    :Visco-Plastic Strain increment                         C
C                                                                      C
C ESTE         :Internal and free energy terms.                        C
C D            :Damage variable .                                      C
C THTA         :Temperature                                            C
C                                                                      C
C       OUTPUT ARGUMENTS-------                                        C
C                                                                      C
C D            :Damage variable .                                      C
C ESTE         :Internal and free energy terms.                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KDAMACAL(STR,EPL,NTENS,PROPS,NPROPS,ESTE,D)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,TWO=2.D0)
C      
      DIMENSION STR(NTENS),EPL(NTENS),PROPS(NPROPS)
C
      AVNUM= PROPS(23)
      BK= PROPS(17)
      AWEIGHT= PROPS(24)
      CONST1 = AWEIGHT/(AVNUM*BK)/10.0
      THTA=PROPS(1)
C      
      DO 10 I =1,3
        ESTE=ESTE+ABS(STR(I)*EPL(I)/THTA)
   10 CONTINUE
C
C     Compute damage
C
      D =ONE-EXP(-CONST1*ESTE)
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE KUHARDNLIN                                            C
C     COMBINED ISOTROPIC/KINEMATIC HARDENING                           C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KUHARDNLIN(SYIELDI,SYIELDK,HMOD,SIG0,SIGSAT,HRDRATE,
     1                      EHARDI,EHARDK,EQPLAS)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0)
C
C     Compute yield stress and hardening modulus corresponding to EQPLAS
C     Isotropic hardenining==> Nonlinear Chaboche model
C     Kinematic hardening  ==> Linear
C
      SYIELDI=DSQRT(TWO/THREE)*(SIG0+SIGSAT*(ONE-EXP(-HRDRATE*EQPLAS)))
      EHARDI=SIGSAT*HRDRATE*(EXP(-HRDRATE*EQPLAS))
      SYIELDK=HMOD*EQPLAS
      EHARDK=HMOD
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE KSPECTRAL                                             C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,E,ENU,EK,D)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0)
C
      DIMENSION Q(6,6),DP(6,6),DC(6,6),DIAG(6,6),GDIA(6,6)
C
      CALL KCLEAR(Q,6,6)
      CALL KCLEAR(DP,6,6)
      CALL KCLEAR(DC,6,6)
      CALL KCLEAR(DIAG,6,6)
      CALL KCLEAR(GDIA,6,6)
C
      EG2=E/(ONE+ENU)
      EG=EG2/TWO
      ELAM=EG2*ENU/(ONE-TWO*ENU)
      TETA1=ONE+(TWO/THREE)*GAM_PAR*(1-D)*EK
      CBETA2=EG2+(TWO/THREE)*(1-D)*EK
      CALFA2=(EG2)+(TWO/THREE)*(1-D)*EK
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
      Q(5,5)=ONE
      Q(6,6)=ONE
C
      DP(1,1)=ONE
      DP(2,2)=ONE
      DP(3,3)=ZERO
      DP(4,4)=TWO
      DP(5,5)=ONE
      DP(6,6)=ONE
C
      DC(1,1)=EG2
      DC(2,2)=EG2
      DC(3,3)=THREE*ELAM+EG2
      DC(4,4)=EG
      DC(5,5)=EG2
      DC(6,6)=EG2
C
      DIAG(1,1)=ONE/(ONE+CBETA2*GAM_PAR)
      DIAG(2,2)=ONE/(ONE+CBETA2*GAM_PAR)
      DIAG(3,3)=ONE/TETA1
      DIAG(4,4)=ONE/(ONE+CBETA2*GAM_PAR)
      DIAG(5,5)=ONE/(ONE+CALFA2*GAM_PAR)
      DIAG(6,6)=ONE/(ONE+CALFA2*GAM_PAR)
C
      GDIA(1,1)=-(CBETA2/((ONE+CBETA2*GAM_PAR)**2))
      GDIA(2,2)=-(CBETA2/((ONE+CBETA2*GAM_PAR)**2))
      GDIA(3,3)=-((TWO/THREE)*EK)/(TETA1**2)
      GDIA(4,4)=-(CBETA2/((ONE+CBETA2*GAM_PAR)**2))
      GDIA(5,5)=-(CALFA2/((ONE+CALFA2*GAM_PAR)**2))
      GDIA(6,6)=-(CALFA2/((ONE+CALFA2*GAM_PAR)**2))
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
      DIMENSION P(6,6)
C
      CALL KCLEAR(P,6,6)
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
      P(5,5)=ONE
      P(6,6)=ONE
C
      RETURN
C
      END
C