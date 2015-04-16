CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     UMAT_VPDCL.for   Visco Plasticity Damage Classical Theory        C
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C     PLANE STRAIN ANALYSIS                                            C
C     DAMAGE MECHANICS THERMOMECHANICAL FATIGUE                        C
C     COMBINED ISOTROPIC/KINEMATIC HARDENING VISCO-PLASYICITY          C
C     RADIAL RETURN INTEGRATIONS SCHEME                                C
C     CREATED BY JUAN GOMEZ                                            C
C     UNIVERSIDAD EAFIT                                                C
C     LABORATORIO DE MECANICA APLICADA                                 C
C     BLOQUE 14-PISO 2                                                 C
C     MEDELLIN, COLOMBIA                                               C
C                                                                      C
C     LAST UPDATED APRIL 16/2015                                       C
C                                                                      C
C     PLANE STRAIN AND 3D ANALYSIS                                     C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C                                                                      C
C     STATE VARIABLES DEFINITION                                       C
C                                                                      C
C      1-4: Elastic strain vector                                      C
C      5-8: Plastic strain vector                                      C
C        9: Equivalent plastic strain.                                 C
C    10-13: Back stress vector.                                        C
C       14: Damage Parameter                                           C
C       15: Helmhlotz free energy term                                 C
C       16: Equivalent stress                                          C
C                                                                      C
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
      DIMENSION STRESS(NTENS), STATEV(NSTATV),
     1DDSDDE(NTENS, NTENS),DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),
     2DSTRAN(NTENS),TIME(2),PREDEF(1), DPRED(1), PROPS(NPROPS),
     3COORDS(3),DROT(3, 3), DFGRD0(3, 3), DFGRD1(3, 3)
C
      DIMENSION EELAS(NTENS), EPLAS(NTENS), FLOW(NTENS)
C
C
      DIMENSION AUX1(NTENS,NTENS),AUX2(NTENS,NTENS),AUX3(NTENS,NTENS),
     1AUX4(NTENS,NTENS),P(NTENS,NTENS),SINVAR(1,1),SDEV(NTENS,1),
     2XBACK(NTENS),STSREL(NTENS),DSTRTHER(NTENS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=10)
C
C**********************************************************************
C     I S O T R O P I C  M I S E S  E L A S T O P L A S T I C I T Y
C     P L A N E  S T R A I N  A N A L Y S I S
C**********************************************************************
C
C     Recover equivalent plastic strain, elastic strains plastic
C     strains,back stress,damage and Helmholtz free energy term.
C
      EQPLAS=STATEV(1+2*NTENS)
      DO K1=1, NTENS
        EELAS(K1)=STATEV(K1)
        EPLAS(K1)=STATEV(K1+NTENS)
        FLOW(K1)=ZERO
      END DO
C
      J=1
      DO K1=10,2*NTENS+5
        XBACK(J)=STATEV(K1)
        J=J+1
      END DO
C
      D=STATEV(14)
      ESTE=STATEV(15)
      SEQUIV=STATEV(16)
C
C     Process material properties
C
      CALL VSPRATE(PROPS,NPROPS,EMOD,ENU,SNETA,SIG0,SIGSAT,HMOD,
     1             HRDRATE)
C
C     Initializes temperature
C
      IF(PROPS(1).NE.0.0)THEN
        THTA=PROPS(1)
      ELSE
        THTA=(TEMP+0.5*DTEMP)+273.15
      END IF
C
C     Thermal strain increment
C
      CALL KCLEAR(DSTRTHER,NTENS,1)
      DO K1=1,NDI
        DSTRTHER(K1)=PROPS(8)*DTEMP
      END DO
C
      SN=PROPS(20)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
C
C     Elastic stiffness(degraded by damage).
C
      DO K1=1, NDI
        DO K2=1, NDI
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
      END DO
C
      CALL KSMULT(DDSDDE,NTENS,NTENS,(ONE-D))
C
C     Calculate predictor stress and elastic strains
C
      DO K1=1, NTENS
        DO K2=1, NTENS
          STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*(DSTRAN(K1)-DSTRTHER(K1))
        END DO
        EELAS(K1)=EELAS(K1)+DSTRAN(K1)-DSTRTHER(K1)
      END DO
C
      SHYDRO=(STRESS(1)+STRESS(2)+STRESS(3))/THREE
C
C     Calculate relative trial stress
C
      CALL KPROYECTOR(P)
      CALL KMMULT(P,NTENS,NTENS,STRESS,NTENS,1,SDEV)
      DO K1=1,NTENS
        STSREL(K1)=SDEV(K1,1)-XBACK(K1)
      END DO
C
C     Calculate equivalent Von Mises stress
C
      SMISES=(STSREL(1)-STSREL(2))**2
     1+(STSREL(2)-STSREL(3))**2+(STSREL(3)-STSREL(1))**2
      DO K1=NDI+1, NTENS
        SMISES=SMISES+SIX*STSREL(K1)**2
      END DO
      SMISES=DSQRT(SMISES/TWO)
      FBAR=DSQRT(TWO/THREE)*SMISES
C
C     Get the yield stress from the specifid hardening function.
C
      CALL KUHARDNLIN(SYIEL0,SYIELDK,HMOD,SIG0,SIGSAT,HRDRATE,
     1                 EHARDI,EHARDK,EQPLAS)
C
C     Determine if actively yielding
C
      SYIELD=SYIEL0
      SYIELDK0=SYIELDK
      IF(FBAR.GT.(ONE+TOLER)*(ONE-D)*SYIEL0) THEN
C
C       Actively yielding-Perform local Newton iterations
C       to find consistncy parameter and equivalent plastic
C       strain
C
        DO K1=1, NTENS
          FLOW(K1)=STSREL(K1)/FBAR
        END DO
C
C       Starts iterations
C
        ITER=1
        GAM_PAR=ZERO+1.0d-25
        IFLAG=0
        PHIINV=0.0D0
        DETDG=0.D0
        DO
          IF(SNETA.NE.0.0) then
             PHIINV=((SNETA/DTIME)**(ONE/SN))*(GAM_PAR**(ONE/SN))
             DETDG=(ONE/SN)*((SNETA/DTIME)**(ONE/SN))
     $             *(GAM_PAR**((ONE/SN)-ONE))
          END IF
          FGAM=-(ONE-D)*SYIELD+FBAR-(ONE-D)*(EG2*GAM_PAR
     &    +DSQRT(TWO/THREE)*(SYIELDK-SYIELDK0))-PHIINV
          TETA2=ONE+((EHARDI+EHARDK)/(EG3))
          FJAC=-(ONE-D)*EG2*TETA2-DETDG
C
C         Updates
C
          GAM_PAR=GAM_PAR-FGAM/FJAC
          EQPLAS1=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR
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
C       Update back stresses, plastic strains, equivalent plastic
C       strains, stresses
C
        DO K1=1,NDI
          XBACK(K1)=XBACK(K1)+DSQRT(TWO/THREE)*(SYIELDK-SYIELDK0)
     &    *FLOW(K1)*(ONE-D)
          EPLAS(K1)=EPLAS(K1)+GAM_PAR*FLOW(K1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
          STRESS(K1)=SDEV(K1,1)+SHYDRO-(ONE-D)*EG2*GAM_PAR*FLOW(K1)
        END DO
C
        DO K1=NDI+1,NTENS
          XBACK(K1)=XBACK(K1)+DSQRT(TWO/THREE)*(SYIELDK-SYIELDK0)
     &    *FLOW(K1)*(ONE-D)
          EPLAS(K1)=EPLAS(K1)+TWO*GAM_PAR*FLOW(K1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
          STRESS(K1)=SDEV(K1,1)-(ONE-D)*EG2*GAM_PAR*FLOW(K1)
        END DO
C
        CALL KDAMACAL(STRESS,GAM_PAR,FLOW,PROPS,NPROPS,NTENS,ESTE,
     1               THTA,D)
C
        EQPLAS=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR
C
C       Formulate the Jacobian (Material tangent)
C
        CALL KCLEAR(DDSDDE,NTENS,NTENS)
C
C       Calculate effective properties 
C
        BETA1=ONE-(ONE-D)*EG2*GAM_PAR/FBAR
        C1=ONE+((EHARDK+EHARDI)/EG3)
        C2=(DETDG/EG2/(ONE-D))
        BETABAR=(ONE/(C1+C2))-((ONE-D)*EG2*GAM_PAR/FBAR)
        EFFG=EG*BETA1
        EFFG2=TWO*EFFG
        EFFLAM=ELAM+(EG2/THREE)*(ONE-BETA1)
        EFFHRD=-EG2*BETABAR
C
        DO K1=1, NDI
          DO K2=1, NDI
            DDSDDE(K2, K1)=(ONE-D)*EFFLAM
          END DO
          DDSDDE(K1, K1)=(ONE-D)*(EFFG2+EFFLAM)
        END DO
        DO K1=NDI+1, NTENS
          DDSDDE(K1, K1)=(ONE-D)*EFFG
        END DO
        DO K1=1, NTENS
          DO K2=1, NTENS
            DDSDDE(K2, K1)=DDSDDE(K2, K1)+(ONE-D)*EFFHRD
     &                     *FLOW(K2)*FLOW(K1)
          END DO
        END DO
      END IF
C
C     Store updated state variables
C
      DO K1=1, NTENS
        STATEV(      K1)=EELAS(K1)
        STATEV(NTENS+K1)=EPLAS(K1)
      END DO
      STATEV(2*NTENS+1)=EQPLAS
C
      J=1
      DO K1=10,2*NTENS+5
        STATEV(K1)=XBACK(J)
        J=J+1
      END DO
C
      STATEV(14)=D
      STATEV(15)=ESTE
      SEQUIV=FBAR
      STATEV(16)=SEQUIV
C
  802 IF (IFLAG.EQ.1) THEN
         WRITE(*,*)
         WRITE(*,*) 'LOCAL PLASTICITY ALGORITHM DID NOT CONVREGED'
         WRITE(*,*) 'AT GAUSS POINT=',NPT, 'ELEMENT=',NOEL
         WRITE(*,*) 'AFTER=',ITER,' ITERATIONS'
         WRITE(*,*) 'LAST CORRECTION=',-fgam/fjac
         CALL XIT
      END IF
C
      RETURN
C
      END
C
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE VSPRATE                                               C
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
C
C
C
      SNETA=0.D0
      POWER=EXP(-CQ/(RGAS*THETA))
      FLUIDITY=((A*D0*EMOD*B)/(BK*THETA)*(EMOD**SN))*((B/DS)**PS)*POWER
      IF(FLUIDITY.NE.0.0) SNETA=ONE/FLUIDITY
C
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
      SUBROUTINE KDAMACAL(STR,GAM_PAR,FLOW,PROPS,NPROPS,NTENS,ESTE,
     1                    THTA,D)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,TWO=2.D0)
C      
      DIMENSION STR(NTENS),STRAVPI(NTENS),PROPS(NPROPS),FLOW(NTENS)
C
      AVNUM= PROPS(23)
      BK= PROPS(17)
      AWEIGHT= PROPS(24)
      CONST1 = AWEIGHT/(AVNUM*BK)/10.0
C      
      DO 10 I =1,3
        STRAVPI(I)=GAM_PAR*FLOW(I)
        ESTE=ESTE+ABS(STR(I)*STRAVPI(I)/THTA)
   10 CONTINUE
C
      STRAVPI(4)=TWO*GAM_PAR*FLOW(I)
      ESTE=ESTE+ABS(STR(4)*STRAVPI(4)/THTA)
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
C     SUBROUTINE KUHARDKIN                                             C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KUHARDKIN(SYIELDI,SYIELDK,EHARD,EQPLAS,NVALUE,TABLE)
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
      SYIELDI=DSQRT(TWO/THREE)*(SYIEL0+EHARD*EQPLAS)
      SYIELDK=EHARD*EQPLAS
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
      SUBROUTINE UHARD(SYIELD, HARD, EQPLAS, EQPLASRT, TIME, DTIME, 
     1TEMP, DTEMP, NOEL, NPT, LAYER, KSPT, KSTEP, KINC, CMNAME, 
     2NSTATV, STATEV, NUMFIELDV, PREDEF, DPRED, NVALUE, TABLE)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION HARD(3), STATEV(NSTATV), TIME(*), PREDEF(NUMFIELDV),
     1DPRED(*)
C
      DIMENSION TABLE(2, NVALUE)
C
      PARAMETER(ZERO=0.D0, TOLER=1.D-6)
C
C    SET YIELD STRESS TO LAST VALUE OF TABLE, HARDENING TO ZERO
C
      SYIELD=TABLE(1, NVALUE)
      HARD(1)=ZERO
C
C    IF MORE THAN ONE ENTRY, SEARCH TABLE
C
      IF(NVALUE.GT.1) THEN
         DO K1=1, NVALUE-1
           EQPL1=TABLE(2, K1+1)
           IF(EQPLAS.LT.EQPL1) THEN
              EQPL0=TABLE(2, K1)
              IF(EQPL1.LE.EQPL0) THEN
                 WRITE(7, 1)
    1            FORMAT(//, 30X, '***ERROR - PLASTIC STRAIN MUST',
     1                           'BE ENTERED IN ASSCENDING ORDER')
                 CALL XIT
              END IF
C
C             CURRENT YIELD STRESS AND HARDENING
C
              IF(K1.EQ.1) THEN
                 SYIEL0=TABLE(1, 1)
                 SYIEL1=TABLE(1, 2)
                 DSYIEL=SYIEL1-SYIEL0
                 HARD(1)=DSYIEL/EQPL1
                 SYIELD=SYIEL0+EQPLAS*HARD(1)
                 GOTO 10
              END IF
C
              DEQPL=EQPL1-EQPL0
              SYIEL0=TABLE(1, K1)
              SYIEL1=TABLE(1, K1+1)
              DSYIEL=SYIEL1-SYIEL0
              HARD(1)=DSYIEL/DEQPL
              SYIELD=SYIEL0+(EQPLAS-EQPL0)*HARD(1)
              GOTO 10
           END IF
         END DO
   10    CONTINUE
      END IF
C
      HARD(2)=ZERO
      HARD(3)=ZERO
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE PROYECTOR                                             C
C     PROJECTS THE STRESS TENSOR INTO THE DEVIATORIC SPACE             C
C     FOR A PLANE STRAIN PROBLEM                                       C
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
      P(1,1)=TWO/THREE
      P(1,2)=-ONE/THREE
      P(1,3)=-ONE/THREE
      P(2,1)=-ONE/THREE
      P(2,2)=TWO/THREE
      P(2,3)=-ONE/THREE
      P(3,1)=-ONE/THREE
      P(3,2)=-ONE/THREE
      P(3,3)=TWO/THREE
      P(4,4)=ONE
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