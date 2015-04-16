C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     GENERAL FINITE ELEMENT LIBRARY TO BE USED WITHIN ABAQUS          C
C     UNIVERSIDAD EAFIT                                                C
C     LABORATORIO DE MECANICA APLICADA                                 C
C     BLOQUE 14-PISO 2                                                 C
C     MEDELLIN, COLOMBIA                                               C
C                                                                      C
C     LAST UPDATED APRIL 16/2015                                       C
C                                                                      C
C     CURRENTLY:                                                       C
C                                                                      C
C     9 NODED ISOPARAMETRIC PURE DISPLACEMENT COSSERAT SOLID ELEMENT   C
C     USING TRANSLATIONAL DEGREES OF FREEDOM ONLY                      C
C     (NON-CONFORMING ELEMENT)                                         C
C----------------------------------------------------------------------C
C                                                                      C
C----------------------------------------------------------------------C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C       --U S E R   E L E M E N T    S U B R O U T I N E S---          C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  UEL9_PCON.for 9-Noded Plastic Cosserat Non-conforming               C
C  USER ELEMENT SUBROUTINE-COSSERAT ELASTO-PLASTIC MATERIAL BEHAVIOR   C
C  AND ISOTROPIC ELASTICITY                                            C
C                                                                      C
C  9-NODED PURE DISPLACEMENT ISOPARAMETRIC ELEMENT                     C
C  3X3 GAUSS INTEGRATION                                               C
C  CREATED BY JUAN GOMEZ                                               C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4PERIOD)
C     
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C  
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.D0,NTENS=6,TWO=2.D0,
     1           NGPTS=9, THREE=3.D0)
C
C     Parameters required by UMAT.f
C
      PARAMETER (NDI=3,NSTATV=14,SSE=0.D0,SCD=0.D0,
     1           RPL=0.D0,DRPLDT=0.D0,TEMP=0.D0,DTEMP=0.D0,NSHR=1,
     2           CELENT=2.D0,LAYER=1,KSPT=1)
C
C     Parameter arrays from UEL.f
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BD(NTENS,NDOFEL),BV(1,NDOFEL),HP(1,1),
     1BT(NDOFEL,NTENS),BDT(NDOFEL,NTENS),BVT(NDOFEL,1),HPT(1,1),
     2FRST1(NDOFEL,NDOFEL),FRST2(NDOFEL,1),XX(2,NNODE),XW(NGPTS),
     3XP(2,NGPTS),AUX1(NTENS,NDOFEL),AUX2(1,NDOFEL),AUX3(NDOFEL,NDOFEL),
     4AKPP(1,1),AKPPI(1,1),AKUP(NDOFEL,1),AKPU(1,NDOFEL),
     5AKUU(NDOFEL,NDOFEL),STRESS(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     6DSTRAND(NTENS),DSTRANV(1),EELAS(NTENS),EPLAS(NTENS),DPRESS(1)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION STATEV(NSTATV),DDSDDT(NTENS),DRPLDE(NTENS),
     1UPREDEF(1), DPRED(1),UCOORDS(3),DROT(3, 3), DFGRD0(3, 3),
     2DFGRD1(3, 3),DDSDDE(NTENS, NTENS)
C
C     Initializes parameters required by UMAT.f
C
      CALL KCLEARV(STRESS,NTENS)
      CALL KCLEARV(STRAN,NTENS)
      CALL KCLEARV(DSTRAN,NTENS)
      CALL KCLEARV(DDSDDT,NTENS)
      CALL KCLEARV(DRPLDE,NTENS)
      CALL KCLEAR(DROT,3,3)
      CALL KCLEAR(DFGRD0,3,3)
      CALL KCLEAR(DFGRD1,3,3)
      CALL KCLEAR(DDSDDE,NTENS,NTENS)
      UPREDEF(1)=ZERO
      DPRED(1)=ZERO
C
C     Clears RHS vector and Stiffness matrix
C
      CALL KCLEARV(RHS,NDOFEL)
      CALL KCLEAR(AMATRX,NDOFEL,NDOFEL)
C
      CMNAME='COSSERAT'
C
C     Generates Gauss points and weights according to JTYPE
C
      CALL KGPOINTS3X3(XP,XW)
      NGPT=9
C
C     Loops around all Gauss points
C
      DO NN=1,NGPT
C
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN)
C
C       Compute State variable index corresponding
C       to current Gauss point and load stress,total strain
C       elastic strain,plastic strain and equivalent plastic
C       strain from state variables as defined in USER ELEMENT. 
C       Different variables are required for different constitutive models.
C
        ISVINT=1+(NN-1)*NSVARS/NGPT
        JJ=1
        DO II=ISVINT,ISVINT+5
          STRESS(JJ)=SVARS(II)
          STRAN(JJ )=SVARS(II+6)
          EELAS(JJ )=SVARS(II+12)
          EPLAS(JJ )=SVARS(II+18)
          JJ=JJ+1
        END DO
        EQPLAS=SVARS(ISVINT+24)
        SMISES= SVARS(ISVINT+25)
C
C       Starts STATEV array definition as required by UMAT.f
C       Different variables for different constitutive models.
C
        DO II=1,NTENS
          STATEV(II)=EELAS(II)
          STATEV(II+NTENS)=EPLAS(II)
        END DO
        STATEV(2*NTENS+1)=EQPLAS
        STATEV(2*NTENS+2)=SMISES
C
C       Ends STATEV array definition as required by UMAT.f
C
C       Assembles B matrix
C
        CALL KSTDM(JELEM,JTYPE,NNODE,NDOFEL,NTENS,COORDS,B,DDET,RII,
     1  SII,XBAR)
C
        PLS=PROPS(3)
        DO K1=1,NDOFEL
          B(5,K1)=B(5,K1)*PLS
          B(6,K1)=B(6,K1)*PLS
        END DO
C
        CALL KMAVEC(B,NTENS,NDOFEL,DU,DSTRAN)
        CALL KUPDVEC(STRAN,NTENS,DSTRAN)
C
C       Assembles material matrix and updates state variables
C
        CALL UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1  DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,
     2  DTEMP, UPREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
     3  PROPS, NPROPS, UCOORDS, DROT, PNEWDT, CELENT, DFGRD0,
     4  DFGRD1, JELEM, NN, LAYER, KSPT, KSTEP, KINC)
C
C       Assembles STIFFNESS matrix.
C
        CALL KMMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL KMTRAN(B,NTENS,NDOFEL,BT)
        CALL KMMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
C
C       Assembles RHS vector.
C
        CALL KMAVEC(BT,NDOFEL,NTENS,STRESS,FRST2)
C
C       Considers Gauss weight and Jacobian determinant representing
C       volume differential.
C
        CALL KSMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
        CALL KSMULT(FRST2,NDOFEL,1,-ALF*DDET*XBAR)
C
C       Updates Stiffness matrix and RHS vector
C       
        CALL KUPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
        CALL KUPDVEC(RHS,NDOFEL,FRST2)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point 
C
        CALL KCLEAR(DDSDDE,NTENS,NTENS)
        CALL KCLEAR(FRST1,NDOFEL,NDOFEL)
C
C       Starts updating of state variables with updated values from UMAT.f
C
        DO II=1,NTENS
          EELAS(II)=STATEV(II)
          EPLAS(II)=STATEV(II+NTENS)
        END DO
        EQPLAS=STATEV(2*NTENS+1)
        SMISES=STATEV(2*NTENS+2)
C
        JJ=1
        DO II=ISVINT,ISVINT+5
          SVARS(II)=STRESS(JJ)
          SVARS(II+6 )=STRAN(JJ)
          SVARS(II+12 )=EELAS(JJ)
          SVARS(II+18)=EPLAS(JJ)
          JJ=JJ+1
        END DO
        SVARS(ISVINT+24)=EQPLAS
        SVARS(ISVINT+25)=SMISES
C
C       Ends updating of state variables with updated values from UMAT.f
C
      END DO
ccc
C
C     EXTRAPOLATE STRAIN TO THE NODES
C
c      RII=-1.0
c      SII=1.0
c      CALL KSTDM(JELEM,JTYPE,NNODE,NDOFEL,NTENS,COORDS,B,BD,BV,HP,
c     1           DDET,RII,SII,XBAR)
c      CALL KMAVEC(B,NTENS,NDOFEL,U,STRAN)
c      svars(235)=stran(1)
c      SVARS(236)=STRAN(5)
ccc
C
      RETURN
C      
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C--C O N S T I T U T I V E  M A T E R I A L S  S U B R O U T I N E S---C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C     ISOTROPIC HARDENING PLASTICITY                                   C
C     CANNOT BE USED FOR PLANE STRESS                                  C
C     TAKEN FROM ABAQUS SAMPLE SUBROUTINES                             C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C     JULY 14/2003                                                     C
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
C     PROPS(3..) - SYIELD AN HARDENING DATA                            C
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
     3 DFGRD0(3, 3), DFGRD1(3, 3)
C
      DIMENSION EELAS(NTENS), EPLAS(NTENS),DS(NTENS),DSTRESS(NTENS)
C
C
      DIMENSION AUX1(1,6),AUX2(6,6),AUX3(6,6),AUX4(6,6),AUX5(6,6),
     1AUX6(6,1),AUX7(1,6),AUX8(1,1),AUX9(1,6),AUX10(1,1),AUX11(6,6),
     2AUX12(6,1),AUX13(1,6),AUX14(1,1),AUX15(1,6),AUX16(6,6),AUX17(6,6),
     3STRESST(1,6),P(6,6),SINVAR(1,1),BI(6,6),STRESSUPD(6,1),
     4SDEV(6,1),EM(6,6),B(6,6),Q(6,6),QT(6,6),DP(6,6),DC(6,6),
     5DEL(6,6),BT(6,6),DIAG(6,6),GDIA(6,6),PM(6,6)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=30,
     2           FOUR=4.D0)
C
C**********************************************************************
C               C O S S E R A T  E L A S T I C I T Y
C**********************************************************************
C
C     Inititalize arrays
C
      CALL KCLEAR(AUX1,1,6)
      CALL KCLEAR(AUX2,6,6)
      CALL KCLEAR(AUX3,6,6)
      CALL KCLEAR(AUX4,6,6)
      CALL KCLEAR(AUX5,6,6)
      CALL KCLEAR(AUX6,6,1)
      CALL KCLEAR(AUX7,1,6)
      CALL KCLEAR(AUX8,1,1)
      CALL KCLEAR(AUX9,1,6)
      CALL KCLEAR(AUX10,1,1)
      CALL KCLEAR(AUX11,6,6)
      CALL KCLEAR(AUX12,6,1)
      CALL KCLEAR(AUX13,1,6)
      CALL KCLEAR(AUX14,1,1)
      CALL KCLEAR(AUX15,1,6)
      CALL KCLEAR(AUX16,6,6)
      CALL KCLEAR(AUX17,6,6)
C
      CALL KCLEAR(P,6,6)
      CALL KCLEAR(Q,6,6)
      CALL KCLEAR(DP,6,6)
      CALL KCLEAR(DC,6,6)
      CALL KCLEAR(QT,6,6)
      CALL KCLEAR(DEL,6,6)
C
      CALL KCLEAR(STRESST,1,6)
      CALL KCLEAR(SINVAR,1,1)
C
      CALL KCLEAR(STRESSUPD,6,1)
      CALL KCLEAR(AUX17,6,6)
      CALL KCLEAR(B,6,6)
      CALL KCLEAR(BT,6,6)
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
      DDSDDE(5,5)=EG2
      DDSDDE(6,6)=EG2
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
      CALL KMTRAN(STRESS,6,1,STRESST)
      CALL KMMULT(P,6,6,STRESS,6,1,SDEV)
      CALL KMMULT(STRESST,1,6,SDEV,6,1,SINVAR)
      FBAR=DSQRT(SINVAR(1,1))
C
C     Get the yield stress from the specifid hardening function.
C
      CALL KUHARD(SYIEL0,EHARD,EQPLAS,2,PROPS(4))
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
          CALL KMTRAN(Q,6,6,QT)
          CALL KMMULT(STRESST,1,6,Q,6,6,AUX1)
          CALL KMMULT(QT,6,6,P,6,6,AUX2)
          CALL KMMULT(AUX2,6,6,Q,6,6,AUX3)
          CALL KMMULT(AUX3,6,6,DIAG,6,6,AUX4)
          CALL KMMULT(AUX4,6,6,QT,6,6,AUX5)
          CALL KMMULT(AUX5,6,6,STRESS,6,1,AUX6)
          CALL KMMULT(AUX1,1,6,DIAG,6,6,AUX7)
          CALL KMMULT(AUX7,1,6,AUX6,6,1,AUX8)
          CALL KMMULT(AUX1,1,6,GDIA,6,6,AUX9)
          CALL KMMULT(AUX9,1,6,AUX6,6,1,AUX10)
          FBAR=DSQRT(AUX8(1,1))
          TETA2=ONE-(TWO/THREE)*EHARD*GAM_PAR
          FJAC=TETA2*AUX10(1,1)/FBAR-(TWO/THREE)*EHARD*FBAR
          FGAM=FBAR-SYIELD
C
C         Updates
C
          GAM_PAR=GAM_PAR-FGAM/FJAC
          EQPLAS1=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR*FBAR
          CALL KUHARD(SYIELD,EHARD,EQPLAS1,2,PROPS(4))
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
        CALL KMMULT(Q,6,6,DIAG,6,6,AUX17)
        CALL KMMULT(AUX17,6,6,QT,6,6,B)
        CALL KMMULT(B,6,6,STRESS,6,1,STRESSUPD)
C
        CALL KCLEAR(STRESS,NTENS,1)
        DO K1=1,NTENS
          STRESS(K1)=STRESSUPD(K1,1)
        END DO
C
        CALL KCLEAR(STRESST,1,6)
        CALL KCLEAR(SDEV,6,1)
        CALL KMTRAN(STRESS,6,1,STRESST)
        CALL KMMULT(P,6,6,STRESS,6,1,SDEV)
        CALL KMMULT(STRESST,1,6,SDEV,6,1,SINVAR)
        FBAR=DSQRT(SINVAR(1,1))
C
        DO K1=1,6
          EPLAS(K1)=EPLAS(K1)+GAM_PAR*SDEV(K1,1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
        END DO
C
        EQPLAS=EQPLAS1
C
C       Formulate the consistent material Jacobian (tangent)
C
        CALL KCLEAR(EM,6,6)
        CALL KMMULT(B,6,6,DDSDDE,6,6,EM)
        CALL KMMULT(EM,6,6,P,6,6,AUX11)
        CALL KMMULT(AUX11,6,6,STRESS,6,1,AUX12)
        CALL KMMULT(STRESST,1,6,P,6,6,AUX13)
        CALL KMMULT(AUX13,1,6,AUX12,6,1,AUX14)
        CALL KMTRAN(AUX12,6,1,AUX15)
        CALL KMMULT(AUX12,6,1,AUX15,1,6,AUX16)
        SCALAR1=ONE/AUX14(1,1)
        CALL KSMULT(AUX16,6,6,SCALAR1)
        TETA2=ONE-(TWO/THREE)*EHARD*GAM_PAR
        CBETA=(TWO/THREE/TETA2/AUX14(1,1))*FBAR*FBAR*EHARD
        SCALAR2=ONE/(ONE+CBETA)
        CALL KSMULT(AUX16,6,6,SCALAR2)
        CALL KCLEAR(DDSDDE,6,6)
        CALL KMATSUB(EM,6,6,AUX16,DDSDDE,0)
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
      CBETA2=EG2
c      CALFA2=TWO*EG2*(RLS**2)
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
      DIAG(3,3)=ONE
      DIAG(4,4)=ONE/(ONE+CBETA2*GAM_PAR)
      DIAG(5,5)=ONE/(ONE+CALFA2*GAM_PAR)
      DIAG(6,6)=ONE/(ONE+CALFA2*GAM_PAR)
      GDIA(1,1)=-(CBETA2/((ONE+CBETA2*GAM_PAR)**2))
      GDIA(2,2)=-(CBETA2/((ONE+CBETA2*GAM_PAR)**2))
      GDIA(3,3)=ZERO
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C         --G E N E R A L  F E M  S U B R O U T I N E S--              C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE KSTDM:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN N-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSTDM(IDEL,IELT,NNE,NDOFEL,NTENS,XX,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),P(2,NNE),PP(3,NNE),XJ1(2,2),
     2          XJ1I(2,2),AUX1(2,NNE),XJ2(3,2),SDOP(3,NNE)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL KCLEAR(B,NTENS,NDOFEL)
      CALL KCLEAR(XJ1,2,2)
      CALL KCLEAR(XJ2,3,2)
      CALL KCLEAR(P,2,NNE)
      CALL KCLEAR(PP,3,NNE)
      CALL KCLEAR(SDOP,3,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL KSFDER(IELT,NDOFEL,NNE,R,S,P,PP)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL KJACOPER(NNE,XJ1,XJ2,XX,P,PP)
      DDET=XJ1(1,1)*XJ1(2,2)-XJ1(1,2)*XJ1(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL KJACINVE(XJ1,DDET,XJ1I)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL KMMULT(XJ1I,2,2,P,2,NNE,AUX1)
C
C     Calculates second degree operator (i.e., der. w.r.t xx yy and xy)
C
      CALL KSOPMAT(XJ1,XJ2,DDET,P,PP,SDOP,NNE)
C
C     Assembles B matrix
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
        B(5,II)=-HALF*SDOP(3,I)
        B(5,II+1)=HALF*SDOP(1,I)
        B(6,II)=-HALF*SDOP(2,I)
        B(6,II+1)=HALF*SDOP(3,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE KSFDER:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSFDER(IELT,NDOFEL,NNE,R,S,P,PP)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE),PP(3,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     9-NODED PURE DISPLACEMENT COSSERAT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C
      P(1,9)=-TWO*R*SMS
      P(1,8)=-HALF*SMS-HALF*P(1,9)
      P(1,7)=-R*SP-HALF*P(1,9)
      P(1,6)= HALF*SMS-HALF*P(1,9)
      P(1,5)=-R*SM-HALF*P(1,9)
      P(1,4)=-QUARTER*SP-HALF*P(1,7)-HALF*P(1,8)-QUARTER*P(1,9)
      P(1,3)= QUARTER*SP-HALF*P(1,7)-HALF*P(1,6)-QUARTER*P(1,9)
      P(1,2)= QUARTER*SM-HALF*P(1,5)-HALF*P(1,6)-QUARTER*P(1,9)
      P(1,1)=-QUARTER*SM-HALF*P(1,8)-HALF*P(1,5)-QUARTER*P(1,9)
C
C     w.r.t.s
C
      P(2,9)=-TWO*S*RMS
      P(2,8)=-S*RM-HALF*P(2,9)
      P(2,7)= HALF*RMS-HALF*P(2,9)
      P(2,6)=-S*RP-HALF*P(2,9)
      P(2,5)=-HALF*RMS-HALF*P(2,9)
      P(2,4)= QUARTER*RM-HALF*P(2,7)-HALF*P(2,8)-QUARTER*P(2,9)
      P(2,3)= QUARTER*RP-HALF*P(2,7)-HALF*P(2,6)-QUARTER*P(2,9)
      P(2,2)=-QUARTER*RP-HALF*P(2,5)-HALF*P(2,6)-QUARTER*P(2,9)
      P(2,1)=-QUARTER*RM-HALF*P(2,8)-HALF*P(2,5)-QUARTER*P(2,9)
C
C     w.r.t.r.r
C
      PP(1,9)=-TWO*SMS
      PP(1,8)=-HALF*PP(1,9)
      PP(1,7)=-SP-HALF*PP(1,9)
      PP(1,6)=-HALF*PP(1,9)
      PP(1,5)=-SM-HALF*PP(1,9)
      PP(1,4)=-HALF*PP(1,7)-HALF*PP(1,8)-QUARTER*PP(1,9)
      PP(1,3)=-HALF*PP(1,7)-HALF*PP(1,6)-QUARTER*PP(1,9)
      PP(1,2)=-HALF*PP(1,5)-HALF*PP(1,6)-QUARTER*PP(1,9)
      PP(1,1)=-HALF*PP(1,8)-HALF*PP(1,5)-QUARTER*PP(1,9)
C
C     w.r.t.s.s
C
      PP(2,9)=-TWO*RMS
      PP(2,8)=-RM-HALF*PP(2,9)
      PP(2,7)=-HALF*PP(2,9)
      PP(2,6)=-RP-HALF*PP(2,9)
      PP(2,5)=-HALF*PP(2,9)
      PP(2,4)=-HALF*PP(2,7)-HALF*PP(2,8)-QUARTER*PP(2,9)
      PP(2,3)=-HALF*PP(2,7)-HALF*PP(2,6)-QUARTER*PP(2,9)
      PP(2,2)=-HALF*PP(2,5)-HALF*PP(2,6)-QUARTER*PP(2,9)
      PP(2,1)=-HALF*PP(2,8)-HALF*PP(2,5)-QUARTER*PP(2,9)
C
C     mixed derivatives
C
      PP(3,9)=FOUR*R*S
      PP(3,8)=S-HALF*PP(3,9)
      PP(3,7)=-R-HALF*PP(3,9)
      PP(3,6)=-S-HALF*PP(3,9)
      PP(3,5)=R-HALF*PP(3,9)
      PP(3,4)=-QUARTER-HALF*PP(3,7)-HALF*PP(3,8)-QUARTER*PP(3,9)
      PP(3,3)= QUARTER-HALF*PP(3,7)-HALF*PP(3,6)-QUARTER*PP(3,9)
      PP(3,2)=-QUARTER-HALF*PP(3,5)-HALF*PP(3,6)-QUARTER*PP(3,9)
      PP(3,1)= QUARTER-HALF*PP(3,8)-HALF*PP(3,5)-QUARTER*PP(3,9)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C    SUBROUTINE KSOPMAT                                                C
C    computes second order operator matrix                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSOPMAT(XJ1,XJ2,DDET,P,PP,SDOP,NNE)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0)
C
      DIMENSION XJ1(2,2),XJ2(3,2),P(2,NNE),PP(3,NNE),SDOP(3,NNE),C(5,3)
C
      CALL KCLEAR(SDOP,3,NNE)
      CALL KCLEAR(C,5,3)
C
      DETSQ=DDET*DDET
      DDETI=ONE/DDET
      DETSQI=ONE/DETSQ
C
C     Second order derivatives w.r.t.x
C
      DJDR=XJ2(3,1)*XJ1(1,2)+XJ2(1,2)*XJ1(2,1)-XJ2(1,1)*XJ1(2,2)-
     &XJ2(3,1)*XJ1(1,1)
      DJDR=DJDR/DETSQ
      DJDS=XJ2(2,1)*XJ1(1,2)+XJ2(3,2)*XJ1(2,1)-XJ2(3,1)*XJ1(2,2)-
     &XJ2(2,2)*XJ1(1,1)
      DJDS=DJDS*DETSQ
C
      C(1,1)=DDETI*XJ1(2,2)*XJ1(2,2)*DJDR+DETSQI*XJ1(2,2)*XJ2(3,2)-
     &DDETI*XJ1(2,2)*XJ1(1,2)*DJDS-DETSQI*XJ1(1,2)*XJ2(2,2)
      C(2,1)=DDETI*XJ1(1,2)*XJ1(1,2)*DJDS+DETSQI*XJ1(1,2)*XJ2(3,2)-
     &DDETI*XJ1(2,2)*XJ1(1,2)*DJDR-DETSQI*XJ1(2,2)*XJ2(1,2)
      C(3,1)=DETSQI*XJ1(2,2)*XJ1(2,2)
      C(4,1)=DETSQI*XJ1(1,2)*XJ1(1,2)
      C(5,1)=-(TWO/DETSQI)*XJ1(2,2)*XJ1(1,2)
C
      C(1,2)=DDETI*XJ1(2,1)*XJ1(2,1)*DJDR+DETSQI*XJ1(2,1)*XJ2(3,1)-
     &DDETI*XJ1(1,1)*XJ1(2,1)*DJDS-DETSQI*XJ1(1,1)*XJ2(2,1)
      C(2,2)=DDETI*XJ1(1,1)*XJ1(1,1)*DJDS+DETSQI*XJ1(1,1)*XJ2(3,1)-
     &DDETI*XJ1(1,1)*XJ1(2,1)*DJDR-DETSQI*XJ1(2,1)*XJ2(1,1)
      C(3,2)=DETSQI*XJ1(2,1)*XJ1(2,1)
      C(4,2)=DETSQI*XJ1(1,1)*XJ1(1,1)
      C(5,2)=-(TWO/DETSQI)*XJ1(2,1)*XJ1(1,1)
C
      C(1,3)=DDETI*XJ1(1,2)*XJ1(2,1)*DJDS+DETSQI*XJ1(1,2)*XJ2(2,1)-
     &DDETI*XJ1(2,2)*XJ1(2,1)*DJDR-DETSQI*XJ1(2,2)*XJ2(3,1)
      C(2,3)=DDETI*XJ1(1,1)*XJ1(2,2)*DJDR+DETSQI*XJ1(2,2)*XJ2(1,1)-
     &DDETI*XJ1(1,1)*XJ1(1,2)*DJDS-DETSQI*XJ1(1,2)*XJ2(3,1)
      C(3,3)=-DETSQI*XJ1(2,2)*XJ1(2,1)
      C(4,3)=-DETSQI*XJ1(1,2)*XJ1(1,1)
      C(5,3)=DETSQI*(XJ1(2,2)*XJ1(1,1)+XJ1(1,2)*XJ1(2,1))
C
      DO K1=1,3
        DO K2=1,NNE
          SDOP(K1,K2)=C(1,K1)*P(1,K2)+C(2,K1)*P(2,K2)+C(3,K1)*PP(1,K2)+
     &    C(4,K1)*PP(2,K2)+C(5,K1)*PP(3,K2)
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
C 3X3 GAUSS POINTS GENERATION                                          C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KGPOINTS3X3(XP,XW)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      DIMENSION XP(2,9),XW(9),RLGP(2,9)
C
      RLGP(1,1)=-ONE
      RLGP(1,2)=ZERO
      RLGP(1,3)= ONE
      RLGP(1,4)=-ONE
      RLGP(1,5)=ZERO
      RLGP(1,6)= ONE
      RLGP(1,7)=-ONE
      RLGP(1,8)=ZERO
      RLGP(1,9)= ONE
C
      RLGP(2,1)=-ONE
      RLGP(2,2)=-ONE
      RLGP(2,3)=-ONE
      RLGP(2,4)=ZERO
      RLGP(2,5)=ZERO
      RLGP(2,6)=ZERO
      RLGP(2,7)= ONE
      RLGP(2,8)= ONE
      RLGP(2,9)= ONE
C
      XW(1)=0.555555555555556D0**2
      XW(2)=0.555555555555556D0*0.888888888888889D0
      XW(3)=0.555555555555556D0**2
C
      XW(4)=0.555555555555556D0*0.888888888888889D0
      XW(5)=0.888888888888889D0**2
      XW(6)=0.555555555555556D0*0.888888888888889D0
C
      XW(7)=0.555555555555556D0**2
      XW(8)=0.555555555555556D0*0.888888888888889D0
      XW(9)=0.555555555555556D0**2
C
      G=DSQRT(0.60D0)
      DO I=1,2
        DO J=1,9 
          XP(I,J)=G*RLGP(I,J)
        END DO
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
C   SUBROUTINE KGPOINTS2X2                                             C
C                                                                      C
C   2X2 GAUSS POINTS GENERATION                                        C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KGPOINTS2X2(XP,XW)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,HALF=0.5d0)
C
      DIMENSION XP(2,9),XW(9)
C
      CALL KCLEARV(XW,9)
      CALL KCLEAR(XP,2,9)
C
c      DO K1=1,4
c        XW(K1)=HALF
c      END DO
C
c      XP(1,1)=-0.288675134595
c      XP(2,1)=-0.288675134595
c      XP(1,2)= 0.288675134595
c      XP(2,2)=-0.288675134595
c      XP(1,3)=-0.288675134595
c      XP(2,3)= 0.288675134595
c      XP(1,4)= 0.288675134595
c      XP(2,4)= 0.288675134595
cc
      DO K1=1,4
        XW(K1)=ONE
      END DO
C
      XP(1,1)=-0.577350269189626
      XP(2,1)=-0.577350269189626
      XP(1,2)= 0.577350269189626
      XP(2,2)=-0.577350269189626
      XP(1,3)=-0.577350269189626
      XP(2,3)= 0.577350269189626
      XP(1,4)= 0.577350269189626
      XP(2,4)= 0.577350269189626
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KJACOPER                                                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KJACOPER(NNE,XJA1,XJA2,XCORD,RDER1,RDER2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
C
      DIMENSION XJA1(2,2),XJA2(3,2),XCORD(2,NNE),RDER1(2,NNE),
     1RDER2(3,NNE)
C
      CALL KCLEAR(XJA1,2,2)
      CALL KCLEAR(XJA2,3,2)
C
      DUM=ZERO
      DO K=1,2
        DO J=1,2
          DO I=1,NNE
            DUM=DUM+RDER1(K,I)*XCORD(J,I)
          END DO
          XJA1(K,J)=DUM
          DUM=ZERO
        END DO
      END DO
C
      DUM=ZERO
      DO K=1,3
        DO J=1,2
          DO I=1,NNE
            DUM=DUM+RDER2(K,I)*XCORD(J,I)
          END DO
          XJA2(K,J)=DUM
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
C      2X2 Jacobian inverse                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KJACINVE(XJA,DD,XJAI)
C
      IMPLICIT REAL*8 (A-H,O-Z)      
C
      DIMENSION XJA(2,2),XJAI(2,2),XADJ(2,2)
C
      XADJ(1,1)=XJA(2,2)
      XADJ(1,2)=XJA(1,2)
      XADJ(2,1)=XJA(2,1)
      XADJ(2,2)=XJA(1,1)
      DO J=1,2
        DO K=1,2
          COFA=((-1)**(J+K))*XADJ(J,K)
          XJAI(J,K)=COFA/DD
        END DO
      END DO
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