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
C     8 NODED ISOPARAMETRIC REDUCED INTEGRATION PENALTHY METHOD        C
C     REDUCED COSSERAT COUPLE STRESS SOLID ELEMENT                     C
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
C  UEL8_VPDCO.for  8-Noded Viscoplasticity Damage Couple Stress Theory C
C  USER ELEMENT SUBROUTINE-COSSERAT ELASTIC MATERIAL BEHAVIOR          C
C                                                                      C
C  8 NODED ISOPARAMETRIC PENALTHY METHOD  COSSERAT SOLID ELEMENT       C
C  3X3 GAUSS INTEGRATION FOR STRESS AND COUPLE STRESS TERMS            C
C  2X2 GAUSS INTEGRATION FOR THE ASSYMETRIC STRESS COMPONENT           C
C                                                                      C
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
     1           THREE=3.D0)
C
C     Parameters required by UMAT.f
C
      PARAMETER (NDI=3,NSTATV=22,SSE=0.D0,SCD=0.D0,
     1           RPL=0.D0,DRPLDT=0.D0,TEMP=0.D0,DTEMP=0.D0,NSHR=1,
     2           CELENT=2.D0,LAYER=1,KSPT=1)
C
C     Parameter arrays from UEL.f
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5JPROPS(*)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),BLM(1,NDOFEL),
     1BLMT(NDOFEL,1),FRST1(NDOFEL,NDOFEL),FRST2(NDOFEL,1),
     2XX(2,NNODE),XWF(9),XPF(2,9),AUX1(NTENS,NDOFEL),
     3AUX2(1,NDOFEL),AUX3(NDOFEL,NDOFEL),AUX4(1,NDOFEL),
     4AUX5(NDOFEL,NDOFEL),STRESS(NTENS),STRAN(NTENS),
     5DSTRAN(NTENS),RELROT(1),DRELROT(1),DDBAR(1,1),STRASYM(1),
     6FRST3(NDOFEL,NDOFEL),FRST4(NDOFEL,1),XWR(4),XPR(2,4),
     7EELAS(NTENS),EPLAS(NTENS),XBACK(NTENS)
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
      CALL KCLEARV(RELROT,1)
      CALL KCLEARV(DRELROT,1)
      CALL KCLEARV(DSTRAN,NTENS)
      CALL KCLEARV(DDSDDT,NTENS)
      CALL KCLEARV(DRPLDE,NTENS)
      CALL KCLEARV(XBACK,NTENS)
      CALL KCLEAR(DROT,3,3)
      CALL KCLEAR(DFGRD0,3,3)
      CALL KCLEAR(DFGRD1,3,3)
      CALL KCLEAR(DDSDDE,NTENS,NTENS)
      CALL KCLEAR(DDBAR,1,1)
      UPREDEF(1)=ZERO
      DPRED(1)=ZERO
C
C     Clears RHS vector and Stiffness matrix
C
      CALL KCLEARV(RHS,NDOFEL)
      CALL KCLEAR(AMATRX,NDOFEL,NDOFEL)
      CALL KCLEAR(AUX1,NTENS,NDOFEL)
      CALL KCLEAR(AUX2,1,NDOFEL)
      CALL KCLEAR(AUX3,NDOFEL,NDOFEL)
      CALL KCLEAR(AUX4,2,NDOFEL)
      CALL KCLEAR(AUX5,NDOFEL,NDOFEL)
C
      NSVARS_R=8
      NSVARS_F=NSVARS-NSVARS_R
C
      CMNAME='COSSERAT'
      PLS=PROPS(25)
      PENNUM=PROPS(26)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               Starts full integration part.                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL KGPOINTS3X3(XPF,XWF)
      NGPT_F=9
C
C     Loops around all Gauss points
C
      DO NN=1,NGPT_F
C
        RII=XPF(1,NN)
        SII=XPF(2,NN)
        ALF=XWF(NN)
C
C       Compute State variable index corresponding
C       to current Gauss point and load stress,total strain
C       elastic strain,plastic strain and equivalent plastic
C       strain from state variables as defined in USER ELEMENT. 
C       Different variables are required for different constitutive models.
C
        ISVINT_F=1+(NN-1)*NSVARS_F/NGPT_F
        JJ=1
        DO II=ISVINT_F,ISVINT_F+5
          STRESS(JJ)=SVARS(II)
          XBACK(JJ )=SVARS(II+6)
          STRAN(JJ )=SVARS(II+12)
          EELAS(JJ )=SVARS(II+18)
          EPLAS(JJ )=SVARS(II+24)
          JJ=JJ+1
        END DO
        EQPLAS=SVARS(ISVINT_F+30)
        D=SVARS(ISVINT_F+31)
        ESTE=SVARS(ISVINT_F+32)
        SEQUIV=SVARS(ISVINT_F+33)
C
C       Starts STATEV array definition as required by UMAT.f
C       Different variables for different constitutive models.
C
        DO II=1,NTENS
          STATEV(II)=EELAS(II)
          STATEV(II+NTENS)=EPLAS(II)
        END DO
        STATEV(2*NTENS+1)=EQPLAS
        JJ=1
        DO II=14,2*NTENS+7
          STATEV(II)=XBACK(JJ)
          JJ=JJ+1
        END DO
        STATEV(20)=D
        STATEV(21)=ESTE
        STATEV(22)=SEQUIV
C
C       Ends STATEV array definition as required by UMAT.f
C
C       Assembles strain-displacement and relative rotation-displacement
C       matrices B and BLM=(L-M)
C
        CALL KSTDM(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,BLM,DDET,RII,
     1  SII,XBAR)
C
C       Incorporte length scale.
C
        DO K1=1,NDOFEL
          B(5,K1)=B(5,K1)*PLS
          B(6,K1)=B(6,K1)*PLS
        END DO
C
C       Computes strain increment and rel.rotation increment
C       and updates strain and rel.rotation vectors.
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
C       Assembles STIFFNESS matrix
c                                        T
C       Generalized stress component.===B M B (stress and couple stress)
C
        CALL KMMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL KMTRAN(B,NTENS,NDOFEL,BT)
        CALL KMMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
C
C       Assembles RHS vector contribution
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
        CALL KCLEAR(FRST2,NDOFEL,1)
C
C       Starts updating of state variables with updated values from UMAT.f
C
        DO II=1,NTENS
          EELAS(II)=STATEV(II)
          EPLAS(II)=STATEV(II+NTENS)
        END DO
        EQPLAS=STATEV(2*NTENS+1)
        JJ=1
        DO II=14,2*NTENS+7
          XBACK(JJ)=STATEV(II)
          JJ=JJ+1
        END DO
        D=STATEV(20)
        ESTE=STATEV(21)
        SEQUIV=STATEV(22)
C
        JJ=1
        DO II=ISVINT_F,ISVINT_F+5
          SVARS(II)=STRESS(JJ)
          SVARS(II+6 )=XBACK(JJ)
          SVARS(II+12 )=STRAN(JJ)
          SVARS(II+18 )=EELAS(JJ)
          SVARS(II+24)=EPLAS(JJ)
          JJ=JJ+1
        END DO
        SVARS(ISVINT_F+30)=EQPLAS
        SVARS(ISVINT_F+31)=D
        SVARS(ISVINT_F+32)=ESTE
        SVARS(ISVINT_F+33)=SEQUIV
C
C       Ends updating of state variables with updated values from UMAT.f
C
      END DO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               Starts reduced integration part.                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL KGPOINTS2X2(XPR,XWR)
      NGPT_R=4
C
C     Loops around all Gauss points
C
      DO NN=1,NGPT_R
C
        RII=XPR(1,NN)
        SII=XPR(2,NN)
        ALF=XWR(NN)
C
C       Compute State variable index
C       Different variables are required for different constitutive models.
C
        ISVINT_R=(ISVINT_F+33)+1+(NN-1)*NSVARS_R/NGPT_R
        JJ=1
        II=ISVINT_R
        STRASYM(JJ)=SVARS(II)
        RELROT(JJ) =SVARS(II+1)
        JJ=JJ+1
C
C       Assembles strain-displacement and relative rotation-displacement
C       matrices B and BLM=(L-M)
C
        CALL KSTDM(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,BLM,DDET,RII,
     1  SII,XBAR)
C
C       Computes relative rotation increment
C       and updates rel.rotation vector.
C
        CALL KMAVEC(BLM,1,NDOFEL,DU,DRELROT)
        CALL KUPDVEC(RELROT,1,DRELROT)
C
C       Updates asymmetric stress component.
C
        STRASYM(1)=STRASYM(1)+PENNUM*DRELROT(1)
C
C       Assembles STIFFNESS matrix
C                                                                  T_
C       Relative rotation component.(i.e constrained part) === (L-M)D(L-M)
C
        DDBAR(1,1)=PENNUM
        CALL KMMULT(DDBAR,1,1,BLM,1,NDOFEL,AUX4)
        CALL KMTRAN(BLM,1,NDOFEL,BLMT)
        CALL KMMULT(BLMT,NDOFEL,1,AUX4,1,NDOFEL,FRST3)
C
C       Assembles RHS vector contribution
C
        DO KK=1,NDOFEL
           FRST4(KK,1)=BLMT(KK,1)*STRASYM(1)
        END DO
C
C       Considers Gauss weight and Jacobian determinant representing
C       volume differential.
C
        CALL KSMULT(FRST3,NDOFEL,NDOFEL,ALF*DDET*XBAR)
        CALL KSMULT(FRST4,NDOFEL,1,-ALF*DDET*XBAR)
C
C       Updates Stiffness matrix and RHS vector
C       
        CALL KUPDMAT(AMATRX,NDOFEL,NDOFEL,FRST3)
        CALL KUPDVEC(RHS,NDOFEL,FRST4)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point 
C
        CALL KCLEAR(FRST3,NDOFEL,NDOFEL)
        CALL KCLEAR(FRST4,NDOFEL,1)
C
C       Starts updating of state variables with updated values from UMAT.f
C
        JJ=1
        II=ISVINT_R
        SVARS(II)=STRASYM(JJ)
        SVARS(II+1)=RELROT(JJ)
        JJ=JJ+1
C
      END DO
C
c
C
C     EXTRAPOLATE STRAIN TO THE NODES
C
c      RII=-1.0
c      SII=1.0
c      CALL KSTDM(JELEM,JTYPE,NNODE,NDOFEL,NTENS,COORDS,B,BD,BV,HP,
c     1           DDET,RII,SII,XBAR)
c      CALL KMAVEC(B,NTENS,NDOFEL,U,STRAN)
c      svars(109)=stress(2)
c      SVARS(110)=STRAN(2)
c
      RETURN
C      
      END
C
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
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     UMAT_VPDCO.for    Viscoplasticity Damage Couple Stress Theory    C
C     Requires a UEL.for subroutine                                    C
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C     DAMAGE MECHANICS COUPLED RATE DEPENDENT HARDENING PLASTICITY     C
C     COMBINED ISOTROPIC/KINEMATIC HARDENING                           C
C     RETURN MAPPING ALGORITHM                                         C
C                                                                      C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C     Last updated OCTOBER 13/2003                                     C
C                                                                      C
C                                                                      C
C                                                                      C
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
      CALL KCLEAR(STSRELT,1,6)
      CALL KCLEAR(SINVAR,1,1)
      CALL KMTRAN(STSREL,6,1,STSRELT)
      CALL KMMULT(P,6,6,STSREL,6,1,SDEV)
      CALL KMMULT(STSRELT,1,6,SDEV,6,1,SINVAR)
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
          CALL KMTRAN(Q,6,6,QT)
          CALL KMMULT(STSRELT,1,6,Q,6,6,AUX1)
          CALL KMMULT(QT,6,6,P,6,6,AUX2)
          CALL KMMULT(AUX2,6,6,Q,6,6,AUX3)
          CALL KMMULT(AUX3,6,6,DIAG,6,6,AUX4)
          CALL KMMULT(AUX4,6,6,QT,6,6,AUX5)
          CALL KMMULT(AUX5,6,6,STSREL,6,1,AUX6)
          CALL KMMULT(AUX1,1,6,DIAG,6,6,AUX7)
          CALL KMMULT(AUX7,1,6,AUX6,6,1,AUX8)
          CALL KMMULT(AUX1,1,6,GDIA,6,6,AUX9)
          CALL KMMULT(AUX9,1,6,AUX6,6,1,AUX10)
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
        CALL KMMULT(STSRELT,1,3,SDEV,3,1,SINVAR)
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
C     SUBROUTINE KSTDM:    GENERATES THE STRAIN-DISPLACEMENT MATRIX B  C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 9-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     BLM  =RELATIVE STRAIN DISPLACEMENT MATRIX                        C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSTDM(IDEL,NNE,NDOFEL,NTENS,XX,B,BLM,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),BLM(1,NDOFEL),
     1P(2,NNE),AUX1(2,NNE),XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL KCLEAR(B,NTENS,NDOFEL)
      CALL KCLEAR(XJ,2,2)
      CALL KCLEAR(XJI,2,2)
      CALL KCLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL KSFDER(IELT,NDOFEL,NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL KJACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL KJACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL KMMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=3*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
        B(5,II+2)=AUX1(2,I)
        B(6,II+2)=AUX1(1,I)
      END DO
C
C     Assembles (L-M) matrix to compute
C     relative strain vector.
C
      DO I=1,NNE
        II=3*(I-1)+1
        BLM(1,II)=HALF*AUX1(2,I)
        BLM(1,II+1)=-HALF*AUX1(1,I)
        BLM(1,II+2)=ONE
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
      SUBROUTINE KSFDER(IELT,NDOFEL,NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     8-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C      
      P(1,6)= HALF*SMS
      P(1,5)=-R*SM
      P(1,8)=-HALF*SMS
      P(1,7)=-R*SP
      P(1,2)= QUARTER*SM-HALF*P(1,5)-HALF*P(1,6)
      P(1,1)=-QUARTER*SM-HALF*P(1,8)-HALF*P(1,5)
      P(1,4)=-QUARTER*SP-HALF*P(1,7)-HALF*P(1,8)
      P(1,3)= QUARTER*SP-HALF*P(1,7)-HALF*P(1,6)
C
C     w.r.t.s
C
      P(2,6)=-S*RP
      P(2,5)=-HALF*RMS
      P(2,8)=-S*RM
      P(2,7)= HALF*RMS
      P(2,2)=-QUARTER*RP-HALF*P(2,5)-HALF*P(2,6)
      P(2,1)=-QUARTER*RM-HALF*P(2,8)-HALF*P(2,5)
      P(2,4)= QUARTER*RM-HALF*P(2,7)-HALF*P(2,8)
      P(2,3)= QUARTER*RP-HALF*P(2,7)-HALF*P(2,6)
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
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      DIMENSION XP(2,4),XW(4)
C
      CALL KCLEARV(XW,4)
      CALL KCLEAR(XP,2,4)
C
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
      SUBROUTINE KJACOPER(NNE,XJA,XCORD,RDER)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
C
      DIMENSION XJA(2,2),XCORD(2,NNE),RDER(2,NNE)
C
      CALL KCLEAR(XJA,2,2)
C
      DUM=ZERO
      DO K=1,2
        DO J=1,2
          DO I=1,NNE
            DUM=DUM+RDER(K,I)*XCORD(J,I)
          END DO
          XJA(K,J)=DUM
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