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
C     9 NODED ISOPARAMETRIC REDUCED INTEGRATION PENALTHY METHOD        C
C     COSSERAT SOLID ELEMENT                                           C
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
C  UEL9_ECO.for   9-Noded Couple stress Elasticity                     C
C  USER ELEMENT SUBROUTINE-COSSERAT ELASTIC MATERIAL BEHAVIOR          C
C                                                                      C
C  9 NODED ISOPARAMETRIC PENALTHY METHOD  COSSERAT SOLID ELEMENT       C
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
      PARAMETER (NDI=3,NSTATV=14,SSE=0.D0,SCD=0.D0,
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
     7EELAS(NTENS),EPLAS(NTENS),dum(18),um(18),cordnp(2,9),
     8rotdisp(9)
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
      PLS=PROPS(3)
ccccc
c      PENNUM=100*pls*props(1)
ccccc
      PENNUM=PROPS(4)
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
          STRAN(JJ )=SVARS(II+6)
          JJ=JJ+1
        END DO
C
C       Starts STATEV array definition as required by UMAT.f
C       Different variables for different constitutive models.
C
C
C       Ends STATEV array definition as required by UMAT.f
C
C       Assembles strain-displacement and relative rotation-displacement
C       matrices B and BLM=(L-M)
C
        CALL KSTDM(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,BLM,DDET,RII,
     1  SII,XBAR)
C
C       Captures length scale and penalty parameter.
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
        JJ=1
        DO II=ISVINT_F,ISVINT_F+5
          SVARS(II)=STRESS(JJ)
          SVARS(II+6 )=STRAN(JJ)
          JJ=JJ+1
        END DO
C
C       Ends updating of state variables with updated values from UMAT.f
C
      END DO
cccc
C
C     Compute nodal relative rotations out of displacments
C
      call knodcord(cordnp)
C
      do kk=1,9
        rii=cordnp(1,kk)
        sii=cordnp(2,kk)
        CALL KSTDM(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,BLM,DDET,RII,
     1  SII,XBAR)
        CALL KMAVEC(BLM,1,NDOFEL,DU,DRELROT)
        rotation=DRELROT(1)
        ibandera=rotation
      end do
cccc
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
        ISVINT_R=(ISVINT_F+11)+1+(NN-1)*NSVARS_R/NGPT_R
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
C
C
ccc
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
ccc
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
C23456789012345678901234567890123456789012345678901234567890123456789012
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C     ISOTROPIC ELASTICITY                                             C
C     CANNOT BE USED FOR PLANE STRESS                                  C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C     JANUARY 26/2004                                                  C
C                                                                      C
C     LOCAL ARRAYS                                                     C
C                                                                      C
C     PROPS(1) - E                                                     C
C     PROPS(2) - NU                                                    C
C                                                                      C
C     DDSDDE() - MATERIAL JACOBIAN                                     C
C     STRESS() - UPDATED STRESS VECTOR                                 C
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
      DIMENSION DS(NTENS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=30,
     2           FOUR=4.D0)
C
C**********************************************************************
C               C O S S E R A T  E L A S T I C I T Y
C**********************************************************************
C
C     Elastic properties
C
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2),ENUMAX)
      PLS=PROPS(3)
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
C     Calculate the stress increment and
C     updates the stress vector.
C
      CALL KMAVEC(DDSDDE,NTENS,NTENS,DSTRAN,DS)
      CALL KUPDVEC(STRESS,NTENS,DS)
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
C     FOR AN 9-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     BLM  =RELATIVE STRAIN DISPLACEMENT MATRIX                        C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KSTDM(IDEL,NNE,NDOFEL,NTENS,XX,B,BLM,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0,ZERO=0.D0)
C
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),BLM(1,NDOFEL),
     1P(2,NNE),AUX1(2,NNE),XJ(2,2),XJI(2,2),SP(NNE)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL KCLEAR(B,NTENS,NDOFEL)
      CALL KCLEAR(XJ,2,2)
      CALL KCLEAR(XJI,2,2)
      CALL KCLEAR(P,2,NNE)
      CALL KCLEAR(SP,NNE,1)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL KSFDER(IELT,NDOFEL,NNE,R,S,P)
      CALL KSHAPEFUN(NNE,R,S,SP)
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
        BLM(1,II+2)=SP(I)
c        BLM(1,II+2)=zero
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
C     9-NODED ELEMENT
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
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE KSHAPEFUN:GENERATES THE SHAPE FUNCTION DERIVATIVES    C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSHAPEFUN(NNE,R,S,SF)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION SF(NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     9-NODED ELEMENT
C
      SF(9)=RMS*SMS
      SF(8)=HALF*SMS*RM-HALF*SF(9)
      SF(7)=HALF*RMS*SP-HALF*SF(9)
      SF(6)=HALF*SMS*RP-HALF*SF(9)
      SF(5)=HALF*RMS*SM-HALF*SF(9)
      SF(4)=QUARTER*RM*SP-HALF*SF(7)-HALF*SF(8)-QUARTER*SF(9)
      SF(3)=QUARTER*RP*SP-HALF*SF(7)-HALF*SF(6)-QUARTER*SF(9)
      SF(2)=QUARTER*RP*SP-HALF*SF(5)-HALF*SF(6)-QUARTER*SF(9)
      SF(1)=QUARTER*RM*SP-HALF*SF(8)-HALF*SF(5)-QUARTER*SF(9)
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
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C 3X3 GAUSS POINTS GENERATION                                          C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KNODCORD(XP)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,HALF=0.5D0)
C
      DIMENSION XP(2,9)
C
      XP(1,1)=-ONE
      XP(1,2)= ONE
      XP(1,3)= ONE
      XP(1,4)=-ONE
      XP(1,5)=ZERO
      XP(1,6)= ONE
      XP(1,7)=ZERO
      XP(1,8)=-ONE
      XP(1,9)= ZERO
C
      XP(2,1)=-ONE
      XP(2,2)=-ONE
      XP(2,3)= ONE
      XP(2,4)= ONE
      XP(2,5)=-ONE
      XP(2,6)=ZERO
      XP(2,7)= ONE
      XP(2,8)=ZERO
      XP(2,9)=ZERO
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