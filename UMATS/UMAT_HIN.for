CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     UMAT_HIN.for  Damaged Coupled Rate Dependent Plasticity          C
C                   Hinton and Owen Integration Algorithm              C
C                                                                      C
C     COMBINED ISOTROPIC/KINEMATIC HARDENING VISCO-PLASYICITY          C
C     CREATED BY HONG TANG                                             C
C     MODIFIED BY JUAN GOMEZ                                           C
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
C     LAST UPDATED March 16/2003                                       C
C                                                                      C
C      INPUT ARGUMENTS-----                                            C
C                                                                      C
C STRESS(NTENS)       :Stress vector per ABAQUS                        C
C STRAN(NTENS)        :Strain vector per ABAQUS                        C
C DSTRAN(NTENS)       :Strain increment vector per ABAQUS              C
C STATEV(NSTATV)      :State variable vector per ABAQUS                C
C PROPS(NPROPS)       :Model properties per ABAQUS                     C
C                      will become PROP() in local notation            C
C                                                                      C
C DTIME               :Time step per ABAQUS                            C
C TEMP                :Temperature per ABAQUS                          C
C DTEMP               :Temperature increment per ABAQUS                C
C NSTATV              :Number of State Variables per ABAQUS            C
C NTENS               :Stress vector dimension                         C
C NPRPOPS             :Number of material properties                   C
C                      specified in ABAQUS input file                  C
C                                                                      C
C      OUTPUT ARGUMENTS-----                                           C
C                                                                      C
C DDSDDE(NTENS,NTENS) :Thermo-Viscoplastic constitutive                C
C                      relationship (2-D)                              C
C STRESS(NTENS)       :Stress vector per ABAQUS                        C
C STATEV(NSTATV)      :State variable vector per ABAQUS                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1                DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,
     2                DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,
     3                PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,
     4                DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C

      INCLUDE 'ABA_PARAM.INC'
C
C                                               Global variables
C

      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C                                               Local Variables
C

      PARAMETER (PI = 3.1415926535)
      DIMENSION STRAINC(4),STRAIN(4),STR(4),STRSINC(4),RST(4)
      DIMENSION CEE(4,4),DEE(4,4),XLEVP(4,4),CEVPI(4,4)
      DIMENSION RATE(4),PROP(21),G2(4)
      DIMENSION XBACK(4),DXBACK(4)
C
C                                               Initializes temperature
C

      IF(PROPS(1).NE.0.0)THEN
        THTA=PROPS(1)
      ELSE
        THTA=(TEMP+0.5*DTEMP)+273.15
      END IF
C
cccc
      if(kinc.eq.8) then
        iflag=1
      end if
cccc
C
C
C					Young's modulus
C					1000 converts to MPa
      
c
      PROP(1)=(PROPS(2)+PROPS(3)*THTA)*1000.0
C					Shear Modulus
      G=(PROPS(4)+PROPS(5)*THTA)*1000.0
C                                       Posisson's ratio.      
      PROP(2)  = PROP(1)/(2.0*G)-1.0
C                                       Termal expansion coefficient.      
      PROP(3)  = PROPS(8)
C                                       Yield strength.(Already in Mpa)       
      PROP(5)= PROPS(6)+PROPS(7)*THTA
C                                       Material constant(used in 
C                                       Isotrop. Hardening evol)      
      PROP(6)= PROPS(11)+PROPS(12)*THTA
C                                       Material constant(used in 
C                                       Isotrop. Hardening evol)         
      PROP(7)= PROPS(13)
C                                       Material constant(used in 
C				        Kinema. Hardening evol)     
      PROP(8)= PROPS(9)
C                                       Material constant(used in
C					Kinema. Hardening evol)
      PROP(9)= PROPS(10)
C
C					Flow rule parameters.
C
      PROP(10)=PROPS(14)
      PROP(11)=PROPS(15)
      PROP(12)=PROPS(16)
      PROP(13)=PROPS(18)
      PROP(15)=PROPS(17)
      PROP(17)=PROPS(21)
      PROP(18)=PROPS(22)
      PROP(19)=PROPS(20)
      PROP(20)=PROPS(19)
C
C					Damage comput.parameters
C
C                                       Avogadro's Number
      PROP(14)=PROPS(23)
C					Atomic Weight(ms/ro see Basaran and Yan 1998)
      PROP(16)=PROPS(24)                                 
C
C
C
      
C
C                                               Initializes parameters
C

      CALL KINITIAL(XKAPPA,TTIME,THTA,DTM,ROOMTEMP,STRAINC,
     &               STRAIN,STR,XBACK,STRSINC)
C
C                                               
C
C                                               Transfer stress vector
C                                               from ABAQUS to loc.var

      STR(1)=STRESS(1)
      STR(2)=STRESS(2)
      STR(4)=STRESS(3)
      STR(3)=STRESS(4)
C
      IF(STRAN(1).EQ.0.0.AND.STRAN(2).EQ.0.0.AND.STRAN(3).EQ.0.0)THEN
        ICYCLE=1
        ISTEP=1
      ELSE
        ICYCLE=2
        ISTEP=2
      END IF
C
C                                              Transfer strain increment
C                                              and strain tensor from ABA.
C                                              to local variable.


      STRAINC(1)=DSTRAN(1)
      STRAINC(2)=DSTRAN(2)
      STRAINC(3)=DSTRAN(4)
C

      STRAIN(1)=STRAN(1)
      STRAIN(2)=STRAN(2)
      STRAIN(3)=STRAN(4)
C                                               Retrieves state variables
C
C                                               AEL effective viscoplastic
C                                               strain.
      AEL=STATEV(1)
C                                               
C                                               XBACK() Back stress tensor.                                                     
      DO 5 I=1,4
        J=I+1
	XBACK(I)=STATEV(J)
    5 CONTINUE
C
C                                               D:       Damage variable
C                                               ESTE:    Increments of
C                                                        internal and 
C                                                        free energy.
C                                               DD:      Damage var. increment.
C                                               DAEL:    Effective visco-plastic
C                                                         strain rate.
C                                               R_HARDEN:Isotropic Hardening parameter.                                        
C    
      D=STATEV(6)
      ESTE=STATEV(7)
      DD=STATEV(8)
      DAEL=STATEV(9)
      R_HARDEN=STATEV(10)
      TOTIME=0.0
C
C                                              Starts iterations.
C                                              Strain not zero only for the
C                                              first iteration.
C


      DO 10 ITER  = 1,10000
        IF(TOTIME.GT.DTIME) GO TO 100
	IF (ITER.NE.1) THEN
	   STRAINC(1) = 0.0
	   STRAINC(2) = 0.0
	   STRAINC(3) = 0.0
	END IF
C
C                                              RATE:Visco-plastic strain rate.
C
        DO 70 I = 1,4
          RATE(I) = 0.0
   70   CONTINUE
C
C                                             Computes constitutive
C                                             relationship.
C   
        CALL KCONSTITU(XLEVP,CEVPI,RATE,G2,AEL,XKAPPA,RST,D,
     &  ITER,ISTEP,ICYCLE,STR,PROP,CONST1,TTIME,RSTE,THTA,DTM,
     &  ROOMTEMP,XBACK,DXBACK,DEE,IFLAG,DAEL,R_HARDEN,STRAIN)
C
C                                            Computes stresses and
C                                            performs increment.
C     
        CALL KSTRESSES(CEE,PROP,CEVPI,ESTE,D,RATE,TTIME,G2,STRAIN,AEL,
     &  STR,XBACK,TMPTIM,DSTNNOR,CTIME,TIMEH,STRAINC,STRSINC,
     &  THTA,ROOMTEMP,DTM,DXBACK,DEE,DD,DAEL)
C     
	IF(IFLAG.EQ.1) GOTO 100
	TOTIME=TOTIME+TTIME
   10 CONTINUE
C
C                                            Updates stresses
C                                             with final values.
C


  100 STRESS(1)=STR(1)
      STRESS(2)=STR(2)
      STRESS(3)=STR(4)
      STRESS(4)=STR(3)
C      
C                                            Updates state variables
C

      STATEV(1)=AEL
      DO 80 I=1,4
	J=I+1
	STATEV(J)=XBACK(I)
 80   CONTINUE
      STATEV(6)=D
      STATEV(7)=ESTE
      STATEV(8)=DD
      STATEV(9)=DAEL
      STATEV(10)=R_HARDEN
C      
C                                              Updates constitutive
C                                              relation.
C
      DDSDDE(1,1)=(1.0-D)*CEVPI(1,1)
      DDSDDE(1,2)=(1.0-D)*CEVPI(1,2)
      DDSDDE(1,3)=(1.0-D)*CEVPI(1,4)
      DDSDDE(1,4)=(1.0-D)*CEVPI(1,3)
      DDSDDE(2,1)=(1.0-D)*CEVPI(2,1)
      DDSDDE(2,2)=(1.0-D)*CEVPI(2,2)
      DDSDDE(2,3)=(1.0-D)*CEVPI(2,4)
      DDSDDE(2,4)=(1.0-D)*CEVPI(2,3)
      DDSDDE(3,1)=(1.0-D)*CEVPI(4,1)
      DDSDDE(3,2)=(1.0-D)*CEVPI(4,2)
      DDSDDE(3,3)=(1.0-D)*CEVPI(4,4)
      DDSDDE(3,4)=(1.0-D)*CEVPI(4,3)
      DDSDDE(4,1)=(1.0-D)*CEVPI(3,1)
      DDSDDE(4,2)=(1.0-D)*CEVPI(3,2)
      DDSDDE(4,3)=(1.0-D)*CEVPI(3,4)
      DDSDDE(4,4)=(1.0-D)*CEVPI(3,3)
C
      RETURN
C      
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C SUBROUTINE KINITIAL                                                  C
C Initializes material parameters                                      C
C                                                                      C
C    OUTPUT ARGUMENTS----                                              C
C                                                                      C
C STRAINC(4)    :Strain increment                                      C
C STRAIN(4)     :Strain vector(local form)                             C
C STR(4)        :Stress vector                                         C
C STRSINC(4)    :Stress increment                                      C
C XBACK(4)      :Back stress vector                                    C
C XKAPPA        :Integration constant.                                 C
C TTIME         :LocalTime step (changes with every iteration(max=.005)C
C DTM           :                                                      C
C ROOMTEMP      :Room Temperature                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KINITIAL(XKAPPA,TTIME,THTA,DTM,ROOMTEMP,STRAINC,
     &                   STRAIN,STR,XBACK,STRSINC)
C     
      IMPLICIT REAL*8 (A-H,O-Z)
C      
      DIMENSION STRAINC(4),STRAIN(4),STR(4),STRSINC(4),XBACK(4)
C

      DO 10 I=1,4
        STRAINC(I)=0.0
        STRAIN(I)=0.0
        STR(I)=0.0
        STRSINC(I)=0.0
        XBACK(I)=0.0
   10 CONTINUE
C
      XKAPPA=0.5
      DTM=0.0
      TTIME=0.005
      ROOMTEMP=293.15
C      
      RETURN
C      
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C SUBROUTINE KCONSTITU                                                 C
C                                                                      C
C Thermo-Viscoplastic constitutive relationship (2-D)                  C
C                                                                      C
C         INPUT ARGUMENTS----                                          C
C                                                                      C
C RATE(4)      :Visco-plastic strain rate.                             C
C G2(4)        :Derivative of visco-plastic strain rate w.r.t temper.  C
C PROP(21)     :Material properties                                    C
C XBACK(4)     :Back stress vector                                     C
C DXBACK(4)    :Increment in back stress vector                        C
C DEE(4,4)     :Elastic constitutive tensor(compliance form)           C 
C STRAIN(4)    :Strain vector                                          C
C STR(4)       :Stress vector                                          C
C                                                                      C
C XKAPPA       :Integration constant (0.5 for Trapezoidal rule)        C
C AEL          :Effective visco-plastic strain                         C
C D            :Damage variable                                        C
C DAEL         :Effective visco-plastic strain rate                    C
C TTIME        :Time step(local)                                       C
C R_HARDEN     :Isotropic Hardening stress component                   C
C THTA         :Temperature                                            C
C DTM          :Temperature increment                                  C
C ROOMTEMP     :Room temperature                                       C
C                                                                      C
C ITER         :Iteration number                                       C
C ISTEP        :Step Flag                                              C
C ICYCLE       :Cycle Flag                                             C
C IFLAG        :Elastic or inelastic return flag                       C
C               (0) Inelastic. (1) Elastic.                            C
C                                                                      C
C       OUTPUT ARGUMENTS---                                            C
C                                                                      C
C CEVPI(4,4)   :Constitutive Relationship Matrix                       C
C RATE(4)      :Visco-plastic strain rate.                             C
C G2(4)        :Derivative of visco-plastic strain rate w.r.t temper.  C
C XLEVP(4,4)   :                                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KCONSTITU(XLEVP,CEVPI,RATE,G2,AEL,XKAPPA,RST,D,
     &           ITER,ISTEP,ICYCLE,STR,PROP,CONST1,TTIME,RSTE,THTA,DTM,
     &           ROOMTEMP,XBACK, DXBACK,DEE,IFLAG,DAEL,R_HARDEN,STRAIN)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      
      DIMENSION  STR(4),SS_KINE(4),DFDS(4),G1(4,4),G2(4),STRAIN(4)
      DIMENSION  CEE(4,4),DEE(4,4),DEVP(4,4),PROP(21),DELTA(4),RST(3)
      DIMENSION  XLEVP(4,4),CEVPI(4,4),RATE(4),CEVPC(4,4)
      DIMENSION  XBACK(4),DXBACK(4)
C      

      IS=4
C                                                Clears storage      
      CALL KDLT2(DELTA,4)
      CALL KNULL1(XLEVP,4)
      CALL KNULL1(CEVPI,4)
      CALL KNULL1(CEVPC,4)
      CALL KNULL1(DEVP,4)
      CALL KNULL1(G1,4)
      CALL KNULL(G2,4)
      CALL KFORMC(CEE,DEE,PROP,THTA,ROOMTEMP)
      IF((ICYCLE.EQ.1.AND.ISTEP.EQ.1.AND.ITER.EQ.1)) THEN
C
C                                                Elastic behav.
C                  
        IFLAG=1
        DO 20 I = 1,4
          DO 10 J = 1,4
            CEVPI(I,J) = CEE(I,J)
	    XLEVP(I,J) = CEE(I,J)
   10     CONTINUE
   20   CONTINUE
        RETURN
      ENDIF
C
C                                               Computes Hardening parameter.
C

      CALL KHARDEN(R_HARDEN,AEL,PROP,XBACK,DXBACK,DAEL,D)
      IFLAG=0
C
C                                                Evaluates Yield Function
C      
      CALL KCHKLOAD(STR,XBACK,SS_KINE,PROP,RJ2D,STR_EFF,IFLAG,
     &             ICONV,IS, R_HARDEN, F,D)
      IF (IFLAG.EQ.0) THEN
C
C                                               Yield function derivative
C                                               w.r.t  J
C                                                       2D
C      
         CALL KDFDJ(RJ2D,FRJ2D)
C
C                                               Yield function derivative
C                                               w.r.t stress tensor.
C         
         CALL KFORMDFDS(DFDS,SS_KINE,FRJ2D)
C
C                                               Derivateives of visco-plastic
C                                               strain rate (G1 and G2)
C         
         CALL KFORMG(FRJ2D,STR_EFF,DFDS,SS_KINE,IS,G1,G2,THTA,PROP)
C
C                                               Visco-plastic strain rate.
C                                                         
	 CALL KVPRATE(RATE,DFDS,IS,THTA,PROP,F)
C
C                                                Effective viscoplastic strain
C                                                and effec. vplas. strain rate.
C

         EPBAR=SQRT((2.0*(RATE(1)*RATE(1)+RATE(2)*RATE(2)+
     $         RATE(4)*RATE(4))+RATE(3)*RATE(3))/3.0)
C     
         TSBAR=SQRT((2.0*(STRAIN(1)*STRAIN(1)+STRAIN(2)*STRAIN(2)+
     $         STRAIN(4)*STRAIN(4))+STRAIN(3)*STRAIN(3))/3.0)
C
C                                                Computes new local time step.
C     
         DOLD=TTIME
         TAU=0.1
         TTIME=TAU*TSBAR/EPBAR
         IF (ITER.GT.1) THEN
           DTMAX=(3.0/2.0)*DOLD
           IF(TTIME.GT.DTMAX) TTIME=DTMAX
         END IF
C
         IF(TTIME.GT.0.005) TTIME=0.005
C
C                                                Updates constitutive tensors.
C
         DO 40 K=1,IS
           DO 30 L=1,IS
             DEVP(K,L)=DEE(K,L)+TTIME*XKAPPA*G1(K,L)
   30      CONTINUE
   40    CONTINUE
         CALL KINVER(DEVP,IS,IS,CEVPI)
      ELSE
         DO 60 I = 1,4
	   DO 50 J = 1,4
	     CEVPI(I,J) = CEE(I,J)
  50	   CONTINUE
  60     CONTINUE
      ENDIF
C
      RETURN
C      

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KFORMC						       C
C								       C
C  Computes the elastic constitutive matrix in stiffness and           C
C  compliance form.                                                    C
C								       C
C       INPUT ARGUMENTS------                                          C
C								       C
C PROP          :Material properties array                             C
C THTA          :Temperature                                           C
C ROOMTEMP      :Room temperature                                      C
C								       C
C     OUTPUT ARGUMENTS-----   					       C   
C                        					       C   
C CEE           :Elastic constitutive matrix in stiffness form         C
C DEE           :Elastic constitutive matrix in compliance form        C
C                        					       C   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KFORMC(CEE,DEE,PROP,THTA,ROOMTEMP)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      
      DIMENSION CEE(4,4),DEE(4,4),PROP(21)
      DATA ONE,TWO/1.E0,2.E0/
C      
      E=PROP(1)
      V=PROP(2)
      TERM1 = ONE-V
      TERM2 = ONE-V**2
      TERM3 = V/TERM1
      TERM5 = TWO*(ONE+V)
      TERM6=ONE/E
      TERM7=TERM2/E
      IS=4
C     

      CEE(1,1) = ONE
      CEE(1,2) = TERM3
      CEE(1,3) = 0.0
      CEE(1,4) = TERM3
      CEE(2,1) = TERM3
      CEE(2,2) = ONE
      CEE(2,3) = 0.0
      CEE(2,4) = TERM3
      CEE(3,1) = 0.0
      CEE(3,2) = 0.0
      CEE(3,3) = (ONE-TWO*V)/(TWO*TERM1)
      CEE(3,4) = 0.0
      CEE(4,1) = TERM3
      CEE(4,2) = TERM3
      CEE(4,3) = 0.0
      CEE(4,4) = ONE
      DEE(1,1) = ONE 
      DEE(1,2) = -V
      DEE(1,3) = 0.0
      DEE(1,4) = -V
      DEE(2,1) = -V 
      DEE(2,2) = ONE 
      DEE(2,3) = 0.0
      DEE(2,4) = -V
      DEE(3,1) = 0.0
      DEE(3,2) = 0.0
      DEE(3,3) = 2.0+2.0*V
      DEE(3,4) = 0.0
      DEE(4,1) = -V
      DEE(4,2) = -V
      DEE(4,3) = 0.0
      DEE(4,4) = ONE
C

      CONST2=E*TERM1/((ONE+V)*(ONE-TWO*V))
C

      DO 30 I=1,IS
        DO 20 J=1,IS
          CEE(I,J)=CONST2*CEE(I,J)
          DEE(I,J)=TERM6*DEE(I,J)
  20    CONTINUE
  30  CONTINUE
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C SUBROUTINE KSTRESSES                                                 C
C								       C
C Computes Element Stresses                                            C
C								       C
C      INPUT ARGUMENTS----					       C
C							               C
C CEE(4,4)        :Elastic constitutive relation matrix(stiffness form)C
C PROP(21)        :Material properties array                           C
C CEVPI(4,4)      :Thermo-viscoplastic constitutive relation matrix    C
C RATE(4)         :Visco-plastic strain rate                           C
C G2(4)           :Derivative of VPSR with respect to temperature      C
C STRAIN(4)       :Strain vector                                       C
C STRAINC(4)      :Strain increment vector                             C
C STRSINC(4)      :Stress increment                                    C
C STR(4)          :Stress vector                                       C
C XBACK(4)        :Back stress vector                                  C
C DXBACK(4)       :Increment of Back stress vector                     C
C DEE(4,4)        :Elastic constitutive relation matrix(complian. form)C
C                                                                      C
C                                                                      C
C ESTE            :Internal and free energy terms.                     C 
C D               :Damage variable.                                    C
C DD              :Increment of damage variable.                       C
C DAEL            :Effective visco-plastic strain rate                 C
C TTIME           :Time step                                           C
C AEL             :Effective visco-plastic strain.                     C
C THTA            :Temperature                                         C
C ROOMTEMP        :Room temperature                                    C
C DTM             :Local step size.                                    C
C                                                                      C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSTRESSES(CEE,PROP,CEVPI,ESTE,D,RATE,TTIME,G2,STRAIN,
     &       AEL,STR,XBACK,TMPTIM,DSTNNOR,CTIME,TIMEH,STRAINC,STRSINC,
     &       THTA,ROOMTEMP,DTM,DXBACK,DEE,DD,DAEL)
C     

      IMPLICIT REAL*8 (A-H,O-Z)
C      

      DIMENSION CEVPI(4,4),PROP(21),RATE(4),G2(4),STR(4),STRAVP(4)
      DIMENSION STRAINC(4),STRAIN(4),STRSINC(4),XBACK(4),DXBACK(4)
      DIMENSION CEE(4,4),DEE(4,4)
      DATA THREE /3.E0/
C

      IS=4
      XKAPA  = 0.5
      ALPHAT = PROP(3)
C
C                                    Computes stress vector increment
C

      CALL KSTRINC(CEVPI,STRAINC,RATE,G2,TTIME,STR,XKAPA,DTM,
     &             STRSINC,ALPHAT,D,DD)
C
C                                    Computes plastic strain increment
C     
      CALL KALSTRN(CEE,STRAINC,STRSINC,DTM,ALPHAT,STRAVP,PROP,THTA,
     &             ROOMTEMP,IS,DEE)
C
C                                    Computes effective visco-plastic 
C                                    strain increment. 
C     
      DAEL=(STRAVP(1)**2+STRAVP(2)**2+STRAVP(4)**2 +
     &     2*(STRAVP(3)/2)**2)*2./3.
C     
      DAEL = SQRT(DAEL)
C
C                                     Computes Effective visco-plastic strain.
C      
      AEL = AEL + DAEL
C
C                                    Updates Backstress vector.
C
      DO 15 I=1,4
        DXBACK(I)=(2.0/3.0)*PROP(9)*PROP(8)*STRAVP(I)
     &            -PROP(9)*XBACK(I)*DAEL
        XBACK(I)=XBACK(I)+DXBACK(I)*(1.0-D)
   15 CONTINUE

c      CALL KHARDEN(R_HARDEN,AEL,PROP,DAEL,D)
                  
C
C                                     Computes damage variable.
C
      CALL KDAMACAL(STR,STRAVP,ESTE,D,PROP,THTA,IS,DD)
C

      RETURN
C      

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C
CSUBROUTINE KSTRINC						       C    
C      Computes the stress vector increment                            C
C      and performs an increment                                       C
C								       C
C      INPUT ARGUMENTS----					       C
C							               C
C CEVPI(4,4)      :Thermo-viscoplastic constitutive relation matrix    C
C RATE(4)         :Visco-plastic strain rate.                          C 
C STRSINC(4)      :Stress increment.                                   C
C G2(4)           :Deriv. of viscoplastic strain rate w.r.t. temperat. C
C STR(4)          :Stress tensor                                       C
C TTIME           :Time step(local)                                    C
C 								       C
C XKAPPA          :Integration constant.                               C
C DTM             :Temperature increment                               C
C ALPHAT          :Termal expansion coefficient                        C
C D               :Damage variable                                     C
C DD              :Damage variable increment                           C
C                                                                      C
C       OUTPUT ARGUMENTS-----                                          C
C                                                                      C 
C STRAINC(4)      :Strain increment vector                             C
C STR(4)          :Stress tensor                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSTRINC(CEVPI,STRAINC,RATE,G2,TTIME,STR,
     &                  XKAPA,DTM,STRSINC,ALPHAT,D,DD)
C     
      IMPLICIT REAL*8 (A-H,O-Z)
C      
      DIMENSION CEVPI(4,4),STRAINC(4),RATE(4),G2(4),STR(4),
     &          STRSINC(4), DELTA(4)
C

      EPSINCX =STRAINC(1)
      EPSINCY =STRAINC(2)
      EPSINCS =STRAINC(3)
C
      DAMAGE=1.0-D
C      


      CALL KDLT2(DELTA,IS)
C      

      STRSINC(1) =CEVPI(1,1)*(EPSINCX-
     &          TTIME*RATE(1)-TTIME*XKAPA*DTM*G2(1)-
     &          ALPHAT*DTM*DELTA(1))+ CEVPI(1,2)*(EPSINCY-TTIME*RATE(2)
     &          -TTIME*XKAPA*DTM*G2(2)-ALPHAT*DTM*DELTA(2))
     &		+ CEVPI(1,3)*(EPSINCS-TTIME*RATE(3)-
     &		TTIME*XKAPA*DTM*G2(3)-ALPHAT*DTM*DELTA(3))
C     

      STRSINC(2) =CEVPI(2,1)*(EPSINCX-
     &          TTIME*RATE(1)-TTIME*XKAPA*DTM*G2(1)-
     &          ALPHAT*DTM*DELTA(1))+ CEVPI(2,2)*(EPSINCY-TTIME*RATE(2)
     &          -TTIME*XKAPA*DTM*G2(2)-ALPHAT*DTM*DELTA(2))
     &          + CEVPI(2,3)*(EPSINCS-TTIME*RATE(3)-
     &          TTIME*XKAPA*DTM*G2(3)-ALPHAT*DTM*DELTA(3))
C     

      STRSINC(3) =CEVPI(3,1)*(EPSINCX-
     &          TTIME*RATE(1)-TTIME*XKAPA*DTM*G2(1)-
     &          ALPHAT*DTM*DELTA(1))+ CEVPI(3,2)*(EPSINCY-TTIME*RATE(2)
     &          -TTIME*XKAPA*DTM*G2(2)-ALPHAT*DTM*DELTA(2))
     &          + CEVPI(3,3)*(EPSINCS-TTIME*RATE(3)-
     &          TTIME*XKAPA*DTM*G2(3)-ALPHAT*DTM*DELTA(3))
C     

      STRSINC(4) =CEVPI(4,1)*(EPSINCX-
     &          TTIME*RATE(1)-TTIME*XKAPA*DTM*G2(1)-
     &          ALPHAT*DTM*DELTA(1))+ CEVPI(4,2)*(EPSINCY-TTIME*RATE(2)
     &          -TTIME*XKAPA*DTM*G2(2)-ALPHAT*DTM*DELTA(2))
     &          + CEVPI(4,3)*(EPSINCS-TTIME*RATE(3)-
     &          TTIME*XKAPA*DTM*G2(3)-ALPHAT*DTM*DELTA(3))
C
       STR(1) = STR(1) + DAMAGE*STRSINC(1) - DD*STR(1)
       STR(2) = STR(2) + DAMAGE*STRSINC(2) - DD*STR(2)
       STR(3) = STR(3) + DAMAGE*STRSINC(3) - DD*STR(3)
       STR(4) = STR(4) + DAMAGE*STRSINC(4) - DD*STR(4)
C

      RETURN
C      

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C SUBROUTINE KALSTRN                                                   C
C Copmputes visco plastic strain increment.(2)                         C
C                                                                      C
C       INPUT ARGUMENTS------                                          C
C                                                                      C
C STRAINC(IS)       :Strain increment vecor                            C
C STRSINC(IS)       :Stress Increment vector                           C
C PROP(21)          :Material properties                               C
C DEE(4,4)          :Elastic constitutive tensor.                      C
C 								       C
C DTM               :Temperature increment.                            C
C ALPHAT            :Termal expansion coefficient.                     C
C THTA              :Temperature.                                      C
C ROOMTEMP          :Room temperature                                  C
C IS                :Identity matrix dimension.                        C
C                                                                      C
C       OUTPUT ARGUMENTS-------                                        C
C                                                                      C
C STRAVP(IS)        :Visco-plastic strain increment vector.            C
C                                                                      C         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C*


       SUBROUTINE KALSTRN(CEE,STRAINC,STRSINC,DTM,ALPHAT,STRAVP,
     &                   PROP,THTA,ROOMTEMP,IS,DEE)
C
       IMPLICIT REAL*8 (A-H,O-Z)
C        

       DIMENSION  STRAINC(IS),STRSINC(IS),DELTA(4)
       DIMENSION  STRAVP(IS),CEE(4,4),DEE(4,4),PROP(21)
C
C                                       Generates Kroneker's delta
C

       CALL KDLT2(DELTA,IS)


       STRAVP(1)= STRAINC(1)-DEE(1,1)*STRSINC(1)
     & -DEE(1,2)*STRSINC(2)- DEE(1,3)*STRSINC(3)-DEE(1,4)*STRSINC(4)
     & -DTM*ALPHAT*DELTA(1)
C     

       STRAVP(2)= STRAINC(2)-DEE(2,1)*STRSINC(1)
     & -DEE(2,2)*STRSINC(2)- DEE(2,3)*STRSINC(3)-DEE(2,4)*STRSINC(4)
     & -DTM*ALPHAT*DELTA(2)
C    

       STRAVP(3)= STRAINC(3) -DEE(3,1)*STRSINC(1)
     & -DEE(3,2)*STRSINC(2)- DEE(3,3)*STRSINC(3)-DEE(3,4)*STRSINC(4)
C     

       STRAVP(4)= 0.0
C            

       RETURN
C       

       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KDAMACAL						       C
C								       C
C  Calculates the energy and damage of constitutive model used         C
C  for viscoplastic damage.					       C
C								       C
C       INPUT ARGUMENTS------                                          C
C								       C
C STR(4)       :Stress vector                                          C
C STRAVP(4)    :Visco-Plastic Strain increment                         C
C PROP(21)     :Material properties                                    C
C								       C
C ESTE         :Internal and free energy terms.                        C
C D            :Damage variable .       			       C
C THTA         :Temperature					       C
C IS           :Identity matrix dimension.      		       C
C								       C
C       OUTPUT ARGUMENTS-------                                        C
C								       C
C DD           :Damage variable increment                              C
C D            :Damage variable .       			       C
C ESTE         :Internal and free energy terms.                        C
C								       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KDAMACAL(STR,STRAVP,ESTE,D,PROP,THTA,IS,DD)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      

      DIMENSION STR(4),STRAVP(4),PROP(21)
C      

      R0    = PROP(14)
      BOLZ_K= PROP(15)
      AV_MOL= PROP(16)
      CONST1 = AV_MOL/(R0*BOLZ_K)/10.0
C      

      DO 10 I =1,IS
        ESTE=ESTE+ABS(STR(I)*STRAVP(I)/THTA)
   10 CONTINUE
C   
      ESTE=ESTE+ABS(STR(3)*STRAVP(3)/THTA)
      D0=D
      D = 1.0 - EXP(-CONST1*ESTE)
      DD=D-D0
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KFORMG						       C
C								       C
C   THIS SUBROUTINE CALCULTES new G1 and G2 MATRICES                   C
C   CONTAINING DERIVATIVES OF VISCOPLASTIC STRAIN RATE W.R.T.STRESS  . C
C   AND TEMP RESPECTIVELY                                              C
C                                                                      C
C       INPUT ARGUMENTS------                                          C
C								       C
C PROP(21)       :Material properties array			       C
C SS_KINE(4)     :Relative stress deviator                             C
C DFDFS(4)       :Yield funcn. deriv. w.r.t stresses                   C
C								       C	
C THTA           :Temperature					       C
C								       C
C        OUTPUT ARGUMENTS-----					       C
C								       C
C G1(4,4)        :Derivative of VPSR w.r.t stress                      C
C G2(4,4)        :Derivative of VPSR w.r.t temper                      C
C								       C			
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KFORMG(FRJ2D,STR_EFF,DFDS,SS_KINE,IS,G1,G2,THTA,PROP)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C

      DIMENSION G1(4,4),G2(4),PROP(21),SS_SIGM(4,4),DFDS(4),SS_KINE(4)
C

      CONST_A=PROP(10)
      D0=PROP(11)
      b=PROP(12)
      boltzk=PROP(15)
      d=PROP(13)
      CONST_p=PROP(20)
      CONST_n=PROP(19)
      Q=PROP(17)
      R=PROP(18)
      CONST1=CONST_A*D0*PROP(1)*b*(b/d)**CONST_p*EXP(-Q/R/THTA)/boltzk/
     &       THTA/(PROP(1))**CONST_n
C     

      CONST4=CONST_A*D0*PROP(1)*b/boltzk*(b/d)**CONST_p*
     &(STR_EFF/PROP(1))**CONST_n/THTA**2*EXP(-Q/R/THTA)*(-1+Q/R/THTA)
C     

      SS_SIGM(1,1)=2./3.
      SS_SIGM(1,2)=-1./3.
      SS_SIGM(1,3)=0.
      SS_SIGM(1,4)=-1./3.
      SS_SIGM(2,1)=-1./3.
      SS_SIGM(2,2)=2./3.
      SS_SIGM(2,3)=0.
      SS_SIGM(2,4)=-1./3.
      SS_SIGM(3,1)=0.
      SS_SIGM(3,2)=0.
      SS_SIGM(3,3)=2.
      SS_SIGM(3,4)=0.
      SS_SIGM(4,1)=-1./3.
      SS_SIGM(4,2)=-1./3.
      SS_SIGM(4,3)=0.
      SS_SIGM(4,4)=2./3.
      DO 10 I=1,IS
	DO 10 J=1,IS
	   G1(I,J) = 3./2.*CONST1*((CONST_n-1)*STR_EFF**(CONST_n-3)
     $              *SS_KINE(I)*SS_KINE(J)*3./2.+STR_EFF**(CONST_n-1)
     $              *SS_SIGM(I,J))
   10 CONTINUE
C   

      DO 20 I=1,IS

        G2(I) = CONST4*DFDS(I)
   20 CONTINUE
C      

      RETURN
C   

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KVPRATE						       C
C FLOW RULE                                                            C
C								       C
C Calculates viscoplastic strain increment (1)                         C
C        							       C
C								       C
C        INPUT ARGUMETNS-----					       C
C								       C
C DFDS(4)       :Yield funcn. deriv. w.r.t stresses                    C
C PROP(21)      :Material properties array			       C
C								       C
C THTA          :Temperature					       C
C F             :Yield function.         			       C
C IS            :Identity matrix dimension.      		       C       
C								       C
C        OUTPUT ARGUMENTS----					       C
C								       C
C RATE(4)       :Viscoplastic strain rate vector		       C
C								       C	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KVPRATE(RATE,DFDS,IS,THTA,PROP,F)
C      

      IMPLICIT REAL*8 (A-H,O-Z)
C      

      DIMENSION DFDS(4),RATE(4),PROP(21)
C
      CONST_A=PROP(10)
      D0=PROP(11)
      b=PROP(12)
      boltzk=PROP(15)
      d=PROP(13)
      CONST_p=PROP(20)
      CONST_n=PROP(19)
      Q=PROP(17)
      R=PROP(18)
      CONST=CONST_A*D0*PROP(1)*b/boltzk/THTA*(b/d)**CONST_p*
     &     (F/PROP(1))**CONST_n*exp(-Q/R/THTA)
      DO 30 I=1,IS
        RATE(I) = CONST*DFDS(I)
   30 CONTINUE
C     
      RETURN
C     	

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KFORMDFDS						       C	  
C								       C
C Calculates derivative of the yield function F w.r.t stresses         C
C								       C
C       INPUT ARGUMENTS------                                          C
C								       C
C SS_KINE(4)    :Relative stress deviator                              C
C FRJ2D         :Second invariant of the relative stress tensor        C
C								       C
C       OUTPUT ARGUMENTS-----					       C
C								       C
C DFDS          :Yield function derivative w.r.t stress tensor.        C
C								       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KFORMDFDS(DFDS,SS_KINE,FRJ2D)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION SS_KINE(4),DFDS(4)
C

      DFDS(1) = FRJ2D*SS_KINE(1)
      DFDS(2) = FRJ2D*SS_KINE(2)
      DFDS(3) = FRJ2D*SS_KINE(3)
      DFDS(4) = FRJ2D*SS_KINE(4)
C

      RETURN
C      
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KHARDEN 						       C
C								       C
C  Computes the hardening parameter				       C
C								       C
C       INPUT ARGUMENTS------                                          C
C PROP(21)      :Material properties array                             C
C XBACK(4)      :Back stress vector                                    C
C DXBACK(4)     :Increment in back stress vector                       C
C								       C
C AEL           :Effective visco-plastic strain.		       C
C DAEL          :Increment of Effective visco-plastic strain.	       C
C D             :Damage variable.				       C 
C								       C
C     OUTPUT ARGUMENTS-----   					       C   
C								       C
C R_HARDEN      :Hardening parameter				       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KHARDEN(R_HARDEN,AEL,PROP,XBACK,DXBACK,DAEL,D)

C
      IMPLICIT REAL*8 (A-H,O-Z)
C      

      DIMENSION PROP(21),XBACK(4),DXBACK(4)

C

      Q_HARDEN = PROP(6)

      B_HARDEN = PROP(7)
      R_HARDEN=R_HARDEN+B_HARDEN*(Q_HARDEN-R_HARDEN)*DAEL*(1.0-D)

C

      RETURN
C      

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KDFDJ						       C
C								       C
C Computes derivative of the yield function F w.r.t J2D		       C
C								       C
C       INPUT ARGUMENTS------                                          C
C								       C
C RJ2D      :Second invariant of the relative stress tensor            C
C								       C
C       OUTPUT ARGUMENTS-----					       C
C								       C
C FRJ2D     :Yield function derivative w.r.t RJ2D                      C
C								       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KDFDJ(RJ2D,FRJ2D)
C
      IMPLICIT REAL*8 (A-H,O-Z)

C
      FRJ2D = 0.5*SQRT(3./RJ2D)
C      
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KCHKLOAD						       C
C								       C
C Evaluates yield function and sets value of IFLAG depending on the    C
C result.							       C
C								       C
C       INPUT ARGUMENTS------                                          C
C								       C
C STR(4)      :Stress vector					       C
C XBACK(4)    :Back stress vector				       C
C PROP(21)    :Material properties array			       C
C								       C
C R_HARDEN    :Isotropic Hardening stress        		       C
C D           :Damage variable					       C
C								       C
c IS          :Identity matrix diemension.			       C
C								       C
C     OUTPUT ARGUMENTS-----   					       C
C								       C
C F           :Yield function value				       C
C STR_EFF     :Effective equivalent stress   			       C
C SS_KINE(4)  :Deviator of relative stress vector.       	       C
C RJ2D        :Second invariant of relative stress tensor              C
C IFLAG       :Viscoplasticity Flag				       C 
C      
C							       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C

      SUBROUTINE KCHKLOAD(STR,XBACK,SS_KINE,PROP,RJ2D,STR_EFF,IFLAG,
     &                   ICONV,IS, R_HARDEN, F,D)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      

      DIMENSION STR(4),XBACK(4),SS_KINE(4),PROP(21)
C      
      DATA TOLCL/1.0e-7/
C

      RK0 = PROP(5)*(1.0-D)
C
C				Computes effective equivalent stress
C

      CALL KINVAR(STR,XBACK,SS_KINE,IS,RJ2D,PROP,STR_EFF)
      F=STR_EFF-RK0-R_HARDEN
C

      IF(F.LT.0.0) IFLAG=1
      IF(F.GT.TOLCL) ICOUN=1
C

      RETURN
C      

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C								       C
C SUBROUTINE KINVAR						       C
C								       C
C Computes 2nd invariant of the relative stress vector and             C  
C effective equivalent stresses                                        C
C								       C
C       INPUT ARGUMENTS------                                          C
C								       C
C STR(4)      :Stress vector					       C
C XBACK(4)    :Back stress vector				       C	
C PROP(21)    :Material properties array			       C
C							               C
C IS          : UNUSED ARGUMENT					       C	
C								       C
C     OUTPUT ARGUMENTS-----   					       C
C								       C
C SS_KINE     :Stress deviator from relative stress vector             C
C STR_EFF     :Effective equivalent stress			       C
C RJ2D        :Second invariant of the relative stress vector          C
C								       C	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KINVAR(STR,XBACK,SS_KINE,IS,RJ2D,PROP,STR_EFF)
C

      IMPLICIT REAL*8 (A-H,O-Z)
C      

      DIMENSION  STR(4),XBACK(4),SS_KINE(4),PROP(21)
C

      V = PROP(2)
C 
C			  Computes relative stress vector
C
      STR_KINE1 = STR(1)-XBACK(1)
      STR_KINE2 = STR(2)-XBACK(2)
      STR_KINE3 = STR(3)-XBACK(3)
      STR_KINE4 = STR(4)-XBACK(4)
C
C                          Computes relative mean stress
C                          and relative stress deviator.
C	

      PKINE = (STR_KINE1+STR_KINE2+STR_KINE4)/3.0
      SS_KINE(1) = STR_KINE1-PKINE
      SS_KINE(2) = STR_KINE2-PKINE
      SS_KINE(4) = STR_KINE4-PKINE
      SS_KINE(3) = STR_KINE3
C
C                          Computes 2nd invariant of relative
C                          stresses.
C	

      RJ2D=0.5*(SS_KINE(1)**2+SS_KINE(2)**2+2*SS_KINE(3)**2
     &         +SS_KINE(4)**2)
C
C                          Computes effective equivalent stress
C     



      STR_EFF=SQRT(3.*RJ2D)
C      	

      RETURN
C    
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C                                                                      C 
C                  U T I L I T I E S   B L O C K                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

        SUBROUTINE KDLT2(DELTA,IS)
C***********************************************************************
C
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTA(4)
        DATA ZERO,ONE/0.E0,1.E0/
        DELTA(1)=ONE
        DELTA(2)=ONE
        DELTA(3)=ZERO
        DELTA(4)=ZERO
        IF (IS.EQ.4) THEN
           DELTA(4)=ONE
        ENDIF

        RETURN
        END

        SUBROUTINE KINVER(A,N,IS,Y)
C***********************************************************************
C
C
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(4,4),Y(4,4),INDX(4)

      CALL KNULL1(Y,4)
      DO 10 I=1,4
   10 Y(I,I)=1.
C
      CALL KLUDCMP(A,N,IS,INDX,D)
c
      DO 20 J=1,N
      CALL KLUBKSB(A,N,IS,INDX,Y(1,J))
   20 CONTINUE
      RETURN
      END
C
C

        SUBROUTINE KLUBKSB(A,N,NP,INDX,B)
C***********************************************************************
C
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(4,4),INDX(4),B(4)
        II=0
        DO 10 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF(II.NE.0)THEN
        DO 20 J=II,I-1
        SUM=SUM-A(I,J)*B(J)
20      CONTINUE
        ELSE IF(SUM.NE.0.0D0) THEN
        II=I
        ENDIF
        B(I)=SUM
10      CONTINUE
        DO 30 I=N,1,-1
        SUM=B(I)
        DO 40 J=I+1,N
        SUM=SUM-A(I,J)*B(J)
40      CONTINUE
        B(I)=SUM/A(I,I)
30      CONTINUE
        RETURN
        END
C
C*
        SUBROUTINE KLUDCMP(A,N,NP,INDX,D)
C***********************************************************************
C
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(4,4),INDX(4),VV(4)
C
        D=0.1E1
        D1=0.1E1
        DO 10 I=1,N
        AAMAX=0.0E0
        DO 20 J=1,N
        IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
20      CONTINUE
        IF(AAMAX.EQ.0.0D0)THEN
	WRITE(*,*) 'SINGULAR MATRIX'
	STOP 1111
	ENDIF 
        VV(I)=D1/AAMAX
10      CONTINUE
        DO 30 J=1,N
        DO 40 I=1,J-1
        SUM=A(I,J)
        DO 50 K=1,I-1
        SUM=SUM-A(I,K)*A(K,J)
50      CONTINUE
        A(I,J)=SUM
40      CONTINUE
        AAMAX=0.0D0
        DO 60 I=J,N
        SUM=A(I,J)
        DO 70 K=1,J-1
        SUM=SUM-A(I,K)*A(K,J)
70      CONTINUE
        A(I,J)=SUM
        DUM=VV(I)*ABS(SUM)
        IF(DUM.GE.AAMAX) THEN
        IMAX=I
        AAMAX=DUM
        ENDIF
60      CONTINUE
        IF(J.NE.IMAX) THEN
        DO 80 K=1,N
        DUM=A(IMAX,K)
        A(IMAX,K)=A(J,K)
        A(J,K)=DUM
80      CONTINUE
        D=-D1*D
        VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.0D0) A(J,J)=TINY
        IF(J.NE.N) THEN
        DUM=D1/A(J,J)
        DO 90 I=J+1,N
        A(I,J)=A(I,J)*DUM
90      CONTINUE
        ENDIF
30      CONTINUE
C
        RETURN
        END
C
C*
        SUBROUTINE KMATADD(A,B,IS,C)
C***********************************************************************
C
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(4,4),B(4,4),C(4,4)
        DO 80 I=1,IS
        DO 80 J=1,IS
80      C(I,J)=A(I,J)+B(I,J)
        RETURN
        END
C
C*
        SUBROUTINE KMULBTC(XLEVP,BM,BTC,NSIZB,IS)
C***********************************************************************
C
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION XLEVP(IS,IS),BM(IS,NSIZB),BTC(NSIZB,IS)
C
        DO 10  J=1,NSIZB
        DO 20  I=1,IS
        SUM=0.0E0
        DO 30  K=1,IS
        SUM=SUM+BM(K,J)*XLEVP(K,I)
30      CONTINUE
        BTC(J,I)=SUM
20      CONTINUE
10      CONTINUE
C
        RETURN
        END
C
C*
        SUBROUTINE KMULBTV(BM,S,NSIZB,IS,C)
C***********************************************************************
C*
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION BM(IS,NSIZB),S(IS),C(NSIZB)
C
        DO 10 J = 1, NSIZB
        X = 0.0
        DO 20 I = 1, IS
 20     X = X + BM(I,J)*S(I)
        C(J) = X
 10     CONTINUE
        RETURN
        END
C*
C*
        SUBROUTINE KMULTB(BC,BM,NSIZB,IS,BCB)
C***********************************************************************
C
C THIS SUBROUTINE MULTIPLIES TWO RECTANGULAR MATRICES
C
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION BC(NSIZB,IS),BM(IS,NSIZB),BCB(NSIZB,NSIZB)
C
        DO 10 J=1,NSIZB
        DO 20 I=1,NSIZB
        SUM=0.0E0
        DO 30 K=1,IS
        SUM=SUM+BC(I,K)*BM(K,J)
 30     CONTINUE
        BCB(I,J)=SUM
 20     CONTINUE
 10     CONTINUE
        RETURN
        END
C*
C*
        SUBROUTINE KMVMLT1(A,B,IS,C)
C***********************************************************************
C
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION  A(4,4),B(4),C(4)
        DO 40 I=1,IS
        X=0.0E0
        DO 50 J=1,IS
50      X=X+A(I,J)*B(J)
        C(I)=X
40      CONTINUE
        RETURN
        END
C
C*
        SUBROUTINE KMVMLT2(A,B,NSIZB,IS,C)
C***********************************************************************
C*
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(NSIZB,IS),B(IS),C(NSIZB)
        DO 40 I=1,NSIZB
        X=0.0E0
        DO 50 J=1,IS
50      X=X+A(I,J)*B(J)
        C(I)=X
40      CONTINUE
        RETURN
        END
C
C*
        SUBROUTINE KNULL(A,N)
C****************************************************************
C*
C****************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(N)
        DO 10 I=1,N
        A(I) = 0.0
10      CONTINUE
        RETURN
        END
C
C*
C***********************************************************************
C
C***********************************************************************
C
      SUBROUTINE KNULL1(A,IS)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(4,4)
      DO 20 I=1,IS
      DO 20 J=1,IS
   20 A(I,J)=0.0E0
      RETURN
      END
        SUBROUTINE KFORMMATX(A,B,IS,C)
C***********************************************************************
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(4),B(4),C(4,4)
        DO 90 I=1,IS
        DO 90 J=1,IS
90      C(I,J)=A(I)*B(J)
        RETURN
        END
        SUBROUTINE KDMEMCP(A1,A2,LEN)
C*******************************************************************
C
C*******************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A1(LEN),A2(LEN)
        DO 10 I =1, LEN
        A2(I) = A1(I)
10      CONTINUE
        RETURN
        END                                
