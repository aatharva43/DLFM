      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 CMNAME
C
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3), TIME(2)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C     
      PARAMETER (ZERO =0.0D0, HALF = 0.5D0, ONE =1.0D0, TWO=2.0D0,
     1           THREE=3.0D0, FOUR = 4.0D0, FIVE=5.0D0, SIX=6.0D0, 
     2           TOLER = 1.0E-8)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
      INTEGER K1
      INTEGER K2
C
      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER L
C
      INTEGER TAG
C
      REAL*8 AnIso1
      REAL*8 AnIso2
      REAL*8 AnIso3
      REAL*8 AnIsoCri
      REAL*8 AlphaD
      REAL*8 AlphaV
      REAL*8 Beta
      REAL*8 EmodAdhesion
      REAL*8 EmodCell1
      REAL*8 EmodCell2
      REAL*8 EmodCF
      REAL*8 EmodCF1
      REAL*8 EmodCF2
      REAL*8 EmodCF3
      REAL*8 EmodF
      REAL*8 EmodLamin
      REAL*8 EnuAdhesion
      REAL*8 EnuCell1
      REAL*8 EnuCell2
      REAL*8 EnuLamin
      REAL*8 Gmod
      REAL*8 GmodC1
      REAL*8 GmodC2
      REAL*8 GpBarMod
      REAL*8 Kmod
      REAL*8 KmodC1
      REAL*8 KmodC2
      REAL*8 KNorm
      REAL*8 KpBarMod
      REAL*8 KRhoMag
      REAL*8 KStressPrinMax
      REAL*8 KStressPrinMid
      REAL*8 KStressPrinMin
      REAL*8 MM
      REAL*8 MMCell
      REAL*8 NN
      REAL*8 NNCell
      REAL*8 PrinEmodVal1
      REAL*8 PrinEmodVal2
      REAL*8 PrinEmodVal3
      REAL*8 PrinStrainVal1
      REAL*8 PrinStrainVal2
      REAL*8 PrinStrainVal3
      REAL*8 PrinStressVal1
      REAL*8 PrinStressVal2
      REAL*8 PrinStressVal3
      REAL*8 Ratio
      REAL*8 RatioI
      REAL*8 RatioCell
      REAL*8 Rho0
      REAL*8 Rho0Bar
      REAL*8 StranC
      REAL*8 StranCCell
      REAL*8 StressRatio
C
      REAL*8 , DIMENSION(3,3,3,3)         :: Dtensor
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CContMatC
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CIsoMat
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CIsoMatC1
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CIsoMatC2
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CFibMat
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CMatC
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CMatC1
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CMatC2
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CMatCInv
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CMatC1Inv
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CMatC2Inv
      REAL*8 , DIMENSION(3,3)             :: EigenProj1
      REAL*8 , DIMENSION(3,3)             :: EigenProj2
      REAL*8 , DIMENSION(3,3)             :: EigenProj3
      REAL*8 , DIMENSION(1:20)            :: HDACcytoplasm
      REAL*8 , DIMENSION(1:20)            :: HDACnucleus
      REAL*8 , DIMENSION(3,3,3,3)         :: IdSym
      REAL*8 , DIMENSION(3,3)             :: Imat
      REAL*8 , DIMENSION(1:2*NTENS)       :: KDStran
      REAL*8 , DIMENSION(NTENS,NTENS)     :: KEyeMat
      REAL*8 , DIMENSION(2*NTENS,2*NTENS) :: KJacob
      REAL*8 , DIMENSION(1:2*NTENS)       :: KRes
      REAL*8 , DIMENSION(1:2*NTENS)       :: KResNeg
      REAL*8 , DIMENSION(1:NTENS)         :: KRho
      REAL*8 , DIMENSION(1:3)             :: KRhoVal
      REAL*8 , DIMENSION(3,3)             :: KRhoVec
      REAL*8 , DIMENSION(1:2*NTENS)       :: KStran
      REAL*8 , DIMENSION(1:NTENS)         :: KStressC1
      REAL*8 , DIMENSION(1:NTENS)         :: KStressC2
      REAL*8 , DIMENSION(1:NTENS)         :: KStressCompr
      REAL*8 , DIMENSION(1:NTENS)         :: KStressH
      REAL*8 , DIMENSION(1:NTENS)         :: KStressF
      REAL*8 , DIMENSION(1:3)             :: KStressVal
      REAL*8 , DIMENSION(3,3)             :: KStressVec
      REAL*8 , DIMENSION(1:3)             :: KStressComprVal
      REAL*8 , DIMENSION(3,3)             :: KStressComprVec
      REAL*8 , DIMENSION(1:NTENS)         :: KTotStran
      REAL*8 , DIMENSION(1:NTENS)         :: KTotStrain
      REAL*8 , DIMENSION(1:3)             :: KVal
      REAL*8 , DIMENSION(3,3)             :: KVec
      REAL*8 , DIMENSION(1:21)            :: MKLcytoplasm
      REAL*8 , DIMENSION(1:21)            :: MKLnucleus
      REAL*8 , DIMENSION(1:3)             :: PrinStrainVec1
      REAL*8 , DIMENSION(1:3)             :: PrinStrainVec2
      REAL*8 , DIMENSION(1:3)             :: PrinStrainVec3
      REAL*8 , DIMENSION(3,3)             :: StrainMat
      REAL*8 , DIMENSION(3,3)             :: StressMat
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C properties from input file
C
      ! cell parameters (cytoskeleton)
C
      TAG            = PROPS(1)    ! switch parameter
      EmodCell1      = PROPS(2)    ! elastic modulus of contractile part (KPa)
      EmodCell2      = PROPS(3)    ! elastic modulus of fibrous part     (KPa)
      EnuCell1       = PROPS(4)    ! poisson's ratio of contractile part
      EnuCell2       = PROPS(5)    ! poisson's ratio of fibrous part
      Beta           = PROPS(6)    ! chemical stiffness (1/KPa)
      AlphaV         = PROPS(7)    ! volumetric chemo-mechanical feedback parameter (1/KPa)
      AlphaD         = PROPS(8)    ! deviatoric chemo-mechanical feedback parameter (1/KPa) 
      Rho0I          = PROPS(9)    ! initial motor density in the quiescent state (KPa)
C
      ! lamin parameters
C
      EmodLamin      = PROPS(10)   ! elastic modulus (KPa)
      EnuLamin       = PROPS(11)   ! poisson's ratio
C
      ! adhesion layer parameters
C
      EmodAdhesion   = PROPS(12)   ! elastic modulus (KPa)
      EnuMdhesion    = PROPS(13)   ! poisson's ratio
C
      ! chromatin parameters
C
      EmodChromatin = PROPS(14)   ! elastic modulus (KPa)
      EnuChromatin  = PROPS(15)   ! poisson's ratio
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C cell or lamin or adhesion layer or chromatin
C
      IF (TAG.EQ.2) THEN
          GO TO 555         ! lamin
      ELSE IF (TAG.EQ.3) THEN
          GO TO 666         ! adhesion layer
      ELSE IF (TAG.EQ.4) THEN
          GO TO 777         ! chromatin
      END IF
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ------------------------    cell    ----------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C cell contractility increases with anisotropy in tension
C  
      Rho0 = Rho0I + 2.0D0 * STATEV(10+4*NTENS+4)
C      
C       Rho0 = Rho0I + 5.0D0 * STATEV(10+4*NTENS+4)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C      
C fiber elastic modulus (second component)
C
c      EmodCF = 1000.0D0 * EmodCell2 * (ONE-EnuCell2) / 
c     1                              ((ONE+EnuCell2)*(ONE-TWO*EnuCell2))
C      
      EmodCF = 10.0D0 * EmodCell2 * (ONE-EnuCell2) / 
     1                              ((ONE+EnuCell2)*(ONE-TWO*EnuCell2))
C
C critical strain (second component)
C
c      StranCCell = 0.060D0
C      
      StranCCell = 0.015D0
C
C other parameters for strain stiffening (second component)
C 
      RatioCell = 0.25D0
c      MMCell    = 300
c      NNCell    = 300
C      
      MMCell    = 10
      NNCell    = 10
C
C bulk modulus (first and second components)
C
      KmodC1 = EmodCell1 / ( THREE * (ONE-TWO*EnuCell1) )
      KmodC2 = EmodCell2 / ( THREE * (ONE-TWO*EnuCell2) )
C
C shear modulus (first and second components)
C
      GmodC1 = EmodCell1 / ( TWO * (ONE+EnuCell1) )
      GmodC2 = EmodCell2 / ( TWO * (ONE+EnuCell2) )
C
C effective bulk modulus for motor density (first component)
C
      KpBarMod = (ONE/THREE) * ( THREE * KmodC1 * AlphaV - ONE ) /
     1                         ( Beta  - AlphaV )
C
C effective shear modulus for motor density (first component)
C     
      GpBarMod = (ONE/TWO)  *  ( TWO   * GmodC1 * AlphaD - ONE ) /
     1                         ( Beta  - AlphaD )
C      
C effective contractility (first component)
C
      Rho0Bar = Beta * Rho0 * (TEMP+DTEMP) / ( Beta - AlphaV )     
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C identity matrix (3 by 3)
C
      DO K1 = 1,3
          DO K2 = 1,3
              Imat(K1,K2) = ZERO
          END DO
      END DO
C
      Imat(1,1) = ONE
      Imat(2,2) = ONE
      Imat(3,3) = ONE
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C symmetric identity fourth order tensor
C
      CALL IdentitySymmetric(IdSym)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C  
C identity matrix (NTENS by NTENS)
C     
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              KEyeMat(K1,K2) = ZERO
          END DO
      END DO
C
      DO K1 = 1,NTENS
          KEyeMat(K1,K1) = ONE
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C isotropic elastic stiffness matrix (first component)
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMatC1(K1,K2) = ZERO
          END DO
      END DO
      !
      CIsoMatC1(1,1) = KmodC1 + (FOUR/THREE)*GmodC1 
      CIsoMatC1(2,2) = KmodC1 + (FOUR/THREE)*GmodC1
      CIsoMatC1(3,3) = KmodC1 + (FOUR/THREE)*GmodC1
      CIsoMatC1(4,4) = GmodC1
      CIsoMatC1(1,2) = KmodC1 - (TWO/THREE) *GmodC1
      CIsoMatC1(1,3) = KmodC1 - (TWO/THREE) *GmodC1
      CIsoMatC1(2,3) = KmodC1 - (TWO/THREE) *GmodC1
      !
      IF (NSHR.EQ.3) THEN
          CIsoMatC1(5,5) = GmodC1
          CIsoMatC1(6,6) = GmodC1
      END IF
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMatC1(K2,K1) = CIsoMatC1(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C isotropic elastic stiffness matrix (second component)
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMatC2(K1,K2) = ZERO
          END DO
      END DO
      !
      CIsoMatC2(1,1) = KmodC2 + (FOUR/THREE)*GmodC2 
      CIsoMatC2(2,2) = KmodC2 + (FOUR/THREE)*GmodC2
      CIsoMatC2(3,3) = KmodC2 + (FOUR/THREE)*GmodC2
      CIsoMatC2(4,4) = GmodC2
      CIsoMatC2(1,2) = KmodC2 - (TWO/THREE) *GmodC2
      CIsoMatC2(1,3) = KmodC2 - (TWO/THREE) *GmodC2
      CIsoMatC2(2,3) = KmodC2 - (TWO/THREE) *GmodC2
      !
      IF (NSHR.EQ.3) THEN
          CIsoMatC2(5,5) = GmodC2
          CIsoMatC2(6,6) = GmodC2
      END IF
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMatC2(K2,K1) = CIsoMatC2(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C contractile elastic stiffness matrix (first component)
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CContMatC(K1,K2) = ZERO
          END DO
      END DO
      !
      CContMatC(1,1) = KpBarMod + (FOUR/THREE)*GpBarMod 
      CContMatC(2,2) = KpBarMod + (FOUR/THREE)*GpBarMod
      CContMatC(3,3) = KpBarMod + (FOUR/THREE)*GpBarMod
      CContMatC(4,4) = GpBarMod
      CContMatC(1,2) = KpBarMod - (TWO/THREE) *GpBarMod
      CContMatC(1,3) = KpBarMod - (TWO/THREE) *GpBarMod
      CContMatC(2,3) = KpBarMod - (TWO/THREE) *GpBarMod
      !
      IF (NSHR.EQ.3) THEN
          CContMatC(5,5) = GpBarMod
          CContMatC(6,6) = GpBarMod
      END IF
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CContMatC(K2,K1) = CContMatC(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C total elastic stiffness (first component)
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CMatC1(K1,K2) = ZERO
          END DO
      END DO
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CMatC1(K1,K2) = CIsoMatC1(K1,K2) + CContMatC(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C    
C total strain
C
      DO K1 = 1,NTENS
          KTotStran(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          KTotStran(K1) = STRAN(K1) + DSTRAN(K1)
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C unknown strains (initial guess)
C
      DO K1 = 1,2*NTENS
          KStran(K1) = ZERO
      END DO
C
      DO K1 = 1,NTENS
          KStran(K1)       =  HALF * KTotStran(K1)      ! first  component
          KStran(K1+NTENS) = -HALF * KTotStran(K1)      ! second component
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C loop to solve the set of nonlinear equations
C
C intial norm
      KNorm = ONE
C
      DO WHILE (KNorm.GT.TOLER)
          !
          ! call residual function
          CALL KResidual(CMatC1 , CIsoMatC2 , KStran , KTotStran ,
     1 IdSym , Imat , KEyeMat ,
     2 NDI, NSHR, NTENS , KINC , NOEL , NPT ,
     3 EmodCF , StranCCell , RatioCell , MMCell , NNCell , Rho0Bar ,
     4 KStressC1 , KStressC2 , CMatC2 , KJacob , KRes ,
     5 PrinEmodVal1C , PrinEmodVal2C , PrinEmodVal3C)
          !
          ! norm of rediusal
          IF (NSHR.EQ.1) THEN
              KNorm = SQRT( KRes(1)  * KRes(1) + KRes(2)  * KRes(2)  +
     1                      KRes(3)  * KRes(3) + KRes(4)  * KRes(4)  +
     2                      KRes(5)  * KRes(5) + KRes(6)  * KRes(6)  +
     3                      KRes(7)  * KRes(7) + KRes(8)  * KRes(8)  )
         ELSE
             KNorm = SQRT( KRes(1)  * KRes(1)  + KRes(2)  * KRes(2)  +
     1                     KRes(3)  * KRes(3)  + KRes(4)  * KRes(4)  +
     2                     KRes(5)  * KRes(5)  + KRes(6)  * KRes(6)  +
     3                     KRes(7)  * KRes(7)  + KRes(8)  * KRes(8)  +
     4                     KRes(9)  * KRes(9)  + KRes(10) * KRes(10) +
     5                     KRes(11) * KRes(11) + KRes(12) * KRes(12) )
         END IF
         !
         ! adding negative sign to the residual vector
         DO K1 = 1,2*NTENS
             KResNeg(K1) = ZERO
         END DO
         !
         DO K1 = 1,2*NTENS
             KResNeg(K1) = -KRes(K1)
         END DO
         !
         ! solve nonlinear equation 
         CALL LU(2*NTENS,KJacob,KResNeg,KDStran)
         !
         ! update strains
         DO K1 = 1,2*NTENS
             KStran(K1) = KStran(K1) + KDStran(K1)
         END DO
         !
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C     
C inverse of stiffness matrix (first component)
C
      CALL KINV( NTENS , CMatC1 , CMatC1Inv )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C     
C inverse of stiffness matrix (second component)
C
      CALL KINV( NTENS , CMatC2 , CMatC2Inv )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C     
C total stiffness matrix (first and second components)
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CMatCInv(K1,K2) = ZERO
          END DO
      END DO
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CMatCInv(K1,K2) = CMatC1Inv(K1,K2) + CMatC2Inv(K1,K2)
          END DO
      END DO
C
      CALL KINV( NTENS , CMatCInv , CMatC )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C motor density
C
      DO K1 = 1,NTENS
          KRho(K1) = ZERO
      END DO
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              KRho(K1) = KRho(K1) + CContMatC(K1,K2) * KStran(K2)
          END DO
      END DO
      !
      KRho(1) = KRho(1) + Rho0Bar
      KRho(2) = KRho(2) + Rho0Bar
      KRho(3) = KRho(3) + Rho0Bar 
C     
C eigen values and eigen vectors of motor density 
C
      CALL SPRIND(KRho, KRhoVal, KRhoVec, 1, NDI, NSHR)
C
C motor polarization magnitude
C
      KRhoMag = HALF * ( MAX(KRhoVal(1),KRhoVal(2),KRhoVal(3)) - 
     1                   MIN(KRhoVal(1),KRhoVal(2),KRhoVal(3)) )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stress on the microtubule network
C
      DO K1 = 1,NTENS
          KStressCompr(K1) = ZERO
      END DO
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              KStressCompr(K1) = KStressCompr(K1) + 
     1                           CIsoMatC1(K1,K2) * KStran(K2)
          END DO
      END DO 
C
C eigen values and eigen vectors of stress on the microtubule network
C
      CALL SPRIND(KStressCompr, KStressComprVal, KStressComprVec, 1,
     1            NDI, NSHR)
C     
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C eigen values and eigen vectors of stress
C
      CALL SPRIND(KStressC1, KStressVal, KStressVec, 1, NDI, NSHR)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C HDAC in cytoplasm
C
C rectangular substrate geometry      
      HDACcytoplasm(1)  = 1.000000000000000
      HDACcytoplasm(2)  = 1.000414356930562
      HDACcytoplasm(3)  = 1.002234602104853
      HDACcytoplasm(4)  = 1.006543624773674
      HDACcytoplasm(5)  = 1.014560935117034
      HDACcytoplasm(6)  = 1.027308066787507
      HDACcytoplasm(7)  = 1.045890374916914
      HDACcytoplasm(8)  = 1.071195151069322
      HDACcytoplasm(9)  = 1.104180838179775
      HDACcytoplasm(10) = 1.145787273270653 
      HDACcytoplasm(11) = 1.187122752477759
      HDACcytoplasm(12) = 1.218201624195104
      HDACcytoplasm(13) = 1.241300164441747
      HDACcytoplasm(14) = 1.257844635512917
      HDACcytoplasm(15) = 1.269319046689393
      HDACcytoplasm(16) = 1.276992804690958 
      HDACcytoplasm(17) = 1.281911177333282 
      HDACcytoplasm(18) = 1.284891240225833 
      HDACcytoplasm(19) = 1.286554546773318 
      HDACcytoplasm(20) = 1.287358311261062 
      HDACcytoplasm(21) = 1.288358311261062
C
C circular substrate geometry      
C      HDACcytoplasm(1)  = 1.000000000000000
C      HDACcytoplasm(2)  = 0.999585823428886
C      HDACcytoplasm(3)  = 0.997770945465787
C      HDACcytoplasm(4)  = 0.993504113712500
C      HDACcytoplasm(5)  = 0.985680265281714
C      HDACcytoplasm(6)  = 0.973553468859974
C      HDACcytoplasm(7)  = 0.956558373536667
C      HDACcytoplasm(8)  = 0.934684930064159
C      HDACcytoplasm(9)  = 0.908247116633826
C      HDACcytoplasm(10) = 0.878003678456830
C      HDACcytoplasm(11) = 0.851681255235354
C      HDACcytoplasm(12) = 0.835309254883147
C      HDACcytoplasm(13) = 0.825417275357930
C      HDACcytoplasm(14) = 0.819913176383588
C      HDACcytoplasm(15) = 0.817173584718320
C      HDACcytoplasm(16) = 0.816069591812066
C      HDACcytoplasm(17) = 0.815865907151123
C      HDACcytoplasm(18) = 0.816095491783992
C      HDACcytoplasm(19) = 0.816493567703002
C      HDACcytoplasm(20) = 0.816909695165019
C      HDACcytoplasm(21) = 0.816409695165019
C     
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C MKL in cytoplasm
C
C rectangulat substrate geometry 
      MKLcytoplasm(1)   = 1.000000000000000
      MKLcytoplasm(2)   = 0.999729095500685
      MKLcytoplasm(3)   = 0.998532166370780
      MKLcytoplasm(4)   = 0.995681961002079
      MKLcytoplasm(5)   = 0.990360053226437
      MKLcytoplasm(6)   = 0.981908416344205
      MKLcytoplasm(7)   = 0.969684743485849
      MKLcytoplasm(8)   = 0.953319050836224
      MKLcytoplasm(9)   = 0.932568486007019
      MKLcytoplasm(10)  = 0.907442772228997
      MKLcytoplasm(11)  = 0.883406796634960
      MKLcytoplasm(12)  = 0.865657756976614
      MKLcytoplasm(13)  = 0.852509966481167
      MKLcytoplasm(14)  = 0.843048027212708
      MKLcytoplasm(15)  = 0.836403129582466
      MKLcytoplasm(16)  = 0.831875699529139
      MKLcytoplasm(17)  = 0.828882106879157
      MKLcytoplasm(18)  = 0.826981258912319
      MKLcytoplasm(19)  = 0.825822511302355
      MKLcytoplasm(20)  = 0.825157406023544
      MKLcytoplasm(21)  = 0.825057406023544
C
C circular substrate geometry 
C      MKLcytoplasm(1)   = 1.000000000000000
C      MKLcytoplasm(2)   = 1.000270901542544
C      MKLcytoplasm(3)   = 1.001467537402632
C      MKLcytoplasm(4)   = 1.004318869298658
C      MKLcytoplasm(5)   = 1.009645355501232
C      MKLcytoplasm(6)   = 1.018121215905421
C      MKLcytoplasm(7)   = 1.030417143733998
C      MKLcytoplasm(8)   = 1.046970799420638
C      MKLcytoplasm(9)   = 1.068131876777390
C      MKLcytoplasm(10)  = 1.094054687731112
C      MKLcytoplasm(11)  = 1.119393969316233
C      MKLcytoplasm(12)  = 1.138905762874322
C      MKLcytoplasm(13)  = 1.154224767725093
C      MKLcytoplasm(14)  = 1.166230156402216
C      MKLcytoplasm(15)  = 1.175685889142103
C      MKLcytoplasm(16)  = 1.183177038904576
C      MKLcytoplasm(17)  = 1.189162431009368
C      MKLcytoplasm(18)  = 1.193983556241833
C      MKLcytoplasm(19)  = 1.197908893618100
C      MKLcytoplasm(20)  = 1.201134470982385
C      MKLcytoplasm(21)  = 1.202134470982385
C     
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C state variables
C
C motor density
C
      DO K1 = 1,NTENS
          STATEV(K1) = KRho(K1)
      END DO
C 
C total stiffness matrix
C
      DO K1 = 1,NTENS
          STATEV(K1+NTENS) = CMatC(K1,K1)
      END DO
C
C eigen values of motor density
C
      STATEV(2*NTENS+1) = KRhoVal(1)
      STATEV(2*NTENS+2) = KRhoVal(2)
      STATEV(2*NTENS+3) = KRhoVal(3)
      STATEV(2*NTENS+4) = KRhoMag
C
C HDAC in cytoplasm (and also HDAC in nucleus)
C
      IF (KSTEP.EQ.1) THEN
          STATEV(1+4*NTENS+4) = HDACcytoplasm(KINC) 
      ELSE
C          STATEV(1+4*NTENS+4) = HDACcytoplasm(KINC+10)         ! low contractility
          STATEV(1+4*NTENS+4) = HDACcytoplasm(KINC+11)         ! high contractility
      END IF
C
C MKL in cytoplasm (and also MKL in nucleus)
C
      IF (KSTEP.EQ.1) THEN
          STATEV(2+4*NTENS+4) = MKLcytoplasm(KINC)
      ELSE
C          STATEV(2+4*NTENS+4) = MKLcytoplasm(KINC+10)           ! low contractility
          STATEV(2+4*NTENS+4) = MKLcytoplasm(KINC+11)           ! high contractility
      END IF
C
C maximum total stiffness 
C
      STATEV(3+4*NTENS+4) = MAX(CMatC(1,1),CMatC(2,2),CMatC(3,3))
C
C maximum stiffness of the second component  
C      
      STATEV(4+4*NTENS+4) = MAX(PrinEmodVal1C , PrinEmodVal2C , 
     1                          PrinEmodVal3C) + EmodCell2
C
      STATEV(5+4*NTENS+4) = PrinEmodVal1C + EmodCell2
      STATEV(6+4*NTENS+4) = PrinEmodVal2C + EmodCell2
      STATEV(7+4*NTENS+4) = PrinEmodVal3C + EmodCell2
C
C von Mises stress
C
      STATEV(8+4*NTENS+4) = ( HALF*((KStressC1(1)-KStressC1(2))**TWO) +
     1                        HALF*((KStressC1(2)-KStressC1(3))**TWO) +
     2                        HALF*((KStressC1(3)-KStressC1(1))**TWO) +
     3                        THREE*KStressC1(4)*KStressC1(4)         + 
     4                        THREE*KStressC1(5)*KStressC1(5)         +
     5                        THREE*KStressC1(6)*KStressC1(6) ) ** HALF
C
C polarization criterion
C
      KStressPrinMax = MAX(KStressVal(1),KStressVal(2),KStressVal(3))
      KStressPrinMin = MIN(KStressVal(1),KStressVal(2),KStressVal(3))
      DO K1 = 1,3
          IF ( (KStressVal(K1).GT.KStressPrinMin).AND.
     1         (KStressVal(K1).LT.KStressPrinMax) ) THEN
              KStressPrinMid = KStressVal(K1)
          END IF
      END DO
C
      IF (KStressPrinMax.LE.0.0D0) THEN
          STATEV(9+4*NTENS+4) = 0.0D0
      ELSE IF (KStressPrinMid.GT.0.0D0) THEN 
          StressRatio = HALF * (KStressPrinMax / KStressPrinMid) - HALF
          STATEV(9+4*NTENS+4) = TANH(StressRatio)
      ELSE
          STATEV(9+4*NTENS+4) = 1.0D0
      END IF    
C 
C aspect ratio of tension anisotropy
C
      STATEV(10+4*NTENS+4) = STATEV(9+4*NTENS+4) * Rho0I
C
C contractile elastic stiffness matrix (first component)
C     
      STATEV(12+4*NTENS+4) = CMatC1(1,1)
      STATEV(12+4*NTENS+4) = CMatC1(2,2)
      STATEV(13+4*NTENS+4) = CMatC1(3,3)
C
C fibrous elastic stiffness matrix (second component)
C
      STATEV(14+4*NTENS+4) = CMatC2(1,1)
      STATEV(15+4*NTENS+4) = CMatC2(2,2)
      STATEV(16+4*NTENS+4) = CMatC2(3,3)
C
C principal stresses
C
      STATEV(17+4*NTENS+4) = KStressPrinMax
      STATEV(18+4*NTENS+4) = KStressPrinMid
      STATEV(19+4*NTENS+4) = KStressPrinMin  
C
C mean contractility
C
      STATEV(20+4*NTENS+4) = ( KRho(1) + KRho(2) + KRho(3) ) / THREE
      STATEV(72) = ( KStressPrinMax + KStressPrinMid + KStressPrinMin ) / THREE
      STATEV(73) = KStressPrinMax
      STATEV(74) = ( KRhoVal(1) + KRhoVal(2) + KRhoVal(3) ) / THREE
      STATEV(75) = MAX( KRhoVal(1),KRhoVal(2),KRhoVal(3) )
C   
C
C strain in tensile and contractile components
C
      STATEV(21+4*NTENS+4) = KStran(1)
      STATEV(22+4*NTENS+4) = KStran(2)
      STATEV(23+4*NTENS+4) = KStran(3)
      STATEV(24+4*NTENS+4) = KStran(4)
      STATEV(25+4*NTENS+4) = KStran(5)
      STATEV(26+4*NTENS+4) = KStran(6)
      STATEV(27+4*NTENS+4) = KStran(7)
      STATEV(28+4*NTENS+4) = KStran(8)
      STATEV(29+4*NTENS+4) = KStran(9)
      STATEV(30+4*NTENS+4) = KStran(10)
      STATEV(31+4*NTENS+4) = KStran(11)
      STATEV(32+4*NTENS+4) = KStran(12)
C
C stress on the microtubule network
C
      DO K1 = 1,NTENS
          STATEV(32+4*NTENS+4+K1) = KStressCompr(K1)
      END DO
C
      STATEV(32+4*NTENS+4+NTENS+1) = KStressComprVal(1)
      STATEV(32+4*NTENS+4+NTENS+2) = KStressComprVal(2)
      STATEV(32+4*NTENS+4+NTENS+3) = KStressComprVal(3)
C
      STATEV(32+4*NTENS+4+NTENS+4) = -STATEV(32+4*NTENS+4+NTENS+2) 
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stress
C
      DO K1 = 1,NTENS
          STRESS(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          STRESS(K1) = KStressC1(K1)
      END DO
C
C jacobean
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              DDSDDE(K1,K2) = ZERO
          END DO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              DDSDDE(K1,K2) = CMatC(K1,K2)
          END DO
      END DO
C      
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C end of the cell part
C 
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
      GO TO 888
C      
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C lamin (beginning of the lamin part)
C
555   CONTINUE
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ------------------------    lamin    --------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C fiber elastic modulus
C
      EmodF = 250.0D0 * EmodLamin * (ONE-EnuLamin) 
     1                          / ((ONE+EnuLamin)*(ONE-TWO*EnuLamin))
C
C critical strain
C
      StranC = 0.020D0
C
C other parameters for strain stiffening
C 
      Ratio = 0.25D0
      MM    = 100
      NN    = 5
C
C bulk modulus
C
      Kmod = EmodLamin / ( THREE * (ONE-TWO*EnuLamin) )
C
C shear modulus
C
      Gmod = EmodLamin / ( TWO * (ONE+EnuLamin) )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C   
C isotropic elastic stiffness matrix
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMat(K1,K2) = ZERO
          END DO
      END DO
      !
      CIsoMat(1,1) = Kmod + (FOUR/THREE)*Gmod 
      CIsoMat(2,2) = Kmod + (FOUR/THREE)*Gmod
      CIsoMat(3,3) = Kmod + (FOUR/THREE)*Gmod
      CIsoMat(4,4) = Gmod
      CIsoMat(1,2) = Kmod - (TWO/THREE) *Gmod
      CIsoMat(1,3) = Kmod - (TWO/THREE) *Gmod
      CIsoMat(2,3) = Kmod - (TWO/THREE) *Gmod
      !
      IF (NSHR.EQ.3) THEN
          CIsoMat(5,5) = Gmod
          CIsoMat(6,6) = Gmod
      END IF
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMat(K2,K1) = CIsoMat(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C total strain
C
      DO K1 = 1,NTENS
          KTotStrain(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          KTotStrain(K1) = STRAN(K1) + DSTRAN(K1)
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C hookean isotropic stress
C
      DO K1 = 1,NTENS
          KStressH(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              KStressH(K1) = KStressH(K1)   +
     1                       CIsoMat(K1,K2) * KTotStrain(K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C strain matrix
C
      StrainMat(1,1) = KTotStrain(1)
      StrainMat(2,2) = KTotStrain(2)
      StrainMat(3,3) = KTotStrain(3)
      StrainMat(1,2) = KTotStrain(4) * HALF
      !
      IF (NSHR.EQ.3) THEN
          StrainMat(1,3) = KTotStrain(5) * HALF
          StrainMat(2,3) = KTotStrain(6) * HALF
      END IF
      !
      DO K1 = 1,3
          DO K2 = 1,3
              StrainMat(K2,K1) = StrainMat(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C eigen values and eigen vectors of total strain 
C
      CALL SPRIND(KTotStrain, KVal, KVec, 2, NDI, NSHR)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C principal strain
C
      PrinStrainVal1 = KVal(1)
      PrinStrainVal2 = KVal(2)
      PrinStrainVal3 = KVal(3)
C
C principal vectors
C
      PrinStrainVec1(1) = KVec(1,1)
      PrinStrainVec1(2) = KVec(1,2)
      PrinStrainVec1(3) = KVec(1,3)
      !
      PrinStrainVec2(1) = KVec(2,1)
      PrinStrainVec2(2) = KVec(2,2)
      PrinStrainVec2(3) = KVec(2,3)
      !
      PrinStrainVec3(1) = KVec(3,1)
      PrinStrainVec3(2) = KVec(3,2)
      PrinStrainVec3(3) = KVec(3,3)  
C
C eigen projection
C
      CALL DyadVect(PrinStrainVec1 , PrinStrainVec1 , 3 , EigenProj1)
      CALL DyadVect(PrinStrainVec2 , PrinStrainVec2 , 3 , EigenProj2)
      CALL DyadVect(PrinStrainVec3 , PrinStrainVec3 , 3 , EigenProj3)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C principal stress
C     
      CALL dWdPrinStrainTension(PrinStrainVal1 , EmodF , StranC , 
     1                   MM , NN , Ratio, PrinStressVal1)
C
      CALL dWdPrinStrainTension(PrinStrainVal2 , EmodF , StranC , 
     1                   MM , NN , Ratio, PrinStressVal2)
C
      CALL dWdPrinStrainTension(PrinStrainVal3 , EmodF , StranC , 
     1                   MM , NN , Ratio, PrinStressVal3)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C derivative of principal stress with respect to principal strain
C     
      CALL d2WdPrinStrainTension(PrinStrainVal1 , EmodF , StranC , 
     1                   MM , NN , Ratio, PrinEmodVal1)
C
      CALL d2WdPrinStrainTension(PrinStrainVal2 , EmodF , StranC ,  
     1                   MM , NN , Ratio, PrinEmodVal2)
C
      CALL d2WdPrinStrainTension(PrinStrainVal3 , EmodF , StranC , 
     1                   MM , NN , Ratio, PrinEmodVal3)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stress matrix
C
      DO K1 = 1,3
          DO K2 = 1,3
              StressMat(K1,K2) = ZERO
          END DO
      END DO
C
      DO K1 = 1,3
          DO K2 = 1,3
              StressMat(K1,K2) = PrinStressVal1 * EigenProj1(K1,K2) +
     1                           PrinStressVal2 * EigenProj2(K1,K2) +
     2                           PrinStressVal3 * EigenProj3(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C fiber stress
C
      DO K1 = 1,NTENS
          KStressF(K1) = ZERO
      END DO
C         
      KStressF(1) = StressMat(1,1)
      KStressF(2) = StressMat(2,2)
      KStressF(3) = StressMat(3,3)
      KStressF(4) = StressMat(1,2)
      !
      IF (NSHR.EQ.3) THEN
          KStressF(5) = StressMat(1,3)
          KStressF(6) = StressMat(2,3)
      END IF
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C identity matrix
C
      DO K1 = 1,3
          DO K2 = 1,3
              Imat(K1,K2) = ZERO
          END DO
      END DO
C
      Imat(1,1) = ONE
      Imat(2,2) = ONE
      Imat(3,3) = ONE
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C symmetric identity fourth order tensor
C
      CALL IdentitySymmetric(IdSym)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stiffness tensor
C
      CALL StiffTens(PrinStrainVal1,PrinStrainVal2,PrinStrainVal3,
     1  PrinStressVal1 , PrinStressVal2 , PrinStressVal3 ,    
     2  PrinEmodVal1   , PrinEmodVal2   , PrinEmodVal3   ,    
     3  EigenProj1     , EigenProj2     , EigenProj3     ,
     4  IdSym          , Imat           , StrainMat      , Dtensor )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stiffness matrix
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CFibMat(K1,K2) = ZERO
          END DO
      END DO
C
      CFibMat(1,1) = Dtensor(1,1,1,1)
      CFibMat(1,2) = Dtensor(1,1,2,2)
      CFibMat(1,3) = Dtensor(1,1,3,3)
      CFibMat(1,4) = Dtensor(1,1,1,2)
      !
      CFibMat(2,1) = Dtensor(2,2,1,1)
      CFibMat(2,2) = Dtensor(2,2,2,2)
      CFibMat(2,3) = Dtensor(2,2,3,3)
      CFibMat(2,4) = Dtensor(2,2,1,2)
      !
      CFibMat(3,1) = Dtensor(3,3,1,1)
      CFibMat(3,2) = Dtensor(3,3,2,2)
      CFibMat(3,3) = Dtensor(3,3,3,3)
      CFibMat(3,4) = Dtensor(3,3,1,2)
      !
      CFibMat(4,1) = Dtensor(1,2,1,1)
      CFibMat(4,2) = Dtensor(1,2,2,2)
      CFibMat(4,3) = Dtensor(1,2,3,3)
      CFibMat(4,4) = Dtensor(1,2,1,2)
      !
      IF (NSHR.EQ.3) THEN
          !
          CFibMat(1,5) = Dtensor(1,1,1,3)
          CFibMat(1,6) = Dtensor(1,1,2,3)
          !
          CFibMat(2,5) = Dtensor(2,2,1,3)
          CFibMat(2,6) = Dtensor(2,2,2,3)
          !
          CFibMat(3,5) = Dtensor(3,3,1,3)
          CFibMat(3,6) = Dtensor(3,3,2,3)
          !
          CFibMat(4,5) = Dtensor(1,2,1,3)
          CFibMat(4,6) = Dtensor(1,2,2,3)
          !
          CFibMat(5,1) = Dtensor(1,3,1,1)
          CFibMat(5,2) = Dtensor(1,3,2,2)
          CFibMat(5,3) = Dtensor(1,3,3,3)
          CFibMat(5,4) = Dtensor(1,3,1,2)
          CFibMat(5,5) = Dtensor(1,3,1,3)
          CFibMat(5,6) = Dtensor(1,3,2,3)
          !
          CFibMat(6,1) = Dtensor(2,3,1,1)
          CFibMat(6,2) = Dtensor(2,3,2,2)
          CFibMat(6,3) = Dtensor(2,3,3,3)
          CFibMat(6,4) = Dtensor(2,3,1,2)
          CFibMat(6,5) = Dtensor(2,3,1,3)
          CFibMat(6,6) = Dtensor(2,3,2,3)
          !
          END IF
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C state variable
C
C total stiffness matrix
      DO K1 = 1,NTENS
          STATEV(K1+2*NTENS+4) = CIsoMat(K1,K1) + CFibMat(K1,K1)
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stress
C
      DO K1 = 1,NTENS
          STRESS(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          STRESS(K1) = KStressH(K1) + KStressF(K1)
      END DO
C
C jacobean
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              DDSDDE(K1,K2) = ZERO
          END DO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              DDSDDE(K1,K2) = CIsoMat(K1,K2) + CFibMat(K1,K2)
          END DO
      END DO
C  
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C end of the lamin part
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
      GO TO 888
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C adhesion layer (beginning of the adhesion part)
C
666   CONTINUE
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ------------------------    adhesion    ------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C fiber elastic modulus
C
      EmodF = 10.0D0 * EmodAdhesion * (ONE-EnuAdhesion) 
     1         / ((ONE+EnuAdhesion) * (ONE-TWO*EnuAdhesion))
C
C critical strain
C
      StranC = 0.02D0
C
C other parameters for strain stiffening
C 
      Ratio = 0.25D0
      MM    = 10
      NN    = 2
C
C bulk modulus
C
      Kmod = EmodAdhesion / ( THREE * (ONE-TWO*EnuAdhesion) )
C
C shear modulus
C
      Gmod = EmodAdhesion / ( TWO * (ONE+EnuAdhesion) )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C   
C isotropic elastic stiffness matrix
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMat(K1,K2) = ZERO
          END DO
      END DO
      !
      CIsoMat(1,1) = Kmod + (FOUR/THREE)*Gmod 
      CIsoMat(2,2) = Kmod + (FOUR/THREE)*Gmod
      CIsoMat(3,3) = Kmod + (FOUR/THREE)*Gmod
      CIsoMat(4,4) = Gmod
      CIsoMat(1,2) = Kmod - (TWO/THREE) *Gmod
      CIsoMat(1,3) = Kmod - (TWO/THREE) *Gmod
      CIsoMat(2,3) = Kmod - (TWO/THREE) *Gmod
      !
      IF (NSHR.EQ.3) THEN
          CIsoMat(5,5) = Gmod
          CIsoMat(6,6) = Gmod
      END IF
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMat(K2,K1) = CIsoMat(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C total strain
C
      DO K1 = 1,NTENS
          KTotStrain(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          KTotStrain(K1) = STRAN(K1) + DSTRAN(K1)
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C hookean isotropic stress
C
      DO K1 = 1,NTENS
          KStressH(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              KStressH(K1) = KStressH(K1)   +
     1                       CIsoMat(K1,K2) * KTotStrain(K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C strain matrix
C
      StrainMat(1,1) = KTotStrain(1)
      StrainMat(2,2) = KTotStrain(2)
      StrainMat(3,3) = KTotStrain(3)
      StrainMat(1,2) = KTotStrain(4) * HALF
      !
      IF (NSHR.EQ.3) THEN
          StrainMat(1,3) = KTotStrain(5) * HALF
          StrainMat(2,3) = KTotStrain(6) * HALF
      END IF
      !
      DO K1 = 1,3
          DO K2 = 1,3
              StrainMat(K2,K1) = StrainMat(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C eigen values and eigen vectors of total strain 
C
      CALL SPRIND(KTotStrain, KVal, KVec, 2, NDI, NSHR)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C principal strain
C
      PrinStrainVal1 = KVal(1)
      PrinStrainVal2 = KVal(2)
      PrinStrainVal3 = KVal(3)
C
C principal vectors
C
      PrinStrainVec1(1) = KVec(1,1)
      PrinStrainVec1(2) = KVec(1,2)
      PrinStrainVec1(3) = KVec(1,3)
      !
      PrinStrainVec2(1) = KVec(2,1)
      PrinStrainVec2(2) = KVec(2,2)
      PrinStrainVec2(3) = KVec(2,3)
      !
      PrinStrainVec3(1) = KVec(3,1)
      PrinStrainVec3(2) = KVec(3,2)
      PrinStrainVec3(3) = KVec(3,3)  
C
C eigen projection
C
      CALL DyadVect(PrinStrainVec1 , PrinStrainVec1 , 3 , EigenProj1)
      CALL DyadVect(PrinStrainVec2 , PrinStrainVec2 , 3 , EigenProj2)
      CALL DyadVect(PrinStrainVec3 , PrinStrainVec3 , 3 , EigenProj3)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C principal stress
C     
      CALL dWdPrinStrain(PrinStrainVal1 , EmodF , StranC , MM , NN , 
     1                   Ratio, PrinStressVal1)
C
      CALL dWdPrinStrain(PrinStrainVal2 , EmodF , StranC , MM , NN , 
     1                   Ratio, PrinStressVal2)
C
      CALL dWdPrinStrain(PrinStrainVal3 , EmodF , StranC , MM , NN , 
     1                   Ratio, PrinStressVal3)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C derivative of principal stress with respect to principal strain
C     
      CALL d2WdPrinStrain(PrinStrainVal1 , EmodF , StranC , MM , NN , 
     1                   Ratio, PrinEmodVal1)
C
      CALL d2WdPrinStrain(PrinStrainVal2 , EmodF , StranC , MM , NN , 
     1                   Ratio, PrinEmodVal2)
C
      CALL d2WdPrinStrain(PrinStrainVal3 , EmodF , StranC , MM , NN , 
     1                   Ratio, PrinEmodVal3)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stress matrix
C
      DO K1 = 1,3
          DO K2 = 1,3
              StressMat(K1,K2) = ZERO
          END DO
      END DO
C
      DO K1 = 1,3
          DO K2 = 1,3
              StressMat(K1,K2) = PrinStressVal1 * EigenProj1(K1,K2) +
     1                           PrinStressVal2 * EigenProj2(K1,K2) +
     2                           PrinStressVal3 * EigenProj3(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C fiber stress
C
      DO K1 = 1,NTENS
          KStressF(K1) = ZERO
      END DO
C         
      KStressF(1) = StressMat(1,1)
      KStressF(2) = StressMat(2,2)
      KStressF(3) = StressMat(3,3)
      KStressF(4) = StressMat(1,2)
      !
      IF (NSHR.EQ.3) THEN
          KStressF(5) = StressMat(1,3)
          KStressF(6) = StressMat(2,3)
      END IF
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C identity matrix
C
      DO K1 = 1,3
          DO K2 = 1,3
              Imat(K1,K2) = ZERO
          END DO
      END DO
C
      Imat(1,1) = ONE
      Imat(2,2) = ONE
      Imat(3,3) = ONE
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C symmetric identity fourth order tensor
C
      CALL IdentitySymmetric(IdSym)
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stiffness tensor
C
      CALL StiffTens(PrinStrainVal1,PrinStrainVal2,PrinStrainVal3,
     1  PrinStressVal1 , PrinStressVal2 , PrinStressVal3 ,    
     2  PrinEmodVal1   , PrinEmodVal2   , PrinEmodVal3   ,    
     3  EigenProj1     , EigenProj2     , EigenProj3     ,
     4  IdSym          , Imat           , StrainMat      , Dtensor )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stiffness matrix
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CFibMat(K1,K2) = ZERO
          END DO
      END DO
C
      CFibMat(1,1) = Dtensor(1,1,1,1)
      CFibMat(1,2) = Dtensor(1,1,2,2)
      CFibMat(1,3) = Dtensor(1,1,3,3)
      CFibMat(1,4) = Dtensor(1,1,1,2)
      !
      CFibMat(2,1) = Dtensor(2,2,1,1)
      CFibMat(2,2) = Dtensor(2,2,2,2)
      CFibMat(2,3) = Dtensor(2,2,3,3)
      CFibMat(2,4) = Dtensor(2,2,1,2)
      !
      CFibMat(3,1) = Dtensor(3,3,1,1)
      CFibMat(3,2) = Dtensor(3,3,2,2)
      CFibMat(3,3) = Dtensor(3,3,3,3)
      CFibMat(3,4) = Dtensor(3,3,1,2)
      !
      CFibMat(4,1) = Dtensor(1,2,1,1)
      CFibMat(4,2) = Dtensor(1,2,2,2)
      CFibMat(4,3) = Dtensor(1,2,3,3)
      CFibMat(4,4) = Dtensor(1,2,1,2)
      !
      IF (NSHR.EQ.3) THEN
          !
          CFibMat(1,5) = Dtensor(1,1,1,3)
          CFibMat(1,6) = Dtensor(1,1,2,3)
          !
          CFibMat(2,5) = Dtensor(2,2,1,3)
          CFibMat(2,6) = Dtensor(2,2,2,3)
          !
          CFibMat(3,5) = Dtensor(3,3,1,3)
          CFibMat(3,6) = Dtensor(3,3,2,3)
          !
          CFibMat(4,5) = Dtensor(1,2,1,3)
          CFibMat(4,6) = Dtensor(1,2,2,3)
          !
          CFibMat(5,1) = Dtensor(1,3,1,1)
          CFibMat(5,2) = Dtensor(1,3,2,2)
          CFibMat(5,3) = Dtensor(1,3,3,3)
          CFibMat(5,4) = Dtensor(1,3,1,2)
          CFibMat(5,5) = Dtensor(1,3,1,3)
          CFibMat(5,6) = Dtensor(1,3,2,3)
          !
          CFibMat(6,1) = Dtensor(2,3,1,1)
          CFibMat(6,2) = Dtensor(2,3,2,2)
          CFibMat(6,3) = Dtensor(2,3,3,3)
          CFibMat(6,4) = Dtensor(2,3,1,2)
          CFibMat(6,5) = Dtensor(2,3,1,3)
          CFibMat(6,6) = Dtensor(2,3,2,3)
          !
          END IF
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C state variable
C  
C total stiffness matrix
      DO K1 = 1,NTENS
          STATEV(K1+3*NTENS+4) = CIsoMat(K1,K1) + CFibMat(K1,K1)
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stress
C
      DO K1 = 1,NTENS
          STRESS(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          STRESS(K1) = KStressH(K1) + KStressF(K1)
      END DO
C
C jacobean
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              DDSDDE(K1,K2) = ZERO
          END DO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              DDSDDE(K1,K2) = CIsoMat(K1,K2) + CFibMat(K1,K2)
          END DO
      END DO
C      
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C end of the adhesion part
C      
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
      GO TO 888
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C chromatin (beginning of the chromatin part)
C
777   CONTINUE
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ------------------------   chromatin    ------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C bulk modulus
C
      Kmod = EmodChromatin / ( THREE * (ONE-TWO*EnuChromatin) )
C
C shear modulus
C
      Gmod = EmodChromatin / ( TWO * (ONE+EnuChromatin) )
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C   
C isotropic elastic stiffness matrix
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMat(K1,K2) = ZERO
          END DO
      END DO
      !
      CIsoMat(1,1) = Kmod + (FOUR/THREE)*Gmod 
      CIsoMat(2,2) = Kmod + (FOUR/THREE)*Gmod
      CIsoMat(3,3) = Kmod + (FOUR/THREE)*Gmod
      CIsoMat(4,4) = Gmod
      CIsoMat(1,2) = Kmod - (TWO/THREE) *Gmod
      CIsoMat(1,3) = Kmod - (TWO/THREE) *Gmod
      CIsoMat(2,3) = Kmod - (TWO/THREE) *Gmod
      !
      IF (NSHR.EQ.3) THEN
          CIsoMat(5,5) = Gmod
          CIsoMat(6,6) = Gmod
      END IF
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CIsoMat(K2,K1) = CIsoMat(K1,K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C total strain
C
      DO K1 = 1,NTENS
          KTotStrain(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          KTotStrain(K1) = STRAN(K1) + DSTRAN(K1)
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C hookean isotropic stress
C
      DO K1 = 1,NTENS
          KStressH(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              KStressH(K1) = KStressH(K1)   +
     1                       CIsoMat(K1,K2) * KTotStrain(K2)
          END DO
      END DO
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C HDAC in nucleus
C
C rectangular substrate geometry 
      HDACnucleus(1)  = 1.000000000000000
      HDACnucleus(2)  = 0.999792821534719
      HDACnucleus(3)  = 0.998882698947573
      HDACnucleus(4)  = 0.996728187613163
      HDACnucleus(5)  = 0.992719532441483
      HDACnucleus(6)  = 0.986345966606247
      HDACnucleus(7)  = 0.977054812541543
      HDACnucleus(8)  = 0.964402424465339
      HDACnucleus(9)  = 0.947909580910113
      HDACnucleus(10) = 0.927106363364673
      HDACnucleus(11) = 0.906438623761120
      HDACnucleus(12) = 0.890899187902448
      HDACnucleus(13) = 0.879349917779127
      HDACnucleus(14) = 0.871077682243541
      HDACnucleus(15) = 0.865340476655304
      HDACnucleus(16) = 0.861503597654521
      HDACnucleus(17) = 0.859044411333359
      HDACnucleus(18) = 0.857554379887084
      HDACnucleus(19) = 0.856722726613341
      HDACnucleus(20) = 0.856320844369469
      HDACnucleus(21) = 0.856120844369469
C
C circular substrate geometry 
C      HDACnucleus(1)  = 1.000000000000000
C      HDACnucleus(2)  = 1.000207088285557
C      HDACnucleus(3)  = 1.001114527267106
C      HDACnucleus(4)  = 1.003247943143750
C      HDACnucleus(5)  = 1.007159867359143
C      HDACnucleus(6)  = 1.013223265570013
C      HDACnucleus(7)  = 1.021720813231666
C      HDACnucleus(8)  = 1.032657534967920
C      HDACnucleus(9)  = 1.045876441683087
C      HDACnucleus(10) = 1.060998160771585
C      HDACnucleus(11) = 1.074159372382323
C      HDACnucleus(12) = 1.082345372558426
C      HDACnucleus(13) = 1.087291362321035
C      HDACnucleus(14) = 1.090043411808206
C      HDACnucleus(15) = 1.091413207640840
C      HDACnucleus(16) = 1.091965204093967
C      HDACnucleus(17) = 1.092067046424438
C      HDACnucleus(18) = 1.091952254108004
C      HDACnucleus(19) = 1.091753216148499
C      HDACnucleus(20) = 1.091545152417491
C      HDACnucleus(21) = 1.091645152417491
C     
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C MKL in nucleus
C
C rectangular substrate geometry
      MKLnucleus(1)   = 1.000000000000000
      MKLnucleus(2)   = 1.000541808998631
      MKLnucleus(3)   = 1.002935667258440
      MKLnucleus(4)   = 1.008636077995842
      MKLnucleus(5)   = 1.019279893547127
      MKLnucleus(6)   = 1.036183167311590
      MKLnucleus(7)   = 1.060630513028301
      MKLnucleus(8)   = 1.093361898327552
      MKLnucleus(9)   = 1.134863027985962
      MKLnucleus(10)  = 1.185114455542005
      MKLnucleus(11)  = 1.233186406730080
      MKLnucleus(12)  = 1.268684486046772
      MKLnucleus(13)  = 1.294980067037667
      MKLnucleus(14)  = 1.313903945574584
      MKLnucleus(15)  = 1.327193740835068
      MKLnucleus(16)  = 1.336248600941721
      MKLnucleus(17)  = 1.342235786241686
      MKLnucleus(18)  = 1.346037482175363
      MKLnucleus(19)  = 1.348354977395289
      MKLnucleus(20)  = 1.349685187952911
      MKLnucleus(21)  = 1.349885187952911
C
C circular substrate geometry
C      MKLnucleus(1)   = 1.000000000000000
C      MKLnucleus(2)   = 0.999458196914912
C      MKLnucleus(3)   = 0.997064925194736
C      MKLnucleus(4)   = 0.991362261402683
C      MKLnucleus(5)   = 0.980709288997536
C      MKLnucleus(6)   = 0.963757568189158
C      MKLnucleus(7)   = 0.939165712532004
C      MKLnucleus(8)   = 0.906058401158724
C      MKLnucleus(9)   = 0.863736246445220
C      MKLnucleus(10)  = 0.811890624537776 
C      MKLnucleus(11)  = 0.761212061367534
C      MKLnucleus(12)  = 0.722188474251356
C      MKLnucleus(13)  = 0.691550464549814
C      MKLnucleus(14)  = 0.667539687195569
C      MKLnucleus(15)  = 0.648628221715794
C      MKLnucleus(16)  = 0.633645922190848
C      MKLnucleus(17)  = 0.621675137981263
C      MKLnucleus(18)  = 0.612032887516333
C      MKLnucleus(19)  = 0.604182212763799
C      MKLnucleus(20)  = 0.597731058035230
C      MKLnucleus(21)  = 0.596731058035230
C     
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C state variables
C
C HDAC in nucleus (and also HDAC in cytoplasm)
      IF (KSTEP.EQ.1) THEN
          STATEV(1+4*NTENS+4) = HDACnucleus(KINC)
      ELSE
C          STATEV(1+4*NTENS+4) = HDACnucleus(KINC+10)   ! low contractility  
          STATEV(1+4*NTENS+4) = HDACnucleus(KINC+11)   ! high contractility
      END IF
C
C MKL in nucleus (and also MKL in cytoplasm)
      IF (KSTEP.EQ.1) THEN
          STATEV(2+4*NTENS+4) = MKLnucleus(KINC)
      ELSE  
C          STATEV(2+4*NTENS+4) = MKLnucleus(KINC+10)      ! low contractility 
          STATEV(2+4*NTENS+4) = MKLnucleus(KINC+11)      ! high contractility
      END IF
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C      
C stress
C
      DO K1 = 1,NTENS
          STRESS(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          STRESS(K1) = KStressH(K1)
      END DO
C
C jacobean
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              DDSDDE(K1,K2) = ZERO
          END DO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              DDSDDE(K1,K2) = CIsoMat(K1,K2) 
          END DO
      END DO 
C      
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C end of the chromatin part
C      
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
888   CONTINUE
C      
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C end of UMAT subroutine
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------- subroutines ----------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C      
C solving set of linear equations using LU factorization
C
      SUBROUTINE LU(neq,A,b,x)
C
      INCLUDE 'ABA_PARAM.INC'
C  
      INTEGER :: neq
      INTEGER :: I, J, K
C
      REAL*8                     :: Eta
C
      REAL*8, DIMENSION(neq,neq) :: A
      REAL*8, DIMENSION(neq)     :: b
      REAL*8, DIMENSION(neq)     :: x
      REAL*8, DIMENSION(neq,neq) :: M
      REAL*8, DIMENSION(neq)     :: FP
C  
      DO I = 1,neq
        DO J = 1,neq
          M(I,J) = A(I,J)
        END DO
        X(I)  = 0.0D0
        FP(I) = 0.0D0
      END DO
C
C make L+U-I
C
      DO I = 1,neq
        DO J = I+1,neq
          Eta    = M(J,I)/M(I,I)
          M(J,I) = Eta
          DO K = I+1,neq
             M(J,K) = M(J,K) - Eta * M(I,K)
          END DO
        END DO
      END DO
C
C forward substitution for f' from f and L
C
      DO J = 1,neq
        DO I = 1,J-1
          FP(J) = FP(J) - M(J,I) * FP(I)
        END DO
        FP(J)   = FP(J) + b(J)
      END DO
C
C backward substitution for u from f' and U
C
      J = neq
      DO WHILE (J>=1)
        I = neq
        DO WHILE (I>=J+1)
          X(J) = X(J) -M(J,I) * X(I)
          I    = I - 1
        END DO
        X(J) = ( X(J) + FP(J) ) / M(J,J)
        J    = J - 1
      END DO
C      
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C inverse of a matrix using LU factorization
C
      SUBROUTINE KINV( neq , A , AINV )
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER :: neq
      INTEGER :: K1 , K2 , K3
C
      REAL*8 , DIMENSION(neq,neq) :: A
      REAL*8 , DIMENSION(neq,neq) :: AINV
      REAL*8 , DIMENSION(1:neq)   :: B
      REAL*8 , DIMENSION(1:neq)   :: X
C
      DO K1 = 1,neq
          DO K2 = 1,neq
              AINV(K1,K2) = 0.0D0
          END DO
      END DO
C
      DO K1 = 1,neq
          ! set to zero
          DO K2 = 1,neq
              B(K2) = 0.0D0
          END DO
          !
          B(K1) = 1.0D0
          ! solve set of equations
          CALL LU(neq,A,B,X)
          !
          DO K3 = 1,neq
              AINV(K3,K1) = X(K3)
          END DO
          !
      END DO
C      
      RETURN
      END      
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C Kronecker delta function
C
      SUBROUTINE KronDelFunc(I,J,D)
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER I
      INTEGER J
      INTEGER D
C
      IF (I.EQ.J) THEN
          D = 1
      ELSE
          D = 0
      END IF
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C symmetric identity fourth order tensor
C
      SUBROUTINE IdentitySymmetric(IdSym)
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER Delik
      INTEGER Deljl
      INTEGER Delil
      INTEGER Deljk
C
      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER L
C
      REAL*8, DIMENSION(3,3,3,3) :: IdSym
C
      DO I = 1,3
          DO J = 1,3
              DO K = 1,3
                  DO L = 1,3
                      !
                      CALL KronDelFunc(I,K,Delik)
                      CALL KronDelFunc(J,L,Deljl)
                      CALL KronDelFunc(I,L,Delil)
                      CALL KronDelFunc(J,K,Deljk)
                      !
                      IdSym(I,J,K,L) = 0.5D0 * ( Delik*Deljl +
     1                                           Delil*Deljk )
                      !
                  END DO
              END DO
          END DO
      END DO
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C dyadic product of two vectors
C
      SUBROUTINE DyadVect(A , B , N , C)
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER K1
      INTEGER K2
C
      INTEGER N
C
      REAL*8 , DIMENSION(1:N) :: A
      REAL*8 , DIMENSION(1:N) :: B
      REAL*8 , DIMENSION(N,N) :: C
C
      DO K1 = 1,N
          DO K2 = 1,N
              C(K1,K2) = A(K1) * B(K2)
          END DO
      END DO
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C dyadic product of two matrices
C
      SUBROUTINE DyadMat(A , B , N , C)
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER N
C
      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER L
C
      REAL*8 , DIMENSION(N,N)     :: A
      REAL*8 , DIMENSION(N,N)     :: B
      REAL*8 , DIMENSION(N,N,N,N) :: C    
C
      DO I = 1,N
          DO J = 1,N
              DO K = 1,N
                  DO L = 1,N
                      C(I,J,K,L) = A(I,J) * B(K,L)
                  END DO
              END DO
          END DO
      END DO
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C first derivative of the energy function with respect to principal strain
C
      SUBROUTINE dWdPrinStrain(X , EF , StranC , MM , NN , Ratio , Y)
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (ZERO =0.0D0, HALF = 0.5D0, ONE =1.0D0, TWO=2.0D0,
     1           THREE=3.0D0, FOUR = 4.0D0, FIVE=5.0D0, SIX=6.0D0, 
     2           TOLER = 1.0E-12)
C
      REAL*8 MM
      REAL*8 NN
      REAL*8 X
      REAL*8 Y
      REAL*8 EF
      REAL*8 Ratio
      REAL*8 StranC
      REAL*8 StranT
C
      StranT = StranC * Ratio
C          
      IF (X .LT. StranC-HALF*StranT) THEN
          Y = ZERO
      ELSE IF (X .LT. StranC+HALF*StranT) THEN
          Y=EF*((X-(StranC-HALF*StranT))/StranT)**(NN+TWO)*StranT*StranT
     1         /((NN+ONE)*(NN+TWO))
      ELSE
          Y=EF*StranT*StranT/((NN+ONE)*(NN+TWO)) +
     1      EF*((ONE+X-(StranC+HALF*StranT))**(MM+TWO)-ONE)/
     2         ((MM+ONE)*(MM+TWO)) +
     3      EF*((StranC+HALF*StranT)-X)/(MM+ONE) +
     4      EF*(X-(StranC+HALF*StranT))*StranT/(NN+ONE)
      END IF
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C second derivative of the energy function with respect to principal strain
C
      SUBROUTINE d2WdPrinStrain(X , EF , StranC , MM , NN , Ratio , Y)
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (ZERO =0.0D0, HALF = 0.5D0, ONE =1.0D0, TWO=2.0D0,
     1           THREE=3.0D0, FOUR = 4.0D0, FIVE=5.0D0, SIX=6.0D0, 
     2           TOLER = 1.0E-12)
C
      REAL*8 MM
      REAL*8 NN
      REAL*8 X
      REAL*8 Y
      REAL*8 EF
      REAL*8 Ratio
      REAL*8 StranC
      REAL*8 StranT
C
      StranT = StranC * Ratio
C  
      IF (X .LT. StranC-HALF*StranT) THEN
          Y = ZERO
      ELSE IF (X .LT. StranC+HALF*StranT) THEN
          Y=EF*((X-(StranC-HALF*StranT))/StranT)**(NN+ONE)*StranT
     1      /(NN+ONE)
      ELSE
          Y=EF*StranT/(NN+ONE)+EF*((ONE+X-(StranC+HALF*StranT))**
     1                (MM+ONE)-ONE) / (MM+ONE)
      END IF
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C stiffness tensor: derivative of stress matrix with respect to strain matrix
C
      SUBROUTINE StiffTens(PrinStrainVal1,PrinStrainVal2,PrinStrainVal3,
     1  PrinStressVal1 , PrinStressVal2 , PrinStressVal3 ,    
     2  PrinEmodVal1   , PrinEmodVal2   , PrinEmodVal3   ,    
     3  EigenProj1     , EigenProj2     , EigenProj3     ,
     4  IdSym          , Imat           , Xmat           , Dtensor )
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (ZERO =0.0D0, HALF = 0.5D0, ONE =1.0D0, TWO=2.0D0,
     1           THREE=3.0D0, FOUR = 4.0D0, FIVE=5.0D0, SIX=6.0D0, 
     2           TOLER = 1.0E-12)
C
      PARAMETER (TolP = 1.0E-12)
C
      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER L
C
      INTEGER Delik
      INTEGER Delil
      INTEGER Deljl
      INTEGER Delkj
C
      REAL*8 Check12
      REAL*8 Check13
      REAL*8 Check23
      REAL*8 Fac1
      REAL*8 Fac2
      REAL*8 Fac3
      REAL*8 PrinEmodVal1
      REAL*8 PrinEmodVal2
      REAL*8 PrinEmodVal3
      REAL*8 PrinStrainVal1
      REAL*8 PrinStrainVal2
      REAL*8 PrinStrainVal3
      REAL*8 PrinStressVal1
      REAL*8 PrinStressVal2
      REAL*8 PrinStressVal3
      REAL*8 s1
      REAL*8 s2
      REAL*8 s3
      REAL*8 s4
      REAL*8 s5
      REAL*8 s6
      REAL*8 x12
      REAL*8 x13
      REAL*8 x23
      REAL*8 x21
      REAL*8 x31
      REAL*8 x32
      REAL*8 x23P
      REAL*8 x31P
      REAL*8 x12P
      REAL*8 y13
      REAL*8 y21
      REAL*8 y32
C
      REAL*8 , DIMENSION(3,3,3,3)     :: Dtensor
      REAL*8 , DIMENSION(3,3,3,3)     :: d2XdX
      REAL*8 , DIMENSION(3,3,3,3)     :: E1E1
      REAL*8 , DIMENSION(3,3,3,3)     :: E2E2
      REAL*8 , DIMENSION(3,3,3,3)     :: E3E3
      REAL*8 , DIMENSION(3,3)         :: EigenProj1
      REAL*8 , DIMENSION(3,3)         :: EigenProj2
      REAL*8 , DIMENSION(3,3)         :: EigenProj3
      REAL*8 , DIMENSION(3,3,3,3)     :: IdSym
      REAL*8 , DIMENSION(3,3,3,3)     :: II
      REAL*8 , DIMENSION(3,3)         :: Imat
      REAL*8 , DIMENSION(3,3,3,3)     :: IX
      REAL*8 , DIMENSION(3,3,3,3)     :: Term1
      REAL*8 , DIMENSION(3,3,3,3)     :: Term2
      REAL*8 , DIMENSION(3,3,3,3)     :: Term3
      REAL*8 , DIMENSION(3,3)         :: Xmat
      REAL*8 , DIMENSION(3,3,3,3)     :: XI
      REAL*8 , DIMENSION(3,3,3,3)     :: XX
C
      DO I = 1,3
          DO J = 1,3
              DO K = 1,3
                  DO L = 1,3
                      E1E1(I,J,K,L)    = ZERO
                      E2E2(I,J,K,L)    = ZERO
                      E3E3(I,J,K,L)    = ZERO
                      XX(I,J,K,L)      = ZERO
                      XI(I,J,K,L)      = ZERO
                      IX(I,J,K,L)      = ZERO
                      II(I,J,K,L)      = ZERO
                      d2XdX(I,J,K,L)   = ZERO
                      Term1(I,J,K,L)   = ZERO
                      Term2(I,J,K,L)   = ZERO
                      Term3(I,J,K,L)   = ZERO
                      Dtensor(I,J,K,L) = ZERO
                  END DO
              END DO
          END DO
      END DO
C      
      CALL DyadMat(EigenProj1 , EigenProj1 , 3 , E1E1)
      CALL DyadMat(EigenProj2 , EigenProj2 , 3 , E2E2)
      CALL DyadMat(EigenProj3 , EigenProj3 , 3 , E3E3)
C
      CALL DyadMat(Xmat , Xmat , 3 , XX)
      CALL DyadMat(Xmat , Imat , 3 , XI)
      CALL DyadMat(Imat , Xmat , 3 , IX)
      CALL DyadMat(Imat , Imat , 3 , II)
C
      DO I = 1,3
          DO J = 1,3
              DO K =1,3
                  DO L = 1,3
                      !
                      CALL KronDelFunc(I,K,Delik)
                      CALL KronDelFunc(I,L,Delil)
                      CALL KronDelFunc(J,L,Deljl)
                      CALL KronDelFunc(K,J,Delkj)
                      !
                      d2XdX(I,J,K,L) = 0.5D0 * ( Delik * Xmat(L,J) +
     1                                           Delil * Xmat(K,J) +
     2                                           Deljl * Xmat(I,K) +
     3                                           Delkj * Xmat(I,L) )
                      !
                  END DO
              END DO
          END DO
      END DO
C
      Check12 = ABS( (PrinStrainVal1-PrinStrainVal2) / PrinStrainVal1 )
      Check13 = ABS( (PrinStrainVal1-PrinStrainVal3) / PrinStrainVal1 )
      Check23 = ABS( (PrinStrainVal2-PrinStrainVal3) / PrinStrainVal2 )
C
      x12  = PrinStrainVal1 - PrinStrainVal2
      x13  = PrinStrainVal1 - PrinStrainVal3
      x23  = PrinStrainVal2 - PrinStrainVal3
      x21  = PrinStrainVal2 - PrinStrainVal1
      x31  = PrinStrainVal3 - PrinStrainVal1
      x32  = PrinStrainVal3 - PrinStrainVal2
C
      x23P = PrinStrainVal2 + PrinStrainVal3
      x31P = PrinStrainVal3 + PrinStrainVal1
      x12P = PrinStrainVal1 + PrinStrainVal2
C     
      y13  = PrinStressVal1 - PrinStressVal3
      y21  = PrinStressVal2 - PrinStressVal1
      y32  = PrinStressVal3 - PrinStressVal2
C
      Fac1 = PrinStressVal1 / ( x12 * x13 )
      Fac2 = PrinStressVal2 / ( x23 * x21 )
      Fac3 = PrinStressVal3 / ( x31 * x32 )                   
C
      IF ( (Check12.GT.TolP).AND.
     1     (Check13.GT.TolP).AND.
     2     (Check23.GT.TolP) ) THEN
          ! x1 =~ x2 =~ x3
          DO I = 1,3
              DO J = 1,3
                  DO K = 1,3
                      DO L = 1,3
                          
                          Term1(I,J,K,L) = Fac1 * ( d2XdX(I,J,K,L) -
     1                    x23P      * IdSym(I,J,K,L) - 
     2                    (x12+x13) *  E1E1(I,J,K,L) - 
     2                    x23       * (E2E2(I,J,K,L)-E3E3(I,J,K,L)) ) +
     3                    PrinEmodVal1 * E1E1(I,J,K,L)
                          
                          Term2(I,J,K,L) = Fac2 * ( d2XdX(I,J,K,L) -
     1                    x31P      * IdSym(I,J,K,L) - 
     2                    (x23+x21) *  E2E2(I,J,K,L) - 
     2                    x31       * (E3E3(I,J,K,L)-E1E1(I,J,K,L)) ) +
     3                    PrinEmodVal2 * E2E2(I,J,K,L)
                          
                          Term3(I,J,K,L) = Fac3 * ( d2XdX(I,J,K,L) -
     1                    x12P      * IdSym(I,J,K,L) - 
     2                    (x31+x32) *  E3E3(I,J,K,L) - 
     2                    x12       * (E1E1(I,J,K,L)-E2E2(I,J,K,L)) ) +
     3                    PrinEmodVal3 * E3E3(I,J,K,L)
                          
                          Dtensor(I,J,K,L) = Term1(I,J,K,L) +
     1                                       Term2(I,J,K,L) +
     2                                       Term3(I,J,K,L)
                          
                      END DO
                  END DO
              END DO
         END DO
      !
      ELSE IF ( (Check12.GT.TolP).AND.
     1          (Check13.GT.TolP).AND.
     2          (Check23.LT.TolP) ) THEN
               ! x1 =~ x2 = x3
               s1 = ( y13                         / (x13**TWO)   ) -
     1              (  PrinEmodVal3               /  x13         )
               s2 = ( y13 * TWO*PrinStrainVal3    / (x13**TWO)   ) -
     1              (  PrinEmodVal3 * x31P        /  x13         )
               s3 = ( y13 * TWO                   / (x13**THREE) ) -
     1              ( (PrinEmodVal1+PrinEmodVal3) / (x13**TWO)   )
               s4 = s3 * PrinStrainVal3
               s5 = s3 * PrinStrainVal3
               s6 = s3 * PrinStrainVal3 * PrinStrainVal3
               DO I = 1,3
                   DO J = 1,3
                       DO K = 1,3
                           DO L = 1,3
                               Dtensor(I,J,K,L) = s1 * d2XdX(I,J,K,L) -
     1                                            s2 * IdSym(I,J,K,L) -
     2                                            s3 * XX(I,J,K,L)    +
     3                                            s4 * XI(I,J,K,L)    +
     4                                            s5 * IX(I,J,K,L)    -
     5                                            s6 * II(I,J,K,L)
                           END DO
                       END DO
                   END DO
               END DO
      !
      ELSE IF ( (Check12.GT.TolP).AND.
     1          (Check13.LT.TolP).AND.
     2          (Check23.GT.TolP) ) THEN
               ! x2 =~ x3 = x1
               s1 = ( y21                         / (x21**TWO)   ) -
     1              (  PrinEmodVal1               /  x21         )
               s2 = ( y21 * TWO*PrinStrainVal1    / (x21**TWO)   ) -
     1              (  PrinEmodVal1 * x12P        /  x21         )
               s3 = ( y21 * TWO                   / (x21**THREE) ) -
     1              ( (PrinEmodVal2+PrinEmodVal1) / (x21**TWO)   )
               s4 = s3 * PrinStrainVal1
               s5 = s3 * PrinStrainVal1
               s6 = s3 * PrinStrainVal1 * PrinStrainVal1
               DO I = 1,3
                   DO J = 1,3
                       DO K = 1,3
                           DO L = 1,3
                               Dtensor(I,J,K,L) = s1 * d2XdX(I,J,K,L) -
     1                                            s2 * IdSym(I,J,K,L) -
     2                                            s3 * XX(I,J,K,L)    +
     3                                            s4 * XI(I,J,K,L)    +
     4                                            s5 * IX(I,J,K,L)    -
     5                                            s6 * II(I,J,K,L)
                           END DO
                       END DO
                   END DO
               END DO
      !
      ELSE IF ( (Check12.LT.TolP).AND.
     1          (Check13.GT.TolP).AND.
     2          (Check23.GT.TolP) ) THEN
               ! x3 =~ x1 = x2
               s1 = ( y32                         / (x32**TWO)   ) -
     1              (  PrinEmodVal2               /  x32         )
               s2 = ( y32 * TWO*PrinStrainVal2    / (x32**TWO)   ) -
     1              (  PrinEmodVal1 * x23P        /  x32         )
               s3 = ( y32 * TWO                   / (x32**THREE) ) -
     1              ( (PrinEmodVal3+PrinEmodVal2) / (x32**TWO)   )
               s4 = s3 * PrinStrainVal2
               s5 = s3 * PrinStrainVal2
               s6 = s3 * PrinStrainVal2 * PrinStrainVal2
               DO I = 1,3
                   DO J = 1,3
                       DO K = 1,3
                           DO L = 1,3
                               Dtensor(I,J,K,L) = s1 * d2XdX(I,J,K,L) -
     1                                            s2 * IdSym(I,J,K,L) -
     2                                            s3 * XX(I,J,K,L)    +
     3                                            s4 * XI(I,J,K,L)    +
     4                                            s5 * IX(I,J,K,L)    -
     5                                            s6 * II(I,J,K,L)
                           END DO
                       END DO
                   END DO
               END DO
      !
      ELSE IF ( (Check12.LT.TolP).AND.
     1          (Check13.LT.TolP).AND.
     2          (Check23.LT.TolP) ) THEN
               ! x1 = x2 = x3
               DO I = 1,3
                   DO J = 1,3
                       DO K = 1,3
                           DO L = 1,3
                               Dtensor(I,J,K,L) = PrinEmodVal1 *
     1                                            IdSym(I,J,K,L)
                           END DO
                       END DO
                   END DO
               END DO
      !
      END IF
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C 
C newton raphson method to solve set of nonlinear equations of cell
C
      SUBROUTINE KResidual(CMatC1 , CIsoMatC2 , KStran , KTotStran ,
     1 IdSym , Imat , KEyeMat ,
     2 NDI, NSHR, NTENS , KINC , NOEL , NPT ,
     3 EmodCF , StranCCell , RatioCell , MMCell , NNCell , Rho0Bar ,
     4 KStressC1 , KStressC2 , CMatC2 , KJacob , KRes ,
     5 PrinEmodVal1C , PrinEmodVal2C , PrinEmodVal3C)
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (ZERO =0.0D0, HALF = 0.5D0, ONE =1.0D0, TWO=2.0D0,
     1           THREE=3.0D0, FOUR = 4.0D0, FIVE=5.0D0, SIX=6.0D0, 
     2           TOLER = 1.0E-12)
C
      INTEGER K1
      INTEGER K2 
C
      REAL*8 EmodCF
      REAL*8 MMCell
      REAL*8 NNCell
      REAL*8 PrinEmodVal1C
      REAL*8 PrinEmodVal2C
      REAL*8 PrinEmodVal3C
      REAL*8 PrinStrainVal1C
      REAL*8 PrinStrainVal2C
      REAL*8 PrinStrainVal3C
      REAL*8 PrinStressVal1C
      REAL*8 PrinStressVal2C
      REAL*8 PrinStressVal3C
      REAL*8 RatioCell
      REAL*8 Rho0Bar
      REAL*8 StranCCell
C
      REAL*8 , DIMENSION(3,3,3,3)         :: DtensorC
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CIsoMatC2
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CFibMatC2
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CMatC1
      REAL*8 , DIMENSION(NTENS,NTENS)     :: CMatC2
      REAL*8 , DIMENSION(3,3)             :: EigenProj1C
      REAL*8 , DIMENSION(3,3)             :: EigenProj2C
      REAL*8 , DIMENSION(3,3)             :: EigenProj3C
      REAL*8 , DIMENSION(3,3,3,3)         :: IdSym
      REAL*8 , DIMENSION(3,3)             :: Imat
      REAL*8 , DIMENSION(NTENS,NTENS)     :: KEyeMat
      REAL*8 , DIMENSION(2*NTENS,2*NTENS) :: KJacob
      REAL*8 , DIMENSION(1:2*NTENS)       :: KRes
      REAL*8 , DIMENSION(1:2*NTENS)       :: KStran
      REAL*8 , DIMENSION(1:NTENS)         :: KStressC1
      REAL*8 , DIMENSION(1:NTENS)         :: KStressC2
      REAL*8 , DIMENSION(1:NTENS)         :: KStressC2H
      REAL*8 , DIMENSION(1:NTENS)         :: KStressC2F
      REAL*8 , DIMENSION(1:NTENS)         :: KTotStran
      REAL*8 , DIMENSION(1:3)             :: KValC
      REAL*8 , DIMENSION(3,3)             :: KVecC
      REAL*8 , DIMENSION(1:3)             :: PrinStrainVec1C
      REAL*8 , DIMENSION(1:3)             :: PrinStrainVec2C
      REAL*8 , DIMENSION(1:3)             :: PrinStrainVec3C
      REAL*8 , DIMENSION(3,3)             :: StranMat
      REAL*8 , DIMENSION(3,3)             :: StrssMat
C
C
C
C total stress (first component)
C
      DO K1 = 1,NTENS
          KStressC1(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              KStressC1(K1) = KStressC1(K1)   +
     1                        CMatC1(K1,K2)   * KStran(K2)
          END DO
      END DO
      !
      KStressC1(1) = KStressC1(1) + Rho0Bar
      KStressC1(2) = KStressC1(2) + Rho0Bar
      KStressC1(3) = KStressC1(3) + Rho0Bar
C
C
C
C hookean isotropic stress (second component)
C
      DO K1 = 1,NTENS
          KStressC2H(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              KStressC2H(K1) = KStressC2H(K1)   +
     1                         CIsoMatC2(K1,K2) * KStran(K2+NTENS)
          END DO
      END DO
C
C
C
C strain matrix (second component)
C
      DO K1 = 1,3
          DO K2 = 1,3
              StranMat(K1,K2) = ZERO
          END DO
      END DO
      !        
      StranMat(1,1) = KStran(1+NTENS)
      StranMat(2,2) = KStran(2+NTENS)
      StranMat(3,3) = KStran(3+NTENS)
      StranMat(1,2) = KStran(4+NTENS) * HALF
      !
      IF (NSHR.EQ.3) THEN
          StranMat(1,3) = KStran(5+NTENS) * HALF
          StranMat(2,3) = KStran(6+NTENS) * HALF
      END IF
      !
      DO K1 = 1,3
          DO K2 = 1,3
              StranMat(K2,K1) = StranMat(K1,K2)
          END DO
      END DO
C
C
C
C eigen values and eigen vectors of strain (second component) 
C    
      CALL SPRIND( KStran(1+NTENS:2*NTENS), KValC, KVecC, 2, NDI, NSHR )
C
C
C
C principal strain (second component)
C
      PrinStrainVal1C = KValC(1)
      PrinStrainVal2C = KValC(2)
      PrinStrainVal3C = KValC(3)
C
C
C
C principal vectors (second component)
C
      PrinStrainVec1C(1) = KVecC(1,1)
      PrinStrainVec1C(2) = KVecC(1,2)
      PrinStrainVec1C(3) = KVecC(1,3)
      !
      PrinStrainVec2C(1) = KVecC(2,1)
      PrinStrainVec2C(2) = KVecC(2,2)
      PrinStrainVec2C(3) = KVecC(2,3)
      !
      PrinStrainVec3C(1) = KVecC(3,1)
      PrinStrainVec3C(2) = KVecC(3,2)
      PrinStrainVec3C(3) = KVecC(3,3)  
C   
C
C
C eigen projection (second component)
C
      CALL DyadVect(PrinStrainVec1C , PrinStrainVec1C , 3 , EigenProj1C)
      CALL DyadVect(PrinStrainVec2C , PrinStrainVec2C , 3 , EigenProj2C)
      CALL DyadVect(PrinStrainVec3C , PrinStrainVec3C , 3 , EigenProj3C)
C
C
C
C principal stress (second component)
C     
C      EmodCF1 = EmodCF
C      IF((PrinStrainVal1C-PrinStrainVal2C-PrinStrainVal3C).GT.ZERO) THEN
C          EmodCF1 = EmodCF + EmodCF * ( PrinStrainVal1C -
C     1                                  PrinStrainVal2C -
C     2                                  PrinStrainVal3C ) * 1.75D0
C      END IF
      
C      EmodCF2 = EmodCF
C      IF((PrinStrainVal2C-PrinStrainVal1C-PrinStrainVal3C).GT.ZERO) THEN
C          EmodCF2 = EmodCF + EmodCF * ( PrinStrainVal2C -
C     1                                  PrinStrainVal1C -
C     2                                  PrinStrainVal3C ) * 1.75D0
C      END IF
      
C      EmodCF3 = EmodCF
C      IF((PrinStrainVal3C-PrinStrainVal1C-PrinStrainVal2C).GT.ZERO) THEN
C          EmodCF3 = EmodCF + EmodCF * ( PrinStrainVal3C -
C     1                                  PrinStrainVal1C -
C     2                                  PrinStrainVal2C ) * 1.75D0
C      END IF

           
      CALL dWdPrinStrain(PrinStrainVal1C , EmodCF    , StranCCell ,  
     1                   MMCell , NNCell , RatioCell , PrinStressVal1C)
C
      CALL dWdPrinStrain(PrinStrainVal2C , EmodCF    , StranCCell ,  
     1                   MMCell , NNCell , RatioCell , PrinStressVal2C)
C
      CALL dWdPrinStrain(PrinStrainVal3C , EmodCF    , StranCCell ,  
     1                   MMCell , NNCell , RatioCell , PrinStressVal3C)
C  
C
C
C derivative of principal stress with respect to principal strain (second component)
C     
      CALL d2WdPrinStrain(PrinStrainVal1C , EmodCF    , StranCCell ,  
     1                    MMCell , NNCell , RatioCell , PrinEmodVal1C)
C
      CALL d2WdPrinStrain(PrinStrainVal2C , EmodCF    , StranCCell ,  
     1                    MMCell , NNCell , RatioCell , PrinEmodVal2C)
C
      CALL d2WdPrinStrain(PrinStrainVal3C , EmodCF    , StranCCell ,  
     1                    MMCell , NNCell , RatioCell , PrinEmodVal3C)
C
C stress matrix (second component)
C
      DO K1 = 1,3
          DO K2 = 1,3
              StrssMat(K1,K2) = ZERO
          END DO
      END DO
      !
      DO K1 = 1,3
          DO K2 = 1,3
              StrssMat(K1,K2) = PrinStressVal1C * EigenProj1C(K1,K2) +
     1                          PrinStressVal2C * EigenProj2C(K1,K2) +
     2                          PrinStressVal3C * EigenProj3C(K1,K2)
          END DO
      END DO
C
C
C
C fiber stress (second component)
C
      DO K1 = 1,NTENS
          KStressC2F(K1) = ZERO
      END DO
      !    
      KStressC2F(1) = StrssMat(1,1)
      KStressC2F(2) = StrssMat(2,2)
      KStressC2F(3) = StrssMat(3,3)
      KStressC2F(4) = StrssMat(1,2)
      !
      IF (NSHR.EQ.3) THEN
          KStressC2F(5) = StrssMat(1,3)
          KStressC2F(6) = StrssMat(2,3)
      END IF
C
C
C
C stiffness tensor (second component)
C
      CALL StiffTens(PrinStrainVal1C,PrinStrainVal2C,PrinStrainVal3C,
     1  PrinStressVal1C , PrinStressVal2C , PrinStressVal3C ,    
     2  PrinEmodVal1C   , PrinEmodVal2C   , PrinEmodVal3C   ,    
     3  EigenProj1C     , EigenProj2C     , EigenProj3C     ,
     4  IdSym           , Imat            , StranMat        , DtensorC )
C
C
C
C stiffness matrix (second component)
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CFibMatC2(K1,K2) = ZERO
          END DO
      END DO
      !
      CFibMatC2(1,1) = DtensorC(1,1,1,1)
      CFibMatC2(1,2) = DtensorC(1,1,2,2)
      CFibMatC2(1,3) = DtensorC(1,1,3,3)
      CFibMatC2(1,4) = DtensorC(1,1,1,2)
      !
      CFibMatC2(2,1) = DtensorC(2,2,1,1)
      CFibMatC2(2,2) = DtensorC(2,2,2,2)
      CFibMatC2(2,3) = DtensorC(2,2,3,3)
      CFibMatC2(2,4) = DtensorC(2,2,1,2)
      !
      CFibMatC2(3,1) = DtensorC(3,3,1,1)
      CFibMatC2(3,2) = DtensorC(3,3,2,2)
      CFibMatC2(3,3) = DtensorC(3,3,3,3)
      CFibMatC2(3,4) = DtensorC(3,3,1,2)
      !
      CFibMatC2(4,1) = DtensorC(1,2,1,1)
      CFibMatC2(4,2) = DtensorC(1,2,2,2)
      CFibMatC2(4,3) = DtensorC(1,2,3,3)
      CFibMatC2(4,4) = DtensorC(1,2,1,2)
      !
      IF (NSHR.EQ.3) THEN
          !
          CFibMatC2(1,5) = DtensorC(1,1,1,3)
          CFibMatC2(1,6) = DtensorC(1,1,2,3)
          !
          CFibMatC2(2,5) = DtensorC(2,2,1,3)
          CFibMatC2(2,6) = DtensorC(2,2,2,3)
          !
          CFibMatC2(3,5) = DtensorC(3,3,1,3)
          CFibMatC2(3,6) = DtensorC(3,3,2,3)
          !
          CFibMatC2(4,5) = DtensorC(1,2,1,3)
          CFibMatC2(4,6) = DtensorC(1,2,2,3)
          !
          CFibMatC2(5,1) = DtensorC(1,3,1,1)
          CFibMatC2(5,2) = DtensorC(1,3,2,2)
          CFibMatC2(5,3) = DtensorC(1,3,3,3)
          CFibMatC2(5,4) = DtensorC(1,3,1,2)
          CFibMatC2(5,5) = DtensorC(1,3,1,3)
          CFibMatC2(5,6) = DtensorC(1,3,2,3)
          !
          CFibMatC2(6,1) = DtensorC(2,3,1,1)
          CFibMatC2(6,2) = DtensorC(2,3,2,2)
          CFibMatC2(6,3) = DtensorC(2,3,3,3)
          CFibMatC2(6,4) = DtensorC(2,3,1,2)
          CFibMatC2(6,5) = DtensorC(2,3,1,3)
          CFibMatC2(6,6) = DtensorC(2,3,2,3)
          !
          END IF
C
C
C
C total elastic stiffness (second component)
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CMatC2(K1,K2) = ZERO
          END DO
      END DO
C
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              CMatC2(K1,K2) = CIsoMatC2(K1,K2) + CFibMatC2(K1,K2)
          END DO
      END DO
C
C
C
C total stress (second component)
C
      DO K1 = 1,NTENS
          KStressC2(K1) = ZERO
      END DO
      !
      DO K1 = 1,NTENS
          KStressC2(K1) = KStressC2H(K1) + KStressC2F(K1)
      END DO
C
C
C
C jacobian matrix for nonlinear equations
C
      DO K1 = 1,2*NTENS
          DO K2 = 1,2*NTENS
              KJacob(K1,K2) = ZERO
          END DO
      END DO
      !
      DO K1 = 1,NTENS
          DO K2 = 1,NTENS
              !
              ! (1:NTENS,1:NTENS)
              KJacob(K1,K2)             =  CMatC1(K1,K2)
              !
              ! (1:NTENS,NTENS+1:2*NTENS)
              KJacob(K1,NTENS+K2)       = -CMatC2(K1,K2)
              !
              ! (NTENS+1:2*NTENS,1:NTENS)
              KJacob(NTENS+K1,K2)       = -KEyeMat(K1,K2)
              !
              ! (NTENS+1:2*NTENS,NTENS+1:2*NTENS)
              KJacob(NTENS+K1,NTENS+K2) = -KEyeMat(K1,K2)
              !
          END DO
      END DO
C
C
C
C set of nonlinear equations
C
      DO K1 = 1,2*NTENS
          KRes(K1) = ZERO
      END DO
      !
      IF (NSHR.EQ.1) THEN
          !
          KRes(1)   = KStressC1(1) - KStressC2(1)
          KRes(2)   = KStressC1(2) - KStressC2(2)
          KRes(3)   = KStressC1(3) - KStressC2(3)
          KRes(4)   = KStressC1(4) - KStressC2(4)
          !
          KRes(5)   = KTotStran(1) - KStran(1) - KStran(5)
          KRes(6)   = KTotStran(2) - KStran(2) - KStran(6)
          KRes(7)   = KTotStran(3) - KStran(3) - KStran(7)
          KRes(8)   = KTotStran(4) - KStran(4) - KStran(8)
          !
      ELSE
          !
          KRes(1)   = KStressC1(1) - KStressC2(1)
          KRes(2)   = KStressC1(2) - KStressC2(2)
          KRes(3)   = KStressC1(3) - KStressC2(3)
          KRes(4)   = KStressC1(4) - KStressC2(4)
          KRes(5)   = KStressC1(5) - KStressC2(5)
          KRes(6)   = KStressC1(6) - KStressC2(6)
          !
          KRes(7)   = KTotStran(1) - KStran(1) - KStran(7)
          KRes(8)   = KTotStran(2) - KStran(2) - KStran(8)
          KRes(9)   = KTotStran(3) - KStran(3) - KStran(9)
          KRes(10)  = KTotStran(4) - KStran(4) - KStran(10)
          KRes(11)  = KTotStran(5) - KStran(5) - KStran(11)
          KRes(12)  = KTotStran(6) - KStran(6) - KStran(12)
          !
      END IF
C     
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C first derivative of the energy function with respect to principal strain
C
      SUBROUTINE dWdPrinStrainTension(X , EF , StranC , MM , NN , 
     1                                    Ratio , Y)
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (ZERO =0.0D0, HALF = 0.5D0, ONE =1.0D0, TWO=2.0D0,
     1           THREE=3.0D0, FOUR = 4.0D0, FIVE=5.0D0, SIX=6.0D0, 
     2           TOLER = 1.0E-12)
C
      REAL*8 MM
      REAL*8 NN
      REAL*8 X
      REAL*8 XX
      REAL*8 Y
      REAL*8 EF
      REAL*8 Ratio
      REAL*8 StranC
      REAL*8 StranT
C
      StranT = StranC * Ratio
C        
C      X = ABS(XX)
C
      IF (X .LT. StranC-HALF*StranT) THEN
          Y = ZERO
      ELSE IF (X .LT. StranC+HALF*StranT) THEN
          Y=EF*((X-(StranC-HALF*StranT))/StranT)**(NN+TWO)*StranT*
     1          StranT/((NN+ONE)*(NN+TWO))
      ELSE
          Y=EF*StranT*StranT/((NN+ONE)*(NN+TWO)) +
     1      EF*((ONE+X-(StranC+HALF*StranT))**(MM+TWO)-ONE)/
     2         ((MM+ONE)*(MM+TWO)) +
     3      EF*((StranC+HALF*StranT)-X)/(MM+ONE) +
     4      EF*(X-(StranC+HALF*StranT))*StranT/(NN+ONE)
      END IF
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C
C second derivative of the energy function with respect to principal strain
C
      SUBROUTINE d2WdPrinStrainTension(X , EF , StranC , MM , NN , 
     1                                     Ratio , Y)
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (ZERO =0.0D0, HALF = 0.5D0, ONE =1.0D0, TWO=2.0D0,
     1           THREE=3.0D0, FOUR = 4.0D0, FIVE=5.0D0, SIX=6.0D0, 
     2           TOLER = 1.0E-12)
C
      REAL*8 MM
      REAL*8 NN
      REAL*8 X
      REAL*8 XX
      REAL*8 Y
      REAL*8 EF
      REAL*8 Ratio
      REAL*8 StranC
      REAL*8 StranT
C
      StranT = StranC * Ratio
C  
C      X = ABS(XX)
C
      IF (X .LT. StranC-HALF*StranT) THEN
          Y = ZERO
      ELSE IF (X .LT. StranC+HALF*StranT) THEN
          Y=EF*((X-(StranC-HALF*StranT))/StranT)**(NN+ONE)*StranT
     1      /(NN+ONE)
      ELSE
          Y=EF*StranT/(NN+ONE)+EF*((ONE+X-(StranC+HALF*StranT))**
     1                (MM+ONE)-ONE) / (MM+ONE)
      END IF
C
      RETURN
      END
C
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
      
      

      
              
      
              
      