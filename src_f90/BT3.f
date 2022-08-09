      SUBROUTINE BT3( SIZE, NLYR, UMU, AZIM, LAMBERT, FISOT,
     $          BTEMP2, EMIS2, BTEMP, LYRTEMP,LYRHGT, KEXT,   
     $          SALB, ASYM, EMIS, STEP, NXCL, NYCL,FILEOUT,
     $          HENYEY, NANG, ANG, PHADI, EMI, EMI2, NEMI,
     $          MEANFP,B0,B1,MXLYR,MXCL,MXAN,MXEMI,NPHO, tsub) 
C
C     Computes the microwave brightness temperature for 
C     horizontally and vertically structured clouds
C     using a backward Monte Carlo method.
C
C
      LOGICAL    LAMBERT, HENYEY 
      REAL  LYRTEMP(0:MXCL,0:MXCL,0:MXLYR),EMIS,EMIS2,AZIM, 
     $     EMI(0:MXEMI), EMI2(0:MXEMI),
     $     LYRHGT(0:MXCL,0:MXCL,0:MXLYR), KEXT(0:MXCL,0:MXCL,MXLYR),
     $     SALB(0:MXCL,0:MXCL,MXLYR), ASYM(0:MXCL,0:MXCL,MXLYR),  
     $     B0(0:MXCL,0:MXCL,MXLYR), B1(0:MXCL,0:MXCL,MXLYR) 
      REAL  MEANFP(0:MXCL,0:MXCL,MXLYR), SINTH, COSTH ,TANTH, 
     $     POS(3),COSS(3),POSBOUND(3),DICO,TB,DIBOMIN,
     $     PHADI(0:MXCL,0:MXCL,MXLYR,0:MXAN),ANG(0:MXAN)
      REAL BTEMP(0:MXCL,0:MXCL)
      REAL TSUB,TSUBUF
      INTEGER EXIT, ISEED,LLL,MMM,NLYR(0:MXCL,0:MXCL)	
      REAL  X, ran1, R
      CHARACTER*50 FILEOUT
      PI = 4.*ATAN(1.)
      ISEED = 3   
C     Initializing the random generating function
      X = ran1(ISEED)
C     
C     
!      OPEN(UNIT=64,FILE=FILEOUT,STATUS='UNKNOWN')
C
C
C                     Calculate some frequently used constants
C
      NX=0
      NY=0
      DO 19  J = 1,NLYR(NX,NY)
       B0(NX,NY,J)  = LYRTEMP(NX,NY,J-1)
       B1(NX,NY,J)  = (LYRTEMP(NX,NY,J) - LYRTEMP(NX,NY,J-1))
     $      /(LYRHGT(NX,NY,J) - LYRHGT(NX,NY,J-1))
       IF (KEXT(NX,NY,J).NE.0) THEN
          MEANFP(NX,NY,J) = 1./KEXT(NX,NY,J)
       ELSE
          MEANFP(NX,NY,J) = 10000.
       ENDIF
 19   CONTINUE
      DO 111 NX = 1, NXCL
         DO 112 NY = 1, NYCL
            DO 20  J = 1,NLYR(NX,NY)
               B0(NX,NY,J)  = LYRTEMP(NX,NY,J-1)
               B1(NX,NY,J)  = (LYRTEMP(NX,NY,J) - LYRTEMP(NX,NY,J-1))
     $              /(LYRHGT(NX,NY,J) - LYRHGT(NX,NY,J-1))
               IF (KEXT(NX,NY,J).NE.0) THEN
                  MEANFP(NX,NY,J) = 1./KEXT(NX,NY,J)
               ELSE
                  MEANFP(NX,NY,J) = 10000.
               ENDIF
 20         CONTINUE
 112     CONTINUE
 111  CONTINUE
C
C
C     Definition of the source parameters
C     
      azim=0
      PI    = 4.*ATAN(1.)
      AZIM  = AZIM *PI/180.
      SINTH = SQRT(1. - UMU**2)	
      TANTH = SINTH/UMU
      NSTEP = NINT(SIZE/STEP)   ! Max number of grid blocks in y
c      write(*,*) azim,sinth,tanth,nstep
c      stop
      DO 121 LLL = 1, NXCL
         DO 122 MMM = 1, NYCL
            TSUB=0.
            DO 33 I = 1,NSTEP
               DO 44 J = 1,NSTEP 
                  TSUBUF=0.
                  DO 300 N =1,NPHO	
C
C          Position and direction of the starting photons
C
        NX = LLL
        NY = MMM
        POS(3) = LYRHGT(NX,NY,NLYR(NX,NY))
        POS(2) = FLOAT(MMM-1) * SIZE + (FLOAT(J)*2. - 1.)  * STEP/2.
        POS(1) = FLOAT(LLL-1) * SIZE + (FLOAT(I)*2. - 1.)  * STEP/2.
        COSS(3) = -UMU
        COSS(2) = SINTH*SIN(PI+AZIM) 
        COSS(1) = SINTH*COS(PI+AZIM)
        NZ = NLYR(NX,NY)
C
C      Computation of the distance to collision
C
c        write(*,*) nx, ny, nlyr(nx,ny)
  110  CALL  DICOLLI(NX,NY, MEANFP,DICO, NZ,MXLYR,MXCL,ISEED)
c       write(*,*) 'DICO=', dico
c       stop
C
C      Computation of the distance to the nearest cell boundary
C
  120  CALL  DISBOUNDS(NX,NY, NZ, LYRHGT, POS, DIBOMIN,COSS, 
     $ EXIT,POSBOUND, SIZE,NEWNX,NEWNY,NEWNZ,NXCL,NYCL,MXLYR,MXCL,ISEED)
        IF (DIBOMIN.GT.DICO) THEN
C
C    If the distance to collision is shorter than the distance to the nearest 
C    boundary the photon is advanced to the point of collision and the nature
C    of the collision (scattering or absorption) is determined
C 
         R = ran1(ISEED)
         POS(3) = POS(3) + DICO * COSS(3)
         POS(2) = POS(2) + DICO * COSS(2)
         POS(1) = POS(1) + DICO * COSS(1)
          IF (R.LT.SALB(NX,NY,NZ)) THEN
C
C    In case of a scattering event new direction parameters are determined
C
           CALL SCATCLOUDS(NX,NY,COSS,ASYM,NZ,HENYEY,PHADI,ANG,
     $           	   NANG,MXLYR,MXCL,MXAN,ISEED)
           GOTO 110
           ELSE
C
C   If the photon is absorbed then is considered as being emitted in the 
C   point of absorption with the physical temperature in that point
C
CJLH added next 3 lines
           IF (NZ.GE.(NLYR(NX,NY)+1)) THEN
             TB =  FISOT
           ELSE
             CALL ABSORBT(NX,NY,NZ,POS,LYRHGT,B0,B1,TB,MXLYR,MXCL,ISEED)
           ENDIF
           GOTO 299
          ENDIF
	ELSE
C
C    If the distance to collision is longer than the distance to the nearest 
C    boundary the photon is advanced to the boundary and a new distance 
C    of collision is determined taking into account the distance already
C    travelled
C
       POS(EXIT) = POSBOUND(EXIT)
       DO 199 M=1,3
          IF (M.NE.EXIT) THEN
             POS(M) = POS(M) + DIBOMIN * COSS(M)
          ENDIF 
 199   CONTINUE
       IF (EXIT.NE.3) THEN
       ENDIF
       IF (EXIT.EQ.3) THEN
          NXOLD = NX
          NYOLD = NY
          NZOLD = NZ
          NZ = NEWNZ
        IF ((NZ.GT.0).AND.(NZ.LT.(NLYR(NX,NY)+1))) THEN
C
C   The photon crosses the interface betwee two layers
C   and a new dinstance to collision is determined
C
        CALL NEWDICOLLI(NX,NY,NZ,NXOLD,NYOLD,NZOLD,DICO,
     $              MEANFP,  DIBOMIN, MXLYR, MXCL,ISEED)
         GOTO 120
         ELSE
C
C    If the boundary is the bottom of a subloud  the 
C    interaction with the surface is considered.
C
            IF (NZ.EQ.0) THEN
            R = ran1(ISEED)
	    IF ((NX.EQ.0).AND.(NY.EQ.0))THEN
C
C    The surface emissivity is determined
C
	    IF (NEMI.NE.1) THEN
	    REMIS=EMIS2
	    ELSE 
	    IF (-COSS(3).GE.(COS(89.*PI/180.))) THEN
	    REMIS=EMI2(NINT(ACOS(-COSS(3))*180./PI))
	    ELSE
	    REMIS=emi2(89)
	    ENDIF	
	    ENDIF
	    ELSE
	    IF ((NX.NE.0).AND.(NY.NE.0)) THEN
	    IF (NEMI.NE.1) THEN
	    REMIS=EMIS
	    ELSE
	    IF (-COSS(3).GE.(COS(89.*PI/180.))) THEN
	    REMIS=EMI(NINT(ACOS(-COSS(3))*180./PI))
	    ELSE
	    REMIS=EMI(89)
	    ENDIF 
	    ENDIF
	    ENDIF
	    ENDIF
C
C    The type of interaction is determined
C
            IF (R.GT.REMIS) THEN
C    The photon is scattered and a new direction is determined
            CALL SCATGROUND(LAMBERT, COSS, ISEED)
	    NZ = NZ + 1
            GOTO 110
            ELSE
C    The photon is absorbed 
	    IF ((NX.NE.0).AND.(NY.NE.0)) THEN
            TB = BTEMP(NX,NY)
	    ELSE
	    TB = BTEMP(NX,NY)
	    ENDIF
            GOTO 299
            ENDIF
           ELSE
C
C    If the boundary is the top of a subcloud  the 
C    photon escapes the subcloud and is considered as being emitted
C    with the cosmic background temperature.
C
CJLH        IF (NZ.EQ.(NLYR(NX,NY)+1)) THEN
CJLH added next line
            IF (NZ.GE.(NLYR(NX,NY)+1)) THEN
            TB =  FISOT
            GOTO 299
            ENDIF
         ENDIF
      ENDIF
      ELSE
         IF (DIBOMIN.EQ.(10000.) ) THEN
            TB = 0.
            GOTO 299
         ENDIF
C     The photons crosses the vertical boundaries between two subclouds
C     and a new distance to collision is determined
         NXOLD = NX
         NYOLD = NY
         NZOLD = NZ 
         NX = NEWNX
         NY = NEWNY
      CALL GETNEWNZ(NX,NY,NEWNZ,POS(3),COSS(3),NLYR,LYRHGT,MXLYR,MXCL,
     $        ISEED)
      NZ = NEWNZ
      IF (NZOLD.NE.NEWNZ) THEN
      ENDIF
      CALL NEWDICOLLI(NX,NY,NZ,NXOLD,NYOLD,NZOLD,DICO,
     $        MEANFP,DIBOMIN,MXLYR,MXCL,ISEED)
      GOTO 120
      ENDIF
      ENDIF
 299  CONTINUE
      TSUBUF=TSUBUF+TB
      !write(*,*) tb, nx, ny, nz
 300  CONTINUE
      TSUBUF=TSUBUF / FLOAT(NPHO)
      TSUB=TSUB+TSUBUF
 44   CONTINUE
 33   CONTINUE
      TSUB=TSUB/NSTEP/NSTEP
C     
C     The brightness temperature of each subcloud is written in a file
C     
      !WRITE(64,*)LLL,MMM,TSUB
      !WRITE(*,*)LLL,MMM,TSUB, nstep
 122  CONTINUE
 121  CONTINUE
      RETURN
      END
C
c23456
      SUBROUTINE GETNEWNZ(NX,NY,NEWNZ,POS3,COS3,NLYR,LYRHGT,MXLYR,MXCL,
     $     ISEED)
      INTEGER NLYR(0:MXCL,0:MXCL)
      REAL LYRHGT (0:MXCL,0:MXCL,0:MXLYR),POS3,COS3
        IF ((COS3.GE.0.).AND.(POS3.EQ.0.)) THEN
        NEWNZ=1
        GOTO 1
        ENDIF
        DO 118 I=1,NLYR(NX,NY)
        IF (POS3.EQ.LYRHGT(NX,NY,I)) THEN
        IF (COS3.GT.0.) THEN
        NEWNZ=I+1
        ELSE
        NEWNZ=I
        ENDIF
        GOTO 1
        ELSE
        IF ((POS3.GT.LYRHGT(NX,NY,I-1)).AND.
     $  (POS3.LT.LYRHGT(NX,NY,I))) THEN
        NEWNZ = I
        GOTO 1
        ENDIF
        ENDIF
  118   CONTINUE
   1    CONTINUE
        RETURN
        END
C
      SUBROUTINE DICOLLI(NX,NY,MEANFP,DICO,NL,MXLYR,MXCL,ISEED)
      REAL DICO,  MEANFP(0:MXCL,0:MXCL,MXLYR) 
      REAL  ran1, R
      INTEGER ISEED
      R = ran1(ISEED)
c      write(*,*) r, nx, ny, nl
      IF (R.NE.(0.)) THEN
         DICO = - MEANFP(NX,NY,NL) *  LOG(R)
      ELSE
         DICO = 10000.
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE DISBOUNDS(NX,NY,NL,LYRHGT,POS,DIBOMIN,COSS, 
     $     EXIT,POSBOUND,SIZE,NEWNX,NEWNY,NEWNZ,NXCL,NYCL,MXLYR,MXCL,
     $     ISEED)
      REAL DIBO(3), LYRHGT(0:MXCL,0:MXCL,0:MXLYR), POS(3), 
     $     COSS(3), POSBOUND(3), 
     $     DIBOMIN, SIZE, TESTDIBO(3), NEWPOS(3)
      INTEGER  EXIT, NEWX(2), NEWY(2)
      LOGICAL LOG1,LOG2
C     
C     
      IF (COSS(3).GT.0) THEN
         DIBO(3) = (LYRHGT(NX,NY,NL) - POS(3)) / COSS(3)
         POSBOUND(3) = LYRHGT(NX,NY,NL)
         NEWNZ = NL+1
      ELSE
         IF (COSS(3).LT.0) THEN
            DIBO(3) = (LYRHGT(NX,NY,NL-1) - POS(3)) / COSS(3)
            POSBOUND(3) = LYRHGT(NX,NY,NL-1)
            NEWNZ = NL - 1
         ELSE 
            DIBO(3) = 10000.
         ENDIF
      ENDIF
      IF ((NX.NE.0).AND.(NY.NE.0)) THEN
         IF (COSS(2).GT.0) THEN
            DIBO(2) = (FLOAT(NY) * SIZE - POS(2)) / COSS(2)
            POSBOUND(2) = FLOAT(NY) * SIZE 
            NEWY(2) = NY+1 
            NEWX(2) = NX
            IF (NEWY(2).GT.NYCL) THEN 
               NEWX(2) = 0
               NEWY(2) = 0
            ENDIF
         ELSE
            IF (COSS(2).LT.0) THEN
               DIBO(2) = ( (FLOAT(NY-1)) * SIZE - POS(2)) / COSS(2)
               POSBOUND(2) =FLOAT((NY-1)) * SIZE 
               NEWY(2) = NY-1
               NEWX(2) = NX
               IF (NEWY(2).LT.1) THEN
                  NEWX(2) = 0
                  NEWY(2) = 0
               ENDIF
            ELSE
               DIBO(2) = 10000.
            ENDIF
         ENDIF
         IF (COSS(1).GT.0) THEN
            DIBO(1) = (FLOAT(NX) * SIZE - POS(1)) / COSS(1)
            POSBOUND(1) =FLOAT(NX) * SIZE
            NEWX(1) = NX+1
            NEWY(1) = NY
            IF (NEWX(1).GT.NXCL) THEN
               NEWX(1) = 0
               NEWY(1) = 0
            ENDIF
         ELSE
            IF (COSS(1).LT.0) THEN
               DIBO(1) = (FLOAT(NX-1)*SIZE - POS(1)) / COSS(1)
               POSBOUND(1) =FLOAT(NX-1) *  SIZE
               NEWX(1) = NX-1
               NEWY(1) = NY
               IF (NEWX(1).LT.1) THEN
                  NEWX(1) = 0
                  NEWY(1) = 0
               ENDIF
            ELSE 
               DIBO(1) = 10000.
            ENDIF
         ENDIF
      ELSE
         IF ((NX.EQ.0).AND.(NY.EQ.0)) THEN
            IF ((COSS(2).GT.0).AND.(POS(2).LT.0)) THEN
               TESTDIBO(2) = -POS(2)/COSS(2)
               NEWPOS(1) = POS(1)+TESTDIBO(2)*COSS(1)
               LOG1=( NEWPOS(1).GT.(0.))
               LOG2=( NEWPOS(1).LT.(FLOAT(NXCL)*SIZE))
               IF (LOG1.AND.LOG2) THEN
                  DIBO(2) = TESTDIBO(2) 
                  POSBOUND(2) = 0. 
                  NEWY(2) = 1
                  CALL GETNEWN(NEWX(2),NEWPOS(1),NXCL,SIZE,ISEED) 
               ELSE
                  DIBO(2) =10000.
               ENDIF
            ELSE
               IF((COSS(2).LT.0).AND.
     $              (POS(2).GT.(FLOAT(NYCL)*SIZE))) THEN
                  TESTDIBO(2) = ((FLOAT(NYCL)*SIZE)-POS(2))/COSS(2)
                  NEWPOS(1) = POS(1)+TESTDIBO(2)*COSS(1)
                  LOG1=( NEWPOS(1).GT.(0.))
                  LOG2=( NEWPOS(1).LT.(FLOAT(NXCL)*SIZE))
                  IF (LOG1.AND.LOG2) THEN
                     DIBO(2) = TESTDIBO(2)
                     POSBOUND(2) = FLOAT(NYCL) * SIZE   
                     NEWY(2) = NYCL 
                     CALL GETNEWN(NEWX(2),NEWPOS(1),NXCL,SIZE,ISEED)
                  ELSE
                     DIBO(2) =10000.
                  ENDIF
               ELSE
                  DIBO(2) = 10000.
               ENDIF
            ENDIF
            IF ((COSS(1).GT.0).AND.(POS(1).LT.0)) THEN
               TESTDIBO(1) = -POS(1)/COSS(1)
               NEWPOS(2) = POS(2)+TESTDIBO(1)*COSS(2)
               LOG1=( NEWPOS(2).GT.(0.))
               LOG2=( NEWPOS(2).LT.(FLOAT(NYCL)*SIZE))
               IF (LOG1.AND.LOG2) THEN
                  DIBO(1) = TESTDIBO(1)
                  POSBOUND(1) = 0.
                  NEWX(1) = 1
                  CALL GETNEWN(NEWY(1),NEWPOS(2),NYCL,SIZE,ISEED)
               ELSE
                  DIBO(1) =10000.
               ENDIF
            ELSE
               IF ((COSS(1).LT.0).AND.
     $              (POS(1).GT.(FLOAT(NXCL)*SIZE))) THEN
                  TESTDIBO(1) = ((FLOAT(NXCL)*SIZE)-POS(1))/COSS(1)
                  NEWPOS(2) = POS(2)+TESTDIBO(1)*COSS(2)
                  LOG1=( NEWPOS(2).GT.(0.))
                  LOG2=( NEWPOS(2).LT.(FLOAT(NYCL)*SIZE))
                  IF (LOG1.AND.LOG2) THEN
                     DIBO(1) = TESTDIBO(1)
                     POSBOUND(1) = FLOAT(NXCL) * SIZE
                     NEWX(1) = NXCL
                     CALL GETNEWN(NEWY(1),NEWPOS(2),NYCL,SIZE,ISEED)
                  ELSE
                     DIBO(1) =10000.
                  ENDIF
               ELSE
                  DIBO(1) = 10000.
               ENDIF
            ENDIF
         ELSE
            WRITE(*,*) 'ERR'
         ENDIF
      ENDIF
      DIBOMIN = DIBO(3)
      EXIT = 3
      DO 66 I=1,2
         IF (DIBO(I).LT.DIBOMIN) THEN
            DIBOMIN = DIBO(I)
            EXIT = I
         ENDIF
 66   CONTINUE
      IF (EXIT.NE.3) THEN
         NEWNX = NEWX(EXIT)
         NEWNY = NEWY(EXIT)
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE GETNEWN(NEWN,NEWPOS,NCL,SIZE,ISEED)
      REAL NEWPOS,SIZE
      DO 11 I = 1, NCL
         IF ((NEWPOS.GE.(I-1)*SIZE).AND.(NEWPOS.LE.(I*SIZE))) THEN
            NEWN = I
            GOTO 91
         ENDIF
 11   CONTINUE
 91   RETURN
      END
C     
C     
      SUBROUTINE SCATCLOUDS(NX,NY,COSS,ASYM,NL,HENYEY,PHADI,ANG,
     $     NANG,MXLYR,MXCL,MXAN, ISEED)
      REAL COSS(3),ASYM(0:MXCL,0:MXCL,MXLYR),G,COST,SINS,PSI,A,B,C,D,
     $     ANG(0:MXAN), PHADI(0:MXCL,0:MXCL,MXLYR,0:MXAN)
      REAL VAR,VAR1,VAR2,VAR3
      REAL ran1,R
      LOGICAL HENYEY
      INTEGER ISEED
C     
      PI = 4.*ATAN(1.)
      R = ran1(ISEED)
      IF (HENYEY)THEN
         
         G = ASYM(NX,NY,NL)
         IF (G.NE.0) THEN
            COST = (1. + G**2 - 
     $           ((1. - G**2)/(G *(2.*R - 1.) + 1.))**2)/2./G
         ELSE
            COST = 1.- 2. * R
         ENDIF
         
      ELSE 
         
         I=0
 22      IF ((R-PHADI(NX,NY,NL,I)).GE.0) THEN
            IF (I.NE.NANG) THEN
               I=I+1
            ENDIF
            GO TO 22
         ENDIF
         RANG1 = ANG(I)*PI/180.
         RANG2 = ANG(I-1)*PI/180.
         RC1=COS(RANG1)
         RC2=COS(RANG2)
         RP1= PHADI(NX,NY,NL,I)
         RP2= PHADI(NX,NY,NL,I-1)
C     Suitable when the phase function is a step function
         COST = RC1-(RP1-R)/(RP1-RP2)*(RC1-RC2) 
C     Suitable when the phase function is concave up
C     COST = (RC1**2-(RP1-R)/(RP1-RP2)*(RC1**2-RC2**2))**0.5 
C     Suitable when the phase function is concave down
C     COST = (RC1**2-(-RP2+R)/(RP1-RP2)*(RC1**2-RC2**2))**0.5 
         
      ENDIF	
c     JLH	IF (COST.NE.(0.)) THEN
c     JLH    SINS = SQRT(1. - COST**2)
c     JLH	ELSE
c     JLH	SINS = 1.
c     JLH	ENDIF
c     JLH added next line
      SINS = SQRT(MAX(0.,1.-COST**2))
      
c     JLH    IF (ABS(ABS(COST)-(1.)).LT.(0.0001)) THEN
c     IF ((COST.EQ.1.).OR.(COST.EQ.-1.)) THEN
c     JLH	SINS = 0.
c     JLH	ENDIF
      A = COST
      B = SINS
      R = ran1(ISEED)
      PSI = (R * 2-1.) *4* ATAN(1.)
      C = COS(PSI)
      D = SIN(PSI) 
      IF (ABS(ABS(COSS(3))-(1.)).LT.(0.0001)) THEN
         COSS(3) = A * COSS(3)
         COSS(2) = B * D 
         COSS(1) = B * C
      ELSE
c     JLH	IF (COSS(3).NE.(0.)) THEN 
c     JLH	VAR = SQRT(1.- COSS(3)**2)
c     JLH	ELSE
c     JLH	VAR = 1.
c     JLH added next line
         VAR = SQRT(MAX(0.,1.-COSS(3)**2))
c     JLH	ENDIF
         VAR1=COSS(1)
         VAR2=COSS(2)
         VAR3=COSS(3)
         COSS(3) = - B * C * VAR + A * VAR3
         COSS(2) = (B*C*VAR3*VAR2 + B*D*VAR1) / VAR + A * VAR2
         COSS(1) = (B*C*VAR3*VAR1 - B*D*VAR2) / VAR + A * VAR1
      ENDIF
c      write(*,*) coss(3), var3, NL
      RETURN
      END
C     
C     
      SUBROUTINE ABSORBT(NX,NY,NL, POS,LYRHGT,B0,B1,TBAIR,
     $     MXLYR,MXCL,ISEED)
      REAL POS(3), LYRHGT(0:MXCL,0:MXCL,0:MXLYR), 
     $     B0(0:MXCL,0:MXCL,MXLYR), B1(0:MXCL,0:MXCL,MXLYR)
      REAL  TBAIR
      TBAIR = B0(NX,NY,NL) + B1(NX,NY,NL) 
     $     * (POS(3) - LYRHGT(NX,NY,NL-1))
      RETURN
      END
C     
C     
C     
C
      SUBROUTINE NEWDICOLLI(NX,NY,NZ,NXOLD,NYOLD,NZOLD,DICO,
     $     MEANFP,  DIBOMIN, MXLYR, MXCL, ISEED)
      REAL DICO,  DIBOMIN 
      REAL MEANFP(0:MXCL,0:MXCL,MXLYR)
      DICO=MEANFP(NX,NY,NZ)/MEANFP(NXOLD,NYOLD,NZOLD)*(DICO-DIBOMIN)
      RETURN
      END
C     
      SUBROUTINE SCATGROUND(LAMBERT, COSS, ISEED)
      LOGICAL LAMBERT
      REAL COSS(3) 
      REAL  ran1,R
      INTEGER ISEED
      IF (LAMBERT) THEN
         R = ran1(ISEED)
         IF (R.NE.(1.)) THEN 
            COSS(3) = SQRT(R)
         ELSE 
            COSS(3)  = 1.
         ENDIF
         R = ran1(ISEED)
         PSI = 2 * 4. * ATAN(1.) * R
         COSS(2) = SIN(PSI) * SQRT(1. - COSS(3)**2)
         COSS(1) = COS(PSI) * SQRT(1. - COSS(3)**2)
      ELSE
         COSS(3) = - COSS(3)
         COSS(2) =  COSS(2)
         COSS(1) =  COSS(1)
      ENDIF
      RETURN
      END
C     
C     
      FUNCTION URAND(IY)
C     $ (HMW)
C     TAKEN FROM FORSYTHE,MALCOLM,MOLER: "COMPUTER METHODS
C     FOR MATHEMATICAL COMPUTATIONS," PRENTICE HALL,1977,259PP.
C     THE ROUTINE CAN BE FOUND IN LAST CHAPTER, P246.
C     PROGRAMMER: LOUIS GARAND, JULY 1983
C     
C     URAND IS A UNIFORM RANDOM NUMBER GENERATOR.  THE 'SEED' IY SHOULD
C     BE INITIALIZED TO AN ARBITRARY INTEGER PRIOR TO THE FIRST CALL
C     TO URAND.  THE CALLING PROGRAM SHOULD NOT ALTER THE VALUE OF IY
C     WHICH IS CHANGED INTERNALLY BY URAND ITSELF.
C     VALUES OF URAND ARE RETURNED IN THE RANGE 0. TO 1.
C     
      REAL HALFM,DATAN,DSQRT
      DATA M2/0/,ITWO/2/
      IF(M2.NE.0) GO TO 120
C     
C     * ON FIRST ENTRY, COMPUTE MACHINE INTEGER WORD LENGTH
      M=1
 110  M2=M
      M=ITWO*M2
      IF(M.GT.M2) GO TO 110
      HALFM=M2
C     
C     * COMPUTE MULTIPLIER AND INCREMENT FOR LINEAR CONGRUENTIAL METHOD
      IA=8*IDINT(HALFM*ATAN(1.D0)/8.D0) + 5
      IC=2*IDINT(HALFM*(0.5D0-SQRT(3.D0)/6.D0)) +1
      MIC=(M2-IC)+M2
C     * OBTAIN SCALE FACTOR FOR CONVERTING TO FLOATING POINT
      S=0.5/HALFM
C     
C     * COMPUTE NEXT RANDOM NUMBER
 120  IY=IY*IA
C     * FOR COMPUTERS NOT ALLOWING INTEGER OVERFLOW ON ADDITION:
      IF(IY.GT.MIC) IY=(IY-M2)-M2
      IY=IY+IC
C     * FOR COMPUTERS WITH WORD LENGTH FOR ADDITION GREATER THAN FOR
C     MULTIPLICATION:
      IF(IY/2.GT.M2) IY=(IY-M2)-M2
C     * FOR COMPUTERS WHERE INTEGER OVERFLOW AFFECTS THE SIGN BIT:
      IF(IY.LT.0) IY=(IY+M2)+M2
      URAND=FLOAT(IY)*S
      RETURN
      END

