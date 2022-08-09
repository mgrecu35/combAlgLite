C
      PROGRAM MC3
C-----------------------------------------------------------------------
C
C     Program to calculate microwave brightness temperatures
C     from horizontally and vertically structured clouds.
C     MC3D is the calling program.  It reads the cloud  
C     parameters from a file and
C     calls the radiative transfer subroutine (BRITEMP).
C     Internal subroutines: 
C     PHANG:     Computes the integral of the phase function
C     PMOM:      Computes the phase function given its moments
C     EMIT       Computes the surface emissivity
C
C     AUTHORS:   Laura Roberti 
C                Dpt. of Elettronica, Politecnico di Torino
C                Corso Duca degli Abruzzi, 24
C                10129, Torino, Italy
C                e-mail: PERONA@polito.it; Fax: ++39-11-564-4015
C
C                Christian Kummerow
C                Laboratory for Atmospheres, NASA/Goddard Space Flight Center
C                Code 912
C                20716 Greenbelt, Maryland, USA 
C                e-mail: kummerow@audry.gsfc.nasa.gov
                
C     First version: September 1993
C     Revised:       August 1994
C
C   Argument list:
C   VAR         I/O   TYPE        DESCRIPTION
C   ---         ---   ----        -----------
C   MXCL         I    I           MAX NUM OF SUBCLOUDS IN X AND Y DIRECTION 
C   MXLYR        I    I           MAX NUM OF  LAYERS 
C   MXMO         I    I           MAX NUM OF PHASE FUNC. MOMENTS
C   MXAN         I    I           MAX NUM OF ANGLES IN DEG. FOR THE PHASE FUNC. 
C   MXEMI        I    I           NUM OF ANGLES IN DEG. FOR THE GROUND
C 			          EMISSIVITY COMPUTATION.
C   STEP         I    R           DISTANCE IN THE X AND Y DIRECTION BETWEEN THE 
C                                 COMPUTATIONAL POINTS ON THE SUBCLOUD TOP 
C                                 IN WHICH  THE BRIGHT. TEMP. IS COMPUTED
C   AZIM         I    R           AZIMUTHAL VIEWING ANGLE IN DEGREE
C   UMU          I    R           COSINE OF POLAR VIEWING ANGLE 
C   FISOT        I    R           BACKGROUND TEMPERATURE (2.7 K)
C   LAMBERT      I   LOG          SURFACE TYPE
C   HENYEY       I   LOG          PHASE FUNCTION INPUT INFORMATION
C   PHAMOM       I   LOG          PHASE FUNCTION INPUT INFORMATION
C   NXCL         I    I           NUMBER OF CLOUDS IN X DIRECTIONS
C   NYCL         I    I           NUMBER OF CLOUDS IN Y DIRECTIONS
C   SIZE         I    R           SIDE OF EACH SQUARE CLOUD (KM)
C   BTEMP        I    R           SURFACE TEMPERATURE [K]
C   BTEMP2       I    R           BACKGROUND SURFACE TEMPERATURE [K]
C   EMIS         I    R           SURFACE EMISSIVITY 
C   EMIS2        I    R           BKGROUND SURFACE EMISSIVITY 
C   EMI          I    R           SURFACE EMISSIVITY FOR DIFFERENT ANGLES
C   EMI2         I    R           BKGROUND SURFACE EMISSIVITY FOR DIFF. ANGLES
C   LYRTEMP      I    R           TEMPERATURE AT THE TOP OF THE LAYERS [K]
C   KEXTTOT      I    R           EXTINCTION COEFF. (1/KM) 
C   SALBTOT      I    R           SINGLE SCATTERING ALBEDO 
C   ASYMTOT      I    R           ASYMMETRY  COEFF.   
C   MOM          I    R           PHASE FUNCTION MOMENTS 
C   ANG          I    R           ANGLES  IN DEG. FOR THE TABULATION 
C                                 OF THE PHASE FUNCTION
C   PHA          I    R           PHASE FUNCTION VALUES IN THE ABOVE ANGLES 
C   NLYR         I    R           NUMBER OF LAYERS 
C   FILEPHA      I   A*50         INPUT FILE CONTAINING THE PHASE FUNCTION
C   FILEINP      I   A*50         IPUT FILE CONTAINING THE ATM. AND CLOUD PARAM.
C   FILEOUT      O   A*50         OUTPUT FILE WITH THE BRIGHT. TEMPERATURES 
C   WIND         I    R           WIND SPEED OVER WATER [m/s]
C   FREQY        I    R           FREQUENCY [GHz]
C   POLN         I   A*1          POLARIZATION  ('V'/'H')                       
C   NEMI         I    I           NEMI=1 SURFACE EMISSIVITY IS COMPUTED
C                                 FOR DIFFERENT POLAR ANGLES
C				  NEMI=0 SURFACE EMISSIVITY IS READ FROM
C                                 THE INPUT FILE
C 
C                 SUBROUTINES CALLED
C   
C  EMIT : computes the surface emissivity as a function of the polar
C         angle, frequency, wind, surface temperature and polarization
C
C  PHANG: Computes the integral of the phase function with respect to
C         the polar angle.
C
C  PHAMO: Computes the phase functions in a set of points, given its
C         moments            
C
C  BT3  : Is the main subroutine that  performs the radiative transfer
C         computations
C
      PARAMETER(MXLYR=49,MXCL=1,MXMO=1,MXAN=180,MXEMI=89,NPHO=100000)
      LOGICAL  LAMBERT,PHAMOM,HENYEY
      REAL  LYRTEMP(0:MXCL,0:MXCL,0:MXLYR),EMI(0:MXEMI),EMI2(0:MXEMI),  
     $     EMIS,EMIS2, tskin(0:MXCL,0:MXCL),
     $     LYRHGT(0:MXCL,0:MXCL,0:MXLYR),
     $     MOM(0:MXCL,0:MXCL,MXLYR,0:MXMO),
     $     ANG(0:MXAN),PHA(0:MXCL,0:MXCL,MXLYR,0:MXAN), 
     $     KEXTTOT(0:MXCL,0:MXCL,MXLYR),
     $     PHADI(0:MXCL,0:MXCL,MXLYR,0:MXAN), 
     $     ASYMTOT(0:MXCL,0:MXCL,MXLYR), SALBTOT(0:MXCL,0:MXCL,MXLYR) 
      REAL  MEANFP(0:MXCL,0:MXCL,MXLYR), WIND, FREQY, EBAR,
     $     B0(0:MXCL,0:MXCL,MXLYR), B1(0:MXCL,0:MXCL,MXLYR),
     $     BTEMP(0:MXCL,0:MXCL)
      CHARACTER*50 FILEINP,FILEOUT,FILEPHA
      CHARACTER POLN*1
      REAL RANG(0:MXAN),RCOS(0:MXAN), RPOL(0:MXMO,0:MXAN),
     $     P1(0:MXAN),P2(0:MXAN),P3(0:MXMO,0:MXAN)
      INTEGER NLYR(0:MXCL,0:MXCL)
      real inc_angle_mwr 
      inc_angle_mwr= 52.8
C     
C     
C     Hard coded input variables that define the type of radiative 
C     transfer computationn
C     


   
      AZIM = 180.
      UMU = 0.6
      umu = cos(inc_angle_mwr*3.141593/180.)
      size=2.2
      FISOT =  2.7
      LAMBERT = .FALSE.
      STEP =  500.          
      PHAMOM =.FALSE.
      HENYEY = .TRUE.
      FILEINP='INPUT.DAT'
      FILEOUT='OUT.DAT'
      FREQY=85.6
      POLN='H'
      WIND=0.
c     JLH
      NEMI = 1
      GOTO 314
c     

 314  CONTINUE
      
C     
C     READ THE CLOUD PARAMETERS 
C     
c   em   OPEN(UNIT=1,FILE=FILEINP,STATUS='OLD')
c      READ(1,*)NXCL,NYCL,SIZE,BTEMP2,BTEMP,EMIS2,EMIS
      NXCL=1
      NYCL=1

      

C     ******  Do the phase function if not HENYEY-GREENSTEIN *****
C
      IF (.NOT.HENYEY) THEN
C     
C     ************************************************************
C     ********  Phase function tabulated by points   **************
C     
      NANG = 180

C     
      ENDIF

C
C      ******************************************************************
C      ------------  Background atmospheric parameters  ----------------- 
C
c        READ(1,*)NX,NY,NLYR(NX,NY),LYRHGT(NX,NY,0),LYRTEMP(NX,NY,0)
c        DO  50  I = 1,NLYR(NX,NY)
c           READ(1,*) ILYR, LYRHGT(NX,NY,I), LYRTEMP(NX,NY,I), 
c     $          KEXTTOT(NX,NY,I),SALBTOT(NX,NY,I),ASYMTOT(NX,NY,I) 
c  50    CONTINUE
C
C       ---------------      Cloud parameters    ------------------------ 

C     
      
C
C       ************ GO TO RADIATIVE TRANSFER SUBROUTINE **************
C     

      call openfile()
      do irec=1,1000
         call readkext(kexttot(1,1,:), salbtot(1,1,:), asymtot(1,1,:), 
     &        emis, lyrtemp(1,1,:), nz, tbe)
         !print*, nz, tbe,emis,asymtot(1,1,1:nz)

         emi2=emi
      !print*, lyrtemp(1,1,1:nz)
      
         emis2=emis
         tskin(1,1)=lyrtemp(1,1,0)
         do i=1,nz
            kexttot(0,0,i)=kexttot(1,1,i)
            asymtot(0,0,i)=asymtot(1,1,i)
            salbtot(0,0,i)=salbtot(1,1,i)
         enddo
         do i=0,nz
            lyrhgt(1,1,i)=i*0.5
         enddo
         
         do i=0,nz
            lyrhgt(0,0,i)=lyrhgt(1,1,i)
            lyrtemp(0,0,i)=lyrtemp(1,1,i)
            tskin(0,0)=tskin(1,1)
         enddo
         NLYR=nz
         BTEMP=tskin
         nxcl=1
         nycl=1
         SIZE=1000.
         step=500.
         !print*, tbe, emis
         if (nz.eq.37 .and. irec/2*2==irec) then
            CALL BT3( SIZE, NLYR, UMU, AZIM, LAMBERT, FISOT,
     $           BTEMP2, EMIS2, BTEMP, LYRTEMP, LYRHGT, KEXTTOT,
     $           SALBTOT, ASYMTOT,EMIS, STEP, NXCL, NYCL,FILEOUT,
     $           HENYEY,NANG,ANG,PHADI,EMI,EMI2,NEMI,
     $           MEANFP,B0,B1,MXLYR,MXCL,MXAN,MXEMI,NPHO,tsub)
            print*, tbe, tsub
         endif
      enddo
 1001 FORMAT(I3,2X,I3,2X,F8.4,2X,F8.4,2X,F8.4,2X,F8.4,2X,F8.4)
 1011 FORMAT(I3,2X,I3,1X,I3,1X,F8.4,2X,F8.4)
 1021 FORMAT(5X,I3,5X,F8.4,2X,F8.4,2X,F8.4,2X,F8.4,2X,F8.4)
      STOP
      END
C
C       
      SUBROUTINE  PHANG(NXCL,NYCL,NLYR,NANG,ANG,PHA,PHADI,
     $              MXLYR, MXCL, MXAN)
      REAL ANG(0:MXAN),PHA(0:MXCL,0:MXCL,MXLYR,0:MXAN),
     $     PHADI(0:MXCL,0:MXCL,MXLYR,0:MXAN)
      INTEGER NLYR(0:MXCL,0:MXCL)
C     
      PI=4.*ATAN(1.)
      
      NX=0
      NY=0
      DO  NL = 1,NLYR(NX,NY)
         PHADI(NX,NY,NL,0)=0.
         DO  I=1,NANG
            RC1=COS(ANG(I)*PI/180.)
            RC2=COS(ANG(I-1)*PI/180.)
            PHADI(NX,NY,NL,I)=PHADI(NX,NY,NL,I-1)+
     $           (PHA(NX,NY,NL,I-1)+PHA(NX,NY,NL,I))*(RC1-RC2)/2.
         enddo
         !write(*,*) phadi(nx,ny,nl,nang)
      enddo
      
      PI=4.*ATAN(1.)
      DO  NX = 1,NXCL 
         DO  NY = 1,NYCL 
            DO  NL = 1,NLYR(NX,NY)
               PHADI(NX,NY,NL,0)=0.
               DO  I=1,NANG
                  RC1=COS(ANG(I)*PI/180.)
                  RC2=COS(ANG(I-1)*PI/180.)
                  PHADI(NX,NY,NL,I)=PHADI(NX,NY,NL,I-1)+
     $                 (PHA(NX,NY,NL,I-1)+PHA(NX,NY,NL,I))*(RC1-RC2)/2.
               enddo
               write(*,*) phadi(nx,ny,nl,nang)
            enddo
         enddo

      enddo
      stop
      RETURN
      END
C     
C     
      SUBROUTINE PHAMO(NXCL,NYCL,NLYR,NMOM,ANG,MOM,PHADI,
     $     RCOS,RPOL,PHA,RANG,P1,P2,P3,
     $     MXLYR, MXCL, MXMO,MXAN)
      REAL ANG(0:MXAN),RANG(0:MXAN),RCOS(0:MXAN),RPOL(0:MXMO,0:MXAN),
     $     PHA(0:MXCL,0:MXCL,MXLYR,0:MXAN),
     $     PHADI(0:MXCL,0:MXCL,MXLYR,0:MXAN),
     $     P1(0:MXAN),P2(0:MXAN),P3(0:MXMO,0:MXAN),
     $     MOM(0:MXCL,0:MXCL,MXLYR,0:MXMO)
      INTEGER NLYR(0:MXCL,0:MXCL), RI
C     
      DO 54 NX=0,NXCL
         DO 54 NY=0,NXCL
            DO 54 NL=1,NLYR(NX,NY)
               DO 54 J=0,NMOM
 54    CONTINUE
       PI = 4.*ATAN(1.)
C     
       DO 555 I=0,127
          P1(I) =0.
          P2(I)=0.
          DO II=0,180
             RPOL(I,II)=0.
          END DO
 555   CONTINUE
C
       NANG=180
C     
       DO 565 RI=180,0,-1
          ANG(180-RI)=RI+0.0
          RANG(180-RI)=(RI+0.0)*PI/180.	
          RCOS(180-RI) = COS (RANG(180-RI))
  565   CONTINUE
C
        P1(0) = 1.         	
        P2(0) = 0. 
        P2(1) = 1.
C
        DO 111 I=2,NMOM
           DO 113 J=0,I
              PP1=-FLOAT(I-1)/FLOAT(I)*P1(J)
              P1(J)=P2(J)
              IF (J.GE.1) THEN
                 PP2=P2(J-1)*(2*FLOAT(I)-1)/FLOAT(I)
              ELSE
                 PP2 = 0.
              ENDIF
              P3(I,J)=PP1+PP2
              IF (J.GE.1) THEN
                 P2(J-1)=P3(I,J-1)
              ENDIF
 113       CONTINUE
           P2(I)=P3(I,I)
 111    CONTINUE
C     
        DO 112 II=0,180
           DO 112 I=2,NMOM
              RPOL(I,II)=0.
              DO 115 J=0,I 
                 RBU= RPOL(I,II)
C     RPOL(I,II) = RPOL(I,II)+P3(I,J)*P3(I,J)*(RCOS(II)**FLOAT(J))	 
    	RPOL(I,II) = RPOL(I,II)+P3(I,J)*(RCOS(II)**FLOAT(J))	 
 115    CONTINUE
 112	CONTINUE
        NX=0
        NY=0
        DO 11 NL = 1,NLYR(NX,NY)
           DO 11 II=0,NANG
              PHA(NX,NY,NL,II)=MOM(NX,NY,NL,0)+MOM(NX,NY,NL,1)*RCOS(II)	
              DO 11 I =2,NMOM
                 PHA(NX,NY,NL,II) = PHA(NX,NY,NL,II) +
     $                MOM(NX,NY,NL,I)*RPOL(I,II)	
 11     CONTINUE
        DO 118 NX = 1,NXCL 
           DO 118 NY = 1,NYCL 
              DO 118 NL = 1,NLYR(NX,NY)
                 DO 118 II=0,NANG
                    PHA(NX,NY,NL,II)=MOM(NX,NY,NL,0)+
     $                   MOM(NX,NY,NL,1)*RCOS(II)	
                    DO 118 I =2,NMOM
        PHA(NX,NY,NL,II) = PHA(NX,NY,NL,II) + MOM(NX,NY,NL,I)*RPOL(I,II)	
 118    CONTINUE
        DO 18 NX = 1,NXCL 
           DO 18 NY = 1,NYCL 
              DO 18 NL = 1,NLYR(NX,NY)
                 DO 18 II=0,NANG
        PHA(NX,NY,NL,II)=PHA(NX,NY,NL,II)/2.
 18    CONTINUE
        CALL PHANG(NXCL,NYCL,NLYR,NANG,ANG,PHA,PHADI,
     $             MXLYR, MXCL, MXAN)
        RETURN
        END



      SUBROUTINE  EMIT( F, POLN, TS, W, UMU, EMIS, EBAR )

      PARAMETER ( NANG = 21 )
      CHARACTER POLN*1
      REAL  ANG(NANG), MU(NANG), ESUM(NANG)

      PI = 2.*ASIN(1.0)
      S = 35.                                    ! SALINITY IN PPM
      ANGLE = ACOS(UMU)*180./PI
C
C     CALCULATE EMIS AT GIVEN ANGLE
      CALL EMISS (F,ANGLE,S,TS,W,EV,EH,EMISH,EMISV)
      IF ( POLN .EQ. 'H' ) EMIS = EH
      IF ( POLN .EQ. 'V' ) EMIS = EV
C
C     CALCULATE EMIS AT VARIOUS ANGLES
      DO 58  I = 1,NANG
       ANG(I) = 4.*( I - 1 )
       MU(I) = COS(ANG(I)*PI/180.)
       ANGLES = ANG(I)
       CALL EMISS (F,ANGLES,S,TS,W,EV,EH,EMISH,EMISV)
       ESUM(I) =  EV + EH
  58  CONTINUE
C
C     CALCULATE EBAR
     SUM = 0.0
      DO 59  I = 1,NANG-1
       EAVG = 0.5*( ESUM(I) + ESUM(I+1) )
       DMU = MU(I) - MU(I+1)
       AVMU = 0.5*( MU(I) + MU(I+1) )
       SUM = SUM + EAVG*AVMU*DMU
  59  CONTINUE
      EBAR = SUM
C
  104 FORMAT(1X,F8.2,F10.2,F9.3,A6,F7.2,2F11.5)
      RETURN
      END
C
C
      SUBROUTINE  DIECON(S,T,FREQ,E1,E2)
      DATA TWOPI /6.283185307/ , EINF /4.9/ , SOLD /0.0/ , TOLD /-99./
      IF (S .EQ. SOLD .AND. T .EQ. TOLD) GO TO 10
      ST = S*T
      S2 = S*S
      T2 = T*T
      SST = S2*T
      STT = T2*S
      SSTT = S2*T2
      ES = 88.-4.339E-01*S+1.71E-03*S2-4.035E-01*T+8.065E-04*T2+6.170
     $  E-03 * ST-8.910E-05*SST-6.934E-05*STT+1.439E-06*SSTT
C
      TAU = (18.70-7.924E-02*S+6.35E-04*S2-5.489E-01*T+5.758E-03*T2+
     $1.889E-03*ST-7.209E-06*SST-5.299E-07*STT-2.101E-07*SSTT)*1.0E-12
C
      SIGMA = (7.788E-03*S-1.672E-06*S2-8.570E-15*T+2.996E-16*T2+4.059E
     $     -04 * ST-3.215E-06*SST-1.423E-06*STT+3.229E-08*SSTT)*1.0E11
C
   10 ZNU = FREQ*1.E09
      OMEGA = TWOPI*ZNU
      DEN = 1. + (OMEGA * TAU) ** 2
      E1 = (ES-EINF)/DEN+EINF
      E2 = (ES-EINF)*OMEGA*TAU/DEN+2.*SIGMA/ZNU
      SOLD = S
      TOLD = T
C
      RETURN
      END


      SUBROUTINE  EMISS(F,ANGLE,S,TS,W,EV,EH,EMISH,EMISV)
      DATA DTR / 0.01745329252 /
      T = TS-273.16
      THETA = ANGLE*DTR
C
      CALL DIECON(S,T,F,E1,E2)
C
      CMHU = COS (THETA)
      CMHU2 = CMHU*CMHU
      FACT1 = E1+CMHU2-1.
      B = (FACT1*FACT1+E2*E2)**0.25
      B2 = B*B
      ARG = E2/FACT1
      PSIZ = 0.5*ATAN(ARG)
      COSPZ = COS(PSIZ)*CMHU*B
      F11 = CMHU2+B2+2.*COSPZ
      PSIY = ATAN(E2/E1)
      G = CMHU*SQRT(E1*E1+E2*E2)
      COSPY = COS(PSIY-PSIZ)*G*B
      F22 = G*G+B2+2.*COSPY
C
C     FOR SPECULAR SURFACES THE HORIZONTAL AND VERTICAL EMISSIVITY
C     CALCULATION
C
      EMISH = 4.*COSPZ/F11
      EMISV = 4.*COSPY/F22
C
C     FROM HERE THE EFFECT OF SURFACE ROUGHNESS AND FOAM ARE INCLUDED
C     BASED ON HOLLINGER MODEL AND FOAM MODEL OF STOGRYN
C
      RSH = 1.-EMISH
      RSV = 1.-EMISV
C
      FOAM = 7.751E-06 * W ** 3.231
C
      GH = 1.-1.748E-3*ANGLE-7.336E-5*ANGLE**2+1.044E-7*ANGLE**3
      GV = 1.-9.946E-4*ANGLE+3.218E-5*ANGLE**2-1.187E-6*ANGLE**3
     $   +7.0E-20*ANGLE**10
C
      A1 = (208.0+1.29*F)/TS
C
      RFV = 1. - A1 * GV
      RFH = 1. - A1 * GH
C
      Y = 7.32E-02*ANGLE
C
C     TS SURFACE TEMP IS IN DEGREE KELVIN
C
      SQRTF = SQRT(F)
C
      CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)
      CORRH = (W*(1.15E-01+3.80E-05*ANGLE**2)*SQRTF/TS)
C
      RRV = RSV-CORRV
      RRH = RSH-CORRH
C
      RV = RRV*(1.-FOAM)+RFV*FOAM
      RH = RRH*(1.-FOAM)+RFH*FOAM
C
      EH = 1.-RH
      EV = 1.-RV
C
      RETURN
      END


