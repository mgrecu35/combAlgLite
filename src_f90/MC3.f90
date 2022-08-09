!C
      
module MC3
!C-----------------------------------------------------------------------
!C
!C     Program to calculate microwave brightness temperatures
!C     from horizontally and vertically structured clouds.
!C     MC3D is the calling program.  It reads the cloud  
!C     parameters from a file and
!C     calls the radiative transfer subroutine (BRITEMP).
!C     Internal subroutines: 
!C     PHANG:     Computes the integral of the phase function
!C     PMOM:      Computes the phase function given its moments
!C     EMIT       Computes the surface emissivity
!C
!C     AUTHORS:   Laura Roberti 
!C                Dpt. of Elettronica, Politecnico di Torino
!C                Corso Duca degli Abruzzi, 24
!C                10129, Torino, Italy
!C                e-mail: PERONA@polito.it; Fax: ++39-11-564-4015
!C
!C                Christian Kummerow
!C                Laboratory for Atmospheres, NASA/Goddard Space Flight Center
!C                Code 912
!C                20716 Greenbelt, Maryland, USA 
!C                e-mail: kummerow@audry.gsfc.nasa.gov
                
!C     First version: September 1993
!C     Revised:       August 1994
!C
!C   Argument list:
!C   VAR         I/O   TYPE        DESCRIPTION
!C   ---         ---   ----        -----------
!C   MXCL         I    I           MAX NUM OF SUBCLOUDS IN X AND Y DIRECTION 
!C   MXLYR        I    I           MAX NUM OF  LAYERS 
!C   MXMO         I    I           MAX NUM OF PHASE FUNC. MOMENTS
!C   MXAN         I    I           MAX NUM OF ANGLES IN DEG. FOR THE PHASE FUNC. 
!C   MXEMI        I    I           NUM OF ANGLES IN DEG. FOR THE GROUND
!C 			          EMISSIVITY COMPUTATION.
!C   STEP         I    R           DISTANCE IN THE X AND Y DIRECTION BETWEEN THE 
!C                                 COMPUTATIONAL POINTS ON THE SUBCLOUD TOP 
!C                                 IN WHICH  THE BRIGHT. TEMP. IS COMPUTED
!C   AZIM         I    R           AZIMUTHAL VIEWING ANGLE IN DEGREE
!C   UMU          I    R           COSINE OF POLAR VIEWING ANGLE 
!C   FISOT        I    R           BACKGROUND TEMPERATURE (2.7 K)
!C   LAMBERT      I   LOG          SURFACE TYPE
!C   HENYEY       I   LOG          PHASE FUNCTION INPUT INFORMATION
!C   PHAMOM       I   LOG          PHASE FUNCTION INPUT INFORMATION
!C   NXCL         I    I           NUMBER OF CLOUDS IN X DIRECTIONS
!C   NYCL         I    I           NUMBER OF CLOUDS IN Y DIRECTIONS
!C   SIZE         I    R           SIDE OF EACH SQUARE CLOUD (KM)
!C   BTEMP        I    R           SURFACE TEMPERATURE [K]
!C   BTEMP2       I    R           BACKGROUND SURFACE TEMPERATURE [K]
!C   EMIS         I    R           SURFACE EMISSIVITY 
!C   EMIS2        I    R           BKGROUND SURFACE EMISSIVITY 
!C   EMI          I    R           SURFACE EMISSIVITY FOR DIFFERENT ANGLES
!C   EMI2         I    R           BKGROUND SURFACE EMISSIVITY FOR DIFF. ANGLES
!C   LYRTEMP      I    R           TEMPERATURE AT THE TOP OF THE LAYERS [K]
!C   KEXTTOT      I    R           EXTINCTION COEFF. (1/KM) 
!C   SALBTOT      I    R           SINGLE SCATTERING ALBEDO 
!C   ASYMTOT      I    R           ASYMMETRY  COEFF.   
!C   MOM          I    R           PHASE FUNCTION MOMENTS 
!C   ANG          I    R           ANGLES  IN DEG. FOR THE TABULATION 
!C                                 OF THE PHASE FUNCTION
!C   PHA          I    R           PHASE FUNCTION VALUES IN THE ABOVE ANGLES 
!C   NLYR         I    R           NUMBER OF LAYERS 
!C   FILEPHA      I   A*50         INPUT FILE CONTAINING THE PHASE FUNCTION
!C   FILEINP      I   A*50         IPUT FILE CONTAINING THE ATM. AND CLOUD PARAM.
!C   FILEOUT      O   A*50         OUTPUT FILE WITH THE BRIGHT. TEMPERATURES 
!C   WIND         I    R           WIND SPEED OVER WATER [m/s]
!C   FREQY        I    R           FREQUENCY [GHz]
!C   POLN         I   A*1          POLARIZATION  ('V'/'H')                       
!C   NEMI         I    I           NEMI=1 SURFACE EMISSIVITY IS COMPUTED
!C                                 FOR DIFFERENT POLAR ANGLES
!C				  NEMI=0 SURFACE EMISSIVITY IS READ FROM
!C                                 THE INPUT FILE
!C 
!C                 SUBROUTINES CALLED
!C   
!C  EMIT : computes the surface emissivity as a function of the polar
!C         angle, frequency, wind, surface temperature and polarization
!C
!C  PHANG: Computes the integral of the phase function with respect to
!C         the polar angle.
!C
!C  PHAMO: Computes the phase functions in a set of points, given its
!C         moments            
!C
!C  BT3  : Is the main subroutine that  performs the radiative transfer
!C         computations
!C
PARAMETER(MXLYR=180,MXCL=1,MXMO=1,MXAN=180,MXEMI=89,NPHO=100000)
LOGICAL  LAMBERT,PHAMOM,HENYEY
REAL  LYRTEMP(0:MXCL,0:MXCL,0:MXLYR),EMI(0:MXEMI),EMI2(0:MXEMI),  &
     EMIS,EMIS2, tskin(0:MXCL,0:MXCL), &
     LYRHGT(0:MXCL,0:MXCL,0:MXLYR),   &
     MOM(0:MXCL,0:MXCL,MXLYR,0:MXMO), &
     ANG(0:MXAN),PHA(0:MXCL,0:MXCL,MXLYR,0:MXAN), &
     KEXTTOT(0:MXCL,0:MXCL,MXLYR), &
     PHADI(0:MXCL,0:MXCL,MXLYR,0:MXAN), &
     ASYMTOT(0:MXCL,0:MXCL,MXLYR), SALBTOT(0:MXCL,0:MXCL,MXLYR) 
REAL  MEANFP(0:MXCL,0:MXCL,MXLYR), WIND, FREQY, EBAR, &
          B0(0:MXCL,0:MXCL,MXLYR), B1(0:MXCL,0:MXCL,MXLYR), &
          BTEMP(0:MXCL,0:MXCL)
CHARACTER*50 FILEINP,FILEOUT,FILEPHA
CHARACTER POLN*1
REAL RANG(0:MXAN),RCOS(0:MXAN), RPOL(0:MXMO,0:MXAN), &
     P1(0:MXAN),P2(0:MXAN),P3(0:MXMO,0:MXAN)
INTEGER NLYR(0:MXCL,0:MXCL)
real inc_angle_mwr, umu 
!inc_angle_mwr= 52.8
!C     
!C     
!C     Hard coded input variables that define the type of radiative 
!C     transfer computationn
!C     

contains
  subroutine init(freq)
    real :: freq
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
    FREQY=freq
    POLN='H'
    WIND=0.
    !c     JLH
    NEMI = 1
    
    
    !C     
    !C     READ THE CLOUD PARAMETERS 
    !C     
    !c   em   OPEN(UNIT=1,FILE=FILEINP,STATUS='OLD')
    !c      READ(1,*)NXCL,NYCL,SIZE,BTEMP2,BTEMP,EMIS2,EMIS
    NXCL=1
    NYCL=1
    
    
    
    !C     ******  Do the phase function if not HENYEY-GREENSTEIN *****
    !C
    IF (.NOT.HENYEY) THEN
       !C     
       !C     ************************************************************
       !C     ********  Phase function tabulated by points   **************
       !C     
       NANG = 180
       
       !C     
    ENDIF
    
    !C
    !C      ******************************************************************
    !C      ------------  Background atmospheric parameters  ----------------- 
    !C
    !c        READ(1,*)NX,NY,NLYR(NX,NY),LYRHGT(NX,NY,0),LYRTEMP(NX,NY,0)
    !c        DO  50  I = 1,NLYR(NX,NY)
    !c           READ(1,*) ILYR, LYRHGT(NX,NY,I), LYRTEMP(NX,NY,I), 
    !c     $          KEXTTOT(NX,NY,I),SALBTOT(NX,NY,I),ASYMTOT(NX,NY,I) 
    !c  50    CONTINUE
    !C
    !C       ---------------      Cloud parameters    ------------------------ 
    
    !C     
    
    !C
    !C       ************ GO TO RADIATIVE TRANSFER SUBROUTINE **************
    !C
  end subroutine init
end module MC3



