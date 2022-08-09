subroutine mcrt(nz1,kext1d,salb1d,asym1d,lyrtemp1d,&
     lyrhgt1d,emis_in,tskin_in,tb_out,inc_angle)
  use mc3
  implicit none
  integer :: nz1
  real :: kext1d(nz1),salb1d(nz1),asym1d(nz1),lyrtemp1d(0:nz1),&
       lyrhgt1d(0:nz1)
  real :: emis_in, tskin_in
  real,intent(out) :: tb_out
  integer :: nxcl, nycl, nz, nemi, nang, i
  real :: tsub, step, size, fisot, btemp2, azim
  real :: inc_angle
  
  AZIM = 180.
  umu = cos(inc_angle*3.141593/180.)
  FISOT =  2.7
  LAMBERT = .FALSE.
  STEP =  500.          
  PHAMOM =.FALSE.
  HENYEY = .TRUE.
  POLN='H'
  NEMI = 1

  emi=emis_in
  emis=emis_in
  emi2=emi
  !print*, lyrtemp(1,1,1:nz)
  
  emis2=emis
  tskin(1,1)=tskin_in
  nz=nz1
  kexttot(1,1,1:nz1)=kext1d
  salbtot(1,1,1:nz1)=salb1d
  asymtot(1,1,1:nz)=asym1d
  do i=1,nz
     kexttot(0,0,i)=kexttot(1,1,i)
     asymtot(0,0,i)=asymtot(1,1,i)
     salbtot(0,0,i)=salbtot(1,1,i)
  enddo
  do i=0,nz
     lyrhgt(1,1,i)=lyrhgt1d(i)
     lyrtemp(1,1,i)=lyrtemp1d(i)
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
  CALL BT3( SIZE, NLYR, UMU, AZIM, LAMBERT, FISOT, &
       BTEMP2, EMIS2, BTEMP, LYRTEMP, LYRHGT, KEXTTOT, &
       SALBTOT, ASYMTOT,EMIS, STEP, NXCL, NYCL,FILEOUT, &
       HENYEY,NANG,ANG,PHADI,EMI,EMI2,NEMI, &
       MEANFP,B0,B1,MXLYR,MXCL,MXAN,MXEMI,NPHO,tsub)
  tb_out=tsub
end subroutine mcrt
