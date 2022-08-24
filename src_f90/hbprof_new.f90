module tableP2
  real :: zKuR(300)
  real :: zKaR(300)
  real :: zKaS(300)
  real :: zKaG(300)
  real :: zKaH(300)
  real :: zKuS(300)
  real :: zKuG(300)
  real :: zKuH(300)
  real :: dmR(300)
  real :: dmG(300)
  real :: dmS(300)
  real :: dmH(300)
  real :: dmBB(300)
  real :: zKuBB(300)
  real :: zKaBB(300)
  real :: attKaBB(300)
  real :: attKuBB(300)
  real :: pRateBB(300)
  real :: rainRate(300)
  real :: graupRate(300)
  real :: snowRate(300)
  real :: hailRate(300)
  real :: rwc(300)
  real :: pwcbb(300)
  real :: gwc(300)
  real :: swc(300)
  real :: hwc(300)
  real :: attKaR(300)
  real :: attKuR(300)
  real :: attKuS(300)
  real :: attKaS(300)
  real :: attKaG(300)
  real :: attKuG(300)
  real :: attKaH(300)
  real :: attKuH(300)
  real :: asymH(300,8)
  real :: asymR(300,8)
  real :: asymG(300,8)
  real :: asymS2(300,8)
  real :: asymBB(300,8)
  
  real :: kextR(300,8)
  real :: kextH(300,8)
  real :: kextS2(300,8)
  real :: kextG(300,8)
  real :: kextBB(300,8)
  
  real :: salbS2(300,8)
  real :: salbG(300,8)
  real :: salbR(300,8)
  real :: salbH(300,8)
  real :: salbBB(300,8)
  
  integer :: nJ
end module tableP2
subroutine initP2
  use tableP2
  use tables2
  implicit none
  integer :: i
  real :: f
  f=1
  kextH(1:289,1:8)=kextTableH(1:289,1:8,3)
  salbH(1:289,1:8)=salbTableH(1:289,1:8,3)
  asymH(1:289,1:8)=asymTableH(1:289,1:8,3)
  kextS2(1:289,1:8)=kextTableS2(1:289,1:8,3)
  salbS2(1:289,1:8)=salbTableS2(1:289,1:8,3)
  asymS2(1:289,1:8)=asymTableS2(1:289,1:8,3)
  kextR(1:289,1:8)=kextTable(1:289,1:8,3)
  salbR(1:289,1:8)=salbTable(1:289,1:8,3)
  asymR(1:289,1:8)=asymTable(1:289,1:8,3)
  kextG(1:289,1:8)=kextTableG(1:289,1:8,3)
  salbG(1:289,1:8)=salbTableG(1:289,1:8,3)
  asymG(1:289,1:8)=asymTableG(1:289,1:8,3)

  kextBB(1:289,1:8)=kextTableBB(1:289,1:8,3)
  salbBB(1:289,1:8)=salbTableBB(1:289,1:8,3)
  asymBB(1:289,1:8)=asymTableBB(1:289,1:8,3)

  print*, kextTable(80,5,1)
  do i=1,nbins
     zKuR(i)=zmin+(i-1)*dzbin
     zKaR(i)=z35Table(i,1)
     zKaG(i)=z35TableG(i,1)
     zKaS(i)=z35TableS2(i,1)
     zKuG(i)=z13TableG(i,1)
     zKuH(i)=z13TableH(i,1)
     zKaH(i)=z35TableH(i,1)
     zKuS(i)=z13TableS2(i,1)
     dmR(i)=d013Table(i,1)
     dmG(i)=d013TableG(i,1)
     dmH(i)=d013TableH(i,1)
     dmS(i)=d013TableS2(i,1)
     dmBB(i)=d013TableBB(i,1)
     attKaR(i)=att35Table(i,1)
     attKuR(i)=att13Table(i,1)
     attKaBB(i)=att35TableBB(i,1)
     attKuBB(i)=att13TableBB(i,1)
     attKaG(i)=att35TableG(i,1)
     attKaG(i)=att35TableG(i,1)
     attKuG(i)=att13TableG(i,1)
     attKuS(i)=att13TableS2(i,1)
     attKaS(i)=att35TableS2(i,1)
     attKuH(i)=att13TableH(i,1)
     attKaH(i)=att35TableH(i,1)
     rainRate(i)=pr13Table(i,1)
     graupRate(i)=pr13TableG(i,1)
     snowRate(i)=pr13TableS2(i,1)
     pRateBB(i)=pr13TableBB(i,1)
     hailRate(i)=pr13TableG(i,1)
     rwc(i)=10**pwc13Table(i,1)
     pwcbb(i)=10**pwc13TableBB(i,1)
     gwc(i)=10**pwc13TableG(i,1)
     swc(i)=10**pwc13TableS2(i,1)
     hwc(i)=10**pwc13TableH(i,1)
  end do
  !print*,f
  !stop
  nJ=nbins
  !print*, d013Table(200:240,1),nbins,nbinS2,nbinH
end subroutine initP2

