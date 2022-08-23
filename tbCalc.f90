subroutine test_n1(n1,binNodes,pwc,dm,sfcBin)
  integer :: n1, binNodes(n1,49,5), sfcBin(n1,49)
  real :: pwc(n1,49,5),dm(n1,49,5)
  print*, n1
end subroutine test_n1

subroutine interpol(xn,yn,fint,binC)
  implicit none
  integer :: xn(4)
  real :: yn(4), fint(88)
  integer :: binC
  integer :: i,n1
  real :: f
  fint=0
  do i=1,3
     !print*, xn(i),min(xn(i+1),binC-1),binC
     do n1=xn(i),min(xn(i+1),binC-1)
        !print*,n1
        f=(n1-xn(i))/(xn(i+1)-xn(i)+1e-3)
        fint(n1+1)=(1-f)*yn(i)+f*yn(i+1)
     end do
  end do
end subroutine interpol

subroutine interpol_kext(xn,yn,kext_int,binC)
  implicit none
  integer :: xn(10)
  real :: yn(10), kext_int(73)
  integer :: binC
  integer :: i,n1
  real :: f
  kext_int=0
  do i=1,9
     do n1=xn(i),min(xn(i+1),binC-1)
        f=(n1-xn(i))/(xn(i+1)-xn(i)+1e-3)
        kext_int(n1+1-xn(1))=(1-f)*yn(i)+f*yn(i+1)
     end do
  end do
end subroutine interpol_kext

     !bins=np.arange(envNode[i,j,0],envNode[i,j,-1]+1)
!kextInt=np.interp(bins,envNode[i,j,:],kextL)

subroutine calc_tb_f90(n1,binNodes,pwc,dm,sfcBin,pType,envNodes,qv,airTemp,press,&
     cldw,umu,sfcTemp,emiss,tbout)
  implicit none
  integer :: n1,binNodes(n1,49,5),sfcBin(n1,49)
  integer :: envNodes(n1,49,10)
  real :: qv(n1,49,10),airTemp(n1,49,10),press(n1,49,10),emiss(n1,49,13),sfcTemp(n1,49),umu
  real :: pwc(n1,49,88),dm(n1,49,88), cldw(n1,49,88)
  real :: tbout(n1,49,13)
  integer :: xn(4)
  integer :: pType(n1,49)
  integer :: i0,j0
  real :: yn(4), fint(88)
  real :: kextH1(88,8),salbH1(88,8), asymH1(88,8), kext_int(73)
  !print*, minval(airTemp)
  !print*, maxval(qv)
  !print*, minval(press)
  !return
  do i0=1,n1
     do j0=1,49
        xn=(/0,binNodes(i0,j0,2),binNodes(i0,j0,4),&
             binNodes(i0,j0,5)/)
        !print*, xn
        yn=(/0,0,1,1/)
        call interpol(xn,yn,fint,sfcBin(i0,j0))
        !fint=np.interp(range(88),xn,[0,0,1,1])
        if(pType(i0,j0)>0) then
           call calcz(n1,pwc,dm,i0,j0,fint,sfcBin,binNodes,&
                kextH1,salbH1,asymH1)
        else
           kextH1=0
           salbH1=0
           asymH1=0
        endif
        call rte(n1,binNodes,sfcBin,dm,kextH1,asymH1,salbH1,i0,j0,emiss,qv,airTemp, &
             press,envNodes,sfcTemp,cldw,umu,tbout)
        !    tb13=rte(binNodes,dm,pwc,kextH1,asymH1,salbH1,i0,j0,sfcEmiss,qv,airTemp,\
        !                 press,envNodes,skTemp,cldw,umu,cAlg)
        !    swp=((1-fint)*pwc[i0,j0,:]).sum()
        !    tb2d[i0,j0,:]=tb13
     end do
  enddo
end subroutine calc_tb_f90

subroutine calcz(n1,pwc_l,dm,i0,j0,fint,sfcBin,binNodes,&
     kextH1,salbH1,asymH1)
  use tablep2
  implicit none
  integer :: i0,j0,k, n1
  real :: fint(88), pwc_l(n1,49,88), dm(n1,49,88)
  integer :: sfcbin(n1,49), binNodes(n1,49,5)
  real :: zKus_l(88), zKas_l(88), zKur_l(88), zKar_l(88)
  real :: zKu1(88),zKa(88), pRate(88), att_Ka(88), att_Ku(88)
  real :: swc_l(88), rwc_l(88)
  real :: swc_bin, rwc_bin, dn, dr, pia_ka, pia_ku
  real :: kextH1(88,8),salbH1(88,8),asymH1(88,8)
  integer :: ibin, n, j
  
  swc_l=pwc_l(i0,j0,:)*(1-fint)
  rwc_l=pwc_l(i0,j0,:)*fint
  
  zKus_l=-99
  zKas_l=-99
  zKur_l=-99
  zKar_l=-99
  pRate=0

  att_Ka=0
  att_Ku=0
  kextH1=0
  salbH1=0
  asymH1=0
  
  dr=0.25
  n=88
  pia_ka=0
  do k=binNodes(i0,j0,1),binNodes(i0,j0,5)
     if (swc_l(k)>1e-3) then
        call bisection2(dms(1:253),253,dm(i0,j0,k),ibin)
        swc_bin=swc(ibin)
        dn=log10(swc_l(k)/swc_bin)
        zKus_l(k)=zkus(ibin)+10*dn
        zKas_l(k)=zkas(ibin)+10*dn
        pRate(k)=snowrate(ibin)*10**dn
        att_ka(k)=attkas(ibin)*10**dn
        kextH1(k,:)=kextH1(k,:)+kexts2(ibin,:)*10**dn
        salbH1(k,:)=salbH1(k,:)+kexts2(ibin,:)*10**dn*&
             salbs2(ibin,:)
        asymH1(k,:)=asymH1(k,:)+kexts2(ibin,:)*10**dn*&
             salbs2(ibin,:)*asyms2(ibin,:)
     end if
  enddo
  zKur_l=-99
  zKar_l=-99
  do k=binNodes(i0,j0,1),binNodes(i0,j0,5)
     if (rwc_l(k)>1e-3) then
        call bisection2(dmr(1:289),289,dm(i0,j0,k),ibin)
        rwc_bin=rwc(ibin)
        dn=log10(rwc_l(k)/rwc_bin)
        zKur_l(k)=zkur(ibin)+10*dn
        zKar_l(k)=zkar(ibin)+10*dn
        pRate(k)=pRate(k)+rainrate(ibin)*10**dn
        att_ka(k)=att_ka(k)+attkar(ibin)*10**dn
        kextH1(k,:)=kextH1(k,:)+kextr(ibin,:)*10**dn
        salbH1(k,:)=salbH1(k,:)+kextr(ibin,:)*10**dn*&
             salbr(ibin,:)
        asymH1(k,:)=asymH1(k,:)+kextr(ibin,:)*10**dn*&
             salbr(ibin,:)*asymr(ibin,:)
     endif
  end do
  do k=binNodes(i0,j0,1),binNodes(i0,j0,5)
     do j=1,8
        if (kextH1(k,j)>1e-5) then
           asymH1(k,j)=asymH1(k,j)/(salbH1(k,j)+1e-7)
           salbH1(k,j)=salbH1(k,j)/kextH1(k,j)
        endif
     enddo
 enddo
 
 do k=1,8
    kextH1(binNodes(i0,j0,5):sfcBin(i0,j0),k)=kextH1(binNodes(i0,j0,5),k)
    asymH1(binNodes(i0,j0,5):sfcBin(i0,j0),k)=asymH1(binNodes(i0,j0,5),k)
    salbH1(binNodes(i0,j0,5):sfcBin(i0,j0),k)=salbH1(binNodes(i0,j0,5),k)
 enddo
 !zKu1D=np.log10(10**(0.1*zKur)+10**(0.1*zKus))*10
 !   zKa1D=np.log10(10**(0.1*zKar)+10**(0.1*zKas))*10
 !   calc_pia(att_ka,pia_ka,dr,n)
 !   zKa1D-=pia_ka

 !return zKu1D,zKa1D,pRate,kextH1,asymH1,salbH1
end subroutine calcz


subroutine rte(n1,binNodes,sfcBin,dm,kextH1,asymH1,salbH1,i,j,emiss,qv,airTemp, &
  press,envNode,sfcTemp,cldw,umu,tbout)
  implicit none
  real ::  freqs(13)
  integer :: nfreq,n1, i, j, k, pol, i_freq
  integer :: npol(13), iFreq(15)
  logical :: lambert
  integer :: envNode(n1,49,10), binNodes(n1,49,5), sfcBin(n1,49)
  real :: dm(n1,49,88), pRate(n1,49,88), kextH1(88,8), asymH1(88,8), salbH1(88,8), &
       airTemp(n1,49,10), qv(n1,49,10), press(n1,49,10), &
       cldw(n1,49,88), emiss(n1,49,13),sfcTemp(n1,49)
  real :: q, t, pa, cld, z_clw, absair, abswv
  real :: kextL(10), emis, ebar, umu, tk, f, kext_int(73)
  real :: kext(73),salb(73),asym(73)
  real :: tbout(n1,49,13), ireturn
  
  freqs=(/10.6,10.6,18.7,18.7,23.,37.,37.,89.,89.,166.,166.,186.3,190.3/)
  npol=(/1,0,1,0,1,1,0,1,0,1,0,1,1/)
  iFreq=(/0,0,1,1,2,3,3,4,4,5,5,6,6,7,7/)
  nfreq=8

  do i_freq=1,13
     pol=npol(i_freq)
     f=freqs(i_freq)
     do k=1,10
        q=qv(i,j,k)
        tk=airTemp(i,j,k)
        pa=press(i,j,k)
        cld=cldw(i,j,k)
        if (q<1e-3) q=1e-3
        call gasabsr98(f,tk,q*1e-3,pa*1e2,absair,abswv,ireturn)
        if(cld>0.001) then
           call gcloud(f,tk,cld,z_clw)
        else
           z_clw=0.
        endif
        kextL(k)=(absair+abswv+z_clw)
     enddo

     call interpol_kext(envNode,kextL,kext_int,sfcBin(i,j))
     kext=kext_int
     kext=kext+kextH1(envNode(i,j,1):envNode(i,j,5),i_freq)
     salb=salbH1(envNode(i,j,1):envNode(i,j,5),i_freq)* &
          (1-kext_int/kext_int)
     asym=asymH1(envNode(i,j,1):envNode(i,j,5),i_freq)
     !tLayer=list(np.interp(bins,envNode[i,j,:],airTemp[i,j,:]))
     !tLayer.append(sfcTemp[i,j])
     !iEnum+=1
     emis=emiss(i,j,i_freq)
     ebar=emis
     !nL=kextInt.shape[0]
     lambert=.False.

     !tb1=cAlg.radtran(umu,sfcTemp[i,j],tLayer[::-1],\
     !                        np.arange(nL+1)*0.25,kextInt[::-1],\
     !                        salb[::-1],asym[::-1],2.7,emis,ebar,lambert)
     !except:
  end do
end subroutine rte
