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
     do n1=xn(i),min(xn(i+1),87)
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
     do n1=xn(i),min(xn(i+1),87)
        f=(n1-xn(i))/(xn(i+1)-xn(i)+1e-3)
        kext_int(n1+1-xn(1))=(1-f)*yn(i)+f*yn(i+1)
     end do
  end do
end subroutine interpol_kext

     !bins=np.arange(envNode[i,j,0],envNode[i,j,-1]+1)
!kextInt=np.interp(bins,envNode[i,j,:],kextL)

subroutine calc_tb_f90(n1,binNodes,pwc,dm,sfcBin,pType,envNodes,qv,airTemp,press,&
     cldw,umu,sfcTemp,emiss,bbPeak,tbout)
  implicit none
  integer :: n1,binNodes(n1,49,5),sfcBin(n1,49)
  integer :: envNodes(n1,49,10)
  real :: qv(n1,49,10),airTemp(n1,49,10),press(n1,49,10),emiss(n1,49,13),sfcTemp(n1,49),umu
  real :: pwc(n1,49,88),dm(n1,49,88), cldw(n1,49,88)
  real,intent(out) :: tbout(n1,49,13)
  integer :: xn(4)
  integer :: pType(n1,49),bbPeak(n1,49)
  integer :: i0,j0
  real :: yn(4), fint(88)
  real :: kextH1(88,8),salbH1(88,8), asymH1(88,8), kext_int(73)
  real :: zKu(n1,49,88), zKa(n1,49,88)
  real :: tb13(13)
  real :: pwc1(88),dm1(88),emiss1(13),qv1(10),airTemp1(10),&
       press1(10),sfcTemp1,cldw1(88),zKu1(88),zKa1(88)
  integer :: sfcBin1,binNodes1(5),envNodes1(10),pType1,bbPeak1

  !print*, minval(airTemp)
  !print*, maxval(qv)
  !print*, minval(press)
  !return
  integer omp_get_max_threads, nthreads
  !nthreads=omp_get_max_threads()
  !print*, nthreads
 
  do i0=1,n1
     !$OMP PARALLEL DO PRIVATE(pwc1,dm1,sfcBin1,binNodes1,envNodes1,pType1,emiss1,&
     !$OMP& qv1,airTemp1,press1,sfcTemp1,cldw1,bbPeak1,tb13,j0,zKu1,zKa1)
     do j0=1,49
        pwc1=pwc(i0,j0,:)
        dm1=dm(i0,j0,:)
        sfcBin1=sfcBin(i0,j0)
        binNodes1=binNodes(i0,j0,:)
        envNodes1=envNodes(i0,j0,:)
        pType1=pType(i0,j0)
        emiss1=emiss(i0,j0,:)
        qv1=qv(i0,j0,:)
        airTemp1=airTemp(i0,j0,:)
        press1=press(i0,j0,:)
        sfcTemp1=sfcTemp(i0,j0)
        cldw1=cldw(i0,j0,:)
        bbPeak1=bbPeak(i0,j0)
        call rte1d(pwc1,dm1,sfcBin1,&
             binNodes1,envNodes1,pType1,&
             emiss1,qv1,airTemp1,&
             press1,sfcTemp1,cldw1,umu,bbPeak1,&
             tb13,zKu1,zKa1)
        tbout(i0,j0,:)=tb13
     end do
     !$OMP END PARALLEL DO
  enddo

end subroutine calc_tb_f90


subroutine calc_tb_f90_flattened(n1,z_obs,binNodes,pwc,dm,sfcBin,pType,envNodes,qv,airTemp,press,&
     cldw,umu,sfcTemp,emiss,bbPeak,tbout,tbout_dfr,pRate_out, pwc_out, swp)
  implicit none
  integer :: n1,binNodes(n1,5),sfcBin(n1)
  real :: z_obs(n1,88,2)
  integer :: envNodes(n1,10),bbPeak(n1)
  real :: qv(n1,10),airTemp(n1,10),press(n1,10),emiss(n1,13),sfcTemp(n1),umu
  real :: pwc(n1,88),dm(n1,88), cldw(n1,88)
  real,intent(out) :: tbout(n1,13),tbout_dfr(n1,13),pRate_out(n1,88), pwc_out(n1,88), swp(n1,2)
  integer :: xn(4)
  integer :: pType(n1)
  integer :: i0,j0
  real :: yn(4), fint(88)
  real :: kextH1(88,8),salbH1(88,8), asymH1(88,8), kext_int(73)
  real :: zKu1(n1,88), zKa1(n1,88)
  !print*, minval(airTemp)
  !print*, maxval(qv)
  !print*, minval(press)
  !return
  !integer omp_get_max_threads, nthreads
  !nthreads=omp_get_max_threads()
  !print*, nthreads
  !$OMP PARALLEL
  !$OMP DO
  do i0=1,n1
     call rte1d(pwc(i0,:),dm(i0,:),sfcBin(i0),&
          binNodes(i0,:),envNodes(i0,:),pType(i0),&
          emiss(i0,:),qv(i0,:),airTemp(i0,:),&
          press(i0,:),sfcTemp(i0),cldw(i0,:),umu,bbPeak(i0),tbout(i0,:),zKu1(i0,:),zKa1(i0,:))
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !$OMP PARALLEL
  !$OMP DO
  do i0=1,n1
     call rte1d_dfr(z_obs(i0,:,:),pwc(i0,:),dm(i0,:),sfcBin(i0),&
          binNodes(i0,:),envNodes(i0,:),pType(i0),&
          emiss(i0,:),qv(i0,:),airTemp(i0,:),&
          press(i0,:),sfcTemp(i0),cldw(i0,:),umu,bbPeak(i0),tbout_dfr(i0,:),zKu1(i0,:),zKa1(i0,:),&
          pRate_out(i0,:),swp(i0,:))
     
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  pwc_out=pwc
end subroutine calc_tb_f90_flattened

subroutine rte1d(pwc1,dm1,sfcBin1,binNodes1,envNodes1,pType1,emiss1,qv1,airTemp1,&
     press1,sfcTemp1,cldw1,umu,bbPeak,tbout_13,zKu1,zKa1)
  implicit none
  real ::  pwc1(88),dm1(88), airTemp1(10), qv1(10), press1(10),sfcTemp1,cldw1(88)
  real :: yn(4), fint(88)
  integer :: xn(4), bbPeak
  real :: kextH1(88,8),salbH1(88,8), asymH1(88,8), kext_int(73)
  integer :: envNodes1(10), binNodes1(5), sfcBin1, pType1
  real :: tbout_13(13), zKu1(88), zKa1(88)
  real :: emiss1(13), umu
  
  xn=(/0,binNodes1(2),binNodes1(4),binNodes1(5)/)
  yn=(/0,0,1,1/)
  call interpol(xn,yn,fint,sfcBin1)
  !fint=0
  if(pType1>0) then
     call calcz1(pwc1,dm1,fint,sfcBin1,binNodes1,bbPeak,&
          kextH1,salbH1,asymH1,zKu1,zKa1)
  else
     kextH1=0
     salbH1=0
     asymH1=0
  endif
  call rte1(binNodes1,sfcBin1,dm1,kextH1,asymH1,salbH1,emiss1,qv1,airTemp1, &
       press1,envNodes1,sfcTemp1,cldw1,umu,tbout_13)
end subroutine rte1d


subroutine rte1d_dfr(z_obs,pwc1,dm1,sfcBin1,binNodes1,envNodes1,pType1,emiss1,qv1,airTemp1,&
     press1,sfcTemp1,cldw1,umu,bbPeak,tbout_13,zKu1,zKa1,prate_out,swp)
  implicit none
  real ::  pwc1(88),dm1(88), airTemp1(10), qv1(10), press1(10),sfcTemp1,cldw1(88)
  real :: yn(4), fint(88), prate_out(88)
  integer :: xn(4), bbPeak
  real :: kextH1(88,8),salbH1(88,8), asymH1(88,8), kext_int(73)
  integer :: envNodes1(10), binNodes1(5), sfcBin1, pType1
  real :: tbout_13(13), zKu1(88), zKa1(88)
  real :: emiss1(13), umu, z_obs(88,2),swp(2)
  
  xn=(/0,binNodes1(2),binNodes1(4),binNodes1(5)/)
  yn=(/0,0,1,1/)
  call interpol(xn,yn,fint,sfcBin1)
  !fint=0
  if(pType1>0) then
     call calcz1_dfr(z_obs,pwc1,dm1,fint,sfcBin1,binNodes1,bbPeak,&
          kextH1,salbH1,asymH1,zKu1,zKa1,prate_out,swp)
  else
     kextH1=0
     salbH1=0
     asymH1=0
  endif
  call rte1(binNodes1,sfcBin1,dm1,kextH1,asymH1,salbH1,emiss1,qv1,airTemp1, &
       press1,envNodes1,sfcTemp1,cldw1,umu,tbout_13)
end subroutine rte1d_dfr


subroutine rte1d_o(pwc1,dm1,sfcBin1,binNodes1,envNodes1,pType1,emiss1,qv1,airTemp1,&
     press1,sfcTemp1,cldw1,umu,bbPeak,tbout_13,zKu1,zKa1)
  implicit none
  real ::  pwc1(88),dm1(88), airTemp1(10), qv1(10), press1(10),sfcTemp1,cldw1(88)
  real :: yn(4), fint(88)
  integer :: xn(4), bbPeak
  real :: kextH1(88,8),salbH1(88,8), asymH1(88,8), kext_int(73)
  integer :: envNodes1(10), binNodes1(5), sfcBin1, pType1
  real,intent(out) :: tbout_13(13), zKu1(88), zKa1(88)
  real :: emiss1(13), umu
  
  xn=(/0,binNodes1(2),binNodes1(4),binNodes1(5)/)
  yn=(/0,0,1,1/)
  call interpol(xn,yn,fint,sfcBin1)
  if(pType1>0) then
     call calcz1(pwc1,dm1,fint,sfcBin1,binNodes1,bbPeak,&
          kextH1,salbH1,asymH1,zKu1,zKa1)
  else
     kextH1=0
     salbH1=0
     asymH1=0
  endif
  call rte1(binNodes1,sfcBin1,dm1,kextH1,asymH1,salbH1,emiss1,qv1,airTemp1, &
       press1,envNodes1,sfcTemp1,cldw1,umu,tbout_13)
end subroutine rte1d_o


subroutine calcz1(pwc1,dm1,fint,sfcBin1,binNodes1,bbPeak,&
     kextH1,salbH1,asymH1,zKu1,zKa1)
  use tablep2
  implicit none
  integer :: k
  real :: fint(88), pwc1(88), dm1(88)
  integer :: sfcbin1, binNodes1(5),bbPeak
  real :: zKus_l(88), zKas_l(88), zKur_l(88), zKar_l(88)
  real :: pRate(88), att_Ka(88), att_Ku(88)
  real :: swc_l(88), rwc_l(88)
  real :: swc_bin, rwc_bin, dn, dr, pia_ka, pia_ku
  real,intent(out) :: kextH1(88,8),salbH1(88,8),asymH1(88,8),zKu1(88),zKa1(88)
  integer :: ibin, n, j, ik
  real :: fract, pwc_bb
  
  swc_l=pwc1(:)*(1-fint)
  rwc_l=pwc1(:)*fint
  
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

  do k=max(binNodes1(1)+1,1),min(1+binNodes1(5),88)
     if (swc_l(k)>1e-3) then
        if(dm1(k)<dms(1)) then
           ibin=1
        else
           call bisection2(dms(1:253),253,dm1(k),ibin)
        end if
        ibin=min(253,ibin)
        ibin=max(1,ibin)
        !print*, ibin, 253, k
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
  !kextH1=0
  !salbH1=0
  !asymH1=0
  !return
  zKur_l=-99
  zKar_l=-99
  
  do k=max(binNodes1(1)+1,1),min(1+binNodes1(5),88)
     if (rwc_l(k)>1e-3) then
        if(dm1(k)<dmr(1)) then
           ibin=1
        else
           call bisection2(dmr(1:289),289,dm1(k),ibin)
        endif
        ibin=min(289,ibin)
        ibin=max(1,ibin)
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
  do k=max(1+binNodes1(1),1),min(1+binNodes1(5),88)
     do j=1,8
        if (kextH1(k,j)>1e-9) then
           asymH1(k,j)=asymH1(k,j)/(salbH1(k,j)+1e-9)
           salbH1(k,j)=salbH1(k,j)/kextH1(k,j)
        endif
     enddo
  enddo
  if (binNodes1(5)+1>88) return
  do ik=max(binNodes1(5)+1,1),min(sfcBin1+1,88)
     do k=1,8
        kextH1(ik,k)=kextH1(binNodes1(5)+1,k)
        asymH1(ik,k)=asymH1(binNodes1(5)+1,k)
        salbH1(ik,k)=salbH1(binNodes1(5)+1,k)
     enddo
  enddo
  if(bbPeak>0.and.bbPeak<88) then
     do k=bbPeak-3,bbPeak
        fract=(bbPeak-k)/4.0
        call bisection2(dmbb(1:289),289,dm1(k),ibin)
        pwc_bb=pwcbb(ibin)
        if(pwc1(k)>0) then
           dn=log10(pwc1(k)/pwc_bb)
           kextH1(k+1,:)=kextH1(k+1,:)*(1-fract)+fract*kextbb(ibin,:)*10**dn
           salbH1(k+1,:)=salbH1(k+1,:)*(1-fract)+fract*salbbb(ibin,:)
           asymH1(k+1,:)=asymH1(k+1,:)*(1-fract)+fract*asymbb(ibin,:)
        endif
     enddo
     do k=bbPeak+1,bbPeak+3
        if(k<88) then
           fract=1-(k-bbPeak)/3.
           call bisection2(dmbb(1:289),289,dm1(k),ibin)
           pwc_bb=pwcbb(ibin)
           if(pwc1(k)>0) then
              dn=log10(pwc1(k)/pwc_bb)
              kextH1(k+1,:)=kextH1(k+1,:)*(1-fract)+fract*kextbb(ibin,:)*10**dn
              salbH1(k+1,:)=salbH1(k+1,:)*(1-fract)+fract*salbbb(ibin,:)
              asymH1(k+1,:)=asymH1(k+1,:)*(1-fract)+fract*asymbb(ibin,:)
           endif
        endif
     enddo
  endif
  do k=1,88
     zKu1(k)=10*log10(10.0**(0.1*zKur_l(k))+10.0**(0.1*zKus_l(k)))
  enddo
end subroutine calcz1


subroutine calcz1_dfr(z_obs,pwc1,dm1,fint,sfcBin1,binNodes1,bbPeak,&
     kextH1,salbH1,asymH1,zKu1,zKa1,pRate_out,swp)
  use tablep2
  implicit none
  integer :: k
  real :: fint(88), pwc1(88), dm1(88), z_obs(88,2)
  integer :: sfcbin1, binNodes1(5),bbPeak
  real :: zKus_l(88), zKas_l(88), zKur_l(88), zKar_l(88)
  real :: pRate(88), att_Ka(88), att_Ku(88)
  real :: swc_l(88), rwc_l(88)
  real :: swc_bin, rwc_bin, dn, dr, pia_ka, pia_ku
  real,intent(out) :: kextH1(88,8),salbH1(88,8),asymH1(88,8),zKu1(88),zKa1(88), pRate_out(88)
  integer :: ibin, n, j, ik
  real :: fract, pwc_bb, dwr1
  real, intent(out) :: swp(2)
  
  swc_l=pwc1(:)*(1-fint)
  rwc_l=pwc1(:)*fint
  
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
  swp=0
  do k=max(binNodes1(1)+1,1),min(1+binNodes1(5),88)
     if (swc_l(k)>1e-3) then
        swp(1)=swp(1)+swc_l(k)
        if(k<binNodes1(3)+1-3)then
           dwr1=z_obs(k,1)-z_obs(k,2)
           !print*, dwr1
           if(dwr1>0.5 .and. dwr1<9 .and. z_obs(k,1)>12) then
              !dwr=(cAlg.tablep2.zkus-cAlg.tablep2.zkas)[:253]
              call bisection2(dwr(1:253),253,dwr1,ibin)
              dn=(z_obs(k,1)-zkus(ibin))/10.
              !print*, swc(ibin)*10**dn, swc_l(k)
              if(swc_l(k)>0.98*pwc1(k))then
                 pwc1(k)=swc(ibin)*10**dn
              end if
              swc_l(k)=swc(ibin)*10**dn
              !print*, dm1(k), dms(ibin)
              dm1(k)=dms(ibin)
           else
              !print*, dwr1, 'out of range'
           endif
        endif
        swp(2)=swp(2)+swc_l(k)
        if(dm1(k)<dms(1)) then
           ibin=1
        else
           call bisection2(dms(1:253),253,dm1(k),ibin)
        end if
        ibin=min(253,ibin)
        ibin=max(1,ibin)
        !print*, ibin, 253, k
        swc_bin=swc(ibin)
        dn=log10(swc_l(k)/swc_bin)
        zKus_l(k)=zkus(ibin)+10*dn
        zKas_l(k)=zkas(ibin)+10*dn
        pRate_out(k)=snowrate(ibin)*10**dn
        att_ka(k)=attkas(ibin)*10**dn
        kextH1(k,:)=kextH1(k,:)+kexts2(ibin,:)*10**dn
        salbH1(k,:)=salbH1(k,:)+kexts2(ibin,:)*10**dn*&
             salbs2(ibin,:)
        asymH1(k,:)=asymH1(k,:)+kexts2(ibin,:)*10**dn*&
             salbs2(ibin,:)*asyms2(ibin,:)
     end if
  enddo
  !kextH1=0
  !salbH1=0
  !asymH1=0
  !return
  zKur_l=-99
  zKar_l=-99
  
  do k=max(binNodes1(1)+1,1),min(1+binNodes1(5),88)
     if (rwc_l(k)>1e-3) then
        if(dm1(k)<dmr(1)) then
           ibin=1
        else
           call bisection2(dmr(1:289),289,dm1(k),ibin)
        endif
        ibin=min(289,ibin)
        ibin=max(1,ibin)
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
  do k=max(1+binNodes1(1),1),min(1+binNodes1(5),88)
     do j=1,8
        if (kextH1(k,j)>1e-9) then
           asymH1(k,j)=asymH1(k,j)/(salbH1(k,j)+1e-9)
           salbH1(k,j)=salbH1(k,j)/kextH1(k,j)
        endif
     enddo
  enddo
  if (binNodes1(5)+1>88) return
  do ik=max(binNodes1(5)+1,1),min(sfcBin1+1,88)
     do k=1,8
        kextH1(ik,k)=kextH1(binNodes1(5)+1,k)
        asymH1(ik,k)=asymH1(binNodes1(5)+1,k)
        salbH1(ik,k)=salbH1(binNodes1(5)+1,k)
     enddo
  enddo
  if(bbPeak>0.and.bbPeak<88) then
     do k=bbPeak-3,bbPeak
        fract=1-(bbPeak-k)/4.0
        call bisection2(dmbb(1:289),289,dm1(k),ibin)
        pwc_bb=pwcbb(ibin)
        if(pwc1(k)>0) then
           dn=log10(pwc1(k)/pwc_bb)
           kextH1(k+1,:)=kextH1(k+1,:)*(1-fract)+fract*kextbb(ibin,:)*10**dn
           salbH1(k+1,:)=salbH1(k+1,:)*(1-fract)+fract*salbbb(ibin,:)
           asymH1(k+1,:)=asymH1(k+1,:)*(1-fract)+fract*asymbb(ibin,:)
        endif
     enddo
     do k=bbPeak+1,bbPeak+3
        if(k<88) then
           fract=1-(k-bbPeak)/3.
           call bisection2(dmbb(1:289),289,dm1(k),ibin)
           pwc_bb=pwcbb(ibin)
           if(pwc1(k)>0) then
              dn=log10(pwc1(k)/pwc_bb)
              kextH1(k+1,:)=kextH1(k+1,:)*(1-fract)+fract*kextbb(ibin,:)*10**dn
              salbH1(k+1,:)=salbH1(k+1,:)*(1-fract)+fract*salbbb(ibin,:)
              asymH1(k+1,:)=asymH1(k+1,:)*(1-fract)+fract*asymbb(ibin,:)
           endif
        endif
     enddo
  endif
  do k=1,88
     zKu1(k)=10*log10(10.0**(0.1*zKur_l(k))+10.0**(0.1*zKus_l(k)))
  enddo
end subroutine calcz1_dfr


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
  integer :: ibin, n, j, ik
  
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
  do k=binNodes(i0,j0,1)+1,min(1+binNodes(i0,j0,5),88)
     if (swc_l(k)>1e-3) then
        if(dm(i0,j0,k)<dms(1)) then
           ibin=1
        else
           call bisection2(dms(1:253),253,dm(i0,j0,k),ibin)
        end if
        !print*, ibin, 253, k
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
  do k=binNodes(i0,j0,1)+1,min(1+binNodes(i0,j0,5),88)
     if (rwc_l(k)>1e-3) then
        if(dm(i0,j0,k)<dmr(1)) then
           ibin=1
        else
           call bisection2(dmr(1:289),289,dm(i0,j0,k),ibin)
        endif
        !print*, ibin, 289, k
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
  do k=1+binNodes(i0,j0,1),min(1+binNodes(i0,j0,5),88)
     do j=1,8
        if (kextH1(k,j)>1e-5) then
           asymH1(k,j)=asymH1(k,j)/(salbH1(k,j)+1e-7)
           salbH1(k,j)=salbH1(k,j)/kextH1(k,j)
        endif
     enddo
  enddo
  if (binNodes(i0,j0,5)+1>88) return
  do ik=binNodes(i0,j0,5)+1,min(sfcBin(i0,j0)+1,88)
     do k=1,8
        kextH1(ik,k)=kextH1(binNodes(i0,j0,5)+1,k)
        asymH1(ik,k)=asymH1(binNodes(i0,j0,5)+1,k)
        salbH1(ik,k)=salbH1(binNodes(i0,j0,5)+1,k)
     enddo
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
  integer :: npol(13), iFreq(13)
  logical :: lambert
  integer :: envNode(n1,49,10), binNodes(n1,49,5), sfcBin(n1,49)
  real :: dm(n1,49,88), pRate(n1,49,88), kextH1(88,8), asymH1(88,8), salbH1(88,8), &
       airTemp(n1,49,10), qv(n1,49,10), press(n1,49,10), &
       cldw(n1,49,88), emiss(n1,49,13),sfcTemp(n1,49)
  real :: q, t, pa, cld, z_clw, absair, abswv
  real :: kextL(10), emis, ebar, umu, tk, f, kext_int(73)
  real :: kext(73),salb(73),asym(73)
  real,intent(out) :: tbout(n1,49,13)
  integer :: ireturn
  real :: lyrhgt(74), lyrtemp(74), kext_rev(73), asym_rev(73), salb_rev(73), tb, fisot, temp_int(73)
  real :: dr
  integer :: nlyr
  
  freqs=(/10.6,10.6,18.7,18.7,23.,37.,37.,89.,89.,166.,166.,186.3,190.3/)
  npol=(/1,0,1,0,1,1,0,1,0,1,0,1,1/)
  iFreq=(/0,0,1,1,2,3,3,4,4,5,5,6,7/)
  nfreq=8
  dr=0.25
  !print*, envNode(i,j,:)
  !print*, airTemp(i,j,:)
  call interpol_kext(envNode(i,j,:),airTemp(i,j,:),temp_int,sfcBin(i,j))
  !print*, temp_int
  salb=0
  asym=0
  !print*, temp_int
  !print*, envNode(i,j,:)
  !print*, airTemp(i,j,:)
  !print*, sfcBin(i,j), '*************'
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

     call interpol_kext(envNode(i,j,:),kextL,kext_int,sfcBin(i,j))
     kext=kext_int
     do k=1,73
        if(k+envNode(i,j,1)<88) then
           kext(k)=kext(k)+kextH1(k+envNode(i,j,1),iFreq(i_freq)+1)
           salb(k)=salbH1(k+envNode(i,j,1),iFreq(i_freq)+1)* &
                (1-kext_int(k)/kext(k))
           asym(k)=asymH1(k+envNode(i,j,1),iFreq(i_freq)+1)
        endif
     enddo
     !tLayer=list(np.interp(bins,envNode[i,j,:],airTemp[i,j,:]))
     !tLayer.append(sfcTemp[i,j])
     !iEnum+=1
     emis=emiss(i,j,i_freq)
     ebar=emis
     !nL=kextInt.shape[0]
     lambert=.False.
     if (sfcBin(i,j)>87) sfcBin(i,j)=87
     if (sfcBin(i,j)>envNode(i,j,5)) sfcBin(i,j)=envNode(i,j,5)
     nlyr=sfcBin(i,j)-envNode(i,j,1)+1
     lyrhgt(1)=0
     !print*, nlyr
     do k=sfcBin(i,j),envNode(i,j,1),-1
        !print*,(k+1-envNode(i,j,1))
        kext_rev(sfcBin(i,j)+1-k)=kext(k+1-envNode(i,j,1))
        asym_rev(sfcBin(i,j)+1-k)=asym(k+1-envNode(i,j,1))
        salb_rev(sfcBin(i,j)+1-k)=salb(k+1-envNode(i,j,1))
        lyrhgt(sfcBin(i,j)+1-k+1)=(sfcBin(i,j)+1-k+1)*dr
        lyrtemp(sfcBin(i,j)+1-k+1)=temp_int(k+1-envNode(i,j,1))
     enddo
     lyrtemp(1)=sfcTemp(i,j)
     fisot=2.7
     !print*, lyrtemp(1:nlyr+1)
     !print*, lyrhgt(1:nlyr+1)
     !print*, sfcTemp(i,j)
     
     call radtran(umu, nlyr, tb, sfcTemp(i,j), lyrtemp, lyrhgt, kext_rev, &
          salb_rev, asym_rev, fisot, emis, ebar, lambert)
     if(isnan(tb)) then
        print*, lyrtemp(1:nlyr+1)
        print*, lyrhgt(1:nlyr+1)
        print*, sfcTemp(i,j)
        print*, nlyr
        print*, kext_rev
        print*, kext_int
        print*, kextL
        call interpol_kext(envNode(i,j,:),kextL,kext_int,sfcBin(i,j))
        print*, kext_int
        print*, envNode(i,j,:), sfcBin(i,j)
        stop
     endif
     !tb1=cAlg.radtran(umu,sfcTemp[i,j],tLayer[::-1],\
     !                        np.arange(nL+1)*0.25,kextInt[::-1],\
     !                        salb[::-1],asym[::-1],2.7,emis,ebar,lambert)
     !except:
     tbout(i,j,i_freq)=tb
     !return
  end do
end subroutine rte



subroutine rte1(binNodes1,sfcBin1,dm1,kextH1,asymH1,salbH1,emiss1,qv1,airTemp1, &
  press1,envNode1,sfcTemp1,cldw1,umu,tbout_13)
  implicit none
  real ::  freqs(13)
  integer :: nfreq,n1, i, j, k, pol, i_freq
  integer :: npol(13), iFreq(15)
  logical :: lambert
  integer :: envNode1(10), binNodes1(5), sfcBin1
  real :: dm1(88), pRate1(88), kextH1(88,8), asymH1(88,8), salbH1(88,8), &
       airTemp1(10), qv1(10), press1(10), &
       cldw1(88), emiss1(13),sfcTemp1
  real :: q, t, pa, cld, z_clw, absair, abswv
  real :: kextL(10), emis, ebar, umu, tk, f, kext_int(73)
  real :: kext(73),salb(73),asym(73)
  real,intent(out) :: tbout_13(13)
  integer :: ireturn
  real :: lyrhgt(74), lyrtemp(74), kext_rev(73), asym_rev(73), salb_rev(73), tb, fisot, temp_int(73)
  real :: dr
  integer :: nlyr
  !print*, sfcBin1
  freqs=(/10.6,10.6,18.7,18.7,23.,37.,37.,89.,89.,166.,166.,186.3,190.3/)
  npol=(/1,0,1,0,1,1,0,1,0,1,0,1,1/)
  iFreq=(/0,0,1,1,2,3,3,4,4,5,5,6,6,7,7/)
  nfreq=8
  dr=0.25
  call interpol_kext(envNode1,airTemp1,temp_int,sfcBin1)
  !print*, sfcBin1
  salb=0
  asym=0

  do i_freq=1,13
     pol=npol(i_freq)
     f=freqs(i_freq)
     do k=1,10
        q=qv1(k)
        tk=airTemp1(k)
        pa=press1(k)
        cld=cldw1(k)
        if (q<1e-3) q=1e-3
        call gasabsr98(f,tk,q*1e-3,pa*1e2,absair,abswv,ireturn)
        if(cld>0.001) then
           call gcloud(f,tk,cld,z_clw)
        else
           z_clw=0.
        endif
        kextL(k)=(absair+abswv+z_clw)
     enddo
     !print*, i_freq, kextL
     call interpol_kext(envNode1,kextL,kext_int,sfcBin1)
     kext=kext_int
     do k=1,73
        if(k+envNode1(1)>0.and.k+envNode1(1)<=88) then
           kext(k)=kext_int(k)+kextH1(k+envNode1(1),iFreq(i_freq)+1)
           salb(k)=salbH1(k+envNode1(1),iFreq(i_freq)+1)*(1-kext_int(k)/kext(k))
           if(salb(k)<0) salb(k)=0
           if(salb(k)>1) salb(k)=1
           asym(k)=asymH1(k+envNode1(1),iFreq(i_freq)+1)
        endif
     enddo
     !kext=kext_int
     !kext=kext+kextH1(envNode1(1):envNode1(5),iFreq(i_freq)+1)
     !salb=salbH1(envNode1(1):envNode1(5),iFreq(i_freq)+1)* &
     !     (kext_int/kext)
     !asym=asymH1(envNode1(1):envNode1(5),iFreq(i_freq)+1)
     emis=emiss1(i_freq)
     ebar=emis
 
     lambert=.False.
     if (sfcBin1>87) sfcBin1=87
     if (sfcBin1>envNode1(10)) sfcBin1=envNode1(10)
     nlyr=sfcBin1-envNode1(1)+1
     lyrhgt(1)=0

     do k=sfcBin1,envNode1(1),-1
        kext_rev(sfcBin1+1-k)=kext(k+1-envNode1(1))
        asym_rev(sfcBin1+1-k)=asym(k+1-envNode1(1))
        salb_rev(sfcBin1+1-k)=salb(k+1-envNode1(1))
        lyrhgt(sfcBin1+1-k+1)=(sfcBin1+1-k+1)*dr
        lyrtemp(sfcBin1+1-k+1)=temp_int(k+1-envNode1(1))
     enddo
     lyrtemp(1)=sfcTemp1
     fisot=2.7
     !print*, kext_rev(1:nlyr)
     !print*, sfcBin1, sfcTemp1, emis, ebar
     nlyr=70
     do while(abs(kext_rev(nlyr))<1e-7.and.nlyr>30)
        nlyr=nlyr-1
     end do
        
     call radtran(umu, nlyr, tb, sfcTemp1, lyrtemp, lyrhgt, kext_rev, &
          salb_rev, asym_rev, fisot, emis, ebar, lambert)


     if(tb>350.or.tb<0.or.isnan(tb)) then
        do k=1,-nlyr
           print*, kext_rev(k),salb_rev(k),asym_rev(k),lyrtemp(k),lyrhgt(k)
        enddo
        !print*, emis, ebar, sfcTemp1, tb, nlyr, umu, lyrtemp(nlyr+1), &
        !     lyrhgt(nlyr+1)
        !stop
        print*,tb,freqs(i_freq)
        tb=-999
     endif
     tbout_13(i_freq)=tb
     !return
  end do
end subroutine rte1


subroutine get_fint(sfcBin1,binNodes1,fint)
  implicit none
  real :: yn(4)
  real, intent(out) :: fint(88)
  integer :: xn(4)
  integer :: binNodes1(5), sfcBin1
  
  xn=(/0,binNodes1(2),binNodes1(4),binNodes1(5)/)
  yn=(/0,0,1,1/)
  call interpol(xn,yn,fint,sfcBin1)

end subroutine get_fint
