module zku_dm_mod
  real:: zku_dm(161,2)
end module zku_dm_mod
subroutine read_zku_dm()
  use zku_dm_mod
  open(10,file="zku_dm.txt")
  do i=1,161
     read(10,*) zku_dm(i,1),zku_dm(i,2)
  enddo
end subroutine read_zku_dm
subroutine iwc_from_dfr(n,z_obs,temp,iwc,dm)
  use tablep2
  implicit none
  integer :: n,k,ibin
  real :: z_obs(n,2), temp(n)
  real, intent(out) :: iwc(n), dm(n)
  real :: dwr1,dn
  do k=1,n
     dwr1=z_obs(k,1)-z_obs(k,2)
     if(dwr1<0.5) dwr1=0.5
     if(dwr1>10.) dwr1=10.
     call bisection2(dwr(1:253),253,dwr1,ibin)
     dn=(z_obs(k,1)-zkus(ibin))/10.
     iwc(k)=swc(ibin)*10**dn
     dm(k)=dms(ibin)
  enddo

end subroutine iwc_from_dfr

subroutine iwc_from_dfr2(n,z_obs,zsim_out,temp,iwc,dm,prate_out)
  use tablep2
  use zku_dm_mod
  implicit none
  integer :: n,k,ibin
  real :: z_obs(n,2), temp(n)
  real, intent(out) :: iwc(n), dm(n), zsim_out(n,3), prate_out(n)
  real :: dwr1,dn
  integer :: ibin2
  real :: zka1,zka2,dm1,dm2, dzka_dm, dn2
  zsim_out=0
  do k=1,n
     if(z_obs(k,1)>40) then
        dwr1=z_obs(k,1)-z_obs(k,2)
        if(dwr1<0.45) dwr1=0.45
        if(dwr1>10.) dwr1=10.
        call bisection2(dwr(1:253),253,dwr1,ibin)
        dn=(z_obs(k,1)-zkus(ibin))/10.
        iwc(k)=swc(ibin)*10**dn
        dm(k)=dms(ibin)
        prate_out(k)=snowrate(ibin)*10**dn
     else
        ibin2=int((z_obs(k,1)-1)/0.25)+1
        ibin2=int((z_obs(k,1)+12)/0.25)+1
        if(ibin2<1) ibin2=1
        !if(ibin2>161) ibin2=161
        if(ibin2>253) ibin2=253
        !call bisection2(dms(1:253),253,zku_dm(ibin2,2),ibin)
        dn=(z_obs(k,1)-zkus(ibin))/10.
        iwc(k)=swc(ibin)*10**dn
        dm(k)=dms(ibin)
        zka1=zkas(ibin)+10*dn
        dm1=dm(k)
        dm2=dm1+0.1
        call bisection2(dms(1:253),253,dm2,ibin2)
        prate_out(k)=snowrate(ibin)*10**dn
        dn2=(z_obs(k,1)-zkus(ibin2))/10.
        zka2=zkas(ibin2)+10*dn2
        if(z_obs(k,2)>12 .and. z_obs(k,2)<35) then
           dzka_dm=(zka2-zka1)/0.1
           dm1=dm1+(z_obs(k,2)-zka1)*dzka_dm/(dzka_dm**2+16)
           call bisection2(dms(1:253),253,dm1,ibin)
           dn=(z_obs(k,1)-zkus(ibin))/10.
           !iwc(k)=swc(ibin)*10**dn
           !dm(k)=dms(ibin)
           zka2=zkas(ibin)+10*dn
           !prate_out(k)=snowrate(ibin)*10**dn
           zsim_out(k,:)=(/z_obs(k,2),zka2,zka1/)! z_obs(k,1)
        endif
        !print*, zku_dm(ibin2,2),dm(k)
     endif
  enddo

end subroutine iwc_from_dfr2
