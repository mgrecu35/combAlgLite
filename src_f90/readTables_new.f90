subroutine read_tables()
    integer :: nmfreq, nmu

    nmfreq=8
    nmu=5
  !call init_random_seed(rseed1,rseed2)
  !print *, 'Random Seeds : ',rseed1,rseed2
  !print*, ialg

    call readtablesLiang2(nmu,nmfreq)  

end subroutine read_tables

!program test
!
!    call read_tables()
!
!end program