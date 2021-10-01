
subroutine initialize()
    use data_params
    use iso_fortran_env
    implicit none
    
    real(kind=real32) :: r
    
    do i=1,lx
      do j=1,ly
        call random_number(r)
        if (r <= 0.5) then
          spin(i,j) = 1
        else
          spin(i,j) = -1    
        end if
      end do
    end do
    
    return
end subroutine initialize

!===========================================================================!

subroutine print_config()
    use data_params
    use iso_fortran_env
    implicit none
    
    print*, new_line('A')
!    print*,'*********************************************************************'
!    print*, new_line('A')
!    print*,'=============================================================='
 !   print*, new_line('A')
    
    print*," -- system configuration -- "
    print*,"-------------------------->>  lx"
    do j=1,ly
      do i=1,lx
        write(*,'(i3)',advance='no') spin(i,j)
!         print*, spin(i,j), ", "
      end do
      print*, new_line('A')
!       print *, spin(i, 1:ly)
    end do
    
!    print*,'*********************************************************************'
!    print*, new_line('A')
    print*,'=============================================================='
!    print*, new_line('A')
    
    return
end subroutine print_config

!===========================================================================!
subroutine put_all_ones()
    use data_params
    use iso_fortran_env
    implicit none
    
    spin = 1
    
    return
end subroutine put_all_ones
!===========================================================================!
subroutine put_all_zeros()
    use data_params
    use iso_fortran_env
    implicit none
    
    spin = 0
    
    return
end subroutine put_all_zeros

!===========================================================================!
subroutine put_header(fout)
    use data_params
    use iso_fortran_env
    implicit none
    
    integer(kind=int32) :: fout
    
    write(fout, *)"**********************************************************************"
    write(fout, *)"*                                                                    *"
    write(fout, *)"*  This code is to simulate a 2-D Ising model, with and without the  *"
    write(fout, *)"*  external magnetic field. Also the exchange interaction parameter  *"
    write(fout, *)"*  could be treated randomly. Physical meaning of this is yet to be  *"
    write(fout, *)"*  understood.                                                       *"
    write(fout, *)"*  This code is written by: Dr. S. Siva Nasarayya Chari,             *"
    write(fout, *)"*  Last modified on: Nov 20, 2020.                                   *"
    write(fout, *)"*                                                                    *"
    write(fout, *)"**********************************************************************"    
    
    return 
end subroutine put_header
!===========================================================================!
subroutine flip_a_random_spin()
    use data_params
    use iso_fortran_env
    implicit none
    
    real(kind=real32) :: r1, r2
    
    call random_number(r1)
    call random_number(r2)
    
    i = mod(int(r1 * lx), lx)+1
    j = mod(int(r2 * ly), ly)+1
    
!    print '(/,a,i0,a,i0,a )',"flipping the spin at the site: (",i,", ",j,")."
    
    spin(i,j) = -spin(i,j)
    
    return
end subroutine flip_a_random_spin
!===========================================================================!
subroutine save_magn_profile()
    use data_params
    use iso_fortran_env
    implicit none
    
!    integer(kind=int32),intent(in) :: fout
    
    do i=1,lx
      write(fout_m_prof, '(i0, 2x, i0)') i, sum(spin(i,:))
    end do
    write(fout_m_prof, *)new_line('A')
    
    return
end subroutine save_magn_profile
!===========================================================================!
subroutine open_all_files()
    use data_params
    use iso_fortran_env
    implicit none
        
    open(unit=fout_m_prof, file=fname_m_prof, status='unknown')
    open(unit=fout_log, file=fname_log, status='unknown')
    open(unit=fout_config, file=fname_config, status='unknown')
    open(unit=fout_thermo, file=fname_thermo, status='unknown')
    
    return
end subroutine open_all_files
!===========================================================================!
subroutine close_all_files()
    use data_params
    use iso_fortran_env
    implicit none
    
    close(fout_m_prof)
    close(fout_log)
    close(fout_config)
    close(fout_thermo)
    
    return
end subroutine close_all_files
!===========================================================================!
subroutine check_plots()
    use data_params
    use iso_fortran_env
    implicit none
    
!    character(*),parameter:: cmd = "xmgrace -free -block "//trim(fname_m_prof)//" -bxy 1:2 &"
    character(100) :: cmd1,cmd2,cmd3
    
    cmd1 = "xmgrace -free -block "//trim(fname_m_prof)//" -bxy 1:2 &"
    cmd2 = "xmgrace -free -block "//trim(fname_thermo)//" -bxy 1:2 -bxy 1:3 &"
    cmd3 = "xmgrace -free -block "//trim(fname_thermo)//" -bxy 1:4 -bxy 1:5 &"

    !call execute_command_line(cmd1)
    call execute_command_line(cmd2)
    call execute_command_line(cmd3)
    
    return
end subroutine check_plots
!===========================================================================!
subroutine save_config_trj()
    use data_params
    use iso_fortran_env
    implicit none
    
    do j=1,ly
      do i=1,lx
        write(fout_config,'(i3)',advance='no') spin(i,j)
      end do
      write(fout_config, *)new_line('A')
    end do
    
    return    
end subroutine save_config_trj
!===========================================================================!
subroutine get_config_energy()
    use data_params
    use iso_fortran_env
    implicit none
    
    integer(kind=int32) :: p,q,l,m
    real(kind=real64) :: esum, nnsum, r1
    
    esum = 0.0d0
    do i=1,lx
      do j=1,ly
        ! get the nearest neighbors (all n.n. pairs for a given spin at site i,j).
        call random_number(r1)
        hVal = 0.005    ! r1*0.01
        
        p = mod(i, lx)+1
        q = mod(i+lx-2, lx)+1
        l = mod(j, ly)+1
        m = mod(j+ly-2, ly)+1
        
        !                     (i,l)
        !                       |
        !                       |
        !        (q,j) <----- (i,j) -----> (p,j)
        !                       |
        !                       |
        !                     (i,m)
        
        nnsum = spin(p,j) + spin(q,j) + spin(i,l) + spin(i,m)
        esum = esum - Jij*spin(i,j)*nnsum*0.5 -hVal*spin(i,j)
      end do
    end do
    
    E = esum
        
    return
end subroutine get_config_energy
!===========================================================================!
subroutine write_log_file()
    use data_params
    use iso_fortran_env
    implicit none
    
    call put_header(fout_log)
    write(fout_log, *)new_line('A')
    
    write(fout_log, '(a, i0, a, i0)') "2D - square lattice, of size lx = ", lx, ", ly = ", ly
    write(fout_log, '(a, f0.2)') "Exchange interaction parameter, J = ", Jij
    write(fout_log, '(a, f0.4)') "Applied external magnetic field value, h = ", hVal
    write(fout_log, '(a, f0.4, a, f0.4, a, f0.4)') "Initial Temperature: ", Tinit, ", final temperature: ", Tfinal, ", dT: ", dT
    write(fout_log, '(a, i0, a)') "At this rate it requires ", TloopSteps, & 
    & " no. of steps to reach the final value of Temperature."
    write(fout_log, '(a, i0)') "Number of MC sweaps = ", nMCSweaps
    
end subroutine write_log_file
!===========================================================================!
subroutine accept_metropolis()
    use data_params
    use iso_fortran_env
    implicit none
    
    real(kind=real64) :: r1, dE
    
    !metro_flag = 0; nAccept = 0; nReject = 0;
    
    !print*,'current value of bta = ', bta
    
    dE = Et - Eold
    if (dE <= 0) then
        spin(i1, j1) = spin(i1, j1)   ! accept the spin-flip
        E = Et
    else
        call random_number(r1)  ! else, accept as per the Boltzmann prob.
        if (r1 < exp(-bta * dE)) then
            spin(i1, j1) = spin(i1, j1)
            E = Et
     !       metro_flag = 1      ! this flag will be 1 if accepted, else 0.
     !       nAccept = nAccept + 1
        else
            spin(i1, j1) = -spin(i1, j1)
            E = Eold
      !      metro_flag = 0
      !      nReject = nReject + 1
        end if
    end if
    
end subroutine accept_metropolis
!===========================================================================!
