program ising_2d
    use data_params
    use iso_fortran_env
    implicit none
    
    integer(kind=int32) :: itr, itr_temp, avgOver, iflip
    real(kind=real64) :: magn, e_sum, e2_sum, m_sum, m2_sum, sp_heat, chi, r1, r2, eavg1, eavg2, mavg1, mavg2, m1, dE, &
                        & varE, varM
!    integer(kind=int32), parameter :: freq_thermo = 100
    
    call initialize()
    call open_all_files()
    
    TloopSteps = int( (Tfinal - Tinit)/dT )+1
    !Tcurr = Tinit
    !bta = 1.0/Tcurr
    
    call write_log_file()
    
    do itr_temp = 1,TloopSteps  ! temperature loop.
    
    Tcurr = Tinit + (itr_temp-1)*1.0 * dT
    bta = 1.0/Tcurr
    
    m_sum = 0.0     ! to calculate < magn >
    e_sum = 0.0     ! to calculate < energy >
    e2_sum = 0.0    ! to obtain < energy^2 >
    m2_sum = 0.0    ! to obtain < magn^2 >
    
    do itr = 1, nMCSweaps
        
        do iflip = 1, lx*ly
        
            call get_config_energy()
            Eold = E
            
            ! site to flip the spin ==> i1, j1
            call random_number(r1)
            call random_number(r2)
            
            i1 = int(r1*lx)+1
            j1 = int(r2*ly)+1
            ! flip it.
            spin(i1, j1) = -spin(i1, j1)
            call get_config_energy()
            Et = E
            call accept_Metropolis()
            
        end do  ! iflip-loop end
        
        if (itr > runEquil) then
            e_sum = e_sum + E               ! + eavg1/real(lx*ly)              !+ E
            e2_sum = e2_sum + E*E           ! + eavg2/real(lx*ly)            !+ E*E
            magn = sum(spin)
            m_sum = m_sum + abs(magn)            ! + mavg1/real(lx*ly)
            m2_sum = m2_sum + magn * magn   ! + mavg2/real(lx*ly)
        end if
        
    end do  ! end-mc-sweaps loop.
    
    avgOver = nMCSweaps - runEquil
    e_sum = e_sum/real(avgOver)
    e_sum = e_sum/real(lx*ly)
    m_sum = m_sum/real(avgOver)
    m_sum = m_sum/real(lx*ly)
    e2_sum = e2_sum/real(avgOver)
    e2_sum = e2_sum/real(lx*ly * lx*ly)
    m2_sum = m2_sum/real(avgOver)
    m2_sum = m2_sum/real(lx*ly * lx*ly)
    varE = e2_sum - e_sum**2
    varM = m2_sum - m_sum**2
    sp_heat = varE * bta**2
    chi = varM * bta
    
    write(fout_thermo,*)Tcurr, e_sum, m_sum, sp_heat, chi
    
    end do  ! end-temperature loop
    
    call save_magn_profile()
    call close_all_files()
    call check_plots()
    
end program ising_2d
