! define the data module.
module data_params
    use iso_fortran_env
    implicit none
    
    !
    !   system parameters.
    !---------------------------------------------------------------------------------
    integer(kind=int32),parameter :: lx = 16, ly = 8
    real(kind=real64),parameter :: Jij = 1.0
    real(kind=real64) :: hVal
    integer(kind=int32), dimension(lx, ly) :: spin
    
    !
    !   simulation parameters.
    !---------------------------------------------------------------------------------
    real(kind=real32),parameter :: Tinit = 1.0, Tfinal = 4.0, dT = 0.1
    integer(kind=int32),parameter :: nMCSweaps = 20000, runEquil = 10000
    integer(kind=int32) :: i,j,k,metro_flag, nAccept, nReject, TloopSteps, i1,j1
    real(kind=real64) :: E, Et, Eold, Tcurr, bta
    
    
    !
    !   output parameters.
    !----------------------------------------------------------------------------------
    integer(kind=int32),parameter :: freq_thermo = 10, freq_trj = 10
                                   
    
    !   files - units
    !----------------------------------------------------------------------------------
    ! input and output file units range.
    !
    ! input files  ==> 10 -- 20.
    ! output files ==> 20 -- 30.
    !
    integer(kind=int32),parameter :: fout_log = 20, &
                                   & fout_m_prof = 21, &
                                   & fout_config = 22, &
                                   & fout_thermo = 23
    
    !   files - names
    !----------------------------------------------------------------------------------
    ! define the file names.
    character(len=*),parameter :: fname_log = "out_log.txt", &
                                & fname_m_prof = "out_m_prof.dat", &
                                & fname_config = "out_config.dat", &
                                & fname_thermo = "out_thermo.dat"

        
end module data_params

!===========================================================================!

