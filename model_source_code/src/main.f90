
program main

  use mo_sas_global_variables,        only: sas, n_balance                                  
  use mo_sas_master_equation,         only: master_equation
  use mo_kind,                        only: i4, dp
  use mo_sas_utils,                   only: n_soil_balance, rmse, read_input, read_param

  implicit none

  !Local variables
  integer(i4)                             :: max_age, sas_function, nTimeStep, median_tt, median_rt, i, j
  integer(i4), dimension(:), allocatable  :: year
  real(dp)                                :: mean_tt, mean_rt, denitri_amount, Nleach, rootMSE, meanTTD
  real(dp)                                :: max_old_fraction, S0, C0, N0, Qmax, cQin
  real(dp), dimension(:), allocatable     :: Qin, Nsurplus, Cobs, cQout, alpha, beta, a, b, half_life, f_instream, cQoutlet
  real(dp), dimension(:), allocatable     :: subNstore, age_rank_discharge, Nwwtp
  character(50)                           :: fName

  !Variable for checking mass balance (optional, delete in final version)
  real(dp) :: temp1, temp2, temp3

  !Read input time series and model setting
  fName = './input/data.txt'
  call read_input(fName, max_age, max_old_fraction, sas_function, S0, C0, N0, Qmax, year, Qin, Nsurplus, Cobs, Nwwtp)

  !Read parameter file and initial condition
  fName = './input/parameter.txt'
  call read_param(fName, alpha, beta, a, b, half_life, f_instream)

  !Allocate output variable
  nTimeStep = size(Qin)
  allocate(cQout(nTimeStep))
  allocate(cQoutlet(nTimeStep))
  allocate(subNstore(nTimeStep))

  allocate(n_balance%n_surplus(nTimeStep))
  allocate(n_balance%soil_N_store(nTimeStep))
  allocate(n_balance%soil_N_remove(nTimeStep))
  allocate(n_balance%soil_N_deni(nTimeStep))
  allocate(n_balance%soil_N_leach(nTimeStep))
  allocate(n_balance%sub_N_export(nTimeStep))
  allocate(n_balance%sub_N_store(nTimeStep))
  allocate(n_balance%sub_N_deni(nTimeStep))

  !open output files
  open(10, file = "./output/RMSE.txt")
  open(11, file = "./output/SimC.txt")
  open(12, file = "./output/SimC_outSub.txt")
  open(13, file = "./output/Nbalance.txt")
  open(14, file = "./output/TTD.txt")


  !iter over number of parameter sets in the './input/parameter.txt' file
  do i = 1, size(alpha)

    !print to screen every 200 iteration
    if (mod(i, 200) == 0) print*, "Iteration: ", i

    !allocate variable at the begining of each iteration
    allocate(sas%stor_age(1))
    allocate(sas%conc_age(1))
    allocate(sas%concage(1)) 

    !assign parameter to sas variable
    sas%ka = a(i)
    sas%b = b(i)
    sas%half_life = half_life(i)

    !initial condition
    sas%stor_age(1) = S0
    sas%conc_age(1) = C0
  
    !mean TTD
    meanTTD = 0.0_dp

    !Iterate over number of time step

    do j = 1, nTimeStep

      !calculate n_soil_balance
      call n_soil_balance(N0, Nsurplus(j), Nleach, alpha(i) * Qin(j)/Qmax)

      cQin = beta(i) * 100.0_dp * Nleach/Qin(j)

      !call master equation: Qin, cQin -> cQout
      call master_equation(sas,                & ! sas_para variables Qin,   
                           Qin(j),             & ! Qin to the subsurface (mm)       
                           cQin,               & ! concentration in the Qin (mg/L) 
                           Qin(j),             & ! Qout from the subsurface (mm) = Qin flow at yearly time step
                           cQout(j),           & ! concentration in the out (mg/L)   
                           max_age,            & ! maximum allowable age in storage
                           max_old_fraction,   & ! maximum allowable olest water fraction in storage
                           sas_function,       & ! StorAge selection function (1=powerlaw, 2=beta)
                           median_tt,          & ! median TT
                           median_rt,          & ! median residence time
                           mean_rt,            & ! mean residence time  
                           mean_tt,            & ! median transit time
                           denitri_amount,     & !
                           subNstore(j),       & !
                           age_rank_discharge)

      n_balance%n_surplus(j) = Nsurplus(j)                        !kg/ha
      n_balance%soil_N_store(j) = N0                              !kg/ha
      n_balance%soil_N_remove(j) = Nleach                         !kg/ha
      n_balance%soil_N_deni(j) = Nleach * (1.0_dp - beta(i))      !kg/ha
      n_balance%soil_N_leach(j) = Nleach * beta(i)                !kg/ha
      n_balance%sub_N_export(j) = cQout(j)* Qin(j)/100.00_dp      !kg/ha
      n_balance%sub_N_store(j) = subNstore(j)                     !kg/ha
      n_balance%sub_N_deni(j) = denitri_amount                    !kg/ha


      !instream nitrate removal - C at outlet (assuming that Qwwtp << Qin)
      cQoutlet(j) = (1 - f_instream(i)) * (cQout(j)  + Nwwtp(j)* 100.0_dp/Qin(j))

      !mean TTD

       if (j > 150) then
         meanTTD = meanTTD + mean_tt * Qin(j)
       end if

    end do

    deallocate(sas%stor_age)
    deallocate(sas%conc_age)
    deallocate(sas%concage)
   
    !Caculate RMSE
    rootMSE  =  rmse(Cobs, cQoutlet)
    write(10, '(f15.5)') rootMSE              !Model performance rmse
    write(11, '(999f15.3)') cQoutlet          !Simulated instream N-NO3
    write(12, '(999f15.3)') cQout             !Simulated instream N-NO3

    write(13, '(999f15.3)') n_balance%n_surplus(151:215)
    write(13, '(999f15.3)') n_balance%soil_N_store(151:215)
    write(13, '(999f15.3)') n_balance%soil_N_remove(151:215)
    write(13, '(999f15.3)') n_balance%soil_N_deni(151:215)
    write(13, '(999f15.3)') n_balance%soil_N_leach(151:215)
    write(13, '(999f15.3)') n_balance%sub_N_export(151:215)
    write(13, '(999f15.3)') n_balance%sub_N_store(151:215)
    write(13, '(999f15.3)') n_balance%sub_N_deni(151:215)

    write(14, '(999f15.3)') meanTTD/(sum(Qin(151:215)))   !mean TT from 1950-2014 (65 years)

  end do

  close(10)
  close(11)
  close(12)
  close(13)
  close(14)



end program main















