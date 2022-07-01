!>       \file mo_sas_master_equation.f90

!>       \brief solve master equation

!>       \details solve master equation to get RTD, TTD, concentration in the outflow

!>       \authors Tam Nguyen

!>       \date June 2020

module mo_sas_master_equation

  use mo_kind,                           only : i4, dp
  use mo_sas_global_variables,           only : sas_para
  use mo_sas_utils,                      only : cumsum,        &
                                                cdfbeta,       &
                                                eval_sas

  implicit none

contains

  ! ------------------------------------------------------------------

  !    NAME
  !        master_equation

  !    PURPOSE
  !>       \brief caculates solute export using the StorAge Selection (SAS) function

  !>       \details caculates solute export using the StorAge Selection (SAS) function

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: L1_inflow"  inflow to the subsurface [mm]
  !>       \param[in] "real(dp), dimension(:) :: L1_outflow" outflow from the subsurface [mm]
  !>       \param[in] "real(dp), dimension(:) :: L1_inc"   solute concentration in the inflow [mg/L]
  !>       \param[in] "real(dp), dimension(:) :: max_age"   allowable maximum age in storage [timestep]
  !>       \param[in] "real(dp), dimension(:) :: max_old_fraction"   allowable maximum old water fraction in storage [-]
  !>       \param[in] "real(dp), dimension(:) :: sas_function"   sas function (1 = powerlaw, 2 = beta)

  !    INTENT(INOUT), OPTIONAL
  !>       \param[in] "real(dp), dimension(:) :: L1_outc"  solute concentration in the outflow [mg/L]
  !>       \param[in] "real(dp), dimension(:) :: L1_sas"  transit time distribution, residence time distribution [-]
  !>       \param[in] "real(dp), dimension(:) :: denitri_amount"  amount of nitrate removed by denitrification [kg N /ha]
  !    HISTORY
  !>       \authors Tam Nguyen

  !>       \date June 2020

  subroutine master_equation(L1_sas,             & ! sas_para variables inflow,   
                             inflow,             & ! inflow to the subsurface (mm)       
                             in_conc,            & ! concentration in the inflow (mg/L) 
                             outflow,            & ! outflow from the subsurface (mm)
                             out_conc,           & ! concentration in the out (mg/L)   
                             max_age,            & ! maximum allowable age in storage
                             max_old_fraction,   & ! maximum allowable olest water fraction in storage
                             sas_function,       & ! StorAge selection function (1=powerlaw, 2=beta)
                             median_tt,          & ! median TT
                             median_rt,          & ! median residence time
                             mean_rt,            & ! mean residence time  
                             mean_tt,            & ! median transit time
                             denitri_amount,     & !
                             sub_n_stor,         & ! subsurface N storage (kg/ha)
                             age_rank_discharge)   ! age of water when loading is 50 percent


    implicit none

    type(sas_para), intent(inout)  :: L1_sas
    real(dp),       intent(in)     :: inflow
    real(dp),       intent(in)     :: in_conc
    real(dp),       intent(in)     :: outflow
    real(dp),       intent(out)    :: out_conc
    integer(i4),    intent(in)     :: max_age
    real(dp),       intent(in)     :: max_old_fraction
    integer(i4),    intent(in)     :: sas_function
    integer(i4),    intent(out)    :: median_tt          ! median transit time
    integer(i4),    intent(out)    :: median_rt          ! median residence time
    real(dp),       intent(out)    :: mean_rt            ! median residence time
    real(dp),       intent(out)    :: mean_tt            ! median residence time
    real(dp),       intent(out)    :: denitri_amount     ! denitrification (kg N ha-1)
    real(dp),       intent(out)    :: sub_n_stor         ! subsurface N storage
    real(dp), dimension(:), allocatable, intent(out) :: age_rank_discharge

    !local variables
    integer(i4)                           :: i, j                                              !counter for loop
    integer(i4)                           :: tmax                                              !maximum age in storage or discharge
    real(dp), dimension(:), allocatable   :: sas, deriv_sas                                    !sas (derivative) function evaluate over the range of normalized age-ranked storage 
    real(dp)                              :: residual_discharge, discharge_age_i, store_age_i  !subsurface denitrification rate
    real(dp)                              :: ka, b, temp, eps_dp                                       !parameter of the powerlaw or beta function
    real(dp)                              :: oldest_water_fraction, age_1                      !volume of water with tmax               
    real(dp), dimension(:), allocatable   :: norm_age_rank_stor, stor_age, conc_age

    eps_dp = epsilon(1.0_dp)

    !get paramter of the powerlaw or beta function
    if (sas_function == 1) then
      ka = L1_sas%ka
      b  = 0.0_dp
    else if (sas_function == 2) then
      ka = L1_sas%ka
      b  = L1_sas%b
    else
      stop
    end if

    !check maximum age in storage
    tmax = size(L1_sas%stor_age)

    !Denitrification amount (regarless of how much N is exported in discharge)
    denitri_amount = 0.01_dp *                                                                      &
                     sum((/inflow, L1_sas%stor_age(:) - (/0.0_dp, L1_sas%stor_age(1:(tmax - 1))/)/)*  &
                         (/in_conc, L1_sas%conc_age/) * (1.0_dp - exp(-L1_sas%half_life)))

    !youngest water in storage
    allocate(stor_age(1))
    stor_age(:) = eval_sas((/L1_sas%stor_age(1)/L1_sas%stor_age(tmax)/), sas_function, ka, b)
    age_1 = max(0.0_dp, inflow - outflow * stor_age(1))

    !update normalized age-ranked storage [0,1]
    allocate(norm_age_rank_stor(tmax + 1))
    norm_age_rank_stor(:) = (/0 + age_1, L1_sas%stor_age(:) + age_1/)/ (L1_sas%stor_age(tmax) + age_1)

    !evaluate the sas function over the normalized age-ranked storage
    allocate(sas(tmax + 1))
    sas(:) = eval_sas(norm_age_rank_stor, sas_function, ka, b)

    !*************************************************************Master equation
    !solve the master equation using Mehthod of Lines and Forward Euler
    deallocate(stor_age)
    allocate(stor_age(tmax + 1))

    !stor_age(:) = (/0.0_dp, L1_sas%stor_age(:)/) + inflow - outflow * sas(:)
    !this solution of the master equation does not ensure that stor_age(:)is a monotonically increasing function

    !update age-ranked storage with inflow (the amount of storage with age < Ti)
    stor_age(:) = (/0.0_dp, L1_sas%stor_age(:)/) + inflow

    !calculate age_rank discharge (the amount of discharge with age < Ti)
    allocate(age_rank_discharge(tmax + 1))
    age_rank_discharge(:) = outflow * sas(:)

    !derivative of stor_age and age_rank_discharge (the volume of storage or discharge with age Ti)
    stor_age(:) = stor_age(:) - (/0.0_dp, stor_age(1:tmax)/)
    age_rank_discharge(:) =  age_rank_discharge(:) - (/0.0_dp, age_rank_discharge(1:tmax)/)

    !initialize residual discharge
    residual_discharge = 0.0_dp

    do i = 1, tmax + 1

      !takes water of this age if residual_discharge > 0.0
      age_rank_discharge(i) = age_rank_discharge(i) + residual_discharge

      !update residual discharge
      residual_discharge = 0.0_dp

      if (age_rank_discharge(i) .le. stor_age(i)) then

        !update storage
        stor_age(i) = stor_age(i) - age_rank_discharge(i)

      else

        !remaining discharge that needs to be taken from older ages
        residual_discharge = age_rank_discharge(i) - stor_age(i)

        !update age rank discharge
        age_rank_discharge(i) = stor_age(i)

        !update stor_age
        stor_age(i) = 0.0_dp

      end if

    end do

    !If with oldest discharge there is no water for outlow, take from all ages (volume weighted)
    age_rank_discharge(:) = age_rank_discharge(:) + residual_discharge * stor_age(:)/sum(stor_age) 
    stor_age(:) =  stor_age(:) - residual_discharge * stor_age(:)/sum(stor_age)

    !Subsurface N storage
    sub_n_stor = 0.01_dp * sum(exp(-L1_sas%half_life) * (/in_conc, L1_sas%conc_age(:)/) * stor_age(:))

    !convert back to cummulative sum of stor_age and age_rank_discharge
    stor_age(:) = cumsum(stor_age)
    age_rank_discharge(:) =  cumsum(age_rank_discharge)

    !update sas
    sas(:) = age_rank_discharge(:)/age_rank_discharge(tmax + 1)
    !*************************************************************end Master Equation

    !update solute concentration in each parcel
    allocate(conc_age(tmax + 1))
    conc_age(:) = (/in_conc, L1_sas%conc_age(:)/) * exp(-L1_sas%half_life)

    deallocate(L1_sas%concage) 
    allocate(L1_sas%concage(tmax + 1)) 
    L1_sas%concage(:) = (/in_conc, L1_sas%concage(1:tmax)/)

    !calcuate pQ * dT
    allocate(deriv_sas(tmax + 1))
    deriv_sas(:) = sas(:) - (/ 0.0_dp, sas(1:tmax)/)

    !solute concentration in the outflow
    out_conc = sum(conc_age(:) * deriv_sas(:))

    !initialized output (median transit time and residence time)
    median_tt = 1_i4
    median_rt = 1_i4
    mean_tt = 0.0_dp

     
    !calcualte median TT, RT50
    do i = 1, size(sas)

      if (sas(i) < 0.5) median_tt = i

      if (i == 1) then
        mean_tt = mean_tt + sas(i)
      else
        mean_tt = mean_tt + (sas(i) - sas(i-1)) * i
      end if

      if (stor_age(i)/stor_age(tmax + 1) < 0.5) median_rt = i
    end do

    ! calculate mean RT
    mean_rt = stor_age(1)

    do i = 2, tmax + 1
      temp = (stor_age(i) - stor_age(i-1))
      if (temp .gt. eps_dp) mean_rt = mean_rt + (stor_age(i) - stor_age(i-1)) * i
    end do

    mean_rt = mean_rt/stor_age(tmax + 1)

    !update TTD and CTD
    deallocate(L1_sas%stor_age)
    deallocate(L1_sas%conc_age)
  
    allocate(L1_sas%stor_age(tmax + 1))
    allocate(L1_sas%conc_age(tmax + 1))   
    L1_sas%stor_age(:) = stor_age(:)
    L1_sas%conc_age(:) = conc_age(:)

  end subroutine master_equation

end module mo_sas_master_equation


























