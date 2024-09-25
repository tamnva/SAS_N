!> \file mo_sas_global_variables.f90

!> \brief Global variables used for water quality modelling based on TTDs.

!> \details
!> \note: time variables (e.g, max_age, deni) have a unit of timetep, e.g., 
!>        if the master equation is call at daily, then the unit is day
!>        ..................................monthly (30 days), then the unit is month (x 30 days)

!> \authors Tam Nguyen
!> \date June 2020


module mo_sas_global_variables

  
  use mo_kind,                                   only : i4, dp

  implicit none

  type sas_para
     real(dp), dimension(:), allocatable              :: stor_age       ! Residence time distribution
     real(dp), dimension(:), allocatable              :: conc_age       ! Concentration distribution
     real(dp)                                         :: ka             ! k (or a) parameter of the powerlaw (or beta)
     real(dp)                                         :: b              ! beta parameter (b)
     real(dp)                                         :: half_life      ! half_life of nitrate [timestep] !change to ln(2)/half_life
     real(dp), dimension(:), allocatable              :: concage       ! Concentration distribution
  end type sas_para

  type(sas_para), public                              :: sas

  type nbalance
    real(dp), dimension(:), allocatable              :: n_surplus
    real(dp), dimension(:), allocatable              :: soil_N_store
    real(dp), dimension(:), allocatable              :: soil_N_remove
    real(dp), dimension(:), allocatable              :: soil_N_deni
    real(dp), dimension(:), allocatable              :: soil_N_leach
    real(dp), dimension(:), allocatable              :: sub_N_export
    real(dp), dimension(:), allocatable              :: sub_N_store
    real(dp), dimension(:), allocatable              :: sub_N_deni   
  end type nbalance

  type(nbalance), public                             :: n_balance

end module mo_sas_global_variables





