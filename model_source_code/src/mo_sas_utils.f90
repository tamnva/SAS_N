!>       \file mo_sas_utils.f90

!>       \brief some ultilities functions for solving the storage selection function

!>       \details some ultilities functions for solving the storage selection function

!>       \authors Tam Nguyen

!>       \date June 2020

MODULE mo_sas_utils

  use mo_kind,                                only : i4, dp

  implicit none

  public :: cdfbeta
  public :: read_param
  public :: read_input
  public :: rmse
  public :: n_soil_balance

contains

  function cdfbeta (x, a, b)
  ! ------------------------------------------------------------------

  !    NAME
  !        cdfbeta

  !    PURPOSE
  !>       \brief calculate value of the cummulative beta distribution function at x

  !>       \details calculate value of the cummulative beta distribution function at x
  !>        cdfbeta(x,a,b)   = beta(x,a,b)/beta(a,b) = incomplete beta function / complete beta function

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: x" input array
  !>       \param[in] "real(dp) :: a" parameter of the beta function
  !>       \param[in] "real(dp) :: b" parameter of the beta function

  !    RETURN
  !>       \return real(dp) :: value of the beta function at x

  !    HISTORY
  !>       \authors Majumder, K.L., Bhattacharjee, G.P. 

  !    REFERENCE
  !>       Majumder, K.L., Bhattacharjee, G.P. (1973).Algorithm AS 63: The incomplete Beta Integral, 
  !>       Applied Statistics, 22(3), 409-411

  !    MODIFIED
  !>       Tam Nguyen (June, 2020)

  implicit none

  real(dp), parameter :: acu = 0.1e-14_dp
  real(dp)            :: ai
  real(dp)            :: beta
  real(dp)            :: cdfbeta
  real(dp)            :: cx
  logical             :: indx
  integer(i4)         :: ns
  real(dp)            :: a
  real(dp)            :: pp
  real(dp)            :: psq
  real(dp)            :: b
  real(dp)            :: qq
  real(dp)            :: rx
  real(dp)            :: temp
  real(dp)            :: term
  real(dp)            :: x
  real(dp)            :: xx

  beta = log(gamma ( a ) ) + log(gamma ( b ) ) - log(gamma ( a + b ))  
  cdfbeta = x

  !Special cases.
  if ( x == 0.0_dp .or. x == 1.0_dp ) then
    return
  end if

  !Change tail if necessary and determine S.
  psq = a + b
  cx = 1.0_dp - x

  if ( a < psq * x ) then
    xx = cx
    cx = x
    pp = b
    qq = a
    indx = .true.
  else
    xx = x
    pp = a
    qq = b
    indx = .false.
  end if

  term = 1.0_dp
  ai = 1.0_dp
  cdfbeta = 1.0_dp
  ns = int ( qq + cx * psq )

  !Use Soper's reduction formula.
  rx = xx / cx
  temp = qq - ai
  if ( ns == 0_i4 ) then
    rx = xx
  end if

  do

    term = term * temp * rx / ( pp + ai )
    cdfbeta = cdfbeta + term
    temp = abs ( term )

    if ( temp <= acu .and. temp <= acu * cdfbeta ) then

      cdfbeta = cdfbeta * exp ( pp * log ( xx ) &
        + ( qq - 1.0_dp ) * log ( cx ) - beta ) / pp

      if ( indx ) then
        cdfbeta = 1.0_dp - cdfbeta
      end if

      exit

    end if

    ai = ai + 1.0_dp
    ns = ns - 1_i4

    if ( 0_i4 <= ns ) then
      temp = qq - ai
      if ( ns == 0_i4 ) then
        rx = xx
      end if
    else
      temp = psq
      psq = psq + 1.0_dp
    end if

  end do

  return

  end function cdfbeta
  ! ------------------------------------------------------------------

  !    NAME
  !        cumsum

  !    PURPOSE
  !>       \brief calculate cummulative sum of an array

  !>       \details calculate cummulative sum of an array
  !>       \cumsum[x]   = /x1, x1 + x2,..., x1 + ... + xn/

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: x" input array

  !    RETURN
  !>       \return real(dp) :: array of cummulative summation of x

  !    HISTORY
  !>       \authors Tam Nguyen

  !>       \date June 2020

  function cumsum(x)
    implicit none
    
    real(dp), dimension(:), intent(in)       :: x
    real(dp), dimension(size(x))             :: cumsum
    integer(i4)                              :: i

    !initialize result
    cumsum = 0.0_dp

    !first element of the output array
    cumsum(1) = x(1)

    !calculate cumulative summation
    if(size(x) > 1) then
       do i = 2,size(x)
          cumsum(i)  = cumsum(i-1) + x(i)
       end do
    end if
   
  end function cumsum

  ! ------------------------------------------------------------------

  !    NAME
  !        eval_sas

  !    PURPOSE
  !>       \brief caculate value of the sas function (beta or powerlaw)

  !>       \details calculate the powerlaw or beta distribution with given the funciton parameters
  !>       \eval_sas[x,ka]   = powerlaw(x,ka) = x ** ka
  !>       \eval_sas[x,ka,b] = beta(x,ka,b) = see the cdfbeta function

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: x" input array

  !    RETURN
  !>       \return real(dp) :: values of the powerlaw or beta at x

  !    HISTORY
  !>       \authors Tam Nguyen

  !>       \date June 2020

  function eval_sas(x, sas_function, ka, b)
    implicit none

    real(dp),     dimension(:),     intent(in)    :: x   
    integer(i4),                    intent(in)    :: sas_function  
    real(dp),                       intent(in)    :: ka      !parameter k (or a) of the powerlaw (or beta) function
    real(dp),                       intent(in)    :: b       !parameter b of the beta function

    !local variables
    real(dp), dimension(size(x))                  :: eval_sas
    integer(i4)                                   :: i, ifault

    !check type of the SAS function
    if (sas_function == 1) then        !powerlaw
       goto 10
    else if (sas_function == 2) then   !beta function
       go to 20
    else                               !unknown funciton
       stop
    end if

    !cummulative powerlaw distribution function
    !if the SAS function is the powerlaw, then ka is the k parameter
10  continue
    eval_sas = x(:) ** ka  !cummulative SAS function
    go to 30

    !cummulative beta distribution function
    !if SAS function is the beta function, ka and b are a, b
20  continue
    do i = 1, size(x)
       eval_sas(i) = cdfbeta(x(i), ka, b)
    end do 

30  continue 

  end function eval_sas

!***************************
  subroutine add_vector(a, b, temp)

    implicit none
    real(dp), dimension(:), intent(inout) :: a, b
    real(dp), dimension(:), allocatable, intent(out) :: temp

    integer(i4) :: n, m

    n = max(size(a), size(b))
    m = min(size(a), size(b))

    allocate(temp(n))

    temp(:) = 0.0_dp
 
    temp(1:m) = a(1:m) + b(1:m)

    if (size(a) < n) temp((m + 1):n) = a(size(a)) + b((m + 1):n)
    if (size(b) < n) temp((m + 1):n) = b(size(b)) + a((m + 1):n)

  end subroutine 

!******************** N mass balance in the soil zone
  subroutine n_soil_balance(N0, Nin, Nout, alpha)

    implicit none
    
    real(dp), intent(in)    :: Nin, alpha
    real(dp), intent(out)   :: Nout
    real(dp), intent(inout) :: N0
 
    !local variable
    real(dp) :: Nt  !N mass at the end of the timestep
    
   !Mass balance equation
   !     dN(t)/dt = N(t-1) + Nin(t) - Nout(t)
   !     Nout(t) = alpha * N(t)

   !Analytical solution with delta_t = 1 [unit of 1/alpha])
   !     N(t) = (N(t-1) - Nin/alpha)*exp(-alpha) + Nin/alpha
   !     Nout(t) = N(t) - N(0) + Nin

   if (alpha .le. 1.0e-10_dp) then
     Nt = N0 + Nin 
     Nout = 0.0_dp
   else
     Nt = (N0 - Nin/alpha) * exp(-alpha) + Nin/alpha
     Nout = N0 + Nin - Nt  
   end if

   N0 = Nt

  end subroutine n_soil_balance

!-----------------------------------------------------------------------------------------------
! Subroutine for reading time series output and model setting in data.txt file
!-----------------------------------------------------------------------------------------------
subroutine read_input(fName, max_age, max_old_fraction, sas_function, &
                      S0, C0, N0, Qmax,                               &
                      year, Qin, Nsurplus, Cobs, Nwwtp)
  implicit none

  integer(i4),                            intent(out) :: max_age, sas_function
  real(dp),                               intent(out) :: max_old_fraction, S0, C0, N0, Qmax
  integer(i4), dimension(:), allocatable, intent(out) :: year
  real(dp),    dimension(:), allocatable, intent(out) :: Qin, Nsurplus, Cobs, Nwwtp
  character(50),                          intent(in)  :: fName                       !Including path

  !Local variables
  integer(i4) :: i, counter

  !Count number of lines
  counter = 0
  open(10, file = trim(fName))

    do 
      read(10, * , end = 100)
      counter = counter + 1
    end do

  100 continue  
  close(10)

  !Length of observed time series data
  counter = counter - 8 

  !Allocate time series data
  allocate(year(counter))
  allocate(Qin(counter))
  allocate(Nsurplus(counter))
  allocate(Cobs(counter))
  allocate(Nwwtp(counter))

  !Read data 
  open(10, file = fName)
    !Read initial condition and model setting

    read(10, *) max_age
    read(10, *) max_old_fraction
    read(10, *) sas_function
    read(10, *) S0
    read(10, *) C0
    read(10, *) N0
    read(10, *) Qmax

    !Read header of time series data
    read(10, *)

    !Start reading time series data
    do i = 1, counter
      read(10, *) year(i), Qin(i), Nsurplus(i), Cobs(i), Nwwtp(i)
    end do

    Qmax = Qmax * maxval(Qin)

  close(10)

end subroutine read_input


!-----------------------------------------------------------------------------------------------
! Function to calculate root mean square error
!-----------------------------------------------------------------------------------------------

function rmse(obs, sim)

  implicit none

  integer(i4)            :: i, counter
  real(dp)               :: rmse
  real(dp), dimension(:) :: sim, obs
  
  rmse = 0.0_dp

  counter = 0

  do i = 1, size(obs)
    if(obs(i) .ge. 0) then
      counter = counter + 1
      rmse = rmse + (obs(i) - sim(i))**2
    end if
  end do 
  
  rmse = sqrt(rmse/counter)

end function rmse

!-----------------------------------------------------------------------------------------------
! Subroutine to read parameter file
!-----------------------------------------------------------------------------------------------
subroutine read_param(fName, alpha, beta, a, b, half_life, f_instream)

  implicit none

  character(50), intent(in)                        :: fName
  real(dp), dimension(:), allocatable, intent(out) :: alpha, beta, a, b, half_life, f_instream

  integer(i4) :: i, counter

  !Count number of lines
  counter = 0
  open(10, file = trim(fName))

    do 
      read(10, * , end = 100)
      counter = counter + 1
    end do

  100 continue  
  close(10)

  !Allocate variable
  allocate(alpha(counter))
  allocate(beta(counter)) 
  allocate(a(counter)) 
  allocate(b(counter)) 
  allocate(half_life(counter))
  allocate(f_instream(counter))

  open(10, file = trim(fName))
  do i = 1, counter
    read(10, *) alpha(i), beta(i), a(i), b(i), half_life(i), f_instream(i)
  end do 
  close(10)


end subroutine read_param

END MODULE mo_sas_utils






