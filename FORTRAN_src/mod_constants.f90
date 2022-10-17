module mod_constants

  use, intrinsic :: iso_fortran_env
  use fftw3
  implicit none

  ! kind parameters

  integer, parameter :: dp = real64

  real, parameter :: pi = 4*atan(1.0)
  real, parameter :: twopi = 2*pi
  real(dp), parameter :: third_dp = 1.0_dp / 3.0_dp
  real(dp), parameter :: sixth_dp = 1.0_dp / 6.0_dp
  real(dp), parameter :: pi_dp = 4*atan(1.0_dp)
  real(dp), parameter :: third_pi_dp = pi_dp / 3.0_dp

end module mod_constants
