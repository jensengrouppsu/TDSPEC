subroutine StokesShift(freq, delta, stokes, stokessolv, freqcut, &
                  nmodes, nexci, fluorstate)

  use constants

  implicit none

  integer, intent(in) :: nmodes, nexci
  integer, intent(in) :: fluorstate
  integer :: i 
 
  real(kindr), intent(in) :: freq(nmodes) 
  real(kindr), intent(in) :: delta(nexci,nmodes) 
  real(kindr), intent(in) :: freqcut, stokessolv
  real(kindr) :: stokesvib

  real(kindr), intent(out) :: stokes

! ======================================================================
! purpose: Evaluate the Stokes shift for fluorescence
!
! input  : freq - vibrational frequencies
!          delta - dimensionless displacements
!          freqcut - low-frequency cut-off
!          fluorstate - state used for fluorescence
!
! out : stokes - Stokes shift (IMDHO method)
!
! ======================================================================

! This routine evaluates the Stokes shift based on the IMDHO method.  
! The equation for the vibrational contribution to this is (when 
! Duschinsky effects are absent):
!
! lambda_v = sum_f 2*s_f*omega_f = sum_f delta_f^2*omega_f
!
! For the "solvent" contribution to the Stokes shift, the code 
! currently treats this as an adjustable parameter related to the 
! inhomogeneous broadening parameter, such that the total Stokes shift
! including "vibrational" and "solvent" contribtutions is
!
! lambda = lambda_v + lambda_s
!        = sum_f delta_f^2*omega_f + theta
!
! For details on evaluating the Stokes shift, see:
!
! T. Petrenko and F. Neese.  J. Chem. Phys., 127, 164319 (2007).
! T. Petrenko, O. Krylova, F. Neese, and M. Sokolowski.  New J. Phys.,
!    11, 015001 (2009).

  ! Stokes shift due to vibrations
  stokesvib = zero

  ! The variable freqcut is chosen such that it is equal to:
  ! 2.0*kb*T = 417 cm^{-1}
  ! 
  ! This can be adjusted in the input file.
  
  do i=0,nmodes
    if (freq(i) < freqcut) then
      stokesvib = stokesvib + delta(fluorstate,i)**2*freq(i)
    end if  
  end do
   
  ! Calculate the total Stokes shift
  stokes = stokesvib + stokessolv

end subroutine StokesShift
