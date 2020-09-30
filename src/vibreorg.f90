subroutine VibrationalReorg(freq, delta, vibreorg, nmodes, nexci)

  use constants

  implicit none

  character(80) :: frmt

  integer :: nmodes, i, j, nexci

  real(kindr), intent(in) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr), intent(out) :: vibreorg(nexci)

! =============================================================================
! purpose: Determine the total vibrational (internal) reorganization energy
!
! input  : freq - vibrational frequencies
!          delta - dimensionless displacements
!          nmodes - number of normal modes
!          nexci - number of excited states
!
! output : vibreorg - total vibrational reorganization energy for a particular
!                   excited state
!
! =============================================================================

  ! The total vibrational reorganization energy is calculated as:
  ! lambda_v = sum_{i=1}^{nmodes} freq*Delta^2/2

  vibreorg = zero

  do i=1,nexci
    do j=1,nmodes
      vibreorg(i) = vibreorg(i) + f12*freq(j)*delta(i,j)**2
    end do
  end do

  ! Print the vibrational reorganization energy

  write(*,*) "# Excited State       Vib. Reorganization Energy (cm-1)"
  frmt = "(1X,1A,10X,I2,10X,F8.2)"
  do i=1,nexci
    write(*,frmt) "#", i, vibreorg(i)
  end do

end subroutine VibrationalReorg
