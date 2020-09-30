subroutine UnMassWeight(tdipder,freq,nmodes,nexci,ndim)

  use constants

  implicit none

  integer, intent(in) :: nexci, nmodes, ndim
  integer :: iex,j,k

  real(kindr), intent(inout) :: tdipder(nexci,nmodes,ndim)
  real(kindr), intent(in) :: freq(nmodes)

! =============================================================================
! purpose: Remove mass-weighting of transition dipole moment derivatives
!
! input  : freq - vibrational frequencies
!          nexci - number of excited states
!          nmodes - number of normal modes
!          ndim - dimension of the vector
!
! in/out : tdipder - un-mass-weighted transition dipole moment derivatives
!
! =============================================================================

  do iex=1,nexci
    do j=1,nmodes
      do k=1,ndim
         tdipder(iex,j,k) = (htmassweight/dsqrt(freq(j)))*tdipder(iex,j,k)
!      tdipder(iex,j,1) = (htmassweight/dsqrt(freq(j)))*tdipder(iex,j,1)
!      tdipder(iex,j,2) = (htmassweight/dsqrt(freq(j)))*tdipder(iex,j,2)
!      tdipder(iex,j,3) = (htmassweight/dsqrt(freq(j)))*tdipder(iex,j,3)
      end do
    end do
  end do

end subroutine UnMassWeight

subroutine UnMassWeightB2(stpmder,freq,nmodes,nexci,ndim)

  use constants

  implicit none

  integer :: nexci, nmodes, ndim
  integer :: iex,j

  real(kindr), intent(inout) :: stpmder(nexci,nmodes,ndim)
  real(kindr), intent(in) :: freq(nmodes)

! =============================================================================
! purpose: Remove mass-weighting of two-photon transition moment derivatives
!
! input  : freq - vibrational frequencies
!          nexci - number of excited states
!          nmodes - number of normal modes
!          ndim - dimension of the tensor (equals 6)
!
! in/out : stpmder - un-mass-weighted two-photon transition moment derivatives
!
! =============================================================================

  do iex=1,nexci
    do j=1,nmodes
      stpmder(iex,j,1) = (htmassweight/dsqrt(freq(j)))*stpmder(iex,j,1)
      stpmder(iex,j,2) = (htmassweight/dsqrt(freq(j)))*stpmder(iex,j,2)
      stpmder(iex,j,3) = (htmassweight/dsqrt(freq(j)))*stpmder(iex,j,3)
      stpmder(iex,j,4) = (htmassweight/dsqrt(freq(j)))*stpmder(iex,j,4)
      stpmder(iex,j,5) = (htmassweight/dsqrt(freq(j)))*stpmder(iex,j,5)
      stpmder(iex,j,6) = (htmassweight/dsqrt(freq(j)))*stpmder(iex,j,6)
    end do
  end do

end subroutine UnMassWeightB2

subroutine UnMassWeightSFG(gdipder,freq,nmodes,nexci,ndim)

  use constants

  implicit none

  integer :: nexci, nmodes, ndim
  integer :: iex,j

  real(kindr), intent(inout) :: gdipder(nexci,nmodes,ndim)
  real(kindr), intent(in) :: freq(nmodes)

! =============================================================================
! purpose: Remove mass-weighting of ground-state dipole moment derivatives
!
! input  : freq - vibrational frequencies
!          nexci - number of excited states
!          nmodes - number of normal modes
!          ndim - dimension of the vector (equals 3)
!
! in/out : gdipder - un-mass-weighted ground-state dipole moment derivatives
!
! =============================================================================

  do iex=1,nexci
    do j=1,nmodes
      gdipder(iex,j,1) = (htmassweight/dsqrt(freq(j)))*gdipder(iex,j,1)
      gdipder(iex,j,2) = (htmassweight/dsqrt(freq(j)))*gdipder(iex,j,2)
      gdipder(iex,j,3) = (htmassweight/dsqrt(freq(j)))*gdipder(iex,j,3)
    end do
  end do

end subroutine UnMassWeightSFG
