subroutine StokesLorentzianConv(npoints, nmodes, excnum, numcb, &
                           freq_foc, rcrs_foc, width, maxfrq)
 
  use constants

  implicit none

  integer :: i, j, npoints, nmodes, excnum, numcb

  real(kindr) :: frqmin, frqmax, dfrq, afac, bfac
  real(kindr) :: conv, frq, factor
  real(kindr), intent(in) :: width
  real(kindr), intent(in) :: maxfrq
  real(kindr), intent(in) :: freq_foc(excnum*nmodes+numcb)
  real(kindr), intent(in) :: rcrs_foc(excnum*nmodes+numcb) 

! =============================================================================
! purpose: Convolute a Stokes spectrum (RRS, RHRS, SFG) with a Lorentzian
!
! input  : nmodes - number of normal modes
!          excnum - excitation number for overtones
!          numcb - number of combination bands
!          freq_foc - vibrational frequencies
!          rcrs_foc - cross sections
!          width - width of the Lorentzian
!          maxfrq - maximum vibrational frequency
!
! =============================================================================

  ! Depending on the type of calculation, frqmax may be in the
  ! fingerprint region or very high in wavenumbers.

  frqmin = zero
  frqmax = maxfrq + hndrtwnty

  ! Energy steps determined by NPOINTS

  dfrq = (frqmax-frqmin)
  dfrq = dfrq/(excnum*npoints-1)

  ! Calculate the factors used in the convolution, using the width 
  ! defined above.  The full width at half maximum is user defined.

  afac = width/(two*pi)
  bfac = width/two

  ! Loop over npoints
  do i=0,npoints*excnum
    conv = zero
    frq = frqmin + i*dfrq
    do j=1,excnum*nmodes+numcb
      factor = afac/( (frq-freq_foc(j))**2 + bfac**2)
      conv = conv + factor*rcrs_foc(j)
    end do
    ! The factor of pi that normalizes the Lorentzian is built into
    ! afac.  In the definition of a normalized Lorentzian, afac is
    ! actually width/2, but the whole expression is multiplied by
    ! 1/pi.  Since we multiply 1/pi into afac, we don't need to do
    ! this again.
    write(*,*) frq, conv
  end do

end subroutine StokesLorentzianConv

subroutine AntiStokesLorentzianConv(npoints, nmodes, freq_foc, &
                               rcrs_foc, width, minfrq, maxfrq)

  use constants

  implicit none

  integer :: i, j, npoints, nmodes

  real(kindr) :: frqmin, frqmax, dfrq, afac, bfac
  real(kindr) :: conv, frq, factor
  real(kindr), intent(in) :: width
  real(kindr), intent(in) :: minfrq
  real(kindr), intent(in) :: maxfrq
  real(kindr), intent(in) :: freq_foc(nmodes)
  real(kindr), intent(in) :: rcrs_foc(nmodes)

! =============================================================================
! purpose: Convolute a anti-Stokes spectrum (AS-RRS, AS-RHRS) with a Lorentzian
!
! input  : nmodes - number of normal modes
!          excnum - excitation number for overtones
!          numcb - number of combination bands
!          freq_foc - vibrational frequencies
!          rcrs_foc - cross sections
!          width - width of the Lorentzian
!          minfrq - minimum vibrational frequency (most negative)
!          maxfrq - maximum vibrational frequency (least negative)
!
! =============================================================================

  ! For anti-Stokes calculations, frqmax should be in the fingerprint 
  ! region due to signal weakness.

  frqmin = minfrq - hndrtwnty
  frqmax = maxfrq + hndrtwnty

  ! Energy steps determined by NPOINTS

  dfrq = (frqmax-frqmin)
  dfrq = dfrq/(npoints-1)

  ! Calculate the factors used in the convolution, using the width 
  ! defined above.  The full width at half maximum is user defined. 

  afac = width/(two*pi)
  bfac = width/two

  ! Loop over npoints
  do i=0,npoints
    conv = zero
    frq = frqmin + i*dfrq
    do j=1,nmodes
      factor = afac/( (frq-freq_foc(j))**2 + bfac**2)
      conv = conv + factor*rcrs_foc(j)
    end do
    ! The factor of pi that normalizes the Lorentzian is built into
    ! afac.  In the definition of a normalized Lorentzian, afac is
    ! actually width/2, but the whole expression is multiplied by
    ! 1/pi.  Since we multiply 1/pi into afac, we don't need to do
    ! this again.
    write(*,*) frq, conv
  end do

end subroutine AntiStokesLorentzianConv

subroutine NormStokesLorentzianConv(npoints, nmodes, excnum, numcb, &
                               freq_foc, rcrs_foc, width, maxfrq, &
                               norm_frq, norm_crs, rmax)

  use constants

  implicit none

  integer :: i, j, npoints, nmodes, excnum, numcb

  real(kindr) :: frqmin, frqmax, dfrq, afac, bfac
  real(kindr) :: conv, frq, factor, rmax
  real(kindr), intent(in) :: width
  real(kindr), intent(in) :: maxfrq
  real(kindr), intent(in) :: freq_foc(excnum*nmodes+numcb)
  real(kindr), intent(in) :: rcrs_foc(excnum*nmodes+numcb)
  real(kindr) :: norm_frq(npoints*excnum+1)
  real(kindr) :: norm_crs(npoints*excnum+1)

! =============================================================================
! purpose: Convolute a Stokes spectrum (RRS, RHRS, SFG) with a Lorentzian, 
!          then normalize it
!
! input  : nmodes - number of normal modes
!          excnum - excitation number for overtones
!          numcb - number of combination bands
!          freq_foc - vibrational frequencies
!          rcrs_foc - cross sections
!          width - width of the Lorentzian
!          maxfrq - maximum vibrational frequency
!          norm_frq - vibrational frequencies for plotting normalized
!                     spectra
!          norm_crs - normalized cross sections
!
! output : rmax - maximum value of the cross section for normalization
!
! =============================================================================

  ! Depending on the type of calculation, frqmax may be in the
  ! fingerprint region or very high in wavenumbers.

  frqmin = zero
  frqmax = maxfrq + hndrtwnty

  ! Energy steps determined by NPOINTS

  dfrq = (frqmax-frqmin)
  dfrq = dfrq/(excnum*npoints-1)

  ! Calculate the factors used in the convolution, using the width 
  ! defined above.

  afac = width/(two*pi)
  bfac = width/two

  ! Loop over npoints
  do i=0,npoints*excnum
    conv = zero
    frq = frqmin + i*dfrq
    do j=1,excnum*nmodes+numcb
      factor = afac/( (frq-freq_foc(j))**2 + bfac**2)
      conv = conv + factor*rcrs_foc(j)
    end do
    norm_frq(i+1) = frq
    norm_crs(i+1) = conv
  end do
  rmax = MAXVAL(norm_crs)
  do i=0,npoints*excnum
    norm_crs(i+1) = norm_crs(i+1)/rmax
    write(*,*) norm_frq(i+1), norm_crs(i+1)
  end do

end subroutine NormStokesLorentzianConv
