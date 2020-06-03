! **********************************************************
! Information on calculating Raman scattering cross sections
! **********************************************************
!
! For the differential Raman cross section equation, see:
!
! Reiher, M.; Neugebauer, J.; and Hess, B. A.  Z. Phys. Chem., 217, 91 
!   (2003).
! Neugebauer, J.; Reiher, M.; Kind, C.; and Hess, B. A.  
!   J. Comp. Chem., 23, 895, (2002).
! Long, Derek A. "The Raman Effect".  John Wiley & Sons: London.
!
! The equation for the differential cross section is:
!
! d(sigma)/d(Omega) = (pi^2/epsilon_0^2)*(lfreq-freq(i))^4*rcrs(i)*N_{p^I}
!
! This assumes that data is collected at a 90 degree angle from the 
! incident radiation path.  The equation has:
!
! epsilon_0 is the vacuum permittivity
! lfreq is the frequency of the incident radiation (laser freq.)
! freq(i) is the frequency of the normal mode
! b = h/(8*pi^2*c*freq(i)) <- Conversion between mass-weighted and 
! dimensionless coordinates, which is not included if you are dealing 
! with general polarizability invariants.
! rcrs(i) = (45*a_i^2+7*gamma_i^2)/45
! N_{p^I} = g_{p^I}*exp[-epsilon_{p^I}/k_BT] / R_{p}
! R_{p} = \sum_{J} g_{p^J}*exp[-epsilon_{p^J}/k_BT]
! 
! In other words, most terms are constants in this equation.  N_{p^I} is the
! population of molecules modeled as harmonic oscillators in the 
! vibrational state I for mode Q_p.  The temperature is initially chosen to
! be 300 K.  According to Long's book, the degeneracy is needed when 
! dealing with differential polarizability invariants, but those are not 
! what is calculated in this code.

subroutine RamanCrossSection(nmodes, boltzpop, lfreq, freq, rcrs)

  use constants
 
  implicit none

  integer i, nmodes

  real(kindr), intent(in) :: boltzpop(nmodes)
  real(kindr), intent(in) :: lfreq
  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(inout) :: rcrs(nmodes) 
  real(kindr) :: frq4th, tempfc

! =============================================================================
! purpose: Evaluate differential Raman cross section for fundamentals
!
! input  : nmodes - number of normal modes
!          boltzpop - Boltzmann population for the normal modes
!          lfreq - laser frequency
!          freq - vibrational frequencies
!
! in/out : rcrs - Differential Raman cross section
!
! =============================================================================

  do i=1,nmodes
    tempfc=boltzpop(i)
    frq4th = (lfreq-freq(i))**4

    ! The factor of 45 is put into the equation for the cross section
    ! when taking the average of the isotropic and anisotropic 
    ! polarizability.  The factor of 1E+12 (100^6) comes from 
    ! converting the cm from wavenumbers to meters for cancellation 
    ! purposes, and putting the differential cross section into units 
    ! of cm^2/sr. 

    ! rrscrsfac = pi^2*(100 cm/1 m)^6/[45*(8.8541878E-12 C V^{-1} m^{-1})^2]

    rcrs(i) = rrscrsfac*tempfc*frq4th*rcrs(i)
  end do

end subroutine RamanCrossSection

subroutine RamanCrossSectionOverCB(excnum, numcb, nmodes, freq, freqcb, &
                     lfreq, boltzpopot, rcrsot, boltzpopcb, rcrscb)

  use constants

  implicit none

  integer :: excnum, numcb, nmodes, i, k

  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(in) :: freqcb(numcb)
  real(kindr), intent(in) :: lfreq
  real(kindr), intent(in) :: boltzpopot(excnum,nmodes)
  real(kindr), intent(in) :: boltzpopcb(numcb)
  real(kindr), intent(inout) :: rcrsot(excnum, nmodes)
  real(kindr), intent(inout) :: rcrscb(numcb)
  real(kindr) :: frq4th, tempfc

! =============================================================================
! purpose: Evaluate differential Raman cross section for overtones/comb. bands
!
! input  : excnum - excitation number
!          numcb - number of combination bands
!          nmodes - number of normal modes
!          freq - vibrational frequencies
!          freqcb - vibrational frequencies for combination bands
!          lfreq - laser frequency
!          boltzpopot - Boltzmann population for the overtones 
!          boltzpopcb - Boltzmann population for the comb. bands 
!
! in/out : rcrsot - Differential Raman cross section for overtones
!          rcrscb - Differential Raman cross section for comb. bands
!
! =============================================================================

  ! Overtone part
  do k=2,excnum
    do i=1,nmodes
      tempfc=boltzpopot(k,i)
      frq4th = (lfreq-k*freq(i))**4
      rcrsot(k,i) = rrscrsfac*tempfc*frq4th*rcrsot(k,i)
    end do
  end do
  ! Combination band part part
  do i=1,numcb
    tempfc=boltzpopcb(i)
    frq4th = (lfreq-freqcb(i))**4
    rcrscb(i) = rrscrsfac*tempfc*frq4th*rcrscb(i)
  end do

end subroutine RamanCrossSectionOverCB

! **********************************************************************
! Information on calculating anti-Stokes Raman scattering cross sections
! **********************************************************************
!
! For the differential cross section equation, see:
!
! Long, Derek A. "The Raman Effect".  John Wiley & Sons: London.
!
! The equation for the differential cross section is:
!
! d(sigma)/d(Omega) = (pi^2/epsilon_0^2)*(lfreq+freq(i))^4*rcrs(i)*N_{p^I}
!
! This assumes that data is collected at a 90 degree angle from the 
! incident radiation path.  The equation has:
!
! epsilon_0 is the vacuum permittivity
! lfreq is the frequency of the incident radiation (laser freq.)
! freq(i) is the frequency of the normal mode
! b = h/(8*pi^2*c*freq(i)) <- Conversion between mass-weighted and 
! dimensionless coordinates, which is not included if you are dealing 
! with general polarizability invariants.
! rcrs(i) = (45*a_i^2+7*gamma_i^2)/45
! N_{p^I} = g_{p^I}*exp[-epsilon_{p^I}/k_BT] / R_{p}
! R_{p} = \sum_{J} g_{p^J}*exp[-epsilon_{p^J}/k_BT]
! 
! In other words, most terms are constants in this equation.  N_{p^I} is the
! population of molecules modeled as harmonic oscillators in the 
! vibrational state I for mode Q_p.  The temperature is initially chosen to
! be 300 K.  According to Long's book, the degeneracy is needed when 
! dealing with differential polarizability invariants, but those are not 
! what is calculated in this code.

subroutine ASRamanCrossSection(nmodes, boltzpopas, lfreq, freq, rcrs)

  use constants

  implicit none

  integer :: nmodes, i

  real(kindr), intent(in) :: boltzpopas(nmodes)
  real(kindr), intent(in) :: lfreq
  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(inout) :: rcrs(nmodes)
  real(kindr) :: frq4th, tempfc

  do i=1,nmodes
    tempfc=boltzpopas(i)
    frq4th = (lfreq+freq(i))**4

    ! The factor of 45 is put into the equation for the cross section
    ! when taking the average of the isotropic and anisotropic 
    ! polarizability.  The factor of 1E+12 (100^6) comes from 
    ! converting the cm from wavenumbers to meters for cancellation 
    ! purposes, and putting the differential cross section into units 
    ! of cm^2/sr. 

    ! rrscrsfac = pi^2*(100 cm/1 m)^6/[45*(8.8541878E-12 C V^{-1} m^{-1})^2]

    rcrs(i) = rrscrsfac*tempfc*frq4th*rcrs(i)
  end do

end subroutine ASRamanCrossSection
