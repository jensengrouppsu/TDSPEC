! ****************************************************************
! Information on calculating hyper-Raman scattering cross sections
! ****************************************************************
!
! For the hyper-Raman differential scattering cross section, see
!
! Craig, D. P. and Thirunamachandran, T.  "Molecular Quantum 
!   Electrodynamics: An Introduction to Radiation Molecule 
!   Interactions".  John Wiley & Sons: West Sussex, England, 2002.
! Kelley, A. M.  J. Phys. Chem. A, 112, 11975, (2008).
! Ziegler, L. D.  J. Raman. Spec., 21, 769 (1990).
!
! Note that the latter two references write the expressions with the
! transition length (i.e. a0), not the transition dipole moment 
! (i.e. e*a0).  Formally, their equations are correct, but they are
! confusing if you think in terms of most quantum chemistry programs.
! The equation for the hyper-Raman differential scattering cross 
! section is (according to Kelley and Ziegler)
!
! d(sigma)/d(Omega) = 
!             (16*pi^2*alphaf^3/N*h*c^2)*N_{p^I}*(lfreq-freq(i))^4*rcrs(i)
!
! This expresion, when written with transition dipole moments yields
!
! d(sigma)/d(Omega) = 
!  (16*pi^2*planck^3*alphaf^3/N*c^2*e^6)*N_{p^I}*(lfreq-freq(i))^4*rcrs(i)
!
! Depending on the setup, rcrs(i) is averaged differently so that it 
! matches the angle that data is collected at with respect to the 
! incident radiation.
!
! alphaf is the fine structure constant
! N is the number of scatterers per unit volume
! h is Planck's constant (6.626E-34 J s)
! c is the speed of light (2.9979E8 m/s)
! lfreq is the frequency of the incident radiation (laser freq.) in Hz
! freq(i) is the frequency of the normal mode in Hz
! rcrs(i) is the averaged square of the hyperpolarizability
! N_{p^I} = g_{p^I}*exp[-epsilon_{p^I}/k_BT] / R_{p}
! R_{p} = \sum_{J} g_{p^J}*exp[-epsilon_{p^J}/k_BT]

subroutine HypRamanCrossSection(nmodes, boltzpop, lfreq, freq, rcrs)

  use constants

  implicit none

  integer :: nmodes, i
 
  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(in) :: lfreq
  real(kindr), intent(in) :: boltzpop(nmodes)
  real(kindr), intent(inout) :: rcrs(nmodes)
  real(kindr) :: frq4th, tempfc 

! =============================================================================
! purpose: Evaluate differential hyper-Raman cross section for fundamentals
!
! input  : nmodes - number of normal modes
!          boltzpop - Boltzmann population for the normal modes
!          lfreq - laser frequency
!          freq - vibrational frequencies
!
! in/out : rcrs - Differential hyper-Raman cross section
!
! =============================================================================

  do i=1,nmodes
    frq4th = (hndr*speed*(lfreq-freq(i)))**4

    tempfc = boltzpop(i)

    ! Calculate hyper-Raman cross sections using SI units, which 
    ! give the final value in units of: cm^4 s / J.  

    ! Note that the extra factor of 100^4 (1E+8) is required for 
    ! converting from m^4 s / J to cm^4 s / J.  This can be converted
    ! to photon normalized units (cm^4 s / photon) by multiplying the 
    ! expression by 2x the photon energy (E = hv = hcv', v' is the 
    ! frequency in wavenumbers). 

    ! rhrscrstmp = 16*pi^2*h^3*alpha^3*(100 cm/1 m)^4/(c^2*e^6)
    ! photonnorm = 2*h*c*(100 cm/1 m)
    ! rhrscrsfac = rhrscrstmp*photonnorm

    rcrs(i) = rhrscrsfac*tempfc*frq4th*rcrs(i)*lfreq
  end do

end subroutine HypRamanCrossSection

subroutine HypRamanCrossSectionOverCB(excnum, numcb, nmodes, freq, freqcb, &
                     lfreq, boltzpopot, rcrsot, boltzpopcb, rcrscb)

  use constants

  implicit none

  integer :: excnum, numcb, nmodes, k, i

  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(in) :: freqcb(numcb)
  real(kindr), intent(in) :: lfreq
  real(kindr), intent(in) :: boltzpopot(excnum,nmodes)
  real(kindr), intent(inout) :: rcrsot(excnum,nmodes)
  real(kindr), intent(in) :: boltzpopcb(numcb)
  real(kindr), intent(inout) :: rcrscb(numcb)
  real(kindr) :: frq4th, tempfc

! =============================================================================
! purpose: Evaluate differential hyper-Raman cross section for overtones/comb. 
!          bands
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
! in/out : rcrsot - Differential hyper-Raman cross section for overtones
!          rcrscb - Differential hyper-Raman cross section for comb. bands
!
! =============================================================================

  ! Overtone part
  do k=2,excnum
    do i=1,nmodes
      frq4th = (hndr*speed*(lfreq-k*freq(i)))**4
      tempfc=boltzpopot(k,i)
      rcrsot(k,i) = rhrscrsfac*tempfc*frq4th*rcrsot(k,i)*lfreq
    end do
  end do
  ! Combination band part
  do i=1,numcb
    frq4th = (hndr*speed*(lfreq-freqcb(i)))**4
    tempfc = boltzpopcb(i)
    rcrscb(i) = rhrscrsfac*tempfc*frq4th*rcrscb(i)*lfreq
  end do

end subroutine HypRamanCrossSectionOverCB

! ****************************************************************************
! Information on calculating anti-Stokes hyper-Raman scattering cross sections
! ****************************************************************************
!
! For the hyper-Raman differential scattering cross section, see
!
! Craig, D. P. and Thirunamachandran, T.  "Molecular Quantum 
!   Electrodynamics: An Introduction to Radiation Molecule 
!   Interactions".  John Wiley & Sons: West Sussex, England, 2002.
! Kelley, A. M.  J. Phys. Chem. A, 112, 11975, (2008).
! Ziegler, L. D.  J. Raman. Spec., 21, 769 (1990).
!
! Note that the latter two references write the expressions with the
! transition length (i.e. a0), not the transition dipole moment 
! (i.e. e*a0).  Formally, their equations are correct, but they are
! confusing if you think in terms of most quantum chemistry programs.
! The equation for the hyper-Raman differential scattering cross 
! section is (according to Kelley and Ziegler)
!
! d(sigma)/d(Omega) = 
!             (16*pi^2*alphaf^3/N*h*c^2)*N_{p^I}*(lfreq-freq(i))^4*rcrs(i)
!
! This expresion, when written with transition dipole moments yields
!
! d(sigma)/d(Omega) = 
!  (16*pi^2*planck^3*alphaf^3/N*c^2*e^6)*N_{p^I}*(lfreq-freq(i))^4*rcrs(i)
!
! Depending on the setup, rcrs(i) is averaged differently so that it 
! matches the angle that data is collected at with respect to the 
! incident radiation.
!
! alphaf is the fine structure constant
! N is the number of scatterers per unit volume
! h is Planck's constant (6.626E-34 J s)
! c is the speed of light (2.9979E8 m/s)
! lfreq is the frequency of the incident radiation (laser freq.) in Hz
! freq(i) is the frequency of the normal mode in Hz
! rcrs(i) is the averaged square of the hyperpolarizability
! N_{p^I} = g_{p^I}*exp[-epsilon_{p^I}/k_BT] / R_{p}
! R_{p} = \sum_{J} g_{p^J}*exp[-epsilon_{p^J}/k_BT]

subroutine ASHypRamanCrossSection(nmodes, boltzpopas, lfreq, freq, rcrs)

  use constants

  implicit none

  integer :: i, nmodes

  real(kindr), intent(in) :: boltzpopas(nmodes)
  real(kindr), intent(in) :: lfreq
  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(inout) :: rcrs(nmodes)
  real(kindr) :: frq4th, tempfc

  do i=1,nmodes
    frq4th = (hndr*speed*(lfreq+freq(i)))**4

    tempfc = boltzpopas(i)

    ! Calculate hyper-Raman cross sections using SI units, which 
    ! give the final value in units of: cm^4 s / J.  

    ! Note that the extra factor of 100^4 (1E+8) is required for 
    ! converting from m^4 s / J to cm^4 s / J.  This can be converted 
    ! to photon normalized units (cm^4 s / photon) by multiplying the 
    ! expression by 2x the photon energy (E = hv = hcv', v' is the 
    ! frequency in wavenumbers). 

    ! rhrscrstmp = 16*pi^2*h^3*alpha^3*(100 cm/1 m)^4/(c^2*e^6)
    ! photonnorm = 2*h*c*(100 cm/1 m)
    ! rhrscrsfac = rhrscrstmp*photonnorm

    rcrs(i) = rhrscrsfac*tempfc*frq4th*rcrs(i)*lfreq
  end do

end subroutine ASHypRamanCrossSection
