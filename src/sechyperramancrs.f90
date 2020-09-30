! ***********************************************************************
! Information on calculating second hyper-Raman scattering cross sections
! ***********************************************************************
!
! Second hyper-Raman differential scattering cross sections can be 
! evaluated from their relationship with intensity.  This intensity is
! given by:
!
! I^{SHRS} = (pi^2*c*E_0^6 / 2*epsilon_0)*<gamma^2>
!
! Here, c is the speed of light, E_0 is the incident field strength,
! and epsilon_0 is the vacuum permittivity.  The intensity can instead 
! be expressed by the incident irradiance, I_0,
!
! I_0 = 0.5*epsilon_0*c*E_0^2
!
! Such that the second hyper-Raman scattering intensity is:
!
! I^{SHRS} = (pi^2*c  / 2*epsilon_0)*(8 / epsilon_0^3*c^3)*I_0^3*<gamma^2>
!          = (4*pi^2  / epsilon_0^4*c^2)*I_0^3*<gamma^2>
!
! We can now insert the definition of the vacuum permittivity:
!
! epsilon_0 = e^2/(2*alphaf*h*c)
!
! This gives:
!
! I^{SHRS} = (4*pi^2  / c^2)*(16*alphaf^4*h^4*c^4/e^8)*I_0^3*<gamma^2>
!          = (64*pi^2*alphaf^4*h^4*c^2 / e^8)*I_0^3*<gamma^2>
!
! We can then use that the intensity of second hyper-Raman scattering is
! given as:
!
! I^{SHRS} = d(sigma)/d(Omega)*I_0^3
!
! So that the differential second hyper-Raman scattering cross section is
!
! d(sigma)/d(Omega) = (64*pi^2*alphaf^4*h^4*c^2 / e^8)*<gamma^2>
!
! This can also be written:
!
! d(sigma)/d(Omega) = 
!  (64*pi^2*planck^4*alphaf^4/c^2*e^8)*N_{p^I}*(lfreq-freq(i))^4*rcrs(i)
!
! alphaf is the fine structure constant
! h is Planck's constant (6.626E-34 J s)
! c is the speed of light (2.9979E8 m/s)
! lfreq is the frequency of the incident radiation (laser freq.) in Hz
! freq(i) is the frequency of the normal mode in Hz
! rcrs(i) is the averaged square of the hyperpolarizability
! N_{p^I} = g_{p^I}*exp[-epsilon_{p^I}/k_BT] / R_{p}
! R_{p} = \sum_{J} g_{p^J}*exp[-epsilon_{p^J}/k_BT]

subroutine SecHypRamanCrossSection(nmodes, boltzpop, lfreq, freq, rcrs)

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

    ! Calculate second hyper-Raman cross sections using SI units, which 
    ! give the final value in units of: cm^6 s^2 / J^2.  

    ! Note that the extra factor of 100^6 (1E+12) is required for 
    ! converting from m^6 s^2 / J^2 to cm^6 s^2 / J^2.  This can be converted
    ! to photon normalized units (cm^6 s^2 / photon^2) by multiplying the 
    ! expression by 3x the photon energy squared (E = hv = hcv', v' is the 
    ! frequency in wavenumbers). 

    ! r2ndhrscrstmp = 64*pi^2*planck^4*alphaf^4*(100 cm/1 m)^6/(c^2*e^8)
    ! photonnorm2 = 3*h*c*(100 cm/1 m)
    ! r2ndhrscrsfac = r2ndhrscrstmp*photonnorm2*photonnorm2

    rcrs(i) = r2ndhrscrsfac*tempfc*frq4th*rcrs(i)*lfreq*lfreq
  end do

end subroutine SecHypRamanCrossSection
