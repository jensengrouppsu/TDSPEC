subroutine Boltzmann(freq, boltzpop, nmodes, temper)

  use constants

  implicit none

  integer :: i
  integer, intent(in) :: nmodes

  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(in) :: temper
  real(kindr), intent(inout) :: boltzpop(nmodes)
  real(kindr) :: denom(nmodes)
  real(kindr) :: pop0(nmodes), pop1(nmodes), pop2(nmodes)
  real(kindr) :: pop3(nmodes), pop4(nmodes), pop5(nmodes)

! =============================================================================
! purpose: Calculate the Boltzmann factors for fundamentals
!
! input  : freq - normal mode frequencies
!          nmodes - number of normal modes
!          temper - temperature of the simulation
!
! output : boltzpop - Boltzmann population for each normal mode (fundamentals)
!
! =============================================================================

! The general Boltzmann population for normal modes is given as:
!
! N_{p^I} = g_{p^I}*exp[-E_{p^I}/k_B*T] / Q_p
! Q_p = sum_j g_{p^j}*exp[-E_{p^j}/k_B*T]
!
! Here, Q_p is the vibrational partition function of the pth normal 
! mode, and g_{p^j} is the vibrational degeneracy of the jth 
! vibrational state.  This expression is taken with small changes
! in notation from:
!
! Long, Derek A. "The Raman Effect".  John Wiley & Sons: London.

! Calculate the Boltzmann population of each normal mode for given 
! vibrational states.  Quantum numbers are given off to the side.  

  do i = 1,nmodes
    pop0(i) = exp(-f12*exparg*freq(i)/temper)  ! v = 0
    pop1(i) = exp(-f32*exparg*freq(i)/temper)  ! v = 1
    pop2(i) = exp(-f52*exparg*freq(i)/temper)  ! v = 2
    pop3(i) = exp(-f72*exparg*freq(i)/temper)  ! v = 3
    pop4(i) = exp(-f92*exparg*freq(i)/temper)  ! v = 4
    pop5(i) = exp(-f112*exparg*freq(i)/temper) ! v = 5
  end do

! Calculate the denominator (vibrational partition function) and 
! numerator and assume that the ground vibrational state is most 
! populated (~100%).  We truncate the infinite sum representing
! the vibrational partition function at the fifth excited state
! of the harmonic oscillator.

  do i = 1,nmodes
    denom(i) = pop0(i) + pop1(i) + pop2(i) + pop3(i) + pop4(i) + &
               pop5(i)
    boltzpop(i) = pop0(i) / denom(i)
  end do

end subroutine Boltzmann

subroutine BoltzmannOvrt(boltzpopot, boltzpop, nmodes, excnum)

  use constants

  implicit none

  integer :: i, j
  integer, intent(in) :: nmodes, excnum

  real(kindr), intent(in) :: boltzpop(nmodes)
  real(kindr), intent(inout) :: boltzpopot(excnum,nmodes)

! =============================================================================
! purpose: Calculate the Boltzmann factors for overtones
!
! input  : boltzpop - Boltzmann population of fundamentals
!          nmodes - number of normal modes
!          excnum - highest excitation level requested
!
! output : boltzpopot - Boltzmann population for each normal mode (overtones)
!
! =============================================================================

! Since the initial state is identical for overtones and fundamentals, we only
! need set indices of different arrays equal to each other.

  do i = 1,excnum
    do j = 1,nmodes
      boltzpopot(i,j) = boltzpop(j)
    end do
  end do

end subroutine BoltzmannOvrt

subroutine BoltzmannComb(freq, boltzpopcb, excnums, nmodes, numcb, temper)

  use constants

  implicit none

  integer :: i,j
  integer, intent(in) :: nmodes, numcb
  integer, intent(in) :: excnums(numcb,nmodes)

  real(kindr), intent(in) :: temper
  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(inout) :: boltzpopcb(numcb)
  real(kindr) :: denom(nmodes)
  real(kindr) :: pop0(nmodes), pop1(nmodes), pop2(nmodes)
  real(kindr) :: pop3(nmodes), pop4(nmodes), pop5(nmodes)
  real(kindr) :: temppop(nmodes)

! =============================================================================
! purpose: Calculate the Boltzmann factors for combination bands
!
! input  : freq - normal mode frequencies
!          excnums - excitation numbers for the combination bands
!          nmodes - number of normal modes
!          numcb - number of combination bands
!          temper - temperature of the simulation
!
! output : boltzpopcb - Boltzmann population for each combination band
!
! =============================================================================

! Calculate the Boltzmann population of each normal mode for given 
! vibrational states.  Quantum numbers are given off to the side.  

  do i = 1,nmodes
    pop0(i) = exp(-f12*exparg*freq(i)/temper)  ! v = 0
    pop1(i) = exp(-f32*exparg*freq(i)/temper)  ! v = 1
    pop2(i) = exp(-f52*exparg*freq(i)/temper)  ! v = 2
    pop3(i) = exp(-f72*exparg*freq(i)/temper)  ! v = 3
    pop4(i) = exp(-f92*exparg*freq(i)/temper)  ! v = 4
    pop5(i) = exp(-f112*exparg*freq(i)/temper) ! v = 5
  end do

! Calculate the denominator (vibrational partition function) and 
! numerator and assume that the ground vibrational state is most 
! populated (~100%).  We truncate the infinite sum representing
! the vibrational partition function at the fifth excited state
! of the harmonic oscillator.

  do i = 1,nmodes
    denom(i) = pop0(i) + pop1(i) + pop2(i) + pop3(i) + pop4(i) + &
               pop5(i)
    temppop(i) = pop0(i) / denom(i)
  end do

! Evaluate the Boltzmann population for the combination bands, assumed
! to be a product of the populations of each contribution normal mode.

  do i = 1,numcb
    boltzpopcb(i) = one
    do j = 1,nmodes
      if (excnums(i,j) /= 0) then
        boltzpopcb(i) = boltzpopcb(i)*temppop(j)
      end if
    end do
  end do

end subroutine BoltzmannComb

subroutine BoltzmannAS(freq, boltzpopas, nmodes, temper)

  use constants

  implicit none

  integer :: i
  integer, intent(in) :: nmodes

  real(kindr), intent(in) :: temper
  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(inout) :: boltzpopas(nmodes)
  real(kindr) :: denom(nmodes)
  real(kindr) :: pop0(nmodes), pop1(nmodes), pop2(nmodes)
  real(kindr) :: pop3(nmodes), pop4(nmodes), pop5(nmodes)

! =============================================================================
! purpose: Calculate the Boltzmann factors for anti-Stokes spectra
!
! input  : freq - normal mode frequencies
!          nmodes - number of normal modes
!          temper - temperature of the simulation
!
! output : boltzpopas - Boltzmann population for each normal mode (anti-Stokes)
!
! =============================================================================

! Calculate the Boltzmann population of each normal mode for given 
! vibrational states.  Quantum numbers are given off to the side.  

  do i = 1,nmodes
    pop0(i) = exp(-f12*exparg*freq(i)/temper)  ! v = 0
    pop1(i) = exp(-f32*exparg*freq(i)/temper)  ! v = 1
    pop2(i) = exp(-f52*exparg*freq(i)/temper)  ! v = 2
    pop3(i) = exp(-f72*exparg*freq(i)/temper)  ! v = 3
    pop4(i) = exp(-f92*exparg*freq(i)/temper)  ! v = 4
    pop5(i) = exp(-f112*exparg*freq(i)/temper) ! v = 5
  end do

! Calculate the denominator (vibrational partition function) and 
! numerator.  We truncate the infinite sum representing
! the vibrational partition function at the fifth excited state
! of the harmonic oscillator.  Since anti-Stokes spectra start in the
! first excited state, the Boltzmann population reflects that.

  do i = 1,nmodes
    denom(i) = pop0(i) + pop1(i) + pop2(i) + pop3(i) + pop4(i) + &
               pop5(i)
    boltzpopas(i) = pop1(i) / denom(i)
  end do

end subroutine BoltzmannAS
