Subroutine ASRRSFCFund(latermoff, x, w, freq, delta, osc, aspol, & 
                     freq0, lfreq, g, nexci, nmodes, tstep)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep

  complex(kindr), intent(inout) :: aspol(nmodes,3,3)
  complex(kindr) :: Ls, rsum, dsum, z
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep) 
  real(kindr) :: osc(nexci,3) 

  complex(kindr) :: zzero=(0d0,0d0) 

  logical :: latermoff

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate A term polarizability anti-Stokes RRS spectra
!
! input  : latermoff - determines if only B terms are calculated
!          x - abcissa points for numerical integration
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          osc - transition dipole moment components
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!
! in/out : aspol - anti-Stokes RRS polarizability tensor
!
! =============================================================================

  ! This is needed because some compilers don't default arrays to being zero.
  aspol = zzero

  if (latermoff.eqv..false.) then
    do j = 1, nmodes
      do iex=1,nexci
        Ls = zzero
          do t=1,tstep
            dsum = zzero
            do i=1,nmodes

              ! Calculate the sum term in the integral (redundant 
              ! part that applies both to Raman and absorption 
              ! spectra). 

              dsum = dsum + &
                    (delta(iex,i)**2/two)*(one-exp(-z*freq(i)*x(t)))
            end do

            ! Calculate the part in brackets { } of equation 19 in 
            ! the Neese and Petrenko paper, for k = 1.        

            ! Negative sign? 
            !rsum = -delta(iex,j)/sqrt(two)*(one-exp(-z*freq(j)*x(t)))
            ! Positive sign?
            rsum = delta(iex,j)/sqrt(two)*(one-exp(-z*freq(j)*x(t)))
            Ls = Ls + z*w(t)*rsum*exp(x(t)* &
              z*(lfreq-freq0(iex)+freq(j)) - g(iex,t) - dsum)
            ! z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - dsum)
          end do

        ! Without converting the units of the lineshape function the 
        ! polarizability has units of a.u.^2 cm since the transition 
        ! dipoles are in a.u. and the weighting function for the 
        ! integration is in cm.  For this part of the calculation, 
        ! the transition dipoles have units of charge times distance 
        ! (in a.u.).  The integral evaluated above has units of 
        ! inverse energy (1/cm^{-1} in this case).  What can be done 
        ! is a conversion of the units of the transition dipoles in
        ! a.u. to SI units.

        ! First, take the transition dipoles and convert them to SI
        ! units using: 1 Cm = 1.1794744E+29 a.u.  Then, it turns out 
        ! that the units of the integral are inverse energy (in this 
        ! case, 1/cm^{-1}).  These can be converted to Joules using 
        ! the fact that 1 cm^{-1} = 1.9864456E-23 J.   

        ! Convert transition dipoles to SI units (C*m) and the units 
        ! of the integral (1/cm^{-1}) to inverse Joules.  The reason 
        ! for doing this is that the SI unit of transition 
        ! polarizability is C^2*m^2/J, which can be written C*m^2/V 
        ! because 1 J = 1 C*V.

        ! pol2si = (1 C*m/1.1794744E+29 a.u.)^2*(1 cm^{-1}/1.9864456E-23 J) 

        aspol(j,1,1) = aspol(j,1,1) + osc(iex,1)*osc(iex,1)*Ls* &
                     pol2si
        aspol(j,1,2) = aspol(j,1,2) + osc(iex,1)*osc(iex,2)*Ls* &
                     pol2si
        aspol(j,1,3) = aspol(j,1,3) + osc(iex,1)*osc(iex,3)*Ls* &
                     pol2si
        aspol(j,2,1) = aspol(j,2,1) + osc(iex,2)*osc(iex,1)*Ls* &
                     pol2si
        aspol(j,2,2) = aspol(j,2,2) + osc(iex,2)*osc(iex,2)*Ls* &
                     pol2si
        aspol(j,2,3) = aspol(j,2,3) + osc(iex,2)*osc(iex,3)*Ls* &
                     pol2si
        aspol(j,3,1) = aspol(j,3,1) + osc(iex,3)*osc(iex,1)*Ls* &
                     pol2si
        aspol(j,3,2) = aspol(j,3,2) + osc(iex,3)*osc(iex,2)*Ls* &
                     pol2si
        aspol(j,3,3) = aspol(j,3,3) + osc(iex,3)*osc(iex,3)*Ls* &
                     pol2si
      end do
    end do
  end if

end subroutine ASRRSFCFund

Subroutine ASRRSHTFund(x, w, freq, delta, osc, tdipder, aspol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb, Ls_htb2, Ls_htb3)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep, r, k, n

  complex(kindr), intent(inout) :: aspol(nmodes,3,3)
  complex(kindr) :: rsum, rsum1, dsum, z
  complex(kindr) :: Ls_htb(nexci,nmodes), Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3),tdipder(nexci,nmodes,3)

  complex(kindr) :: zzero=(0d0,0d0)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B term polarizability for anti-Stokes RRS spectra
!
! input  : x - abcissa points for numerical integration 
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          osc - transition dipole moment components
!          tdipder - derivative of transition dipole moment components
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!
! in/out : aspol - anti-Stokes RRS polarizability tensor
!
! =============================================================================

! Calculations of the Herzberg-Teller terms for anti-Stokes resonance 
! Raman scattering are based on equations for time dependent 
! integrals implemented previously for Stokes scattering.  The only
! significant difference is that the sum of initial normal mode 
! frequencies (denoted epsilon_{I_0}) includes an excited normal
! mode (the index "j" is used in the code).  This contrasts Stokes
! scattering where the initial state energy is simply the sum of
! zero-point energies for every oscillator.  An additional frequency
! is included in nearly every lineshape as a result.

! Keep in mind that the initial and final states of the system
! are actually vectors containing the state of each harmonic 
! oscillator.  Since the overlap integral looks like sum_a<F|Qa|I>,
! each normal mode is raised and lowered, not just the mode in
! its excited state.

  do j = 1, nmodes
    do iex=1,nexci
      Ls_htb(1:nexci,1:nmodes) = zzero
      Ls_htb2(1:nexci,1:nmodes) = zzero
      Ls_htb3(1:nexci,1:nmodes) = zzero
      do t=1,tstep
        dsum = zzero
        do i=1,nmodes
          dsum = dsum + &
            (delta(iex,i)**2/two)*(one-exp(-z*freq(i)*x(t)))
        end do
        ! Calculate the contribution from the excited normal mode.
        ! This is either raised (rsum) or has nothing performed
        ! (rsum1).
        rsum = delta(iex,j)**2/(sqrt(eight))* &
               (one-exp(-z*freq(j)*x(t)))**2
        rsum1 = delta(iex,j)/(sqrt(two))* &
                (one-exp(-z*freq(j)*x(t)))
        ! Calculate the contribution from other normal modes.  This
        ! is either raised or has nothing performed, if the mode is
        ! the current mode.
        do r=1,nmodes
          if (r /= j) then
            expsum(r) = delta(iex,r)/sqrt(two)* &
                        (one-exp(-z*freq(r)*x(t)))
          else
            expsum(r) = one
          end if
        end do
        ! Calculate lineshape functions
        do n=1,nmodes
          ! Integrals of the form: <(F+1)|I(t)>
          if (n == j) then
            ! Raise the initial state of the current mode,
            ! so that I = F+1 = 1.  A factor of 1-delta^2 is
            ! applied for the 1-1 overlap of a normal mode.
            ! In anti-Stokes, the extra frequency from the
            ! lineshape function cancels the frequency of the
            ! initial state.
            Ls_htb(iex,n) = Ls_htb(iex,n) + &
                     z*(one-delta(iex,n)**2+ &
                     delta(iex,n)**2*cos(freq(n)*x(t)))* &
                     w(t)*exp(x(t)* &
                     z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
                   !  z*(lfreq-freq0(iex)) - conjg(g(iex,t)) - dsum)
          else
            ! Raise the final state of a different mode,
            ! so I = 0 and F+1 = 1.  Here, there's a product
            ! accounting for multiple modes in an excited state.
             Ls_htb(iex,n) = Ls_htb(iex,n) + z*w(t)*rsum1*expsum(n)* &
                      exp(x(t)* &
                      z*(lfreq-freq0(iex)+freq(j)) - g(iex,t) - &
                    !  z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
                      dsum)
          end if
          ! Integrals of the form: <F|(I-1)(t)>
          if (n == j) then
            ! Lower the initial state of the excited mode.
            Ls_htb2(iex,n) = Ls_htb2(iex,n) + z*w(t)*exp(x(t)* &
                      z*(lfreq-freq0(iex)+freq(j)) - g(iex,t) - &
                    !  z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
                      dsum)
          else
            ! Can't lower a harmonic oscillator past its ground
            ! state.
            Ls_htb2(iex,n) = zzero
          end if
          ! Integrals of the form: <F|(I+1)(t)> 
          if (n == j) then
            ! Raise the initial state of the current mode,
            ! so that F = 0 and I+1 = 2
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      z*w(t)*root2*rsum*exp(x(t)* &
                      z*(lfreq-freq0(iex)+freq(j)) - g(iex,t) - &
                    !  z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
                      dsum)
          else
            ! Raise the initial state of a different mode,
            ! so that F = 0 and I = 1.
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      z*w(t)*expsum(n)*rsum1*exp(x(t)* &
                      z*(lfreq-freq0(iex)+freq(j)) - g(iex,t) - &
                    !  z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
                      dsum)
          end if
        end do
      end do
      do k = 1, nmodes
        aspol(j,1,1) = aspol(j,1,1) + &
                     (osc(iex,1)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,1)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        aspol(j,1,2) = aspol(j,1,2) + &
                     (osc(iex,1)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,1)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        aspol(j,1,3) = aspol(j,1,3) + &
                     (osc(iex,1)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,1)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        aspol(j,2,1) = aspol(j,2,1) + &
                     (osc(iex,2)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,2)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        aspol(j,2,2) = aspol(j,2,2) + &
                     (osc(iex,2)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,2)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        aspol(j,2,3) = aspol(j,2,3) + &
                     (osc(iex,2)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,2)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        aspol(j,3,1) = aspol(j,3,1) + &
                     (osc(iex,3)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,3)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        aspol(j,3,2) = aspol(j,3,2) + &
                     (osc(iex,3)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,3)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        aspol(j,3,3) = aspol(j,3,3) + &
                     (osc(iex,3)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,3)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
      end do
    end do
  end do

end subroutine ASRRSHTFund
