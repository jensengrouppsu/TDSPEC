subroutine ASRHRSFCFund(latermoff, x, w, freq, delta, osc, stpm, ashpol, &
                   freq0, lfreq, g, nexci, nmodes, tstep)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep

  complex(kindr), intent(inout) :: ashpol(nmodes,3,3,3)
  complex(kindr) :: Ls, rsum, dsum, z
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3), stpm(nexci,6)

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: latermoff

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate A term hyperpolarizability for anti-Stokes RHRS spectra 
!          (fundamentals)
!
! input  : latermoff - determines if only B terms are calculated
!          x - abcissa points for numerical integration
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          osc - transition dipole moment components
!          stpm - two-photon transition moment components
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!
! in/out : ashpol - anti-Stokes RHRS hyperpolarizability tensor
!
! =============================================================================

  ! This is needed because some compilers don't default arrays to being zero.
  ashpol = zzero

  ! Integration is carried out for this part very similarly to how 
  ! it was done for normal Raman scattering.

  if (latermoff.eqv..false.) then
    do j = 1, nmodes
      do iex=1,nexci
        Ls = zzero
        do t=1,tstep
          dsum = zzero
          do i=1,nmodes
            dsum = dsum + &
                   (delta(iex,i)**2/two)*(one-exp(-z*freq(i)*x(t)))
          end do
          rsum = delta(iex,j)/(sqrt(two))* &
                 (one-exp(-z*freq(j)*x(t)))
          Ls = Ls + z*w(t)*rsum*exp(x(t)* &
               z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
               dsum)
        end do

        ! The lineshape function is still in units of cm in this 
        ! code.  One-photon and two-photon transition dipole moments
        ! are most likely in atomic units, where the two-photon 
        ! transition dipole moments are actually in units of 
        ! a.u.^2 / Hartree.  Without performing conversions, the 
        ! units of the hyperpolarizability are: a.u.^3*cm/Hartree 
        ! or alternatively, a.u.^3/Hartree*cm^{-1}.

        ! Again, this part of the code utilizes that atomic units are
        ! related to SI units for dipoles as: 
        ! 1 Cm = 1.1794744E+29 a.u.  Also, energy units are easily 
        ! converted to SI units by using 1 cm^{-1} = 1.9864456E-23 J 
        ! and 1 Hartree = 219474.63068 cm^{-1}.  Performing those 
        ! conversions allows hyperpolarizability to be stated in SI 
        ! units: Cm^3/V^2 (since 1 V = 1 J/C). 

        ! To alleviate future confusion on my part:
        ! S_{xx} = stpm(#,1)
        ! S_{yy} = stpm(#,2)
        ! S_{zz} = stpm(#,3)
        ! S_{xy} = stpm(#,4) = S_{yx}
        ! S_{xz} = stpm(#,5) = S_{zx}
        ! S_{yz} = stpm(#,6) = S_{zy}

        ! Constants used:
        ! dipconv = 1 C*m / 1.1794744E+29 a.u. , 1 a.u. = 1 e*a_0
        ! hartree2cm = 219474.63068 cm^{-1} / 1 Hartree
        ! hartree2joules = 4.359743357E-18 J / 1 Hartree
    
        ! hpol2si = (1 C*m/1.1794744E+29 a.u.)^3*
        !           (219474.63068 cm^{-1}/1 Hartree)/
        !           (4.359743357E-18 J/1 Hartree)^2


        ashpol(j,1,1,1) = ashpol(j,1,1,1) + osc(iex,1)*stpm(iex,1)*Ls* &
          hpol2si
        ashpol(j,1,1,2) = ashpol(j,1,1,2) + osc(iex,1)*stpm(iex,4)*Ls* &
          hpol2si
        ashpol(j,1,1,3) = ashpol(j,1,1,3) + osc(iex,1)*stpm(iex,5)*Ls* &
          hpol2si

        ashpol(j,1,2,1) = ashpol(j,1,2,1) + osc(iex,1)*stpm(iex,4)*Ls* &
          hpol2si
        ashpol(j,1,2,2) = ashpol(j,1,2,2) + osc(iex,1)*stpm(iex,2)*Ls* &
          hpol2si
        ashpol(j,1,2,3) = ashpol(j,1,2,3) + osc(iex,1)*stpm(iex,6)*Ls* &
          hpol2si

        ashpol(j,1,3,1) = ashpol(j,1,3,1) + osc(iex,1)*stpm(iex,5)*Ls* &
          hpol2si
        ashpol(j,1,3,2) = ashpol(j,1,3,2) + osc(iex,1)*stpm(iex,6)*Ls* &
          hpol2si
        ashpol(j,1,3,3) = ashpol(j,1,3,3) + osc(iex,1)*stpm(iex,3)*Ls* &
          hpol2si

        ashpol(j,2,1,1) = ashpol(j,2,1,1) + osc(iex,2)*stpm(iex,1)*Ls* &
          hpol2si
        ashpol(j,2,1,2) = ashpol(j,2,1,2) + osc(iex,2)*stpm(iex,4)*Ls* &
          hpol2si
        ashpol(j,2,1,3) = ashpol(j,2,1,3) + osc(iex,2)*stpm(iex,5)*Ls* &
          hpol2si

        ashpol(j,2,2,1) = ashpol(j,2,2,1) + osc(iex,2)*stpm(iex,4)*Ls* &
          hpol2si
        ashpol(j,2,2,2) = ashpol(j,2,2,2) + osc(iex,2)*stpm(iex,2)*Ls* &
          hpol2si
        ashpol(j,2,2,3) = ashpol(j,2,2,3) + osc(iex,2)*stpm(iex,6)*Ls* &
          hpol2si

        ashpol(j,2,3,1) = ashpol(j,2,3,1) + osc(iex,2)*stpm(iex,5)*Ls* &
          hpol2si
        ashpol(j,2,3,2) = ashpol(j,2,3,2) + osc(iex,2)*stpm(iex,6)*Ls* &
          hpol2si
        ashpol(j,2,3,3) = ashpol(j,2,3,3) + osc(iex,2)*stpm(iex,3)*Ls* &
          hpol2si

        ashpol(j,3,1,1) = ashpol(j,3,1,1) + osc(iex,3)*stpm(iex,1)*Ls* &
          hpol2si
        ashpol(j,3,1,2) = ashpol(j,3,1,2) + osc(iex,3)*stpm(iex,4)*Ls* &
          hpol2si
        ashpol(j,3,1,3) = ashpol(j,3,1,3) + osc(iex,3)*stpm(iex,5)*Ls* &
          hpol2si

        ashpol(j,3,2,1) = ashpol(j,3,2,1) + osc(iex,3)*stpm(iex,4)*Ls* &
          hpol2si
        ashpol(j,3,2,2) = ashpol(j,3,2,2) + osc(iex,3)*stpm(iex,2)*Ls* &
          hpol2si
        ashpol(j,3,2,3) = ashpol(j,3,2,3) + osc(iex,3)*stpm(iex,6)*Ls* &
          hpol2si

        ashpol(j,3,3,1) = ashpol(j,3,3,1) + osc(iex,3)*stpm(iex,5)*Ls* &
          hpol2si
        ashpol(j,3,3,2) = ashpol(j,3,3,2) + osc(iex,3)*stpm(iex,6)*Ls* &
          hpol2si
        ashpol(j,3,3,3) = ashpol(j,3,3,3) + osc(iex,3)*stpm(iex,3)*Ls* &
          hpol2si
      end do
    end do
  end if

end subroutine ASRHRSFCFund

subroutine ASRHRSB2Fund(x, w, freq, delta, osc, stpmder, ashpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb2, Ls_htb3)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep, r, k, n

  complex(kindr), intent(inout) :: ashpol(nmodes,3,3,3)
  complex(kindr) :: rsum, rsum1, dsum, z
  complex(kindr) :: Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3),stpmder(nexci,nmodes,6)

  complex(kindr) :: zzero=(0d0,0d0)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B2 term hyperpolarizability for anti-Stokes RHRS spectra 
!          (fundamentals)
!
! input  : x - abcissa points for numerical integration 
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          osc - transition moment components
!          stpmder - derivative of two-photon transition moment components
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!
! in/out : ashpol - anti-Stokes RHRS hyperpolarizability tensor
!
! =============================================================================

! Evaluates the the Herzberg-Teller terms for each normal mode. 
! Equations are based on the conversion from the vibronic theory 
! to the time-dependent formalism.
!
! Y. C. Chung and L. D. Ziegler.  J. Chem. Phys. 88, 7287 (1988).
!
! This part evaluates B2 of the Herzberg-Teller terms,  which involves 
! derivatives of the two-photon transition moment between the initial 
! and final electronic states.  For anti-Stokes spectra, the lineshape
! function of the B2 term is the same as the B1 term lineshape function
! for Stokes spectra.

  do j = 1, nmodes
    do iex=1,nexci
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
          ! Integrals of the form: <F|(I-1)(t)>
          if (n == j) then
            ! Lower the initial state of the current mode,
            ! so that F = I-1 = 0
            Ls_htb2(iex,n) = Ls_htb2(iex,n) + z*w(t)*exp(x(t)* &
                      z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
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
                      z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
                      dsum)
          else
            ! Raise the final state of a different mode,
            ! so that I = 0 and F = 1.
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      z*w(t)*expsum(n)*rsum1*exp(x(t)* &
                      z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
                      dsum)
          end if
        end do
      end do

      ! Add the A, B1, and B2 terms together.
      ! For future reference:
      !
      ! stpmder(#,#,1) = dS_{xx}/dQ
      ! stpmder(#,#,2) = dS_{yy}/dQ
      ! stpmder(#,#,3) = dS_{zz}/dQ
      ! stpmder(#,#,4) = dS_{xy}/dQ = dS_{yx}/dQ
      ! stpmder(#,#,5) = dS_{xz}/dQ = dS_{zx}/dQ
      ! stpmder(#,#,6) = dS_{yz}/dQ = dS_{zy}/dQ
      !
      ! For the single-residue quadratic response calculations, 
      ! it turns out that S is a symmetric matrix, and so is the
      ! matrix composed of derivatives of S.

      do k = 1, nmodes

        ashpol(j,1,1,1) = ashpol(j,1,1,1) + &
          osc(iex,1)*stpmder(iex,k,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,1,1,2) = ashpol(j,1,1,2) + &
          osc(iex,1)*stpmder(iex,k,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,1,1,3) = ashpol(j,1,1,3) + &
          osc(iex,1)*stpmder(iex,k,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        ashpol(j,1,2,1) = ashpol(j,1,2,1) + &
          osc(iex,1)*stpmder(iex,k,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,1,2,2) = ashpol(j,1,2,2) + &
          osc(iex,1)*stpmder(iex,k,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,1,2,3) = ashpol(j,1,2,3) + &
          osc(iex,1)*stpmder(iex,k,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        ashpol(j,1,3,1) = ashpol(j,1,3,1) + &
          osc(iex,1)*stpmder(iex,k,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,1,3,2) = ashpol(j,1,3,2) + &
          osc(iex,1)*stpmder(iex,k,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,1,3,3) = ashpol(j,1,3,3) + &
          osc(iex,1)*stpmder(iex,k,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        ashpol(j,2,1,1) = ashpol(j,2,1,1) + &
          osc(iex,2)*stpmder(iex,k,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,2,1,2) = ashpol(j,2,1,2) + &
          osc(iex,2)*stpmder(iex,k,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,2,1,3) = ashpol(j,2,1,3) + &
          osc(iex,2)*stpmder(iex,k,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        ashpol(j,2,2,1) = ashpol(j,2,2,1) + &
          osc(iex,2)*stpmder(iex,k,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,2,2,2) = ashpol(j,2,2,2) + &
          osc(iex,2)*stpmder(iex,k,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,2,2,3) = ashpol(j,2,2,3) + &
          osc(iex,2)*stpmder(iex,k,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        ashpol(j,2,3,1) = ashpol(j,2,3,1) + &
          osc(iex,2)*stpmder(iex,k,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,2,3,2) = ashpol(j,2,3,2) + &
          osc(iex,2)*stpmder(iex,k,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,2,3,3) = ashpol(j,2,3,3) + &
          osc(iex,2)*stpmder(iex,k,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        ashpol(j,3,1,1) = ashpol(j,3,1,1) + &
          osc(iex,3)*stpmder(iex,k,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,3,1,2) = ashpol(j,3,1,2) + &
          osc(iex,3)*stpmder(iex,k,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,3,1,3) = ashpol(j,3,1,3) + &
          osc(iex,3)*stpmder(iex,k,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        ashpol(j,3,2,1) = ashpol(j,3,2,1) + &
          osc(iex,3)*stpmder(iex,k,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,3,2,2) = ashpol(j,3,2,2) + &
          osc(iex,3)*stpmder(iex,k,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,3,2,3) = ashpol(j,3,2,3) + &
          osc(iex,3)*stpmder(iex,k,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        ashpol(j,3,3,1) = ashpol(j,3,3,1) + &
          osc(iex,3)*stpmder(iex,k,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,3,3,2) = ashpol(j,3,3,2) + &
          osc(iex,3)*stpmder(iex,k,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        ashpol(j,3,3,3) = ashpol(j,3,3,3) + &
          osc(iex,3)*stpmder(iex,k,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
      end do
    end do
  end do

end subroutine ASRHRSB2Fund

subroutine ASRHRSB1Fund(x, w, freq, delta, stpm, tdipder, ashpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep, r, k, n

  complex(kindr), intent(inout) :: ashpol(nmodes,3,3,3)
  complex(kindr) :: rsum, rsum1, dsum, z
  complex(kindr) :: Ls_htb(nexci,nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: stpm(nexci,6),tdipder(nexci,nmodes,3)

  complex(kindr) :: zzero=(0d0,0d0)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B1 term hyperpolarizability for anti-Stokes RHRS spectra 
!          (fundamentals)
!
! input  : x - abcissa points for numerical integration 
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          stpm - two-photon transition moment components
!          tdipder - derivative of transition dipole moment components
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!
! in/out : ashpol - anti-Stokes RHRS hyperpolarizability tensor
!
! =============================================================================

! This part evaluates the remaining part (B2) of the Herzberg-Teller 
! expansion for resonance hyper-Raman scattering.  It is only used if
! the user requests this part of the code, because the program requires 
! derivatives of the two-photon moments along each normal mode to 
! evaluate these contributions. 

  do j = 1, nmodes
    do iex=1,nexci
      Ls_htb(1:nexci,1:nmodes) = zzero
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
          ! Integrals of the form: <F|(I+1)(t)>
          if (n == j) then
            ! Raise the final state of the current mode,
            ! so that I = F+1 = 1, while leaving every other
            ! mode in their ground states. 
            ! In anti-Stokes, the extra frequency from the
            ! lineshape function cancels the frequency of the
            ! initial state.
            Ls_htb(iex,n) = Ls_htb(iex,n) + z*w(t)* &
                     (one-delta(iex,n)**2+ &
                     delta(iex,n)**2*cos(freq(n)*x(t)))*exp(x(t)* &
                     z*(lfreq-freq0(iex)) - conjg(g(iex,t)) - dsum)
          else
            ! Raise the final state of a different mode,
            ! so F+1 = 1 and I = 0.  Here, there's a product
            ! accounting for multiple modes in an excited state.
            Ls_htb(iex,n) = Ls_htb(iex,n) + z*w(t)*rsum1*expsum(n)* &
                     exp(x(t)* &
                     z*(lfreq-freq0(iex)+freq(j)) - conjg(g(iex,t)) - &
                     dsum)
          end if
        end do
      end do

      ! Add the A term to the B1 term to get the total 
      ! hyperpolarizability.

      do k = 1, nmodes

        ashpol(j,1,1,1) = ashpol(j,1,1,1) + &
          tdipder(iex,k,1)*stpm(iex,1)*Ls_htb(iex,k)*hpol2si
        ashpol(j,1,1,2) = ashpol(j,1,1,2) + &
          tdipder(iex,k,1)*stpm(iex,4)*Ls_htb(iex,k)*hpol2si
        ashpol(j,1,1,3) = ashpol(j,1,1,3) + &
          tdipder(iex,k,1)*stpm(iex,5)*Ls_htb(iex,k)*hpol2si

        ashpol(j,1,2,1) = ashpol(j,1,2,1) + &
          tdipder(iex,k,1)*stpm(iex,4)*Ls_htb(iex,k)*hpol2si
        ashpol(j,1,2,2) = ashpol(j,1,2,2) + &
          tdipder(iex,k,1)*stpm(iex,2)*Ls_htb(iex,k)*hpol2si
        ashpol(j,1,2,3) = ashpol(j,1,2,3) + &
          tdipder(iex,k,1)*stpm(iex,6)*Ls_htb(iex,k)*hpol2si

        ashpol(j,1,3,1) = ashpol(j,1,3,1) + &
          tdipder(iex,k,1)*stpm(iex,5)*Ls_htb(iex,k)*hpol2si
        ashpol(j,1,3,2) = ashpol(j,1,3,2) + &
          tdipder(iex,k,1)*stpm(iex,6)*Ls_htb(iex,k)*hpol2si
        ashpol(j,1,3,3) = ashpol(j,1,3,3) + &
          tdipder(iex,k,1)*stpm(iex,3)*Ls_htb(iex,k)*hpol2si

        ashpol(j,2,1,1) = ashpol(j,2,1,1) + &
          tdipder(iex,k,2)*stpm(iex,1)*Ls_htb(iex,k)*hpol2si
        ashpol(j,2,1,2) = ashpol(j,2,1,2) + &
          tdipder(iex,k,2)*stpm(iex,4)*Ls_htb(iex,k)*hpol2si
        ashpol(j,2,1,3) = ashpol(j,2,1,3) + &
          tdipder(iex,k,2)*stpm(iex,5)*Ls_htb(iex,k)*hpol2si

        ashpol(j,2,2,1) = ashpol(j,2,2,1) + &
          tdipder(iex,k,2)*stpm(iex,4)*Ls_htb(iex,k)*hpol2si
        ashpol(j,2,2,2) = ashpol(j,2,2,2) + &
          tdipder(iex,k,2)*stpm(iex,2)*Ls_htb(iex,k)*hpol2si
        ashpol(j,2,2,3) = ashpol(j,2,2,3) + &
          tdipder(iex,k,2)*stpm(iex,6)*Ls_htb(iex,k)*hpol2si

        ashpol(j,2,3,1) = ashpol(j,2,3,1) + &
          tdipder(iex,k,2)*stpm(iex,5)*Ls_htb(iex,k)*hpol2si
        ashpol(j,2,3,2) = ashpol(j,2,3,2) + &
          tdipder(iex,k,2)*stpm(iex,6)*Ls_htb(iex,k)*hpol2si
        ashpol(j,2,3,3) = ashpol(j,2,3,3) + &
          tdipder(iex,k,2)*stpm(iex,3)*Ls_htb(iex,k)*hpol2si

        ashpol(j,3,1,1) = ashpol(j,3,1,1) + &
          tdipder(iex,k,3)*stpm(iex,1)*Ls_htb(iex,k)*hpol2si
        ashpol(j,3,1,2) = ashpol(j,3,1,2) + &
          tdipder(iex,k,3)*stpm(iex,4)*Ls_htb(iex,k)*hpol2si
        ashpol(j,3,1,3) = ashpol(j,3,1,3) + &
          tdipder(iex,k,3)*stpm(iex,5)*Ls_htb(iex,k)*hpol2si

        ashpol(j,3,2,1) = ashpol(j,3,2,1) + &
          tdipder(iex,k,3)*stpm(iex,4)*Ls_htb(iex,k)*hpol2si
        ashpol(j,3,2,2) = ashpol(j,3,2,2) + &
          tdipder(iex,k,3)*stpm(iex,2)*Ls_htb(iex,k)*hpol2si
        ashpol(j,3,2,3) = ashpol(j,3,2,3) + &
          tdipder(iex,k,3)*stpm(iex,6)*Ls_htb(iex,k)*hpol2si

        ashpol(j,3,3,1) = ashpol(j,3,3,1) + &
          tdipder(iex,k,3)*stpm(iex,5)*Ls_htb(iex,k)*hpol2si
        ashpol(j,3,3,2) = ashpol(j,3,3,2) + &
          tdipder(iex,k,3)*stpm(iex,6)*Ls_htb(iex,k)*hpol2si
        ashpol(j,3,3,3) = ashpol(j,3,3,3) + &
          tdipder(iex,k,3)*stpm(iex,3)*Ls_htb(iex,k)*hpol2si
      end do
    end do
  end do

end subroutine ASRHRSB1Fund
