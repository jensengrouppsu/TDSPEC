subroutine RHRSIRCoupling(x, w, freq, delta, osc, tdipder, hpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htir)

  use constants
 
  implicit none

  integer :: nmodes, nexci, tstep
  integer :: j, iex, t, i, r, n, k

  complex(kindr), intent(inout) :: hpol(nmodes,3,3,3)
  complex(kindr) :: rsum1, dsum, z
  complex(kindr) :: Ls_htir(nexci,nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3),tdipder(nexci,nmodes,3)

  complex(kindr) :: zzero=(0d0,0d0)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate A' term hyperpolarizability for RHRS spectra (fundamentals)
!
! input  : x - abcissa points for numerical integration 
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          osc - transition dipole moment components
!          tdipder - ground state dipole moment derivatives
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!
! in/out : hpol - RHRS hyperpolarizability tensor
!
! =============================================================================

! Code for the coupling of IR active modes to hyper-Raman scattering, 
! based on the concept outlined in
!
! L. C. T. Shoute, M. Blanchard-Desce, and A. M. Kelley.  
!   J. Phys. Chem. A 109, 10503 (2005).
!
! Note that doing this correctly requires you to account for the sum over normal
! modes.

  do j = 1, nmodes
    do iex=1,nexci
      Ls_htir(1:nexci,1:nmodes) = zzero
      do t=1,tstep
        dsum = zzero
        do i=1,nmodes
          dsum = dsum + &
                 (delta(iex,i)**2/two)*(one-exp(-z*freq(i)*x(t)))
        end do
        ! Calculate the contribution from the excited normal mode.
        ! If the current mode is raised, the integral is a modified
        ! <0|0(t)> overlap.  If a different mode is raised, use the
        ! <1|0(t)> overlap (rsum1).
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
            ! Raise the initial state of the current mode,
            ! so that I+1 = F = 1, while leaving every other
            ! mode in their ground states.
            Ls_htir(iex,n) = Ls_htir(iex,n) + z*w(t)* &
                     (one-delta(iex,n)**2+ &
                     delta(iex,n)**2*cos(freq(n)*x(t)))*exp(x(t)* &
                     z*(lfreq-freq0(iex)-freq(n)) - g(iex,t) - dsum)
          else
            ! Raise the initial state of a different mode,
            ! so I+1 = 1 and F = 0.  Here, there's a product
            ! accounting for multiple modes in an excited state.
            Ls_htir(iex,n) = Ls_htir(iex,n) + z*w(t)*rsum1*expsum(n)* &
                     exp(x(t)* &
                     z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
          end if
        end do
      end do

      ! Add the A term and the A' term to get the total 
      ! hyperpolarizability.

      do k = 1, nmodes

        hpol(j,1,1,1) = hpol(j,1,1,1) + osc(iex,1)* &
          ((osc(iex,1)*tdipder(iex,k,1) + &
          osc(iex,1)*tdipder(iex,k,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,1,1,2) = hpol(j,1,1,2) + osc(iex,1)* &
          ((osc(iex,1)*tdipder(iex,k,2) + &
          tdipder(iex,k,1)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,1,1,3) = hpol(j,1,1,3) + osc(iex,1)* &
          ((osc(iex,1)*tdipder(iex,k,3) + &
          tdipder(iex,k,1)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si

        hpol(j,1,2,1) = hpol(j,1,2,1) + osc(iex,1)* &
          ((osc(iex,2)*tdipder(iex,k,1) + &
          tdipder(iex,k,2)*osc(iex,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,1,2,2) = hpol(j,1,2,2) + osc(iex,1)* &
          ((osc(iex,2)*tdipder(iex,k,2) + &
          tdipder(iex,k,2)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,1,2,3) = hpol(j,1,2,3) + osc(iex,1)* &
          ((osc(iex,2)*tdipder(iex,k,3) + &
          tdipder(iex,k,2)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si

        hpol(j,1,3,1) = hpol(j,1,3,1) + osc(iex,1)* &
          ((osc(iex,3)*tdipder(iex,k,1) + &
          tdipder(iex,k,3)*osc(iex,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,1,3,2) = hpol(j,1,3,2) + osc(iex,1)* &
          ((osc(iex,3)*tdipder(iex,k,2) + &
          tdipder(iex,k,3)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,1,3,3) = hpol(j,1,3,3) + osc(iex,1)* &
          ((osc(iex,3)*tdipder(iex,k,3) + &
          tdipder(iex,k,3)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si

        hpol(j,2,1,1) = hpol(j,2,1,1) + osc(iex,2)* &
          ((osc(iex,1)*tdipder(iex,k,1) + &
          tdipder(iex,k,1)*osc(iex,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,2,1,2) = hpol(j,2,1,2) + osc(iex,2)* &
          ((osc(iex,1)*tdipder(iex,k,2) + &
          tdipder(iex,k,1)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,2,1,3) = hpol(j,2,1,3) + osc(iex,2)* &
          ((osc(iex,1)*tdipder(iex,k,3) + &
          tdipder(iex,k,1)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si

        hpol(j,2,2,1) = hpol(j,2,2,1) + osc(iex,2)* &
          ((osc(iex,2)*tdipder(iex,k,1) + &
          tdipder(iex,k,2)*osc(iex,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,2,2,2) = hpol(j,2,2,2) + osc(iex,2)* &
          ((osc(iex,2)*tdipder(iex,k,2) + &
          tdipder(iex,k,2)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,2,2,3) = hpol(j,2,2,3) + osc(iex,2)* &
          ((osc(iex,2)*tdipder(iex,k,3) + &
          tdipder(iex,k,2)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si

        hpol(j,2,3,1) = hpol(j,2,3,1) + osc(iex,2)* &
          ((osc(iex,3)*tdipder(iex,k,1) + &
          tdipder(iex,k,3)*osc(iex,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,2,3,2) = hpol(j,2,3,2) + osc(iex,2)* &
          ((osc(iex,3)*tdipder(iex,k,2) + &
          tdipder(iex,k,3)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,2,3,3) = hpol(j,2,3,3) + osc(iex,2)* &
          ((osc(iex,3)*tdipder(iex,k,3) + &
          tdipder(iex,k,3)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si

        hpol(j,3,1,1) = hpol(j,3,1,1) + osc(iex,3)* &
          ((osc(iex,1)*tdipder(iex,k,1) + &
          tdipder(iex,k,1)*osc(iex,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,3,1,2) = hpol(j,3,1,2) + osc(iex,3)* &
          ((osc(iex,1)*tdipder(iex,k,2) + &
          tdipder(iex,k,1)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,3,1,3) = hpol(j,3,1,3) + osc(iex,3)* &
          ((osc(iex,1)*tdipder(iex,k,3) + &
          tdipder(iex,k,1)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si

        hpol(j,3,2,1) = hpol(j,3,2,1) + osc(iex,3)* &
          ((osc(iex,2)*tdipder(iex,k,1) + &
          tdipder(iex,k,2)*osc(iex,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,3,2,2) = hpol(j,3,2,2) + osc(iex,3)* &
          ((osc(iex,2)*tdipder(iex,k,2) + &
          tdipder(iex,k,2)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,3,2,3) = hpol(j,3,2,3) + osc(iex,3)* &
          ((osc(iex,2)*tdipder(iex,k,3) + &
          tdipder(iex,k,2)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si

        hpol(j,3,3,1) = hpol(j,3,3,1) + osc(iex,3)* &
          ((osc(iex,3)*tdipder(iex,k,1) + &
          tdipder(iex,k,3)*osc(iex,1))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,3,3,2) = hpol(j,3,3,2) + osc(iex,3)* &
          ((osc(iex,3)*tdipder(iex,k,2) + &
          tdipder(iex,k,3)*osc(iex,2))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
        hpol(j,3,3,3) = hpol(j,3,3,3) + osc(iex,3)* &
          ((osc(iex,3)*tdipder(iex,k,3) + &
          tdipder(iex,k,3)*osc(iex,3))/ &
          (freq(k)-f12*lfreq))* &
          Ls_htir(iex,k)*hpol2si
      end do
    end do
  end do

end subroutine RHRSIRCoupling
