subroutine SFGFCFund(latermoff, x, w, freq, delta, osc, gdipder, hpol, &
                   freq0, lfreq, g, nexci, nmodes, tstep, IRUV)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep

  complex(kindr), intent(inout) :: hpol(nmodes,3,3,3)
  complex(kindr) :: Ls, rsum, dsum, z
  complex(kindr),intent(in) :: g(nexci,tstep)
  
  real(kindr) :: freq(nmodes), delta(nexci,nmodes), lfreqadj(nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3), gdipder(nexci,nmodes,3)

  complex(kindr) :: zzero=(0d0,0d0)
  
  logical :: latermoff, IRUV

  z = csqrt((-1,0))
  if (IRUV.eqv..true.) then
    lfreqadj = freq
  else 
    lfreqadj = 0 
  end if

! =============================================================================
! purpose: Evaluate A term hyperpolarizability for SFG spectra (fundamentals)
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
!          IRUV - Whether or not IR-UV SFG is being calculated
!          lfreqadj - Adjust the laser freqency for IR-UV SFG
!
! in/out : hpol - SFG hyperpolarizability tensor
!
! =============================================================================

  ! This is needed because some compilers don't default arrays to being zero.
  hpol = zzero

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
               z*(lfreq+lfreqadj(j)-freq0(iex)) - g(iex,t) - dsum)
        end do

        ! Units are a little strange at this point: (e*a0)^3 / cm^{-1}^2
        ! These are converted to "alternative SI units" in the 
        ! averaging routine according to:
        ! Shelton, D. P. and Rice, J. E.  Chem. Rev. 94, 3 (1994).

        hpol(j,1,1,1) = hpol(j,1,1,1) + osc(iex,1)*osc(iex,1)*gdipder(iex,j,1)*Ls 
        hpol(j,1,1,2) = hpol(j,1,1,2) + osc(iex,1)*osc(iex,1)*gdipder(iex,j,2)*Ls 
        hpol(j,1,1,3) = hpol(j,1,1,3) + osc(iex,1)*osc(iex,1)*gdipder(iex,j,3)*Ls 

        hpol(j,1,2,1) = hpol(j,1,2,1) + osc(iex,1)*osc(iex,2)*gdipder(iex,j,1)*Ls 
        hpol(j,1,2,2) = hpol(j,1,2,2) + osc(iex,1)*osc(iex,2)*gdipder(iex,j,2)*Ls 
        hpol(j,1,2,3) = hpol(j,1,2,3) + osc(iex,1)*osc(iex,2)*gdipder(iex,j,3)*Ls 

        hpol(j,1,3,1) = hpol(j,1,3,1) + osc(iex,1)*osc(iex,3)*gdipder(iex,j,1)*Ls 
        hpol(j,1,3,2) = hpol(j,1,3,2) + osc(iex,1)*osc(iex,3)*gdipder(iex,j,2)*Ls 
        hpol(j,1,3,3) = hpol(j,1,3,3) + osc(iex,1)*osc(iex,3)*gdipder(iex,j,3)*Ls 

        hpol(j,2,1,1) = hpol(j,2,1,1) + osc(iex,2)*osc(iex,1)*gdipder(iex,j,1)*Ls 
        hpol(j,2,1,2) = hpol(j,2,1,2) + osc(iex,2)*osc(iex,1)*gdipder(iex,j,2)*Ls 
        hpol(j,2,1,3) = hpol(j,2,1,3) + osc(iex,2)*osc(iex,1)*gdipder(iex,j,3)*Ls 

        hpol(j,2,2,1) = hpol(j,2,2,1) + osc(iex,2)*osc(iex,2)*gdipder(iex,j,1)*Ls 
        hpol(j,2,2,2) = hpol(j,2,2,2) + osc(iex,2)*osc(iex,2)*gdipder(iex,j,2)*Ls 
        hpol(j,2,2,3) = hpol(j,2,2,3) + osc(iex,2)*osc(iex,2)*gdipder(iex,j,3)*Ls 

        hpol(j,2,3,1) = hpol(j,2,3,1) + osc(iex,2)*osc(iex,3)*gdipder(iex,j,1)*Ls 
        hpol(j,2,3,2) = hpol(j,2,3,2) + osc(iex,2)*osc(iex,3)*gdipder(iex,j,2)*Ls 
        hpol(j,2,3,3) = hpol(j,2,3,3) + osc(iex,2)*osc(iex,3)*gdipder(iex,j,3)*Ls 

        hpol(j,3,1,1) = hpol(j,3,1,1) + osc(iex,3)*osc(iex,1)*gdipder(iex,j,1)*Ls 
        hpol(j,3,1,2) = hpol(j,3,1,2) + osc(iex,3)*osc(iex,1)*gdipder(iex,j,2)*Ls 
        hpol(j,3,1,3) = hpol(j,3,1,3) + osc(iex,3)*osc(iex,1)*gdipder(iex,j,3)*Ls 

        hpol(j,3,2,1) = hpol(j,3,2,1) + osc(iex,3)*osc(iex,2)*gdipder(iex,j,1)*Ls 
        hpol(j,3,2,2) = hpol(j,3,2,2) + osc(iex,3)*osc(iex,2)*gdipder(iex,j,2)*Ls 
        hpol(j,3,2,3) = hpol(j,3,2,3) + osc(iex,3)*osc(iex,2)*gdipder(iex,j,3)*Ls 

        hpol(j,3,3,1) = hpol(j,3,3,1) + osc(iex,3)*osc(iex,3)*gdipder(iex,j,1)*Ls 
        hpol(j,3,3,2) = hpol(j,3,3,2) + osc(iex,3)*osc(iex,3)*gdipder(iex,j,2)*Ls 
        hpol(j,3,3,3) = hpol(j,3,3,3) + osc(iex,3)*osc(iex,3)*gdipder(iex,j,3)*Ls 

      end do
    end do
  end if

end subroutine SFGFCFund

subroutine SFGHTFund(x, w, freq, delta, osc, tdipder, gdipder,hpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb, Ls_htb2, Ls_htb3, IRUV)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep, r, k, n

  complex(kindr), intent(inout) :: hpol(nmodes,3,3,3)
  complex(kindr) :: ht_hpol(nmodes,nmodes,3,3,3)
  complex(kindr) :: rsum, rsum1, dsum, z
  complex(kindr) :: Ls_htb(nexci,nmodes)
  complex(kindr) :: Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq, lfreqadj(nmodes)
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3),tdipder(nexci,nmodes,3)
  real(kindr) :: gdipder(nexci,nmodes,3)
  real(kindr) :: htlfreqdiff(nexci)
  real(kindr) :: ht_rcrs(nmodes,nmodes)
  real(kindr) :: ht_rcrs_sum(nmodes)
  real(kindr) :: htexcist

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: IRUV

  z = csqrt((-1,0))
  if (IRUV.eqv..true.) then 
    lfreqadj = freq
  else
    lfreqadj = 0
  end if

! =============================================================================
! purpose: Evaluate B term hyperpolarizability for SFG spectra (fundamentals)
!
! input  : x - abcissa points for numerical integration 
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          tdipder - derivative of transition dipole moment components
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!          IRUV - Whether or not IR-UV SFG is being calculated
!          lfreqadj - Adjust the laser freqency for IR-UV SFG
!
! in/out : hpol - SFG hyperpolarizability tensor
!
! =============================================================================

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
          ! Integrals of the form: <F|(I+1)(t)>
          if (n == j) then
            ! Raise the initial state of the current mode,
            ! so that I+1 = F = 1.  A factor of 1-delta^2 is
            ! applied for the 1-1 overlap of a normal mode.
            Ls_htb(iex,n) = Ls_htb(iex,n) + &
                     z*(one-delta(iex,n)**2+ &
                     delta(iex,n)**2*cos(freq(n)*x(t)))* &
                     w(t)*exp(x(t)* &
                     z*(lfreq+lfreqadj(n)-freq0(iex)-freq(n)) - g(iex,t) - dsum)
          else
            ! Raise the initial state of a different mode,
            ! so I+1 = 1 and F = 0.  Here, there's a product
            ! accounting for multiple modes in an excited state.
             Ls_htb(iex,n) = Ls_htb(iex,n) + z*w(t)*rsum1*expsum(n)* &
                      exp(x(t)* &
                      z*(lfreq+lfreqadj(n)-freq0(iex)) - g(iex,t) - dsum)
          end if
          ! Integrals of the form: <(F-1)|I(t)>
          if (n == j) then
            ! Lower the final state of the current mode,
            ! so that F-1 = I = 0
            Ls_htb2(iex,n) = Ls_htb2(iex,n) + z*w(t)*exp(x(t)* &
                      z*(lfreq+lfreqadj(n)-freq0(iex)) - g(iex,t) - dsum)
          else
            ! Can't lower a harmonic oscillator past its ground
            ! state.
            Ls_htb2(iex,n) = zzero
          end if
          ! Integrals of the form: <(F+1)|I(t)> 
          if (n == j) then
            ! Raise the final state of the current mode,
            ! so that I = 0 and F+1 = 2
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      z*w(t)*root2*rsum*exp(x(t)* &
                      z*(lfreq+lfreqadj(n)-freq0(iex)) - g(iex,t) - dsum)
          else
            ! Raise the final state of a different mode,
            ! so that I = 0 and F = 1.
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      z*w(t)*expsum(n)*rsum1*exp(x(t)* &
                      z*(lfreq+lfreqadj(n)-freq0(iex)) - g(iex,t) - dsum)
          end if
        end do
      end do

      ! Add the A term to the B term to get the total 
      ! hyperpolarizability.

      ! Units are a little strange at this point: (e*a0)^3 / cm^{-1}
      ! These are converted to "alternative SI units" in the 
      ! averaging routine according to:
      ! Shelton, D. P. and Rice, J. E.  Chem. Rev. 94, 3 (1994).

      do k = 1, nmodes

!===============================
        hpol(j,1,1,1) = hpol(j,1,1,1) + &
                        (osc(iex,1)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,1,1,2) = hpol(j,1,1,2) + &
                        (osc(iex,1)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,1,1,3) = hpol(j,1,1,3) + &
                        (osc(iex,1)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
        hpol(j,1,2,1) = hpol(j,1,2,1) + &
                        (osc(iex,1)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,1,2,2) = hpol(j,1,2,2) + &
                        (osc(iex,1)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,1,2,3) = hpol(j,1,2,3) + &
                        (osc(iex,1)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
        hpol(j,1,3,1) = hpol(j,1,3,1) + &
                        (osc(iex,1)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,1,3,2) = hpol(j,1,3,2) + &
                        (osc(iex,1)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,1,3,3) = hpol(j,1,3,3) + &
                        (osc(iex,1)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,1)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
        hpol(j,2,1,1) = hpol(j,2,1,1) + &
                        (osc(iex,2)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,2,1,2) = hpol(j,2,1,2) + &
                        (osc(iex,2)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,2,1,3) = hpol(j,2,1,3) + &
                        (osc(iex,2)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
        hpol(j,2,2,1) = hpol(j,2,2,1) + &
                        (osc(iex,2)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,2,2,2) = hpol(j,2,2,2) + &
                        (osc(iex,2)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,2,2,3) = hpol(j,2,2,3) + &
                        (osc(iex,2)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
        hpol(j,2,3,1) = hpol(j,2,3,1) + &
                        (osc(iex,2)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,2,3,2) = hpol(j,2,3,2) + &
                        (osc(iex,2)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,2,3,3) = hpol(j,2,3,3) + &
                        (osc(iex,2)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,2)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
        hpol(j,3,1,1) = hpol(j,3,1,1) + &
                        (osc(iex,3)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,3,1,2) = hpol(j,3,1,2) + &
                        (osc(iex,3)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,3,1,3) = hpol(j,3,1,3) + &
                        (osc(iex,3)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
        hpol(j,3,2,1) = hpol(j,3,2,1) + &
                        (osc(iex,3)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,3,2,2) = hpol(j,3,2,2) + &
                        (osc(iex,3)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,3,2,3) = hpol(j,3,2,3) + &
                        (osc(iex,3)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
        hpol(j,3,3,1) = hpol(j,3,3,1) + &
                        (osc(iex,3)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,1)

        hpol(j,3,3,2) = hpol(j,3,3,2) + &
                        (osc(iex,3)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,2)

        hpol(j,3,3,3) = hpol(j,3,3,3) + &
                        (osc(iex,3)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                        tdipder(iex,k,3)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                        Ls_htb3(iex,k)))*gdipder(iex,k,3)
!===============================
      end do
    end do
  end do

end subroutine SFGHTFund

subroutine SFGFullCalc(freq, fhpol, nmodes, width, freqstrt, freqend, res, nhpol, lnonRes, shpol)

  use constants
  implicit none

  integer :: i,j, nmodes, freqend,a,b,c

  real(kindr) :: freq(nmodes),width, freqstrt
  complex(kindr) :: nhpol(freqend,3,3,3), z, shpol(3,3,3)
  complex(kindr) :: fhpol(nmodes, 3,3,3)
  logical :: lnonRes

  real(kindr) :: res, IR

! =============================================================================
! purpose: Calculate each hyperpolariazability for each wavelength for a given
!          range without doing a convolution
!
! input  : freq - vibrational frequencies
!          fhpol - SFG hyperpolarizability tensor
!          nmodes - number of normal modes
!          width - lifetime of the IR transition in wavenumbers
!          freqstrt - adjusted starting wavelength
!          freqend - number of points scanned
!          res - resolution or step size in wavenumbers
!           
! in/out : nhpol - SFG hyperpolarizabily tensor for full spectrum
!
! =============================================================================




  z = csqrt((-1,0))    
  do i = 1 , nmodes
    do j = 1 , freqend
      !adjust the wavelength since the IR trasntion is first
      IR = res*j + freqstrt

      nhpol(j,1,1,1) = nhpol(j,1,1,1) + fhpol(i,1,1,1)*(1/((freq(i)-IR-(z*width)))) 
      nhpol(j,1,1,2) = nhpol(j,1,1,2) + fhpol(i,1,1,2)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,1,1,3) = nhpol(j,1,1,3) + fhpol(i,1,1,3)*(1/((freq(i)-IR-(z*width))))

      nhpol(j,1,2,1) = nhpol(j,1,2,1) + fhpol(i,1,2,1)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,1,2,2) = nhpol(j,1,2,2) + fhpol(i,1,2,2)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,1,2,3) = nhpol(j,1,2,3) + fhpol(i,1,2,3)*(1/((freq(i)-IR-(z*width))))

      nhpol(j,1,3,1) = nhpol(j,1,3,1) + fhpol(i,1,3,1)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,1,3,2) = nhpol(j,1,3,2) + fhpol(i,1,3,2)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,1,3,3) = nhpol(j,1,3,3) + fhpol(i,1,3,3)*(1/((freq(i)-IR-(z*width))))

      nhpol(j,2,1,1) = nhpol(j,2,1,1) + fhpol(i,2,1,1)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,2,1,2) = nhpol(j,2,1,2) + fhpol(i,2,1,2)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,2,1,3) = nhpol(j,2,1,3) + fhpol(i,2,1,3)*(1/((freq(i)-IR-(z*width))))

      nhpol(j,2,2,1) = nhpol(j,2,2,1) + fhpol(i,2,2,1)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,2,2,2) = nhpol(j,2,2,2) + fhpol(i,2,2,2)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,2,2,3) = nhpol(j,2,2,3) + fhpol(i,2,2,3)*(1/((freq(i)-IR-(z*width))))

      nhpol(j,2,3,1) = nhpol(j,2,3,1) + fhpol(i,2,3,1)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,2,3,2) = nhpol(j,2,3,2) + fhpol(i,2,3,2)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,2,3,3) = nhpol(j,2,3,3) + fhpol(i,2,3,3)*(1/((freq(i)-IR-(z*width))))

      nhpol(j,3,1,1) = nhpol(j,3,1,1) + fhpol(i,3,1,1)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,3,1,2) = nhpol(j,3,1,2) + fhpol(i,3,1,2)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,3,1,3) = nhpol(j,3,1,3) + fhpol(i,3,1,3)*(1/((freq(i)-IR-(z*width))))

      nhpol(j,3,2,1) = nhpol(j,3,2,1) + fhpol(i,3,2,1)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,3,2,2) = nhpol(j,3,2,2) + fhpol(i,3,2,2)*(1/((freq(i)-IR-(z*width))))
      nhpol(j,3,2,3) = nhpol(j,3,2,3) + fhpol(i,3,2,3)*(1/((freq(i)-IR-(z*width))))
  
      nhpol(j,3,3,1) = nhpol(j,3,3,1) + fhpol(i,3,3,1)*(1/((freq(j)-IR-(z*width))))
      nhpol(j,3,3,2) = nhpol(j,3,3,2) + fhpol(i,3,3,2)*(1/((freq(j)-IR-(z*width))))
      nhpol(j,3,3,3) = nhpol(j,3,3,3) + fhpol(i,3,3,3)*(1/((freq(j)-IR-(z*width))))

      end do
  end do    

  if (lnonRes) then
    shpol = shpol * (1/(hartree2cm**2))
    do j = 1,freqend
      do a=1,3
        do b=1,3
          do c=1,3
            nhpol(j,a,b,c) = nhpol(j,a,b,c) + shpol(a,b,c)
          end do
        end do 
      end do
    end do
  end if
end subroutine SFGFullCalc



subroutine SFGnonRes(nmodes, hpol, width, shpol)



! =============================================================================
! purpose: To add in the non-resonant hyperpolarizability
!
! input  : nmodes - number of normal modes
!          hpol - overall hyperpolarizability
!          width - sets the width of the Lorentzian line function used for
!                  convolution. Here to correct for when the whole
!                  hyperpolarizability is divided by the width in sfgaverage.f90 
!          shpol - static hyperpolarizability
!
! in/out : hpol - SFG hyperpolarizability tensor
!
! =============================================================================
    
    use constants 
    implicit none
    complex(kindr), intent(inout) :: hpol(nmodes,3,3,3) 
    integer :: nmodes,j,a,b,c
    real(kindr) :: width
    complex(kindr) :: z
    complex(kindr) :: htown
    complex(kindr), intent(inout) :: shpol(3,3,3)
    
    z = csqrt((-1,0))
    !Correct the static hyperpolarizability so it has proper units when run
    !through SFGFundAvg
    htown = ((width*z)/(hartree2cm**2))

    do  j=1,nmodes

        hpol(j,1,1,1) = hpol(j,1,1,1) + (shpol(1,1,1)* htown)
        hpol(j,1,1,2) = hpol(j,1,1,2) + (shpol(1,1,2)* htown)
        hpol(j,1,1,3) = hpol(j,1,1,3) + (shpol(1,1,3)* htown)
        
        hpol(j,1,2,1) = hpol(j,1,2,1) + (shpol(1,2,1)* htown)
        hpol(j,1,2,2) = hpol(j,1,2,2) + (shpol(1,2,2)* htown)
        hpol(j,1,2,3) = hpol(j,1,2,3) + (shpol(1,2,3)* htown)
        
        hpol(j,1,3,1) = hpol(j,1,3,1) + (shpol(1,3,1)* htown)
        hpol(j,1,3,2) = hpol(j,1,3,2) + (shpol(1,3,2)* htown)
        hpol(j,1,3,3) = hpol(j,1,3,3) + (shpol(1,3,3)* htown)
        
        hpol(j,2,1,1) = hpol(j,2,1,1) + (shpol(2,1,1)* htown)
        hpol(j,2,1,2) = hpol(j,2,1,2) + (shpol(2,1,2)* htown)
        hpol(j,2,1,3) = hpol(j,2,1,3) + (shpol(2,1,3)* htown)
        
        hpol(j,2,2,1) = hpol(j,2,2,1) + (shpol(2,2,1)* htown)
        hpol(j,2,2,2) = hpol(j,2,2,2) + (shpol(2,2,2)* htown)
        hpol(j,2,2,3) = hpol(j,2,2,3) + (shpol(2,2,3)* htown)

        hpol(j,2,3,1) = hpol(j,2,3,1) + (shpol(2,3,1)* htown)
        hpol(j,2,3,2) = hpol(j,2,3,2) + (shpol(2,3,2)* htown)
        hpol(j,2,3,3) = hpol(j,2,3,3) + (shpol(2,3,3)* htown)

        hpol(j,3,1,1) = hpol(j,3,1,1) + (shpol(3,1,1)* htown)
        hpol(j,3,1,2) = hpol(j,3,1,2) + (shpol(3,1,2)* htown)
        hpol(j,3,1,3) = hpol(j,3,1,3) + (shpol(3,1,3)* htown)

        hpol(j,3,2,1) = hpol(j,3,2,1) + (shpol(3,2,1)* htown)
        hpol(j,3,2,2) = hpol(j,3,2,2) + (shpol(3,2,2)* htown)
        hpol(j,3,2,3) = hpol(j,3,2,3) + (shpol(3,2,2)* htown)

        hpol(j,3,3,1) = hpol(j,3,3,1) + (shpol(3,3,1)* htown)
        hpol(j,3,3,2) = hpol(j,3,3,2) + (shpol(3,3,2)* htown)
        hpol(j,3,3,3) = hpol(j,3,3,3) + (shpol(3,3,3)* htown)
    end do
end subroutine SFGnonRes
