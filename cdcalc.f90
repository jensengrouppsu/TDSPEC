Subroutine CDPlot(lherztell, latermoff, x, w, tstep, osc, mtdip, &
                    freq, delta, lfreq, freq0, g, nexci, &
                    nmodes, nabspts)

  use constants

  implicit none

  integer :: j, i, t, n, r, k, p
  integer :: nexci, nabspts, nmodes, tstep

  complex(kindr) :: Ls, Lshtint, Lsht, Lsht1
  complex(kindr) :: dsum, rsum, lstmp, rsum1
!  complex(kindr) :: Ls_htabs(nexci,nmodes,nmodes)
!  complex(kindr) :: Ls_htb3(nexci,nmodes)

  real(kindr) :: freq(nmodes), freq0(nexci)
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: lfreq
  real(kindr) :: osc(nexci,3), mtdip(nexci,3)
  real(kindr) :: delta(nexci,nmodes)

  logical :: lherztell, latermoff

  complex(kindr) :: zzero = (0d0,0d0) 
  complex(kindr) :: z
  complex(kindr),intent(in) :: g(nexci,tstep)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Simulate circular dichroism spectra
!
! input  : lherztell - determines if Herzberg-Teller terms are calculated
!          latermoff - determines if only Herzberg-Teller terms are calculated
!          x - abcissa points for numerical integration
!          w - quadrature weights for numerical integration
!          tstep - number of abcissa points
!          osc - transition dipole moment components
!          mtdip - magnetic transition dipole moment components
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          lfreq - laser frequency
!          freq0 - vertical excitation energy
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          nabspts - number of points used to make an absorbance spectrum 
!
! =============================================================================

  do j = 1,nabspts
    Ls = zzero
    Lshtint = zzero
    Lsht = zzero
    Lsht1 = zzero
    do i=1,nexci
      lstmp = zero
!      if (lherztell) then
!        Ls_htabs(1:nexci,1:nmodes,1:nmodes) = zero
!        Ls_htb3(1:nexci,1:nmodes) = zero
!      end if
      do t=1,tstep
        dsum = zzero
        do n=1,nmodes
          dsum = dsum + (delta(i,n)**2/two)*(one-exp(-z*freq(n)*x(t)))
        end do

        ! Evaluate the circular dichroism intensity using equation 31 from 
        ! Neese and Petrenko.   

        lstmp = lstmp + w(t)*exp(x(t)* &
               z*(lfreq-freq0(i)) - g(i,t) - dsum)

!        ! Calculate the Herzberg-Teller contribution to the circular dichroism
!        ! spectrum, if requested.
!
!        if (lherztell) then
!
!          ! Calculate the interference or  mixed FC/HT term and the
!          ! HT term.  The interference term looks a lot like the A term for 
!          ! resonance Raman scattering, while the HT term looks a lot
!          ! like a combination band.
!          do r=1,nmodes
!            ! Calculate the interference term
!            ! With negative sign?
!!            rsum = -delta(i,r)/(sqrt(two))*(one-exp(-z*freq(r)*x(t)))
!            ! No negative sign?
!            rsum = delta(i,r)/(sqrt(two))*(one-exp(-z*freq(r)*x(t)))
!            Ls_htb3(i,r) = Ls_htb3(i,r) + w(t)*rsum*exp(x(t)* &
!                z*(lfreq-freq0(i)) - g(i,t) - dsum)
!            ! Calculate the diagonal contribution to the HT term
!            Ls_htabs(i,r,r) = Ls_htabs(i,r,r) + &
!                   (one-delta(i,r)**2+delta(i,r)**2*cos(freq(r)*x(t)))* &
!                           w(t)*exp(x(t)* &
!                           z*(lfreq-freq0(i)-freq(r)) - g(i,t) - dsum)
!            do p=1,nmodes
!              ! Calculate the off-diagonal contribution to the HT term
!              if (r /= p) then
!                ! With negative sign?
!!                rsum1 = -delta(i,p)/(sqrt(two))*(one-exp(-z*freq(p)*x(t)))
!                ! No negative sign?
!                rsum1 = delta(i,p)/(sqrt(two))*(one-exp(-z*freq(p)*x(t)))
!                Ls_htabs(i,r,p) = Ls_htabs(i,r,p) + &
!                  w(t)*rsum*rsum1*exp(x(t)* &
!                  z*(lfreq-freq0(i)) - g(i,t) - dsum)
!              end if
!            end do
!          end do
!        end if
      end do

      if (latermoff.eqv..false.) then
        Ls = Ls + &
            (osc(i,1)*mtdip(i,1) + osc(i,2)*mtdip(i,2) + &
             osc(i,3)*mtdip(i,3))*lstmp
      end if

      ! Add in the Herzberg-Teller contribution if requested.

!      if (lherztell) then
!        do k=1,nmodes
!          ! Contribution from the mixed FC/HT term (interference term).
!          Lshtint = Lshtint + two* &
!               (osc(i,1)*tdipder(i,k,1)+ &
!               osc(i,2)*tdipder(i,k,2)+osc(i,3)*tdipder(i,k,3))* &
!               Ls_htb3(i,k)
!          ! Contribution from the diagonal part of the HT term.
!          Lsht = Lsht + &
!             (tdipder(i,k,1)*tdipder(i,k,1)+ &
!             tdipder(i,k,2)*tdipder(i,k,2)+ &
!             tdipder(i,k,3)*tdipder(i,k,3))* &
!             Ls_htabs(i,k,k)
!          do p=1,nmodes
!            ! Contribution from the off-diagonal part of the HT term.
!            if (k /= p) then
!              Lsht1 = Lsht1 + &
!                 (tdipder(i,k,1)*tdipder(i,p,1)+ &
!                 tdipder(i,k,2)*tdipder(i,p,2)+ &
!                 tdipder(i,k,3)*tdipder(i,p,3))* &
!                 Ls_htabs(i,k,p)
!            end if
!          end do
!        end do
!      end if
    end do

    ! Output the excitation energy and absorbance cross section. 
    ! If you want wavelength (in nm) in the output, you need to 
    ! take 1.0E+7 divided by the laser frequency.

    ! The terms in the sum are:
    ! Ls - FC term
    ! Lshtint - Interference term (FC/HT term)
    ! Lsht - HT term (diagonal, where the modes are identical)
    ! Lsht1 - HT term (off-diagonal, where modes are different)

    write(*,*) 1.0E+7_kindr/lfreq,lfreq*Real(Ls)
!    write(*,*) 1.0E+7_kindr/lfreq,lfreq*absconv*(Real(Ls)+ &
!                                  Real(Lshtint)+Real(Lsht)+Real(Lsht1))
    lfreq = lfreq + 2.0
  end do

End Subroutine CDPlot
