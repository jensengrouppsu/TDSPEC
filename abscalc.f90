Subroutine AbsPlot(lherztell, latermoff, x, w, tstep, osc, tdipder, &
                    freq, delta, lfreq, freq0, g, nexci, &
                    nmodes, nabspts)

  use constants

  implicit none

  integer :: j, i, t, n, r, k, p
  integer :: nexci, nabspts, nmodes, tstep

  complex(kindr) :: Ls, Lshtint, Lsht, Lsht1
  complex(kindr) :: lstmp
  complex(kindr) :: dsum(tstep,nexci)
  complex(kindr) :: rsum(nmodes,tstep,nexci)
  complex(kindr) :: Ls_htabs(nmodes,nmodes)
  complex(kindr) :: Ls_htb3(nmodes)
  complex(kindr) :: expfac
  complex(kindr) :: rsumexpfac

  real(kindr) :: freq(nmodes), freq0(nexci)
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: lfreq, lfreqt
  real(kindr) :: osc(nexci,3), tdipder(nexci,nmodes,3)
  real(kindr) :: delta(nexci,nmodes)
  real(kindr) :: osctdipder(nmodes,nexci)
  real(kindr) :: tdipdersq(nmodes,nexci)
  real(kindr) :: tdipdersqoffdiag(nmodes,nmodes,nexci)
  real(kindr) :: deltairsq(nmodes,nexci)
!  real(kindr) :: freqs(nabspts), absorption(nabspts)

  logical :: lherztell, latermoff

  complex(kindr) :: zzero = (0d0,0d0) 
  complex(kindr) :: z
  complex(kindr),intent(in) :: g(nexci,tstep)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Simulate one-photon absorbance spectra
!
! input  : lherztell - determines if Herzberg-Teller terms are calculated
!          latermoff - determines if only Herzberg-Teller terms are calculated
!          x - abcissa points for numerical integration
!          w - quadrature weights for numerical integration
!          tstep - number of abcissa points
!          osc - transition dipole moment components
!          tdipder - derivatives of transition dipole moments along each mode
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

  ! Evaluate the sum appearing in the integral equations of the 
  ! line shape function (see equation 20 from Neese/Petrenko).
  dsum = zzero
  do i=1,nexci
    do t=1,tstep
      do n=1,nmodes
        dsum(t,i) = dsum(t,i) + (delta(i,n)**2/two)*(one-exp(-z*freq(n)*x(t)))
        if (lherztell) then
          ! The rsum term may be computed with either a positive or negative term
          rsum(n,t,i) = delta(i,n)/(root2)*(one-exp(-z*freq(n)*x(t)))
        end if
      end do
    end do

    ! Calculate the square of the transition dipole derivatives
    if (lherztell) then
      do n=1,nmodes
        deltairsq(n,i)  = delta(i,n)**2
        osctdipder(n,i) = osc(i,1)*tdipder(i,n,1) + osc(i,2)*tdipder(i,n,2) &
                        + osc(i,3)*tdipder(i,n,3)
        tdipdersq(n,i)  = tdipder(i,n,1)*tdipder(i,n,1) + tdipder(i,n,2)*tdipder(i,n,2) &
                        + tdipder(i,n,3)*tdipder(i,n,3)
        do k=1,nmodes
          tdipdersqoffdiag(k,n,i) = tdipder(i,n,1)*tdipder(i,k,1) + tdipder(i,n,2)*tdipder(i,k,2) &
                                  + tdipder(i,n,3)*tdipder(i,k,3)
        end do
      end do
    end if

  end do

  ! Cycle over each absorption point
  do j = 1,nabspts
    lfreqt = lfreq + (j-1)*20.0_kindr ! 20 cm-1 ~ 1 nm at 700 nm.
    Ls = zzero
    Lshtint = zzero
    Lsht = zzero
    Lsht1 = zzero
    do i=1,nexci
      lstmp = zero
      if (lherztell) then
        Ls_htabs = zzero
        Ls_htb3 = zzero
      end if

      do t=1,tstep
        ! Evaluate the absorption cross section using equation 31 from 
        ! Neese and Petrenko.   

        expfac = w(t)*exp(x(t)*z*(lfreqt-freq0(i)) - g(i,t) - dsum(t,i))

        lstmp = lstmp + expfac

        ! Calculate the Herzberg-Teller contribution to the absorbance
        ! spectrum, if requested.

        if (lherztell) then

          ! Calculate the interference or  mixed FC/HT term and the
          ! HT term.  The interference term looks a lot like the A term for 
          ! resonance Raman scattering, while the HT term looks a lot
          ! like a combination band.

          !$OMP PARALLEL DO PRIVATE(p,rsumexpfac)
          do r=1,nmodes
            rsumexpfac = rsum(r,t,i) * expfac
            ! Calculate the interference term
            Ls_htb3(r) = Ls_htb3(r) + rsumexpfac
            ! Calculate the diagonal contribution to the HT term
            Ls_htabs(r,r) = Ls_htabs(r,r) + &
                   (one-deltairsq(r,i)+deltairsq(r,i)*cos(freq(r)*x(t)))* &
                           w(t)*exp(x(t)* &
                           z*(lfreqt-freq0(i)-freq(r)) - g(i,t) - dsum(t,i))
            do p=1,nmodes
              ! Calculate the off-diagonal contribution to the HT term
              if (r /= p) then
                Ls_htabs(p,r) = Ls_htabs(p,r) + rsum(p,t,i)*rsumexpfac
              end if
            end do
          end do
          !$OMP END PARALLEL DO

        end if ! lherztell
      end do ! t (tstep)

      ! The units of the absorption cross section should be those of 
      ! area (i.e. Angstrom^2).  The transition dipole has atomic 
      ! units, while the lineshape function is in units of cm.  There
      ! is also a constant in front of the integral:
      !
      ! (4*pi*lfreqt*alphaf*bohr2ang^2)/(3)
      !
      ! For this to work, the units of the laser frequency should be
      ! wavenumbers since that is what they need to be for the 
      ! exponential.  The final result is put into units of 
      ! Angstrom^2/molecule.  The wavenumbers on the horizontal axis 
      ! are scaled for plotting purposes.
      !
      ! Constants used (units are shown, when appropriate):
      ! bohr2ang => 0.5291772 Angstrom / Bohr
      ! abscrfac => (four*pi*alphaf)/(three) = 0.030567081 

      ! absconv = abscrfac*bohr2ang^2

      if (latermoff.eqv..false.) then
        Ls = Ls + &
            (osc(i,1)**2 + osc(i,2)**2 + osc(i,3)**2)*lstmp
      end if

      ! Add in the Herzberg-Teller contribution if requested.

      if (lherztell) then

        ! Contribution from the mixed FC/HT term (interference term).
        Lshtint = Lshtint + two*SUM(osctdipder(:,i)*Ls_htb3(:))

        do k=1,nmodes
          ! Contribution from the diagonal part of the HT term.
          Lsht = Lsht + tdipdersq(k,i)*Ls_htabs(k,k)
          ! Contribution from the off-diagonal part of the HT term.
          do p=1,nmodes
            if (k /= p) then
              Lsht1 = Lsht1 + tdipdersqoffdiag(p,k,i)*Ls_htabs(p,k)
            end if
          end do
        end do

      end if ! lherztell
    end do ! i (nexci)

    ! Output the excitation energy and absorbance cross section. 
    ! If you want wavelength (in nm) in the output, you need to 
    ! take 1.0E+7 divided by the laser frequency.

    ! The terms in the sum are:
    ! Ls - FC term
    ! Lshtint - Interference term (FC/HT term)
    ! Lsht - HT term (diagonal, where the modes are identical)
    ! Lsht1 - HT term (off-diagonal, where modes are different)

!    freqs(j) = 1.0E+7_kindr/lfreqt
!    absorption(j) = lfreqt*absconv*(Real(Ls)+ Real(Lshtint)+Real(Lsht)+Real(Lsht1))
    write(*,"(F15.10,10X,F22.18)") 1.0E+7_kindr/lfreqt, lfreqt*absconv*(Real(Ls)+ Real(Lshtint)+Real(Lsht)+Real(Lsht1))
  end do ! j (nabpts)

!  do j = 1,nabspts
!    write(*,*) freqs(j), absorption(j)
!  end do

End Subroutine AbsPlot
