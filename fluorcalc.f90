Subroutine FluorPlot(fluorstate, x, w, tstep, osc, freq, delta, lfreq, &
               freq0, g, nexci, nmodes, nabspts, stokes)

  use constants

  implicit none

  integer :: j, t, n
  integer :: nexci, nabspts, nmodes, tstep
  integer :: fluorstate

  complex(kindr) :: Ls
  complex(kindr) :: lstmp
  complex(kindr) :: dsum(tstep)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), freq0(nexci)
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: lfreq, stokes
  real(kindr) :: osc(nexci,3), delta(nexci,nmodes)

  complex(kindr) :: zzero = (0d0,0d0) 
  complex(kindr) :: z

  z = csqrt((-1,0))

! =============================================================================
! purpose: Simulate fluorescence spectra
!
! input  : fluorstate - excited state fluorescence is based on
!          x - abcissa points for numerical integration
!          w - quadrature weights for numerical integration
!          tstep - number of abcissa points
!          osc - transition dipole moment components
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          lfreq - laser frequency
!          freq0 - vertical excitation energy
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          nabspts - number of points used to make a fluorescence spectrum 
!          stokes - Stokes shift 
!
! =============================================================================

! For details on evaluating the fluorescence, see:
!
! T. Petrenko and F. Neese.  J. Chem. Phys., 127, 164319 (2007).
! T. Petrenko, O. Krylova, F. Neese, and M. Sokolowski.  New J. Phys.,
!    11, 015001 (2009).
! D. Tannor.  Introduction to Quantum Mechanics: A Time-Dependent 
!    Perspective.  (University Science Books, Sausalito, CA, 2007).
!    Chapter 14.  pp. 441-445.

  dsum = zzero
  do t=1,tstep
    do n=1,nmodes
      ! Evaluate the sum appearing in the integral equations of the 
      ! line shape function (see equation 20 from Neese/Petrenko).
      dsum(t) = dsum(t) + (delta(fluorstate,n)**2/two)*(one-exp(-z*freq(n)*x(t)))
    end do
  end do

  do j = 1,nabspts
    Ls = zzero
    lstmp = zero
    do t=1,tstep
      ! Evaluate the fluorescence efficiency using equation 4 from 
      ! Neese and Petrenko. 

      ! The New J. Phys. paper above yields a total Stokes shift by 
      ! giving 50% of the shift to absorbance and 50% of it to 
      ! fluorescence.  Since I give no Stokes shift to the absorbance
      ! spectrum, the present method involves a total shift defined
      ! by moving the placement of the fluorescence spectrum only.  To
      ! use Neese and Petrenko's method, divide the stokes parameter
      ! by 2 and also place an equivalent shift in the absorbance 
      ! routine.

      lstmp = lstmp + w(t)*exp(x(t)* &
             z*(freq0(fluorstate)-lfreq-stokes)-conjg(g(fluorstate,t)) &
             -dsum(t))

    end do

    ! Based on the book by Tannor, the quantity we can directly 
    ! calculate for fluorescence is the second derivative of the 
    ! emission rate (k) with respect to the emitted (scattered) photon
    ! energy and solid angle.  Thus, we have:
    !
    ! d2k/(dEs*dOmega) = (1/(3*pi^2*hbar^2*e^2))*(omega/c)^3*tdip^2*Ls
    !
    ! In comparison to Eq. 14.35 in the book by Tannor, we include an 
    ! extra factor of e^2.  This is required to cancel the charge units 
    ! from the transition dipole moment, otherwise the derivative does 
    ! not have sensible units.  Examining the left side of the equation,
    ! the rate has units of s^{-1}, energy Es is in Joules, and the 
    ! solid angle is effectively unitless (steradians).  So the right
    ! side of the expression must have units of s^{-1} J^{-1} sr^{-1}.
    ! Alternatively we could integrate over the solid angle, obtain
    ! an additional factor of 4*pi sr, and multiply this into the
    ! expression:
    !
    ! dk/dEs = (4/(3*pi*hbar^2*e^2))*(omega/c)^3*tdip^2*Ls 
    ! 
    ! This seems better so we use this convention.  Temporarily I'll
    ! call this the differential fluorescence emission rate. 
    
    ! planckbar  => 6.6260755E-34/(two*pi) J*s
    ! speed      => 2.99792458E+8 m / s
    ! bohr2ang   => 0.5291772 Angstrom / Bohr
    ! wvn2joules => 1.9864456E-23 J / cm^{-1}
    ! joule2ev   => 6.24152862E+18 eV / J

    ! fl1 = 4*c^3*(100 cm/1 m)^3
    ! fl2 = 3*pi*hbar^2*c^3*(100 cm/1 m)^3
    ! fl3 = bohr2ang^2*(1 m/100 cm)^8*wvn2joules
    ! fl4 = c*(1 m/100 cm)*joule2ev
    ! flourconv = (fl1*fl3)/(fl2*fl4)

    !Ls = Ls + (four/(three*pi*planckbar**2*speed**3*hndr**3))* &
    !    (osc(fluorstate,1)**2 + osc(fluorstate,2)**2 + &
    !    osc(fluorstate,3)**2)*lstmp* &
    !    (bohr2ang**2*(1E-8_kindr)**2*wvn2joules)/ &
    !    (speed*hndr*joule2ev)
    Ls = Ls + &
        (osc(fluorstate,1)**2 + osc(fluorstate,2)**2 + &
        osc(fluorstate,3)**2)*lstmp

    ! Output the excitation energy and absorbance cross section. 
    ! If you want wavelength (in nm) in the output, you need to 
    ! take 1.0E+7 divided by the laser frequency.

    ! The terms in the sum are:
    ! Ls - FC term
    !
    ! We convert the emitted frequency "lfreq" from wavenumbers
    ! to Hertz as the final step.  We also normalize the units
    ! into per eV instead of per Joule

    !write(*,*) 1.0E+7_kindr/lfreq, (lfreq*speed*hndr)**3*Real(Ls)
    write(*,*) 1.0E+7_kindr/lfreq, fluorconv*lfreq**3*Real(Ls)
    lfreq = lfreq - two
  end do

End Subroutine FluorPlot
