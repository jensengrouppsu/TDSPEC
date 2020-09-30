Subroutine TwoAbsPlot2(lherztell, latermoff, x, w, tstep, &
             stpm, stpmder, freq, delta, lfreq, &
             freq0, g, nexci, nmodes, nabspts)

  use constants

  implicit none

  integer :: j, i, t, n, r, k, p, iex
  integer :: nexci, nabspts, nmodes, tstep

  complex(kindr) :: Ls, Lshtint, Lsht, Lsht1
  complex(kindr) :: dsum, rsum, lstmp, rsum1
  complex(kindr) :: Ls_htabs(nexci,nmodes,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), freq0(nexci)
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: lfreq
  real(kindr) :: stpm(nexci,6), stpmder(nexci,nmodes,6)
  real(kindr) :: delta(nexci,nmodes)
  real(kindr) :: df(nexci), dg(nexci)
  real(kindr) :: dfht(nexci,nmodes), dght(nexci,nmodes)
  real(kindr) :: dfht2(nexci,nmodes), dght2(nexci,nmodes)
  real(kindr) :: dfht3(nexci,nmodes,nmodes), dght3(nexci,nmodes,nmodes)

  logical :: lherztell, latermoff

  complex(kindr) :: zzero = (0d0,0d0) 
  complex(kindr) :: z

  z = csqrt((-1,0))

! =============================================================================
! purpose: Simulate two-photon absorbance spectra
!
! input  : lherztell - determines if Herzberg-Teller terms are calculated
!          latermoff - determines if only Herzberg-Teller terms are calculated
!          x - abcissa points for numerical integration
!          w - quadrature weights for numerical integration
!          tstep - number of abcissa points
!          stpm - two-photon transition moment components
!          stpmder - derivatives of two-photon transition moments along each mode
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

  ! Calculate the averages of the two-photon transition dipole 
  ! moments.  Referring to
  !
  ! McClain, W. M.  J. Chem. Phys., 55, 2789, (1971).
  !
  ! the elements df and dg are molecular parameters of the two-photon
  ! absorption probability, which are independent of the incident 
  ! radiation.  Based on the Macak et. al. paper, the two-photon 
  ! absorption probability is given by (for linear polarized light):
  !
  ! dtp = 2*sum(a,b){Saa*Sbb + Sab*Sab + Sab*Sba}
  !     = sum(a,b){2*Saa*Sbb + 4*Sab*Sab}
  !
  ! Where the sum is over the x, y, and z directions of the 
  ! two-photon transition dipole moment tensor.  Note that there is 
  ! no restriction on the sum, so both a and b can be equal.  It turns
  ! out in slightly different notation the two-photon absorption 
  ! probability is:
  !
  ! dtp = 2*df + 4*dg
  ! df = sum(a,b){Saa*Sbb}
  ! dg = sum(a,b){Sab*Sab} = sum(a,b){Sab*Sba}
  ! 
  ! In DALTON, the orientational averaging is included by using the 
  ! factor of 1/30, which is built in to the df and dg coefficients.

  do iex=1,nexci
    df(iex) = stpm(iex,1)**2 + stpm(iex,2)**2 + stpm(iex,3)**2
    dg(iex) = df(iex)

    df(iex) = df(iex) + two*stpm(iex,1)*stpm(iex,2)
    df(iex) = df(iex) + two*stpm(iex,1)*stpm(iex,3)
    df(iex) = df(iex) + two*stpm(iex,2)*stpm(iex,3)
    df(iex) = df(iex)/thrty

    dg(iex) = dg(iex) + two*stpm(iex,4)*stpm(iex,4)
    dg(iex) = dg(iex) + two*stpm(iex,5)*stpm(iex,5)
    dg(iex) = dg(iex) + two*stpm(iex,6)*stpm(iex,6)
    dg(iex) = dg(iex)/thrty
  end do

  ! Do the averaging for Herzberg-Teller coupling if it's requested. 
  ! This involves multiplying the two-photon moment times its 
  ! derivative, similar to what was done above.  This part is mode 
  ! dependent, however (thus the extra loop).

  if (lherztell) then
    do iex=1,nexci
      do k=1,nmodes

        ! Contribution to the mixed FC/HT term.

        dfht(iex,k) = stpm(iex,1)*stpmder(iex,k,1) + &
                      stpm(iex,2)*stpmder(iex,k,2) + &
                      stpm(iex,3)*stpmder(iex,k,3)
        dght(iex,k) = dfht(iex,k)

        dfht(iex,k) = dfht(iex,k) + stpmder(iex,k,1)*stpm(iex,2)
        dfht(iex,k) = dfht(iex,k) + stpmder(iex,k,1)*stpm(iex,3)
        dfht(iex,k) = dfht(iex,k) + stpmder(iex,k,2)*stpm(iex,1)
        dfht(iex,k) = dfht(iex,k) + stpmder(iex,k,2)*stpm(iex,3)
        dfht(iex,k) = dfht(iex,k) + stpmder(iex,k,3)*stpm(iex,1)
        dfht(iex,k) = dfht(iex,k) + stpmder(iex,k,3)*stpm(iex,2)
        dfht(iex,k) = dfht(iex,k)/thrty

        dght(iex,k) = dght(iex,k) + two*stpmder(iex,k,4)*stpm(iex,4)
        dght(iex,k) = dght(iex,k) + two*stpmder(iex,k,5)*stpm(iex,5)
        dght(iex,k) = dght(iex,k) + two*stpmder(iex,k,6)*stpm(iex,6)
        dght(iex,k) = dght(iex,k)/thrty

        ! Contribution from the diagonal part of the HT term.

        dfht2(iex,k) = stpmder(iex,k,1)*stpmder(iex,k,1) + &
                       stpmder(iex,k,2)*stpmder(iex,k,2) + &
                       stpmder(iex,k,3)*stpmder(iex,k,3)
        dght2(iex,k) = dfht2(iex,k)

        dfht2(iex,k) = dfht2(iex,k) + stpmder(iex,k,1)* &
                                      stpmder(iex,k,2)
        dfht2(iex,k) = dfht2(iex,k) + stpmder(iex,k,1)* &
                                      stpmder(iex,k,3)
        dfht2(iex,k) = dfht2(iex,k) + stpmder(iex,k,2)* &
                                      stpmder(iex,k,1)
        dfht2(iex,k) = dfht2(iex,k) + stpmder(iex,k,2)* &
                                      stpmder(iex,k,3)
        dfht2(iex,k) = dfht2(iex,k) + stpmder(iex,k,3)* &
                                      stpmder(iex,k,1)
        dfht2(iex,k) = dfht2(iex,k) + stpmder(iex,k,3)* &
                                      stpmder(iex,k,2)
        dfht2(iex,k) = dfht2(iex,k)/thrty

        dght2(iex,k) = dght2(iex,k) + two*stpmder(iex,k,4)* &
                                          stpmder(iex,k,4)
        dght2(iex,k) = dght2(iex,k) + two*stpmder(iex,k,5)* &
                                          stpmder(iex,k,5)
        dght2(iex,k) = dght2(iex,k) + two*stpmder(iex,k,6)* &
                                          stpmder(iex,k,6)
        dght2(iex,k) = dght2(iex,k)/thrty

        ! Contribution from the off-diagonal part of the HT term.
        ! We will avoid double counting here, since the contribution
        ! where the modes are equal are held in dght2.

        do p=1,nmodes
          if (p == k) then
            cycle
          else
            dfht3(iex,k,p) = stpmder(iex,k,1)*stpmder(iex,p,1) + &
                             stpmder(iex,k,2)*stpmder(iex,p,2) + &
                             stpmder(iex,k,3)*stpmder(iex,p,3)
            dght3(iex,k,p) = dfht3(iex,k,p)

            dfht3(iex,k,p) = dfht3(iex,k,p) + stpmder(iex,k,1)* &
                                          stpmder(iex,p,2)
            dfht3(iex,k,p) = dfht3(iex,k,p) + stpmder(iex,k,1)* &
                                          stpmder(iex,p,3)
            dfht3(iex,k,p) = dfht3(iex,k,p) + stpmder(iex,k,2)* &
                                          stpmder(iex,p,1)
            dfht3(iex,k,p) = dfht3(iex,k,p) + stpmder(iex,k,2)* &
                                          stpmder(iex,p,3)
            dfht3(iex,k,p) = dfht3(iex,k,p) + stpmder(iex,k,3)* &
                                          stpmder(iex,p,1)
            dfht3(iex,k,p) = dfht3(iex,k,p) + stpmder(iex,k,3)* &
                                          stpmder(iex,p,2)
            dfht3(iex,k,p) = dfht3(iex,k,p)/thrty

            dght3(iex,k,p) = dght3(iex,k,p) + two*stpmder(iex,k,4)* &
                                              stpmder(iex,p,4)
            dght3(iex,k,p) = dght3(iex,k,p) + two*stpmder(iex,k,5)* &
                                              stpmder(iex,p,5)
            dght3(iex,k,p) = dght3(iex,k,p) + two*stpmder(iex,k,6)* &
                                              stpmder(iex,p,6)
            dght3(iex,k,p) = dght3(iex,k,p)/thrty
          end if
        end do
      end do
    end do
  end if

  ! Perform the same integration used to get the lineshape for 
  ! one-photon absorption to get the lineshape function for two-photon
  ! absorption.  The only difference is that an incident frequency of
  ! twice the incident frequency from the laser is used.

  do j = 1,nabspts
    Ls = zzero
    Lshtint = zzero
    Lsht = zzero
    Lsht1 = zzero
    do i=1,nexci
      lstmp = zero
      if (lherztell) then
        Ls_htb3(1:nexci,1:nmodes) = zero
        Ls_htabs(1:nexci,1:nmodes,1:nmodes) = zero
      end if
      do t=1,tstep
        dsum = zzero
        do n=1,nmodes
          dsum = dsum + (delta(i,n)**2/two)*(one-exp(-z*freq(n)*x(t)))
        end do
        lstmp = lstmp + w(t)*exp(x(t)* &
          z*(two*lfreq-freq0(i)) - g(i,t) - dsum)

        ! Calculate the Herzberg-Teller contribution to the absorbance
        ! spectrum, if requested.

        if (lherztell) then
          ! Calculate the interference or  mixed FC/HT term and the
          ! HT term.  The interference term looks a lot like the A term for 
          ! resonance Raman scattering, while the HT term looks a lot
          ! like a combination band.
          do r=1,nmodes
            ! Calculate the interference term
            ! Negative sign?
!            rsum = -delta(i,r)/(sqrt(two))*(one-exp(-z*freq(r)*x(t)))
            ! Positive sign?
            rsum = delta(i,r)/(sqrt(two))*(one-exp(-z*freq(r)*x(t)))
            Ls_htb3(i,r) = Ls_htb3(i,r) + w(t)*rsum*exp(x(t)* &
                z*(two*lfreq-freq0(i)) - g(i,t) - dsum)
            ! Calculate the HT term
            ! Diagonal contribution to the HT term (raising/lowering 
            ! the same mode).
            Ls_htabs(i,r,r) = Ls_htabs(i,r,r) + &
!                           (one-delta(i,r)**2)*w(t)*exp(x(t)* &
                    (one-delta(i,r)**2+delta(i,r)**2*cos(freq(r)*x(t)))* &
                            w(t)*exp(x(t)* &
                            z*(two*lfreq-freq0(i)-freq(r)) - g(i,t) - &
                            dsum)
            do p=1,nmodes
              ! Off-Diagonal contribution to the HT term (raising/
              ! lowering two different modes).
              if (r /= p) then
                ! Negative sign?
!                rsum1 = -delta(i,p)/(sqrt(two))* &
!                        (one-exp(-z*freq(p)*x(t)))
                ! Positive sign?
                rsum1 = delta(i,p)/(sqrt(two))* &
                        (one-exp(-z*freq(p)*x(t)))
                Ls_htabs(i,r,p) = Ls_htabs(i,r,p) + &
                  w(t)*rsum*rsum1*exp(x(t)* &
                  z*(two*lfreq-freq0(i)) - g(i,t) - dsum)
              end if
            end do
          end do
        end if
      end do

      ! Add in the correct units to the two-photon absorption 
      ! cross-section.  According to:
      !
      ! Norman, P.; Luo, Y.; and Agren, H.  J. Chem. Phys., 111, 
      !   7758, (1999).
      ! Norman, P.; Cronstrand, P.; and Ericsson, J.  Chem. Phys., 
      !   285, 207, (2002).
      !
      ! The two-photon absorption cross-section is given by:
      !
      ! \sigma_{TPA} = (8*pi^3*alphaf^2)*lfreq^2*lstmp*dtp
      !
      ! Note that this equation appears different than what is 
      ! presented in Eq. 9 of the second reference, but it actually
      ! is the same.  In this particular equation, every variable 
      ! needs to be placed into atomic units.  In particular, note 
      ! that the laser frequency lfreq needs to be in per time units 
      ! (i.e. hbar*lfreq is the photon energy), so this must
      ! be converted properly.  The latter reference states that if 
      ! this calculation is performed in a.u., a conversion of
      !
      ! au2gm = (5.291772E-9 cm/Bohr)^4*(2.418884E-17 s/a.u.)
      !
      ! will give the correct TPA cross section units of 
      ! (cm^4 s)/photon.
      !
      ! hartree2cm = 219474 cm^{-1}/Hartree
      !
      ! Note that in atomic units, the energy is identical to the
      ! angular frequency (since hbar = 1).
      !
      ! Based on other people's work, it seems like the conversion
      ! above has an extra factor of 2, so it should include a 
      ! 4 instead of an 8.  This is based on counting the number
      ! of photons instead of counting photons as pairs.  We will 
      ! use the conversion with a factor of 4 to be consistent
      ! with most other work.

      if (latermoff.eqv..false.) then
        Ls = Ls + (two*df(i) + four*dg(i))*lstmp
      end if

      ! Calculate the Herzberg-Teller term if requested.

      if (lherztell) then
        do k=1,nmodes
          ! Contribution from the mixed FC/HT term
          Lshtint = Lshtint + &
               (four*dfht(i,k)+eight*dght(i,k))*Ls_htb3(i,k)
          ! Contribution from the diagonal part of the HT term.
          Lsht = Lsht + &
             (two*dfht2(i,k)+four*dght2(i,k))*Ls_htabs(i,k,k)
          do p=1,nmodes
            ! Contribution from the diagonal part of the HT term.
            if (k /= p) then
              Lsht1 = Lsht1 + &
                 (two*dfht3(i,k,p)+four*dght3(i,k,p))* &
                 Ls_htabs(i,k,p)
            end if
          end do
        end do
      end if
    end do

    ! wvn2hz2au = (100 cm/m)*(2.99792458E+8 m/s)*(2.418884E-17 s/a.u.)

    ! tpa1 = 4*pi^3*alpha^2*(219474 cm^{-1}/Hartree)*
    !        (5.291772E-9 cm/Bohr)^4*(2.418884E-17 s/a.u.)
    ! tpa2 = 4*pi^2*(100 cm/m)*(2.99792458E+8 m/s)*(2.418884E-17 s/a.u.)
    ! tpaconv = tpa1*tpa2

    ! Extra factor of (2*pi)^2 results from the conversion of
    ! wavenumbers to Hz to orbital period (a.u.).  This conversion
    ! does not include the factor of 2*pi for the angular frequency
    ! (i.e. E = hbar*omega), so it needs to be added externally.

    ! The terms in the sum are:
    ! Ls - FC term
    ! Lshtint - Interference term (FC/HT term)
    ! Lsht - HT term (diagonal, where the modes are identical)
    ! Lsht1 - HT term (off-diagonal, where modes are different)

    write(*,*) 1.0E+7_kindr/(two*lfreq), &
               tpaconv*lfreq**2*(Real(Ls)+Real(Lshtint)+ &
               Real(Lsht)+Real(Lsht1))
    lfreq = lfreq + two
  end do

End Subroutine TwoAbsPlot2
