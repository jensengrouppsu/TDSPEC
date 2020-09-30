subroutine RHRSFCFund(latermoff, x, w, freq, delta, osc, stpm, hpol, &
                   freq0, lfreq, g, nexci, nmodes, tstep)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep

  complex(kindr), intent(inout) :: hpol(nmodes,3,3,3)
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
! purpose: Evaluate A term hyperpolarizability for RHRS spectra (fundamentals)
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
! in/out : hpol - RHRS hyperpolarizability tensor
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
               z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
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


        hpol(j,1,1,1) = hpol(j,1,1,1) + osc(iex,1)*stpm(iex,1)*Ls* &
          hpol2si
        hpol(j,1,1,2) = hpol(j,1,1,2) + osc(iex,1)*stpm(iex,4)*Ls* &
          hpol2si
        hpol(j,1,1,3) = hpol(j,1,1,3) + osc(iex,1)*stpm(iex,5)*Ls* &
          hpol2si

        hpol(j,1,2,1) = hpol(j,1,2,1) + osc(iex,1)*stpm(iex,4)*Ls* &
          hpol2si
        hpol(j,1,2,2) = hpol(j,1,2,2) + osc(iex,1)*stpm(iex,2)*Ls* &
          hpol2si
        hpol(j,1,2,3) = hpol(j,1,2,3) + osc(iex,1)*stpm(iex,6)*Ls* &
          hpol2si

        hpol(j,1,3,1) = hpol(j,1,3,1) + osc(iex,1)*stpm(iex,5)*Ls* &
          hpol2si
        hpol(j,1,3,2) = hpol(j,1,3,2) + osc(iex,1)*stpm(iex,6)*Ls* &
          hpol2si
        hpol(j,1,3,3) = hpol(j,1,3,3) + osc(iex,1)*stpm(iex,3)*Ls* &
          hpol2si

        hpol(j,2,1,1) = hpol(j,2,1,1) + osc(iex,2)*stpm(iex,1)*Ls* &
          hpol2si
        hpol(j,2,1,2) = hpol(j,2,1,2) + osc(iex,2)*stpm(iex,4)*Ls* &
          hpol2si
        hpol(j,2,1,3) = hpol(j,2,1,3) + osc(iex,2)*stpm(iex,5)*Ls* &
          hpol2si

        hpol(j,2,2,1) = hpol(j,2,2,1) + osc(iex,2)*stpm(iex,4)*Ls* &
          hpol2si
        hpol(j,2,2,2) = hpol(j,2,2,2) + osc(iex,2)*stpm(iex,2)*Ls* &
          hpol2si
        hpol(j,2,2,3) = hpol(j,2,2,3) + osc(iex,2)*stpm(iex,6)*Ls* &
          hpol2si

        hpol(j,2,3,1) = hpol(j,2,3,1) + osc(iex,2)*stpm(iex,5)*Ls* &
          hpol2si
        hpol(j,2,3,2) = hpol(j,2,3,2) + osc(iex,2)*stpm(iex,6)*Ls* &
          hpol2si
        hpol(j,2,3,3) = hpol(j,2,3,3) + osc(iex,2)*stpm(iex,3)*Ls* &
          hpol2si

        hpol(j,3,1,1) = hpol(j,3,1,1) + osc(iex,3)*stpm(iex,1)*Ls* &
          hpol2si
        hpol(j,3,1,2) = hpol(j,3,1,2) + osc(iex,3)*stpm(iex,4)*Ls* &
          hpol2si
        hpol(j,3,1,3) = hpol(j,3,1,3) + osc(iex,3)*stpm(iex,5)*Ls* &
          hpol2si

        hpol(j,3,2,1) = hpol(j,3,2,1) + osc(iex,3)*stpm(iex,4)*Ls* &
          hpol2si
        hpol(j,3,2,2) = hpol(j,3,2,2) + osc(iex,3)*stpm(iex,2)*Ls* &
          hpol2si
        hpol(j,3,2,3) = hpol(j,3,2,3) + osc(iex,3)*stpm(iex,6)*Ls* &
          hpol2si

        hpol(j,3,3,1) = hpol(j,3,3,1) + osc(iex,3)*stpm(iex,5)*Ls* &
          hpol2si
        hpol(j,3,3,2) = hpol(j,3,3,2) + osc(iex,3)*stpm(iex,6)*Ls* &
          hpol2si
        hpol(j,3,3,3) = hpol(j,3,3,3) + osc(iex,3)*stpm(iex,3)*Ls* &
          hpol2si

      end do

      ! Test
      !write(6,*) 'A term hyperpolarizability'
      !write(6,*) freq(j)
      !!------------------------------
      !write(6,*) 'xxx', hpol(j,1,1,1)
      !write(6,*) 'xxy', hpol(j,1,1,2)
      !write(6,*) 'xxz', hpol(j,1,1,3)
      !write(6,*) 'xyx', hpol(j,1,2,1)
      !write(6,*) 'xyy', hpol(j,1,2,2)
      !write(6,*) 'xyz', hpol(j,1,2,3)
      !write(6,*) 'xzx', hpol(j,1,3,1)
      !write(6,*) 'xzy', hpol(j,1,3,2)
      !write(6,*) 'xzz', hpol(j,1,3,3)
      !!------------------------------
      !write(6,*) 'yxx', hpol(j,2,1,1)
      !write(6,*) 'yxy', hpol(j,2,1,2)
      !write(6,*) 'yxz', hpol(j,2,1,3)
      !write(6,*) 'yyx', hpol(j,2,2,1)
      !write(6,*) 'yyy', hpol(j,2,2,2)
      !write(6,*) 'yyz', hpol(j,2,2,3)
      !write(6,*) 'yzx', hpol(j,2,3,1)
      !write(6,*) 'yzy', hpol(j,2,3,2)
      !write(6,*) 'yzz', hpol(j,2,3,3)
      !!------------------------------
      !write(6,*) 'zxx', hpol(j,3,1,1)
      !write(6,*) 'zxy', hpol(j,3,1,2)
      !write(6,*) 'zxz', hpol(j,3,1,3)
      !write(6,*) 'zyx', hpol(j,3,2,1)
      !write(6,*) 'zyy', hpol(j,3,2,2)
      !write(6,*) 'zyz', hpol(j,3,2,3)
      !write(6,*) 'zzx', hpol(j,3,3,1)
      !write(6,*) 'zzy', hpol(j,3,3,2)
      !write(6,*) 'zzz', hpol(j,3,3,3)
      !write(6,*) '------------------------------'

    end do
  end if

end subroutine RHRSFCFund

subroutine RHRSFCOverCB(latermoff, x, w, freq, delta, osc, stpm, &
                     hpolot, hpolcb, freq0, lfreq, g, nexci, &
                     nmodes, tstep, excnum, numcb, excnums, factorial)

  use constants

  implicit none

  integer :: nexci, nmodes, excnum, numcb
  integer :: j, iex, t, i, tstep, r, k, m
  integer :: excnums(numcb,nmodes)

  complex(kindr), intent(inout) :: hpolot(excnum,nmodes,3,3,3)
  complex(kindr), intent(inout) :: hpolcb(numcb,3,3,3)
  complex(kindr) :: Ls, rsum, dsum, z
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3), stpm(nexci,6)
  real(kindr) :: factorial(16)

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: latermoff

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate A term hyperpolarizability for RHRS spectra (overtones and
!          combination bands)
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
! in/out : hpolot - RHRS hyperpolarizability tensor for overtones
!          hpolcb - RHRS hyperpolarizability tensor for combination bands
!
! =============================================================================

  if (excnum > 1) then
    ! Calculate the transition hyperpolarizability for the A term for 
    ! overtones
    if (latermoff.eqv..false.) then
      do k=2,excnum
        do j = 1, nmodes
          do iex=1,nexci
            Ls = zzero
            do t=1,tstep
              dsum = zzero
              do i=1,nmodes
                dsum = dsum + &
                       (delta(iex,i)**2/two)* &
                       (one-exp(-z*freq(i)*x(t)))
              end do
              rsum = (delta(iex,j))**k/ &
                     (sqrt(two**k*factorial(k)))* &
                     (one-exp(-z*freq(j)*x(t)))**k
              Ls = Ls + z*w(t)*rsum*exp(x(t)* &
                   z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
            end do

            hpolot(k,j,1,1,1) = hpolot(k,j,1,1,1) + &
                osc(iex,1)*stpm(iex,1)*Ls* &
                hpol2si
            hpolot(k,j,1,1,2) = hpolot(k,j,1,1,2) + &
                osc(iex,1)*stpm(iex,4)*Ls* &
                hpol2si
            hpolot(k,j,1,1,3) = hpolot(k,j,1,1,3) + &
                osc(iex,1)*stpm(iex,5)*Ls* &
                hpol2si

            hpolot(k,j,1,2,1) = hpolot(k,j,1,2,1) + &
                osc(iex,1)*stpm(iex,4)*Ls* &
                hpol2si
            hpolot(k,j,1,2,2) = hpolot(k,j,1,2,2) + &
                osc(iex,1)*stpm(iex,2)*Ls* &
                hpol2si
            hpolot(k,j,1,2,3) = hpolot(k,j,1,2,3) + &
                osc(iex,1)*stpm(iex,6)*Ls* &
                hpol2si

            hpolot(k,j,1,3,1) = hpolot(k,j,1,3,1) + &
                osc(iex,1)*stpm(iex,5)*Ls* &
                hpol2si
            hpolot(k,j,1,3,2) = hpolot(k,j,1,3,2) + &
                osc(iex,1)*stpm(iex,6)*Ls* &
                hpol2si
            hpolot(k,j,1,3,3) = hpolot(k,j,1,3,3) + &
                osc(iex,1)*stpm(iex,3)*Ls* &
                hpol2si

            hpolot(k,j,2,1,1) = hpolot(k,j,2,1,1) + &
                osc(iex,2)*stpm(iex,1)*Ls* &
                hpol2si
            hpolot(k,j,2,1,2) = hpolot(k,j,2,1,2) + &
                osc(iex,2)*stpm(iex,4)*Ls* &
                hpol2si
            hpolot(k,j,2,1,3) = hpolot(k,j,2,1,3) + &
                osc(iex,2)*stpm(iex,5)*Ls* &
                hpol2si

            hpolot(k,j,2,2,1) = hpolot(k,j,2,2,1) + &
                osc(iex,2)*stpm(iex,4)*Ls* &
                hpol2si
            hpolot(k,j,2,2,2) = hpolot(k,j,2,2,2) + &
                osc(iex,2)*stpm(iex,2)*Ls* &
                hpol2si
            hpolot(k,j,2,2,3) = hpolot(k,j,2,2,3) + &
                osc(iex,2)*stpm(iex,6)*Ls* &
                hpol2si

            hpolot(k,j,2,3,1) = hpolot(k,j,2,3,1) + &
                osc(iex,2)*stpm(iex,5)*Ls* &
                hpol2si
            hpolot(k,j,2,3,2) = hpolot(k,j,2,3,2) + &
                osc(iex,2)*stpm(iex,6)*Ls* &
                hpol2si
            hpolot(k,j,2,3,3) = hpolot(k,j,2,3,3) + &
                osc(iex,2)*stpm(iex,3)*Ls* &
                hpol2si

            hpolot(k,j,3,1,1) = hpolot(k,j,3,1,1) + &
                osc(iex,3)*stpm(iex,1)*Ls* &
                hpol2si
            hpolot(k,j,3,1,2) = hpolot(k,j,3,1,2) + &
                osc(iex,3)*stpm(iex,4)*Ls* &
                hpol2si
            hpolot(k,j,3,1,3) = hpolot(k,j,3,1,3) + &
                osc(iex,3)*stpm(iex,5)*Ls* &
                hpol2si

            hpolot(k,j,3,2,1) = hpolot(k,j,3,2,1) + &
                osc(iex,3)*stpm(iex,4)*Ls* &
                hpol2si
            hpolot(k,j,3,2,2) = hpolot(k,j,3,2,2) + &
                osc(iex,3)*stpm(iex,2)*Ls* &
                hpol2si
            hpolot(k,j,3,2,3) = hpolot(k,j,3,2,3) + &
                osc(iex,3)*stpm(iex,6)*Ls* &
                hpol2si

            hpolot(k,j,3,3,1) = hpolot(k,j,3,3,1) + &
                osc(iex,3)*stpm(iex,5)*Ls* &
                hpol2si
            hpolot(k,j,3,3,2) = hpolot(k,j,3,3,2) + &
                osc(iex,3)*stpm(iex,6)*Ls* &
                hpol2si
            hpolot(k,j,3,3,3) = hpolot(k,j,3,3,3) + &
                osc(iex,3)*stpm(iex,3)*Ls* &
                hpol2si
          end do
        end do
      end do
    end if
    ! Calculate the transition hyperpolarizability for the A term for
    ! combination bands
    if (latermoff.eqv..false.) then
      do r=1,numcb
        do iex=1,nexci
          Ls = zzero
          do t=1,tstep
            dsum = zzero
            do i=1,nmodes
              dsum = dsum + &
                     (delta(iex,i)**2/two)* &
                     (one-exp(-z*freq(i)*x(t)))
            end do
            rsum = one
            do m=1,nmodes
              if (excnums(r,m) /= 0) then
                rsum = rsum* &
                    (delta(iex,m))**excnums(r,m)/ &
                    sqrt(two**excnums(r,m)* &
                    factorial(excnums(r,m)))* &
                    (one-exp(-z*freq(m)*x(t)))**excnums(r,m)
              end if
            end do
            Ls = Ls + z*w(t)*rsum*exp(x(t)* &
                 z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
          end do

          hpolcb(r,1,1,1) = hpolcb(r,1,1,1) + &
              osc(iex,1)*stpm(iex,1)*Ls* &
              hpol2si
          hpolcb(r,1,1,2) = hpolcb(r,1,1,2) + &
              osc(iex,1)*stpm(iex,4)*Ls* &
              hpol2si
          hpolcb(r,1,1,3) = hpolcb(r,1,1,3) + &
              osc(iex,1)*stpm(iex,5)*Ls* &
              hpol2si

          hpolcb(r,1,2,1) = hpolcb(r,1,2,1) + &
              osc(iex,1)*stpm(iex,4)*Ls* &
              hpol2si
          hpolcb(r,1,2,2) = hpolcb(r,1,2,2) + &
              osc(iex,1)*stpm(iex,2)*Ls* &
              hpol2si
          hpolcb(r,1,2,3) = hpolcb(r,1,2,3) + &
              osc(iex,1)*stpm(iex,6)*Ls* &
              hpol2si

          hpolcb(r,1,3,1) = hpolcb(r,1,3,1) + &
              osc(iex,1)*stpm(iex,5)*Ls* &
              hpol2si
          hpolcb(r,1,3,2) = hpolcb(r,1,3,2) + &
              osc(iex,1)*stpm(iex,6)*Ls* &
              hpol2si
          hpolcb(r,1,3,3) = hpolcb(r,1,3,3) + &
              osc(iex,1)*stpm(iex,3)*Ls* &
              hpol2si

          hpolcb(r,2,1,1) = hpolcb(r,2,1,1) + &
              osc(iex,2)*stpm(iex,1)*Ls* &
              hpol2si
          hpolcb(r,2,1,2) = hpolcb(r,2,1,2) + &
              osc(iex,2)*stpm(iex,4)*Ls* &
              hpol2si
          hpolcb(r,2,1,3) = hpolcb(r,2,1,3) + &
              osc(iex,2)*stpm(iex,5)*Ls* &
              hpol2si

          hpolcb(r,2,2,1) = hpolcb(r,2,2,1) + &
              osc(iex,2)*stpm(iex,4)*Ls* &
              hpol2si
          hpolcb(r,2,2,2) = hpolcb(r,2,2,2) + &
              osc(iex,2)*stpm(iex,2)*Ls* &
              hpol2si
          hpolcb(r,2,2,3) = hpolcb(r,2,2,3) + &
              osc(iex,2)*stpm(iex,6)*Ls* &
              hpol2si

          hpolcb(r,2,3,1) = hpolcb(r,2,3,1) + &
              osc(iex,2)*stpm(iex,5)*Ls* &
              hpol2si
          hpolcb(r,2,3,2) = hpolcb(r,2,3,2) + &
              osc(iex,2)*stpm(iex,6)*Ls* &
              hpol2si
          hpolcb(r,2,3,3) = hpolcb(r,2,3,3) + &
              osc(iex,2)*stpm(iex,3)*Ls* &
              hpol2si

          hpolcb(r,3,1,1) = hpolcb(r,3,1,1) + &
              osc(iex,3)*stpm(iex,1)*Ls* &
              hpol2si
          hpolcb(r,3,1,2) = hpolcb(r,3,1,2) + &
              osc(iex,3)*stpm(iex,4)*Ls* &
              hpol2si
          hpolcb(r,3,1,3) = hpolcb(r,3,1,3) + &
              osc(iex,3)*stpm(iex,5)*Ls* &
              hpol2si

          hpolcb(r,3,2,1) = hpolcb(r,3,2,1) + &
              osc(iex,3)*stpm(iex,4)*Ls* &
              hpol2si
          hpolcb(r,3,2,2) = hpolcb(r,3,2,2) + &
              osc(iex,3)*stpm(iex,2)*Ls* &
              hpol2si
          hpolcb(r,3,2,3) = hpolcb(r,3,2,3) + &
              osc(iex,3)*stpm(iex,6)*Ls* &
              hpol2si

          hpolcb(r,3,3,1) = hpolcb(r,3,3,1) + &
              osc(iex,3)*stpm(iex,5)*Ls* &
              hpol2si
          hpolcb(r,3,3,2) = hpolcb(r,3,3,2) + &
              osc(iex,3)*stpm(iex,6)*Ls* &
              hpol2si
          hpolcb(r,3,3,3) = hpolcb(r,3,3,3) + &
              osc(iex,3)*stpm(iex,3)*Ls* &
              hpol2si
        end do
      end do
    end if
  end if

end subroutine RHRSFCOverCB

subroutine RHRSB1Fund(x, w, freq, delta, stpm, tdipder, hpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb2, Ls_htb3, htlfreqdiff, &
                     ht_hpol, ht_rcrs, ht_rcrs_sum, lhtints)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep, r, k, n

  complex(kindr), intent(inout) :: hpol(nmodes,3,3,3)
  complex(kindr) :: ht_hpol(nmodes,nmodes,3,3,3)
  complex(kindr) :: rsum, rsum1, dsum, z
  complex(kindr) :: Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: stpm(nexci,6),tdipder(nexci,nmodes,3)
  real(kindr) :: htlfreqdiff(nexci)
  real(kindr) :: ht_rcrs(nmodes,nmodes)
  real(kindr) :: ht_rcrs_sum(nmodes)
  real(kindr) :: htexcist

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: lhtints

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B1 term hyperpolarizability for RHRS spectra (fundamentals)
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
!          lhtints - HT integrals are calculated?
!
! in/out : hpol - RHRS hyperpolarizability tensor
!
! =============================================================================

! Evaluates the the Herzberg-Teller terms for each normal mode. 
! Equations are based on the conversion from the vibronic theory 
! to the time-dependent formalism.
!
! Y. C. Chung and L. D. Ziegler.  J. Chem. Phys. 88, 7287 (1988).
!
! This part evaluates what I call B1 of the Herzberg-Teller terms, 
! which involves derivatives of the one photon transition dipole 
! moment between the initial and final electronic states.

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
          ! Integrals of the form: <(F-1)|I(t)>
          if (n == j) then
            ! Lower the final state of the current mode,
            ! so that F-1 = I = 0
            Ls_htb2(iex,n) = Ls_htb2(iex,n) + z*w(t)*exp(x(t)* &
                      z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
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
                      z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
          else
            ! Raise the final state of a different mode,
            ! so that I = 0 and F = 1.
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      z*w(t)*expsum(n)*rsum1*exp(x(t)* &
                      z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
          end if
        end do
      end do

      ! Add the A term to the B1 term to get the total 
      ! hyperpolarizability.

      do k = 1, nmodes

        hpol(j,1,1,1) = hpol(j,1,1,1) + &
          tdipder(iex,k,1)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,1,1,2) = hpol(j,1,1,2) + &
          tdipder(iex,k,1)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,1,1,3) = hpol(j,1,1,3) + &
          tdipder(iex,k,1)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        hpol(j,1,2,1) = hpol(j,1,2,1) + &
          tdipder(iex,k,1)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,1,2,2) = hpol(j,1,2,2) + &
          tdipder(iex,k,1)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,1,2,3) = hpol(j,1,2,3) + &
          tdipder(iex,k,1)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        hpol(j,1,3,1) = hpol(j,1,3,1) + &
          tdipder(iex,k,1)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,1,3,2) = hpol(j,1,3,2) + &
          tdipder(iex,k,1)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,1,3,3) = hpol(j,1,3,3) + &
          tdipder(iex,k,1)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        hpol(j,2,1,1) = hpol(j,2,1,1) + &
          tdipder(iex,k,2)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,2,1,2) = hpol(j,2,1,2) + &
          tdipder(iex,k,2)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,2,1,3) = hpol(j,2,1,3) + &
          tdipder(iex,k,2)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        hpol(j,2,2,1) = hpol(j,2,2,1) + &
          tdipder(iex,k,2)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,2,2,2) = hpol(j,2,2,2) + &
          tdipder(iex,k,2)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,2,2,3) = hpol(j,2,2,3) + &
          tdipder(iex,k,2)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        hpol(j,2,3,1) = hpol(j,2,3,1) + &
          tdipder(iex,k,2)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,2,3,2) = hpol(j,2,3,2) + &
          tdipder(iex,k,2)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,2,3,3) = hpol(j,2,3,3) + &
          tdipder(iex,k,2)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        hpol(j,3,1,1) = hpol(j,3,1,1) + &
          tdipder(iex,k,3)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,3,1,2) = hpol(j,3,1,2) + &
          tdipder(iex,k,3)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,3,1,3) = hpol(j,3,1,3) + &
          tdipder(iex,k,3)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        hpol(j,3,2,1) = hpol(j,3,2,1) + &
          tdipder(iex,k,3)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,3,2,2) = hpol(j,3,2,2) + &
          tdipder(iex,k,3)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,3,2,3) = hpol(j,3,2,3) + &
          tdipder(iex,k,3)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 

        hpol(j,3,3,1) = hpol(j,3,3,1) + &
          tdipder(iex,k,3)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,3,3,2) = hpol(j,3,3,2) + &
          tdipder(iex,k,3)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
        hpol(j,3,3,3) = hpol(j,3,3,3) + &
          tdipder(iex,k,3)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si 
      end do

      ! Evaluate the coupling between modes, if requested.
      if (lhtints) then
        ! Find the excited state nearest in energy to the incident
        ! frequency.
        do i = 1, nexci
          htlfreqdiff(i) = lfreq - freq0(i)
          htlfreqdiff(i) = abs(htlfreqdiff(i))
        end do
        htexcist = minval(htlfreqdiff)
        call B1HTInts(Ls_htb2, Ls_htb3, freq, stpm, tdipder, &
                      freq0, htexcist, ht_hpol, ht_rcrs, &
                      ht_rcrs_sum, delta, lfreq, nexci, nmodes, &
                      3, 6, iex, j)
      end if
    end do

    ! Test
    !write(6,*) 'A + B1 term hyperpolarizability'
    !write(6,*) freq(j)
    !!------------------------------
    !write(6,*) 'xxx', hpol(j,1,1,1)
    !write(6,*) 'xxy', hpol(j,1,1,2)
    !write(6,*) 'xxz', hpol(j,1,1,3)
    !write(6,*) 'xyx', hpol(j,1,2,1)
    !write(6,*) 'xyy', hpol(j,1,2,2)
    !write(6,*) 'xyz', hpol(j,1,2,3)
    !write(6,*) 'xzx', hpol(j,1,3,1)
    !write(6,*) 'xzy', hpol(j,1,3,2)
    !write(6,*) 'xzz', hpol(j,1,3,3)
    !!------------------------------
    !write(6,*) 'yxx', hpol(j,2,1,1)
    !write(6,*) 'yxy', hpol(j,2,1,2)
    !write(6,*) 'yxz', hpol(j,2,1,3)
    !write(6,*) 'yyx', hpol(j,2,2,1)
    !write(6,*) 'yyy', hpol(j,2,2,2)
    !write(6,*) 'yyz', hpol(j,2,2,3)
    !write(6,*) 'yzx', hpol(j,2,3,1)
    !write(6,*) 'yzy', hpol(j,2,3,2)
    !write(6,*) 'yzz', hpol(j,2,3,3)
    !!------------------------------
    !write(6,*) 'zxx', hpol(j,3,1,1)
    !write(6,*) 'zxy', hpol(j,3,1,2)
    !write(6,*) 'zxz', hpol(j,3,1,3)
    !write(6,*) 'zyx', hpol(j,3,2,1)
    !write(6,*) 'zyy', hpol(j,3,2,2)
    !write(6,*) 'zyz', hpol(j,3,2,3)
    !write(6,*) 'zzx', hpol(j,3,3,1)
    !write(6,*) 'zzy', hpol(j,3,3,2)
    !write(6,*) 'zzz', hpol(j,3,3,3)
    !write(6,*) '------------------------------'

  end do

end subroutine RHRSB1Fund

subroutine RHRSB2Fund(x, w, freq, delta, osc, stpmder, hpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb, htlfreqdiff, &
                     ht_hpol, ht_rcrs, ht_rcrs_sum, lhtints)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep, r, k, n

  complex(kindr), intent(inout) :: hpol(nmodes,3,3,3)
  complex(kindr) :: ht_hpol(nmodes,nmodes,3,3,3)
  complex(kindr) :: rsum, rsum1, dsum, z
  complex(kindr) :: Ls_htb(nexci,nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3),stpmder(nexci,nmodes,6)
  real(kindr) :: htlfreqdiff(nexci)
  real(kindr) :: ht_rcrs(nmodes,nmodes)
  real(kindr) :: ht_rcrs_sum(nmodes)
  real(kindr) :: htexcist

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: lhtints

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B2 term hyperpolarizability for RHRS spectra (fundamentals)
!
! input  : x - abcissa points for numerical integration 
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          osc - transition dipole moment components
!          stpmder - derivative of two-photon transition moment components
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!          lhtints - HT integrals are calculated?
!
! in/out : hpol - RHRS hyperpolarizability tensor
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
            ! Raise the initial state of the current mode,
            ! so that I+1 = F = 1, while leaving every other
            ! mode in their ground states. 
            Ls_htb(iex,n) = Ls_htb(iex,n) + z*w(t)* &
                     (one-delta(iex,n)**2+ &
                     delta(iex,n)**2*cos(freq(n)*x(t)))*exp(x(t)* &
                     z*(lfreq-freq0(iex)-freq(n)) - g(iex,t) - dsum)
          else
            ! Raise the initial state of a different mode,
            ! so I+1 = 1 and F = 0.  Here, there's a product
            ! accounting for multiple modes in an excited state.
            Ls_htb(iex,n) = Ls_htb(iex,n) + z*w(t)*rsum1*expsum(n)* &
                     exp(x(t)* &
                     z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
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

        hpol(j,1,1,1) = hpol(j,1,1,1) + &
          osc(iex,1)*stpmder(iex,k,1)*Ls_htb(iex,k)*hpol2si
        hpol(j,1,1,2) = hpol(j,1,1,2) + &
          osc(iex,1)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
        hpol(j,1,1,3) = hpol(j,1,1,3) + &
          osc(iex,1)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si

        hpol(j,1,2,1) = hpol(j,1,2,1) + &
          osc(iex,1)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
        hpol(j,1,2,2) = hpol(j,1,2,2) + &
          osc(iex,1)*stpmder(iex,k,2)*Ls_htb(iex,k)*hpol2si
        hpol(j,1,2,3) = hpol(j,1,2,3) + &
          osc(iex,1)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si

        hpol(j,1,3,1) = hpol(j,1,3,1) + &
          osc(iex,1)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si
        hpol(j,1,3,2) = hpol(j,1,3,2) + &
          osc(iex,1)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si
        hpol(j,1,3,3) = hpol(j,1,3,3) + &
          osc(iex,1)*stpmder(iex,k,3)*Ls_htb(iex,k)*hpol2si

        hpol(j,2,1,1) = hpol(j,2,1,1) + &
          osc(iex,2)*stpmder(iex,k,1)*Ls_htb(iex,k)*hpol2si
        hpol(j,2,1,2) = hpol(j,2,1,2) + &
          osc(iex,2)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
        hpol(j,2,1,3) = hpol(j,2,1,3) + &
          osc(iex,2)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si

        hpol(j,2,2,1) = hpol(j,2,2,1) + &
          osc(iex,2)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
        hpol(j,2,2,2) = hpol(j,2,2,2) + &
          osc(iex,2)*stpmder(iex,k,2)*Ls_htb(iex,k)*hpol2si
        hpol(j,2,2,3) = hpol(j,2,2,3) + &
          osc(iex,2)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si

        hpol(j,2,3,1) = hpol(j,2,3,1) + &
          osc(iex,2)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si
        hpol(j,2,3,2) = hpol(j,2,3,2) + &
          osc(iex,2)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si
        hpol(j,2,3,3) = hpol(j,2,3,3) + &
          osc(iex,2)*stpmder(iex,k,3)*Ls_htb(iex,k)*hpol2si

        hpol(j,3,1,1) = hpol(j,3,1,1) + &
          osc(iex,3)*stpmder(iex,k,1)*Ls_htb(iex,k)*hpol2si
        hpol(j,3,1,2) = hpol(j,3,1,2) + &
          osc(iex,3)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
        hpol(j,3,1,3) = hpol(j,3,1,3) + &
          osc(iex,3)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si

        hpol(j,3,2,1) = hpol(j,3,2,1) + &
          osc(iex,3)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
        hpol(j,3,2,2) = hpol(j,3,2,2) + &
          osc(iex,3)*stpmder(iex,k,2)*Ls_htb(iex,k)*hpol2si
        hpol(j,3,2,3) = hpol(j,3,2,3) + &
          osc(iex,3)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si

        hpol(j,3,3,1) = hpol(j,3,3,1) + &
          osc(iex,3)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si
        hpol(j,3,3,2) = hpol(j,3,3,2) + &
          osc(iex,3)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si
        hpol(j,3,3,3) = hpol(j,3,3,3) + &
          osc(iex,3)*stpmder(iex,k,3)*Ls_htb(iex,k)*hpol2si
      end do

      ! Evaluate the coupling between modes, if requested.
      if (lhtints) then
        ! Find the excited state nearest in energy to the incident
        ! frequency.
        do i = 1, nexci
          htlfreqdiff(i) = lfreq - freq0(i)
          htlfreqdiff(i) = abs(htlfreqdiff(i))
        end do
        htexcist = minval(htlfreqdiff)
        call B2HTInts(Ls_htb, freq, stpmder, osc, &
                      freq0, htexcist, ht_hpol, ht_rcrs, &
                      ht_rcrs_sum, delta, lfreq, nexci, nmodes, &
                      3, 6, iex, j)
      end if
    end do

    !! Test
    !write(6,*) 'A + B term hyperpolarizability'
    !write(6,*) freq(j)
    !!------------------------------
    !write(6,*) 'xxx', hpol(j,1,1,1)
    !write(6,*) 'xxy', hpol(j,1,1,2)
    !write(6,*) 'xxz', hpol(j,1,1,3)
    !write(6,*) 'xyx', hpol(j,1,2,1)
    !write(6,*) 'xyy', hpol(j,1,2,2)
    !write(6,*) 'xyz', hpol(j,1,2,3)
    !write(6,*) 'xzx', hpol(j,1,3,1)
    !write(6,*) 'xzy', hpol(j,1,3,2)
    !write(6,*) 'xzz', hpol(j,1,3,3)
    !!------------------------------
    !write(6,*) 'yxx', hpol(j,2,1,1)
    !write(6,*) 'yxy', hpol(j,2,1,2)
    !write(6,*) 'yxz', hpol(j,2,1,3)
    !write(6,*) 'yyx', hpol(j,2,2,1)
    !write(6,*) 'yyy', hpol(j,2,2,2)
    !write(6,*) 'yyz', hpol(j,2,2,3)
    !write(6,*) 'yzx', hpol(j,2,3,1)
    !write(6,*) 'yzy', hpol(j,2,3,2)
    !write(6,*) 'yzz', hpol(j,2,3,3)
    !!------------------------------
    !write(6,*) 'zxx', hpol(j,3,1,1)
    !write(6,*) 'zxy', hpol(j,3,1,2)
    !write(6,*) 'zxz', hpol(j,3,1,3)
    !write(6,*) 'zyx', hpol(j,3,2,1)
    !write(6,*) 'zyy', hpol(j,3,2,2)
    !write(6,*) 'zyz', hpol(j,3,2,3)
    !write(6,*) 'zzx', hpol(j,3,3,1)
    !write(6,*) 'zzy', hpol(j,3,3,2)
    !write(6,*) 'zzz', hpol(j,3,3,3)
    !write(6,*) '------------------------------'

  end do

end subroutine RHRSB2Fund

subroutine RHRSB1OverCB(x, w, freq, delta, stpm, tdipder, hpolot, &
               hpolcb, freq0, lfreq, g, nexci, nmodes, tstep, &
               excnum, numcb, excnums, factorial, rootn, &
               expsum, Ls_htb2, Ls_htb3, sumr, sumr1)

  use constants

  implicit none

  integer :: nexci, nmodes, excnum, numcb
  integer :: j, iex, t, i, tstep, r, k, n, p, m
  integer :: excnums(numcb,nmodes)

  complex(kindr), intent(inout) :: hpolot(excnum,nmodes,3,3,3)
  complex(kindr), intent(inout) :: hpolcb(numcb,3,3,3)
  complex(kindr) :: rsum, dsum, z
  complex(kindr) :: Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr) :: sumr(nmodes), sumr1(nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr) :: rsum1, rsum2
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: stpm(nexci,6),tdipder(nexci,nmodes,3)
  real(kindr) :: factorial(16), rootn(12)

  complex(kindr) :: zzero=(0d0,0d0)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B1 term hyperpolarizability for RHRS spectra (overtones and
!          combination bands)
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
!          excnum - excitation level
!          numcb - number of important combination bands
!          excnums - array of each mode's state
!          factorial - array of factorial values
!          rootn - array of square roots of double precision numbers
!
! in/out : hpolot - RHRS hyperpolarizability tensor for overtones
!          hpolcb - RHRS hyperpolarizability tensor for combination bands
!
! =============================================================================

  ! Calculate the B1 term overtones
  do r=3,excnum+1
    do j = 1, nmodes
      do iex=1,nexci
        Ls_htb2(1:nexci,1:nmodes) = zzero
        Ls_htb3(1:nexci,1:nmodes) = zzero
        do t=1,tstep
          dsum = zzero
          do i=1,nmodes
            dsum = dsum + (delta(iex,i)**2/two)* &
                   (one-exp(-z*freq(i)*x(t)))
          end do
          ! 3 rsums are: current mode's raised overtone (rsum, 
          ! at excitation number r), the current mode's 
          ! overtone (rsum1, at excitation number r-1), and the
          ! current mode's lowered overtone (rsum2, at excitation
          ! number r-2)
          rsum = (delta(iex,j))**r/ &
                 (sqrt(two**r*factorial(r)))* &
                 (one-exp(-z*freq(j)*x(t)))**r
          rsum1 = (delta(iex,j))**(r-1)/ &
                  (sqrt(two**(r-1)*factorial(r-1)))* &
                  (one-exp(-z*freq(j)*x(t)))**(r-1)
          rsum2 = (delta(iex,j))**(r-2)/ &
                  (sqrt(two**(r-2)*factorial(r-2)))* &
                  (one-exp(-z*freq(j)*x(t)))**(r-2)
          ! Calculate the contribution from other normal modes.  This
          ! is either raised or has nothing performed, if the mode is
          ! the current mode.
          do m=1,nmodes
            if (m /= j) then
              expsum(m) = delta(iex,m)/sqrt(two)* &
                          (one-exp(-z*freq(r)*x(t)))
            else
              expsum(m) = one
            end if
          end do
          ! Calculate lineshape functions
          do n=1,nmodes
            ! Integrals of the form: <(F-1)|I(t)>
            if (n == j) then
              ! Lower the final state of the current mode,
              ! so that F-1 = excnum-1, I = 0
              Ls_htb2(iex,n) = Ls_htb2(iex,n) + &
                        z*rootn(r-1)*rsum2*w(t)*exp(x(t)* &
                        z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
            else
              ! Can't lower a harmonic oscillator past its ground
              ! state.
              Ls_htb2(iex,n) = zzero
            end if
            ! Integrals of the form: <(F+1)|I(t)> 
            if (n == j) then
              ! Raise the final state of the current mode,
              ! so that I = 0 and F+1 = excnum+1
              Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                        z*w(t)*rootn(r+1)*rsum*exp(x(t)* &
                        z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
            else
              ! Raise the final state of a different mode,
              ! so that I = 0 and F = 1.
              Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                        z*w(t)*expsum(n)*rsum1*exp(x(t)* &
                        z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
            end if
          end do
        end do

        ! Add the A term to the B1 term to get the total overtone 
        ! hyperpolarizability.

        do k = 1, nmodes

          hpolot(r-1,j,1,1,1) = hpolot(r-1,j,1,1,1) + &
            tdipder(iex,k,1)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,1,1,2) = hpolot(r-1,j,1,1,2) + &
            tdipder(iex,k,1)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,1,1,3) = hpolot(r-1,j,1,1,3) + &
            tdipder(iex,k,1)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si

          hpolot(r-1,j,1,2,1) = hpolot(r-1,j,1,2,1) + &
            tdipder(iex,k,1)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,1,2,2) = hpolot(r-1,j,1,2,2) + &
            tdipder(iex,k,1)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,1,2,3) = hpolot(r-1,j,1,2,3) + &
            tdipder(iex,k,1)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si

          hpolot(r-1,j,1,3,1) = hpolot(r-1,j,1,3,1) + &
            tdipder(iex,k,1)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,1,3,2) = hpolot(r-1,j,1,3,2) + &
            tdipder(iex,k,1)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,1,3,3) = hpolot(r-1,j,1,3,3) + &
            tdipder(iex,k,1)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si

          hpolot(r-1,j,2,1,1) = hpolot(r-1,j,2,1,1) + &
            tdipder(iex,k,2)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,2,1,2) = hpolot(r-1,j,2,1,2) + &
            tdipder(iex,k,2)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,2,1,3) = hpolot(r-1,j,2,1,3) + &
            tdipder(iex,k,2)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si

          hpolot(r-1,j,2,2,1) = hpolot(r-1,j,2,2,1) + &
            tdipder(iex,k,2)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,2,2,2) = hpolot(r-1,j,2,2,2) + &
            tdipder(iex,k,2)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,2,2,3) = hpolot(r-1,j,2,2,3) + &
            tdipder(iex,k,2)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si

          hpolot(r-1,j,2,3,1) = hpolot(r-1,j,2,3,1) + &
            tdipder(iex,k,2)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,2,3,2) = hpolot(r-1,j,2,3,2) + &
            tdipder(iex,k,2)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,2,3,3) = hpolot(r-1,j,2,3,3) + &
            tdipder(iex,k,2)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si

          hpolot(r-1,j,3,1,1) = hpolot(r-1,j,3,1,1) + &
            tdipder(iex,k,3)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,3,1,2) = hpolot(r-1,j,3,1,2) + &
            tdipder(iex,k,3)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,3,1,3) = hpolot(r-1,j,3,1,3) + &
            tdipder(iex,k,3)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si

          hpolot(r-1,j,3,2,1) = hpolot(r-1,j,3,2,1) + &
            tdipder(iex,k,3)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,3,2,2) = hpolot(r-1,j,3,2,2) + &
            tdipder(iex,k,3)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,3,2,3) = hpolot(r-1,j,3,2,3) + &
            tdipder(iex,k,3)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si

          hpolot(r-1,j,3,3,1) = hpolot(r-1,j,3,3,1) + &
            tdipder(iex,k,3)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,3,3,2) = hpolot(r-1,j,3,3,2) + &
            tdipder(iex,k,3)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
          hpolot(r-1,j,3,3,3) = hpolot(r-1,j,3,3,3) + &
            tdipder(iex,k,3)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
            Ls_htb3(iex,k))*hpol2si
        end do
      end do
    end do
  end do
  ! Calculate the B1 term combination bands
  do r=1,numcb
    do iex=1,nexci
      Ls_htb2(1:nexci,1:nmodes) = zzero
      Ls_htb3(1:nexci,1:nmodes) = zzero
      do t=1,tstep
        dsum = zzero
        do i=1,nmodes
          dsum = dsum + (delta(iex,i)**2/two)* &
                 (one-exp(-z*freq(i)*x(t)))
        end do
        ! Results for raising and lowering operators.
        ! sumr -  contribution after raising the final state
        ! sumr1 - contribution after lowering the final state
        sumr(1:nmodes) = one
        sumr1(1:nmodes) = one
        j = 1
        do p=1,nmodes
          ! Contribution from <(F-1)|I(t)>
          do i=1,nmodes
            if (i==j) then
              if (excnums(r,j) /= 0) then
                sumr1(j) = sumr1(j)* &
                    (delta(iex,i))**(excnums(r,i)-1)/ &
                    sqrt(two**(excnums(r,i)-1)* &
                    factorial(excnums(r,i)-1))* &
                    (one-exp(-z*freq(i)*x(t)))**(excnums(r,i)-1)
              else
                sumr1(j) = sumr1(j)*zzero
              end if
            else
              sumr1(j) = sumr1(j)* &
                  (delta(iex,i))**(excnums(r,i))/ &
                  sqrt(two**(excnums(r,i))* &
                  factorial(excnums(r,i)))* &
                  (one-exp(-z*freq(i)*x(t)))**(excnums(r,i))
            end if
          end do
          ! Contribution from <(F+1)|I(t)>
          do i=1,nmodes
            if (i==j) then
              if (excnums(r,j) /= 0) then
                sumr(j) = sumr(j)* &
                    (delta(iex,i))**(excnums(r,i)+1)/ &
                    sqrt(two**(excnums(r,i)+1)* &
                    factorial(excnums(r,i)+1))* &
                    (one-exp(-z*freq(i)*x(t)))**(excnums(r,i)+1)
              else
                sumr(j) = sumr(j)* &
                  (delta(iex,i))**(excnums(r,i)+1)/ &
                  sqrt(two**(excnums(r,i)+1)* &
                  factorial(excnums(r,i)+1))* &
                  (one-exp(-z*freq(i)*x(t)))**(excnums(r,i)+1)
              end if
            else
              sumr(j) = sumr(j)* &
                  (delta(iex,i))**(excnums(r,i))/ &
                  sqrt(two**(excnums(r,i))* &
                  factorial(excnums(r,i)))* &
                  (one-exp(-z*freq(i)*x(t)))**(excnums(r,i))
            end if
          end do
          j = j+1
        end do
        ! Calculate lineshapes
        do i=1,nmodes
          ! Integrals of the form: <(F-1)|I(t)>
          Ls_htb2(iex,i) = Ls_htb2(iex,i) + &
                      z*rootn(excnums(r,i))*sumr1(i)*w(t)* &
                      exp(x(t)* &
                      z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
          ! Integrals of the form: <(F+1)|I(t)>
          Ls_htb3(iex,i) = Ls_htb3(iex,i) + &
                      z*rootn(excnums(r,i)+1)*sumr(i)*w(t)* &
                      exp(x(t)* &
                      z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
        end do
      end do

      ! Add the A term to the B1 term to get the total  
      ! combination band hyperpolarizability.

      do k = 1, nmodes

        hpolcb(r,1,1,1) = hpolcb(r,1,1,1) + &
          tdipder(iex,k,1)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,1,1,2) = hpolcb(r,1,1,2) + &
          tdipder(iex,k,1)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,1,1,3) = hpolcb(r,1,1,3) + &
          tdipder(iex,k,1)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si

        hpolcb(r,1,2,1) = hpolcb(r,1,2,1) + &
          tdipder(iex,k,1)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,1,2,2) = hpolcb(r,1,2,2) + &
          tdipder(iex,k,1)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,1,2,3) = hpolcb(r,1,2,3) + &
          tdipder(iex,k,1)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si

        hpolcb(r,1,3,1) = hpolcb(r,1,3,1) + &
          tdipder(iex,k,1)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,1,3,2) = hpolcb(r,1,3,2) + &
          tdipder(iex,k,1)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,1,3,3) = hpolcb(r,1,3,3) + &
          tdipder(iex,k,1)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si

        hpolcb(r,2,1,1) = hpolcb(r,2,1,1) + &
          tdipder(iex,k,2)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,2,1,2) = hpolcb(r,2,1,2) + &
          tdipder(iex,k,2)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,2,1,3) = hpolcb(r,2,1,3) + &
          tdipder(iex,k,2)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si

        hpolcb(r,2,2,1) = hpolcb(r,2,2,1) + &
          tdipder(iex,k,2)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,2,2,2) = hpolcb(r,2,2,2) + &
          tdipder(iex,k,2)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,2,2,3) = hpolcb(r,2,2,3) + &
          tdipder(iex,k,2)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si

        hpolcb(r,2,3,1) = hpolcb(r,2,3,1) + &
          tdipder(iex,k,2)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,2,3,2) = hpolcb(r,2,3,2) + &
          tdipder(iex,k,2)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,2,3,3) = hpolcb(r,2,3,3) + &
          tdipder(iex,k,2)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si

        hpolcb(r,3,1,1) = hpolcb(r,3,1,1) + &
          tdipder(iex,k,3)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,3,1,2) = hpolcb(r,3,1,2) + &
          tdipder(iex,k,3)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,3,1,3) = hpolcb(r,3,1,3) + &
          tdipder(iex,k,3)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si

        hpolcb(r,3,2,1) = hpolcb(r,3,2,1) + &
          tdipder(iex,k,3)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,3,2,2) = hpolcb(r,3,2,2) + &
          tdipder(iex,k,3)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,3,2,3) = hpolcb(r,3,2,3) + &
          tdipder(iex,k,3)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si

        hpolcb(r,3,3,1) = hpolcb(r,3,3,1) + &
          tdipder(iex,k,3)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,3,3,2) = hpolcb(r,3,3,2) + &
          tdipder(iex,k,3)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
        hpolcb(r,3,3,3) = hpolcb(r,3,3,3) + &
          tdipder(iex,k,3)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
          Ls_htb3(iex,k))*hpol2si
      end do
    end do
  end do

end subroutine RHRSB1OverCB

subroutine RHRSB2OverCB(x, w, freq, delta, osc, stpmder, hpolot, &
               hpolcb, freq0, lfreq, g, nexci, nmodes, tstep, &
               excnum, numcb, excnums, factorial, num, &
               expsum, Ls_htb, sumr2)

  use constants

  implicit none

  integer :: nexci, nmodes, excnum, numcb
  integer :: j, iex, t, i, tstep, r, k, n, p, m
  integer :: excnums(numcb,nmodes)

  complex(kindr), intent(inout) :: hpolot(excnum,nmodes,3,3,3)
  complex(kindr), intent(inout) :: hpolcb(numcb,3,3,3)
  complex(kindr) :: rsum, dsum, z
  complex(kindr) :: Ls_htb(nexci,nmodes)
  complex(kindr) :: sumr2(nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr) :: rsum1, rsum3
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3),stpmder(nexci,nmodes,6)
  real(kindr) :: factorial(16), num(11)

  complex(kindr) :: zzero=(0d0,0d0)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B2 term hyperpolarizability for RHRS spectra (overtones and
!          combination bands)
!
! input  : x - abcissa points for numerical integration 
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          osc - transition dipole moment components
!          stpmder - derivative of two-photon transition moment components
!          freq0 - vertical excitation energy
!          lfreq - laser frequency
!          g - lifetime function g(t)
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!          excnum - excitation level
!          numcb - number of important combination bands
!          excnums - array of each mode's state
!          factorial - array of factorial values
!          num - array of double precision numbers
!
! in/out : hpolot - RHRS hyperpolarizability tensor for overtones
!          hpolcb - RHRS hyperpolarizability tensor for combination bands
!
! =============================================================================

  ! Calculate B2 term overtones
  do r=2,excnum
    do j = 1, nmodes
      do iex=1,nexci
        Ls_htb(1:nexci,1:nmodes) = zzero
        do t=1,tstep
          dsum = zzero
          do i=1,nmodes
            dsum = dsum + (delta(iex,i)**2/two)* &
                   (one-exp(-z*freq(i)*x(t)))
          end do
          ! rsum1 is evaluated by decreasing the excitation
          ! number (r) by 1, for use with the integrals below
          rsum = (delta(iex,j))**r/ &
                 (sqrt(two**r*factorial(r)))* &
                 (one-exp(-z*freq(j)*x(t)))**r
          rsum1 = (delta(iex,j))**(r-1)/ &
                  (sqrt(two**(r-1)*factorial(r-1)))* &
                  (one-exp(-z*freq(j)*x(t)))**(r-1)
          ! The pattern for this is:
          ! <n|1(t)> = delta^(n-1)/sqrt(2^(n-1)n!)* &
          !            (n-delta^2+delta^2*cos(omega*t))* & 
          !            e^(-s(1-e^(-i*omega*t)-3*i*omega*t/2)
          rsum3 = (delta(iex,j))**(r-1)/ &
                  (sqrt(two**(r-1)*factorial(r)))* &
                  (num(r) - delta(iex,j)**2 + &
                  delta(iex,j)**2*cos(freq(j)*x(t)))* &
                  (one-exp(-z*freq(j)*x(t)))**(r-1)
          ! Calculate the contribution from other normal modes.  
          ! This is either raised or has nothing performed if 
          ! the mode is the current mode.
          do m=1,nmodes
            if (m /= j) then
              expsum(m) = delta(iex,m)/sqrt(two)* &
                          (one-exp(-z*freq(m)*x(t)))
            else
              expsum(m) = one
            end if
          end do
          ! Calculate lineshape functions
          do n=1,nmodes
            ! Integrals of the form: <F|(I+1)(t)>
            if (n == j) then
              ! Raise the initial state of the current mode,
              ! so that I+1 = 1, F = excnum.
              Ls_htb(iex,n) = Ls_htb(iex,n) + z*rsum3*w(t)* &
                       exp(x(t)* &
                       z*(lfreq-freq0(iex)-freq(n)) - g(iex,t) - dsum)
            else
              ! Raise the initial state of a different mode,
              ! so I+1 = 1 and F = 0.  Here, there's a product
              ! accounting for multiple modes in an excited state.
              Ls_htb(iex,n) = Ls_htb(iex,n) + &
                              z*w(t)*rsum*expsum(n)* &
                              exp(x(t)* &
                              z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
            end if
          end do
        end do

        do k = 1, nmodes

          hpolot(r,j,1,1,1) = hpolot(r,j,1,1,1) + &
            osc(iex,1)*stpmder(iex,k,1)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,1,1,2) = hpolot(r,j,1,1,2) + &
            osc(iex,1)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,1,1,3) = hpolot(r,j,1,1,3) + &
            osc(iex,1)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
            hpol2si

          hpolot(r,j,1,2,1) = hpolot(r,j,1,2,1) + &
            osc(iex,1)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,1,2,2) = hpolot(r,j,1,2,2) + &
            osc(iex,1)*stpmder(iex,k,2)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,1,2,3) = hpolot(r,j,1,2,3) + &
            osc(iex,1)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
            hpol2si

          hpolot(r,j,1,3,1) = hpolot(r,j,1,3,1) + &
            osc(iex,1)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,1,3,2) = hpolot(r,j,1,3,2) + &
            osc(iex,1)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,1,3,3) = hpolot(r,j,1,3,3) + &
            osc(iex,1)*stpmder(iex,k,3)*Ls_htb(iex,k)* &
            hpol2si

          hpolot(r,j,2,1,1) = hpolot(r,j,2,1,1) + &
            osc(iex,2)*stpmder(iex,k,1)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,2,1,2) = hpolot(r,j,2,1,2) + &
            osc(iex,2)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,2,1,3) = hpolot(r,j,2,1,3) + &
            osc(iex,2)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
            hpol2si

          hpolot(r,j,2,2,1) = hpolot(r,j,2,2,1) + &
            osc(iex,2)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,2,2,2) = hpolot(r,j,2,2,2) + &
            osc(iex,2)*stpmder(iex,k,2)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,2,2,3) = hpolot(r,j,2,2,3) + &
            osc(iex,2)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
            hpol2si

          hpolot(r,j,2,3,1) = hpolot(r,j,2,3,1) + &
            osc(iex,2)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,2,3,2) = hpolot(r,j,2,3,2) + &
            osc(iex,2)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,2,3,3) = hpolot(r,j,2,3,3) + &
            osc(iex,2)*stpmder(iex,k,3)*Ls_htb(iex,k)* &
            hpol2si

          hpolot(r,j,3,1,1) = hpolot(r,j,3,1,1) + &
            osc(iex,3)*stpmder(iex,k,1)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,3,1,2) = hpolot(r,j,3,1,2) + &
            osc(iex,3)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,3,1,3) = hpolot(r,j,3,1,3) + &
            osc(iex,3)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
            hpol2si

          hpolot(r,j,3,2,1) = hpolot(r,j,3,2,1) + &
            osc(iex,3)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,3,2,2) = hpolot(r,j,3,2,2) + &
            osc(iex,3)*stpmder(iex,k,2)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,3,2,3) = hpolot(r,j,3,2,3) + &
            osc(iex,3)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
            hpol2si

          hpolot(r,j,3,3,1) = hpolot(r,j,3,3,1) + &
            osc(iex,3)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,3,3,2) = hpolot(r,j,3,3,2) + &
            osc(iex,3)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
            hpol2si
          hpolot(r,j,3,3,3) = hpolot(r,j,3,3,3) + &
            osc(iex,3)*stpmder(iex,k,3)*Ls_htb(iex,k)* &
            hpol2si
        end do
      end do
    end do
  end do

  ! Combination band contribution to the B2 term.

  do r=1,numcb
    do iex=1,nexci
      Ls_htb(1:nexci,1:nmodes) = zzero
      do t=1,tstep
        dsum = zzero
        do i=1,nmodes
          dsum = dsum + (delta(iex,i)**2/two)* &
                 (one-exp(-z*freq(i)*x(t)))
        end do
        ! Results for raising and lowering operators.
        ! sumr2 - contribution after raising the initial state
        sumr2(1:nmodes) = one
        j = 1
        do p=1,nmodes
          ! Contribution from <F|(I+1)(t)>
          do i=1,nmodes
            if (i==j) then
              if (excnums(r,j) /= 0) then
                sumr2(j) = sumr2(j)* &
                    (delta(iex,i))**(excnums(r,i)-1)/ &
                    sqrt(two**(excnums(r,i)-1)* &
                    factorial(excnums(r,i)-1))* &
                    (num(excnums(r,i)) - delta(iex,i)**2 + &
                    delta(iex,i)**2*cos(freq(i)*x(t)))* &
                    (one-exp(-z*freq(i)*x(t)))**(excnums(r,i)-1)
              else
                sumr2(j) = sumr2(j)* &
                    (delta(iex,i))**(excnums(r,i)+1)/ &
                    sqrt(two**(excnums(r,i)+1)* &
                    factorial(excnums(r,i)+1))* &
                    (one-exp(-z*freq(i)*x(t)))**(excnums(r,i)+1)
              end if
            else
              sumr2(j) = sumr2(j)* &
                  (delta(iex,i))**(excnums(r,i))/ &
                  sqrt(two**(excnums(r,i))* &
                  factorial(excnums(r,i)))* &
                  (one-exp(-z*freq(i)*x(t)))**(excnums(r,i))
            end if
          end do
          j = j+1
        end do
        ! Calculate lineshapes
        do i=1,nmodes
          ! Integrals of the form: <F|(I+1)(t)>.  In effect,
          ! the current mode is lowered by raising its initial
          ! state.
          Ls_htb(iex,i) = Ls_htb(iex,i) + z*sumr2(i)*w(t)* &
                      exp(x(t)* &
                      z*(lfreq-freq0(iex)-freq(i)) - g(iex,t) - dsum)
        end do
      end do

      ! Add the A term to the B2 term to get the total  
      ! combination band hyperpolarizability.

      do k = 1, nmodes

        hpolcb(r,1,1,1) = hpolcb(r,1,1,1) + &
          osc(iex,1)*stpmder(iex,k,1)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,1,1,2) = hpolcb(r,1,1,2) + &
          osc(iex,1)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,1,1,3) = hpolcb(r,1,1,3) + &
          osc(iex,1)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
          hpol2si

        hpolcb(r,1,2,1) = hpolcb(r,1,2,1) + &
          osc(iex,1)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,1,2,2) = hpolcb(r,1,2,2) + &
          osc(iex,1)*stpmder(iex,k,2)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,1,2,3) = hpolcb(r,1,2,3) + &
          osc(iex,1)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
          hpol2si

        hpolcb(r,1,3,1) = hpolcb(r,1,3,1) + &
          osc(iex,1)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,1,3,2) = hpolcb(r,1,3,2) + &
          osc(iex,1)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,1,3,3) = hpolcb(r,1,3,3) + &
          osc(iex,1)*stpmder(iex,k,3)*Ls_htb(iex,k)* &
          hpol2si

        hpolcb(r,2,1,1) = hpolcb(r,2,1,1) + &
          osc(iex,2)*stpmder(iex,k,1)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,2,1,2) = hpolcb(r,2,1,2) + &
          osc(iex,2)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,2,1,3) = hpolcb(r,2,1,3) + &
          osc(iex,2)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
          hpol2si

        hpolcb(r,2,2,1) = hpolcb(r,2,2,1) + &
          osc(iex,2)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,2,2,2) = hpolcb(r,2,2,2) + &
          osc(iex,2)*stpmder(iex,k,2)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,2,2,3) = hpolcb(r,2,2,3) + &
          osc(iex,2)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
          hpol2si

        hpolcb(r,2,3,1) = hpolcb(r,2,3,1) + &
          osc(iex,2)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,2,3,2) = hpolcb(r,2,3,2) + &
          osc(iex,2)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,2,3,3) = hpolcb(r,2,3,3) + &
          osc(iex,2)*stpmder(iex,k,3)*Ls_htb(iex,k)* &
          hpol2si

        hpolcb(r,3,1,1) = hpolcb(r,3,1,1) + &
          osc(iex,3)*stpmder(iex,k,1)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,3,1,2) = hpolcb(r,3,1,2) + &
          osc(iex,3)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,3,1,3) = hpolcb(r,3,1,3) + &
          osc(iex,3)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
          hpol2si

        hpolcb(r,3,2,1) = hpolcb(r,3,2,1) + &
          osc(iex,3)*stpmder(iex,k,4)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,3,2,2) = hpolcb(r,3,2,2) + &
          osc(iex,3)*stpmder(iex,k,2)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,3,2,3) = hpolcb(r,3,2,3) + &
          osc(iex,3)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
          hpol2si

        hpolcb(r,3,3,1) = hpolcb(r,3,3,1) + &
          osc(iex,3)*stpmder(iex,k,5)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,3,3,2) = hpolcb(r,3,3,2) + &
          osc(iex,3)*stpmder(iex,k,6)*Ls_htb(iex,k)* &
          hpol2si
        hpolcb(r,3,3,3) = hpolcb(r,3,3,3) + &
          osc(iex,3)*stpmder(iex,k,3)*Ls_htb(iex,k)* &
          hpol2si
      end do
    end do
  end do

end subroutine RHRSB2OverCB
