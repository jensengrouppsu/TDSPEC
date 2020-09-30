Subroutine RRSFCFund(latermoff, x, w, freq, delta, osc, pol, freq0, &
                     lfreq, g, nexci, nmodes, tstep, polA, tquad, polG, mtdip, &
                     polC, polD, lprintAtensor, lprintGtensor)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep
  ! NB: Tquad is entered in the following order:
  ! Tquad XX, XY, XZ, YY, YZ, ZZ
  integer, parameter :: xx=1, xy=2, xz=3, yy=4, yz=5, zz=6

  complex(kindr), intent(inout)           :: pol(nmodes,3,3)
  complex(kindr)                          :: Ls, rsum, z
  complex(kindr)                          :: dsum(tstep,nexci)
  complex(kindr), intent(in)              :: g(nexci,tstep)
  complex(kindr), intent(inout)           :: polA(nmodes,3,3,3)
  complex(kindr), intent(inout)           :: polG(nmodes,3,3)
  real(kindr), intent(in)                 :: tquad(nexci,6)
  real(kindr), intent(in)                 :: mtdip(nexci,3)
  logical, intent(in)                     :: lprintAtensor
  logical, intent(in)                     :: lprintGtensor

  ! Higher order tensors
  complex(kindr), intent(inout)           :: polC(nmodes,6,6)
  complex(kindr), intent(inout)           :: polD(nmodes,6,3)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep) 
  real(kindr) :: osc(nexci,3)

  complex(kindr) :: zzero=(0d0,0d0) 

  logical :: latermoff

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate A term polarizability for RRS spectra (fundamentals)
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
! in/out : pol - RRS polarizability tensor
!          polA - A-tensor
!          polG - G-tensor
!          polC - C-tensor
!          polD - D-tensor
!
! =============================================================================

  ! This is needed because some compilers don't default arrays to being zero.
  pol = zzero
  polA = zzero
  polG = zzero
  polC = zzero
  polD = zzero

  if (latermoff.eqv..false.) then
    ! Calculate the sum term in the integral (redundant 
    ! part that applies both to Raman and absorption 
    ! spectra).
    dsum = zzero
    !$OMP PARALLEL
    !$OMP DO
    do iex=1,nexci
      do t=1,tstep
        do i=1,nmodes
          dsum(t,iex) = dsum(t,iex) + (delta(iex,i)**2/two)*(one-exp(-z*freq(i)*x(t)))
        end do
      end do
    end do
    !$OMP END DO

    !$OMP DO PRIVATE(Ls,rsum)
    do j = 1, nmodes
      do iex=1,nexci
        Ls = zzero
          do t=1,tstep
            ! Calculate the part in brackets { } of equation 19 in 
            ! the Neese and Petrenko paper, for k = 1.        

            ! Negative sign? 
            !rsum = -delta(iex,j)/sqrt(two)*(one-exp(-z*freq(j)*x(t)))
            ! Positive sign?
            rsum = delta(iex,j)/sqrt(two)*(one-exp(-z*freq(j)*x(t)))
            Ls = Ls + z*w(t)*rsum*exp(x(t)* &
              z*(lfreq-freq0(iex)) - g(iex,t) - dsum(t,iex))
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

        pol(j,1,1) = pol(j,1,1) + osc(iex,1)*osc(iex,1)*Ls* &
                     pol2si
        pol(j,1,2) = pol(j,1,2) + osc(iex,1)*osc(iex,2)*Ls* &
                     pol2si
        pol(j,1,3) = pol(j,1,3) + osc(iex,1)*osc(iex,3)*Ls* &
                     pol2si
        pol(j,2,1) = pol(j,2,1) + osc(iex,2)*osc(iex,1)*Ls* &
                     pol2si
        pol(j,2,2) = pol(j,2,2) + osc(iex,2)*osc(iex,2)*Ls* &
                     pol2si
        pol(j,2,3) = pol(j,2,3) + osc(iex,2)*osc(iex,3)*Ls* &
                     pol2si
        pol(j,3,1) = pol(j,3,1) + osc(iex,3)*osc(iex,1)*Ls* &
                     pol2si
        pol(j,3,2) = pol(j,3,2) + osc(iex,3)*osc(iex,2)*Ls* &
                     pol2si
        pol(j,3,3) = pol(j,3,3) + osc(iex,3)*osc(iex,3)*Ls* &
                     pol2si

        if (lprintAtensor) then
          ! Calculate the A-tensor if needed
          polA(j,1,1,1) = polA(j,1,1,1) + osc(iex,1)*Tquad(iex,xx)*Ls
          polA(j,1,1,2) = polA(j,1,1,2) + osc(iex,1)*Tquad(iex,xy)*Ls
          polA(j,1,1,3) = polA(j,1,1,3) + osc(iex,1)*Tquad(iex,xz)*Ls
          polA(j,1,2,1) = polA(j,1,2,1) + osc(iex,1)*Tquad(iex,xy)*Ls
          polA(j,1,2,2) = polA(j,1,2,2) + osc(iex,1)*Tquad(iex,yy)*Ls
          polA(j,1,2,3) = polA(j,1,2,3) + osc(iex,1)*Tquad(iex,yz)*Ls
          polA(j,1,3,1) = polA(j,1,3,1) + osc(iex,1)*Tquad(iex,xz)*Ls
          polA(j,1,3,2) = polA(j,1,3,2) + osc(iex,1)*Tquad(iex,yz)*Ls
          polA(j,1,3,3) = polA(j,1,3,3) + osc(iex,1)*Tquad(iex,zz)*Ls

          polA(j,2,1,1) = polA(j,2,1,1) + osc(iex,2)*Tquad(iex,xx)*Ls
          polA(j,2,1,2) = polA(j,2,1,2) + osc(iex,2)*Tquad(iex,xy)*Ls
          polA(j,2,1,3) = polA(j,2,1,3) + osc(iex,2)*Tquad(iex,xz)*Ls
          polA(j,2,2,1) = polA(j,2,2,1) + osc(iex,2)*Tquad(iex,xy)*Ls
          polA(j,2,2,2) = polA(j,2,2,2) + osc(iex,2)*Tquad(iex,yy)*Ls
          polA(j,2,2,3) = polA(j,2,2,3) + osc(iex,2)*Tquad(iex,yz)*Ls
          polA(j,2,3,1) = polA(j,2,3,1) + osc(iex,2)*Tquad(iex,xz)*Ls
          polA(j,2,3,2) = polA(j,2,3,2) + osc(iex,2)*Tquad(iex,yz)*Ls
          polA(j,2,3,3) = polA(j,2,3,3) + osc(iex,2)*Tquad(iex,zz)*Ls

          polA(j,3,1,1) = polA(j,3,1,1) + osc(iex,3)*Tquad(iex,xx)*Ls
          polA(j,3,1,2) = polA(j,3,1,2) + osc(iex,3)*Tquad(iex,xy)*Ls
          polA(j,3,1,3) = polA(j,3,1,3) + osc(iex,3)*Tquad(iex,xz)*Ls
          polA(j,3,2,1) = polA(j,3,2,1) + osc(iex,3)*Tquad(iex,xy)*Ls
          polA(j,3,2,2) = polA(j,3,2,2) + osc(iex,3)*Tquad(iex,yy)*Ls
          polA(j,3,2,3) = polA(j,3,2,3) + osc(iex,3)*Tquad(iex,yz)*Ls
          polA(j,3,3,1) = polA(j,3,3,1) + osc(iex,3)*Tquad(iex,xz)*Ls
          polA(j,3,3,2) = polA(j,3,3,2) + osc(iex,3)*Tquad(iex,yz)*Ls
          polA(j,3,3,3) = polA(j,3,3,3) + osc(iex,3)*Tquad(iex,zz)*Ls

          ! Calculate the C-tensor, because we can
          polC(j,1,1) = polC(j,1,1) + Tquad(iex,xx)*Tquad(iex,xx)*Ls
          polC(j,1,2) = polC(j,1,2) + Tquad(iex,xx)*Tquad(iex,xy)*Ls
          polC(j,1,3) = polC(j,1,3) + Tquad(iex,xx)*Tquad(iex,xz)*Ls
          polC(j,1,4) = polC(j,1,4) + Tquad(iex,xx)*Tquad(iex,yy)*Ls
          polC(j,1,5) = polC(j,1,5) + Tquad(iex,xx)*Tquad(iex,yz)*Ls
          polC(j,1,6) = polC(j,1,6) + Tquad(iex,xx)*Tquad(iex,zz)*Ls

          polC(j,2,1) = polC(j,2,1) + Tquad(iex,xy)*Tquad(iex,xx)*Ls
          polC(j,2,2) = polC(j,2,2) + Tquad(iex,xy)*Tquad(iex,xy)*Ls
          polC(j,2,3) = polC(j,2,3) + Tquad(iex,xy)*Tquad(iex,xz)*Ls
          polC(j,2,4) = polC(j,2,4) + Tquad(iex,xy)*Tquad(iex,yy)*Ls
          polC(j,2,5) = polC(j,2,5) + Tquad(iex,xy)*Tquad(iex,yz)*Ls
          polC(j,2,6) = polC(j,2,6) + Tquad(iex,xy)*Tquad(iex,zz)*Ls

          polC(j,3,1) = polC(j,3,1) + Tquad(iex,xz)*Tquad(iex,xx)*Ls
          polC(j,3,2) = polC(j,3,2) + Tquad(iex,xz)*Tquad(iex,xy)*Ls
          polC(j,3,3) = polC(j,3,3) + Tquad(iex,xz)*Tquad(iex,xz)*Ls
          polC(j,3,4) = polC(j,3,4) + Tquad(iex,xz)*Tquad(iex,yy)*Ls
          polC(j,3,5) = polC(j,3,5) + Tquad(iex,xz)*Tquad(iex,yz)*Ls
          polC(j,3,6) = polC(j,3,6) + Tquad(iex,xz)*Tquad(iex,zz)*Ls

          polC(j,4,1) = polC(j,4,1) + Tquad(iex,yy)*Tquad(iex,xx)*Ls
          polC(j,4,2) = polC(j,4,2) + Tquad(iex,yy)*Tquad(iex,xy)*Ls
          polC(j,4,3) = polC(j,4,3) + Tquad(iex,yy)*Tquad(iex,xz)*Ls
          polC(j,4,4) = polC(j,4,4) + Tquad(iex,yy)*Tquad(iex,yy)*Ls
          polC(j,4,5) = polC(j,4,5) + Tquad(iex,yy)*Tquad(iex,yz)*Ls
          polC(j,4,6) = polC(j,4,6) + Tquad(iex,yy)*Tquad(iex,zz)*Ls

          polC(j,5,1) = polC(j,5,1) + Tquad(iex,yz)*Tquad(iex,xx)*Ls
          polC(j,5,2) = polC(j,5,2) + Tquad(iex,yz)*Tquad(iex,xy)*Ls
          polC(j,5,3) = polC(j,5,3) + Tquad(iex,yz)*Tquad(iex,xz)*Ls
          polC(j,5,4) = polC(j,5,4) + Tquad(iex,yz)*Tquad(iex,yy)*Ls
          polC(j,5,5) = polC(j,5,5) + Tquad(iex,yz)*Tquad(iex,yz)*Ls
          polC(j,5,6) = polC(j,5,6) + Tquad(iex,yz)*Tquad(iex,zz)*Ls

          polC(j,6,1) = polC(j,6,1) + Tquad(iex,zz)*Tquad(iex,xx)*Ls
          polC(j,6,2) = polC(j,6,2) + Tquad(iex,zz)*Tquad(iex,xy)*Ls
          polC(j,6,3) = polC(j,6,3) + Tquad(iex,zz)*Tquad(iex,xz)*Ls
          polC(j,6,4) = polC(j,6,4) + Tquad(iex,zz)*Tquad(iex,yy)*Ls
          polC(j,6,5) = polC(j,6,5) + Tquad(iex,zz)*Tquad(iex,yz)*Ls
          polC(j,6,6) = polC(j,6,6) + Tquad(iex,zz)*Tquad(iex,zz)*Ls
       end if

       ! Calcuate the G-tensor
       if (lprintGtensor) then
          polG(j,1,1) = polG(j,1,1) + osc(iex,1)*MTdip(iex,1)*Ls
          polG(j,1,2) = polG(j,1,2) + osc(iex,1)*MTdip(iex,2)*Ls
          polG(j,1,3) = polG(j,1,3) + osc(iex,1)*MTdip(iex,3)*Ls

          polG(j,2,1) = polG(j,2,1) + osc(iex,2)*MTdip(iex,1)*Ls
          polG(j,2,2) = polG(j,2,2) + osc(iex,2)*MTdip(iex,2)*Ls
          polG(j,2,3) = polG(j,2,3) + osc(iex,2)*MTdip(iex,3)*Ls

          polG(j,3,1) = polG(j,3,1) + osc(iex,3)*MTdip(iex,1)*Ls
          polG(j,3,2) = polG(j,3,2) + osc(iex,3)*MTdip(iex,2)*Ls
          polG(j,3,3) = polG(j,3,3) + osc(iex,3)*MTdip(iex,3)*Ls
       end if

       ! Calculate the D-tensor
       if (lprintAtensor .and. lprintGtensor) then
          polD(j,1,1) = polD(j,1,1) + Tquad(iex,xx)*MTdip(iex,1)*Ls
          polD(j,1,2) = polD(j,1,2) + Tquad(iex,xx)*MTdip(iex,2)*Ls
          polD(j,1,3) = polD(j,1,3) + Tquad(iex,xx)*MTdip(iex,3)*Ls

          polD(j,2,1) = polD(j,2,1) + Tquad(iex,xy)*MTdip(iex,1)*Ls
          polD(j,2,2) = polD(j,2,2) + Tquad(iex,xy)*MTdip(iex,2)*Ls
          polD(j,2,3) = polD(j,2,3) + Tquad(iex,xy)*MTdip(iex,3)*Ls

          polD(j,3,1) = polD(j,3,1) + Tquad(iex,xz)*MTdip(iex,1)*Ls
          polD(j,3,2) = polD(j,3,2) + Tquad(iex,xz)*MTdip(iex,2)*Ls
          polD(j,3,3) = polD(j,3,3) + Tquad(iex,xz)*MTdip(iex,3)*Ls

          polD(j,4,1) = polD(j,4,1) + Tquad(iex,yy)*MTdip(iex,1)*Ls
          polD(j,4,2) = polD(j,4,2) + Tquad(iex,yy)*MTdip(iex,2)*Ls
          polD(j,4,3) = polD(j,4,3) + Tquad(iex,yy)*MTdip(iex,3)*Ls

          polD(j,5,1) = polD(j,5,1) + Tquad(iex,yz)*MTdip(iex,1)*Ls
          polD(j,5,2) = polD(j,5,2) + Tquad(iex,yz)*MTdip(iex,2)*Ls
          polD(j,5,3) = polD(j,5,3) + Tquad(iex,yz)*MTdip(iex,3)*Ls

          polD(j,6,1) = polD(j,6,1) + Tquad(iex,zz)*MTdip(iex,1)*Ls
          polD(j,6,2) = polD(j,6,2) + Tquad(iex,zz)*MTdip(iex,2)*Ls
          polD(j,6,3) = polD(j,6,3) + Tquad(iex,zz)*MTdip(iex,3)*Ls
       end if

      end do ! iex

      ! Test
      !write(6,*)
      !write(6,*) 'A term polarizability'
      !write(6,"(1X,A,1X,F7.2)") 'Frequency', freq(j)
      !write(6,*) '--------------------------------'
      !write(6,*) 'Real part (a.u.)'
      !write(6,*) '--------------------------------'
      !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
      !           real(pol(j,1,1))*hartree2cm/pol2si, &
      !           real(pol(j,1,2))*hartree2cm/pol2si, &
      !           real(pol(j,1,3))*hartree2cm/pol2si
      !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
      !           real(pol(j,2,1))*hartree2cm/pol2si, &
      !           real(pol(j,2,2))*hartree2cm/pol2si, &
      !           real(pol(j,2,3))*hartree2cm/pol2si
      !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
      !           real(pol(j,3,1))*hartree2cm/pol2si, &
      !           real(pol(j,3,2))*hartree2cm/pol2si, &
      !           real(pol(j,3,3))*hartree2cm/pol2si
      !write(6,*) '--------------------------------'
      !write(6,*) 'Imaginary part (a.u.)'
      !write(6,*) '--------------------------------'
      !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
      !           aimag(pol(j,1,1))*hartree2cm/pol2si, &
      !           aimag(pol(j,1,2))*hartree2cm/pol2si, &
      !           aimag(pol(j,1,3))*hartree2cm/pol2si
      !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
      !           aimag(pol(j,2,1))*hartree2cm/pol2si, &
      !           aimag(pol(j,2,2))*hartree2cm/pol2si, &
      !           aimag(pol(j,2,3))*hartree2cm/pol2si
      !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
      !           aimag(pol(j,3,1))*hartree2cm/pol2si, &
      !           aimag(pol(j,3,2))*hartree2cm/pol2si, &
      !           aimag(pol(j,3,3))*hartree2cm/pol2si
      !write(6,*) '--------------------------------'
    end do ! j
    !$OMP END DO
    !$OMP END PARALLEL
  end if ! latermoff

end subroutine RRSFCFund

Subroutine RRSFCOverCB(latermoff, x, w, freq, delta, osc, polot, &
                     polcb, freq0, lfreq, g, nexci, nmodes, &
                     tstep, excnum, numcb, excnums, factorial)

  use constants

  implicit none

  integer :: nexci, nmodes, excnum, numcb
  integer :: j, iex, t, i, tstep, r, k
  integer :: excnums(numcb,nmodes)

  complex(kindr), intent(inout) :: polot(excnum,nmodes,3,3)
  complex(kindr), intent(inout) :: polcb(numcb,3,3)
  complex(kindr) :: Ls, rsum, dsum, z
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3)
  real(kindr) :: factorial(16)

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: latermoff

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate A term polarizability for RRS spectra (overtones/combination bands)
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
!          excnum - excitation level
!          numcb - number of important combination bands
!          excnums - array of each mode's state
!          factorial - array of factorial values
!
! in/out : polot - RRS polarizability tensor for overtones
!          polcb - RRS polarizability tensor for combination bands
!
! =============================================================================

  if (excnum > 1) then
    ! Calculate the transition polarizability for the A term for 
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
                     sqrt(two**k*factorial(k))* &
                     (one-exp(-z*freq(j)*x(t)))**k
              Ls = Ls + z*w(t)*rsum*exp(x(t)* &
                  z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
            end do
            polot(k,j,1,1) = polot(k,j,1,1) + osc(iex,1)*osc(iex,1)* &
                         Ls*pol2si
            polot(k,j,1,2) = polot(k,j,1,2) + osc(iex,1)*osc(iex,2)* &
                         Ls*pol2si
            polot(k,j,1,3) = polot(k,j,1,3) + osc(iex,1)*osc(iex,3)* &
                         Ls*pol2si
            polot(k,j,2,1) = polot(k,j,2,1) + osc(iex,2)*osc(iex,1)* &
                         Ls*pol2si
            polot(k,j,2,2) = polot(k,j,2,2) + osc(iex,2)*osc(iex,2)* &
                         Ls*pol2si
            polot(k,j,2,3) = polot(k,j,2,3) + osc(iex,2)*osc(iex,3)* &
                         Ls*pol2si
            polot(k,j,3,1) = polot(k,j,3,1) + osc(iex,3)*osc(iex,1)* &
                         Ls*pol2si
            polot(k,j,3,2) = polot(k,j,3,2) + osc(iex,3)*osc(iex,2)* &
                         Ls*pol2si
            polot(k,j,3,3) = polot(k,j,3,3) + osc(iex,3)*osc(iex,3)* &
                         Ls*pol2si
          end do ! iex (nexci)
        end do ! j (nmodes)
      end do ! k (ExcNum)
    end if
    ! Calculate the transition polarizability for the A term for
    ! combination bands
    if (latermoff.eqv..false.) then
      do k=1,numcb
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
            do r=1,nmodes
              if (excnums(k,r) /= 0) then
                rsum = rsum* &
                    (delta(iex,r))**excnums(k,r)/ &
                    sqrt(two**excnums(k,r)* &
                    factorial(excnums(k,r)))* &
                    (one-exp(-z*freq(r)*x(t)))**excnums(k,r)
              end if
            end do
            Ls = Ls + z*w(t)*rsum*exp(x(t)* &
                z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
          end do
          polcb(k,1,1) = polcb(k,1,1) + osc(iex,1)*osc(iex,1)* &
                   Ls*pol2si
          polcb(k,1,2) = polcb(k,1,2) + osc(iex,1)*osc(iex,2)* &
                   Ls*pol2si
          polcb(k,1,3) = polcb(k,1,3) + osc(iex,1)*osc(iex,3)* &
                   Ls*pol2si
          polcb(k,2,1) = polcb(k,2,1) + osc(iex,2)*osc(iex,1)* &
                   Ls*pol2si
          polcb(k,2,2) = polcb(k,2,2) + osc(iex,2)*osc(iex,2)* &
                   Ls*pol2si
          polcb(k,2,3) = polcb(k,2,3) + osc(iex,2)*osc(iex,3)* &
                   Ls*pol2si
          polcb(k,3,1) = polcb(k,3,1) + osc(iex,3)*osc(iex,1)* &
                   Ls*pol2si
          polcb(k,3,2) = polcb(k,3,2) + osc(iex,3)*osc(iex,2)* &
                   Ls*pol2si
          polcb(k,3,3) = polcb(k,3,3) + osc(iex,3)*osc(iex,3)* &
                   Ls*pol2si
        end do
      end do
    end if
  end if

end subroutine RRSFCOverCB

Subroutine RRSHTFund(x, w, freq, delta, osc, tdipder, pol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     Ls_htb, Ls_htb2, Ls_htb3, htlfreqdiff, &
                     ht_pol, ht_rcrs, ht_rcrs_sum, lhtints, &
                     polA, polG, polC, polD, tquad, tquadder, mtdip, &
                     mtdipder, polAs, polGs, polDs, &
                     lprintAtensor, lprintGtensor)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep, k, n
  ! NB: Tquad is entered in the following order:
  ! Tquad XX, XY, XZ, YY, YZ, ZZ
  integer, parameter :: xx=1, xy=2, xz=3, yy=4, yz=5, zz=6

  complex(kindr), intent(inout) :: pol(nmodes,3,3)
  complex(kindr), intent(inout) :: polA(nmodes,3,3,3)
  complex(kindr), intent(inout) :: polG(nmodes,3,3)
  complex(kindr), intent(inout) :: polC(nmodes,6,6)
  complex(kindr), intent(inout) :: polD(nmodes,6,3)
  complex(kindr), intent(inout) :: polAs(nmodes,3,3,3)
  complex(kindr), intent(inout) :: polGs(nmodes,3,3)
  complex(kindr), intent(inout) :: polDs(nmodes,3,6)
  complex(kindr) :: ht_pol(nmodes,nmodes,3,3)
  complex(kindr) :: rsum, rsum1, z
  complex(kindr) :: dsum(nexci,tstep)
  complex(kindr) :: expfac
  complex(kindr) :: Ls_htb(nexci,nmodes), Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr) :: expsum(nexci,tstep,nmodes)
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3),tdipder(nexci,nmodes,3)
  real(kindr) :: htlfreqdiff(nexci)
  real(kindr) :: ht_rcrs(nmodes,nmodes)
  real(kindr) :: ht_rcrs_sum(nmodes)
  real(kindr) :: htexcist
  real(kindr), intent(in) :: tquad(nexci,6)
  real(kindr), intent(in) :: tquadder(nexci,nmodes,6)
  real(kindr), intent(in) :: mtdip(nexci,3)
  real(kindr), intent(in) :: mtdipder(nexci,nmodes,3)

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: lhtints
  logical, intent(in) :: lprintAtensor
  logical, intent(in) :: lprintGtensor

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B term polarizability for RRS spectra (fundamentals)
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
!          lhtints - HT integrals are calculated?
!
! in/out : pol - RRS polarizability tensor
!
! =============================================================================

! Calculations of the Herzberg-Teller terms for resonance Raman 
! scattering are based on equations for time dependent integrals 
! and ideas in the following sources
!
! Spiro, T. G.  "Biological Applications of Raman Spectrscopy".  
!   John Wiley & Sons, Inc.  1987.
! Markel, F.; Ferris, N. S.; Gould, I. R. and Myers, A. B.  
!   J. Am. Chem. Soc. 114, 6208 (1992).

! Keep in mind that the initial and final states of the system
! are actually vectors containing the state of each harmonic 
! oscillator.  Since the overlap integral looks like sum_a<F|Qa|I>,
! each normal mode is raised and lowered, not just the mode in
! its excited state.

  dsum = zzero
  !$OMP PARALLEL 
  !$OMP DO
  do t=1,tstep
    do iex=1,nexci
      do i=1,nmodes
        ! Generate dsum
        dsum(iex,t) = dsum(iex,t) + (delta(iex,i)**2/two)*(one-exp(-z*freq(i)*x(t)))
        ! Calculate the contribution from other normal modes.  This
        ! is either raised or has nothing performed, if the mode is
        ! the current mode.
        expsum(iex,t,i) = delta(iex,i)/(root2) * (one-exp(-z*freq(i)*x(t)))
      end do
    end do
  end do
  !$OMP END DO

  !$OMP DO PRIVATE(Ls_htb,Ls_htb2,Ls_htb3,rsum,rsum1,expfac)
  do j = 1, nmodes
    do iex=1,nexci
      Ls_htb(1:nexci,1:nmodes) = zzero
      Ls_htb2(1:nexci,1:nmodes) = zzero
      Ls_htb3(1:nexci,1:nmodes) = zzero
      do t=1,tstep
        ! Calculate the contribution from the excited normal mode.
        ! This is either raised (rsum) or has nothing performed (rsum1).
!        rsum = delta(iex,j)**2/(root8)* &
!               (one-exp(-z*freq(j)*x(t)))**2
        rsum1 = delta(iex,j)/(root2)* &
                (one-exp(-z*freq(j)*x(t)))
        rsum = rsum1**2 / root2
        expfac = z*w(t)*exp(x(t)*z*(lfreq-freq0(iex)) - g(iex,t) - dsum(iex,t))
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
                     z*(lfreq-freq0(iex)-freq(n)) - g(iex,t) - dsum(iex,t))
          else
            ! Raise the initial state of a different mode,
            ! so I+1 = 1 and F = 0.  Here, there's a product
            ! accounting for multiple modes in an excited state.
             Ls_htb(iex,n) = Ls_htb(iex,n) + rsum1*expsum(iex,t,n)*expfac
          end if
          ! Integrals of the form: <(F-1)|I(t)>
          if (n == j) then
            ! Lower the final state of the current mode,
            ! so that F-1 = I = 0
            Ls_htb2(iex,n) = Ls_htb2(iex,n) + expfac
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
                      root2*rsum*expfac
          else
            ! Raise the final state of a different mode,
            ! so that I = 0 and F = 1.
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      expsum(iex,t,n)*rsum1*expfac
          end if
        end do ! n (nmodes)
      end do ! t (tstep)
      do k = 1, nmodes
        pol(j,1,1) = pol(j,1,1) + &
                     (osc(iex,1)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,1)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        pol(j,1,2) = pol(j,1,2) + &
                     (osc(iex,1)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,1)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        pol(j,1,3) = pol(j,1,3) + &
                     (osc(iex,1)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,1)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        pol(j,2,1) = pol(j,2,1) + &
                     (osc(iex,2)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,2)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        pol(j,2,2) = pol(j,2,2) + &
                     (osc(iex,2)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,2)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        pol(j,2,3) = pol(j,2,3) + &
                     (osc(iex,2)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,2)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        pol(j,3,1) = pol(j,3,1) + &
                     (osc(iex,3)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,3)*osc(iex,1)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        pol(j,3,2) = pol(j,3,2) + &
                     (osc(iex,3)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,3)*osc(iex,2)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si
        pol(j,3,3) = pol(j,3,3) + &
                     (osc(iex,3)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
                     tdipder(iex,k,3)*osc(iex,3)*(Ls_htb2(iex,k)+ &
                     Ls_htb3(iex,k)))*pol2si

        if (lprintAtensor) then
          ! Do A-tensor (Roman-type)
          polA(j,1,1,1) = polA(j,1,1,1) + osc(iex,1)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,1,1,2) = polA(j,1,1,2) + osc(iex,1)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,1,1,3) = polA(j,1,1,3) + osc(iex,1)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,1,2,1) = polA(j,1,2,1) + osc(iex,1)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,1,2,2) = polA(j,1,2,2) + osc(iex,1)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,1,2,3) = polA(j,1,2,3) + osc(iex,1)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,1,3,1) = polA(j,1,3,1) + osc(iex,1)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,1,3,2) = polA(j,1,3,2) + osc(iex,1)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,1,3,3) = polA(j,1,3,3) + osc(iex,1)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                          tdipder(iex,k,1)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polA(j,2,1,1) = polA(j,2,1,1) + osc(iex,2)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,2,1,2) = polA(j,2,1,2) + osc(iex,2)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,2,1,3) = polA(j,2,1,3) + osc(iex,2)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,2,2,1) = polA(j,2,2,1) + osc(iex,2)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,2,2,2) = polA(j,2,2,2) + osc(iex,2)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,2,2,3) = polA(j,2,2,3) + osc(iex,2)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,2,3,1) = polA(j,2,3,1) + osc(iex,2)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,2,3,2) = polA(j,2,3,2) + osc(iex,2)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,2,3,3) = polA(j,2,3,3) + osc(iex,2)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                          tdipder(iex,k,2)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polA(j,3,1,1) = polA(j,3,1,1) + osc(iex,3)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,3,1,2) = polA(j,3,1,2) + osc(iex,3)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,3,1,3) = polA(j,3,1,3) + osc(iex,3)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,3,2,1) = polA(j,3,2,1) + osc(iex,3)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,3,2,2) = polA(j,3,2,2) + osc(iex,3)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,3,2,3) = polA(j,3,2,3) + osc(iex,3)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,3,3,1) = polA(j,3,3,1) + osc(iex,3)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,3,3,2) = polA(j,3,3,2) + osc(iex,3)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polA(j,3,3,3) = polA(j,3,3,3) + osc(iex,3)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                          tdipder(iex,k,3)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          ! Do A-tensor (script-type)
          polAs(j,1,1,1) = polAs(j,1,1,1) + tdipder(iex,k,1)*Tquad(iex,1)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,1,1,2) = polAs(j,1,1,2) + tdipder(iex,k,1)*Tquad(iex,2)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,1,1,3) = polAs(j,1,1,3) + tdipder(iex,k,1)*Tquad(iex,3)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,1,2,1) = polAs(j,1,2,1) + tdipder(iex,k,1)*Tquad(iex,2)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,1,2,2) = polAs(j,1,2,2) + tdipder(iex,k,1)*Tquad(iex,4)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,1,2,3) = polAs(j,1,2,3) + tdipder(iex,k,1)*Tquad(iex,5)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,1,3,1) = polAs(j,1,3,1) + tdipder(iex,k,1)*Tquad(iex,3)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,1,3,2) = polAs(j,1,3,2) + tdipder(iex,k,1)*Tquad(iex,5)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,1,3,3) = polAs(j,1,3,3) + tdipder(iex,k,1)*Tquad(iex,6)*Ls_htb(iex,k)+&
                          osc(iex,1)*Tquadder(iex,k,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polAs(j,2,1,1) = polAs(j,2,1,1) + tdipder(iex,k,2)*Tquad(iex,1)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,2,1,2) = polAs(j,2,1,2) + tdipder(iex,k,2)*Tquad(iex,2)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,2,1,3) = polAs(j,2,1,3) + tdipder(iex,k,2)*Tquad(iex,3)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,2,2,1) = polAs(j,2,2,1) + tdipder(iex,k,2)*Tquad(iex,2)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,2,2,2) = polAs(j,2,2,2) + tdipder(iex,k,2)*Tquad(iex,4)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,2,2,3) = polAs(j,2,2,3) + tdipder(iex,k,2)*Tquad(iex,5)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,2,3,1) = polAs(j,2,3,1) + tdipder(iex,k,2)*Tquad(iex,3)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,2,3,2) = polAs(j,2,3,2) + tdipder(iex,k,2)*Tquad(iex,5)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,2,3,3) = polAs(j,2,3,3) + tdipder(iex,k,2)*Tquad(iex,6)*Ls_htb(iex,k)+&
                          osc(iex,2)*Tquadder(iex,k,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polAs(j,3,1,1) = polAs(j,3,1,1) + tdipder(iex,k,3)*Tquad(iex,1)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,3,1,2) = polAs(j,3,1,2) + tdipder(iex,k,3)*Tquad(iex,2)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,3,1,3) = polAs(j,3,1,3) + tdipder(iex,k,3)*Tquad(iex,3)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,3,2,1) = polAs(j,3,2,1) + tdipder(iex,k,3)*Tquad(iex,2)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,3,2,2) = polAs(j,3,2,2) + tdipder(iex,k,3)*Tquad(iex,4)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,3,2,3) = polAs(j,3,2,3) + tdipder(iex,k,3)*Tquad(iex,5)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,3,3,1) = polAs(j,3,3,1) + tdipder(iex,k,3)*Tquad(iex,3)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,3,3,2) = polAs(j,3,3,2) + tdipder(iex,k,3)*Tquad(iex,5)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polAs(j,3,3,3) = polAs(j,3,3,3) + tdipder(iex,k,3)*Tquad(iex,6)*Ls_htb(iex,k)+&
                          osc(iex,3)*Tquadder(iex,k,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          ! Do C-tensor
          polC(j,1,1) = polC(j,1,1) + Tquad(iex,1)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,1,2) = polC(j,1,2) + Tquad(iex,1)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,1,3) = polC(j,1,3) + Tquad(iex,1)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,1,4) = polC(j,1,4) + Tquad(iex,1)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,1,5) = polC(j,1,5) + Tquad(iex,1)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,1,6) = polC(j,1,6) + Tquad(iex,1)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polC(j,2,1) = polC(j,2,1) + Tquad(iex,2)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,2,2) = polC(j,2,2) + Tquad(iex,2)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,2,3) = polC(j,2,3) + Tquad(iex,2)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,2,4) = polC(j,2,4) + Tquad(iex,2)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,2,5) = polC(j,2,5) + Tquad(iex,2)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,2,6) = polC(j,2,6) + Tquad(iex,2)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polC(j,3,1) = polC(j,3,1) + Tquad(iex,3)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,3,2) = polC(j,3,2) + Tquad(iex,3)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,3,3) = polC(j,3,3) + Tquad(iex,3)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,3,4) = polC(j,3,4) + Tquad(iex,3)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,3,5) = polC(j,3,5) + Tquad(iex,3)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,3,6) = polC(j,3,6) + Tquad(iex,3)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polC(j,4,1) = polC(j,4,1) + Tquad(iex,4)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,4,2) = polC(j,4,2) + Tquad(iex,4)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,4,3) = polC(j,4,3) + Tquad(iex,4)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,4,4) = polC(j,4,4) + Tquad(iex,4)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,4,5) = polC(j,4,5) + Tquad(iex,4)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,4,6) = polC(j,4,6) + Tquad(iex,4)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polC(j,5,1) = polC(j,5,1) + Tquad(iex,5)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,5,2) = polC(j,5,2) + Tquad(iex,5)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,5,3) = polC(j,5,3) + Tquad(iex,5)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,5,4) = polC(j,5,4) + Tquad(iex,5)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,5,5) = polC(j,5,5) + Tquad(iex,5)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,5,6) = polC(j,5,6) + Tquad(iex,5)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polC(j,6,1) = polC(j,6,1) + Tquad(iex,6)*Tquadder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*Tquad(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,6,2) = polC(j,6,2) + Tquad(iex,6)*Tquadder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*Tquad(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,6,3) = polC(j,6,3) + Tquad(iex,6)*Tquadder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*Tquad(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,6,4) = polC(j,6,4) + Tquad(iex,6)*Tquadder(iex,k,4)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*Tquad(iex,4)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,6,5) = polC(j,6,5) + Tquad(iex,6)*Tquadder(iex,k,5)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*Tquad(iex,5)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polC(j,6,6) = polC(j,6,6) + Tquad(iex,6)*Tquadder(iex,k,6)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*Tquad(iex,6)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
        end if

        if (lprintGtensor) then
          ! DO G-tensor (Roman-type)
          polG(j,1,1) = polG(j,1,1) + osc(iex,1)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        tdipder(iex,k,1)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polG(j,1,2) = polG(j,1,2) + osc(iex,1)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        tdipder(iex,k,1)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polG(j,1,3) = polG(j,1,3) + osc(iex,1)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        tdipder(iex,k,1)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polG(j,2,1) = polG(j,2,1) + osc(iex,2)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        tdipder(iex,k,2)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polG(j,2,2) = polG(j,2,2) + osc(iex,2)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        tdipder(iex,k,2)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polG(j,2,3) = polG(j,2,3) + osc(iex,2)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        tdipder(iex,k,2)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polG(j,3,1) = polG(j,3,1) + osc(iex,3)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        tdipder(iex,k,3)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polG(j,3,2) = polG(j,3,2) + osc(iex,3)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        tdipder(iex,k,3)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polG(j,3,3) = polG(j,3,3) + osc(iex,3)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        tdipder(iex,k,3)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          ! Do G-tensor (script-type)
          polGs(j,1,1) = polGs(j,1,1) - MTdip(iex,1)*tdipder(iex,k,1)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,1)*osc(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polGs(j,1,2) = polGs(j,1,2) - MTdip(iex,1)*tdipder(iex,k,2)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,1)*osc(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polGs(j,1,3) = polGs(j,1,3) - MTdip(iex,1)*tdipder(iex,k,3)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,1)*osc(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polGs(j,2,1) = polGs(j,2,1) - MTdip(iex,2)*tdipder(iex,k,1)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,2)*osc(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polGs(j,2,2) = polGs(j,2,2) - MTdip(iex,2)*tdipder(iex,k,2)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,2)*osc(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polGs(j,2,3) = polGs(j,2,3) - MTdip(iex,2)*tdipder(iex,k,3)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,2)*osc(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polGs(j,3,1) = polGs(j,3,1) - MTdip(iex,3)*tdipder(iex,k,1)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,3)*osc(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polGs(j,3,2) = polGs(j,3,2) - MTdip(iex,3)*tdipder(iex,k,2)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,3)*osc(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polGs(j,3,3) = polGs(j,3,3) - MTdip(iex,3)*tdipder(iex,k,3)*Ls_htb(iex,k)-&
                        MTdipder(iex,k,3)*osc(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
        end if

        if (lprintAtensor .and. lprintGtensor) then
          ! Do D-tensor (Roman-type)
          polD(j,1,1) = polD(j,1,1) + Tquad(iex,1)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,1,2) = polD(j,1,2) + Tquad(iex,1)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,1,3) = polD(j,1,3) + Tquad(iex,1)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,1)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polD(j,2,1) = polD(j,2,1) + Tquad(iex,2)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,2,2) = polD(j,2,2) + Tquad(iex,2)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,2,3) = polD(j,2,3) + Tquad(iex,2)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,2)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polD(j,3,1) = polD(j,3,1) + Tquad(iex,3)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,3,2) = polD(j,3,2) + Tquad(iex,3)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,3,3) = polD(j,3,3) + Tquad(iex,3)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,3)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polD(j,4,1) = polD(j,4,1) + Tquad(iex,4)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,4,2) = polD(j,4,2) + Tquad(iex,4)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,4,3) = polD(j,4,3) + Tquad(iex,4)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,4)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polD(j,5,1) = polD(j,5,1) + Tquad(iex,5)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,5,2) = polD(j,5,2) + Tquad(iex,5)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,5,3) = polD(j,5,3) + Tquad(iex,5)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,5)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polD(j,6,1) = polD(j,6,1) + Tquad(iex,6)*MTdipder(iex,k,1)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*MTdip(iex,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,6,2) = polD(j,6,2) + Tquad(iex,6)*MTdipder(iex,k,2)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*MTdip(iex,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polD(j,6,3) = polD(j,6,3) + Tquad(iex,6)*MTdipder(iex,k,3)*Ls_htb(iex,k)+&
                        Tquadder(iex,k,6)*MTdip(iex,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          ! Do D-tensor (script-type)
          polDs(j,1,1) = polDs(j,1,1) - Tquadder(iex,k,1)*MTdip(iex,1)*Ls_htb(iex,k)-&
                        Tquad(iex,1)*MTdipder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,2,1) = polDs(j,2,1) - Tquadder(iex,k,1)*MTdip(iex,2)*Ls_htb(iex,k)-&
                        Tquad(iex,1)*MTdipder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,3,1) = polDs(j,3,1) - Tquadder(iex,k,1)*MTdip(iex,3)*Ls_htb(iex,k)-&
                        Tquad(iex,1)*MTdipder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polDs(j,1,2) = polDs(j,1,2) - Tquadder(iex,k,2)*MTdip(iex,1)*Ls_htb(iex,k)-&
                        Tquad(iex,2)*MTdipder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,2,2) = polDs(j,2,2) - Tquadder(iex,k,2)*MTdip(iex,2)*Ls_htb(iex,k)-&
                        Tquad(iex,2)*MTdipder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,3,2) = polDs(j,3,2) - Tquadder(iex,k,2)*MTdip(iex,3)*Ls_htb(iex,k)-&
                        Tquad(iex,2)*MTdipder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polDs(j,1,3) = polDs(j,1,3) - Tquadder(iex,k,3)*MTdip(iex,1)*Ls_htb(iex,k)-&
                        Tquad(iex,3)*MTdipder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,2,3) = polDs(j,2,3) - Tquadder(iex,k,3)*MTdip(iex,2)*Ls_htb(iex,k)-&
                        Tquad(iex,3)*MTdipder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,3,3) = polDs(j,3,3) - Tquadder(iex,k,3)*MTdip(iex,3)*Ls_htb(iex,k)-&
                        Tquad(iex,3)*MTdipder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polDs(j,1,4) = polDs(j,1,4) - Tquadder(iex,k,4)*MTdip(iex,1)*Ls_htb(iex,k)-&
                        Tquad(iex,4)*MTdipder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,2,4) = polDs(j,2,4) - Tquadder(iex,k,4)*MTdip(iex,2)*Ls_htb(iex,k)-&
                        Tquad(iex,4)*MTdipder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,3,4) = polDs(j,3,4) - Tquadder(iex,k,4)*MTdip(iex,3)*Ls_htb(iex,k)-&
                        Tquad(iex,4)*MTdipder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polDs(j,1,5) = polDs(j,1,5) - Tquadder(iex,k,5)*MTdip(iex,1)*Ls_htb(iex,k)-&
                        Tquad(iex,5)*MTdipder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,2,5) = polDs(j,2,5) - Tquadder(iex,k,5)*MTdip(iex,2)*Ls_htb(iex,k)-&
                        Tquad(iex,5)*MTdipder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,3,5) = polDs(j,3,5) - Tquadder(iex,k,5)*MTdip(iex,3)*Ls_htb(iex,k)-&
                        Tquad(iex,5)*MTdipder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))

          polDs(j,1,6) = polDs(j,1,6) - Tquadder(iex,k,6)*MTdip(iex,1)*Ls_htb(iex,k)-&
                        Tquad(iex,6)*MTdipder(iex,k,1)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,2,6) = polDs(j,2,6) - Tquadder(iex,k,6)*MTdip(iex,2)*Ls_htb(iex,k)-&
                        Tquad(iex,6)*MTdipder(iex,k,2)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
          polDs(j,3,6) = polDs(j,3,6) - Tquadder(iex,k,6)*MTdip(iex,3)*Ls_htb(iex,k)-&
                        Tquad(iex,6)*MTdipder(iex,k,3)*(Ls_htb2(iex,k)+Ls_htb3(iex,k))
        end if

      end do ! k (nmodes)

      ! Evaluate the coupling between modes, if requested.
      if (lhtints) then
        do i = 1, nexci
          htlfreqdiff(i) = lfreq - freq0(i)
          htlfreqdiff(i) = abs(htlfreqdiff(i))
        end do
        htexcist = minval(htlfreqdiff)

        call RRSHTInts(Ls_htb, Ls_htb2, Ls_htb3, freq, osc, &
                       tdipder, freq0, htexcist, ht_pol, &
                       ht_rcrs, ht_rcrs_sum, delta, lfreq, nexci, &
                       nmodes, 3, iex, j)
      end if
    end do ! iex (nexci)
    ! Test
    !write(6,*)
    !write(6,*) 'A + B term polarizability'
    !write(6,"(1X,A,1X,F7.2)") 'Frequency', freq(j)
    !write(6,*) '--------------------------------'
    !write(6,*) 'Real part (a.u.)'
    !write(6,*) '--------------------------------'
    !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
    !           real(pol(j,1,1))*hartree2cm/pol2si, &
    !           real(pol(j,1,2))*hartree2cm/pol2si, &
    !           real(pol(j,1,3))*hartree2cm/pol2si
    !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
    !           real(pol(j,2,1))*hartree2cm/pol2si, &
    !           real(pol(j,2,2))*hartree2cm/pol2si, &
    !           real(pol(j,2,3))*hartree2cm/pol2si
    !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
    !           real(pol(j,3,1))*hartree2cm/pol2si, &
    !           real(pol(j,3,2))*hartree2cm/pol2si, &
    !           real(pol(j,3,3))*hartree2cm/pol2si
    !write(6,*) '--------------------------------'
    !write(6,*) 'Imaginary part (a.u.)'
    !write(6,*) '--------------------------------'
    !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
    !           aimag(pol(j,1,1))*hartree2cm/pol2si, &
    !           aimag(pol(j,1,2))*hartree2cm/pol2si, &
    !           aimag(pol(j,1,3))*hartree2cm/pol2si
    !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
    !           aimag(pol(j,2,1))*hartree2cm/pol2si, &
    !           aimag(pol(j,2,2))*hartree2cm/pol2si, &
    !           aimag(pol(j,2,3))*hartree2cm/pol2si
    !write(6, "(1X,F10.4,1X,F10.4,1X,F10.4)") &
    !           aimag(pol(j,3,1))*hartree2cm/pol2si, &
    !           aimag(pol(j,3,2))*hartree2cm/pol2si, &
    !           aimag(pol(j,3,3))*hartree2cm/pol2si
    !write(6,*) '--------------------------------'
  end do ! j (nmodes)
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine RRSHTFund

subroutine RRSHTOverCB(x, w, freq, delta, osc, tdipder, polot, polcb, &
               freq0, lfreq, g, nexci, nmodes, tstep, &
               excnum, numcb, excnums, factorial, rootn, num, &
               Ls_htb, Ls_htb2, Ls_htb3, sumr, sumr1, &
               sumr2)

  use constants

  implicit none

  integer :: nexci, nmodes, excnum, numcb
  integer :: j, iex, t, i, tstep, r, k, n, p, m
  integer :: excnums(numcb,nmodes)

  complex(kindr), intent(inout) :: polot(excnum,nmodes,3,3)
  complex(kindr), intent(inout) :: polcb(numcb,3,3)
  complex(kindr) :: rsum, dsum, z
  complex(kindr) :: Ls_htb(nexci,nmodes), Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr) :: sumr(nmodes), sumr1(nmodes), sumr2(nmodes)
  complex(kindr) :: expsum(nmodes)
  complex(kindr) :: rsum1, rsum2, rsum3
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3),tdipder(nexci,nmodes,3)
  real(kindr) :: factorial(16), rootn(12), num(11)

  complex(kindr) :: zzero=(0d0,0d0)

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate B term polarizability for RRS spectra (overtones/combination
!          bands)
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
!          excnum - excitation level
!          numcb - number of important combination bands
!          excnums - array of each mode's state
!          factorial - array of factorial values
!          rootn - array of square roots of double precision numbers
!          num - array of double precision numbers
!
! in/out : polot - RRS polarizability tensor for overtones
!          polcb - RRS polarizability tensor for combination bands
!
! =============================================================================

! Due to the raising operator, the loop involving excitation numbers
! goes from 3 to exnum+1 to be consistent with the Franck-Condon
! term.
!
! One needs to be cautious with the Herzberg-Teller term overtones 
! because the form of each equation is a little different, resulting
! from the raising and lowering operators.  

! For combination bands in particular, it is useful to notice that the
! raising and lowering operators only cause the excitation number to 
! increase or decrease by one value, around the original FC overlap
! integral.  (Clearly, the modes in their ground state are not lowered,
! which is accounted for in this code).

  ! Combine the HT term and FC term for overtones
  do r=3,excnum+1
    do j = 1, nmodes
      do iex=1,nexci
        Ls_htb(1:nexci,1:nmodes) = zzero
        Ls_htb2(1:nexci,1:nmodes) = zzero
        Ls_htb3(1:nexci,1:nmodes) = zzero
        do t=1,tstep
          dsum = zzero
          do i=1,nmodes
            dsum = dsum + (delta(iex,i)**2/two)* &
                   (one-exp(-z*freq(i)*x(t)))
          end do
          ! 4 rsums are: current mode's raised overtone (rsum, 
          ! at excitation number r), the current mode's 
          ! overtone (rsum1, at excitation number r-1), the
          ! current mode's lowered overtone (rsum2, at excitation
          ! number r-2), and the current mode's initial state
          ! raised (rsum3, at excitation number r-2)
          rsum = (delta(iex,j))**r/ &
                 (sqrt(two**r*factorial(r)))* &
                 (one-exp(-z*freq(j)*x(t)))**r
          rsum1 = (delta(iex,j))**(r-1)/ &
                  (sqrt(two**(r-1)*factorial(r-1)))* &
                  (one-exp(-z*freq(j)*x(t)))**(r-1)
          rsum2 = (delta(iex,j))**(r-2)/ &
                  (sqrt(two**(r-2)*factorial(r-2)))* &
                  (one-exp(-z*freq(j)*x(t)))**(r-2)
          ! The pattern for this is (keep in mind r starts at 3):
          ! <n|1(t)> = delta^(n-1)/sqrt(2^(n-1)n!)* &
          !            (n-delta^2+delta^2*cos(omega*t))* & 
          !            e^(-s(1-e^(-i*omega*t)-3*i*omega*t/2)
          rsum3 = (delta(iex,j))**(r-2)/ &
                  (sqrt(two**(r-2)*factorial(r-1)))* &
                  (num(r-1) - delta(iex,j)**2 + &
                  delta(iex,j)**2*cos(freq(j)*x(t)))* &
                  (one-exp(-z*freq(j)*x(t)))**(r-2)
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
          ! Calculate lineshape functions.  These assume that 
          ! the initial state is |0 0 0 ... 0 0>.
          do n=1,nmodes
            ! Integrals of the form: <F|(I+1)(t)>
            if (n == j) then
              ! Raise the initial state of the current mode,
              ! so that I+1 = 1, F = excnum.  Since the loop 
              ! starts with the current mode at excitation 
              ! level 3, you must lower the current mode twice.
              Ls_htb(iex,n) = Ls_htb(iex,n) + z*rsum3*w(t)* &
                       exp(x(t)* &
                       z*(lfreq-freq0(iex)-freq(n)) - g(iex,t) - dsum)
            else
              ! Raise the initial state of a different mode,
              ! so I+1 = 1 and F = 0.  Here, there's a product
              ! accounting for multiple modes in an excited state.
              Ls_htb(iex,n) = Ls_htb(iex,n) + &
                              z*w(t)*rsum1*expsum(n)* &
                              exp(x(t)* &
                              z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
            end if
            ! Integrals of the form: <(F-1)|I(t)>
            if (n == j) then
              ! Lower the final state of the current mode,
              ! so that F-1 = excnum-1, I = 0.  Again, you lower
              ! the current mode twice since the loop starts with 
              ! excitation level 3 (current mode at excitation 
              ! level 2).
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
              ! so that I = 0 and F+1 = excnum+1.
              Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                        z*w(t)*rootn(r)*rsum*exp(x(t)* &
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
        do k = 1, nmodes
          polot(r-1,j,1,1) = polot(r-1,j,1,1) + &
            (osc(iex,1)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
            tdipder(iex,k,1)*osc(iex,1)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
          polot(r-1,j,1,2) = polot(r-1,j,1,2) + &
            (osc(iex,1)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
            tdipder(iex,k,1)*osc(iex,2)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
          polot(r-1,j,1,3) = polot(r-1,j,1,3) + &
            (osc(iex,1)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
            tdipder(iex,k,1)*osc(iex,3)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
          polot(r-1,j,2,1) = polot(r-1,j,2,1) + &
            (osc(iex,2)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
            tdipder(iex,k,2)*osc(iex,1)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
          polot(r-1,j,2,2) = polot(r-1,j,2,2) + &
            (osc(iex,2)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
            tdipder(iex,k,2)*osc(iex,2)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
          polot(r-1,j,2,3) = polot(r-1,j,2,3) + &
            (osc(iex,2)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
            tdipder(iex,k,2)*osc(iex,3)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
          polot(r-1,j,3,1) = polot(r-1,j,3,1) + &
            (osc(iex,3)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
            tdipder(iex,k,3)*osc(iex,1)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
          polot(r-1,j,3,2) = polot(r-1,j,3,2) + &
            (osc(iex,3)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
            tdipder(iex,k,3)*osc(iex,2)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
          polot(r-1,j,3,3) = polot(r-1,j,3,3) + &
            (osc(iex,3)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
            tdipder(iex,k,3)*osc(iex,3)* &
            (Ls_htb2(iex,k)+Ls_htb3(iex,k)))* &
            pol2si
        end do
      end do
    end do
  end do
  ! Combine the HT and FC terms for comb. bands.

  do k=1,numcb
    do iex=1,nexci
      Ls_htb(1:nexci,1:nmodes) = zzero
      Ls_htb2(1:nexci,1:nmodes) = zzero
      Ls_htb3(1:nexci,1:nmodes) = zzero
      do t=1,tstep
        dsum = zzero
        do i=1,nmodes
          dsum = dsum + &
              (delta(iex,i)**2/two)* &
              (one-exp(-z*freq(i)*x(t)))
        end do
        ! Results for raising and lowering operators.
        ! sumr -  contribution after raising the final state
        ! sumr1 - contribution after lowering the final state
        ! sumr2 - contribution after raising the initial state
        sumr(1:nmodes) = one
        sumr1(1:nmodes) = one
        sumr2(1:nmodes) = one
        ! j is a comparison variable, that increases each time
        ! the loop over p goes up by 1.
        j = 1
        do p=1,nmodes
          ! Contribution from <F|(I+1)(t)>
          do r=1,nmodes
            if (r==j) then
              if (excnums(k,j) /= 0) then
                ! The pattern for this is:
                ! <n|1(t)> = delta^(n-1)/sqrt(2^(n-1)n!)* &
                !            (n-delta^2)^(n-1)* & 
                !            e^(-s(1-e^(-i*omega*t)-3*i*omega*t/2)
                sumr2(j) = sumr2(j)* &
                    (delta(iex,r))**(excnums(k,r)-1)/ &
                    sqrt(two**(excnums(k,r)-1)* &
                    factorial(excnums(k,r)))* &
                    (num(excnums(k,r)) - delta(iex,r)**2 + &
                    delta(iex,r)**2*cos(freq(r)*x(t)))* &
                    (one-exp(-z*freq(r)*x(t)))**(excnums(k,r)-1)
              else
                sumr2(j) = sumr2(j)* &
                    (delta(iex,r))**(excnums(k,r)+1)/ &
                    sqrt(two**(excnums(k,r)+1)* &
                    factorial(excnums(k,r)+1))* &
                    (one-exp(-z*freq(r)*x(t)))**(excnums(k,r)+1)
              end if
            else
              sumr2(j) = sumr2(j)* &
                  (delta(iex,r))**(excnums(k,r))/ &
                  sqrt(two**(excnums(k,r))* &
                  factorial(excnums(k,r)))* &
                  (one-exp(-z*freq(r)*x(t)))**(excnums(k,r))
            end if
          end do
          ! Contribution from <(F-1)|I(t)>
          do r=1,nmodes
            if (r==j) then
              if (excnums(k,j) /= 0) then
                sumr1(j) = sumr1(j)* &
                    (delta(iex,r))**(excnums(k,r)-1)/ &
                    sqrt(two**(excnums(k,r)-1)* &
                    factorial(excnums(k,r)-1))* &
                    (one-exp(-z*freq(r)*x(t)))**(excnums(k,r)-1)
              else
                sumr1(j) = sumr1(j)*zzero
              end if
            else
              sumr1(j) = sumr1(j)* &
                  (delta(iex,r))**(excnums(k,r))/ &
                  sqrt(two**(excnums(k,r))* &
                  factorial(excnums(k,r)))* &
                  (one-exp(-z*freq(r)*x(t)))**(excnums(k,r))
            end if
          end do
          ! Contribution from <(F+1)|I(t)>
          do r=1,nmodes
            if (r==j) then
              if (excnums(k,j) /= 0) then
                sumr(j) = sumr(j)* &
                    (delta(iex,r))**(excnums(k,r)+1)/ &
                    sqrt(two**(excnums(k,r)+1)* &
                    factorial(excnums(k,r)+1))* &
                    (one-exp(-z*freq(r)*x(t)))**(excnums(k,r)+1)
              else
                sumr(j) = sumr(j)* &
                  (delta(iex,r))**(excnums(k,r)+1)/ &
                  sqrt(two**(excnums(k,r)+1)* &
                  factorial(excnums(k,r)+1))* &
                  (one-exp(-z*freq(r)*x(t)))**(excnums(k,r)+1)
              end if
            else
              sumr(j) = sumr(j)* &
                  (delta(iex,r))**(excnums(k,r))/ &
                  sqrt(two**(excnums(k,r))* &
                  factorial(excnums(k,r)))* &
                  (one-exp(-z*freq(r)*x(t)))**(excnums(k,r))
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
          ! Integrals of the form: <(F-1)|I(t)>
          Ls_htb2(iex,i) = Ls_htb2(iex,i) + &
                      z*rootn(excnums(k,i))*sumr1(i)*w(t)* &
                      exp(x(t)* &
                      z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
          ! Integrals of the form: <(F+1)|I(t)>
          Ls_htb3(iex,i) = Ls_htb3(iex,i) + &
                      z*rootn(excnums(k,i)+1)*sumr(i)*w(t)* &
                      exp(x(t)* &
                      z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
        end do
      end do
      do r=1,nmodes
        polcb(k,1,1) = polcb(k,1,1) + &
          (osc(iex,1)*tdipder(iex,r,1)*Ls_htb(iex,r)+ &
          tdipder(iex,r,1)*osc(iex,1)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
        polcb(k,1,2) = polcb(k,1,2) + &
          (osc(iex,1)*tdipder(iex,r,2)*Ls_htb(iex,r)+ &
          tdipder(iex,r,1)*osc(iex,2)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
        polcb(k,1,3) = polcb(k,1,3) + &
          (osc(iex,1)*tdipder(iex,r,3)*Ls_htb(iex,r)+ &
          tdipder(iex,r,1)*osc(iex,3)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
        polcb(k,2,1) = polcb(k,2,1) + &
          (osc(iex,2)*tdipder(iex,r,1)*Ls_htb(iex,r)+ &
          tdipder(iex,r,2)*osc(iex,1)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
        polcb(k,2,2) = polcb(k,2,2) + &
          (osc(iex,2)*tdipder(iex,r,2)*Ls_htb(iex,r)+ &
          tdipder(iex,r,2)*osc(iex,2)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
        polcb(k,2,3) = polcb(k,2,3) + &
          (osc(iex,2)*tdipder(iex,r,3)*Ls_htb(iex,r)+ &
          tdipder(iex,r,2)*osc(iex,3)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
        polcb(k,3,1) = polcb(k,3,1) + &
          (osc(iex,3)*tdipder(iex,r,1)*Ls_htb(iex,r)+ &
          tdipder(iex,r,3)*osc(iex,1)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
        polcb(k,3,2) = polcb(k,3,2) + &
          (osc(iex,3)*tdipder(iex,r,2)*Ls_htb(iex,r)+ &
          tdipder(iex,r,3)*osc(iex,2)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
        polcb(k,3,3) = polcb(k,3,3) + &
          (osc(iex,3)*tdipder(iex,r,3)*Ls_htb(iex,r)+ &
          tdipder(iex,r,3)*osc(iex,3)*(Ls_htb2(iex,r)+ &
          Ls_htb3(iex,r)))*pol2si
      end do
    end do
  end do

end subroutine RRSHTOverCB
