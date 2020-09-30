Subroutine RVROAFCFund(latermoff, x, w, freq, delta, osc, mtdip, tquad, &
                     pol, polg, polgs, pola, polas, freq0, lfreq, g, &
                     nexci, nmodes, tstep)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep

  complex(kindr), intent(inout) :: pol(nmodes,3,3)
  complex(kindr), intent(inout) :: polg(nmodes,3,3)
  complex(kindr), intent(inout) :: polgs(nmodes,3,3)
  complex(kindr), intent(inout) :: pola(nmodes,3,3,3)
  complex(kindr), intent(inout) :: polas(nmodes,3,3,3)
  complex(kindr) :: Ls, rsum, dsum, z
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep) 
  real(kindr) :: osc(nexci,3) 
  real(kindr) :: mtdip(nexci,3) 
  real(kindr) :: tquad(nexci,6) 

  complex(kindr) :: zzero=(0d0,0d0) 

  logical :: latermoff

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate A term tensors for RVROA spectra (fundamentals)
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
! in/out : pol - RVROA electric dipole-electric dipole polarizability tensor
!          polg - RVROA electric dipole-magnetic dipole tensor
!          polgs - RVROA magnetic dipole-electric dipole tensor
!          pola - RVROA electric dipole-electric quadrupole tensor
!          polas - RVROA electric quadrupole-electric dipole tensor
!
! =============================================================================

! For information on the equations used, see:
!
! Luber, S.; Neugebauer, J.; and Reiher, M.  J. Chem. Phys. 132, 044113 (2010).

  ! This is needed because some compilers don't default arrays to being zero.
  pol = zzero
  polg = zzero
  polgs = zzero
  pola = zzero
  polas = zzero

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

            rsum = delta(iex,j)/sqrt(two)*(one-exp(-z*freq(j)*x(t)))
            Ls = Ls + z*w(t)*rsum*exp(x(t)* &
              z*(lfreq-freq0(iex)) - g(iex,t) - dsum)
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

        ! Electric dipole-electric dipole polarizability

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

        ! Electric dipole-magnetic dipole polarizability (figure out
        ! how to convert to SI units)
        polg(j,1,1) = polg(j,1,1) + osc(iex,1)*mtdip(iex,1)*Ls
        polg(j,1,2) = polg(j,1,2) + osc(iex,1)*mtdip(iex,2)*Ls
        polg(j,1,3) = polg(j,1,3) + osc(iex,1)*mtdip(iex,3)*Ls
        polg(j,2,1) = polg(j,2,1) + osc(iex,2)*mtdip(iex,1)*Ls
        polg(j,2,2) = polg(j,2,2) + osc(iex,2)*mtdip(iex,2)*Ls
        polg(j,2,3) = polg(j,2,3) + osc(iex,2)*mtdip(iex,3)*Ls
        polg(j,3,1) = polg(j,3,1) + osc(iex,3)*mtdip(iex,1)*Ls
        polg(j,3,2) = polg(j,3,2) + osc(iex,3)*mtdip(iex,2)*Ls
        polg(j,3,3) = polg(j,3,3) + osc(iex,3)*mtdip(iex,3)*Ls

        ! Magnetic dipole-electric dipole polarizability (figure out
        ! how to convert to SI units)
        polgs(j,1,1) = polgs(j,1,1) + mtdip(iex,1)*osc(iex,1)*Ls
        polgs(j,1,2) = polgs(j,1,2) + mtdip(iex,1)*osc(iex,2)*Ls
        polgs(j,1,3) = polgs(j,1,3) + mtdip(iex,1)*osc(iex,3)*Ls
        polgs(j,2,1) = polgs(j,2,1) + mtdip(iex,2)*osc(iex,1)*Ls
        polgs(j,2,2) = polgs(j,2,2) + mtdip(iex,2)*osc(iex,2)*Ls
        polgs(j,2,3) = polgs(j,2,3) + mtdip(iex,2)*osc(iex,3)*Ls
        polgs(j,3,1) = polgs(j,3,1) + mtdip(iex,3)*osc(iex,1)*Ls
        polgs(j,3,2) = polgs(j,3,2) + mtdip(iex,3)*osc(iex,2)*Ls
        polgs(j,3,3) = polgs(j,3,3) + mtdip(iex,3)*osc(iex,3)*Ls

        ! Electric dipole-electric quadupole polarizability (figure out
        ! how to convert to SI units)
        ! xx = 1 yy = 2 zz = 3 xy = 4 xz = 5 yz = 6
        pola(j,1,1,1) = pola(j,1,1,1) + osc(iex,1)*tquad(iex,1)*Ls
        pola(j,1,1,2) = pola(j,1,1,2) + osc(iex,1)*tquad(iex,4)*Ls
        pola(j,1,1,3) = pola(j,1,1,3) + osc(iex,1)*tquad(iex,5)*Ls

        pola(j,1,2,1) = pola(j,1,2,1) + osc(iex,1)*tquad(iex,4)*Ls
        pola(j,1,2,2) = pola(j,1,2,2) + osc(iex,1)*tquad(iex,2)*Ls
        pola(j,1,2,3) = pola(j,1,2,3) + osc(iex,1)*tquad(iex,6)*Ls
        
        pola(j,1,3,1) = pola(j,1,3,1) + osc(iex,1)*tquad(iex,5)*Ls
        pola(j,1,3,2) = pola(j,1,3,2) + osc(iex,1)*tquad(iex,6)*Ls
        pola(j,1,3,3) = pola(j,1,3,3) + osc(iex,1)*tquad(iex,3)*Ls

        pola(j,2,1,1) = pola(j,2,1,1) + osc(iex,2)*tquad(iex,1)*Ls
        pola(j,2,1,2) = pola(j,2,1,2) + osc(iex,2)*tquad(iex,4)*Ls
        pola(j,2,1,3) = pola(j,2,1,3) + osc(iex,2)*tquad(iex,5)*Ls

        pola(j,2,2,1) = pola(j,2,2,1) + osc(iex,2)*tquad(iex,4)*Ls
        pola(j,2,2,2) = pola(j,2,2,2) + osc(iex,2)*tquad(iex,2)*Ls
        pola(j,2,2,3) = pola(j,2,2,3) + osc(iex,2)*tquad(iex,6)*Ls

        pola(j,2,3,1) = pola(j,2,3,1) + osc(iex,2)*tquad(iex,5)*Ls
        pola(j,2,3,2) = pola(j,2,3,2) + osc(iex,2)*tquad(iex,6)*Ls
        pola(j,2,3,3) = pola(j,2,3,3) + osc(iex,2)*tquad(iex,3)*Ls

        pola(j,3,1,1) = pola(j,3,1,1) + osc(iex,3)*tquad(iex,1)*Ls
        pola(j,3,1,2) = pola(j,3,1,2) + osc(iex,3)*tquad(iex,4)*Ls
        pola(j,3,1,3) = pola(j,3,1,3) + osc(iex,3)*tquad(iex,5)*Ls

        pola(j,3,2,1) = pola(j,3,2,1) + osc(iex,3)*tquad(iex,4)*Ls
        pola(j,3,2,2) = pola(j,3,2,2) + osc(iex,3)*tquad(iex,2)*Ls
        pola(j,3,2,3) = pola(j,3,2,3) + osc(iex,3)*tquad(iex,6)*Ls

        pola(j,3,3,1) = pola(j,3,3,1) + osc(iex,3)*tquad(iex,5)*Ls
        pola(j,3,3,2) = pola(j,3,3,2) + osc(iex,3)*tquad(iex,6)*Ls
        pola(j,3,3,3) = pola(j,3,3,3) + osc(iex,3)*tquad(iex,3)*Ls

        ! Electric quadrupole-electric dipole polarizability (figure out
        ! how to convert to SI units)
        polas(j,1,1,1) = polas(j,1,1,1) + osc(iex,1)*tquad(iex,1)*Ls
        polas(j,1,1,2) = polas(j,1,1,2) + osc(iex,1)*tquad(iex,4)*Ls
        polas(j,1,1,3) = polas(j,1,1,3) + osc(iex,1)*tquad(iex,5)*Ls

        polas(j,1,2,1) = polas(j,1,2,1) + osc(iex,1)*tquad(iex,4)*Ls
        polas(j,1,2,2) = polas(j,1,2,2) + osc(iex,1)*tquad(iex,2)*Ls
        polas(j,1,2,3) = polas(j,1,2,3) + osc(iex,1)*tquad(iex,6)*Ls
        
        polas(j,1,3,1) = polas(j,1,3,1) + osc(iex,1)*tquad(iex,5)*Ls
        polas(j,1,3,2) = polas(j,1,3,2) + osc(iex,1)*tquad(iex,6)*Ls
        polas(j,1,3,3) = polas(j,1,3,3) + osc(iex,1)*tquad(iex,3)*Ls

        polas(j,2,1,1) = polas(j,2,1,1) + osc(iex,2)*tquad(iex,1)*Ls
        polas(j,2,1,2) = polas(j,2,1,2) + osc(iex,2)*tquad(iex,4)*Ls
        polas(j,2,1,3) = polas(j,2,1,3) + osc(iex,2)*tquad(iex,5)*Ls

        polas(j,2,2,1) = polas(j,2,2,1) + osc(iex,2)*tquad(iex,4)*Ls
        polas(j,2,2,2) = polas(j,2,2,2) + osc(iex,2)*tquad(iex,2)*Ls
        polas(j,2,2,3) = polas(j,2,2,3) + osc(iex,2)*tquad(iex,6)*Ls

        polas(j,2,3,1) = polas(j,2,3,1) + osc(iex,2)*tquad(iex,5)*Ls
        polas(j,2,3,2) = polas(j,2,3,2) + osc(iex,2)*tquad(iex,6)*Ls
        polas(j,2,3,3) = polas(j,2,3,3) + osc(iex,2)*tquad(iex,3)*Ls

        polas(j,3,1,1) = polas(j,3,1,1) + osc(iex,3)*tquad(iex,1)*Ls
        polas(j,3,1,2) = polas(j,3,1,2) + osc(iex,3)*tquad(iex,4)*Ls
        polas(j,3,1,3) = polas(j,3,1,3) + osc(iex,3)*tquad(iex,5)*Ls

        polas(j,3,2,1) = polas(j,3,2,1) + osc(iex,3)*tquad(iex,4)*Ls
        polas(j,3,2,2) = polas(j,3,2,2) + osc(iex,3)*tquad(iex,2)*Ls
        polas(j,3,2,3) = polas(j,3,2,3) + osc(iex,3)*tquad(iex,6)*Ls

        polas(j,3,3,1) = polas(j,3,3,1) + osc(iex,3)*tquad(iex,5)*Ls
        polas(j,3,3,2) = polas(j,3,3,2) + osc(iex,3)*tquad(iex,6)*Ls
        polas(j,3,3,3) = polas(j,3,3,3) + osc(iex,3)*tquad(iex,3)*Ls
      end do
    end do
  end if

end subroutine RVROAFCFund

