subroutine R2NDHRSFCFund(latermoff, x, w, freq, delta, osc, ttpm, sechpol, &
                   freq0, lfreq, g, nexci, nmodes, tstep)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep

  complex(kindr), intent(inout) :: sechpol(nmodes,3,3,3,3)
  complex(kindr) :: Ls, rsum, dsum, z
  complex(kindr),intent(in) :: g(nexci,tstep)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci), lfreq
  real(kindr) :: x(tstep), w(tstep)
  real(kindr) :: osc(nexci,3), ttpm(nexci,10)

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: latermoff

  z = csqrt((-1,0))

! ====================================================================================
! purpose: Evaluate A term second hyperpolarizability for RSHRS spectra (fundamentals)
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
! in/out : sechpol - RSHRS second hyperpolarizability tensor
!
! ====================================================================================

  ! This is needed because some compilers don't default arrays to being zero.
  sechpol = zzero

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
        ! code.  One-photon and three-photon transition dipole moments
        ! are most likely in atomic units, where the two-photon 
        ! transition dipole moments are actually in units of 
        ! a.u.^3 / Hartree^2.  Without performing conversions, the 
        ! units of the hyperpolarizability are: a.u.^4*cm/Hartree^2 
        ! or alternatively, a.u.^4/Hartree^2*cm^{-1}.

        ! Again, this part of the code utilizes that atomic units are
        ! related to SI units for dipoles as: 
        ! 1 Cm = 1.1794744E+29 a.u.  Also, energy units are easily 
        ! converted to SI units by using 1 cm^{-1} = 1.9864456E-23 J 
        ! and 1 Hartree = 219474.63068 cm^{-1}.  Performing those 
        ! conversions allows second hyperpolarizability to be stated in SI 
        ! units: Cm^4/V^3 (since 1 V = 1 J/C). 

        ! To alleviate future confusion on my part:
        ! T_{xxx} = ttpm(#,1)
        ! T_{yyy} = ttpm(#,2)
        ! T_{zzz} = ttpm(#,3)
        ! T_{xxy} = ttpm(#,4) = T_{xyx} = T_{yxx}
        ! T_{xyy} = ttpm(#,5) = T_{yxy} = T_{yyx}
        ! T_{xxz} = ttpm(#,6) = T_{xzx} = T_{zxx}
        ! T_{xzz} = ttpm(#,7) = T_{zxz} = T_{zzx}
        ! T_{yyz} = ttpm(#,8) = T_{yzy} = T_{zyy}
        ! T_{yzz} = ttpm(#,9) = T_{zyz} = T_{zzy}
        ! T_{xyz} = ttpm(#,10) = T_{yzx} = T_{zxy} = T_{xzy} = T_{zyx} = T_{yxz}

        ! Constants used:
        ! Constants used:
        ! dipconv = 1 C*m / 1.1794744E+29 a.u. , 1 a.u. = 1 e*a_0
        ! hartree2cm = 219474.63068 cm^{-1} / 1 Hartree
        ! hartree2joules = 4.359743357E-18 J / 1 Hartree

        ! sechpol2si = (1 C*m/1.1794744E+29 a.u.)^4*
        !              (219474.63068 cm^{-1}/1 Hartree)/
        !              (4.359743357E-18 J/1 Hartree)^3
    
        ! x-direction
        ! xxxx, xxxy, xxxz
        sechpol(j,1,1,1,1) = sechpol(j,1,1,1,1) + osc(iex,1)*ttpm(iex,1)*Ls* &
          sechpol2si
        sechpol(j,1,1,1,2) = sechpol(j,1,1,1,2) + osc(iex,1)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,1,1,1,3) = sechpol(j,1,1,1,3) + osc(iex,1)*ttpm(iex,6)*Ls* &
          sechpol2si

        ! xxyx, xxyy, xxyz
        sechpol(j,1,1,2,1) = sechpol(j,1,1,2,1) + osc(iex,1)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,1,1,2,2) = sechpol(j,1,1,2,2) + osc(iex,1)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,1,1,2,3) = sechpol(j,1,1,2,3) + osc(iex,1)*ttpm(iex,10)*Ls* &
          sechpol2si

        ! xxzx, xxzy, xxzz
        sechpol(j,1,1,3,1) = sechpol(j,1,1,3,1) + osc(iex,1)*ttpm(iex,6)*Ls* &
          sechpol2si
        sechpol(j,1,1,3,2) = sechpol(j,1,1,3,2) + osc(iex,1)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,1,1,3,3) = sechpol(j,1,1,3,3) + osc(iex,1)*ttpm(iex,7)*Ls* &
          sechpol2si

        ! xyxx, xyxy, xyxz
        sechpol(j,1,2,1,1) = sechpol(j,1,2,1,1) + osc(iex,1)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,1,2,1,2) = sechpol(j,1,2,1,2) + osc(iex,1)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,1,2,1,3) = sechpol(j,1,2,1,3) + osc(iex,1)*ttpm(iex,10)*Ls* &
          sechpol2si

        ! xyyx, xyyy, xyyz
        sechpol(j,1,2,2,1) = sechpol(j,1,2,2,1) + osc(iex,1)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,1,2,2,2) = sechpol(j,1,2,2,2) + osc(iex,1)*ttpm(iex,2)*Ls* &
          sechpol2si
        sechpol(j,1,2,2,3) = sechpol(j,1,2,2,3) + osc(iex,1)*ttpm(iex,8)*Ls* &
          sechpol2si

        ! xyzx, xyzy, xyzz
        sechpol(j,1,2,3,1) = sechpol(j,1,2,3,1) + osc(iex,1)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,1,2,3,2) = sechpol(j,1,2,3,2) + osc(iex,1)*ttpm(iex,8)*Ls* &
          sechpol2si
        sechpol(j,1,2,3,3) = sechpol(j,1,2,3,3) + osc(iex,1)*ttpm(iex,9)*Ls* &
          sechpol2si

        ! xzxx, xzxy, xzxz
        sechpol(j,1,3,1,1) = sechpol(j,1,3,1,1) + osc(iex,1)*ttpm(iex,6)*Ls* &
          sechpol2si
        sechpol(j,1,3,1,2) = sechpol(j,1,3,1,2) + osc(iex,1)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,1,3,1,3) = sechpol(j,1,3,1,3) + osc(iex,1)*ttpm(iex,7)*Ls* &
          sechpol2si

        ! xzyx, xzyy, xzyz
        sechpol(j,1,3,2,1) = sechpol(j,1,3,2,1) + osc(iex,1)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,1,3,2,2) = sechpol(j,1,3,2,2) + osc(iex,1)*ttpm(iex,8)*Ls* &
          sechpol2si
        sechpol(j,1,3,2,3) = sechpol(j,1,3,2,3) + osc(iex,1)*ttpm(iex,9)*Ls* &
          sechpol2si

        ! xzzx, xzzy, xzzz
        sechpol(j,1,3,3,1) = sechpol(j,1,3,3,1) + osc(iex,1)*ttpm(iex,7)*Ls* &
          sechpol2si
        sechpol(j,1,3,3,2) = sechpol(j,1,3,3,2) + osc(iex,1)*ttpm(iex,9)*Ls* &
          sechpol2si
        sechpol(j,1,3,3,3) = sechpol(j,1,3,3,3) + osc(iex,1)*ttpm(iex,3)*Ls* &
          sechpol2si

        ! y-direction
        ! yxxx, yxxy, yxxz
        sechpol(j,2,1,1,1) = sechpol(j,2,1,1,1) + osc(iex,2)*ttpm(iex,1)*Ls* &
          sechpol2si
        sechpol(j,2,1,1,2) = sechpol(j,2,1,1,2) + osc(iex,2)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,2,1,1,3) = sechpol(j,2,1,1,3) + osc(iex,2)*ttpm(iex,6)*Ls* &
          sechpol2si

        ! yxyx, yxyy, yxyz
        sechpol(j,2,1,2,1) = sechpol(j,2,1,2,1) + osc(iex,2)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,2,1,2,2) = sechpol(j,2,1,2,2) + osc(iex,2)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,2,1,2,3) = sechpol(j,2,1,2,3) + osc(iex,2)*ttpm(iex,10)*Ls* &
          sechpol2si

        ! yxzx, yxzy, yxzz
        sechpol(j,2,1,3,1) = sechpol(j,2,1,3,1) + osc(iex,2)*ttpm(iex,6)*Ls* &
          sechpol2si
        sechpol(j,2,1,3,2) = sechpol(j,2,1,3,2) + osc(iex,2)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,2,1,3,3) = sechpol(j,2,1,3,3) + osc(iex,2)*ttpm(iex,7)*Ls* &
          sechpol2si

        ! yyxx, yyxy, yyxz
        sechpol(j,2,2,1,1) = sechpol(j,2,2,1,1) + osc(iex,2)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,2,2,1,2) = sechpol(j,2,2,1,2) + osc(iex,2)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,2,2,1,3) = sechpol(j,2,2,1,3) + osc(iex,2)*ttpm(iex,10)*Ls* &
          sechpol2si

        ! yyyx, yyyy, yyyz
        sechpol(j,2,2,2,1) = sechpol(j,2,2,2,1) + osc(iex,2)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,2,2,2,2) = sechpol(j,2,2,2,2) + osc(iex,2)*ttpm(iex,2)*Ls* &
          sechpol2si
        sechpol(j,2,2,2,3) = sechpol(j,2,2,2,3) + osc(iex,2)*ttpm(iex,8)*Ls* &
          sechpol2si

        ! yyzx, yyzy, yyzz
        sechpol(j,2,2,3,1) = sechpol(j,2,2,3,1) + osc(iex,2)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,2,2,3,2) = sechpol(j,2,2,3,2) + osc(iex,2)*ttpm(iex,8)*Ls* &
          sechpol2si
        sechpol(j,2,2,3,3) = sechpol(j,2,2,3,3) + osc(iex,2)*ttpm(iex,9)*Ls* &
          sechpol2si

        ! yzxx, yzxy, yzxz
        sechpol(j,2,3,1,1) = sechpol(j,2,3,1,1) + osc(iex,2)*ttpm(iex,6)*Ls* &
          sechpol2si
        sechpol(j,2,3,1,2) = sechpol(j,2,3,1,2) + osc(iex,2)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,2,3,1,3) = sechpol(j,2,3,1,3) + osc(iex,2)*ttpm(iex,7)*Ls* &
          sechpol2si

        ! yzyx, yzyy, yzyz
        sechpol(j,2,3,2,1) = sechpol(j,2,3,2,1) + osc(iex,2)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,2,3,2,2) = sechpol(j,2,3,2,2) + osc(iex,2)*ttpm(iex,8)*Ls* &
          sechpol2si
        sechpol(j,2,3,2,3) = sechpol(j,2,3,2,3) + osc(iex,2)*ttpm(iex,9)*Ls* &
          sechpol2si

        ! yzzx, yzzy, yzzz
        sechpol(j,2,3,3,1) = sechpol(j,2,3,3,1) + osc(iex,2)*ttpm(iex,7)*Ls* &
          sechpol2si
        sechpol(j,2,3,3,2) = sechpol(j,2,3,3,2) + osc(iex,2)*ttpm(iex,9)*Ls* &
          sechpol2si
        sechpol(j,2,3,3,3) = sechpol(j,2,3,3,3) + osc(iex,2)*ttpm(iex,3)*Ls* &
          sechpol2si

        ! z-direction
        ! zxxx, zxxy, zxxz
        sechpol(j,3,1,1,1) = sechpol(j,3,1,1,1) + osc(iex,3)*ttpm(iex,1)*Ls* &
          sechpol2si
        sechpol(j,3,1,1,2) = sechpol(j,3,1,1,2) + osc(iex,3)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,3,1,1,3) = sechpol(j,3,1,1,3) + osc(iex,3)*ttpm(iex,6)*Ls* &
          sechpol2si

        ! zxyx, zxyy, zxyz
        sechpol(j,3,1,2,1) = sechpol(j,3,1,2,1) + osc(iex,3)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,3,1,2,2) = sechpol(j,3,1,2,2) + osc(iex,3)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,3,1,2,3) = sechpol(j,3,1,2,3) + osc(iex,3)*ttpm(iex,10)*Ls* &
          sechpol2si

        ! zxzx, zxzy, zxzz
        sechpol(j,3,1,3,1) = sechpol(j,3,1,3,1) + osc(iex,3)*ttpm(iex,6)*Ls* &
          sechpol2si
        sechpol(j,3,1,3,2) = sechpol(j,3,1,3,2) + osc(iex,3)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,3,1,3,3) = sechpol(j,3,1,3,3) + osc(iex,3)*ttpm(iex,7)*Ls* &
          sechpol2si

        ! zyxx, zyxy, zyxz
        sechpol(j,3,2,1,1) = sechpol(j,3,2,1,1) + osc(iex,3)*ttpm(iex,4)*Ls* &
          sechpol2si
        sechpol(j,3,2,1,2) = sechpol(j,3,2,1,2) + osc(iex,3)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,3,2,1,3) = sechpol(j,3,2,1,3) + osc(iex,3)*ttpm(iex,10)*Ls* &
          sechpol2si

        ! zyyx, zyyy, zyyz
        sechpol(j,3,2,2,1) = sechpol(j,3,2,2,1) + osc(iex,3)*ttpm(iex,5)*Ls* &
          sechpol2si
        sechpol(j,3,2,2,2) = sechpol(j,3,2,2,2) + osc(iex,3)*ttpm(iex,2)*Ls* &
          sechpol2si
        sechpol(j,3,2,2,3) = sechpol(j,3,2,2,3) + osc(iex,3)*ttpm(iex,8)*Ls* &
          sechpol2si

        ! zyzx, zyzy, zyzz
        sechpol(j,3,2,3,1) = sechpol(j,3,2,3,1) + osc(iex,3)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,3,2,3,2) = sechpol(j,3,2,3,2) + osc(iex,3)*ttpm(iex,8)*Ls* &
          sechpol2si
        sechpol(j,3,2,3,3) = sechpol(j,3,2,3,3) + osc(iex,3)*ttpm(iex,9)*Ls* &
          sechpol2si

        ! zzxx, zzxy, zzxz
        sechpol(j,3,3,1,1) = sechpol(j,3,3,1,1) + osc(iex,3)*ttpm(iex,6)*Ls* &
          sechpol2si
        sechpol(j,3,3,1,2) = sechpol(j,3,3,1,2) + osc(iex,3)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,3,3,1,3) = sechpol(j,3,3,1,3) + osc(iex,3)*ttpm(iex,7)*Ls* &
          sechpol2si

        ! zzyx, zzyy, zzyz
        sechpol(j,3,3,2,1) = sechpol(j,3,3,2,1) + osc(iex,3)*ttpm(iex,10)*Ls* &
          sechpol2si
        sechpol(j,3,3,2,2) = sechpol(j,3,3,2,2) + osc(iex,3)*ttpm(iex,8)*Ls* &
          sechpol2si
        sechpol(j,3,3,2,3) = sechpol(j,3,3,2,3) + osc(iex,3)*ttpm(iex,9)*Ls* &
          sechpol2si

        ! zzzx, zzzy, zzzz
        sechpol(j,3,3,3,1) = sechpol(j,3,3,3,1) + osc(iex,3)*ttpm(iex,7)*Ls* &
          sechpol2si
        sechpol(j,3,3,3,2) = sechpol(j,3,3,3,2) + osc(iex,3)*ttpm(iex,9)*Ls* &
          sechpol2si
        sechpol(j,3,3,3,3) = sechpol(j,3,3,3,3) + osc(iex,3)*ttpm(iex,3)*Ls* &
          sechpol2si

      end do

    end do
  end if

end subroutine R2NDHRSFCFund
