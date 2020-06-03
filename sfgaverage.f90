subroutine SFGFundAvg(hpol, width, psi, theta, diroption, rcrs, nmodes)

  use constants

  implicit none

  integer :: i, nmodes

  complex(kindr), intent(inout) :: hpol(nmodes,3,3,3)
  real(kindr), intent(in) :: width

  real(kindr), intent(inout) :: rcrs(nmodes)


  complex(kindr) :: tmpcrs(nmodes)
  real(kindr),intent(in) :: theta
  real(kindr),intent(in) :: psi
  integer,intent(in)     :: diroption
  complex(kindr) :: z

! ================================================================================
! purpose: Orientational averaging of the SFG hyperpolarizability for fundamentals
!
! input  : hpol - SFG hyperpolarizability tensor
!          width - lifetime of the IR transition (< 30 fs)
!          psi - Euler angle
!          theta - Euler angle
!          diroption - what polarization is used for the light?
!          nmodes - number of vibrational normal modes
!
! out : rcrs - SFG scattering factor 
!
! ================================================================================


  z = csqrt((-1,0))

  ! Include the lifetime for the IR transition and also convert the
  ! units completely to a.u.  The hyperpolarizability in a.u. is 
  ! defined as: (e*a0)^3/Eh^2, where e is the fundamental charge,
  ! a0 is the Bohr radius, and Eh is a Hartree.  We then convert
  ! to "alternative SI" units, where the hyperpolarizability is
  ! converted to SI units (C^3 m^3/J^2) and divided by the
  ! vacuum permittivity (in C^2 m^{-1} J^{-1}).  This leaves
  ! hyperpolarizabilities in m^4 V^{-1}.

  ! For specifics on units, see:
  ! Shelton, D. P. and Rice, J. E.  Chem. Rev., 94, 3 (1994).

  ! au2altsi = 3.621294e-42 m^4/V / 1 e^3*a0^3/Eh^2
  ! hartree2cm = 219474.63068 cm^{-1} / 1 Eh 

  hpol = hpol*au2altsi*hartree2cm**2 / (z*width)

  do i = 1, nmodes

     if (diroption == 1 ) then
     ! DO CHI2_ZZZ
     ! PPP - polarization

     tmpcrs(i) = (cos(theta)**3) * hpol(i,3,3,3) &
             + sin(theta) * sin(psi) * ( hpol(i,2,3,3) + hpol(i,3,2,3) + hpol(i,3,3,2) ) &
             - sin(theta) * cos(psi) * ( hpol(i,1,3,3) + hpol(i,3,1,3) + hpol(i,3,3,1) )  &
             + (sin(theta)**2) * cos(theta) * (sin(psi)**2) * ( hpol(i,2,2,3) + hpol(i,2,3,2) + hpol(i,3,2,2) ) &
             + (sin(theta)**2) * cos(theta) * (cos(psi)**2) * ( hpol(i,1,1,3) + hpol(i,1,3,1) + hpol(i,3,1,1) )  &
             - (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * ( hpol(i,1,2,3) + hpol(i,1,3,2) &
                                                                    + hpol(i,2,1,3) + hpol(i,2,3,1) &
                                                                    + hpol(i,3,1,2) + hpol(i,3,2,1) ) &
             + (sin(theta)**3) * sin(psi) * ( hpol(i,1,1,2) + hpol(i,1,2,1) + hpol(i,2,1,1) &
                                            - hpol(i,2,3,3) - hpol(i,3,2,3) - hpol(i,3,3,2) ) &
             + (sin(theta)**3) * cos(psi) * ( hpol(i,1,3,3) + hpol(i,3,1,3) + hpol(i,3,3,1) &
                                            - hpol(i,1,2,2) - hpol(i,2,1,2) - hpol(i,2,2,1) ) &
             + (sin(theta)**3) * (sin(psi)**3) * ( hpol(i,2,2,2) - hpol(i,1,1,2) - hpol(i,1,2,1) - hpol(i,2,1,1) ) & 
             + (sin(theta)**3) * (cos(psi)**3) * ( -hpol(i,1,1,1) + hpol(i,1,2,2) + hpol(i,2,1,2) + hpol(i,2,2,1) )

     elseif (diroption == 2) then
     ! DO CHI2_ZXX
     ! PSS - Polarization

     tmpcrs(i) = (sin(theta)**2) * cos(theta)  * hpol(i,3,3,3) & 
             + cos(theta) * ( hpol(i,3,1,1) + hpol(i,3,2,2) ) &
             - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * ( hpol(i,2,2,3) + hpol(i,2,3,2) + hpol(i,3,2,2) ) & 
             - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * ( hpol(i,1,1,3) + hpol(i,1,3,1) + hpol(i,3,1,1) ) &
             + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * ( hpol(i,1,2,3) + hpol(i,1,3,2) &
                                                                    + hpol(i,2,1,3) + hpol(i,2,3,1) &
                                                                    + hpol(i,3,1,2) + hpol(i,3,2,1) ) &
             + sin(theta) * sin(psi) * ( hpol(i,2,2,2) + hpol(i,2,1,1) - hpol(i,3,2,3) - hpol(i,3,3,2) ) &
             + sin(theta) * cos(psi) * ( -hpol(i,1,1,1) - hpol(i,1,2,2) + hpol(i,3,1,3) + hpol(i,3,3,1) ) &
             + (sin(theta)**3) * sin(psi) * ( -hpol(i,1,1,2) - hpol(i,1,2,1) - hpol(i,2,1,1) &
                                            + hpol(i,2,3,3) + hpol(i,3,2,3) + hpol(i,3,3,2) ) &
             + (sin(theta)**3) * cos(psi) * ( hpol(i,1,2,2) + hpol(i,2,1,2) + hpol(i,2,2,1) &
                                            - hpol(i,1,3,3) - hpol(i,3,1,3) - hpol(i,3,3,1) ) &
             + (sin(theta)**3) * (sin(psi)**3) * ( -hpol(i,2,2,2) + hpol(i,1,1,2) + hpol(i,1,2,1) + hpol(i,2,1,1) ) &
             + (sin(theta)**3) * (cos(psi)**3) * ( hpol(i,1,1,1) - hpol(i,1,2,2) - hpol(i,2,1,2) - hpol(i,2,2,1) ) 

     tmpcrs(i) = tmpcrs(i) * f12

     elseif (diroption == 3) then
     ! DO CHI2_XZX
     ! SPS - Polarization

     tmpcrs(i) = (sin(theta)**2) * cos(theta)  * hpol(i,3,3,3) &
             + cos(theta) * ( hpol(i,1,3,1) + hpol(i,2,3,2) ) &
             - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * ( hpol(i,2,2,3) + hpol(i,2,3,2) + hpol(i,3,2,2) ) &
             - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * ( hpol(i,1,1,3) + hpol(i,1,3,1) + hpol(i,3,1,1) ) &
             + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * ( hpol(i,1,2,3) + hpol(i,1,3,2) + hpol(i,2,1,3) &
                                                                    + hpol(i,2,3,1) + hpol(i,3,1,2) + hpol(i,3,2,1) ) &
             + sin(theta) * sin(psi) * ( hpol(i,2,2,2) + hpol(i,1,2,1) - hpol(i,2,3,3) - hpol(i,3,3,2) ) &
             + sin(theta) * cos(psi) * ( -hpol(i,1,1,1) - hpol(i,2,1,2) + hpol(i,1,3,3) + hpol(i,3,3,1) ) &
             + (sin(theta)**3) * sin(psi) * ( -hpol(i,1,1,2) - hpol(i,1,2,1) - hpol(i,2,1,1) &
                                            + hpol(i,2,3,3) + hpol(i,3,2,3) + hpol(i,3,3,2) ) &
             + (sin(theta)**3) * cos(psi) * ( hpol(i,1,2,2) + hpol(i,2,1,2) + hpol(i,2,2,1) &
                                            - hpol(i,1,3,3) - hpol(i,3,1,3) - hpol(i,3,3,1) ) &
             + (sin(theta)**3) * (sin(psi)**3) * ( -hpol(i,2,2,2) + hpol(i,1,1,2) + hpol(i,1,2,1) + hpol(i,2,1,1) ) &
             + (sin(theta)**3) * (cos(psi)**3) * ( hpol(i,1,1,1) - hpol(i,1,2,2) - hpol(i,2,1,2) - hpol(i,2,2,1) )

     tmpcrs(i) = tmpcrs(i) * f12
   
     elseif (diroption == 4) then
     ! DO  CHI2_XXZ
     ! SSP - polarization

     tmpcrs(i) = (sin(theta)**2) * cos(theta)  * hpol(i,3,3,3) &
             + cos(theta) * ( hpol(i,1,1,3) + hpol(i,2,2,3) ) &
             - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * ( hpol(i,2,2,3) + hpol(i,2,3,2) + hpol(i,3,2,2) ) &
             - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * ( hpol(i,1,1,3) + hpol(i,1,3,1) + hpol(i,3,1,1) ) &
             + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * ( hpol(i,1,2,3) + hpol(i,1,3,2) + hpol(i,2,1,3) &
                                                                    + hpol(i,2,3,1) + hpol(i,3,1,2) + hpol(i,3,2,1) ) &
             + sin(theta) * sin(psi) * ( hpol(i,2,2,2) + hpol(i,1,1,2) - hpol(i,2,3,3) - hpol(i,3,2,3) ) &
             + sin(theta) * cos(psi) * ( -hpol(i,1,1,1) - hpol(i,2,2,1) + hpol(i,1,3,3) + hpol(i,3,1,3) ) &
             + (sin(theta)**3) * sin(psi) * ( -hpol(i,1,1,2) - hpol(i,1,2,1) - hpol(i,2,1,1) &
                                            + hpol(i,2,3,3) + hpol(i,3,2,3) + hpol(i,3,3,2) ) &
             + (sin(theta)**3) * cos(psi) * ( hpol(i,1,2,2) + hpol(i,2,1,2) + hpol(i,2,2,1) &
                                            - hpol(i,1,3,3) - hpol(i,3,1,3) - hpol(i,3,3,1)) &
             + (sin(theta)**3) * (sin(psi)**3) * ( -hpol(i,2,2,2) + hpol(i,1,1,2) + hpol(i,1,2,1) + hpol(i,2,1,1) ) &
             + (sin(theta)**3) * (cos(psi)**3) * ( hpol(i,1,1,1) - hpol(i,1,2,2) - hpol(i,2,1,2) - hpol(i,2,2,1) ) 

     tmpcrs(i) = tmpcrs(i) * f12 

     elseif (diroption == 5 ) then
     ! DO CHI2_XYZ
     ! SSP Polarization

     tmpcrs(i) = (cos(theta)**2) * ( hpol(i,1,2,3) - hpol(i,2,1,3) ) &
             + (sin(theta)**2) * (sin(psi)**2) * ( hpol(i,3,1,2) - hpol(i,1,3,2) ) &
             + (sin(theta)**2) * (cos(psi)**2) * ( hpol(i,2,3,1) - hpol(i,3,2,1) ) &
             + (sin(theta)**2) * sin(psi) * cos(psi) * ( hpol(i,1,3,1) - hpol(i,2,3,2) &
                                                       - hpol(i,3,1,1) + hpol(i,3,2,2) ) &
             + sin(theta) * cos(theta) * sin(psi) * ( hpol(i,1,2,2) - hpol(i,1,3,3) &
                                                    - hpol(i,2,1,2) + hpol(i,3,1,3) ) &
             + sin(theta) * cos(theta) * cos(psi) * ( -hpol(i,1,2,1) + hpol(i,2,1,1) &
                                                    - hpol(i,2,3,3) + hpol(i,3,2,3) )

     tmpcrs(i) = tmpcrs(i) * f12 

     elseif (diroption == 6) then
     ! DO CHI2_XZY
     ! SPS Polarization

     tmpcrs(i) = (cos(theta)**2) * ( hpol(i,1,3,2) - hpol(i,2,3,1) ) &
             + (sin(theta)**2) * (sin(psi)**2) * ( hpol(i,3,2,1) - hpol(i,1,2,3) ) &
             + (sin(theta)**2) * (cos(psi)**2) * ( hpol(i,2,1,3) - hpol(i,3,1,2) ) &
             + (sin(theta)**2) * sin(psi) * cos(psi) * ( hpol(i,1,1,3) - hpol(i,2,2,3) &
                                                       - hpol(i,3,1,1) + hpol(i,3,2,2) ) &
             + sin(theta) * cos(theta) * sin(psi) * ( hpol(i,1,2,2) - hpol(i,1,3,3) &
                                                    - hpol(i,2,2,1) + hpol(i,3,3,1) ) & 
             + sin(theta) * cos(theta) * cos(psi) * ( -hpol(i,1,1,2) + hpol(i,2,1,1) &
                                                    - hpol(i,2,3,3) + hpol(i,3,3,2) ) 

     tmpcrs(i) = tmpcrs(i) * f12

     elseif (diroption == 7) then
     ! DO CHI2_ZXY
     ! PSS Polarization

     tmpcrs(i) = (cos(theta)**2) * ( hpol(i,3,1,2) - hpol(i,3,2,1) ) &
             + (sin(theta)**2) * (sin(psi)**2) * ( hpol(i,2,3,1) - hpol(i,2,1,3) ) &
             + (sin(theta)**2) * (cos(psi)**2) * ( hpol(i,1,2,3) - hpol(i,1,3,2) ) &
             + (sin(theta)**2) * sin(psi) * cos(psi) * ( hpol(i,1,1,3) - hpol(i,1,3,1) &
                                                       - hpol(i,2,2,3) + hpol(i,2,3,2) ) &
             + sin(theta) * cos(theta) * sin(psi) * ( hpol(i,2,1,2) - hpol(i,2,2,1) &
                                                    - hpol(i,3,1,3) + hpol(i,3,3,1) ) &
             + sin(theta) * cos(theta) * cos(psi) * ( -hpol(i,1,1,2) + hpol(i,1,2,1) &
                                                    - hpol(i,3,2,3) + hpol(i,3,3,2) )

     tmpcrs(i) = tmpcrs(i) * f12

     else
       stop 
     end if

   end do
      
   rcrs = tmpcrs*conjg(tmpcrs)
    
end subroutine SFGFundAvg

subroutine SFGFullAvg(nhpol, psi, theta, diroption, rcrs_fsfg, freqend)

  use constants

  implicit none

  integer :: i, freqend

  complex(kindr), intent(inout) :: nhpol(freqend,3,3,3)

  real(kindr), intent(out) :: rcrs_fsfg(freqend)


  complex(kindr) :: tmpcrs(freqend)
  real(kindr),intent(in) :: theta
  real(kindr),intent(in) :: psi
  integer,intent(in)     :: diroption
  complex(kindr) :: z

! ================================================================================
! purpose: Orientational averaging of the SFG hyperpolarizability for
!          fundamentals in the SFG full spectrum
!
! input  : nhpol - SFG hyperpolarizability tensor for full range
!          psi - Euler angle
!          theta - Euler angle
!          diroption - what polarization is used for the light?
!          freqend - number of points in the SFG full spectrum
!
! out : rcrs_fsfg - SFG scattering factor 
!
! ================================================================================


  z = csqrt((-1,0))

  ! Include the lifetime for the IR transition and also convert the
  ! units completely to a.u.  The hyperpolarizability in a.u. is 
  ! defined as: (e*a0)^3/Eh^2, where e is the fundamental charge,
  ! a0 is the Bohr radius, and Eh is a Hartree.  We then convert
  ! to "alternative SI" units, where the hyperpolarizability is
  ! converted to SI units (C^3 m^3/J^2) and divided by the
  ! vacuum permittivity (in C^2 m^{-1} J^{-1}).  This leaves
  ! hyperpolarizabilities in m^4 V^{-1}.

  ! For specifics on units, see:
  ! Shelton, D. P. and Rice, J. E.  Chem. Rev., 94, 3 (1994).

  ! au2altsi = 3.621294e-42 m^4/V / 1 e^3*a0^3/Eh^2
  ! hartree2cm = 219474.63068 cm^{-1} / 1 Eh 
  nhpol = nhpol*au2altsi*hartree2cm**2 

  do i = 1, freqend
     if (diroption == 1 ) then
     ! DO CHI2_ZZZ
     ! PPP - polarization

     tmpcrs(i) = (cos(theta)**3) * nhpol(i,3,3,3) &
             + sin(theta) * sin(psi) * ( nhpol(i,2,3,3) + nhpol(i,3,2,3) + nhpol(i,3,3,2) ) &
             - sin(theta) * cos(psi) * ( nhpol(i,1,3,3) + nhpol(i,3,1,3) + nhpol(i,3,3,1) ) &
             + (sin(theta)**2) * cos(theta) * (sin(psi)**2) * ( nhpol(i,2,2,3) + nhpol(i,2,3,2) + nhpol(i,3,2,2) ) &
             + (sin(theta)**2) * cos(theta) * (cos(psi)**2) * ( nhpol(i,1,1,3) + nhpol(i,1,3,1) + nhpol(i,3,1,1) ) &
             - (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * ( nhpol(i,1,2,3) + nhpol(i,1,3,2) + nhpol(i,2,1,3) &
                                                                    + nhpol(i,2,3,1) + nhpol(i,3,1,2) + nhpol(i,3,2,1) ) &
             + (sin(theta)**3) * sin(psi) * ( nhpol(i,1,1,2) + nhpol(i,1,2,1) + nhpol(i,2,1,1) &
                                            - nhpol(i,2,3,3) - nhpol(i,3,2,3) - nhpol(i,3,3,2) ) &
             + (sin(theta)**3) * cos(psi) * ( nhpol(i,1,3,3) + nhpol(i,3,1,3) + nhpol(i,3,3,1) &
                                            - nhpol(i,1,2,2) - nhpol(i,2,1,2) - nhpol(i,2,2,1) ) &
             + (sin(theta)**3) * (sin(psi)**3) * ( nhpol(i,2,2,2) - nhpol(i,1,1,2) - nhpol(i,1,2,1) - nhpol(i,2,1,1) ) & 
             + (sin(theta)**3) * (cos(psi)**3) * ( -nhpol(i,1,1,1) + nhpol(i,1,2,2) + nhpol(i,2,1,2) + nhpol(i,2,2,1) )

     elseif (diroption == 2) then
     ! DO CHI2_ZXX
     ! PSS - Polarization

     tmpcrs(i) = (sin(theta)**2) * cos(theta)  * nhpol(i,3,3,3) & 
             + cos(theta) * ( nhpol(i,3,1,1) + nhpol(i,3,2,2) ) &
             - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * ( nhpol(i,2,2,3) + nhpol(i,2,3,2) + nhpol(i,3,2,2) ) & 
             - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * ( nhpol(i,1,1,3) + nhpol(i,1,3,1) + nhpol(i,3,1,1) ) &
             + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * ( nhpol(i,1,2,3) + nhpol(i,1,3,2) + nhpol(i,2,1,3) &
                                                                    + nhpol(i,2,3,1) + nhpol(i,3,1,2) + nhpol(i,3,2,1) ) &
             + sin(theta) * sin(psi) * ( nhpol(i,2,2,2) + nhpol(i,2,1,1) - nhpol(i,3,2,3) - nhpol(i,3,3,2) ) &
             + sin(theta) * cos(psi) * ( -nhpol(i,1,1,1) - nhpol(i,1,2,2) + nhpol(i,3,1,3) + nhpol(i,3,3,1) ) &
             + (sin(theta)**3) * sin(psi) * ( -nhpol(i,1,1,2) - nhpol(i,1,2,1) - nhpol(i,2,1,1) &
                                            + nhpol(i,2,3,3) + nhpol(i,3,2,3) + nhpol(i,3,3,2) ) &   
             + (sin(theta)**3) * cos(psi) * ( nhpol(i,1,2,2) + nhpol(i,2,1,2) + nhpol(i,2,2,1) &
                                            - nhpol(i,1,3,3) - nhpol(i,3,1,3) - nhpol(i,3,3,1) ) &
             + (sin(theta)**3) * (sin(psi)**3) * ( -nhpol(i,2,2,2) + nhpol(i,1,1,2) + nhpol(i,1,2,1) + nhpol(i,2,1,1) ) &
             + (sin(theta)**3) * (cos(psi)**3) * ( nhpol(i,1,1,1) - nhpol(i,1,2,2) - nhpol(i,2,1,2) - nhpol(i,2,2,1) ) 

     tmpcrs(i) = tmpcrs(i) * f12

     elseif (diroption == 3) then
     ! DO CHI2_XZX
     ! SPS - Polarization

     tmpcrs(i) = (sin(theta)**2) * cos(theta)  * nhpol(i,3,3,3) &
             + cos(theta) * ( nhpol(i,1,3,1) + nhpol(i,2,3,2) ) &
             - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * ( nhpol(i,2,2,3) + nhpol(i,2,3,2) + nhpol(i,3,2,2) ) &
             - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * ( nhpol(i,1,1,3) + nhpol(i,1,3,1) + nhpol(i,3,1,1) ) &
             + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * ( nhpol(i,1,2,3) + nhpol(i,1,3,2) + nhpol(i,2,1,3) &
                                                                    + nhpol(i,2,3,1) + nhpol(i,3,1,2) + nhpol(i,3,2,1) ) &
             + sin(theta) * sin(psi) * ( nhpol(i,2,2,2) + nhpol(i,1,2,1) - nhpol(i,2,3,3) - nhpol(i,3,3,2) ) &
             + sin(theta) * cos(psi) * ( -nhpol(i,1,1,1) - nhpol(i,2,1,2) + nhpol(i,1,3,3) + nhpol(i,3,3,1) ) &
             + (sin(theta)**3) * sin(psi) * ( -nhpol(i,1,1,2) - nhpol(i,1,2,1) - nhpol(i,2,1,1) &
                                            + nhpol(i,2,3,3) + nhpol(i,3,2,3) + nhpol(i,3,3,2) ) &
             + (sin(theta)**3) * cos(psi) * ( nhpol(i,1,2,2) + nhpol(i,2,1,2) + nhpol(i,2,2,1) &
                                            - nhpol(i,1,3,3) - nhpol(i,3,1,3) - nhpol(i,3,3,1) ) &
             + (sin(theta)**3) * (sin(psi)**3) * ( -nhpol(i,2,2,2) + nhpol(i,1,1,2) + nhpol(i,1,2,1) + nhpol(i,2,1,1) ) &
             + (sin(theta)**3) * (cos(psi)**3) * ( nhpol(i,1,1,1) - nhpol(i,1,2,2) - nhpol(i,2,1,2) - nhpol(i,2,2,1) )
 
     tmpcrs(i) = tmpcrs(i) * f12
   
     elseif (diroption == 4) then
     ! DO  CHI2_XXZ
     ! SSP - polarization

     tmpcrs(i) = (sin(theta)**2) * cos(theta)  * nhpol(i,3,3,3) &
             + cos(theta) * ( nhpol(i,1,1,3) + nhpol(i,2,2,3)) &
             - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * ( nhpol(i,2,2,3) + nhpol(i,2,3,2) + nhpol(i,3,2,2) ) &
             - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * ( nhpol(i,1,1,3) + nhpol(i,1,3,1) + nhpol(i,3,1,1) ) &
             + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * ( nhpol(i,1,2,3) + nhpol(i,1,3,2) + nhpol(i,2,1,3) &
                                                                    + nhpol(i,2,3,1) + nhpol(i,3,1,2) + nhpol(i,3,2,1) ) &
             + sin(theta) * sin(psi) * ( nhpol(i,2,2,2) + nhpol(i,1,1,2) - nhpol(i,2,3,3) - nhpol(i,3,2,3) ) &
             + sin(theta) * cos(psi) * ( -nhpol(i,1,1,1) - nhpol(i,2,2,1) + nhpol(i,1,3,3) + nhpol(i,3,1,3) ) &
             + (sin(theta)**3) * sin(psi) * ( -nhpol(i,1,1,2) - nhpol(i,1,2,1) - nhpol(i,2,1,1) &
                                            + nhpol(i,2,3,3) + nhpol(i,3,2,3) + nhpol(i,3,3,2)) &
             + (sin(theta)**3) * cos(psi) * ( nhpol(i,1,2,2) + nhpol(i,2,1,2) + nhpol(i,2,2,1) &
                                            - nhpol(i,1,3,3) - nhpol(i,3,1,3) - nhpol(i,3,3,1)) &
             + (sin(theta)**3) * (sin(psi)**3) * ( -nhpol(i,2,2,2) + nhpol(i,1,1,2) + nhpol(i,1,2,1) + nhpol(i,2,1,1) ) &
             + (sin(theta)**3) * (cos(psi)**3) * ( nhpol(i,1,1,1) - nhpol(i,1,2,2)- nhpol(i,2,1,2) - nhpol(i,2,2,1) ) 

     tmpcrs(i) = tmpcrs(i) * f12 

     elseif (diroption == 5 ) then
     ! DO CHI2_XYZ
     ! SSP Polarization

     tmpcrs(i) = (cos(theta)**2) * ( nhpol(i,1,2,3) - nhpol(i,2,1,3) ) &
             + (sin(theta)**2) * (sin(psi)**2) * ( nhpol(i,3,1,2) - nhpol(i,1,3,2) ) &
             + (sin(theta)**2) * (cos(psi)**2) * ( nhpol(i,2,3,1) - nhpol(i,3,2,1) ) &
             + (sin(theta)**2) * sin(psi) * cos(psi) * ( nhpol(i,1,3,1) - nhpol(i,2,3,2) - nhpol(i,3,1,1) + nhpol(i,3,2,2) ) &
             + sin(theta) * cos(theta) * sin(psi) * ( nhpol(i,1,2,2) - nhpol(i,1,3,3) - nhpol(i,2,1,2) + nhpol(i,3,1,3) ) &
             + sin(theta) * cos(theta) * cos(psi) * ( -nhpol(i,1,2,1) + nhpol(i,2,1,1) - nhpol(i,2,3,3) + nhpol(i,3,2,3) )

     tmpcrs(i) = tmpcrs(i) * f12 

     elseif (diroption == 6) then
     ! DO CHI2_XZY
     ! SPS Polarization

     tmpcrs(i) = (cos(theta)**2) * ( nhpol(i,1,3,2) - nhpol(i,2,3,1) ) &
             + (sin(theta)**2) * (sin(psi)**2) * ( nhpol(i,3,2,1) - nhpol(i,1,2,3) ) &
             + (sin(theta)**2) * (cos(psi)**2) * ( nhpol(i,2,1,3) - nhpol(i,3,1,2) ) &
             + (sin(theta)**2) * sin(psi) * cos(psi) * ( nhpol(i,1,1,3) - nhpol(i,2,2,3) - nhpol(i,3,1,1) + nhpol(i,3,2,2) ) &
             + sin(theta) * cos(theta) * sin(psi) * ( nhpol(i,1,2,2) - nhpol(i,1,3,3) - nhpol(i,2,2,1) + nhpol(i,3,3,1) ) & 
             + sin(theta) * cos(theta) * cos(psi) * ( -nhpol(i,1,1,2) + nhpol(i,2,1,1) - nhpol(i,2,3,3) + nhpol(i,3,3,2) ) 

     tmpcrs(i) = tmpcrs(i) * f12

     elseif (diroption == 7) then
     ! DO CHI2_ZXY
     ! PSS Polarization

     tmpcrs(i) = (cos(theta)**2) * ( nhpol(i,3,1,2) - nhpol(i,3,2,1) ) &
             + (sin(theta)**2) * (sin(psi)**2) * ( nhpol(i,2,3,1) - nhpol(i,2,1,3) ) &
             + (sin(theta)**2) * (cos(psi)**2) * ( nhpol(i,1,2,3) - nhpol(i,1,3,2) ) &
             + (sin(theta)**2) * sin(psi) * cos(psi) * ( nhpol(i,1,1,3) - nhpol(i,1,3,1) - nhpol(i,2,2,3) + nhpol(i,2,3,2) ) &
             + sin(theta) * cos(theta) * sin(psi) * ( nhpol(i,2,1,2) - nhpol(i,2,2,1) - nhpol(i,3,1,3) + nhpol(i,3,3,1) ) &
             + sin(theta) * cos(theta) * cos(psi) * ( -nhpol(i,1,1,2) + nhpol(i,1,2,1) - nhpol(i,3,2,3) + nhpol(i,3,3,2) )

     tmpcrs(i) = tmpcrs(i) * f12

     else
       stop 
     end if

   end do
   
   rcrs_fsfg = tmpcrs*conjg(tmpcrs)

end subroutine SFGFullAvg
