subroutine UnitSpherePol(pol, freq, nmodes, lfreq, lherztell, &
                      latermoff, sphereradius)

  use constants

  implicit none

  character*80 :: fname
  character*80 :: frmt
  character*80 :: charfreq
  character*80 :: charlfreq
  character*80 :: charfcht
  character*80 :: chargridx, chargridy, chargridz
  character*80 :: charpolx, charpoly, charpolz
  character*80 :: charmax, chartheta, charphi
  character*100 :: line
  character*40 :: colorline

  integer :: nmodes, i, j, k
  integer :: intlfreq
  integer :: ind(2)

  complex(kindr) :: pol(nmodes,3,3)
  complex(kindr) :: polcont(nmodes,36,18,3)

  real(kindr) :: freq(nmodes)
  real(kindr) :: lfreq
  real(kindr) :: theta, phi
  real(kindr) :: grid(36,18,3)
  real(kindr) :: maxpolcont(nmodes,36,18)
  real(kindr) :: norm(nmodes,36,18)
  real(kindr) :: maximum
  real(kindr) :: color(nmodes,36,18)
  real(kindr) :: sphereradius

  logical :: lherztell, latermoff

! =============================================================================
! purpose: Build a grid of points and draw vector representations of the 
!         polarizability along a unit sphere.  This is presently meant
!         for interfacing with VMD.
!
! input  : nmodes - number of normal modes
!          pol - RRS polarizability tensor (Stokes or anti-Stokes)
!          freq - Vibrational frequency
!          lfreq - Incident frequency
!          lherztell - Inclusion of B terms
!          latermoff - Neglect the A term
!          sphereradius - Radius of the sphere used in the images
!
! =============================================================================

  ! Theta is the polar angle (between the positive x-axis and the vector)
  ! Phi is the azimuthal angle (between the positive z-axis and the vector)
  ! 0 < Theta < 360, 0 < Phi < 180

  ! Build the spherical grid of unit electric fields.
  do i = 1,36 ! Number of theta points (step by 10 degrees)
    theta = real(i)*ten*pi/hndreghty ! radians
    do j = 1,18 ! Number of phi points (step by 10 degrees)
      phi = real(j)*ten*pi/hndreghty ! radians
      grid(i,j,1) = unitradius*dsin(theta)*dcos(phi)
      grid(i,j,2) = unitradius*dsin(theta)*dsin(phi)
      grid(i,j,3) = unitradius*dcos(theta)
    end do
  end do

  ! Contract the polarizability with the electric field.
  ! 11 = xx, 12 = xy, 13 = xz
  ! 21 = yx, 22 = yy, 23 = yz
  ! 31 = zx, 32 = zy, 33 = zz 
  maxpolcont = zero
  do i = 1,nmodes
    maximum = zero
    do j = 1,36
      do k = 1,18
        ! The contraction of the polarizability with the electric
        ! field can be written:
        ! |\alpha_a^{eff}| = \sum_{b=x,y,z}|\alpha_{ab}|E_b
        polcont(i,j,k,1) = &
          (cdsqrt(pol(i,1,1)*conjg(pol(i,1,1)))*grid(j,k,1) + &
           cdsqrt(pol(i,1,2)*conjg(pol(i,1,2)))*grid(j,k,2) + &
           cdsqrt(pol(i,1,3)*conjg(pol(i,1,3)))*grid(j,k,3)) * &
           hartree2cm / pol2si

        polcont(i,j,k,2) = &
          (cdsqrt(pol(i,2,1)*conjg(pol(i,2,1)))*grid(j,k,1) + &
           cdsqrt(pol(i,2,2)*conjg(pol(i,2,2)))*grid(j,k,2) + &
           cdsqrt(pol(i,2,3)*conjg(pol(i,2,3)))*grid(j,k,3)) * &
           hartree2cm / pol2si

        polcont(i,j,k,3) = &
          (cdsqrt(pol(i,3,1)*conjg(pol(i,3,1)))*grid(j,k,1) + &
           cdsqrt(pol(i,3,2)*conjg(pol(i,3,2)))*grid(j,k,2) + &
           cdsqrt(pol(i,3,3)*conjg(pol(i,3,3)))*grid(j,k,3)) * &
           hartree2cm / pol2si

        ! Store the norm of the polarizability vector.
        maximum = maxval(real(polcont(i,j,k,1:3)))

        norm(i,j,k) = real(polcont(i,j,k,1))**2 + &
               real(polcont(i,j,k,2))**2 + &
               real(polcont(i,j,k,3))**2
        norm(i,j,k) = dsqrt(norm(i,j,k))

        if (maximum > maxpolcont(i,j,k)) then
          maxpolcont(i,j,k) = maximum
        end if
      end do
    end do
  end do

  ! Normalize the polarizability vectors to the longest vector.
  do i = 1,nmodes
    maximum = maxval(norm(i,1:36,1:18))
    do j = 1,36
      do k = 1,18
        polcont(i,j,k,1) = polcont(i,j,k,1) / maximum
        polcont(i,j,k,2) = polcont(i,j,k,2) / maximum
        polcont(i,j,k,3) = polcont(i,j,k,3) / maximum
!        polcont(i,j,k,1) = polcont(i,j,k,1) / norm(i,j,k)
!        polcont(i,j,k,2) = polcont(i,j,k,2) / norm(i,j,k)
!        polcont(i,j,k,3) = polcont(i,j,k,3) / norm(i,j,k)
        ! Find the color value.
        color(i,j,k) = norm(i,j,k) / maximum
      end do
    end do
  end do

  ! Convert the grid to the user input radius.
  do i = 1,36 ! Number of theta points (step by 10 degrees)
    theta = real(i)*ten*pi/hndreghty ! radians
    do j = 1,18 ! Number of phi points (step by 10 degrees)
      phi = real(j)*ten*pi/hndreghty ! radians
      grid(i,j,1) = sphereradius*dsin(theta)*dcos(phi)
      grid(i,j,2) = sphereradius*dsin(theta)*dsin(phi)
      grid(i,j,3) = sphereradius*dcos(theta)
    end do
  end do

  ! Determine what terms are included in the calculation
  if (lherztell.and.latermoff.eqv..false.) then
    charfcht = 'AB'
  else if (lherztell.and.latermoff) then
    charfcht = 'B'
  else
    charfcht = 'A'
  end if

  ! Format for the output
  frmt = '(87A)'

  ! Convert the incident frequency to a string  
  intlfreq = nint(1.0E+7_kindr/lfreq)
  write(charlfreq,*) intlfreq
  charlfreq = adjustl(charlfreq)

  ! Build the unit sphere vector files
  do i=1,nmodes
    ! Convert the normal mode frequency to a string
    write(charfreq,'(F7.2)') freq(i)
    charfreq = adjustl(charfreq)
    fname = 'mode' // trim(charfreq) // '_' // trim(charlfreq) // &
            'exc_' // trim(charfcht) // 'term_unitsphere.tcl'
    open(5001, file=fname, form='formatted')
    ! Write the normalization maximum to the file.  Also note the 
    ! angles (in degrees) that give the maximum.
    write(5001,'(66A)') '# Normalization information (angles of the' // &
                        ' maximum and norm of the'
    write(5001,'(32A)') '# largest polarizability vector)'
    ind = maxloc(norm(i,1:36,1:18))
    write(chartheta,'(F5.1)') real(ind(1))*ten
    write(charphi,'(F5.1)') real(ind(2))*ten
    chartheta = adjustl(chartheta)
    charphi = adjustl(charphi)
    write(5001,'(46A)') '# Theta(polar) = ' // trim(chartheta) // &
                        ', Phi(azimuthal) = ' // trim(charphi)
    maximum = maxval(norm(i,1:36,1:18))
    write(charmax,'(ES12.5)') maximum
    write(5001,'(44A)') '# Normalization Maximum (a.u.): ' // trim(charmax)
    ! Revised coloring scheme due to the inflexible/stupid way colors 
    ! are implemented in VMD 1.8.7 (black is color number 16 so it 
    ! isn't replaced).  We try to avoid replacing commonly used
    ! colors for different atom types.
    write(5001,'(31A)') 'color change rgb 11 0.1 0.1 1.0'
    write(5001,'(31A)') 'color change rgb 12 0.2 0.2 1.0'
    write(5001,'(31A)') 'color change rgb 13 0.3 0.3 1.0'
    write(5001,'(31A)') 'color change rgb 14 0.4 0.4 1.0'
    write(5001,'(31A)') 'color change rgb 15 0.5 0.5 1.0'
    write(5001,'(31A)') 'color change rgb 17 0.6 0.6 1.0'
    write(5001,'(31A)') 'color change rgb 18 0.7 0.7 1.0'
    write(5001,'(31A)') 'color change rgb 19 0.8 0.8 1.0'
    write(5001,'(31A)') 'color change rgb 20 0.9 0.9 1.0'
    write(5001,'(31A)') 'color change rgb 21 1.0 0.1 0.1'
    write(5001,'(31A)') 'color change rgb 22 1.0 0.2 0.2'
    write(5001,'(31A)') 'color change rgb 23 1.0 0.3 0.3'
    write(5001,'(31A)') 'color change rgb 24 1.0 0.4 0.4'
    write(5001,'(31A)') 'color change rgb 25 1.0 0.5 0.5'
    write(5001,'(31A)') 'color change rgb 26 1.0 0.6 0.6'
    write(5001,'(31A)') 'color change rgb 27 1.0 0.7 0.7'
    write(5001,'(31A)') 'color change rgb 28 1.0 0.8 0.8'
    write(5001,'(31A)') 'color change rgb 29 1.0 0.9 0.9'
    do j = 1,36
      do k = 1,18
        ! Convert the grid to a string
        write(chargridx,'(F8.5)') grid(j,k,1)
        write(chargridy,'(F8.5)') grid(j,k,2)
        write(chargridz,'(F8.5)') grid(j,k,3)
        ! Convert the contracted polarizability to a string
        write(charpolx,'(ES12.5)') real(polcont(i,j,k,1))
        write(charpoly,'(ES12.5)') real(polcont(i,j,k,2))
        write(charpolz,'(ES12.5)') real(polcont(i,j,k,3))
        ! Determine the color of the vector (obnoxious)
        if (0.00000<=color(i,j,k).and.color(i,j,k)<0.047619) then
          colorline = 'draw color 0' ! Blue (smallest intensity)
        else if (0.047619<=color(i,j,k).and.color(i,j,k)<0.095238) then
          colorline = 'draw color 11'
        else if (0.095238<=color(i,j,k).and.color(i,j,k)<0.142857) then
          colorline = 'draw color 12'
        else if (0.142857<=color(i,j,k).and.color(i,j,k)<0.190476) then
          colorline = 'draw color 13'
        else if (0.190476<=color(i,j,k).and.color(i,j,k)<0.238095) then
          colorline = 'draw color 14'
        else if (0.238095<=color(i,j,k).and.color(i,j,k)<0.285714) then
          colorline = 'draw color 15'
        else if (0.285714<=color(i,j,k).and.color(i,j,k)<0.333333) then
          colorline = 'draw color 17'
        else if (0.333333<=color(i,j,k).and.color(i,j,k)<0.380952) then
          colorline = 'draw color 18'
        else if (0.380952<=color(i,j,k).and.color(i,j,k)<0.428571) then
          colorline = 'draw color 19'
        else if (0.428571<=color(i,j,k).and.color(i,j,k)<0.476190) then
          colorline = 'draw color 20'
        else if (0.476190<=color(i,j,k).and.color(i,j,k)<0.523809) then
          colorline = 'draw color 8' ! White
        else if (0.523809<=color(i,j,k).and.color(i,j,k)<0.571428) then
          colorline = 'draw color 29'
        else if (0.571428<=color(i,j,k).and.color(i,j,k)<0.619047) then
          colorline = 'draw color 28'
        else if (0.619047<=color(i,j,k).and.color(i,j,k)<0.666666) then
          colorline = 'draw color 27'
        else if (0.666666<=color(i,j,k).and.color(i,j,k)<0.714285) then
          colorline = 'draw color 26'
        else if (0.714285<=color(i,j,k).and.color(i,j,k)<0.761904) then
          colorline = 'draw color 25'
        else if (0.761904<=color(i,j,k).and.color(i,j,k)<0.809523) then
          colorline = 'draw color 24'
        else if (0.809523<=color(i,j,k).and.color(i,j,k)<0.857142) then
          colorline = 'draw color 23'
        else if (0.857142<=color(i,j,k).and.color(i,j,k)<0.904761) then
          colorline = 'draw color 22'
        else if (0.904761<=color(i,j,k).and.color(i,j,k)<0.952380) then
          colorline = 'draw color 21'
        else if (0.952380<=color(i,j,k).and.color(i,j,k)<=1.00000) then
          colorline = 'draw color 1' ! Red (largest intensity)
        end if
        line = 'vmd_draw_vector 0 {' // trim(chargridx) // ' ' // &
               trim(chargridy) // ' ' // trim(chargridz) // '} {' // &
               trim(charpolx) // ' ' // trim(charpoly) // ' ' // &
               trim(charpolz) // '} 1.0 30 0.08'
        write(5001,'(13A)') colorline
        write(5001,frmt) line
      end do
    end do
    close(5001)
  end do

end subroutine UnitSpherePol

subroutine UnitSphereHpol(hpol, freq, nmodes, lfreq, lherztell, &
                      lfullherztell, latermoff, lrhrsb2, sphereradius)

  use constants
 
  implicit none

  character*80 :: fname
  character*80 :: frmt
  character*80 :: charfreq
  character*80 :: charlfreq
  character*80 :: charfcht
  character*80 :: chargridx, chargridy, chargridz
  character*80 :: charhpolx, charhpoly, charhpolz
  character*80 :: charmax, chartheta, charphi
  character*100 :: line
  character*40 :: colorline

  integer :: nmodes, i, j, k
  integer :: intlfreq
  integer :: ind(2)

  complex(kindr) :: hpol(nmodes,3,3,3)
  complex(kindr) :: hpolcont(nmodes,36,18,3)

  real(kindr) :: freq(nmodes)
  real(kindr) :: lfreq
  real(kindr) :: theta, phi
  real(kindr) :: grid(36,18,3)
  real(kindr) :: maxhpolcont(nmodes,36,18)
  real(kindr) :: norm(nmodes,36,18)
  real(kindr) :: maximum
  real(kindr) :: color(nmodes,36,18)
  real(kindr) :: sphereradius
 
  logical :: lherztell, lfullherztell, latermoff, lrhrsb2

! =============================================================================
! purpose: Build a grid of points and draw vector representations of the 
!         hyperpolarizability along a unit sphere.  This is presently meant
!         for interfacing with VMD.
!
! input  : nmodes - number of normal modes
!          hpol - RHRS hyperpolarizability tensor (Stokes or anti-Stokes)
!          freq - Vibrational frequency
!          lfreq - Incident frequency
!          lherztell - Inclusion of B1 terms only
!          lfullherztell - Inclusion of full B terms
!          latermoff - Neglect the A term
!          lrhrsb2 - Inclusion of B2 terms only
!          sphereradius - Radius of the sphere used in the images
!
! =============================================================================

  ! Theta is the polar angle (between the positive x-axis and the vector)
  ! Phi is the azimuthal angle (between the positive z-axis and the vector)
  ! 0 < Theta < 360, 0 < Phi < 180

  ! Build the spherical grid of unit electric fields.
  do i = 1,36 ! Number of theta points (step by 10 degrees)
    theta = real(i)*ten*pi/hndreghty ! radians
    do j = 1,18 ! Number of phi points (step by 10 degrees)
      phi = real(j)*ten*pi/hndreghty ! radians
      grid(i,j,1) = unitradius*dsin(theta)*dcos(phi)
      grid(i,j,2) = unitradius*dsin(theta)*dsin(phi)
      grid(i,j,3) = unitradius*dcos(theta)
    end do
  end do 

  ! Contract the hyperpolarizability with the electric field.
  ! 111 = xxx, 112 = xxy, 113 = xxz
  ! 121 = xyx, 122 = xyy, 123 = xyz
  ! 131 = xzx, 132 = xzy, 133 = xzz 
  !
  ! 211 = yxx, 212 = yxy, 213 = yxz
  ! 221 = yyx, 222 = yyy, 223 = yyz
  ! 231 = yzx, 232 = yzy, 233 = yzz
  !
  ! 311 = zxx, 312 = zxy, 313 = zxz
  ! 321 = zyx, 322 = zyy, 323 = zyz
  ! 331 = zzx, 332 = zzy, 333 = zzz
  maxhpolcont = zero
  do i = 1,nmodes
    maximum = zero
    do j = 1,36
      do k = 1,18
        ! The contraction of the hyperpolarizability with the electric
        ! field can be written:
        ! \beta_a^{eff} = \sum_{b=x,y,z}\sum_{c=x,y,z}\beta_{abc}E_cE_b
        hpolcont(i,j,k,1) = &
          (cdsqrt(hpol(i,1,1,1)*conjg(hpol(i,1,1,1)))*grid(j,k,1)*grid(j,k,1) + &
           cdsqrt(hpol(i,1,1,2)*conjg(hpol(i,1,1,2)))*grid(j,k,1)*grid(j,k,2) + & 
           cdsqrt(hpol(i,1,1,3)*conjg(hpol(i,1,1,3)))*grid(j,k,1)*grid(j,k,3) + & 
           cdsqrt(hpol(i,1,2,1)*conjg(hpol(i,1,2,1)))*grid(j,k,2)*grid(j,k,1) + & 
           cdsqrt(hpol(i,1,2,2)*conjg(hpol(i,1,2,2)))*grid(j,k,2)*grid(j,k,2) + & 
           cdsqrt(hpol(i,1,2,3)*conjg(hpol(i,1,2,3)))*grid(j,k,2)*grid(j,k,3) + & 
           cdsqrt(hpol(i,1,3,1)*conjg(hpol(i,1,3,1)))*grid(j,k,3)*grid(j,k,1) + & 
           cdsqrt(hpol(i,1,3,2)*conjg(hpol(i,1,3,2)))*grid(j,k,3)*grid(j,k,2) + & 
           cdsqrt(hpol(i,1,3,3)*conjg(hpol(i,1,3,3)))*grid(j,k,3)*grid(j,k,3)) / &
           au2si

        hpolcont(i,j,k,2) = &
          (cdsqrt(hpol(i,2,1,1)*conjg(hpol(i,2,1,1)))*grid(j,k,1)*grid(j,k,1) + &
           cdsqrt(hpol(i,2,1,2)*conjg(hpol(i,2,1,2)))*grid(j,k,1)*grid(j,k,2) + & 
           cdsqrt(hpol(i,2,1,3)*conjg(hpol(i,2,1,3)))*grid(j,k,1)*grid(j,k,3) + & 
           cdsqrt(hpol(i,2,2,1)*conjg(hpol(i,2,2,1)))*grid(j,k,2)*grid(j,k,1) + & 
           cdsqrt(hpol(i,2,2,2)*conjg(hpol(i,2,2,2)))*grid(j,k,2)*grid(j,k,2) + & 
           cdsqrt(hpol(i,2,2,3)*conjg(hpol(i,2,2,3)))*grid(j,k,2)*grid(j,k,3) + & 
           cdsqrt(hpol(i,2,3,1)*conjg(hpol(i,2,3,1)))*grid(j,k,3)*grid(j,k,1) + & 
           cdsqrt(hpol(i,2,3,2)*conjg(hpol(i,2,3,2)))*grid(j,k,3)*grid(j,k,2) + & 
           cdsqrt(hpol(i,2,3,3)*conjg(hpol(i,2,3,3)))*grid(j,k,3)*grid(j,k,3)) / &
           au2si

        hpolcont(i,j,k,3) = &
          (cdsqrt(hpol(i,3,1,1)*conjg(hpol(i,3,1,1)))*grid(j,k,1)*grid(j,k,1) + &
           cdsqrt(hpol(i,3,1,2)*conjg(hpol(i,3,1,2)))*grid(j,k,1)*grid(j,k,2) + & 
           cdsqrt(hpol(i,3,1,3)*conjg(hpol(i,3,1,3)))*grid(j,k,1)*grid(j,k,3) + & 
           cdsqrt(hpol(i,3,2,1)*conjg(hpol(i,3,2,1)))*grid(j,k,2)*grid(j,k,1) + & 
           cdsqrt(hpol(i,3,2,2)*conjg(hpol(i,3,2,2)))*grid(j,k,2)*grid(j,k,2) + & 
           cdsqrt(hpol(i,3,2,3)*conjg(hpol(i,3,2,3)))*grid(j,k,2)*grid(j,k,3) + & 
           cdsqrt(hpol(i,3,3,1)*conjg(hpol(i,3,3,1)))*grid(j,k,3)*grid(j,k,1) + & 
           cdsqrt(hpol(i,3,3,2)*conjg(hpol(i,3,3,2)))*grid(j,k,3)*grid(j,k,2) + & 
           cdsqrt(hpol(i,3,3,3)*conjg(hpol(i,3,3,3)))*grid(j,k,3)*grid(j,k,3)) / &
           au2si

        ! Store the norm of the hyperpolarizability vector.
        maximum = maxval(real(hpolcont(i,j,k,1:3)))

        norm(i,j,k) = real(hpolcont(i,j,k,1))**2 + &
               real(hpolcont(i,j,k,2))**2 + &
               real(hpolcont(i,j,k,3))**2
        norm(i,j,k) = dsqrt(norm(i,j,k))

        if (maximum > maxhpolcont(i,j,k)) then
          maxhpolcont(i,j,k) = maximum
        end if
      end do
    end do
  end do

  ! Normalize the hyperpolarizability vectors to the longest vector.
  do i = 1,nmodes
    maximum = maxval(norm(i,1:36,1:18))
    do j = 1,36
      do k = 1,18
        hpolcont(i,j,k,1) = hpolcont(i,j,k,1) / maximum
        hpolcont(i,j,k,2) = hpolcont(i,j,k,2) / maximum
        hpolcont(i,j,k,3) = hpolcont(i,j,k,3) / maximum
!        hpolcont(i,j,k,1) = hpolcont(i,j,k,1) / norm(i,j,k)
!        hpolcont(i,j,k,2) = hpolcont(i,j,k,2) / norm(i,j,k)
!        hpolcont(i,j,k,3) = hpolcont(i,j,k,3) / norm(i,j,k)
        ! Find the color value.
        color(i,j,k) = norm(i,j,k) / maximum
      end do
    end do
  end do

  ! Convert the grid to the user input radius.
  do i = 1,36 ! Number of theta points (step by 10 degrees)
    theta = real(i)*ten*pi/hndreghty ! radians
    do j = 1,18 ! Number of phi points (step by 10 degrees)
      phi = real(j)*ten*pi/hndreghty ! radians
      grid(i,j,1) = sphereradius*dsin(theta)*dcos(phi)
      grid(i,j,2) = sphereradius*dsin(theta)*dsin(phi)
      grid(i,j,3) = sphereradius*dcos(theta)
    end do
  end do

  ! Determine what terms are included in the calculation
  if (lfullherztell.and.(latermoff.eqv..false.).and. &
      (lrhrsb2.eqv..false.)) then
    charfcht = 'AB'
  else if (lherztell.and.latermoff.eqv..false.) then
    charfcht = 'AB1'
  else if (lherztell.and.latermoff) then
    charfcht = 'B1'
  else if (lfullherztell.and.lrhrsb2.and.latermoff.eqv..false.) then
    charfcht = 'AB2'
  else if (lrhrsb2.and.latermoff) then
    charfcht = 'B2'
  else if (lfullherztell.and.latermoff) then
    charfcht = 'B'
  else
    charfcht = 'A'
  end if

  ! Format for the output
  frmt = '(87A)'

  ! Convert the incident frequency to a string  
  intlfreq = nint(2.0E+7_kindr/lfreq)
  write(charlfreq,*) intlfreq
  charlfreq = adjustl(charlfreq)

  ! Build the unit sphere vector files
  do i=1,nmodes
    ! Convert the normal mode frequency to a string
    write(charfreq,'(F7.2)') freq(i)
    charfreq = adjustl(charfreq)
    fname = 'mode' // trim(charfreq) // '_' // trim(charlfreq) // &
            'exc_' // trim(charfcht) // 'term_unitsphere.tcl'
    open(5002, file=fname, form='formatted')
    ! Write the normalization maximum to the file.  Also note the 
    ! angles (in degrees) that give the maximum.
    write(5002,'(66A)') '# Normalization information (angles of the' // &
                        ' maximum and norm of the'
    write(5002,'(37A)') '# largest hyperpolarizability vector)'
    ind = maxloc(norm(i,1:36,1:18))
    write(chartheta,'(F5.1)') real(ind(1))*ten
    write(charphi,'(F5.1)') real(ind(2))*ten
    chartheta = adjustl(chartheta)
    charphi = adjustl(charphi)
    write(5002,'(46A)') '# Theta(polar) = ' // trim(chartheta) // &
                        ', Phi(azimuthal) = ' // trim(charphi)
    maximum = maxval(norm(i,1:36,1:18))
    write(charmax,'(ES12.5)') maximum
    write(5002,'(44A)') '# Normalization Maximum (a.u.): ' // trim(charmax)
    ! Revised coloring scheme due to the inflexible/stupid way colors 
    ! are implemented in VMD 1.8.7 (black is color number 16 so it 
    ! isn't replaced).  We try to avoid replacing colors for commonly
    ! used atom types.
    write(5002,'(31A)') 'color change rgb 11 0.1 0.1 1.0'
    write(5002,'(31A)') 'color change rgb 12 0.2 0.2 1.0'
    write(5002,'(31A)') 'color change rgb 13 0.3 0.3 1.0'
    write(5002,'(31A)') 'color change rgb 14 0.4 0.4 1.0'
    write(5002,'(31A)') 'color change rgb 15 0.5 0.5 1.0'
    write(5002,'(31A)') 'color change rgb 17 0.6 0.6 1.0'
    write(5002,'(31A)') 'color change rgb 18 0.7 0.7 1.0'
    write(5002,'(31A)') 'color change rgb 19 0.8 0.8 1.0'
    write(5002,'(31A)') 'color change rgb 20 0.9 0.9 1.0'
    write(5002,'(31A)') 'color change rgb 21 1.0 0.1 0.1'
    write(5002,'(31A)') 'color change rgb 22 1.0 0.2 0.2'
    write(5002,'(31A)') 'color change rgb 23 1.0 0.3 0.3'
    write(5002,'(31A)') 'color change rgb 24 1.0 0.4 0.4'
    write(5002,'(31A)') 'color change rgb 25 1.0 0.5 0.5'
    write(5002,'(31A)') 'color change rgb 26 1.0 0.6 0.6'
    write(5002,'(31A)') 'color change rgb 27 1.0 0.7 0.7'
    write(5002,'(31A)') 'color change rgb 28 1.0 0.8 0.8'
    write(5002,'(31A)') 'color change rgb 29 1.0 0.9 0.9'
    do j = 1,36
      do k = 1,18
        ! Convert the grid to a string
        write(chargridx,'(F8.5)') grid(j,k,1)
        write(chargridy,'(F8.5)') grid(j,k,2)
        write(chargridz,'(F8.5)') grid(j,k,3)
        ! Convert the contracted hyperpolarizability to a string
        write(charhpolx,'(ES12.5)') real(hpolcont(i,j,k,1))
        write(charhpoly,'(ES12.5)') real(hpolcont(i,j,k,2))
        write(charhpolz,'(ES12.5)') real(hpolcont(i,j,k,3))
        ! Determine the color of the vector (obnoxious)
        if (0.00000<=color(i,j,k).and.color(i,j,k)<0.047619) then
          colorline = 'draw color 0' ! Blue (smallest intensity)
        else if (0.047619<=color(i,j,k).and.color(i,j,k)<0.095238) then
          colorline = 'draw color 11'
        else if (0.095238<=color(i,j,k).and.color(i,j,k)<0.142857) then
          colorline = 'draw color 12'
        else if (0.142857<=color(i,j,k).and.color(i,j,k)<0.190476) then
          colorline = 'draw color 13'
        else if (0.190476<=color(i,j,k).and.color(i,j,k)<0.238095) then
          colorline = 'draw color 14'
        else if (0.238095<=color(i,j,k).and.color(i,j,k)<0.285714) then
          colorline = 'draw color 15'
        else if (0.285714<=color(i,j,k).and.color(i,j,k)<0.333333) then
          colorline = 'draw color 17'
        else if (0.333333<=color(i,j,k).and.color(i,j,k)<0.380952) then
          colorline = 'draw color 18'
        else if (0.380952<=color(i,j,k).and.color(i,j,k)<0.428571) then
          colorline = 'draw color 19'
        else if (0.428571<=color(i,j,k).and.color(i,j,k)<0.476190) then
          colorline = 'draw color 20'
        else if (0.476190<=color(i,j,k).and.color(i,j,k)<0.523809) then
          colorline = 'draw color 8' ! White
        else if (0.523809<=color(i,j,k).and.color(i,j,k)<0.571428) then
          colorline = 'draw color 29'
        else if (0.571428<=color(i,j,k).and.color(i,j,k)<0.619047) then
          colorline = 'draw color 28'
        else if (0.619047<=color(i,j,k).and.color(i,j,k)<0.666666) then
          colorline = 'draw color 27'
        else if (0.666666<=color(i,j,k).and.color(i,j,k)<0.714285) then
          colorline = 'draw color 26'
        else if (0.714285<=color(i,j,k).and.color(i,j,k)<0.761904) then
          colorline = 'draw color 25'
        else if (0.761904<=color(i,j,k).and.color(i,j,k)<0.809523) then
          colorline = 'draw color 24'
        else if (0.809523<=color(i,j,k).and.color(i,j,k)<0.857142) then
          colorline = 'draw color 23'
        else if (0.857142<=color(i,j,k).and.color(i,j,k)<0.904761) then
          colorline = 'draw color 22'
        else if (0.904761<=color(i,j,k).and.color(i,j,k)<0.952380) then
          colorline = 'draw color 21'
        else if (0.952380<=color(i,j,k).and.color(i,j,k)<=1.00000) then
          colorline = 'draw color 1' ! Red (largest intensity)
        end if
        line = 'vmd_draw_vector 0 {' // trim(chargridx) // ' ' // &
               trim(chargridy) // ' ' // trim(chargridz) // '} {' // &
               trim(charhpolx) // ' ' // trim(charhpoly) // ' ' // &
               trim(charhpolz) // '} 1.0 30 0.08'
        write(5002,'(13A)') colorline
        write(5002,frmt) line
      end do
    end do
    close(5002)
  end do

end subroutine UnitSphereHpol
