subroutine OutputHeader()

  implicit none

  character(3) :: version
  character(11) :: vline

! =============================================================================
! purpose: Output the TDSPEC header.
! =============================================================================

  version = '2.0'
  vline = 'TDSPEC v' // version

  write(*,*) repeat('#', 70)
  write(*,'(X,A,68X,A)') '#', '#'
  write(*,'(X,2(A,28X),A)') '#', vline, ' #'
  write(*,'(X,2(A,12X),A)') '#', &
    'Time-Dependent Spectroscopic Simulations of', ' #'
  write(*,'(X,2(A,15X),A)') '#', &
    'Linear and Nonlinear Optical Processes', '#'
  write(*,'(X,A,68X,A)') '#', '#'
  write(*,'(X,2(A,19X),A)') '#', 'Authors: Daniel W. Silverstein', '#'
  write(*,'(X,2(A,25X),A)') '#   ', 'Philip A. Weiss', '#'
  write(*,'(X,A,25X,A,23X,A)') '#   ', 'Dhabih V. Chulhai', '#'
  write(*,'(X,A,25X,A,29X,A)') '#   ', 'Zhongwei Hu', '#'
  write(*,'(X,2(A,28X),A)') '#', 'Lasse Jensen', '#'
  write(*,'(X,A,68X,A)') '#', '#'
  write(*,*) repeat('#', 70)

end subroutine OutputHeader

subroutine OutputCalcInfo(labs, lraman, lasraman, ltpa, lhyperraman, &
                 lashyperraman, lsfg, lfluor, lcd, lrvroa, lsecondhyperraman, &
                 lfreq, sshift, solvmod, excnum, fc_screen, &
                 freq_screen, theta, psi, orientation, lintensity)

  use constants
 
  implicit none

  logical :: labs, lraman, lasraman, ltpa, lhyperraman
  logical :: lashyperraman, lsfg, lfluor, lcd, lrvroa, lintensity
  logical :: lsecondhyperraman

  real(kindr), intent(in) :: lfreq
  real(kindr), intent(in) :: sshift
  real(kindr), intent(in) :: fc_screen
  real(kindr), intent(in) :: freq_screen
  real(kindr), intent(in) :: theta
  real(kindr), intent(in) :: psi

  integer, intent(in) :: solvmod(4)
  integer, intent(in) :: excnum
  integer, intent(in) :: orientation

  character(15) :: calctype
  character(70) :: line

! =============================================================================
! purpose: Output calculation information.
!
! input: labs - OPA calculation
!        lraman - RRS calculation
!        lasraman - AS-RRS calculation
!        ltpa - TPA calculation
!        lhyperraman - RHRS calculation
!        lashyperraman - AS-RHRS calculation
!        lsfg - DR-SFG calculation
!        lfluor - Fluorescence calculation
!        lcd - CD calculation
!        lrvroa - RVROA calculation
!        lsecondhyperraman - RSHRS calculation
!        lfreq - laser frequency
!        sshift - solvent shift
!        solvmod - type of solvent model
!        excnum - excitation number
!        fc_screen - FC integral screening for combination bands
!        freq_screen - frequency screening for combination bands
!        theta - Euler angle theta for DR-SFG
!        psi - Euler angle psi for DR-SFG
!        orientation - polarization of the light for averaging in DR-SFG
!        lintensity - true if a stick spectrum is needed
! 
! =============================================================================

  ! Output information about the calculation type
  if (labs) then
    calctype = 'OPA'
  else if (lraman) then
    calctype = 'RRS'
  else if (lasraman) then
    calctype = 'AS-RRS'
  else if (ltpa) then
    calctype = 'TPA'
  else if (lhyperraman) then
    calctype = 'RHRS'
  else if (lashyperraman) then
    calctype = 'AS-RHRS'
  else if (lsfg) then
    calctype = 'DR-SFG'
  else if (lfluor) then
    calctype = 'Fluorescence'
  else if (lcd) then
    calctype = 'CD'
  else if (lrvroa) then
    calctype = 'RVROA'
  else if (lsecondhyperraman) then
    calctype = 'RSHRS'
  end if
  line = '# TDSPEC ' // trim(calctype) // ' Calculation'
  write(*,*) line

  ! Output the laser frequency and solvent shift
  write(*,*) '# Incident Frequency', lfreq
  write(*,*) '# Solvent shift', sshift
  
  ! Output the solvent model being used
  if (solvmod(1)==1) then
    write(*,*) &
      '# Lifetime Model = Simple Homogeneous and Inhomogenous Broadening'
  else if (solvmod(2)==1) then
    write(*,*) &
      '# Lifetime Model = High Temp. Limit Overdamped Brownian Oscillator'
  else if (solvmod(3)==1) then
    write(*,*) &
      '# Lifetime Model = General Temp. Limit Overdamped Brownian Oscillator'
  else if (solvmod(4)==1) then
    write(*,*) &
      '# Lifetime Model = High Temp. Limit Overdamped Brownian Oscillator'
    write(*,*) &
      '#                  with excitation energy distribution'
  end if

  ! Print information on screening the combination bands when the
  ! excitation number is greater than 1.
  if (excnum > 1) then
    write(*,'(A22,E9.1)') '# FC factor screening', fc_screen
    write(*,'(A32,F10.2)') '# Frequency screening threshold', freq_screen
    ! Automatically print out the stick spectrum if studying 
    ! combination bands and overtones.
    lintensity = .true.
  end if

  ! Write the Euler angles of the DR-SFG calculation (in degrees).  Also
  ! output the polarization.
  if (lsfg) then
    write(*,'(A14,F6.2)') ' # Angle theta ', theta * 180.0_kindr / pi
    write(*,'(A12,2X,F6.2)') ' # Angle psi', psi * 180.0_kindr / pi
    if (orientation==1) then
      write(*,'(A23)') '# Polarization = PPP 1'
    else if (orientation==2) then
      write(*,'(A23)') '# Polarization = PSS 2'
    else if (orientation==3) then
      write(*,'(A23)') '# Polarization = SPS 3'
    else if (orientation==4) then
      write(*,'(A23)') '# Polarization = SSP 4'
    else if (orientation==5) then
      write(*,'(A23)') '# Polarization = SSP 5'
    else if (orientation==6) then
      write(*,'(A23)') '# Polarization = SPS 6'
    else if (orientation==7) then
      write(*,'(A23)') '# Polarization = PSS 7'
    end if
  end if

end subroutine OutputCalcInfo

subroutine OutputLifetime(nexci, solvmod, gam, gam2, kappa, freqbath, deltabath)

  use constants

  implicit none

  integer, intent(in) :: nexci
  integer :: iex
  integer, intent(in) :: solvmod(4)

  real(kindr), intent(in) :: gam(nexci)
  real(kindr), intent(in) :: gam2
  real(kindr), intent(in) :: kappa(nexci)
  real(kindr), intent(in) :: freqbath(nexci) 
  real(kindr), intent(in) :: deltabath(nexci) 

  character(80) :: frmt

! =============================================================================
! purpose: Output information for the lifetime of the excited states.
!
! input: nexci - number of excited states
!        solvmod - type of solvent model
!        gam - homogeneous broadening of the excited state
!        gam2 - inhomogeneous broadening for the simple solvent model
!        kappa - lifetime parameter (solvent model II or III) 
!        freqbath - bath frequency for solvent model III
!        deltabath - bath dimensionless displacements for solvent model III
! 
! =============================================================================

  !write(*,*) 'NEXCI ', nexci
  if (solvmod(1)==1) then
    write(*,*) "# Excited State       Gamma (cm-1)      Theta (cm-1)"
    frmt = "(1X,1A,10X,I2,10X,F7.2,10X,F8.2)"
    do iex=1,nexci
      write(*,frmt) "#", iex, gam(iex), gam2
    end do
  else if (solvmod(2)==1) then
    write(*,*) "# Excited State       FWHM (cm-1)       Kappa"
    frmt = "(1X,1A,10X,I2,11X,F7.2,8X,F6.2)"
    do iex=1,nexci
      write(*,frmt) "#", iex, gam(iex), kappa(iex)
    end do
  else if (solvmod(3)==1) then
    write(*,*) "# Excited State       Omega_B (cm-1)     Delta_B      Kappa"
    do iex=1,nexci
      write(*,*) "#", iex, "      ", freqbath(iex), "      ", deltabath(iex), &
                "      ", kappa(iex)
    end do
  end if  

end subroutine OutputLifetime

subroutine PrintPol(freq, lfreq, pol, nmodes, polA, polG, polC, polD, &
                    polAs, polGs, polDs, lprintAtensor, lprintGtensor)

  use constants

  implicit none

  integer, intent(in)                  :: nmodes
  real(kindr), intent(in)              :: lfreq
  real(kindr), intent(in)              :: freq(nmodes)
  complex(kindr)                       :: pol(nmodes,3,3)
  complex(kindr)                       :: polA(nmodes,3,3,3)
  complex(kindr)                       :: polG(nmodes,3,3)
  complex(kindr)                       :: polC(nmodes,6,6)
  complex(kindr)                       :: polD(nmodes,6,3)
  complex(kindr)                       :: polAs(nmodes,3,3,3)
  complex(kindr)                       :: polGs(nmodes,3,3)
  complex(kindr)                       :: polDs(nmodes,3,6)
  logical, intent(in)                  :: lprintAtensor
  logical, intent(in)                  :: lprintGtensor

  integer                                   :: i, j, k, l, iu, istat
  real(kindr), parameter                    :: threehalf = 1.5_kindr, &
                                               half = 0.5_kindr, &
                                               quart = 0.25_kindr
  real(kindr)                               :: conversion, polconv
  complex(kindr)                            :: trace, trace2(6), trace3
  character(len=1), dimension(3), parameter :: dir = (/'x','y','z'/)
  character(len=2), dimension(6), parameter :: dir2 = (/'xx','xy','xz','yy','yz','zz'/)
  character(3)                              :: version
  character(11)                             :: vline

! =========================================================
! purpose: prints the RRS polarizability derivatives
!          in atomic units (similar to the output
!          in ADF).
!
! input: freq - frequency of the normal modes
!        pol - polarizability derivatives
!        nmodes - number of normal modes
!        polA  - the dipole-quadrupole tensor
!        polAs - the quadrupole-dipole tensor
!        polG  - the dipole-magnetic dipole tensor
!        polGs - the magnetic dipole-electric dipole tensor
!        polC  - the quadrupole-quadrupole tensor
!        polD  - the quadrupole-magnetic dipole tensor
!        polDs  - the magnetic dipole-electric quadrupole
! =========================================================

  ! open file for printing, named 'tdspec_pol.out'
  iu = 101
  open(unit=iu, file='_TDSPEC.pol', status='REPLACE', iostat=istat)

  ! Print header
  version = '1.1'
  vline = 'TDSPEC v' // version

  write(iu,*) repeat('#', 70)
  write(iu,'(X,A,68X,A)') '#', '#'
  write(iu,'(X,2(A,28X),A)') '#', vline, ' #'
  write(iu,'(X,2(A,12X),A)') '#', &
    'Time-Dependent Spectroscopic Simulations of', ' #'
  write(iu,'(X,2(A,15X),A)') '#', &
    'Linear and Nonlinear Optical Processes', '#'
  write(iu,'(X,A,68X,A)') '#', '#'
  write(iu,'(X,2(A,19X),A)') '#', 'Authors: Daniel W. Silverstein', '#'
  write(iu,'(X,2(A,25X),A)') '#   ', 'Philip A. Weiss', '#'
  write(iu,'(X,A,25X,A,23X,A)') '#   ', 'Dhabih V. Chulhai', '#'
  write(iu,'(X,A,25X,A,29X,A)') '#   ', 'Zhongwei Hu', '#'
  write(iu,'(X,2(A,28X),A)') '#', 'Lasse Jensen', '#'
  write(iu,'(X,A,68X,A)') '#', '#'
  write(iu,*) repeat('#', 70)


  ! Write standard information
  write(iu,*) repeat('=', 50)
  write(iu,'(A,F14.2,A)') ' Excitation Frequency: ', lfreq, ' cm-1'
  write(iu,*) repeat('=', 50)
  write(iu,*) "Polarizability derivatives (atomic units)"
  write(iu,*) repeat('-', 50)

  ! Cycle through each mode
  do i=1,nmodes
    ! we convert from unitless polarizability derivatives into
    ! Atomic units of the form outputed by ADF
    conversion = dsqrt(( 8_kindr * pi**2 * speed * freq(i) ) &
        / ( planck )) * dsqrt(( 4e-40_kindr * epsilon0**2 ) / ( amu ))
    polconv = conversion / pol2si
    write(iu,"(A14,F7.2)") " Mode (cm-1)   ", freq(i)

    ! ==============
    ! Polarizability
    ! ==============
    do j=1,3
      do k=1,3
        write(iu,"(' alpha_',2A1,'    ',2E18.8)") dir(j), dir(k), &
              polconv*real(pol(i,j,k)), polconv*aimag(pol(i,j,k))
      end do
    end do

    ! ========
    ! A-tensor
    ! ========
    if (lprintAtensor) then
      ! -----------------------
      ! Roman-type:
      ! El Dipole-El Quadrupole
      ! -----------------------
      do j=1,3
        ! Make traceless
        trace = ( polA(i,j,1,1) + polA(i,j,2,2) + polA(i,j,3,3) ) * half
        do k=1,3
          do l=k,3
            polA(i,j,k,l) = polA(i,j,k,l) * threehalf
            if (k==l) polA(i,j,k,k) = polA(i,j,k,k) - trace
            write(iu,"(' A_',A1,',',2A1,'      ',2E18.8)") dir(j), dir(k), dir(l), &
                  conversion*real(polA(i,j,k,l)), conversion*aimag(polA(i,j,k,l))
          end do
        end do
      end do
      ! -----------------------
      ! Script-type:
      ! El Quadrupole-El Dipole
      ! -----------------------
      do j=1,3
        ! Make traceless
        trace = ( polAs(i,j,1,1) + polAs(i,j,2,2) + polAs(i,j,3,3) ) * half
        do k=1,3
          do l=k,3
            polAs(i,j,k,l) = polAs(i,j,k,l) * threehalf
            if (k==l) polAs(i,j,k,k) = polAs(i,j,k,k) - trace
            write(iu,"(' As_',A1,',',2A1,'     ',2E18.8)") dir(j), dir(k), dir(l), &
                  conversion*real(polAs(i,j,k,l)), conversion*aimag(polAs(i,j,k,l))
          end do
        end do
      end do
    end if

    ! ========
    ! G-tensor
    ! ========
    if (lprintGtensor) then
      ! --------------------
      ! Roman-type:
      ! El Dipole-Mag Dipole
      ! --------------------
      do j=1,3
        do k=1,3
          write(iu,"(' G_',2A1,'        ',2E18.8)") dir(j), dir(k), &
                conversion*real(polG(i,j,k)), conversion*aimag(polG(i,j,k))
        end do
      end do
      ! --------------------
      ! Script-type:
      ! Mag Dipole-El Dipole
      ! --------------------
      do j=1,3
        do k=1,3
          write(iu,"(' Gs_',2A1,'       ',2E18.8)") dir(j), dir(k), &
                conversion*real(polGs(i,j,k)), conversion*aimag(polGs(i,j,k))
        end do
      end do
    end if

    ! ========
    ! C-tensor
    ! ========
    if (lprintAtensor) then
      ! Make traceless
      trace2(1) = ( polC(i,1,1) + polC(i,4,1) + polC(i,6,1) ) * half
      trace2(2) = ( polC(i,1,2) + polC(i,4,2) + polC(i,6,2) ) * half
      trace2(3) = ( polC(i,1,3) + polC(i,4,3) + polC(i,6,3) ) * half
      trace2(4) = ( polC(i,1,4) + polC(i,4,4) + polC(i,6,4) ) * half
      trace2(5) = ( polC(i,1,5) + polC(i,4,5) + polC(i,6,5) ) * half
      trace2(6) = ( polC(i,1,6) + polC(i,4,6) + polC(i,6,6) ) * half
      trace3 = polC(i,1,1) + polC(i,1,4) + polC(i,1,6) &
             + polC(i,4,1) + polC(i,4,4) + polC(i,4,6) &
             + polC(i,6,1) + polC(i,6,4) + polC(i,6,6)
      do j=1,6
        trace = ( polC(i,j,1) + polC(i,j,4) + polC(i,j,6) ) * half
        do k=1,6
          polC(i,j,k) = polC(i,j,k) * threehalf * threehalf
          if (k==1 .or. k==4 .or. k==6) polC(i,j,k) = polC(i,j,k) - trace * threehalf
          if (j==1 .or. j==4 .or. j==6) polC(i,j,k) = polC(i,j,k) - trace2(k) * threehalf
          if ((j==1 .or. j==4 .or. j==6) .and. (k==1 .or. k==4 .or. k==6)) &
            polC(i,j,k) = polC(i,j,k) + trace3 * quart
          write(iu,"(' C_',A2,',',A2,'     ',2E18.8)") dir2(j), dir2(k), &
                conversion*real(polC(i,j,k)), conversion*aimag(polC(i,j,k))
        end do
      end do
    end if

    ! ========
    ! D-tensor
    ! ========
    if (lprintAtensor .and. lprintGtensor) then
      ! ------------------------
      ! Roman-type:
      ! El Quadrupole-Mag Dipole
      ! ------------------------
      do k=1,3
        ! Make tranceless
        trace = ( polD(i,1,k) + polD(i,4,k) + polD(i,6,k) ) * half
        do j=1,6
          polD(i,j,k) = polD(i,j,k) * threehalf
          if (j==1 .or. j==4 .or. j==6) polD(i,j,k) = polD(i,j,k) - trace
          write(iu,"(' D_',A2,','A1,'      ',2E18.8)") dir2(j), dir(k), &
                conversion*real(polD(i,j,k)), conversion*aimag(polD(i,j,k))
        end do
      end do
      ! ------------------------
      ! Script-type:
      ! Mag Dipole-El Quadrupole
      ! ------------------------
      do k=1,3
        ! Make traceless
        trace = ( polDs(i,k,1) + polDs(i,k,4) + polDs(i,k,6) ) * half
        do j=1,6
          polDs(i,k,j) = polDs(i,k,j) * threehalf
          if (j==1 .or. j==4 .or. j==6) polDs(i,k,j) = polDs(i,k,j) - trace
          write(iu,"(' Ds_',A1,','A2,'     ',2E18.8)") dir(k), dir2(j), &
                conversion*real(polDs(i,k,j)), conversion*aimag(polDs(i,k,j))
        end do
      end do
    end if
    write(iu,*) repeat('-', 50)
  end do
  write(iu,*) repeat('=', 50)

  ! close file
  close(iu)

end subroutine PrintPol

subroutine PrintHpol(freq, lfreq, hpol, nmodes)

  use constants

  implicit none

  integer, intent(in)                  :: nmodes
  real(kindr), intent(in)              :: lfreq
  real(kindr), intent(in)              :: freq(nmodes)
  complex(kindr)                       :: hpol(nmodes,3,3,3)

  integer                                   :: i, j, k, l, iu, istat
  real(kindr)                               :: conversion, hpolconv
  character(len=1), dimension(3), parameter :: dir = (/'x','y','z'/)
  character(3)                              :: version
  character(11)                             :: vline

! =========================================================
! purpose: prints the RHRS hyperpolarizability derivatives
!          in atomic units (similar to the output
!          in ADF).
!
! input: freq - frequency of the normal modes
!        hpol - hyperpolarizability derivatives
!        nmodes - number of normal modes
! =========================================================

  ! open file for printing, named 'tdspec_hpol.out'
  iu = 101
  open(unit=iu, file='_TDSPEC.hpol', status='REPLACE', iostat=istat)

  ! Print header
  version = '1.1'
  vline = 'TDSPEC v' // version

  write(iu,*) repeat('#', 70)
  write(iu,'(X,A,68X,A)') '#', '#'
  write(iu,'(X,2(A,28X),A)') '#', vline, ' #'
  write(iu,'(X,2(A,12X),A)') '#', &
    'Time-Dependent Spectroscopic Simulations of', ' #'
  write(iu,'(X,2(A,15X),A)') '#', &
    'Linear and Nonlinear Optical Processes', '#'
  write(iu,'(X,A,68X,A)') '#', '#'
  write(iu,'(X,2(A,19X),A)') '#', 'Authors: Daniel W. Silverstein', '#'
  write(iu,'(X,2(A,25X),A)') '#   ', 'Philip A. Weiss', '#'
  write(iu,'(X,A,25X,A,23X,A)') '#   ', 'Dhabih V. Chulhai', '#'
  write(iu,'(X,A,25X,A,29X,A)') '#   ', 'Zhongwei Hu', '#'
  write(iu,'(X,2(A,28X),A)') '#', 'Lasse Jensen', '#'
  write(iu,'(X,A,68X,A)') '#', '#'
  write(iu,*) repeat('#', 70)


  ! Write standard information
  write(iu,*) repeat('=', 50)
  write(iu,'(A,F14.2,A)') ' Excitation Frequency: ', lfreq, ' cm-1'
  write(iu,*) repeat('=', 50)
  write(iu,*) "Hyperpolarizability derivatives (atomic units)"
  write(iu,*) repeat('-', 50)

  ! Cycle through each mode
  do i=1,nmodes
    ! we convert from unitless hyperpolarizability derivatives into
    ! Atomic units of the form outputed by ADF
    conversion = dsqrt(( 8_kindr * pi**2 * speed * freq(i) ) &
        / ( planck )) * dsqrt(( 4e-40_kindr * epsilon0**2 ) / ( amu )) &
        * pi / 2
    hpolconv = conversion / hpol2si
    write(iu,"(A14,F7.2)") " Mode (cm-1)   ", freq(i)

    ! ===================
    ! Hyperpolarizability
    ! ===================
    do j=1,3
      do k=1,3
         do l=1,3
           write(iu,"(' beta_',3A1,'    ',2E18.8)") dir(j), dir(k), dir(l), &
                 hpolconv*real(hpol(i,j,k,l)), hpolconv*aimag(hpol(i,j,k,l))
         end do
      end do
    end do
    write(iu,*) repeat('-', 50)
  end do
  write(iu,*) repeat('=', 50)

  ! close file
  close(iu)

end subroutine PrintHpol
