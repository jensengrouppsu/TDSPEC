subroutine StickSpectra(lintensity, lraman, lhyperraman, lsfg, &
                lsecondhyperraman, &
                latermoff, lherztell, lfullherztell, lrhrsb2, &
                lnormalize, lfreq, width, freq_foc, rcrs_foc, &
                nmodes, excnum, numcb, rmax, theta, psi, orientation)

  use constants

  implicit none

  character(80) :: charfreq, charpol, charpsi, chartheta, charfcht
  character(80) :: fname

  integer :: nmodes, excnum, numcb, i
  integer :: intfreq, inttheta, intpsi, orientation

  logical :: lintensity, lraman, lhyperraman, lsfg, lnormalize
  logical :: latermoff, lherztell, lfullherztell, lrhrsb2
  logical :: lsecondhyperraman

  real(kindr), intent(in) :: width
  real(kindr), intent(in) :: lfreq
  real(kindr), intent(in) :: freq_foc(excnum*nmodes+numcb)
  real(kindr), intent(in) :: rcrs_foc(excnum*nmodes+numcb)
  real(kindr), intent(in) :: rmax
  real(kindr), intent(in) :: theta
  real(kindr), intent(in) :: psi

  ! The maximum of a Lorentzian is defined as:
  ! 2/(pi*width)
  ! 
  ! We normalize the stick spectrum such that there is correspondance
  ! with the intensity of the Lorentzian broadened spectrum.

  ! Determine if the calculation includes the FC, HT, or both terms.

  if ((lherztell.and.lraman.and.latermoff.eqv..false.).or. &
      (lfullherztell.and.lhyperraman.and.latermoff.eqv..false..and. &
       lrhrsb2.eqv..false.).or. &
      (lherztell.and.lsfg.and.latermoff.eqv..false.)) then
    charfcht = 'AB'
  else if (lherztell.and.lhyperraman.and.latermoff.eqv..false.) then
    charfcht = 'AB1'
  else if (lherztell.and.lhyperraman.and.latermoff) then
    charfcht = 'B1'
  else if (lfullherztell.and.lrhrsb2.and.lhyperraman.and. &
           latermoff.eqv..false.) then
    charfcht = 'AB2'
  else if (lrhrsb2.and.lhyperraman.and.latermoff) then
    charfcht = 'B2'
  else if ((lherztell.and.lraman.and.latermoff).or. &
           (lfullherztell.and.lhyperraman.and.latermoff).or. &
           (lherztell.and.lsfg.and.latermoff)) then
    charfcht = 'B'
  else
    charfcht = 'A'
  end if

  charfcht = adjustl(charfcht)

  ! Stick spectra for RRS

  if (lintensity.and.lraman) then
    ! Convert the laser frequency to a character.
    intfreq = int(1.0E+7_kindr/lfreq)
    write(charfreq,*) intfreq
    charfreq = adjustl(charfreq)
    ! Write the filename (concatenate strings).
    fname = 'RRS_' // trim(charfreq) // 'exc_' // trim(charfcht) // &
            'term_stick.out'
    open(90,file=fname,form='formatted')
    do i=1,excnum*nmodes+numcb
      if (lnormalize.eqv..true.) then
        write(90,*) freq_foc(i),rcrs_foc(i)*two/(pi*width*rmax)
      else
        write(90,*) freq_foc(i),rcrs_foc(i)*two/(pi*width)
      end if
    end do
    close(90)
  end if

  ! Stick spectra for RHRS

  if (lintensity.and.lhyperraman) then
    ! Convert the laser frequency to a character.
    intfreq = int(2.0E+7_kindr/lfreq)
    write(charfreq,*) intfreq
    charfreq = adjustl(charfreq)
    ! Write the filename (concatenate strings).
    fname = 'RHRS_' // trim(charfreq) // 'exc_' // trim(charfcht) // &
            'term_stick.out'
    open(91,file=fname,form='formatted')
    do i=1,excnum*nmodes+numcb
      if (lnormalize.eqv..true.) then
        write(91,*) freq_foc(i),rcrs_foc(i)*two/(pi*width*rmax)
      else
        write(91,*) freq_foc(i),rcrs_foc(i)*two/(pi*width)
      end if
    end do
    close(91)
  end if

  ! Stick spectra for RSHRS

  if (lintensity.and.lsecondhyperraman) then
    ! Convert the laser frequency to a character.
    intfreq = int(3.0E+7_kindr/lfreq)
    write(charfreq,*) intfreq
    charfreq = adjustl(charfreq)
    ! Write the filename (concatenate strings).
    fname = 'RSHRS_' // trim(charfreq) // 'exc_' // trim(charfcht) // &
            'term_stick.out'
    open(93,file=fname,form='formatted')
    do i=1,excnum*nmodes+numcb
      if (lnormalize.eqv..true.) then
        write(93,*) freq_foc(i),rcrs_foc(i)*two/(pi*width*rmax)
      else
        write(93,*) freq_foc(i),rcrs_foc(i)*two/(pi*width)
      end if
    end do
    close(93)
  end if

  ! Stick spectra for DR-SFG

  if (lintensity.and.lsfg) then
    ! Convert the laser frequency and theta/psi to characters.  Also
    ! convert the polarization to a character.
    intfreq = nint(1.0E+7_kindr/lfreq) ! Laser frequency
    write(charfreq,*) intfreq
    charfreq = adjustl(charfreq)
    inttheta = nint(theta*180.0_kindr/pi) ! Angle theta
    write(chartheta,*) inttheta
    chartheta = adjustl(chartheta)
    intpsi = nint(psi*180.0_kindr/pi) ! Angle psi
    write(charpsi,*) intpsi
    charpsi = adjustl(charpsi)
    write(charpol,*) orientation ! Polarization
    charpol = adjustl(charpol)
    ! Write the filename (concatenate strings)
    fname = 'DRSFG_' // trim(charfreq) // 'exc_p' // trim(charpsi) // &
            '_t' // trim(chartheta) // '_pol' // trim(charpol) // &
            '_' // trim(charfcht) // 'term_stick.out'
    open(92,file=fname,form='formatted')
    do i=1,nmodes
      if (lnormalize.eqv..true.) then
        write(92,*) freq_foc(i),rcrs_foc(i)*two/(pi*width*rmax)
      else
        write(92,*) freq_foc(i),rcrs_foc(i)*two/(pi*width)
      end if
    end do
    close(92)
  end if

end subroutine StickSpectra
