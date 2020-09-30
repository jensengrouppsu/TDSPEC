subroutine Error(errortype, runtype, indata)

  implicit none

  integer :: errortype
  
  character(80) :: runtype
  character(80) :: indata

! =============================================================================
! purpose: Output errors from TDSPEC, if they exist. 
!
! input  : errortype - error number 
!          runtype - type of calculation
!          indata - name of the data file
!
! =============================================================================

  if (errortype.eq.1) then
    ! For errors dealing with the input file contents 
    write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/'
    write(*,*) 'There is an error in the input file.'
    write(*,*) 'Please check that the following keys (and corresponding'
    write(*,*) 'values) are present:'
    write(*,*) '1. Runtype'
    write(*,*) '2. Freq'
    write(*,*) '3. Datafile'
    write(*,*) '4. Width (if RRS, RHRS, SFG)'
    write(*,*) '5. End'
    write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/'
  else if (errortype.eq.2) then
    ! For errors dealing with the runtype
    runtype = adjustl(runtype)
    write(*,*) '---------------------------------------------' // &
      '---------------------'
    write(*,*) 'Input File Error!'
    write(*,'(3A)') ' Option given for Runtype, ', trim(runtype), &
      ' is not recognized.'
    write(*,*) 'Try one of the following options (meaning)'
    write(*,*) '1. Abs (for one-photon absorbance)'
    write(*,*) '2. Raman (for resonance Raman scattering)'
    write(*,*) '3. ASRaman (for anti-Stokes resonance Raman ' // &
      'scattering)'
    write(*,*) '4. TwoAbs (for two-photon absorbance)'
    write(*,*) '5. HyperRaman (for resonance hyper-Raman ' // &
      'scattering)'
    write(*,*) '6. ASHyperRaman (for anti-Stokes resonance ' // &
      'hyper-Raman scattering)'
    write(*,*) '7. SFG (for double resonance sum-frequency ' // &
      'generation)'
    write(*,*) '8. Fluor (for one-photon fluorescence)'
    write(*,*) '9. CD (for circular dichroism)'
    write(*,*) '10. RVROA (for resonance vibrational Raman ' // &
      'optical activity)'
    write(*,*) '---------------------------------------------' // &
      '---------------------'
  else if (errortype.eq.3) then
    ! For errors dealing with not finding a data file
    indata = adjustl(indata)
    write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\'
    write(*,'(X,3A)') 'Data file ', trim(indata), ' cannot be found.'
    write(*,'(X,A)')  'Please check that the data file exists in the'
    write(*,'(X,A)')  'current directory.'
    write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\'
  else if (errortype.eq.4) then
    ! For errors in reading the data file
    indata = adjustl(indata)
    write(*,'(A)') '  /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/'
    write(*,'(A)') '  There is an error in the data file.'
    write(*,'(A)') '  Please check that either:'
    write(*,'(A)') '  1. NEXCI is present in the file.'
    write(*,'(3A)') '  2. That ', trim(indata), &
      ' is a TDSPEC data file.'
    write(*,'(A)') '  /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/'
  end if

  ! End the program execution due to errors
  stop

end subroutine error
