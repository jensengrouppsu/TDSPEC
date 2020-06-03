subroutine Timing(labs, lraman, lasraman, ltpa, lhyperraman, &
             lashyperraman, lsfg, lfluor, lcd, lrvroa, lsecondhyperraman, &
             time1, time2)

  implicit none

  real :: time1, time2

  logical :: labs, lraman, lasraman, ltpa, lhyperraman
  logical :: lashyperraman, lsfg, lfluor, lcd, lrvroa
  logical :: lsecondhyperraman

  character(80) :: calctype

! =============================================================================
! purpose: Output TDSPEC timing information.
!
! input  : labs - OPA calculation
!          lraman - RRS calculation
!          lasraman - AS-RRS calculation
!          ltpa - TPA calculation
!          lhyperraman - RHRS calculation
!          lashyperraman - AS-RHRS calculation
!          lsfg - DR-SFG calculation
!          lfluor - Fluorescence calculation
!          lcd - CD calculation
!          lrvroa - RVROA calculation
!          lsecondhyperraman - RSHRS calculation
!          time1 - start time
!          time2 - end time
!
! =============================================================================

  ! Determine the output string
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
  calctype = adjustl(calctype)
  calctype = '# TDSPEC ' // trim(calctype) // ' Spectrum'

  ! Output information
  write(*,*) trim(calctype)
  write(*,*) "# CPU Time in seconds:", time2 - time1

end subroutine Timing
