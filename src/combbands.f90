Subroutine CombBands(freq, freqcb, excnums, nmodes, numcb)

  use constants

  implicit none

  integer :: numcb, i, j, nmodes
  integer :: excnums(numcb,nmodes)

  real(kindr) :: freq(nmodes), freqcb(numcb)

  character(80) :: fname

! =============================================================================
! purpose: Output the number of combination bands and their contributing normal
!          modes.
!
! input  : freq - normal mode frequencies
!          freqcb - combination band frequencies
!          excnums - excitation numbers of the normal modes in a comb. band
!          nmodes - number of normal modes
!          numcb - number of combination bands
!
! output : Data file containing the combination band information.
!
! =============================================================================

  write(fname,'(a)') 'comb_bands.out'
  open(1000,file=fname,form='formatted')
  write(1000,*) 'Number of combination bands:', numcb
  write(1000,*) 'Order of output: combination band number,'
  write(1000,*) 'normal mode number, excitation number,'
  write(1000,*) 'normal mode frequency'
  do i=1,numcb
    write(1000,*) 'Combination band:', freqcb(i)
    do j=1,nmodes
      if (excnums(i,j) /= 0) then
        write(1000,*) i,j,excnums(i,j),freq(j)
      end if
    end do
  end do
  close(1000)

End Subroutine CombBands
