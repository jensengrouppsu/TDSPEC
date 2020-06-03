function position (n, nmodes, l, maxl, ncmb, ntot) result(pos)

   implicit none

   integer             :: pos
   integer, intent(in) :: n(nmodes), nmodes, l, maxl, &
                                ncmb(nmodes + maxl, maxl + 1), ntot(maxl + 1)

! =============================================================================
! purpose: Calculate the position of the current state
!
! input  : n      - vibrational quantum numbers
!          nmodes - number of modes
!          l      - level (sum of n)
!          maxl   - maximum number of levels
!          ncmb   - number of combinations of quantum numbers
!          ntot   - total number of combinations per level
!
! output : pos - position of the current state
! =============================================================================

   integer :: i, j

   ! ---------------------------------------
   ! Procedure POSITION of Ruhoff and Ratner
   ! ---------------------------------------
   pos = ntot(l) -1
   i = 1;
   j = l
   do while (j > 0)
      if (n(i) > 0) j = j - n(i)
      if (j > 0) then
         pos = pos - ncmb(nmodes-i+1, j - 1)
      end if
      i = i + 1
   end do

end function position
