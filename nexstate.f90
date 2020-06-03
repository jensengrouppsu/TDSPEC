subroutine nextstate (n, nmodes, l, first, last)

   implicit none

   integer, intent(in)    :: nmodes
   integer, intent(inout) :: n(nmodes), l, first, last

! =============================================================================
! purpose: Calculate the vibrational quantum numbers of the next state
!
! input  : nmodes - number of modes
!
! in/out : n     - vibrational quantum numbers
!          l     - level (sum of n)
!          first - index of first non-zero element in n
!          last  - index of last non-zero element in n
! =============================================================================

   ! ----------------------------------------
   ! Procedure NEXTSTATE of Ruhoff and Ratner
   ! ----------------------------------------
   if (n(1) == l) then
      first = nmodes
      last = nmodes
      l = l + 1
      n(1) = 0
      n(nmodes) = l
   else
      if (n(first) == l) then
         n(first) = 0
         first = first - 1
         n(first) = 1
         if (l > 1) then
            last = nmodes
            n(last) = l - 1
         else
            last = first
         end if
      else
         n(last - 1) = n(last - 1) + 1
         if ((n(last) > 1) .and. (last < nmodes)) then
            n(nmodes) = n(last) - 1
            n(last) = 0
            last = nmodes
         else
            n(last) = n(last) - 1
            if (n(last) == 0) last = last - 1
         end if
      end if
   end if

end subroutine nextstate

