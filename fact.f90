module fact

  implicit none

  CONTAINS

  Real FUNCTION f(j)
  
    ! This module is used to obtain factorials up to 15, for use in 
    ! calculating FC overlap integrals.  The code is currently tested for 
    ! factorials ranging from 0 to 15, as required by the program fcfactor.
    ! It may cause problems with integer overflows if you try to make
    ! the maximum number too large (i.e. the factorial becomes enormous for
    ! integers between 10 and 20).

    implicit none

    integer, intent(in) :: j
    
    integer i
    real :: Fact

    if (j == 0) then
      Fact = 1
    else
      Fact = 1
      do i=1,j
        Fact = Fact*i
      end do
    end if

    f = Fact

  End FUNCTION f

end module fact
