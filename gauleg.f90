subroutine GAULEG (x1,x2,x,w,n)

!************************************************!
!     Subroutine for Gauss Legendre Quadrature   !
!                                                !
! This program uses Gauss Legendre quadrature to !
! do numerical integration of a function.  What  !
! it does is it uses the roots (y = 0) of the    !
! Legendre polynomials to determine the steps    !
! along the x-axis for integration.  These don't !
! need to be uniformly spaced.  Knowing the      !
! abscissas (x-values), one can then perform the !
! desired integral using appropriate             !
! approximations as sums.                        !
!************************************************!

  Implicit none

!----------------------------------------
! Double precision
!----------------------------------------
  integer, parameter :: kindr = kind(0d0)
!----------------------------------------

  integer, intent(in) :: n
  real(kindr), intent(in) :: x1,x2
  real(kindr), intent(out) :: x(n), w(n)

  real(kindr), parameter :: eps =3.d-14

  integer :: i,j,m
  real(kindr) :: p1, p2, p3, pp, xl, xm, z, z1

!----------------------------------------

! Use the values read in from the tdspec program, where the following
! variables in this subroutine correspond to the variables in tdspec:
!
! n  => tstep
! x1 => tstart
! x2 => tstop
!
! Those values can be changed in the tdspec program to fit the user's 
! needs for the desired accuracy in finding the roots.

  m = (n+1)/2
  xm = 0.5_kindr*(x2+x1)
  xl = 0.5_kindr*(x2-x1)
! Determine the Legendre polynomial (p1) at point z.  
  do i=1,m
     z = cos(3.14592654d0*(i-0.25d0)/(n+0.5d0))
1    continue
     p1 = 1.d0
     p2 = 0.d0
     do j=1,n
        p3=p2
        p2=p1
        p1 = ((2.d0*j-1.d0)*z*p2 - (j-1.d0)*p3)/j
     end do
! Calculate the derivative (pp) of the Legendre polynomial.
     pp = n*(z*p1-p2)/(z*z -1.d0)
     z1 = z
! Use Newton's method for calculating z.
     z = z1 - p1/pp
     if (abs(z-z1).gt.eps) goto 1

! Scale the root (x) to the desired interval of integration, and also the
! symmetric counterpart (integration value at the other end of the data range).

     x(i) = xm-xl*z
     x(n+1-i) = xm+xl*z

! Scale the weighting factors in the same way that the integration values were
! done.

     w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
     w(n+1-i) = w(i)
  end do

end subroutine GAULEG
