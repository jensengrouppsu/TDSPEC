Subroutine RRSPref(rrs_pref, osc, tdipder, delta, freq, nexci, &
                   nmodes, ndim)

  use constants

  implicit none

  integer :: nexci, nmodes, ndim, iex, j
  integer :: a, b

  character(80) :: fname, fname1

  real(kindr), intent(in) :: osc(nexci,ndim)
  real(kindr), intent(out) :: rrs_pref(nexci,nmodes)
  real(kindr), intent(in) :: tdipder(nexci,nmodes,ndim)
  real(kindr) :: betaA, polmn
  real(kindr) :: pref(3,3)
  real(kindr) :: delta(nexci,nmodes), freq(nmodes)

! =============================================================================
! purpose: Orientational averaging of the B term prefactor for RRS
!
! input  : rrs_pref - orientational averages of the B term prefactors
!          osc - transition dipole moment
!          tdipder - derivatives of transition dipole moments along each mode
!          delta - dimensionless displacements
!          freq - vibrational frequencies
!          nexci - number of excited states
!          nmodes - number of normal modes
!          ndim - final dimension of the property arrays (equals 3)
!
! =============================================================================

  !*********************************************
  ! Calculate the RRS prefactors for the B terms
  !*********************************************

  ! This part calculates the average prefactor for a Herzberg-Teller 
  ! term.  This is analogous to finding the scattering factor, where
  ! we sum the square of the anisotropic ("gamma") and isotropic ("a")
  ! parts of the prefactor:
  !
  ! S_k = 45*a_k^2 + 7*gamma_k^2
  !
  ! In the equation below, betaA is the anisotropic part squared and 
  ! polmn is the isotropic part squared.

  write(fname,'(a)') 'RRS_BTerm_Prefactors.out'
  write(fname1,'(a)') 'RRS_BTerm_Prefactor_Comps.out'
  rrs_pref(1:nexci,1:nmodes) = zero
  do iex = 1,nexci
    do j = 1,nmodes
      do a=1,3
        do b=1,3
          pref(a,b) = osc(iex,a)*tdipder(iex,j,b)
        end do
      end do  

      open(80,file=fname1,form='formatted')
      write(80,'(i5,i5)') iex, j
      write(80,'(a14,a14,a14)') 'x','y','z'
      write(80,'(a4,f14.4,f14.4,f14.4)') 'x',(pref(1,b), b=1,3)
      write(80,'(a4,f14.4,f14.4,f14.4)') 'y',(pref(2,b), b=1,3)
      write(80,'(a4,f14.4,f14.4,f14.4)') 'z',(pref(3,b), b=1,3)
      betaA = zero
      polmn = zero
      do a=1,3
        do b=1,3
          betaA = betaA + f12*(three*pref(a,b)*pref(a,b) - &
            pref(a,a)*pref(b,b))
        end do
      end do
    
      do a=1,3
        polmn = polmn + pref(a,a)*pref(a,a)
      end do
      polmn = polmn/nine
      rrs_pref(iex,j) = fortyfv*polmn + seven*betaA
      open(90,file=fname,form='formatted')
      write(90,'(i5,i5,f8.2,f6.2,f14.2)') iex, j, &
                               freq(j), delta(iex,j), &
                               rrs_pref(iex,j)
    end do
    write(90,*)
  end do
  close(80)
  close(90)

End Subroutine RRSPref
