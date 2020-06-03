subroutine HPolAvg(hpol, freq, nmodes, lfreq, lherztell, &
                   lfullherztell, latermoff, lrhrsb2)

  use constants

  implicit none

  character*80 :: frmt
  character*80 :: fname
  character*80 :: charfcht
  character*80 :: charfreq

  integer :: m, i, j, nmodes
  integer :: intfreq

  complex(kindr), intent(in) :: hpol(nmodes,3,3,3)
  complex(kindr) :: betaavg(nmodes,3)
  complex(kindr) :: betanorm(nmodes)

  real(kindr), intent(in) :: freq(nmodes)
  real(kindr), intent(in) :: lfreq

  complex(kindr) :: zzero=(0d0,0d0)

  logical :: lherztell
  logical :: lfullherztell
  logical :: latermoff
  logical :: lrhrsb2

! =============================================================================
! purpose: Build orientational averages of the hyperpolarizability to get 
!         Cartesian axis oriented hyperpolarizabilities.
!
! input  : nmodes - number of normal modes
!          hpol - RHRS hyperpolarizability tensor
!          freq - Vibrational frequency
!          lfreq - Incident frequency
!          lherztell - Inclusion of B1 terms only
!          lfullherztell - Inclusion of full B terms
!          latermoff - Neglect the A term
!          lrhrsb2 - Inclusion of B2 terms only
!
! =============================================================================

! Here we perform a specific orientational average of the hyperpolarizability:
!
! |\beta_{i}^{p}|^2 = 
!       (1/5)\sum_{j=x,y,z}(|\beta_{ijj}^{p}|^2
!                          +|\beta_{jij}^{p}|^2
!                          +|\beta_{jji}^{p}|^2)
!
! i is one of the Cartesian directions (x, y, z).  Some authors also call
! this beta-parallel, but then they choose a specific axis (usually z).
! We must somehow account for the frequency-dependent hyperpolarizablity
! being complex, since this averaging is in general done for the
! static (real) hyperpolarizability.  Since we're using SI units for
! the hyperpolarizability, the numbers are ~10^{-50} so since we're
! printing absolute squares, we convert to a.u. and multiply by 10^{-8} 

  betaavg = zzero

  do m=1,nmodes
    do i=1,3
      do j=1,3
        betaavg(m,i) = betaavg(m,i) + &
                       hpol(m,i,j,j)*conjg(hpol(m,i,j,j)) + &
                       hpol(m,j,i,j)*conjg(hpol(m,j,i,j)) + &
                       hpol(m,j,j,i)*conjg(hpol(m,j,j,i))
      end do
      betaavg(m,i) = f15*tentomeight*betaavg(m,i)/(au2si**2)
    end do
    betanorm(m) = betaavg(m,1) + betaavg(m,2) + betaavg(m,3)
  end do

  if (lfullherztell.and.(latermoff.eqv..false.).and. &
    (lrhrsb2.eqv..false.)) then
    charfcht = 'AB'
  else if (lherztell.and.(latermoff.eqv..false.)) then
    charfcht = 'AB1'
  else if (lherztell.and.latermoff) then
    charfcht = 'B1'
  else if (lfullherztell.and.lrhrsb2.and.(latermoff.eqv..false.)) then
    charfcht = 'AB2'
  else if (lrhrsb2.and.latermoff) then
    charfcht = 'B2'
  else if (lfullherztell.and.latermoff) then
    charfcht = 'B'
  else
    charfcht = 'A'
  end if

  ! Convert the laser frequency to a character.
  intfreq = int(2.0E+7_kindr/lfreq)
  write(charfreq,*) intfreq
  charfreq = adjustl(charfreq)

  charfcht = adjustl(charfcht)

  fname = 'Beta_Cartesian_exc' // trim(charfreq) // '_' // &
          trim(charfcht) // 'term.out'

  frmt = '(F7.2,2X,ES12.6,2X,ES12.6,2X,ES12.6,2X,ES12.6)'

  open(2000,file=fname,form='formatted')

  write(2000,'(51A)') '#Squared hyperpolarizabilities printed with units:'  
  write(2000,'(14A)') '#x10^8 a.u.^2.'  
  write(2000,'(A)')   '#'  
  write(2000,'(20A)') '#Useful conversions:'  
  write(2000,'(58A)') '# 1 a.u. = 8.6392E-33 esu, where 1 esu = 1 cm^5 statC^{-1}'  
  write(2000,'(32A)') '# 1 a.u. = 3.6213E-42 m^4 V^{-1}'  
  write(2000,'(35A)') '# 1 a.u. = 3.206E-53 C^3 m^3 J^{-2}'  
  write(2000,'(A)')   '#'  
  write(2000,'(65A)') '#Freq     |Beta_x|^2    |Beta_y|^2    |Beta_z|^2    |Beta_norm|^2'  

  do m=1,nmodes
    write(2000,frmt) freq(m), real(betaavg(m,1)), real(betaavg(m,2)), & 
      real(betaavg(m,3)), real(betanorm(m))
  end do

  write(2000,*)
  write(2000,'(33A)') '#Normalized hyperpolarizabilities'  
  write(2000,'(65A)') '#Freq     |Beta_x|^2    |Beta_y|^2    |Beta_z|^2    |Beta_norm|^2'  
  
  do m=1,nmodes
    if (real(betanorm(m)).eq.zero) then
      write(2000,frmt) freq(m), zero, zero, zero, real(betanorm(m))
    else 
      write(2000,frmt) freq(m), real(betaavg(m,1))/real(betanorm(m)), &
        real(betaavg(m,2))/real(betanorm(m)), &
        real(betaavg(m,3))/real(betanorm(m)), real(betanorm(m))
    end if
  end do

  close(2000)

end subroutine HPolAvg
