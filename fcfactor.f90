!Subroutine FCfactor(freq, delta, excnum, fc_screen, &
!              freq_screen, freqcb, excnums, k, nexci, &
!              nmodes, numcb)
Subroutine FCfactor(freq, delta, excnum, fc_screen, &
              freq_screen, freqcbtmp, excnumstmp, nexci, &
              nmodes, numcb, cbguess)

  use fact

  Implicit none

!----------------------------------------
! Double precision
!----------------------------------------
  integer, parameter :: kindr = kind(0d0)
!----------------------------------------
  integer,intent(in) :: nmodes
  integer,intent(in) :: nexci
  integer,intent(in) :: excnum
  integer,intent(in) :: cbguess
  integer :: lmax
!  integer :: s
!  integer :: t
  integer :: q(nmodes)
  integer :: imax 
  integer :: i 
!  integer :: iex
  integer :: pmax
  integer :: first,last
  integer :: j,k,l,p,r,rtmp
!  integer :: m 
  integer :: countzero = 0
  integer,intent(out) :: numcb

  real(kindr),intent(in) :: fc_screen 
  real(kindr),intent(in) :: freq_screen 

  real(kindr), parameter :: zero  = 0.0_kindr
  real(kindr) :: freq(nmodes)
  real(kindr) :: delta(nexci,nmodes)
  real(kindr) :: fc(nmodes,0:15)
!  real(kindr) :: fc(nexci,nmodes,0:10)
  real(kindr) :: fctot, fctmp
!  real(kindr) :: fc0
  real(kindr) :: cb_freq
  real(kindr), intent(inout) :: freqcbtmp(cbguess)
!  real(kindr), intent(out) :: freqcb(nexci,100)

  integer, external :: position

  integer :: ncmb(nmodes+excnum,excnum+1)
  integer :: ntot(excnum+1)
!  integer, allocatable :: ncmb(:,:),ntot(:)
  integer :: n(nmodes)
!  integer, allocatable :: n(:)
  integer,intent(inout) :: excnumstmp(cbguess,nmodes)
!  integer, intent(inout) :: excnums(nexci,s,t)

!----------------------------------------

!===========================================================
! purpose: Calculate FC factors for determining the relative
! importance of different combination bands, used for making
! different types of spectra with TDSPEC.
!
! input  : freq - normal mode frequencies
!          delta - dimensionless deltas
!          nexci - number of excited states
!          nmodes - number of normal modes
!          excnum - maximum excitation level requested
!          fc_screen - FC factor screening threshold
!          freq_screen - frequency screening threshold
!
! output : freqcb - frequencies of the combination bands
!          excnums - excitation numbers for comb. bands.
!          numcb - number of combination bands
!===========================================================

  ! Define/intialize variables

  imax = nmodes
  lmax = excnum
  numcb = 0

  ! Allocate arrays storing combination bands and excitation numbers.

  fc = zero

! f(j) comes from the module fact, defining the factorial j!.  This part
! determines the FC factor using the formula for a Poisson distribution:
!
! P(X) = e^{-mu}*mu^x / x!
!
! Translating this to Chemistry terms, mu is related to the intensity
! of the transition (in this case, the dimensionless displacement), and x is
! the excitation number (number of vibrational quanta).

! Fix this to work with multiple excited states!!!!!!!! (third loop)

  do i=1,imax
    do j=0,15
      fc(i,j) = exp(-delta(1,i))*(delta(1,i)**j)/(f(j))   
    end do
  end do

! Include multiple excited states.

!  do k=1,nexci
!    do i=1,imax
!      do j=0,15
!        fc(k,i,j) = exp(-delta(k,i))*(delta(k,i)**j)/(f(j))
!      end do
!    end do
!  end do

! ---------------------
! calculate the 0-0 FCF
! ---------------------
!
! Determine the total FC factor, by multiplying the contribution from all of 
! the normal modes.
!
  fctot = 0
  fctmp = 1.0
  q = 0
  do i=1,imax
    fctmp = fctmp*fc(i,q(i))
  end do 
  fctot = fctot + fctmp

!  do k=1,nexci
!    fctot = 0
!    fctmp = 1.0
!    q = 0
!    do i=1,imax
!      fctmp = fctmp*fc(k,i,q(i))
!    end do
!    fctot = fctot + fctmp
!  end do
 
! ---------------------------------------------------
! Calculate the number of quantum number combintaions
! ---------------------------------------------------
!
! Determine the value of S(alpha,l) from Ruhoff and Ratner's paper.  
! For multiple excited states this does not change.
!
!  allocate(ncmb(imax+lmax,lmax+1))

  ! Define the binomial coefficients, defined by the formula:
  ! B(alpha+l-1,l) = (alpha+l-1)! / (alpha-l)!l!
  !
  ! Here, alpha is the number of normal modes, l is the excitation
  ! state level.

  ! Here, the dummy variable "i" defines the number of modes times the
  ! maximum excitation state level, "j" defines the number of excitation
  ! levels.

  do i=1,imax+lmax
    do j=1,lmax+1
      ncmb(i,j) = 1
      if (j <= i) then
        do k=1,j
          ncmb(i,j) = (ncmb(i,j)*(i-k+1)) / k
        end do
      end if
    end do
  end do

  ! Total number of nodes, N(alpha,l), from Ruhoff and Ratner's paper.

!  allocate(ntot(lmax+1))

  ntot = 1
  do i = 1, lmax
    do j = 1, i
      ntot(i + 1) = ntot(i + 1) + ncmb(imax + j - 1, j)
    end do
  end do

  ! Total number of permutations

  pmax = ntot(lmax+1)

  ! -------------------------------
  ! Calculate the overlap integrals
  ! -------------------------------

  ! Determine the FC overlap integrals, in accordance with the example code
  ! on p. 387 of Ruhoff and Ratner's paper.

  r = 1
  rtmp = r

  ! This part defines the array containing the excitation levels, n(#)

!  allocate (n(imax))
!  do m=1,nexci
  ! Initially n is zero.
  n = 0
  ! Loop over all the states. 
  n(imax) = 1
  l = 1
  ! The variables first and last are changed from imax to different
  ! values after the first iteration of the loop.  This performs the 
  ! loop over every possible permutation
  first = imax 
  last = imax 
  do i =1, pmax -1
  
    ! Subtract 1 from n(#) (initially zero).  Calculate the 
    ! contribution from level l-1
  
    n(first) = n(first) - 1 
  
    ! Calculate the position of the FC index tree.
  
    p = position(n,imax, l-1, lmax, ncmb, ntot)
  
    ! Add 1 to n(#) (set to one initially)
  
    n(first) = n(first) + 1
  
    do k=first+1,imax
      ! This is only performed if a vibrational quantum number is
      ! not zero (can't lower past the ground state). 
      if(n(k) > 0) then
        n(k) = n(k) -1
        p = position(n,imax, l-1, lmax, ncmb, ntot)
        n(k) = n(k) +1
      end if
    end do
    ! Calculate the contribution from level l-2
    if (n(first) > 1) then
      n(first) = n(first) -2
      p = position(n,imax, l-2, lmax, ncmb, ntot)
      n(first) = n(first) + 2
    end if
  
    ! Calculate the FC overlap integrals.
    fctmp = 1
    do j=1,imax
      fctmp = fctmp*fc(j,n(j))
!      fctmp = fctmp*fc(m,j,n(j))
    end do
  
    ! Attempt at storing the valuable information.  First check if the
    ! overlap integral is greater than the threshold value.
  
    if (fctmp > fc_screen) then
      cb_freq=0.0_kindr
      do j=1,imax
        cb_freq = cb_freq + n(j)*freq(j)
      end do
      ! Test if the vibrational frequency corresponds to a fundamental
      ! or overtone.  If not, store the frequency (this only stores 
      ! combination bands).  
      do j=1,imax
        if (n(j) == 0) then
          countzero = countzero + 1
        end if
      end do
      if (countzero /= imax - 1) then
        freqcbtmp(rtmp) = cb_freq
!        write(*,*) freqcbtmp(rtmp)
!        freqcb(m,rtmp) = cb_freq
      end if
      ! Store the excitation numbers for the combination band, if the
      ! frequency is less than the threshold.
      if (cb_freq < freq_screen) then
!        write(*,*) i
!        if (countzero /= imax - 1) then
!          write(*,*) freqcbtmp(rtmp)
!        end if
        do j=1,imax
          if (countzero /= imax - 1) then
            excnumstmp(rtmp,j) = n(j)
!            if (excnumstmp(rtmp,j) /= 0) then
!              write(*,*) j, excnumstmp(rtmp,j)
!            end if
!            write(*,*) j, n(j)
!            excnums(m,rtmp,j) = n(j)
          end if
        end do
!        write(*,*)
        if (countzero /= imax - 1) then
          numcb = numcb + 1
        end if
        if (countzero /= imax - 1) then
          rtmp = rtmp+1
        end if
      end if
    end if
  
    countzero = 0
    fctot = fctot + fctmp
  
    ! The subroutine nextstate changes l, first, and last, defining
    ! excitation level, first nonzeo element of n(#) and last
    ! nonzero element of n(#), respectively.  
  
    call nextstate(n, imax, l, first, last)
  end do
!  end do

end subroutine FCfactor
