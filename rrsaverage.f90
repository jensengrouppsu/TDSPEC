subroutine RRSFundAvg(pol, rcrs, nmodes)

  use constants

  implicit none

  integer :: m, j, i, nmodes

  complex(kindr), intent(in) :: pol(nmodes,3,3)

  real(kindr) :: betaA, polmn, delta
  real(kindr), intent(out) :: rcrs(nmodes)

  real(kindr), parameter :: thresh = 1.0e-90_kindr

! =============================================================================
! purpose: Orientational averaging of the polarizability for fundamentals
!
! input  : nmodes - number of normal modes
!          pol - RRS polarizability tensor
!
! out : rcrs - RRS scattering factor (perpendicular orientation)
!
! =============================================================================

! This part calculates the average polarizability (I've seen it called
! scattering factor or cross section also).  This is done by taking the square 
! of the symmetric anisotropic  polarizability ("gamma"), antisymmetric anistropic
! polarizability ("delta"), and isotropic polarizability ("a") and summing them 
! as:
!
! S_k = 45*a_k^2 + 7*gamma_k^2 + 5*delta_k^2
!
! If A term scattering is done, the polarizability tensor is symmetric in the
! vibronic model and delta_k is zero.  For overlapping excited states, B term 
! scattering can make the polarizability tensor asymmetric, so delta_k is
! nonzero.
!
! In the equation below, betaA is the anisotropic part of the 
! polarizability squared and polmn is the isotropic part of the 
! polarizability squared.
!
! For the expressions here, see Chapter A14, section 7.5 of the Long book called
! "The Raman Effect".  The coefficients for the tensor invariants are grabbed
! from Table 5.2(a) from the same book.

  do m=1,nmodes
    betaA = zero
    polmn = zero
    delta = zero
    ! Symmetric anistropic part
    do i=1,3
      do j=1,3
        betaA = betaA + f12* &
          ( three * pol(m,i,j) * conjg(pol(m,i,j)) &
          - pol(m,i,i) * conjg(pol(m,j,j)) )
      end do
    end do

    ! Antisymmetric anisotropic part
    ! This expression is easier to write explicitly rather than
    ! using summation notation.  You need to use the traceless 
    ! antisymmetric polarizability tensor to build delta in
    ! summation notation.
    delta = delta + f32*f12* &
          ( ( pol(m,1,2) - pol(m,2,1) ) * &
            conjg( pol(m,1,2) - pol(m,2,1) ) &
          + ( pol(m,2,3) - pol(m,3,2) ) * &
            conjg( pol(m,2,3) - pol(m,3,2) ) &
          + ( pol(m,3,1) - pol(m,1,3) ) * &
            conjg( pol(m,2,3) - pol(m,3,2) ) )

    ! Presumably, this can be very small.  We zero delta here if it
    ! is below a threshold to avoid problems with noise.
    if (abs(delta).lt.thresh) delta = zero

    ! Isotropic part
    do i=1,3
      do j=1,3
        polmn = polmn + pol(m,i,i)*conjg(pol(m,j,j))
      end do
    end do
    !polmn = polmn/nine
    rcrs(m) = fortyfv*polmn + seven*betaA + five*delta
    !rcrs(m) = pol(m,3,3)*conjg(pol(m,3,3)) 
  end do

end subroutine RRSFundAvg

subroutine RRSOverCBAvg(polot, polcb, rcrsot, rcrscb, nmodes, excnum, & 
                        numcb)

  use constants

  implicit none

  integer :: k, m, j, i, nmodes, excnum, numcb

  complex(kindr), intent(in) :: polot(excnum,nmodes,3,3)
  complex(kindr), intent(in) :: polcb(numcb,3,3)

  real(kindr) :: betaA, polmn, delta
  real(kindr), intent(out) :: rcrsot(excnum,nmodes)
  real(kindr), intent(out) :: rcrscb(numcb)

  real(kindr), parameter :: thresh = 1.0e-90_kindr

! =============================================================================
! purpose: Orientational averaging of the polarizability for fundamentals
!
! input  : nmodes - number of normal modes
!          polot - RRS polarizability tensor for overtones
!          polcb - RRS polarizability tensor for comb. bands
!          excnum - user input excitation number
!          numcb - number of combination bands from FCFactor routine
!
! out : rcrsot - RRS scattering factor for overtones (perpendicular orientation)
!       rcrscb - RRS scattering factor for comb. bands (perpendicular orientation)
!
! =============================================================================

! This part calculates the average polarizability (I've seen it called
! scattering factor or cross section also).  This is done by taking the
! of the symmetric anisotropic  polarizability ("gamma"), antisymmetric anistropic
! polarizability ("delta"), and isotropic polarizability ("a") and summing them 
! as:
!
! S_k = 45*a_k^2 + 7*gamma_k^2 + 5*delta_k^2
!
! If A term scattering is done, the polarizability tensor is symmetric in the
! vibronic model and delta_k is zero.  For overlapping excited states, B term 
! scattering can make the polarizability tensor asymmetric, so delta_k is
! nonzero.

  ! Averaging for overtones

  do k=2,excnum
    do m =1,nmodes
      betaA = zero
      polmn = zero
      delta = zero
      ! Symmetric anistropic part
      do i=1,3
        do j=1,3
          betaA = betaA + f12* &
            (three*polot(k,m,i,j)*conjg(polot(k,m,i,j)) - &
            polot(k,m,i,i)*conjg(polot(k,m,j,j)))
        end do
      end do

      ! Antisymmetric anisotropic part
      ! This expression is easier to write explicitly rather than
      ! using summation notation.  You need to use the traceless 
      ! antisymmetric polarizability tensor to build delta in
      ! summation notation.
      delta = delta + f32*f12* &
            ( ( polot(k,m,1,2) - polot(k,m,2,1) ) * &
              conjg( polot(k,m,1,2) - polot(k,m,2,1) ) &
            + ( polot(k,m,2,3) - polot(k,m,3,2) ) * &
              conjg( polot(k,m,2,3) - polot(k,m,3,2) ) &
            + ( polot(k,m,3,1) - polot(k,m,1,3) ) * &
              conjg( polot(k,m,2,3) - polot(k,m,3,2) ) )

      ! Presumably, this can be very small.  We zero delta here if it
      ! is below a threshold to avoid problems with noise.
      if (abs(delta).lt.thresh) delta = zero

      ! Isotropic part
      do i=1,3
        do j=1,3
          polmn = polmn + polot(k,m,i,i)*conjg(polot(k,m,j,j))
        end do
      end do
      polmn = polmn/nine
      rcrsot(k,m) = fortyfv*polmn + seven*betaA + five*delta
    end do
  end do

  ! Averaging for combination bands

  do k=1,numcb
    betaA = zero
    polmn = zero
    delta = zero
    ! Symmetric anisotropic part
    do i=1,3
      do j=1,3
        betaA = betaA + f12* &
          (three*polcb(k,i,j)*conjg(polcb(k,i,j)) - &
          polcb(k,i,i)*conjg(polcb(k,j,j)))
      end do
    end do

    ! Asymmetric anisotropic part
    do i=1,3
      do j=1,3
        delta = delta + f32* &
          ( polcb(k,i,j) * conjg(polcb(k,i,j)) )
      end do
    end do

    ! Antisymmetric anisotropic part
    ! This expression is easier to write explicitly rather than
    ! using summation notation.  You need to use the traceless 
    ! antisymmetric polarizability tensor to build delta in
    ! summation notation.
    delta = delta + f32*f12* &
          ( ( polcb(k,1,2) - polcb(k,2,1) ) * &
            conjg( polcb(k,1,2) - polcb(k,2,1) ) &
          + ( polcb(k,2,3) - polcb(k,3,2) ) * &
            conjg( polcb(k,2,3) - polcb(k,3,2) ) &
          + ( polcb(k,3,1) - polcb(k,1,3) ) * &
            conjg( polcb(k,2,3) - polcb(k,3,2) ) )

    ! Presumably, this can be very small.  We zero delta here if it
    ! is below a threshold to avoid problems with noise.
    if (abs(delta).lt.thresh) delta = zero

    ! Isotropic part
    do i=1,3
      do j=1,3
        polmn = polmn + polcb(k,i,i)*conjg(polcb(k,j,j))
      end do
    end do
    polmn = polmn/nine
    rcrscb(k) = fortyfv*polmn + seven*betaA + five*delta
  end do

end subroutine RRSOverCBAvg

subroutine ASRRSFundAvg(aspol, rcrs, nmodes)

  use constants

  implicit none

  integer :: m, j, i, nmodes

  complex(kindr), intent(in) :: aspol(nmodes,3,3)

  real(kindr) :: betaA, polmn, delta
  real(kindr), intent(out) :: rcrs(nmodes)

  real(kindr), parameter :: thresh = 1.0e-90_kindr

! =============================================================================
! purpose: Orientational averaging of the polarizability for fundamentals
!
! input  : nmodes - number of normal modes
!          pol - RRS polarizability tensor
!
! out : rcrs - anti-Stokes RRS scattering factor (perpendicular orientation)
!
! =============================================================================

! This part calculates the average polarizability (I've seen it called
! scattering factor or cross section also).  This is done by taking the
! of the symmetric anisotropic  polarizability ("gamma"), antisymmetric anistropic
! polarizability ("delta"), and isotropic polarizability ("a") and summing them 
! as:
!
! S_k = 45*a_k^2 + 7*gamma_k^2 + 5*delta_k^2
!
! If A term scattering is done, the polarizability tensor is symmetric in the
! vibronic model and delta_k is zero.  For overlapping excited states, B term 
! scattering can make the polarizability tensor asymmetric, so delta_k is
! nonzero.

  do m=1,nmodes
    betaA = zero
    polmn = zero
    delta = zero
    ! Symmetric anistropic part
    do i=1,3
      do j=1,3
        betaA = betaA + f12* &
          (three*aspol(m,i,j)*conjg(aspol(m,i,j)) - &
          aspol(m,i,i)*conjg(aspol(m,j,j)))
      end do
    end do

    ! Antisymmetric anisotropic part
    ! This expression is easier to write explicitly rather than
    ! using summation notation.  You need to use the traceless 
    ! antisymmetric polarizability tensor to build delta in
    ! summation notation.
    delta = delta + f32*f12* &
          ( ( aspol(m,1,2) - aspol(m,2,1) ) * &
            conjg( aspol(m,1,2) - aspol(m,2,1) ) &
          + ( aspol(m,2,3) - aspol(m,3,2) ) * &
            conjg( aspol(m,2,3) - aspol(m,3,2) ) &
          + ( aspol(m,3,1) - aspol(m,1,3) ) * &
            conjg( aspol(m,2,3) - aspol(m,3,2) ) )

    ! Presumably, this can be very small.  We zero delta here if it
    ! is below a threshold to avoid problems with noise.
    if (abs(delta).lt.thresh) delta = zero

    ! Isotropic part
    do i=1,3
      do j=1,3
        polmn = polmn + aspol(m,i,i)*conjg(aspol(m,j,j))
      end do
    end do
    polmn = polmn/nine
    rcrs(m) = fortyfv*polmn + seven*betaA + five*delta
  end do

end subroutine ASRRSFundAvg
