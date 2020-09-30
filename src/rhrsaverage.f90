subroutine RHRSFundAvg(hpol, rcrs, nmodes)

  use constants

  implicit none

  integer :: m, i, j, k, nmodes

  complex(kindr), intent(in) :: hpol(nmodes,3,3,3)

  real(kindr) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13
  real(kindr) :: b14, b15
  real(kindr) :: bmean1, bmean2
  real(kindr), intent(out) :: rcrs(nmodes)

! =============================================================================
! purpose: Orientational averaging of the hyperpolarizability for fundamentals
!
! input  : nmodes - number of normal modes
!          hpol - RHRS hyperpolarizability tensor
!
! out : rcrs - RHRS scattering factor (perpendicular orientation)
!
! =============================================================================

! Much like what is done for the polarizability, the averages of the
! hyperpolarizability are calculated.  These were done using equations
! in the following references for measurements of scattered light 
! perpendicular to the incident radiation direction:
!
! Quinet et al.  Int. J. Quant. Chem. 106, 599 (2005).
! Yang et al.  J. Chem. Phys 97, 3831 (1992).

  do m=1,nmodes
    b1 = zero
    do i=1,3
      b1 = b1 + hpol(m,i,i,i)*conjg(hpol(m,i,i,i))
    end do

    b2 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b2 = b2 + hpol(m,i,i,j)*conjg(hpol(m,i,i,j))
      end do
    end do

    b3 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b3 = b3 + hpol(m,i,i,i)*conjg(hpol(m,i,j,j))
      end do
    end do

    b4 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b4 = b4 + hpol(m,j,i,i)*conjg(hpol(m,i,i,j))
      end do
    end do

    b5 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b5 = b5 + hpol(m,i,i,i)*conjg(hpol(m,j,j,i))
      end do
    end do

    b6 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b6 = b6 + hpol(m,j,i,i)*conjg(hpol(m,j,i,i))
      end do
    end do

    b7=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b7 = b7 + hpol(m,i,i,j)*conjg(hpol(m,j,k,k))
        end do
      end do
    end do

    b8=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b8 = b8 + hpol(m,j,i,i)*conjg(hpol(m,j,k,k))
        end do
      end do
    end do

    b9=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b9 = b9 + hpol(m,i,i,j)*conjg(hpol(m,k,k,j))
        end do
      end do
    end do

    b10=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b10 = b10 + hpol(m,i,j,k)*conjg(hpol(m,i,j,k))
        end do
      end do
    end do

    b11=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b11 = b11 + hpol(m,i,j,k)*conjg(hpol(m,j,i,k))
        end do
      end do
    end do

    b12 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
        b12 = b12 + hpol(m,i,j,j)*conjg(hpol(m,i,j,j))
      end do
    end do

    b13 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
        b13 = b13 + hpol(m,i,i,j)*conjg(hpol(m,j,i,i))
      end do
    end do

    b14 = zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b14 = b14 + hpol(m,i,j,j)*conjg(hpol(m,i,k,k))
        end do
      end do
    end do

    b15 = zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b15 = b15 + hpol(m,i,i,k)*conjg(hpol(m,j,j,k))
        end do
      end do
    end do

    ! The hyperpolarizability averaging was performed based on the 
    ! equations in: Quinet, O.; Champagne, B.; and 
    ! van Gisbergen, S. J. A.  Int. J. Quantum Chem.  106, 599, 
    ! (2006).

    bmean1 = zero
    bmean1 = b1/seven
    bmean1 = bmean1 + four*b2/thrtfv + two*b3/thrtfv + &
             four*b4/thrtfv + four*b5/thrtfv + b6/thrtfv
    bmean1 = bmean1 + four*b7/hndrfv + b8/hndrfv + &
             four*b9/hndrfv + two*b10/hndrfv + four*b11/hndrfv

    bmean2 = zero
    bmean2 = b1/thrtfv
    bmean2 = bmean2 + four*b3/hndrfv - four*b5/seventy + &
             eight*b2/hndrfv + three*b12/thrtfv - four*b13/seventy
    bmean2 = bmean2 + b14/thrtfv -four*b15/twhndrten - &
             four*b7/twhndrten + two*b10/thrtfv -four*b11/twhndrten

    rcrs(m) = bmean1 + bmean2
  end do

end subroutine RHRSFundAvg

subroutine RHRSOverCBAvg(hpolot, hpolcb, rcrsot, rcrscb, nmodes, &
                         excnum, numcb)

  use constants

  implicit none

  integer :: m, i, j, k, r, nmodes, excnum, numcb

  complex(kindr), intent(in) :: hpolot(excnum,nmodes,3,3,3)
  complex(kindr), intent(in) :: hpolcb(numcb,3,3,3)

  real(kindr) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13
  real(kindr) :: b14, b15
  real(kindr) :: bmean1, bmean2
  real(kindr), intent(out) :: rcrsot(excnum,nmodes)
  real(kindr), intent(out) :: rcrscb(numcb)

! =============================================================================
! purpose: Orientational averaging of the hyperpolarizability for overtones and
!          combination bands
!
! input  : nmodes - number of normal modes
!          hpolot - RHRS hyperpolarizability tensor for overtones
!          hpolcb - RHRS hyperpolarizability tensor for comb. bands
!          excnum - excitation number
!          numcb - number of combination bands from FCFactor routine
!
! out : rcrsot - RHRS scattering factor for overtones (perpendicular orientation)
!       rcrscb - RHRS scattering factor for comb. bands (perpendicular orientation)
!
! =============================================================================

! Much like what is done for the polarizability, the averages of the
! hyperpolarizability are calculated.  These were done using equations
! in the following references for measurements of scattered light 
! perpendicular to the incident radiation direction:
!
! Quinet et al.  Int. J. Quant. Chem. 106, 599 (2005).
! Yang et al.  J. Chem. Phys 97, 3831 (1992).

  ! Perform averaging for overtones
  do r=2,excnum
    do m=1,nmodes
      b1 = zero
      do i=1,3
        b1 = b1 + hpolot(r,m,i,i,i)*conjg(hpolot(r,m,i,i,i))
      end do

      b2 = zero
      do i=1,3
        do j=1,3
          if (i==j) cycle
            b2 = b2 + hpolot(r,m,i,i,j)*conjg(hpolot(r,m,i,i,j))
        end do
      end do

      b3 = zero
      do i=1,3
        do j=1,3
          if (i==j) cycle
            b3 = b3 + hpolot(r,m,i,i,i)*conjg(hpolot(r,m,i,j,j))
        end do
      end do

      b4 = zero
      do i=1,3
        do j=1,3
          if (i==j) cycle
            b4 = b4 + hpolot(r,m,j,i,i)*conjg(hpolot(r,m,i,i,j))
        end do
      end do

      b5 = zero
      do i=1,3
        do j=1,3
          if (i==j) cycle
            b5 = b5 + hpolot(r,m,i,i,i)*conjg(hpolot(r,m,j,j,i))
        end do
      end do

      b6 = zero
      do i=1,3
        do j=1,3
          if (i==j) cycle
            b6 = b6 + hpolot(r,m,j,i,i)*conjg(hpolot(r,m,j,i,i))
        end do
      end do

      b7=zero
      do i=1,3
        do j=1,3
          do k=1,3
            if (i==j.or.i==k.or.j==k) cycle
            b7 = b7 + hpolot(r,m,i,i,j)*conjg(hpolot(r,m,j,k,k))
          end do
        end do
      end do

      b8=zero
      do i=1,3
        do j=1,3
          do k=1,3
            if (i==j.or.i==k.or.j==k) cycle
            b8 = b8 + hpolot(r,m,j,i,i)*conjg(hpolot(r,m,j,k,k))
          end do
        end do
      end do

      b9=zero
      do i=1,3
        do j=1,3
          do k=1,3
            if (i==j.or.i==k.or.j==k) cycle
            b9 = b9 + hpolot(r,m,i,i,j)*conjg(hpolot(r,m,k,k,j))
          end do
        end do
      end do

      b10=zero
      do i=1,3
        do j=1,3
          do k=1,3
            if (i==j.or.i==k.or.j==k) cycle
            b10 = b10 + hpolot(r,m,i,j,k)*conjg(hpolot(r,m,i,j,k))
          end do
        end do
      end do

      b11=zero
      do i=1,3
        do j=1,3
          do k=1,3
            if (i==j.or.i==k.or.j==k) cycle
            b11 = b11 + hpolot(r,m,i,j,k)*conjg(hpolot(r,m,j,i,k))
          end do
        end do
      end do

      b12 = zero
      do i=1,3
        do j=1,3
          if (i==j) cycle
          b12 = b12 + hpolot(r,m,i,j,j)*conjg(hpolot(r,m,i,j,j))
        end do
      end do

      b13 = zero
      do i=1,3
        do j=1,3
          if (i==j) cycle
          b13 = b13 + hpolot(r,m,i,i,j)*conjg(hpolot(r,m,j,i,i))
        end do
      end do

      b14 = zero
      do i=1,3
        do j=1,3
          do k=1,3
            if (i==j.or.i==k.or.j==k) cycle
            b14 = b14 + hpolot(r,m,i,j,j)*conjg(hpolot(r,m,i,k,k))
          end do
        end do
      end do

      b15 = zero
      do i=1,3
        do j=1,3
          do k=1,3
            if (i==j.or.i==k.or.j==k) cycle
            b15 = b15 + hpolot(r,m,i,i,k)*conjg(hpolot(r,m,j,j,k))
          end do
        end do
      end do

      bmean1 = zero
      bmean1 = b1/seven
      bmean1 = bmean1 + four*b2/thrtfv + two*b3/thrtfv + &
               four*b4/thrtfv + four*b5/thrtfv + b6/thrtfv
      bmean1 = bmean1 + four*b7/hndrfv + b8/hndrfv + &
               four*b9/hndrfv + two*b10/hndrfv + four*b11/hndrfv

      bmean2 = zero
      bmean2 = b1/thrtfv
      bmean2 = bmean2 + four*b3/hndrfv - four*b5/seventy + &
               eight*b2/hndrfv + three*b12/thrtfv - &
               four*b13/seventy
      bmean2 = bmean2 + b14/thrtfv -four*b15/twhndrten - &
               four*b7/twhndrten + two*b10/thrtfv -four*b11/twhndrten

      rcrsot(r,m) = bmean1 + bmean2
    end do
  end do

  ! Perform averaging for combination bands
  do m=1,numcb
    b1 = zero
    do i=1,3
      b1 = b1 + hpolcb(m,i,i,i)*conjg(hpolcb(m,i,i,i))
    end do

    b2 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b2 = b2 + hpolcb(m,i,i,j)*conjg(hpolcb(m,i,i,j))
      end do
    end do

    b3 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b3 = b3 + hpolcb(m,i,i,i)*conjg(hpolcb(m,i,j,j))
      end do
    end do

    b4 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b4 = b4 + hpolcb(m,j,i,i)*conjg(hpolcb(m,i,i,j))
      end do
    end do

    b5 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b5 = b5 + hpolcb(m,i,i,i)*conjg(hpolcb(m,j,j,i))
      end do
    end do

    b6 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b6 = b6 + hpolcb(m,j,i,i)*conjg(hpolcb(m,j,i,i))
      end do
    end do

    b7=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b7 = b7 + hpolcb(m,i,i,j)*conjg(hpolcb(m,j,k,k))
        end do
      end do
    end do

    b8=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b8 = b8 + hpolcb(m,j,i,i)*conjg(hpolcb(m,j,k,k))
        end do
      end do
    end do

    b9=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b9 = b9 + hpolcb(m,i,i,j)*conjg(hpolcb(m,k,k,j))
        end do
      end do
    end do

    b10=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b10 = b10 + hpolcb(m,i,j,k)*conjg(hpolcb(m,i,j,k))
        end do
      end do
    end do

    b11=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b11 = b11 + hpolcb(m,i,j,k)*conjg(hpolcb(m,j,i,k))
        end do
      end do
    end do

    b12 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
        b12 = b12 + hpolcb(m,i,j,j)*conjg(hpolcb(m,i,j,j))
      end do
    end do

    b13 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
        b13 = b13 + hpolcb(m,i,i,j)*conjg(hpolcb(m,j,i,i))
      end do
    end do

    b14 = zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b14 = b14 + hpolcb(m,i,j,j)*conjg(hpolcb(m,i,k,k))
        end do
      end do
    end do

    b15 = zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b15 = b15 + hpolcb(m,i,i,k)*conjg(hpolcb(m,j,j,k))
        end do
      end do
    end do

    ! The hyperpolarizability averaging was performed based on the 
    ! equations in: Quinet, O.; Champagne, B.; and 
    ! van Gisbergen, S. J. A.  Int. J. Quantum Chem.  106, 599, 
    ! (2006).

    bmean1 = zero
    bmean1 = b1/seven
    bmean1 = bmean1 + four*b2/thrtfv + two*b3/thrtfv + &
             four*b4/thrtfv + four*b5/thrtfv + b6/thrtfv
    bmean1 = bmean1 + four*b7/hndrfv + b8/hndrfv + &
             four*b9/hndrfv + two*b10/hndrfv + four*b11/hndrfv

    bmean2 = zero
    bmean2 = b1/thrtfv
    bmean2 = bmean2 + four*b3/hndrfv - four*b5/seventy + &
             eight*b2/hndrfv + three*b12/thrtfv - four*b13/seventy
    bmean2 = bmean2 + b14/thrtfv -four*b15/twhndrten - &
             four*b7/twhndrten + two*b10/thrtfv -four*b11/twhndrten

    rcrscb(m) = bmean1 + bmean2
  end do

end subroutine RHRSOverCBAvg

subroutine ASRHRSFundAvg(ashpol, rcrs, nmodes)

  use constants

  implicit none

  integer :: m, i, j, k, nmodes

  complex(kindr), intent(in) :: ashpol(nmodes,3,3,3)

  real(kindr) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13
  real(kindr) :: b14, b15
  real(kindr) :: bmean1, bmean2
  real(kindr), intent(out) :: rcrs(nmodes)

! =============================================================================
! purpose: Orientational averaging of the anti-Stokes hyperpolarizability for 
!          fundamentals
!
! input  : nmodes - number of normal modes
!          ashpol - RHRS hyperpolarizability tensor
!
! out : rcrs - RHRS scattering factor (perpendicular orientation)
!
! =============================================================================

! Much like what is done for the polarizability, the averages of the
! hyperpolarizability are calculated.  These were done using equations
! in the following references for measurements of scattered light 
! perpendicular to the incident radiation direction:
!
! Quinet et al.  Int. J. Quant. Chem. 106, 599 (2005).
! Yang et al.  J. Chem. Phys 97, 3831 (1992).

  do m=1,nmodes
    b1 = zero
    do i=1,3
      b1 = b1 + ashpol(m,i,i,i)*conjg(ashpol(m,i,i,i))
    end do

    b2 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b2 = b2 + ashpol(m,i,i,j)*conjg(ashpol(m,i,i,j))
      end do
    end do

    b3 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b3 = b3 + ashpol(m,i,i,i)*conjg(ashpol(m,i,j,j))
      end do
    end do

    b4 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b4 = b4 + ashpol(m,j,i,i)*conjg(ashpol(m,i,i,j))
      end do
    end do

    b5 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b5 = b5 + ashpol(m,i,i,i)*conjg(ashpol(m,j,j,i))
      end do
    end do

    b6 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
          b6 = b6 + ashpol(m,j,i,i)*conjg(ashpol(m,j,i,i))
      end do
    end do

    b7=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b7 = b7 + ashpol(m,i,i,j)*conjg(ashpol(m,j,k,k))
        end do
      end do
    end do

    b8=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b8 = b8 + ashpol(m,j,i,i)*conjg(ashpol(m,j,k,k))
        end do
      end do
    end do

    b9=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b9 = b9 + ashpol(m,i,i,j)*conjg(ashpol(m,k,k,j))
        end do
      end do
    end do

    b10=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b10 = b10 + ashpol(m,i,j,k)*conjg(ashpol(m,i,j,k))
        end do
      end do
    end do

    b11=zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b11 = b11 + ashpol(m,i,j,k)*conjg(ashpol(m,j,i,k))
        end do
      end do
    end do

    b12 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
        b12 = b12 + ashpol(m,i,j,j)*conjg(ashpol(m,i,j,j))
      end do
    end do

    b13 = zero
    do i=1,3
      do j=1,3
        if (i==j) cycle
        b13 = b13 + ashpol(m,i,i,j)*conjg(ashpol(m,j,i,i))
      end do
    end do

    b14 = zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b14 = b14 + ashpol(m,i,j,j)*conjg(ashpol(m,i,k,k))
        end do
      end do
    end do

    b15 = zero
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          b15 = b15 + ashpol(m,i,i,k)*conjg(ashpol(m,j,j,k))
        end do
      end do
    end do

    ! The hyperpolarizability averaging was performed based on the 
    ! equations in: Quinet, O.; Champagne, B.; and 
    ! van Gisbergen, S. J. A.  Int. J. Quantum Chem.  106, 599, 
    ! (2006).

    bmean1 = zero
    bmean1 = b1/seven
    bmean1 = bmean1 + four*b2/thrtfv + two*b3/thrtfv + &
             four*b4/thrtfv + four*b5/thrtfv + b6/thrtfv
    bmean1 = bmean1 + four*b7/hndrfv + b8/hndrfv + &
             four*b9/hndrfv + two*b10/hndrfv + four*b11/hndrfv

    bmean2 = zero
    bmean2 = b1/thrtfv
    bmean2 = bmean2 + four*b3/hndrfv - four*b5/seventy + &
             eight*b2/hndrfv + three*b12/thrtfv - four*b13/seventy
    bmean2 = bmean2 + b14/thrtfv -four*b15/twhndrten - &
             four*b7/twhndrten + two*b10/thrtfv -four*b11/twhndrten

    rcrs(m) = bmean1 + bmean2
  end do

end subroutine ASRHRSFundAvg
