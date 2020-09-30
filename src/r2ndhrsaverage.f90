subroutine R2NDHRSFundAvg(sechpol, rcrs, nmodes)

  use constants

  implicit none

  integer :: m, i, j, k, nmodes

  complex(kindr), intent(in) :: sechpol(nmodes,3,3,3,3)

  real(kindr) :: g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13
  real(kindr) :: g14, g15, g16, g17, g18, g19, g20, g21, g22, g23, g24
  real(kindr) :: g25, g26, g27, g28, g29
  real(kindr) :: giiii1, giiii2, giiii3, giiii4
  real(kindr) :: gijjj1, gijjj2, gijjj3, gijjj4
  real(kindr), intent(out) :: rcrs(nmodes)

! ===================================================================================
! purpose: Orientational averaging of the second hyperpolarizability for fundamentals
!
! input  : nmodes - number of normal modes
!          sechpol - RSHRS second hyperpolarizability tensor
!
! out : rcrs - RSHRS scattering factor (perpendicular orientation)
!
! ===================================================================================

  do m=1,nmodes
    ! Initialize all quantities
    g1 = zero
    g2 = zero
    g3 = zero
    g4 = zero
    g5 = zero
    g6 = zero
    g7 = zero
    g8 = zero
    g9 = zero
    g10 = zero
    g11 = zero
    g12 = zero
    g13 = zero
    g14 = zero
    g15 = zero
    g16 = zero
    g17 = zero
    g18 = zero
    g19 = zero
    g20 = zero
    g21 = zero
    g22 = zero
    g23 = zero
    g24 = zero
    g25 = zero
    g26 = zero
    g27 = zero
    g28 = zero
    g29 = zero

    do i=1,3
      g1 = g1 + sechpol(m,i,i,i,i)*conjg(sechpol(m,i,i,i,i))
    end do

    ! G1 term in my notes
    giiii1 = one*g1/nine 
    gijjj1 = one*g1/sixtythree

    do i=1,3
      do j=1,3
        if (i==j) cycle
        g2 = g2 + sechpol(m,i,i,i,j)*conjg(sechpol(m,i,i,i,j))
        g3 = g3 + sechpol(m,i,i,i,i)*conjg(sechpol(m,i,i,j,j))
        g4 = g4 + sechpol(m,i,i,i,j)*conjg(sechpol(m,j,i,i,i))
        g5 = g5 + sechpol(m,j,i,i,i)*conjg(sechpol(m,j,i,i,i))
        g6 = g6 + sechpol(m,i,i,i,i)*conjg(sechpol(m,j,i,i,j))
      end do
    end do

    ! G2 term in my notes
    giiii2 = one*g2/seven + two*g3/twentyone + two*g4/twentyone + &
             one*g5/sixtythree + two*g6/twentyone
    gijjj2 = two*g2/thrtfv + four*g3/hndrfv - one*g4/twentyone + &
             four*g5/sixtythree - one*g6/twentyone

    do i=1,3
      do j=1,3
        if (i==j) cycle
        g7 = g7 + sechpol(m,i,i,j,j)*conjg(sechpol(m,i,i,j,j))
        g8 = g8 + sechpol(m,i,i,i,j)*conjg(sechpol(m,i,j,j,j))
        g9 = g9 + sechpol(m,i,j,j,j)*conjg(sechpol(m,j,i,i,i))
        g10 = g10 + sechpol(m,i,i,j,j)*conjg(sechpol(m,j,i,i,j))
        g11 = g11 + sechpol(m,j,i,i,j)*conjg(sechpol(m,j,i,i,j))
        g12 = g12 + sechpol(m,j,i,i,i)*conjg(sechpol(m,j,i,j,j))
        g13 = g13 + sechpol(m,i,i,i,j)*conjg(sechpol(m,j,j,j,i))
        g14 = g14 + sechpol(m,i,i,i,i)*conjg(sechpol(m,j,j,j,j))
      end do
    end do

    ! G3 term in my notes
    giiii3 = three*g7/thrtfv + two*g8/thrtfv + two*g9/hndrfv + &
             six*g10/thrtfv + three*g11/thrtfv + two*g12/thrtfv + &
             six*g13/thrtfv + two*g14/hndrfv
    gijjj3 = three*g7/thrtfv + two*g8/thrtfv - one*g9/hndrfv - &
             three*g10/thrtfv + three*g11/thrtfv - three*g12/thrtfv + &
             two*g13/thrtfv - one*g14/hndrfv

    do i=1,3
      do j=1,3
        do k=1,3
          if (i==j.or.i==k.or.j==k) cycle
          g15 = g15 + sechpol(m,i,i,j,k)*conjg(sechpol(m,i,i,j,k))
          g16 = g16 + sechpol(m,i,i,j,j)*conjg(sechpol(m,i,i,k,k))
          g17 = g17 + sechpol(m,i,i,i,k)*conjg(sechpol(m,i,j,j,k))
          g18 = g18 + sechpol(m,i,j,k,k)*conjg(sechpol(m,j,i,i,i))
          g19 = g19 + sechpol(m,i,i,k,k)*conjg(sechpol(m,j,i,i,j))
          g20 = g20 + sechpol(m,i,i,j,k)*conjg(sechpol(m,j,i,i,k))
          g21 = g21 + sechpol(m,j,i,i,k)*conjg(sechpol(m,j,i,i,k))
          g22 = g22 + sechpol(m,i,i,i,k)*conjg(sechpol(m,j,i,j,k))
          g23 = g23 + sechpol(m,i,i,i,j)*conjg(sechpol(m,j,i,k,k))
          g24 = g24 + sechpol(m,j,i,i,i)*conjg(sechpol(m,j,i,k,k))
          g25 = g25 + sechpol(m,i,i,i,i)*conjg(sechpol(m,j,j,k,k))
          g26 = g26 + sechpol(m,j,i,j,k)*conjg(sechpol(m,k,i,i,i))
          g27 = g27 + sechpol(m,j,i,i,k)*conjg(sechpol(m,k,i,i,j))
          g28 = g28 + sechpol(m,j,i,i,j)*conjg(sechpol(m,k,i,i,k))
          g29 = g29 + sechpol(m,j,i,i,i)*conjg(sechpol(m,k,i,j,k))
        end do
      end do
    end do

    ! G4 term in my notes
    giiii4 = two*g15/thrtfv + one*g16/thrtfv + two*g17/thrtfv + &
             two*g18/hndrfv + two*g19/thrtfv + four*g20/thrtfv + &
             one*g21/hndrfv + four*g22/thrtfv + two*g23/thrtfv + &
             two*g24/hndrfv + one*g25/sixtythree + four*g26/hndrfv + &
             one*g27/thrtfv + one*g28/thrtfv + one*g29/thrhndrftn
    gijjj4 = two*g15/thrtfv + one*g16/thrtfv + two*g17/thrtfv - &
             one*g18/hndrfv - one*g19/thrtfv - two*g20/thrtfv + &
             four*g21/hndrfv - two*g22/thrtfv - one*g23/thrtfv + &
             eight*g24/hndrfv - one*g25/hndrfv - one*g26/twhndrten - &
             one*g27/seventy - one*g28/seventy - one*g29/seventy

    rcrs(m) = giiii1 + giiii2 + giiii3 + giiii4 + &
              gijjj1 + gijjj2 + gijjj3 + gijjj4 
  end do

end subroutine R2NDHRSFundAvg
