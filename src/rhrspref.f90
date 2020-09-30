Subroutine RHRSPrefB1(rhrs_pref, stpm, tdipder, delta, freq, nexci, &
                      nmodes, tdim, sdim) 

  use constants

  implicit none

  integer :: nexci, nmodes, tdim, sdim, iex, j
  integer :: a, b, c

  character(80) :: fname, fname1

  real(kindr), intent(in) :: stpm(nexci,sdim)
  real(kindr), intent(in) :: tdipder(nexci,nmodes,tdim)
  real(kindr), intent(inout) :: rhrs_pref(nexci,nmodes)
  real(kindr) :: bmean1, bmean2
  real(kindr) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11
  real(kindr) :: b12, b13, b14, b15
  real(kindr) :: pref(3,3,3) 
  real(kindr) :: delta(nexci,nmodes), freq(nmodes)

! =============================================================================
! purpose: Orientational averaging of the B1 term prefactor for RHRS
!
! input  : rhrs_pref - orientational averages of the B1 term prefactors
!          freq - vibrational frequencies
!          stpm - two-photon transition moment components
!          tdipder - derivatives of transition dipole moments along each mode
!          delta - dimensionless displacements
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tdim - final dimension of the transition dipole moment derivative
!          array (equals 3)
!          sdim - final dimension of the two-photon transition moment 
!          array (equals 6)
!
! =============================================================================


  !******************************************
  ! Calculate RHRS B1 Term Prefactor Averages
  !******************************************

  ! Much like what is done for the prefactors for the RRS HT term, the 
  ! averages of the prefactor for the RHRS B1 term are calculated.  These
  ! were done using equations in the following references for 
  ! measurements of scattered light perpendicular to the incident 
  ! radiation direction:
  !
  ! Quinet et al.  Int. J. Quant. Chem. 106, 599 (2005).
  ! Yang et al.  J. Chem. Phys 97, 3831 (1992).

  ! Calculate every combination of prefactors.

  write(fname,'(a)') 'RHRS_B1Term_Prefactors.out'
  write(fname1,'(a)') 'RHRS_B1Term_Prefactor_Comps.out'
  rhrs_pref(1:nexci,1:nmodes) = zero
  do iex = 1,nexci
    do j = 1,nmodes
      pref(1,1,1) = stpm(iex,1)*tdipder(iex,j,1) ! xxx
      pref(1,1,2) = stpm(iex,1)*tdipder(iex,j,2) ! xxy
      pref(1,1,3) = stpm(iex,1)*tdipder(iex,j,3) ! xxz
      pref(1,2,1) = stpm(iex,4)*tdipder(iex,j,1) ! xyx
      pref(1,2,2) = stpm(iex,4)*tdipder(iex,j,2) ! xyy
      pref(1,2,3) = stpm(iex,4)*tdipder(iex,j,3) ! xyz
      pref(1,3,1) = stpm(iex,5)*tdipder(iex,j,1) ! xzx
      pref(1,3,2) = stpm(iex,5)*tdipder(iex,j,2) ! xzy
      pref(1,3,3) = stpm(iex,5)*tdipder(iex,j,3) ! xzz
      pref(2,1,1) = stpm(iex,4)*tdipder(iex,j,1) ! yxx
      pref(2,1,2) = stpm(iex,4)*tdipder(iex,j,2) ! yxy
      pref(2,1,3) = stpm(iex,4)*tdipder(iex,j,3) ! yxz
      pref(2,2,1) = stpm(iex,2)*tdipder(iex,j,1) ! yyx
      pref(2,2,2) = stpm(iex,2)*tdipder(iex,j,2) ! yyy
      pref(2,2,3) = stpm(iex,2)*tdipder(iex,j,3) ! yyz
      pref(2,3,1) = stpm(iex,6)*tdipder(iex,j,1) ! yzx
      pref(2,3,2) = stpm(iex,6)*tdipder(iex,j,2) ! yzy
      pref(2,3,3) = stpm(iex,6)*tdipder(iex,j,3) ! yzz
      pref(3,1,1) = stpm(iex,5)*tdipder(iex,j,1) ! zxx
      pref(3,1,2) = stpm(iex,5)*tdipder(iex,j,2) ! zxy
      pref(3,1,3) = stpm(iex,5)*tdipder(iex,j,3) ! zxz
      pref(3,2,1) = stpm(iex,6)*tdipder(iex,j,1) ! zyx
      pref(3,2,2) = stpm(iex,6)*tdipder(iex,j,2) ! zyy
      pref(3,2,3) = stpm(iex,6)*tdipder(iex,j,3) ! zyz
      pref(3,3,1) = stpm(iex,3)*tdipder(iex,j,1) ! zzx
      pref(3,3,2) = stpm(iex,3)*tdipder(iex,j,2) ! zzy
      pref(3,3,3) = stpm(iex,3)*tdipder(iex,j,3) ! zzz

      open(81,file=fname1,form='formatted')
      write(81,'(i5,i5,f10.2)') iex, j, freq(j)
!      write(81,'(a5)') 'X'
      write(81,'(a14,a14,a14)') 'x','y','z'
      write(81,'(a4,f14.4,f14.4,f14.4)') 'xx',(pref(1,1,b), b=1,3)
      write(81,'(a4,f14.4,f14.4,f14.4)') 'xy',(pref(1,2,b), b=1,3)
      write(81,'(a4,f14.4,f14.4,f14.4)') 'xz',(pref(1,3,b), b=1,3)
!      write(81,'(a5)') 'Y'
      write(81,'(a14,a14,a14)') 'x','y','z'
      write(81,'(a4,f14.4,f14.4,f14.4)') 'yx',(pref(2,1,b), b=1,3)
      write(81,'(a4,f14.4,f14.4,f14.4)') 'yy',(pref(2,2,b), b=1,3)
      write(81,'(a4,f14.4,f14.4,f14.4)') 'yz',(pref(2,3,b), b=1,3)
!      write(81,'(a5)') 'Z'
      write(81,'(a14,a14,a14)') 'x','y','z'
      write(81,'(a4,f14.4,f14.4,f14.4)') 'zx',(pref(3,1,b), b=1,3)
      write(81,'(a4,f14.4,f14.4,f14.4)') 'zy',(pref(3,2,b), b=1,3)
      write(81,'(a4,f14.4,f14.4,f14.4)') 'zz',(pref(3,3,b), b=1,3)
      b1 = zero 
      do a=1,3
        b1 = b1 + pref(a,a,a)*pref(a,a,a)
      end do
    
      b2 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b2 = b2 + pref(a,a,b)*pref(a,a,b)
        end do
      end do
    
      b3 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b3 = b3 + pref(a,a,a)*pref(a,b,b)
        end do
      end do
    
      b4 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b4 = b4 + pref(b,a,a)*pref(a,a,b)
        end do
      end do
    
      b5 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b5 = b5 + pref(a,a,a)*pref(b,b,a)
        end do
      end do
    
      b6 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b6 = b6 + pref(b,a,a)*pref(b,a,a)
        end do
      end do
    
      b7 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b7 = b7 + pref(a,a,b)*pref(b,c,c)
          end do
        end do
      end do
    
      b8 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b8 = b8 + pref(b,a,a)*pref(b,c,c)
          end do
        end do
      end do
    
      b9 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b9 = b9 + pref(a,a,b)*pref(c,c,b)
          end do
        end do
      end do
    
      b10 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.c==c) cycle
            b10 = b10 + pref(a,b,c)*pref(a,b,c)
          end do
        end do
      end do
    
      b11 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b11 = b11 + pref(a,b,c)*pref(b,a,c)
          end do
        end do
      end do
    
      b12 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
          b12 = b12 + pref(a,b,b)*pref(a,b,b)
        end do
      end do
    
      b13 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
          b13 = b13 + pref(a,a,b)*pref(b,a,a)
        end do
      end do
    
      b14 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b14 = b14 + pref(a,b,b)*pref(a,c,c)
          end do
        end do
      end do
    
      b15 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b15 = b15 + pref(a,a,c)*pref(b,b,c)
          end do
        end do
      end do
    
      bmean1 = b1/seven
      bmean1 = bmean1 + four*b2/thrtfv + two*b3/thrtfv + &
               four*b4/thrtfv + four*b5/thrtfv + b6/thrtfv
      bmean1 = bmean1 + four*b7/hndrfv + b8/hndrfv + &
               four*b9/hndrfv + two*b10/hndrfv + four*b11/hndrfv
    
      bmean2 = b1/thrtfv
      bmean2 = bmean2 + four*b3/hndrfv - four*b5/seventy + &
               eight*b2/hndrfv + three*b12/thrtfv - four*b13/seventy
      bmean2 = bmean2 + b14/thrtfv -four*b15/twhndrten - &
               four*b7/twhndrten + two*b10/thrtfv -four*b11/twhndrten
    
      rhrs_pref(iex,j) = bmean1 + bmean2

      open(82,file=fname,form='formatted')
      write(82,'(i5,i5,f8.2,f6.2,f14.2)') iex, j, &
                               freq(j), delta(iex,j), &
                               rhrs_pref(iex,j)
    end do
    write(82,*)
  end do
  close(81)
  close(82)

End Subroutine RHRSPrefB1

Subroutine RHRSPrefB2(rhrs_pref, osc, stpmder, delta, freq, nexci, &
                      nmodes, tdim, sdim) 

  use constants

  implicit none

  integer :: nexci, nmodes, tdim, sdim, iex, j
  integer :: a, b, c

  character(80) :: fname, fname1

  real(kindr), intent(in) :: osc(nexci,tdim)
  real(kindr), intent(in) :: stpmder(nexci,nmodes,sdim)
  real(kindr), intent(inout) :: rhrs_pref(nexci,nmodes)  
  real(kindr) :: bmean1, bmean2
  real(kindr) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11
  real(kindr) :: b12, b13, b14, b15
  real(kindr) :: pref(3,3,3) 
  real(kindr) :: delta(nexci,nmodes), freq(nmodes)

! =============================================================================
! purpose: Orientational averaging of the B2 term prefactor for RHRS
! 
! input  : rhrs_pref - orientational averages of the B2 term prefactors
!          freq - vibrational frequencies
!          stpmder - derivatives of two-photon transition moment components
!          osc - transition dipole moment
!          delta - dimensionless displacements
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tdim - final dimension of the transition dipole moment derivative
!          array (equals 3)
!          sdim - final dimension of the two-photon transition moment 
!          array (equals 6)
!     
! =============================================================================

  !******************************************
  ! Calculate RHRS B2 Term Prefactor Averages
  !******************************************

  ! Much like what is done for the prefactors for the RRS HT term, the 
  ! averages of the prefactor for the RHRS B1 term are calculated.  These
  ! were done using equations in the following references for 
  ! measurements of scattered light perpendicular to the incident 
  ! radiation direction:
  !
  ! Quinet et al.  Int. J. Quant. Chem. 106, 599 (2005).
  ! Yang et al.  J. Chem. Phys 97, 3831 (1992).

  ! Calculate every combination of prefactors.

  write(fname,'(a)') 'RHRS_B2Term_Prefactors.out'
  write(fname1,'(a)') 'RHRS_B2Term_Prefactor_Comps.out'
  rhrs_pref(1:nexci,1:nmodes) = zero
  do iex = 1,nexci
    do j = 1,nmodes

      pref(1,1,1) = stpmder(iex,j,1)*osc(iex,1) ! xxx
      pref(1,1,2) = stpmder(iex,j,1)*osc(iex,2) ! xxy
      pref(1,1,3) = stpmder(iex,j,1)*osc(iex,3) ! xxz
      pref(1,2,1) = stpmder(iex,j,4)*osc(iex,1) ! xyx
      pref(1,2,2) = stpmder(iex,j,4)*osc(iex,2) ! xyy
      pref(1,2,3) = stpmder(iex,j,4)*osc(iex,3) ! xyz
      pref(1,3,1) = stpmder(iex,j,5)*osc(iex,1) ! xzx
      pref(1,3,2) = stpmder(iex,j,5)*osc(iex,2) ! xzy
      pref(1,3,3) = stpmder(iex,j,5)*osc(iex,3) ! xzz
      pref(2,1,1) = stpmder(iex,j,4)*osc(iex,1) ! yxx
      pref(2,1,2) = stpmder(iex,j,4)*osc(iex,2) ! yxy
      pref(2,1,3) = stpmder(iex,j,4)*osc(iex,3) ! yxz
      pref(2,2,1) = stpmder(iex,j,2)*osc(iex,1) ! yyx
      pref(2,2,2) = stpmder(iex,j,2)*osc(iex,2) ! yyy
      pref(2,2,3) = stpmder(iex,j,2)*osc(iex,3) ! yyz
      pref(2,3,1) = stpmder(iex,j,6)*osc(iex,1) ! yzx
      pref(2,3,2) = stpmder(iex,j,6)*osc(iex,2) ! yzy
      pref(2,3,3) = stpmder(iex,j,6)*osc(iex,3) ! yzz
      pref(3,1,1) = stpmder(iex,j,5)*osc(iex,1) ! zxx
      pref(3,1,2) = stpmder(iex,j,5)*osc(iex,2) ! zxy
      pref(3,1,3) = stpmder(iex,j,5)*osc(iex,3) ! zxz
      pref(3,2,1) = stpmder(iex,j,6)*osc(iex,1) ! zyx
      pref(3,2,2) = stpmder(iex,j,6)*osc(iex,2) ! zyy
      pref(3,2,3) = stpmder(iex,j,6)*osc(iex,3) ! zyz
      pref(3,3,1) = stpmder(iex,j,3)*osc(iex,1) ! zzx
      pref(3,3,2) = stpmder(iex,j,3)*osc(iex,2) ! zzy
      pref(3,3,3) = stpmder(iex,j,3)*osc(iex,3) ! zzz

      open(83,file=fname1,form='formatted')
      write(83,'(i5,i5,f10.2)') iex, j, freq(j)
!      write(83,'(a5)') 'X'
      write(83,'(a14,a14,a14)') 'x','y','z'
      write(83,'(a4,f14.4,f14.4,f14.4)') 'xx',(pref(1,1,b), b=1,3)
      write(83,'(a4,f14.4,f14.4,f14.4)') 'xy',(pref(1,2,b), b=1,3)
      write(83,'(a4,f14.4,f14.4,f14.4)') 'xz',(pref(1,3,b), b=1,3)
!      write(83,'(a5)') 'Y'
      write(83,'(a14,a14,a14)') 'x','y','z'
      write(83,'(a4,f14.4,f14.4,f14.4)') 'yx',(pref(2,1,b), b=1,3)
      write(83,'(a4,f14.4,f14.4,f14.4)') 'yy',(pref(2,2,b), b=1,3)
      write(83,'(a4,f14.4,f14.4,f14.4)') 'yz',(pref(2,3,b), b=1,3)
!      write(83,'(a5)') 'Z'
      write(83,'(a14,a14,a14)') 'x','y','z'
      write(83,'(a4,f14.4,f14.4,f14.4)') 'zx',(pref(3,1,b), b=1,3)
      write(83,'(a4,f14.4,f14.4,f14.4)') 'zy',(pref(3,2,b), b=1,3)
      write(83,'(a4,f14.4,f14.4,f14.4)') 'zz',(pref(3,3,b), b=1,3)
      
      b1 = zero 
      do a=1,3
        b1 = b1 + pref(a,a,a)*pref(a,a,a)
      end do
    
      b2 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b2 = b2 + pref(a,a,b)*pref(a,a,b)
        end do
      end do
    
      b3 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b3 = b3 + pref(a,a,a)*pref(a,b,b)
        end do
      end do
    
      b4 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b4 = b4 + pref(b,a,a)*pref(a,a,b)
        end do
      end do
    
      b5 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b5 = b5 + pref(a,a,a)*pref(b,b,a)
        end do
      end do
    
      b6 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
            b6 = b6 + pref(b,a,a)*pref(b,a,a)
        end do
      end do
    
      b7 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b7 = b7 + pref(a,a,b)*pref(b,c,c)
          end do
        end do
      end do
    
      b8 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b8 = b8 + pref(b,a,a)*pref(b,c,c)
          end do
        end do
      end do
    
      b9 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b9 = b9 + pref(a,a,b)*pref(c,c,b)
          end do
        end do
      end do
    
      b10 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.c==c) cycle
            b10 = b10 + pref(a,b,c)*pref(a,b,c)
          end do
        end do
      end do
    
      b11 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b11 = b11 + pref(a,b,c)*pref(b,a,c)
          end do
        end do
      end do
    
      b12 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
          b12 = b12 + pref(a,b,b)*pref(a,b,b)
        end do
      end do
    
      b13 = zero 
      do a=1,3
        do b=1,3
          if (a==b) cycle
          b13 = b13 + pref(a,a,b)*pref(b,a,a)
        end do
      end do
    
      b14 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b14 = b14 + pref(a,b,b)*pref(a,c,c)
          end do
        end do
      end do
    
      b15 = zero 
      do a=1,3
        do b=1,3
          do c=1,3
            if (a==b.or.a==c.or.b==c) cycle
            b15 = b15 + pref(a,a,c)*pref(b,b,c)
          end do
        end do
      end do
    
      bmean1 = b1/seven
      bmean1 = bmean1 + four*b2/thrtfv + two*b3/thrtfv + &
               four*b4/thrtfv + four*b5/thrtfv + b6/thrtfv
      bmean1 = bmean1 + four*b7/hndrfv + b8/hndrfv + &
               four*b9/hndrfv + two*b10/hndrfv + four*b11/hndrfv
    
      bmean2 = b1/thrtfv
      bmean2 = bmean2 + four*b3/hndrfv - four*b5/seventy + &
               eight*b2/hndrfv + three*b12/thrtfv - four*b13/seventy
      bmean2 = bmean2 + b14/thrtfv -four*b15/twhndrten - &
               four*b7/twhndrten + two*b10/thrtfv -four*b11/twhndrten
    
      rhrs_pref(iex,j) = bmean1 + bmean2

      open(84,file=fname,form='formatted')
      write(84,'(i5,i5,f8.2,f6.2,f14.2)') iex, j, &
                               freq(j), delta(iex,j), &
                               rhrs_pref(iex,j)
    end do
    write(84,*)
  end do
  close(83)
  close(84)

End Subroutine RHRSPrefB2
