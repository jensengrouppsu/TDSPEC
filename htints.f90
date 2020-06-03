Subroutine RRSHTInts(Ls_htb, Ls_htb2, Ls_htb3, freq, osc, &
                    tdipder, freq0, htexcist, ht_pol, &
                    ht_rcrs, ht_rcrs_sum, delta, lfreq, nexci, &
                    nmodes, tdim, iex, j)

  use constants

  implicit none

  integer :: nexci, nmodes, tdim, j, iex, m, i, k

  character(80) :: fname

  real(kindr) :: freq(nmodes), osc(nexci,tdim)
  real(kindr) :: tdipder(nexci,nmodes,tdim), freq0(nexci)
  real(kindr) :: htexcist, lfreq
  real(kindr) :: ht_rcrs(nmodes,nmodes), ht_rcrs_sum(nmodes)
  real(kindr) :: delta(nexci,nmodes)
  real(kindr) :: betaA, polmn

  complex(kindr) :: Ls_htb(nexci,nmodes), Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes), ht_pol(nmodes,nmodes,3,3)

! =============================================================================
! purpose: Evaluate Herzberg-Teller Coupling between modes in RRS
!
! input  : Ls_htb - lineshape from one contribution to the B term
!          Ls_htb2 - lineshape from one contribution to the B term
!          Ls_htb3 - lineshape from one contribution to the B term
!          freq - vibrational frequencies
!          osc - transition dipole moment components
!          tdipder - derivatives of transition dipole moments along each mode
!          freq0 - vertical excitation energy
!          htexcist - excited state of interest
!          ht_pol - polarizability tensor from the HT term
!          ht_rcrs - RRS cross section for a particular mode
!          ht_rcrs_sum - RRS cross section for all modes summed
!          delta - dimensionless displacements
!          lfreq - laser frequency
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tdim - final dimension of the transition dipole moment derivative
!          array (equals 3)
!          iex - number of the current excited state
!          j - index for the current mode 
!
! =============================================================================

  if (abs(lfreq - freq0(iex)) == htexcist) then
    ! Calculate the polarizability related to the HT term.
    ! Keep each normal mode's contribution separate (hence,
    ! the extra dimension to the array).
    do k = 1, nmodes
      ht_pol(j,k,1,1) = ht_pol(j,k,1,1) + &
              (osc(iex,1)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
              tdipder(iex,k,1)*osc(iex,1)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
      ht_pol(j,k,1,2) = ht_pol(j,k,1,2) + &
              (osc(iex,1)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
              tdipder(iex,k,1)*osc(iex,2)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
      ht_pol(j,k,1,3) = ht_pol(j,k,1,3) + &
              (osc(iex,1)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
              tdipder(iex,k,1)*osc(iex,3)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
      ht_pol(j,k,2,1) = ht_pol(j,k,2,1) + &
              (osc(iex,2)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
              tdipder(iex,k,2)*osc(iex,1)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
      ht_pol(j,k,2,2) = ht_pol(j,k,2,2) + &
              (osc(iex,2)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
              tdipder(iex,k,2)*osc(iex,2)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
      ht_pol(j,k,2,3) = ht_pol(j,k,2,3) + &
              (osc(iex,2)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
              tdipder(iex,k,2)*osc(iex,3)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
      ht_pol(j,k,3,1) = ht_pol(j,k,3,1) + &
              (osc(iex,3)*tdipder(iex,k,1)*Ls_htb(iex,k)+ &
              tdipder(iex,k,3)*osc(iex,1)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
      ht_pol(j,k,3,2) = ht_pol(j,k,3,2) + &
              (osc(iex,3)*tdipder(iex,k,2)*Ls_htb(iex,k)+ &
              tdipder(iex,k,3)*osc(iex,2)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
      ht_pol(j,k,3,3) = ht_pol(j,k,3,3) + &
              (osc(iex,3)*tdipder(iex,k,3)*Ls_htb(iex,k)+ &
              tdipder(iex,k,3)*osc(iex,3)*(Ls_htb2(iex,k)+ &
              Ls_htb3(iex,k)))*pol2si
    end do
    ! Perform the averaging of each mode's polarizability, and 
    ! calculate the scattering factor.  Also find the total
    ! contribution to the polarizability from every mode.
    ht_rcrs_sum(j) = zero
    do m=1,nmodes
      betaA = zero
      polmn = zero
      do i=1,3
        do k=1,3
          betaA = betaA + f12* &
            (three*ht_pol(j,m,i,k)*conjg(ht_pol(j,m,i,k)) - &
            ht_pol(j,m,i,i)*conjg(ht_pol(j,m,k,k)))
        end do
      end do

      do i=1,3
        polmn = polmn + ht_pol(j,m,i,i)*conjg(ht_pol(j,m,i,i))
      end do
      polmn = polmn/nine
      ht_rcrs(j,m) = fortyfv*polmn + seven*betaA
      ht_rcrs_sum(j) = ht_rcrs_sum(j) + ht_rcrs(j,m)
    end do
    ! Output the data.
    write(fname,'(a,i3,a,i1,a)') 'RRS_PolBTerm_', &
                int(1.0E+7_kindr/lfreq),'exc_',iex,'state.out'
    open(8,file=fname,form='formatted')
    do k = 1, nmodes
        write(8,'(i5,i5,f8.2,f6.2,e16.4E3,f10.4)') j, k, &
                     freq(k), delta(iex,k), ht_rcrs(j,k), &
                     (ht_rcrs(j,k)/ht_rcrs_sum(j))*100_kindr
    end do
    write(8,*)
    close(8)
  else
    ! Do nothing!
  end if

End Subroutine RRSHTInts

Subroutine B1HTInts(Ls_htb2, Ls_htb3, freq, stpm, tdipder, &
                    freq0, htexcist, ht_hpol, ht_rcrs, &
                    ht_rcrs_sum, delta, lfreq, nexci, nmodes, &
                    tdim, sdim, iex, j)
  
  use constants

  implicit none

  integer :: nexci, nmodes, tdim, sdim, j, iex, m, i, r, k

  character(80) :: fname

  real(kindr) :: freq(nmodes), stpm(nexci,sdim)
  real(kindr) :: tdipder(nexci,nmodes,tdim), freq0(nexci)
  real(kindr) :: htexcist, lfreq
  real(kindr) :: ht_rcrs(nmodes,nmodes), ht_rcrs_sum(nmodes)
  real(kindr) :: delta(nexci,nmodes)
  real(kindr) :: bmean1, bmean2
  real(kindr) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11
  real(kindr) :: b12, b13, b14, b15

  complex(kindr) :: Ls_htb2(nexci,nmodes), Ls_htb3(nexci,nmodes)
  complex(kindr) :: ht_hpol(nmodes,nmodes,3,3,3)  

! =============================================================================
! purpose: Evaluate Herzberg-Teller Coupling between modes in RHRS (B1 Term)
!
! input  : Ls_htb2 - lineshape from one contribution to the B1 term
!          Ls_htb3 - lineshape from one contribution to the B1 term
!          freq - vibrational frequencies
!          stpm - two-photon transition moment components
!          tdipder - derivatives of transition dipole moments along each mode
!          freq0 - vertical excitation energy
!          htexcist - excited state of interest
!          ht_hpol - hyperpolarizability tensor from the HT term
!          ht_rcrs - RHRS cross section for a particular mode
!          ht_rcrs_sum - RHRS cross section for all modes summed
!          delta - dimensionless displacements
!          lfreq - laser frequency
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tdim - final dimension of the transition dipole moment derivative
!          array (equals 3)
!          sdim - final dimension of the two-photon transition moment 
!          array (equals 6)
!          iex - number of the current excited state
!          j - index for the current mode 
!
! =============================================================================

  if (abs(lfreq - freq0(iex)) == htexcist) then
    do k = 1, nmodes

      ht_hpol(j,k,1,1,1) = &
        tdipder(iex,k,1)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,1,1,2) = &
        tdipder(iex,k,1)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,1,1,3) = &
        tdipder(iex,k,1)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si

      ht_hpol(j,k,1,2,1) = &
        tdipder(iex,k,1)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,1,2,2) = &
        tdipder(iex,k,1)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,1,2,3) = &
        tdipder(iex,k,1)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si

      ht_hpol(j,k,1,3,1) = &
        tdipder(iex,k,1)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,1,3,2) = &
        tdipder(iex,k,1)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,1,3,3) = &
        tdipder(iex,k,1)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si

      ht_hpol(j,k,2,1,1) = &
        tdipder(iex,k,2)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,2,1,2) = &
        tdipder(iex,k,2)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,2,1,3) = &
        tdipder(iex,k,2)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si

      ht_hpol(j,k,2,2,1) = &
        tdipder(iex,k,2)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,2,2,2) = &
        tdipder(iex,k,2)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,2,2,3) = &
        tdipder(iex,k,2)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si

      ht_hpol(j,k,2,3,1) = &
        tdipder(iex,k,2)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,2,3,2) = &
        tdipder(iex,k,2)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,2,3,3) = &
        tdipder(iex,k,2)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si

      ht_hpol(j,k,3,1,1) = &
        tdipder(iex,k,3)*stpm(iex,1)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,3,1,2) = &
        tdipder(iex,k,3)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,3,1,3) = &
        tdipder(iex,k,3)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si

      ht_hpol(j,k,3,2,1) = &
        tdipder(iex,k,3)*stpm(iex,4)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,3,2,2) = &
        tdipder(iex,k,3)*stpm(iex,2)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,3,2,3) = &
        tdipder(iex,k,3)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si

      ht_hpol(j,k,3,3,1) = &
        tdipder(iex,k,3)*stpm(iex,5)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,3,3,2) = &
        tdipder(iex,k,3)*stpm(iex,6)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
      ht_hpol(j,k,3,3,3) = &
        tdipder(iex,k,3)*stpm(iex,3)*(Ls_htb2(iex,k)+ &
        Ls_htb3(iex,k))*hpol2si
    end do
    ! Perform the averaging of each mode's hyperpolarizability, and 
    ! calculate the scattering factor.  Also find the total
    ! contribution to the polarizability from every mode.
    ht_rcrs_sum(j) = zero
    do m=1,nmodes
      b1 = zero
      do i=1,3
        b1 = b1 + ht_hpol(j,m,i,i,i)*conjg(ht_hpol(j,m,i,i,i))
      end do

      b2 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b2 = b2 + ht_hpol(j,m,i,i,r)*conjg(ht_hpol(j,m,i,i,r))
        end do
      end do

      b3 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b3 = b3 + ht_hpol(j,m,i,i,i)*conjg(ht_hpol(j,m,i,r,r))
        end do
      end do

      b4 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b4 = b4 + ht_hpol(j,m,r,i,i)*conjg(ht_hpol(j,m,i,i,r))
        end do
      end do

      b5 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b5 = b5 + ht_hpol(j,m,i,i,i)*conjg(ht_hpol(j,m,r,r,i))
        end do
      end do

      b6 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b6 = b6 + ht_hpol(j,m,r,i,i)*conjg(ht_hpol(j,m,r,i,i))
        end do
      end do

      b7=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b7 = b7 + ht_hpol(j,m,i,i,r)*conjg(ht_hpol(j,m,r,k,k))
          end do
        end do
      end do

      b8=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b8 = b8 + ht_hpol(j,m,r,i,i)*conjg(ht_hpol(j,m,r,k,k))
          end do
        end do
      end do

      b9=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b9 = b9 + ht_hpol(j,m,i,i,r)*conjg(ht_hpol(j,m,k,k,r))
          end do
        end do
      end do

      b10=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b10 = b10 + ht_hpol(j,m,i,r,k)*conjg(ht_hpol(j,m,i,r,k))
          end do
        end do
      end do

      b11=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b11 = b11 + ht_hpol(j,m,i,r,k)*conjg(ht_hpol(j,m,r,i,k))
          end do
        end do
      end do

      b12 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
          b12 = b12 + ht_hpol(j,m,i,r,r)*conjg(ht_hpol(j,m,i,r,r))
        end do
      end do

      b13 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
          b13 = b13 + ht_hpol(j,m,i,i,r)*conjg(ht_hpol(j,m,r,i,i))
        end do
      end do

      b14 = zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b14 = b14 + ht_hpol(j,m,i,r,r)*conjg(ht_hpol(j,m,i,k,k))
          end do
        end do
      end do

      b15 = zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b15 = b15 + ht_hpol(j,m,i,i,k)*conjg(ht_hpol(j,m,r,r,k))
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

      ht_rcrs(j,m) = bmean1 + bmean2
      ht_rcrs_sum(j) = ht_rcrs_sum(j) + ht_rcrs(j,m)
    end do

    ! Output the data.
    write(fname,'(a,i3,a,i1,a)') 'RHRS_HpolB1Term_', &
                int(1.0E+7_kindr/lfreq),'exc_',iex,'state.out'
    open(8,file=fname,form='formatted')
    do k = 1, nmodes
        write(8,'(i5,i5,f8.2,f6.2,e14.4E3,f10.4)') j, k, &
                     freq(k), delta(iex,k), ht_rcrs(j,k), &
                     (ht_rcrs(j,k)/ht_rcrs_sum(j))*100_kindr
    end do
    write(8,*)
    close(8)
  else
    ! Do nothing!
  end if

End Subroutine B1HTInts

Subroutine B2HTInts(Ls_htb, freq, stpmder, osc, &
                    freq0, htexcist, ht_hpol, ht_rcrs, &
                    ht_rcrs_sum, delta, lfreq, nexci, nmodes, &
                    tdim, sdim, iex, j)

  use constants

  implicit none

  integer :: nexci, nmodes, tdim, sdim, j, iex, m, i, r, k

  character(80) :: fname

  real(kindr) :: freq(nmodes), stpmder(nexci,nmodes,sdim)
  real(kindr) :: osc(nexci,tdim), freq0(nexci)
  real(kindr) :: htexcist, lfreq
  real(kindr) :: ht_rcrs(nmodes,nmodes), ht_rcrs_sum(nmodes)
  real(kindr) :: delta(nexci,nmodes)
  real(kindr) :: bmean1, bmean2
  real(kindr) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11
  real(kindr) :: b12, b13, b14, b15

  complex(kindr) :: Ls_htb(nexci,nmodes)
  complex(kindr) :: ht_hpol(nmodes,nmodes,3,3,3)

! =============================================================================
! purpose: Evaluate Herzberg-Teller Coupling between modes in RHRS (B2 Term)
!
! input  : Ls_htb - lineshape from one contribution to the B2 term
!          freq - vibrational frequencies
!          stpmder - derivative of two-photon transition moment components
!          osc - transition dipole moments
!          freq0 - vertical excitation energy
!          htexcist - excited state of interest
!          ht_hpol - hyperpolarizability tensor from the HT term
!          ht_rcrs - RHRS cross section for a particular mode
!          ht_rcrs_sum - RHRS cross section for all modes summed
!          delta - dimensionless displacements
!          lfreq - laser frequency
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tdim - final dimension of the transition dipole moment derivative
!          array (equals 3)
!          sdim - final dimension of the two-photon transition moment 
!          array (equals 6)
!          iex - number of the current excited state
!          j - index for the current mode 
!
! =============================================================================

  if (abs(lfreq - freq0(iex)) == htexcist) then
    do k = 1, nmodes

      ht_hpol(j,k,1,1,1) = &
        osc(iex,1)*stpmder(iex,k,1)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,1,1,2) = &
        osc(iex,1)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,1,1,3) = &
        osc(iex,1)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si

      ht_hpol(j,k,1,2,1) = &
        osc(iex,1)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,1,2,2) = &
        osc(iex,1)*stpmder(iex,k,2)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,1,2,3) = &
        osc(iex,1)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si

      ht_hpol(j,k,1,3,1) = &
        osc(iex,1)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,1,3,2) = &
        osc(iex,1)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,1,3,3) = &
        osc(iex,1)*stpmder(iex,k,3)*Ls_htb(iex,k)*hpol2si

      ht_hpol(j,k,2,1,1) = &
        osc(iex,2)*stpmder(iex,k,1)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,2,1,2) = &
        osc(iex,2)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,2,1,3) = &
        osc(iex,2)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si

      ht_hpol(j,k,2,2,1) = &
        osc(iex,2)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,2,2,2) = &
        osc(iex,2)*stpmder(iex,k,2)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,2,2,3) = &
        osc(iex,2)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si

      ht_hpol(j,k,2,3,1) = &
        osc(iex,2)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,2,3,2) = &
        osc(iex,2)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,2,3,3) = &
        osc(iex,2)*stpmder(iex,k,3)*Ls_htb(iex,k)*hpol2si

      ht_hpol(j,k,3,1,1) = &
        osc(iex,3)*stpmder(iex,k,1)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,3,1,2) = &
        osc(iex,3)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,3,1,3) = &
        osc(iex,3)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si

      ht_hpol(j,k,3,2,1) = &
        osc(iex,3)*stpmder(iex,k,4)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,3,2,2) = &
        osc(iex,3)*stpmder(iex,k,2)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,3,2,3) = &
        osc(iex,3)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si

      ht_hpol(j,k,3,3,1) = &
        osc(iex,3)*stpmder(iex,k,5)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,3,3,2) = &
        osc(iex,3)*stpmder(iex,k,6)*Ls_htb(iex,k)*hpol2si
      ht_hpol(j,k,3,3,3) = &
        osc(iex,3)*stpmder(iex,k,3)*Ls_htb(iex,k)*hpol2si
    end do

!    do k=1,nmodes
!      do m=1,3
!        do i=1,3
!          do r=1,3
!            write(*,*) ht_hpol(j,k,m,i,r)
!          end do
!        end do
!      end do
!    end do

    ! Perform the averaging of each mode's hyperpolarizability, and 
    ! calculate the scattering factor.  Also find the total
    ! contribution to the polarizability from every mode.
    ht_rcrs_sum(j) = zero
    do m=1,nmodes
      b1 = zero
      do i=1,3
        b1 = b1 + ht_hpol(j,m,i,i,i)*conjg(ht_hpol(j,m,i,i,i))
      end do

      b2 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b2 = b2 + ht_hpol(j,m,i,i,r)*conjg(ht_hpol(j,m,i,i,r))
        end do
      end do

      b3 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b3 = b3 + ht_hpol(j,m,i,i,i)*conjg(ht_hpol(j,m,i,r,r))
        end do
      end do

      b4 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b4 = b4 + ht_hpol(j,m,r,i,i)*conjg(ht_hpol(j,m,i,i,r))
        end do
      end do

      b5 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b5 = b5 + ht_hpol(j,m,i,i,i)*conjg(ht_hpol(j,m,r,r,i))
        end do
      end do

      b6 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
            b6 = b6 + ht_hpol(j,m,r,i,i)*conjg(ht_hpol(j,m,r,i,i))
        end do
      end do

      b7=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b7 = b7 + ht_hpol(j,m,i,i,r)*conjg(ht_hpol(j,m,r,k,k))
          end do
        end do
      end do

      b8=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b8 = b8 + ht_hpol(j,m,r,i,i)*conjg(ht_hpol(j,m,r,k,k))
          end do
        end do
      end do

      b9=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b9 = b9 + ht_hpol(j,m,i,i,r)*conjg(ht_hpol(j,m,k,k,r))
          end do
        end do
      end do

      b10=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b10 = b10 + ht_hpol(j,m,i,r,k)*conjg(ht_hpol(j,m,i,r,k))
          end do
        end do
      end do

      b11=zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b11 = b11 + ht_hpol(j,m,i,r,k)*conjg(ht_hpol(j,m,r,i,k))
          end do
        end do
      end do

      b12 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
          b12 = b12 + ht_hpol(j,m,i,r,r)*conjg(ht_hpol(j,m,i,r,r))
        end do
      end do

      b13 = zero
      do i=1,3
        do r=1,3
          if (i==r) cycle
          b13 = b13 + ht_hpol(j,m,i,i,r)*conjg(ht_hpol(j,m,r,i,i))
        end do
      end do

      b14 = zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b14 = b14 + ht_hpol(j,m,i,r,r)*conjg(ht_hpol(j,m,i,k,k))
          end do
        end do
      end do

      b15 = zero
      do i=1,3
        do r=1,3
          do k=1,3
            if (i==r.or.i==k.or.r==k) cycle
            b15 = b15 + ht_hpol(j,m,i,i,k)*conjg(ht_hpol(j,m,r,r,k))
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

      ht_rcrs(j,m) = bmean1 + bmean2
      ht_rcrs_sum(j) = ht_rcrs_sum(j) + ht_rcrs(j,m)
    end do
    ! Output the data.
    write(fname,'(a,i3,a,i1,a)') 'RHRS_HpolB2Term_', &
                int(1.0E+7_kindr/lfreq),'exc_',iex,'state.out'
    open(8,file=fname,form='formatted')
    do k = 1, nmodes
        write(8,'(i5,i5,f8.2,f6.2,e14.4E3,f10.4)') j, k, &
                     freq(k), delta(iex,k), ht_rcrs(j,k), &
                     (ht_rcrs(j,k)/ht_rcrs_sum(j))*100_kindr
    end do
    write(8,*)
    close(8)
  else
    ! Do nothing!
  end if

End Subroutine B2HTInts
