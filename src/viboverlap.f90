Subroutine VibOverlap(x, w, freq, delta, voverlap, freq0, nexci, &
                      nmodes, tstep)

  use constants

  implicit none

  integer :: nexci, nmodes
  integer :: j, iex, t, i, tstep, r, k, n

  complex(kindr), intent(inout) :: voverlap(nexci,nmodes,16)
  complex(kindr) :: Ls(nexci), rsum, rsum1, dsum, z
  complex(kindr) :: Ls_htb(nexci,nmodes), Ls_htb2(nexci,nmodes)
  complex(kindr) :: Ls_htb3(nexci,nmodes)
  complex(kindr) :: expsum(nmodes)

  real(kindr) :: freq(nmodes), delta(nexci,nmodes)
  real(kindr) :: freq0(nexci)
  real(kindr) :: x(tstep), w(tstep) 
  real(kindr) :: zeropt 

  complex(kindr) :: zzero=(0d0,0d0) 

  character(80) :: fnameaa, fnamebia, fnamebiia, fnamebiiia
  character(80) :: fnameabi, fnamebibi, fnamebiibi, fnamebiiibi
  character(80) :: fnameabii, fnamebibii, fnamebiibii, fnamebiiibii
  character(80) :: fnameabiii, fnamebibiii, fnamebiibiii, fnamebiiibiii

  z = csqrt((-1,0))

! =============================================================================
! purpose: Evaluate products of vibrational overlap integrals, assuming that 
! the exponential decay (detuning) can be ignored.  This code is identical
! for both RRS and RHRS since the overlap integrals do not distinguish 
! methods.
!
! input  : x - abcissa points for numerical integration
!          w - quadrature weights for numerical integration
!          freq - vibrational frequencies
!          delta - dimensionless displacements
!          freq0 - vertical excitation energy
!          gam - homogeneous broadening
!          gam2 - inhomogeneous broadening
!          nexci - number of excited states
!          nmodes - number of normal modes
!          tstep - number of abcissa points
!
! in/out : pol - RRS polarizability tensor
!
! =============================================================================

! Here is a correpsondence table for the different products of specific
! terms.  The complex conjugate is denoted with a star (*), and 
! different terms are represented by their corresponding overlap
! integrals:
!
! A    => <F|I(t)>
! Bi   => <F|(Ia+1)(t)>
! Bii  => <(Fa-1)|I(t)>
! Biii => <(Fa+1)|I(t)>
!
! 1  voverlap(nexci,nmodes,1)  = A x A*
! 2  voverlap(nexci,nmodes,2)  = Bi x A*
! 3  voverlap(nexci,nmodes,3)  = Bii x A*
! 4  voverlap(nexci,nmodes,4)  = Biii x A*
! 5  voverlap(nexci,nmodes,5)  = A x Bi*
! 6  voverlap(nexci,nmodes,6)  = Bi x Bi*
! 7  voverlap(nexci,nmodes,7)  = Bii x Bi*
! 8  voverlap(nexci,nmodes,8)  = Biii x Bi*
! 9  voverlap(nexci,nmodes,9)  = A x Bii*
! 10 voverlap(nexci,nmodes,10) = Bi x Bii*
! 11 voverlap(nexci,nmodes,11) = Bii x Bii*
! 12 voverlap(nexci,nmodes,12) = Biii x Bii*
! 13 voverlap(nexci,nmodes,13) = A x Biii*
! 14 voverlap(nexci,nmodes,14) = Bi x Biii*
! 15 voverlap(nexci,nmodes,15) = Bii x Biii*
! 16 voverlap(nexci,nmodes,16) = Biii x Biii* 

  ! Find the zero point energies of each normal mode.  Note that
  ! the factor of hbar is not included since it is divided into time
  ! in the code.
  zeropt = zero
  do i=1,nmodes
    zeropt = zeropt + freq(i)/two
  end do

  ! Determine the vibrational overlap for all modes.
  do j = 1, nmodes
    do iex=1,nexci
      Ls(iex) = zzero
      Ls_htb(iex,1:nmodes) = zzero
      Ls_htb2(iex,1:nmodes) = zzero
      Ls_htb3(iex,1:nmodes) = zzero
      do t=1,tstep
        dsum = zzero
        do i=1,nmodes

          ! Calculate the sum term in the integral (redundant 
          ! part that applies both to Raman and absorption 
          ! spectra).

          dsum = dsum + &
                (delta(iex,i)**2/two)*(one-exp(-z*freq(i)*x(t))) 
        end do

        ! rsum contributes to 0-1 overlaps 
        ! rsum1 contributes to 0-2 overlaps        

        rsum = delta(iex,j)/sqrt(two)*(one-exp(-z*freq(j)*x(t)))
        rsum1 = delta(iex,j)**2/(sqrt(eight))* &
               (one-exp(-z*freq(j)*x(t)))**2
        ! Calculate the contribution from other normal modes.  This
        ! is either raised or has nothing performed, if the mode is
        ! the current mode.
        do r=1,nmodes
          if (r /= j) then
            expsum(r) = delta(iex,r)/sqrt(two)* &
                        (one-exp(-z*freq(r)*x(t)))
          else
            expsum(r) = one
          end if
        end do
        ! A term overlap
        Ls(iex) = Ls(iex) + z*w(t)*rsum*exp(x(t)* &
          z*(zeropt+freq0(iex)) - dsum)
        ! B terms
        do n=1,nmodes
          ! Bi term overlap
          if (n == j) then
            Ls_htb(iex,n) = Ls_htb(iex,n) + &
                     z*(one-delta(iex,n)**2+ &
                     delta(iex,n)**2*cos(freq(n)*x(t)))* &
                     w(t)*exp(x(t)* &
                     z*(zeropt+freq0(iex)) - dsum)
          else
             Ls_htb(iex,n) = Ls_htb(iex,n) + z*w(t)*rsum*expsum(n)* &
                      exp(x(t)* &
                      z*(zeropt+freq0(iex)) - dsum)
          end if
          ! Bii term overlap
          if (n == j) then
            Ls_htb2(iex,n) = Ls_htb2(iex,n) + z*w(t)*exp(x(t)* &
                      z*(zeropt+freq0(iex)) - dsum)
          else
            Ls_htb2(iex,n) = zzero
          end if
          ! Biii term overlap 
          if (n == j) then
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      z*w(t)*root2*rsum1*exp(x(t)* &
                      z*(zeropt+freq0(iex)) - dsum)
          else
            Ls_htb3(iex,n) = Ls_htb3(iex,n) + &
                      z*w(t)*expsum(n)*rsum*exp(x(t)* &
                      z*(zeropt+freq0(iex)) - dsum)
          end if
        end do
      end do
    end do
  end do

  ! Store Stuff
  do iex=1,nexci
    do n=1,nmodes
      voverlap(iex,n,1) = Ls(iex)*conjg(Ls(iex))
      voverlap(iex,n,2) = Ls_htb(iex,n)*conjg(Ls(iex))
      voverlap(iex,n,3) = Ls_htb2(iex,n)*conjg(Ls(iex))
      voverlap(iex,n,4) = Ls_htb3(iex,n)*conjg(Ls(iex))
      voverlap(iex,n,5) = Ls(iex)*conjg(Ls_htb(iex,n))
      voverlap(iex,n,6) = Ls_htb(iex,n)*conjg(Ls_htb(iex,n))
      voverlap(iex,n,7) = Ls_htb2(iex,n)*conjg(Ls_htb(iex,n))
      voverlap(iex,n,8) = Ls_htb3(iex,n)*conjg(Ls_htb(iex,n))
      voverlap(iex,n,9) = Ls(iex)*conjg(Ls_htb2(iex,n))
      voverlap(iex,n,10) = Ls_htb(iex,n)*conjg(Ls_htb2(iex,n))
      voverlap(iex,n,11) = Ls_htb2(iex,n)*conjg(Ls_htb2(iex,n))
      voverlap(iex,n,12) = Ls_htb3(iex,n)*conjg(Ls_htb2(iex,n))
      voverlap(iex,n,13) = Ls(iex)*conjg(Ls_htb3(iex,n))
      voverlap(iex,n,14) = Ls_htb(iex,n)*conjg(Ls_htb3(iex,n))
      voverlap(iex,n,15) = Ls_htb2(iex,n)*conjg(Ls_htb3(iex,n))
      voverlap(iex,n,16) = Ls_htb3(iex,n)*conjg(Ls_htb3(iex,n))
    end do
  end do

  ! Output the data.  Only the real part is nonzero.

  ! A x A*
  write(fnameaa,'(a)') 'VibOverlap_AA.out'
  open(7000,file=fnameaa,form='formatted')
  do iex = 1, nexci
    write(7000,*) "Excitation ", iex
    write(7000,*)
    do k = 1, nmodes 
      write(7000,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,1)), aimag(voverlap(iex,k,1))
    end do
  end do
  write(7000,*)
  close(7000)

  ! Bi x A*
  write(fnamebia,'(a)') 'VibOverlap_BiA.out'
  open(7001,file=fnamebia,form='formatted')
  do iex = 1, nexci
    write(7001,*) "Excitation ", iex
    write(7001,*)
    do k = 1, nmodes 
      write(7001,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,2)), aimag(voverlap(iex,k,2))
    end do
  end do
  write(7001,*)
  close(7001)

  ! Bii x A*
  write(fnamebiia,'(a)') 'VibOverlap_BiiA.out'
  open(7002,file=fnamebiia,form='formatted')
  do iex = 1, nexci
    write(7002,*) "Excitation ", iex
    write(7002,*)
    do k = 1, nmodes 
      write(7002,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,3)), aimag(voverlap(iex,k,3))
    end do
  end do
  write(7002,*)
  close(7002)

  ! Biii x A*
  write(fnamebiiia,'(a)') 'VibOverlap_BiiiA.out'
  open(7003,file=fnamebiiia,form='formatted')
  do iex = 1, nexci
    write(7003,*) "Excitation ", iex
    write(7003,*)
    do k = 1, nmodes 
      write(7003,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,4)), aimag(voverlap(iex,k,4))
    end do
  end do
  write(7003,*)
  close(7003)

  ! A x Bi*
  write(fnameabi,'(a)') 'VibOverlap_ABi.out'
  open(7004,file=fnameabi,form='formatted')
  do iex = 1, nexci
    write(7004,*) "Excitation ", iex
    write(7004,*)
    do k = 1, nmodes 
      write(7004,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,5)), aimag(voverlap(iex,k,5))
    end do
  end do
  write(7004,*)
  close(7004)

  ! Bi x Bi*
  write(fnamebibi,'(a)') 'VibOverlap_BiBi.out'
  open(7005,file=fnamebibi,form='formatted')
  do iex = 1, nexci
    write(7005,*) "Excitation ", iex
    write(7005,*)
    do k = 1, nmodes 
      write(7005,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,6)), aimag(voverlap(iex,k,6))
    end do
  end do
  write(7005,*)
  close(7005)

  ! Bii x Bi*
  write(fnamebiibi,'(a)') 'VibOverlap_BiiBi.out'
  open(7006,file=fnamebiibi,form='formatted')
  do iex = 1, nexci
    write(7006,*) "Excitation ", iex
    write(7006,*)
    do k = 1, nmodes 
      write(7006,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,7)), aimag(voverlap(iex,k,7))
    end do
  end do
  write(7006,*)
  close(7006)

  ! Biii x Bi*
  write(fnamebiiibi,'(a)') 'VibOverlap_BiiiBi.out'
  open(7007,file=fnamebiiibi,form='formatted')
  do iex = 1, nexci
    write(7007,*) "Excitation ", iex
    write(7007,*)
    do k = 1, nmodes 
      write(7007,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,8)), aimag(voverlap(iex,k,8))
    end do
  end do
  write(7007,*)
  close(7007)

  ! A x Bii*
  write(fnameabii,'(a)') 'VibOverlap_ABii.out'
  open(7008,file=fnameabii,form='formatted')
  do iex = 1, nexci
    write(7008,*) "Excitation ", iex
    write(7008,*)
    do k = 1, nmodes 
      write(7008,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,9)), aimag(voverlap(iex,k,9))
    end do
  end do
  write(7008,*)
  close(7008)

  ! Bi x Bii*
  write(fnamebibii,'(a)') 'VibOverlap_BiBii.out'
  open(7009,file=fnamebibii,form='formatted')
  do iex = 1, nexci
    write(7009,*) "Excitation ", iex
    write(7009,*)
    do k = 1, nmodes 
      write(7009,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,10)), aimag(voverlap(iex,k,10))
    end do
  end do
  write(7009,*)
  close(7009)

  ! Bii x Bii*
  write(fnamebiibii,'(a)') 'VibOverlap_BiiBii.out'
  open(7010,file=fnamebiibii,form='formatted')
  do iex = 1, nexci
    write(7010,*) "Excitation ", iex
    write(7010,*)
    do k = 1, nmodes 
      write(7010,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,11)), aimag(voverlap(iex,k,11))
    end do
  end do
  write(7010,*)
  close(7010)

  ! Biii x Bii*
  write(fnamebiiibii,'(a)') 'VibOverlap_BiiiBii.out'
  open(7011,file=fnamebiiibii,form='formatted')
  do iex = 1, nexci
    write(7011,*) "Excitation ", iex
    write(7011,*)
    do k = 1, nmodes 
      write(7011,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,12)), aimag(voverlap(iex,k,12))
    end do
  end do
  write(7011,*)
  close(7011)

  ! A x Biii*
  write(fnameabiii,'(a)') 'VibOverlap_ABiii.out'
  open(7012,file=fnameabiii,form='formatted')
  do iex = 1, nexci
    write(7012,*) "Excitation ", iex
    write(7012,*)
    do k = 1, nmodes 
      write(7012,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,13)), aimag(voverlap(iex,k,13))
    end do
  end do
  write(7012,*)
  close(7012)

  ! Bi x Biii*
  write(fnamebibiii,'(a)') 'VibOverlap_BiBiii.out'
  open(7013,file=fnamebibiii,form='formatted')
  do iex = 1, nexci
    write(7013,*) "Excitation ", iex
    write(7013,*)
    do k = 1, nmodes 
      write(7013,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,14)), aimag(voverlap(iex,k,14))
    end do
  end do
  write(7013,*)
  close(7013)

  ! Bii x Biii*
  write(fnamebiibiii,'(a)') 'VibOverlap_BiiBiii.out'
  open(7014,file=fnamebiibiii,form='formatted')
  do iex = 1, nexci
    write(7014,*) "Excitation ", iex
    write(7014,*)
    do k = 1, nmodes 
      write(7014,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,15)), aimag(voverlap(iex,k,15))
    end do
  end do
  write(7014,*)
  close(7014)

  ! Biii x Biii*
  write(fnamebiiibiii,'(a)') 'VibOverlap_BiiiBiii.out'
  open(7015,file=fnamebiiibiii,form='formatted')
  do iex = 1, nexci
    write(7015,*) "Excitation ", iex
    write(7015,*)
    do k = 1, nmodes 
      write(7015,'(f8.2,3X,e12.4E2,3X,e12.4E2)') &
            freq(k), real(voverlap(iex,k,16)), aimag(voverlap(iex,k,16))
    end do
  end do
  write(7015,*)
  close(7015)

end subroutine VibOverlap
