subroutine Lineshape(solvmod, g, x, gam, gam2, kappa, &
                     freqbath, deltabath, tstep, nexci, temper, &
                     labs, ltpa, lfluor, lraman, lhyperraman, lsfg, &
                     lsecondhyperraman) 

  use constants

  implicit none

  character(80) :: frmt

  integer :: solvmod(4), t, tstep, nexci, iex, j
  integer :: sumpts = 50

  complex(kindr),intent(inout) :: g(nexci,tstep)
  complex(kindr) :: z

  real(kindr) :: gr, gi
  real(kindr),intent(in) :: x(tstep)
  real(kindr),intent(in) :: gam(nexci), gam2, kappa(nexci), freqbath(nexci)
  real(kindr),intent(in) :: deltabath(nexci), temper
  real(kindr) :: dsolv(nexci),lsolv(nexci),solvreorg(nexci)
  real(kindr) :: tfactor
  real(kindr) :: nun, numer, denom, sumn
  real(kindr) :: arg

  logical :: labs, ltpa, lfluor
  logical :: lraman, lhyperraman, lsfg
  logical :: lsecondhyperraman

  z = csqrt((-1,0))

! =============================================================================
! purpose: Determine the lifetime function g(t)
!
! input  : solvmod - type of modeling of the solvent
!          x - abcissa points for numerical integration
!          gam - homogeneous broadening
!          gam2 - inhomogeneous broadening
!          kappa - lifetime parameter
!          freqbath - bath frequency
!          deltabath - dimensionless displacement of the bath
!          tstep - number of abcissa points
!          nexci - number of excited states
!          temper - temperature 
!
! output : g - lifetime function
!
! =============================================================================

  ! The functional form of g(t) depends on the model.
  ! 1. For the simple model, where the excited lifetime is Lorentzian 
  !   and the broadening due to solvent is Gaussian, see:
  !
  !   Petrenko, T. and Neese, F.  J. Chem. Phys., 127, 164319 (2007).  
  !   Silverstein, D. W. and Jensen, L.  J. Chem. Theo. Comp., 6, 2845 
  !     (2010). 
  !
  ! 2. For the high-temperature limit of the overdamped Brownian 
  !   oscillator, see:
  !
  !   Mukamel, S.  "Principles of Nonlinear Optical Spectroscopy".
  !     Oxford University Press: New York, 1995.
  !   Myers, A. B.  Ann. Rev. Phys. Chem., 49, 267 (1998).
  !   Li, B.; Johnson, A. E.; Mukamel, S. and Myers, A. B.  
  !     J. Am. Chem. Soc., 116, 11039 (1994).
  !
  ! 3. For the general temperature limit of the overdamped Brownian
  !   oscillator, see:
  !
  !   Mukamel, S.  "Principles of Nonlinear Optical Spectroscopy".
  !     Oxford University Press: New York, 1995.
  !   Li, B.; Johnson, A. E.; Mukamel, S. and Myers, A. B.  
  !     J. Am. Chem. Soc., 116, 11039 (1994).
  !
  !   * Note that there are typos in the JACS reference of 2 and 3 for
  !   the imaginary part of g(t).

  ! Simple model (Lorentzian lifetime and Gaussian solvent 
  ! inhomogeneous broadening).  The Gaussian distribution is only 
  ! included for absorption/fluorescence since it can only be Fourier
  ! transformed in those cases.
  if (solvmod(1)==1) then
    if ((labs.eqv..true.).or.(ltpa.eqv..true.).or.(lfluor.eqv..true.)) then
      do iex=1,nexci
        do t=1,tstep
          gr = gam(iex)*x(t) + f12*(x(t)*gam2)**2
          gi = zero
          g(iex,t) = gr + z*gi 
        end do
      end do
    else if ((lraman.eqv..true.).or.(lhyperraman.eqv..true.).or. &
             (lsfg.eqv..true.).or.(lsecondhyperraman.eqv..true.)) then
      do iex=1,nexci
        do t=1,tstep
          gr = gam(iex)*x(t) 
          gi = zero
          g(iex,t) = gr + z*gi 
        end do
      end do
    end if
  ! High temperature limit of the overdamped Brownian oscillator
  else if (solvmod(2)==1) then
    ! Determine the solvent coupling strength to the electronic 
    ! transition via a Pade approximant.  Also find the solvent
    ! correlation frequency.
    do iex=1,nexci
      dsolv(iex) = gam(iex)* &
        (one+0.85_kindr*kappa(iex)+0.88_kindr*kappa(iex)**2)/ &
        (2.355_kindr + 1.76_kindr*kappa(iex))
      lsolv(iex) = kappa(iex)*dsolv(iex)
      solvreorg(iex) = two*pi*planckbar*speed*100.0_kindr*dsolv(iex)**2/ &
        (two*boltz*temper)
      ! Determine g(t).
      do t=1,tstep
        gr = (dsolv(iex)/lsolv(iex))**2* &
          (exp(-lsolv(iex)*x(t))+lsolv(iex)*x(t)-one)
        gi = (solvreorg(iex)/lsolv(iex))* &
          (one-exp(-lsolv(iex)*x(t))-lsolv(iex)*x(t)) 
        g(iex,t) = gr + z*gi 
      end do
    end do
  ! General temperature limit of the overdamped Brownian oscillator
  else if (solvmod(3)==1) then
    do iex=1,nexci
      ! Determine the solvent coupling strength and solvent correlation
      ! frequency.
      solvreorg(iex) = two*pi*freqbath(iex)*deltabath(iex)**2/two
      tfactor = planckbar*speed*100.0_kindr*freqbath(iex)/(boltz*temper)
      dsolv(iex) = solvreorg(iex)*freqbath(iex)*(two/(exp(tfactor)-one)+one)
      dsolv(iex) = sqrt(dsolv(iex))
      lsolv(iex) = kappa(iex)*dsolv(iex)
      ! Determine g(t).
      do t=1,tstep
        ! Evaluate the sum over the Matsubara frequencies.  This 
        ! appears to be rapidly convergent, even though it is a sum to
        ! infinity.
        sumn = zero
        do j=1,sumpts
          nun = &
            (two*pi*boltz*temper/(planckbar*speed*100.0_kindr))* &
            dfloat(j)
          numer = exp(-nun*x(t)) + nun*x(t) - one
          denom = nun*(nun**2-lsolv(iex)**2)
          sumn = sumn + numer/denom
        end do
        arg = planckbar*lsolv(iex)*100.0_kindr*speed/(two*boltz*temper)
        gr = (solvreorg(iex)/lsolv(iex))*cos(arg)/sin(arg)* &
          (exp(-lsolv(iex)*x(t))+lsolv(iex)*x(t)-one) + &
          four*solvreorg(iex)*lsolv(iex)*boltz*temper*sumn/ &
          (planckbar*speed*100.0_kindr)
        gi = (solvreorg(iex)/lsolv(iex))* &
          (one-exp(-lsolv(iex)*x(t))-lsolv(iex)*x(t)) 
        g(iex,t) = gr + z*gi 
      end do
    end do
  ! High temperature limit of the overdamped Brownian oscillator with 
  ! a Gaussian distribution of excitation energies.
  else if (solvmod(4)==1) then
    ! Determine the solvent coupling strength to the electronic 
    ! transition via a Pade approximant.  Also find the solvent
    ! correlation frequency.
    do iex=1,nexci
      dsolv(iex) = gam(iex)* &
        (one+0.85_kindr*kappa(iex)+0.88_kindr*kappa(iex)**2)/ &
        (2.355_kindr + 1.76_kindr*kappa(iex))
      lsolv(iex) = kappa(iex)*dsolv(iex)
      solvreorg(iex) = two*pi*planckbar*speed*100.0_kindr*dsolv(iex)**2/ &
        (two*boltz*temper)
      ! Determine g(t).
      do t=1,tstep
        gr = (dsolv(iex)/lsolv(iex))**2* &
          (exp(-lsolv(iex)*x(t))+lsolv(iex)*x(t)-one) &
          + f12*(x(t)*gam2)**2
        gi = (solvreorg(iex)/lsolv(iex))* &
          (one-exp(-lsolv(iex)*x(t))-lsolv(iex)*x(t))
        g(iex,t) = gr + z*gi
      end do
    end do
  end if

  ! Print the solvent contribution to the reorganization, if possible.  
  if (solvmod(1)==1) then
    write(*,*) "# Solvent contribution to the reorganization energy"
    write(*,*) "# cannot be calculated in the simple solvent model."
  else if ((solvmod(2)==1).or.(solvmod(3)==1).or.(solvmod(4)==1)) then
    write(*,*) "# Excited State       Solvent Reorganization Energy (cm-1)"
    frmt = "(1X,1A,10X,I2,10X,F8.2)"
    do iex=1,nexci
      write(*,frmt) "#", iex, solvreorg(iex)
    end do
  end if 

end subroutine Lineshape
