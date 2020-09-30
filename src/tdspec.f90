program TDSPEC

!***************
! Running TDSPEC
!***************

! For this version of the code, the input file can contain
!
! Runtype 
! Freq 
! Gamma2 
! Shift 
! Datafile 
! Int 
! HerzTell
! IRCoupling
! Normalize
! Width
! ExcNum
! FullHerzTell
! ATermOff
! FCScreen
! FrqScreen
! Tstop
! Tstep
! NABSpts
! HTInts
! FluorState
! HTPref
! RamanPts
! Frqcut
! StokesSolv
! VOverlap
! RHRSB2
! Polarization
! Angles
! Temp
! NonRes
! SFGFull
! IRUV
! SolvMod
! HPolAvg
! UnitSphere
! SphereRadius
! PrintPol
! PrintAtensor
! PrintGtensor
! PrintAllTensors
! PrintHpol
! End
!
! By default, the options: Runtype, Freq, Gamma2, Shift, Datafile, and
! Width need to be assigned a value in the input file.  The keyword End
! must also be present.  The options are defined as follows
! 
! 1. Runtype - tells the program what type of simulation is being 
! performed (i.e. Abs, Raman, ASRaman, TwoAbs, HyperRaman, ASHyperRaman,
! SFG, Fluor, CD, RVROA, SecondHyperRaman).
!
! 2. Freq - energy of the incident radiation in nanometers
!
! 3. Gamma2 - inhomogeneous broadening parameter in wavenumbers
!
! 4. Shift - solvent shift applied to move each excitation 
! (in Hartrees)
!
! 5. Datafile - file containing the excitation data (i.e. energy, one 
! and two photon transition dipole moments, vibrational frequencies, 
! dimensionless excited state displacements, and dipole derivatives).
!
! 6. Int - the program prints the intensity of each normal mode, without 
! convoluting with a Lorentzian lineshape.  This can be used to obtain 
! stick spectra and compare peak intensities directly.
!
! 7. HerzTell - can be Yes or No.  By default, this is set to No.  If 
! Yes, the program includes Herzberg-Teller terms for the Raman 
! calculations.  This is also used to evaluate the B1 term for 
! hyper-Raman calculations.  OPA and TPA can now be plotted including 
! their Herzberg-Teller contributions also.  Only set this option or 
! RHRSBHT, not both.
!
! 8. IRCoupling - can be Yes or No.  By default, this is set to No.  If
! Yes, the program includes coupling to IR active modes in the 
! resonance hyper-Raman scattering simulation, based on the work of Ann 
! Myers Kelley's group (the paper is referenced in the actual code).
!
! 9. Normalize - can be Yes or No.  By default, this is set to No.  If
! Yes, the program normalizes the resonance Raman or hyper Raman 
! spectrum to the mode with the highest intensity.
!
! 10. Width - sets the width of the Lorentzian line function used for
! convolution.  This must be set to a nonzero, positive value or you 
! will obtain unphysical results.
!
! 11. FullHerzTell - can be Yes or No.  By default, this is set to No.  
! If Yes, the program looks for two-photon moment derivatives along a 
! normal mode.  These are used to evaluate the full B term in 
! hyper-Raman scattering.  Only set this option or HerzTell, not both.
!
! 12. ATermOff - can be Yes or No.  By default, this is set to No.  If
! Yes, the program will not calculate the A term in resonance 
! Raman/hyper-Raman.  Option is primarily for comparing the relative 
! intensity of the Franck-Condon term to the Herzberg-Teller term.
!
! 13. ExcNum - set the "excitation number", which allows you to account
! for overtones and combination bands on Raman/hyper-Raman spectra.  
! The default value is 1, which includes first harmonics only.  This 
! can be set as large as 10, but that number is generally excessive.
! Depending on how the FC factors are screened, the maximum excitation
! number with combination bands is about 4 or 5.  For large molcules
! where 3N-6 is a big number, it is likely that those values are
! excessive.
!
! 14. FCScreen - set the FC overlap integral screening threshold used 
! by the subroutine fcfactor.  Testing is being done on this value, but
! it appears that the default, 1e-4, is sufficient to pick up all 
! important (i.e. perceivable peaks on the Ramam spectrum) combination
! bands.
!
! 15. FrqScreen - set the frequency screening threshold used by the
! subroutine fcfactor.  This is set to 8000 cm-1 by default.
!
! 16. Tstop - stopping point for the numerical integration.  The 
! default value is 0.300 (time units), but this should be set
! larger if a larger range on the absorbance spectrum is required.
!
! 17. Tstep - number of time steps in the range of numerical
! integration (basically the fineness of the integration grid). This
! is set to 500 by default.  Larger numbers are required to remove
! baseline oscillations in the absorption spectrum resulting from
! numerical integration.  It is also required that the value be 
! increased if a larger wavelength range on the absorbance spectrum
! is needed.  Clearly, making this value very large will cause
! the code to run slowly, so trial and error is required.
!
! 18. NABSpts - number of frequencies used to plot the OPA or TPA
! spectrum.  By default, this is set to 10000, where the integration
! is performed starting at the incident frequency and ending at 
! a frequency 20000 cm-1 higher in energy.  Setting this larger is 
! only required if a larger wavelength range is desired.
!
! 19. HTInts - can be Yes or No.  By default this is set to No.  If
! Yes, the program will output the HT term coupling matrices.  These 
! involve the overlap integrals for each contribution to the HT terms.
!
! 20. FluorState - the number of the excited state fluorescence occurs
! from.
!
! 21. HTPref - calculates HT term prefactors.
!
! 22. RamanPts - number of points used in the Lorentzian convolution 
! in Raman and hyper-Raman.  Default value is 2000.
!
! 23. Frqcut - cut off value for the largest frequency used in finding
! the Stokes shift.  Default value is 417 (2*kb*T, T = 300 K).
!
! 24. StokesSolv - Stokes shift due to solvent effects.  Default value 
! is 400.
!
! 25. VOverlap - can be Yes or No.  By default this is set to No.  If
! Yes, the program calculates vibrational overlap integrals in the 
! limit that the exponential detuning term is unity.  This is done in
! addition to the normal calculation. 
!
! 26. RHRSB2 - can be Yes or No.  By default this is set to No.  If 
! Yes, TDSPEC will calculate only the B2 term in RHRS.
!
! 27. Polarization - number from 1 - 7 determing laser polarization in 
! SFG
! 
! 28. Angles - psi and theta - euler angles for SFG
!
! 29. NonRes - can be Yes or No.  By default this is set to No.  If
! Yes, TDSPEC will add in the nonresonant term (static 
! hyperpolarizatbility) to an SFG spectrum
! 
! 30. SFGFull - the syntax for this input is:
! SFGFull freqstart freqend res
! 
! where freqstart is the starting IR frequency, freqend is the ending
! IR frequency, and res is the resolution (stepsize) for stepping along
! the IR frequencies.  With this inserted, the calculation won't
! do a Lorentzian convolution. 
!
! 31. Temp - temperature (in Kelvin) used for simulation.  If this is
! not defined by the user, the calculation defaults to use 300 K.  This
! is useful solely for solvent models II and III (described below).
!
! 32. SolvMod - defines the "solvent" model used.  The following 
! options define how lifetime/solvent effects are incorporated:
! 
!   I. Simple (lifetime is Lorentzian and solvent is Gaussian,
!     giving a Voigt lineshape overall).
!   II. HighTemp (high temperature limit of the overdamped Brownian
!     oscillator model).
!   III. GenTemp (general temperature limit of the overdamped Brownian
!     oscillator model).  This is an expert option.
!
! 33. HPolAvg - defines whether the user wants the hyperpolarizability
! orientationally averaged (see hpolaverage.f90 for details).
!
! 34. UnitSphere - defines whether the user wants plots of the mode
! dependent hyperpolarizability plotted using a unit sphere 
! representation.  
!
! 35. SphereRadius - defines the radius of the sphere for making images
! in the unit sphere representation.  This should be be larger than the 
! longest axis of the molecule.  The default value is 1.0 Angstroms.
!
! 36. PrintPol - prints the polarizability to file _TDSPEC.pol when
! calculating the RRS spectrum. This file may be collected by the
! CHEM package.
!
! 37. PrintAtensor - prints the A-tensor (and C-tensor) to file
! _TDSPEC.pol (see above). Requires the transition quadrupole
! for each excited state being used in the calculation (keyword
! Tquad in Datafile).
! Within the Herzberg-Teller approximation, requires the transition
! quadrupole gradients. These are assumed to the the sixth-eleventh
! numbers in the datafile.
!
! 38. PrintGtensor - prints the dipole-magnetic dipole (G) tensor
! to file _TDSPEC.pol (see # 36.) Also, if present along with
! the PrintAtensor key, will print the quadrupole-magnetic dipole
! (D) tensor. Requres magnetic dipole transitions for each excited
! state being used (keyword mtdip in Datafile).
! Within the Herzberg-Teller approximation, requires the transition
! magnetic dipoles. These are assumed to tbe the sixth-ninth (or
! the twelveth-fourteenth if the PrintAtensor or PrintAllTensors
! key(s) are included) numbers in the datafile.
!
! 39. PrintAllTensors - prints the polarizability, A, G, C and D-tensors.
! See # 36-39 above. Requires Tquad and mtdip in Datafile.
!
! 40. PrintHpol - prints the hyperpolarizability to file _TDSPEC.hpol
! when calculating the RHRS spectrum. This file may be collected by the
! CHEM package.
!
! I have observed segmentation faults for small molecules when
! FCScreen and FrqScreen are poorly chosen, possibly resulting from
! large numbers of combination bands.  These issues are definitely 
! associated with the subroutine FCfactor, which calculates FC overlaps
! using the dimensionless deltas.  Secondly, if you want overtones
! only, this can be accomplished using very strict settings for those
! screening values (i.e. FCScreen 1.0 or FrqScreen 0.0).

  use constants

  Implicit none

  integer, parameter :: iudata    = 30
  integer :: openstat,readstat,stat
  integer :: nexci
  integer, allocatable :: excnums(:,:), excnumstmp(:,:)
  integer :: nmodes
  integer :: numcb = 0 
  integer :: fluorstate = 0
  integer :: nlines = 0
  integer :: freqend, widthrng
  integer :: iend
 
  character(80) :: indata
  character(80) :: string
  character(80) :: stmp,runtype
  character(80) :: herztell
  character(80) :: ircoupling
  character(80) :: normalize
  character(80) :: fullherztell
  character(80) :: atermoff
  character(80) :: htints
  character(80) :: fname
  character(80) :: htpref
  character(80) :: vboverlap
  character(80) :: rhrsb2
  character(80) :: solvation
  character(80) :: charfreq, charpol, charpsi, chartheta, charfcht
  character(80) :: nonRes

  real(kindr), allocatable :: freq0(:)
  real(kindr), allocatable :: osc(:,:)
  real(kindr), allocatable :: mtdip(:,:)
  real(kindr), allocatable :: tquad(:,:)
  real(kindr), allocatable :: stpm(:,:)
  real(kindr), allocatable :: ttpm(:,:)
  real(kindr), allocatable :: tdipder(:,:,:)
  real(kindr), allocatable :: mtdipder(:,:,:) ! Transition mag dipole derivatives
  real(kindr), allocatable :: tquadder(:,:,:) ! Transition quadrupole derivatives
  real(kindr), allocatable :: gdipder(:,:,:)
  real(kindr), allocatable :: stpmder(:,:,:)
  real(kindr), allocatable :: freq(:),delta(:,:), freq_fsfg(:)
  real(kindr), allocatable :: rcrs(:), rcrs_st(:), rcrs_fsfg(:)
  real(kindr), allocatable :: ht_rcrs(:,:)
  real(kindr), allocatable :: ht_rcrs_sum(:)
  real(kindr), allocatable :: rcrsot(:,:),rcrscb(:)
  real(kindr), allocatable :: boltzpop(:), boltzpopot(:,:)
  real(kindr), allocatable :: boltzpopcb(:), boltzpopas(:) 
  real(kindr), allocatable :: gam(:)
  real(kindr), allocatable :: norm_crs(:)
  real(kindr), allocatable :: norm_frq(:)
  real(kindr), allocatable :: factorial(:)
  real(kindr), allocatable :: rootn(:)
  real(kindr), allocatable :: num(:)
  real(kindr) :: gam2 = 0.0_kindr
  real(kindr) :: sshift = 0.0_kindr
  real(kindr) :: fc_screen = 0.0001_kindr
  real(kindr) :: freq_screen = 8.0E+3_kindr
  real(kindr), allocatable :: freqcb(:), freqcbtmp(:)
  real(kindr), allocatable :: freq_foc(:)
  real(kindr), allocatable :: rcrs_foc(:)
  real(kindr), allocatable :: htlfreqdiff(:)
  real(kindr), allocatable :: rrs_pref(:,:)
  real(kindr), allocatable :: rhrs_pref(:,:)
  real(kindr) :: maxfrq, minfrq
  real(kindr) :: freqcut = two*kbT
  real(kindr) :: stokessolv = 400.0_kindr
  real(kindr) :: stokes
  real(kindr) :: temper = 3.0E+2_kindr
  real(kindr), allocatable :: kappa(:)
  real(kindr), allocatable :: freqbath(:)
  real(kindr), allocatable :: deltabath(:)
  real(kindr), allocatable :: vibreorg(:)
  real(kindr) :: sphereradius = one
  real(kindr) :: freqstrt, res

  complex(kindr), allocatable :: pol(:,:,:)
  complex(kindr), allocatable :: polg(:,:,:)
  complex(kindr), allocatable :: polgs(:,:,:)
  complex(kindr), allocatable :: pola(:,:,:,:)
  complex(kindr), allocatable :: polas(:,:,:,:)
  complex(kindr), allocatable :: aspol(:,:,:)
  complex(kindr), allocatable :: ht_pol(:,:,:,:)
  complex(kindr), allocatable :: polot(:,:,:,:)
  complex(kindr), allocatable :: polcb(:,:,:)
  complex(kindr), allocatable :: hpol(:,:,:,:)
  complex(kindr), allocatable :: ashpol(:,:,:,:)
  complex(kindr), allocatable :: ht_hpol(:,:,:,:,:)
  complex(kindr), allocatable :: hpolot(:,:,:,:,:)
  complex(kindr), allocatable :: hpolcb(:,:,:,:)
  complex(kindr), allocatable :: sechpol(:,:,:,:,:)
  complex(kindr), allocatable :: Ls_htb(:,:), Ls_htb2(:,:), Ls_htb3(:,:)
  complex(kindr), allocatable :: Ls_htir(:,:)
  complex(kindr), allocatable :: expsum(:)
  complex(kindr), allocatable :: voverlap(:,:,:)
  complex(kindr), allocatable :: sumr(:), sumr1(:), sumr2(:)
  complex(kindr), allocatable :: g(:,:)
  complex(kindr), allocatable :: nhpol(:,:,:,:) 
  complex(kindr), allocatable :: fhpol(:,:,:,:)
  complex(kindr), allocatable :: polC(:,:,:)
  complex(kindr), allocatable :: polD(:,:,:)
  complex(kindr), allocatable :: polDs(:,:,:)
   
  complex(kindr), allocatable :: shpol(:,:,:)
  real(kindr) :: tstart
  real(kindr) :: tstop = 0.300_kindr
  real(kindr), allocatable :: x(:),w(:)
  real(kindr) :: lfreq = zero
  real(kindr) :: conv_st
  real(kindr) :: width 
  real(kindr) :: rmax
  real :: time1, time2
  real(kindr) :: psi = 0
  real(kindr) :: theta = 0
  
  complex(kindr) :: z

  integer :: i,j,k,n,iex,r
  integer :: tstep = 500
  integer :: excnum = 1
  integer :: npoints = 2000
  integer :: nabspts = 10000
  integer :: orientation = 1
  integer :: cbguess = 9000
  integer :: solvmod(4)
  integer :: intfreq, inttheta, intpsi

  logical :: labs = .false.
  logical :: ltpa = .false.
  logical :: lraman = .false.
  logical :: lasraman = .false.
  logical :: lhyperraman = .false.
  logical :: lashyperraman = .false.
  logical :: lsfg = .false.
  logical :: lfluor = .false.
  logical :: lcd = .false.
  logical :: lrvroa = .false.
  logical :: lsecondhyperraman = .false.
  logical :: lnormalize = .false.
  logical :: lintensity = .false.
  logical :: lherztell = .false.
  logical :: lircoupling = .false.
  logical :: lfullherztell = .false.
  logical :: latermoff = .false.
  logical :: lhtints = .false.
  logical :: lhtpref = .false.
  logical :: lvoverlap = .false.
  logical :: lrhrsb2 = .false.
  logical :: lhpolavg = .false.
  logical :: lunitsphere = .false.
  logical :: lnonRes = .false.
  logical :: lSFGFull = .false.
  logical :: IRUV = .false.
  logical :: lerror = .false.
  logical :: lruntype = .false.
  logical :: lprintPol = .false.
  logical :: lprintAtensor = .false.
  logical :: lprintGtensor = .false.
  ! Zhongwei:
  logical :: lprintHpol = .false.
!----------------------------------------

  z = csqrt((-1,0))
  ! Initialize the solvent model used.
  solvmod = 0
  solvmod(1) = 1

  ! Defines factorials for overtones and combination bands.

  allocate(factorial(0:15))

  factorial(0) = 1.0_kindr
  factorial(1) = 1.0_kindr
  factorial(2) = 2.0_kindr
  factorial(3) = 6.0_kindr
  factorial(4) = 24.0_kindr
  factorial(5) = 120.0_kindr
  factorial(6) = 720.0_kindr
  factorial(7) = 5040.0_kindr
  factorial(8) = 40320.0_kindr
  factorial(9) = 362880.0_kindr
  factorial(10) = 3628800.0_kindr
  factorial(11) = 39916800.0_kindr
  factorial(12) = 479001600.0_kindr
  factorial(13) = 6227020800.0_kindr
  factorial(14) = 87178291200.0_kindr
  factorial(15) = 1307674368000.0_kindr

  ! Allocate normalization factors due to raising and lowering
  ! operators, used for overtones, combination bands, and
  ! Herzberg-Teller terms.

  allocate(rootn(0:11))

  rootn(0) = sqrt(0.0_kindr)
  rootn(1) = sqrt(1.0_kindr)
  rootn(2) = sqrt(2.0_kindr)
  rootn(3) = sqrt(3.0_kindr)
  rootn(4) = sqrt(4.0_kindr)
  rootn(5) = sqrt(5.0_kindr)
  rootn(6) = sqrt(6.0_kindr)
  rootn(7) = sqrt(7.0_kindr)
  rootn(8) = sqrt(8.0_kindr)
  rootn(9) = sqrt(9.0_kindr)
  rootn(10) = sqrt(10.0_kindr)
  rootn(11) = sqrt(11.0_kindr)

  ! Allocate numbers for HT term overtones and combination bands
  ! so that they may be used in double precision instead of
  ! integer values from the loop.

  allocate(num(0:10))

  num(0) = 0.0_kindr
  num(1) = 1.0_kindr
  num(2) = 2.0_kindr
  num(3) = 3.0_kindr
  num(4) = 4.0_kindr
  num(5) = 5.0_kindr
  num(6) = 6.0_kindr
  num(7) = 7.0_kindr
  num(8) = 8.0_kindr
  num(9) = 9.0_kindr
  num(10) = 10.0_kindr


  !****************************
  ! Define the Calculation Type
  !****************************

  do  
    read(*,'(a80)',iostat=stat) string
    if (stat /= 0) then
      call Error(1, runtype, indata)
    end if
    if (string(1:7)=='Runtype') then
      read(string,*,iostat=readstat) stmp, runtype
      if (runtype=='Abs') then
        labs = .true.
      else if (runtype=='Raman') then
        lraman = .true.
      else if (runtype=='ASRaman') then
        lasraman = .true.
      else if (runtype=='TwoAbs') then
        ltpa = .true.
      else if (runtype=='HyperRaman') then
        lhyperraman = .true.
      else if (runtype=='ASHyperRaman') then
        lashyperraman = .true.
      else if (runtype=='SFG') then
        lsfg = .true.
      else if (runtype=='Fluor') then
        lfluor = .true.
      else if (runtype=='CD') then
        lcd = .true.
      else if (runtype=='RVROA') then
        lrvroa = .true.
      else if (runtype=='SecondHyperRaman') then
        lsecondhyperraman = .true.
      else
        call Error(2, runtype, indata)
      end if
      lruntype = .true.
    else if (string(1:4)=='Freq') then
      read(string,*,iostat=readstat) stmp, lfreq
      ! Divide the wavelength by 2 so that the incident frequency does not
      ! need to be multiplied by 2 in the code later for hyper-Raman. 
      if (lhyperraman.or.lashyperraman) lfreq=f12*lfreq
      ! Divide the wavelength by 3 so the incident frequency does not need
      ! to be multiplied by 3 in the code later for second hyper-Raman
      if (lsecondhyperraman) lfreq = f13*lfreq
      ! The incident frequency is given in units of nm.  The units are then
      ! converted to wavenumbers using the conversion involving 
      ! 1.0E+7 nm cm^{-1} make sense.
      lfreq = 1.0E+7_kindr/lfreq
    else if (string(1:5)=='Shift') then
      read(string,*,iostat=readstat) stmp, sshift
      ! Read in the solvent shift (if necessary).
    else if (string(1:7)=='Gamma2') then
      ! Reads in the value of gamma 2, for accounting for inhomogeneous
      ! broadening effects.
      read(string,*,iostat=readstat) stmp, gam2
    else if (string(1:8)=='Datafile') then
      ! Tells the program which file to look through for information such as
      ! the excitation energies, transition dipoles (one- and two-photon),
      ! vibrational frequencies, and dimensionless deltas.
      read(string,*,iostat=readstat) stmp, indata
    else if (string(1:3)=='Int') then
      ! Reads whether or not you want spectral lines broadened with a 
      ! Lorentzian line shape function.
      lintensity = .true.
    else if (string(1:8)=='HerzTell') then
      ! Determines whether Herzberg-Teller terms (normal coordinate 
      ! dependence of the transition dipole moment) are included in 
      ! the calculation.  This can be used for
      !
      ! 1. Resonance Raman to obtain Albrecht's B term in the 
      ! time-dependent formalism.
      ! 2. Resonance hyper Raman to obtain part of Chung and Ziegler's 
      ! B term in the time-dependent formalism.  This keyword only 
      ! includes the derivative of the transition dipole moment between
      ! the initial and final electronic states).
      read(string,*,iostat=readstat) stmp, herztell
      if (herztell=='Yes') then
        lherztell = .true.
      else
        lherztell = .false.
      end if
    else if (string(1:10)=='IRCoupling') then
      ! Determines whether coupling to IR active normal modes based on
      ! A. M. Kelley's work should be applied for hyper Raman.
      read(string,*,iostat=readstat) stmp, ircoupling
      if (ircoupling=='Yes') then
        lircoupling = .true.
      else
        lircoupling = .false.
      end if
    else if (string(1:9)=='Normalize') then
      ! Determines whether or not you normalize the spectrum to the 
      ! maximum intensity normal mode.
      read(string,*,iostat=readstat) stmp, normalize
      if (normalize=='Yes') then
        lnormalize = .true.
      else
        lnormalize = .false.
      end if
    else if (string(1:5)=='Width') then
      ! Define the line width for the Lorentzian function that is 
      ! convoluted with the Raman or hyper Raman data.
      read(string,*,iostat=readstat) stmp, width
    else if (string(1:12)=='FullHerzTell') then
      ! Determines if the full B term is to be calculated for 
      ! hyper-Raman.  This is important because the input requires
      ! more data for a full B term in hyper-Raman scattering.
      read(string,*,iostat=readstat) stmp, fullherztell
      if (fullherztell=='Yes') then
        lfullherztell = .true.
      else
        lfullherztell = .false.
      end if
    else if (string(1:8)=='ATermOff') then
      ! Determines if the user wants the A term (Franck-Condon) to be 
      ! calculated.
      read(string,*,iostat=readstat) stmp, atermoff
      if (atermoff=='Yes') then
        latermoff = .true.
      else
        latermoff = .false.
      end if
    else if (string(1:6)=='ExcNum') then
      ! Determine the excitation number (for overtones and combination
      ! bands).  By default this is set to 1 (fundamentals only).
      read(string,*,iostat=readstat) stmp, excnum
    else if (string(1:8)=='FCScreen') then
      ! Determine the FC factor screening (for combination bands).
      read(string,*,iostat=readstat) stmp, fc_screen
    else if (string(1:10)=='FrqScreen') then
      ! Determine the frequency screening (for combination bands).
      read(string,*,iostat=readstat) stmp, freq_screen
    else if (string(1:5)=='Tstop') then
      ! Determine the end time value for numerical integration.
      read(string,*,iostat=readstat) stmp, tstop
    else if (string(1:5)=='Tstep') then
      ! Determine the number of steps for numerical integration.
      read(string,*,iostat=readstat) stmp, tstep
    else if (string(1:7)=='NABSpts') then
      ! Determine the number of points for plotting the OPA or TPA 
      ! spectrum.
      read(string,*,iostat=readstat) stmp, nabspts
    else if (string(1:6)=='HTInts') then
      ! Determines if the user wants the B term coupling matrix/matrices 
      ! printed
      read(string,*,iostat=readstat) stmp, htints
      if (htints=='Yes') then
        lhtints = .true.
      else
        lhtints = .false.
      end if
    else if (string(1:10)=='FluorState') then
      ! Determines which excited state to monitor for fluorescence. 
      read(string,*,iostat=readstat) stmp, fluorstate
    else if (string(1:6)=='HTPref') then
      ! Determines if the user wants Herzberg-Teller term prefactor
      ! averages
      read(string,*,iostat=readstat) stmp, htpref
      if (htpref=='Yes') then
        lhtpref = .true.
      else
        lhtpref = .false.
      end if
    else if (string(1:8)=='RamanPts') then
      ! Number of points used in the Lorentzian convolution. 
      read(string,*,iostat=readstat) stmp, npoints
    else if (string(1:6)=='Frqcut') then
      ! Frequency cutoff for low frequency normal modes.  This is
      ! used to find the Stokes shift. 
      read(string,*,iostat=readstat) stmp, freqcut
    else if (string(1:10)=='StokesSolv') then
      ! Stokes shift due to solvent effects.
      read(string,*,iostat=readstat) stmp, stokessolv
    else if (string(1:8)=='VOverlap') then
      ! Does the user want vibrational overlap integral products?
      read(string,*,iostat=readstat) stmp, vboverlap
      if (vboverlap=='Yes') then
        lvoverlap = .true.
      else
        lvoverlap = .false.
      end if
    else if (string(1:6)=='RHRSB2') then
      ! Does the user want B2 terms only in RHRS?
      read(string,*,iostat=readstat) stmp, rhrsb2
      if (rhrsb2=='Yes') then
        lrhrsb2 = .true.
      else
        lrhrsb2 = .false.
      end if
    else if (string(1:6)=='Angles') then
      ! Define the Euler angles for DR-SFG.
      read(string,*,iostat=readstat) stmp, psi, theta
      ! convert to radians
      theta = theta * pi / 180.0_kindr
      psi = psi * pi / 180.0_kindr
    else if (string(1:12)=='Polarization') then
      ! Define the polarization of the incident and scattered photons
      ! in DR-SFG.
      read(string,*,iostat=readstat) stmp,orientation
    else if (string(1:4)=='Temp') then
      ! Define a value for the temperature other than the default.
      read(string,*,iostat=readstat) stmp,temper
      ! Does the user want the non-resonant term used?
    else if (string(1:6)=='NonRes') then
       read(string,*,iostat=readstat) stmp, NonRes
       if (NonRes == 'Yes') then
            lnonRes=.true.
       end if
      ! Do a full DR-SFG calculation without the need for Lorentzian convolution
    else if (string(1:7)=='SFGFull') then
      read(string,*,iostat=readstat) stmp, freqstrt, freqend, res
      lSFGFull = .true.
    else if (string(1:4)=='IRUV') then 
      IRUV = .true.
    else if (string(1:7)=='SolvMod') then
      ! Define the solvation model.
      read(string,*,iostat=readstat) stmp,solvation
      if (solvation=='Simple') then
        solvmod(1) = 1
      else if (solvation=='HighTemp') then
        solvmod(2) = 1
        solvmod(1) = 0
      else if (solvation=='GenTemp') then
        solvmod(3) = 1
        solvmod(1) = 0
      ! Combines the Brownian oscillator with a Gaussian distribution
      ! of excitation energies for absorption spectra.
      else if (solvation=='HighTempGauss') then
        solvmod(4) = 1
        solvmod(1) = 0
      end if
    else if (string(1:7)=='HPolAvg') then
      ! Determine if the user wants Cartesian-orientational averages
      ! of the hyperpolarizability.
      lhpolavg = .true.
    else if (string(1:10)=='UnitSphere') then
      ! Determine if the user wants plots for the hyperpolarizability
      ! in the unit sphere representation.
      lunitsphere = .true.
    else if (string(1:12)=='SphereRadius') then
      ! Defines the radius of the sphere used to make images of the 
      ! unit sphere representation.
      read(string,*,iostat=readstat) stmp,sphereradius
    else if (string(1:5)=='Print') then
      if (string(6:12)=='Atensor') then
        lprintAtensor = .true.
        lprintPol = .true.
      else if (string(6:8)=='Pol') then
        lprintPol = .true.
      ! -------------------------------
      ! Zhongwei:
      else if (string(6:9)=='Hpol') then
        lprintHpol = .true.
      ! -------------------------------
      else if (string(6:12)=='Gtensor') then
        lprintGtensor = .true.
        lprintPol = .true.
      else if (string(6:15)=='AllTensors') then
        lprintPol = .true.
        lprintAtensor = .true.
        lprintGtensor = .true.
      end if
    else if (string(1:3)=='End') then
      exit
    end if
  end do

  ! Check that the necessary parts of the input are present.  If not,
  ! exit with an error.
  if (lruntype.eqv..false.) then
    lerror = .true.
  else if (lfreq==zero) then
    lerror = .true.
  else
    indata = adjustl(indata)
    indata = trim(indata)
    if (len(indata)==0) then
      lerror = .true.
    end if
  end if  

  if (lerror) then
    call Error(1, runtype, indata)
  end if

  ! Output header and calculation information.
  call OutputHeader()
  call OutputCalcInfo(labs, lraman, lasraman, ltpa, lhyperraman, &
             lashyperraman, lsfg, lfluor, lcd, lrvroa, lsecondhyperraman, &
             lfreq, sshift, solvmod, excnum, fc_screen, &
             freq_screen, theta, psi, orientation, lintensity)

! ----------------------------
! open the molecular data file
! ----------------------------
  open(unit=iudata,file=indata,status='old',iostat=openstat)
  if (lruntype) then
    if (openstat /= 0) then
      call Error(3, runtype, indata)
    end if
  end if

! -------------
! read the data
! -------------

! Count the number of lines in the file

  do
    read(iudata,*,iostat=stat)
    nlines = nlines + 1
    if (stat /= 0) then
      nlines = nlines - 1 ! Don't want to count EOF as a line
      exit
    end if
  end do

  rewind(iudata)

! Read the static hyperpolarizabilities for the SFG non-resonant term 

  read(iudata,'(a80)') string
  if (string(1:4) == 'SHyp') then
    allocate (shpol(3,3,3))
    do i=1,3
      read(iudata, *) shpol(i,1,1), shpol(i,1,2), shpol(i,1,3), shpol(i,2,1),&
        shpol(i,2,2), shpol(i,2,3), shpol(i,3,1), shpol(i,3,2), shpol(i,3,3)
    end do
  read(iudata,'(a80)') string
  end if
! Now actually do the reading. 

! iex in this code is the indice counting the number of the excited
! states being examined.

  iex = 0
  if (string(1:5) == 'NEXCI') then
    read(string,*,iostat=readstat) stmp, nexci
  else
    ! If the top line isn't NEXCI, search through lines until
    ! that keyword is found.
    do i = 1,nlines
      read(iudata,'(a80)',iostat=stat) string
      if (stat /= 0) then
        call Error(4, runtype, indata)
      end if
      if (string(1:5) == 'NEXCI') then
        read(string,*,iostat=readstat) stmp, nexci
        exit
      end if
    end do
  end if

  allocate(freq0(nexci))
  allocate(osc(nexci,3))
  allocate(stpm(nexci,6))
  allocate(ttpm(nexci,10))
  allocate(gam(nexci))
  allocate(kappa(nexci))
  allocate(freqbath(nexci))
  allocate(deltabath(nexci))
  allocate(mtdip(nexci,3))
  allocate(tquad(nexci,6))

  do i=1,nexci
    do
      read(iudata,'(a80)') string
      ! Read in the excitation number

      if (string(1:10)=="Excitation") then
        read(string,*,iostat=readstat) stmp,iex
        write(*,'(A36,I2)') '# start reading info excited state:',iex
      end if

      ! Read in the energy for the excitation, accounting for solvent
      ! effects by including an optional solvent shift.  This shift 
      ! is necessary for positioning the peak from theory in the same 
      ! place as experiment.  Also note that the program assumes a 
      ! redshift due to the solvent if the shift is positive, but you 
      ! can do blueshifts using negative shifts.

      if (string(1:6)=="Energy") then
        read(string,*,iostat=readstat) stmp,freq0(i)

        ! Account for the "solvent shift" and convert the units from
        ! Hartrees to wavenumbers.  This allows you to translate the 
        ! excitation to whichever energy you want it to be at.  
        ! Inserting a negative solvent shift will blueshift the 
        ! excitation, a positive one will redshift the excitation.

        freq0(i) = (freq0(i)-sshift)*hartree2cm
      end if

      ! Read the value of homogeneous linewidth (due to lifetime of the
      ! excited state).

      if (string(1:5)=="Gamma") then
        read(string,*,iostat=readstat) stmp,gam(i)
        if ((solvmod(2)==1).or.(solvmod(4)==1)) then
          ! For the high temperature limit, convert to the full width
          ! at half maximum.
          gam(i) = two*gam(i)
        end if
      end if

      ! Define the lineshape parameter for solvent models II and III.

      if (string(1:5)=="Kappa") then
        read(string,*,iostat=readstat) stmp,kappa(i)
      end if

      ! Define the bath frequency for solvent model III.

      if (string(1:9)=="OmegaBath") then
        read(string,*,iostat=readstat) stmp,freqbath(i)
      end if

      ! Define the bath dimensionless displacement for solvent model III.

      if (string(1:9)=="DeltaBath") then
        read(string,*,iostat=readstat) stmp,deltabath(i)
      end if

      ! Read the transition dipole components for a vertical 
      ! transition.  These have atomic units from NWChem.

      if (string(1:4)=="Tdip") then
        read(string,*,iostat=readstat) stmp,osc(i,1:3)
      end if

      ! Read the two photon transition moments for a vertical
      ! transition.  These have atomic units from Dalton.

      if (string(1:4)=="Stpm") then
        read(string,*,iostat=readstat) stmp,stpm(i,1:6)
      end if

      ! Read the three photon transition moments for a vertical
      ! transition.  These have atomic units from Dalton.

      if (string(1:4)=="Ttpm") then
        read(string,*,iostat=readstat) stmp,ttpm(i,1:10)
      end if

      ! Read in the magnetic transition dipole components for a 
      ! vertical transition.  These have atomic units from NWChem.

      if (string(1:5)=="MTdip") then
        read(string,*,iostat=readstat) stmp,mtdip(i,1:3)
      end if 

      ! Read in the transition quadrupole moment components for a
      ! vertical transition.  These have atomic units from NWChem.
      
      if (string(1:5)=="Tquad") then
        read(string,*,iostat=readstat) stmp,tquad(i,1:6)
      end if

      ! Read the number of normal modes for the purpose of giving each
      ! array a defined length.  Take the vibrational frequency and 
      ! dimensionless delta from the file and store them in the freq(j)
      ! and delta(i,j) arrays.  The frequencies are in wavenumbers at 
      ! this stage of the program.

      if (string(1:11)=="Frequencies") then
        read(string,*,iostat=readstat) stmp,nmodes
        if (iex==1) allocate(freq(nmodes))
        if (iex==1) allocate(delta(nexci,nmodes))
        if (iex==1) allocate(tdipder(nexci,nmodes,3))
        if (iex==1) allocate(gdipder(nexci,nmodes,3))
        if (iex==1) allocate(stpmder(nexci,nmodes,6))
        if (iex==1) allocate(tquadder(nexci,nmodes,6))
        if (iex==1) allocate(mtdipder(nexci,nmodes,3))

        ! Collects data depending on whether Herzberg-Teller terms are
        ! being calculated.  This also includes the IR coupling model 
        ! of Ann Myers Kelly's group.  If the full B term is being 
        ! determined in hyper-Raman, a lot of data is passed into the 
        ! program.

        if ((lhyperraman.or.lashyperraman).and.lfullherztell) then
          do j=1,nmodes
            read(iudata,*) freq(j),delta(i,j),tdipder(i,j,1),tdipder(i,j,2), &
                           tdipder(i,j,3),stpmder(i,j,1),stpmder(i,j,2),     &
                           stpmder(i,j,3),stpmder(i,j,4),stpmder(i,j,5),     &
                           stpmder(i,j,6)
          end do
        else if (ltpa.and.lherztell) then
          do j=1,nmodes
            read(iudata,*) freq(j),delta(i,j),tdipder(i,j,1),tdipder(i,j,2), &
                           tdipder(i,j,3),stpmder(i,j,1),stpmder(i,j,2),     &
                           stpmder(i,j,3),stpmder(i,j,4),stpmder(i,j,5),     &
                           stpmder(i,j,6)
          end do
        else if (lraman.or.labs.or.lhyperraman.or.lfluor.or.lasraman &
               .or.lashyperraman) then
          if ((lherztell.or.lircoupling).and..not.lfullherztell) then
            do j=1,nmodes
              if (lprintGtensor .and. lprintAtensor) then
                read(iudata,*) freq(j),delta(i,j), &
                               tdipder(i,j,1), &
                               tdipder(i,j,2),tdipder(i,j,3), &
                               tquadder(i,j,1),tquadder(i,j,2),tquadder(i,j,3),&
                               tquadder(i,j,4),tquadder(i,j,5),tquadder(i,j,6),&
                               mtdipder(i,j,1),mtdipder(i,j,2),mtdipder(i,j,3)
              else if (lprintGtensor) then
                read(iudata,*) freq(j),delta(i,j), &
                               tdipder(i,j,1), &
                               tdipder(i,j,2),tdipder(i,j,3), &
                               mtdipder(i,j,1),mtdipder(i,j,2),mtdipder(i,j,3)
              else if (lprintAtensor) then
                read(iudata,*) freq(j),delta(i,j), &
                               tdipder(i,j,1), &
                               tdipder(i,j,2),tdipder(i,j,3), &
                               tquadder(i,j,1),tquadder(i,j,2),tquadder(i,j,3),&
                               tquadder(i,j,4),tquadder(i,j,5),tquadder(i,j,6)
              else
                read(iudata,*) freq(j),delta(i,j), &
                               tdipder(i,j,1), &
                               tdipder(i,j,2),tdipder(i,j,3)
              end if
            end do
          else
            do j=1,nmodes
              read(iudata,*) freq(j),delta(i,j)
            end do
          end if
        else if (lsfg) then
          if (lherztell) then
            do j=1,nmodes
              read(iudata,*) freq(j),delta(i,j),gdipder(i,j,1), &
                             gdipder(i,j,2),gdipder(i,j,3) , tdipder(i,j,1), &
                             tdipder(i,j,2),tdipder(i,j,3)
            end do
          else
            do j=1,nmodes
              read(iudata,*) freq(j),delta(i,j),gdipder(i,j,1), &
                             gdipder(i,j,2),gdipder(i,j,3)
            end do
          end if
        ! For the moment, second hyper-Raman only goes through here (A term
        ! only)
        else
          do j=1,nmodes
            read(iudata,*) freq(j),delta(i,j)
          end do
        end if
      end if
      if (string=='End') then
        write(*,'(A39,I2)') '# done reading info for excited state:',iex
        exit 
      end if
    end do
  end do

  close(iudata)

  ! Write broadening information to the output
  call OutputLifetime(nexci, solvmod, gam, gam2, kappa, freqbath, deltabath)

! --------
! end read
! --------
!****************************
! Evaluate HT Term Prefactors
!****************************

! Based on the user input, if Herzberg-Teller terms are calculated,
! this part determines prefactors depending on the calculation
! type.  These prefactors tend to cause most of the peak intensity
! in Herzberg-Teller coupling with the model we use. 

! Transition dipole moment prefactor for B term resonance Raman 
! scattering.  Because the polarizability is a tensor, orientational
! averaging is performed as if the prefactor is a tensor (with 
! integrals set equal to one).

  if (lraman.and.lherztell.and.lhtpref) then
    allocate(rrs_pref(nexci,nmodes))
    call RRSPref(rrs_pref, osc, tdipder, delta, freq, nexci, &
                     nmodes, 3)
  end if

! Prefactor for the resonance hyper-Raman B1 term only, which involve 
! two-photon moment and the transition dipole moment derivatives.

  if ((lhyperraman.and.lherztell.and.lhtpref)) then
    allocate(rhrs_pref(nexci,nmodes))
    call RHRSPrefB1(rhrs_pref, stpm, tdipder, delta, freq, nexci, &
                    nmodes, 3, 6)
  end if

! Prefactors for the resonance hyper-Raman B1 and B2 terms, involving 
! the products of the transition dipole moment and derivatives of the
! two-photon moment, and products of the two-photon moment and the 
! transition dipole moment derivatives

  if ((lhyperraman.and.lfullherztell.and.lhtpref)) then
    allocate(rhrs_pref(nexci,nmodes))
    call RHRSPrefB1(rhrs_pref, stpm, tdipder, delta, freq, nexci, &
                    nmodes, 3, 6)  
    call RHRSPrefB2(rhrs_pref, osc, stpmder, delta, freq, nexci, &
                    nmodes, 3, 6)
  end if

!**************************************
! Convert Response Property Derivatives
!**************************************

! Derivatives are read by TDSPEC in mass-weighted units, so these
! must be converted using the square root term (resulting from the 
! harmonic model/second quantized form of the normal coordinate).

  if (labs) then
    if (lherztell) then
      call UnMassWeight(tdipder,freq,nmodes,nexci,3)
    end if
  end if
  
  if (ltpa) then
    if (lherztell) then
      call UnMassWeightB2(stpmder,freq,nmodes,nexci,6)
    end if
  end if
  
  if (lraman) then
    if (lherztell) then
      call UnMassWeight(tdipder,freq,nmodes,nexci,3)
      call UnMassWeight(tquadder,freq,nmodes,nexci,6)
      call UnMassWeight(mtdipder,freq,nmodes,nexci,3)
    end if
  end if
  
  if (lhyperraman) then
    if (lherztell.or.lfullherztell) then
      call UnMassWeight(tdipder,freq,nmodes,nexci,3)
    end if
    if (lfullherztell) then
      call UnMassWeightB2(stpmder,freq,nmodes,nexci,6)
    end if
  end if
  
  if (lsfg) then
    if (lherztell) then
      call UnMassWeight(tdipder,freq,nmodes,nexci,3)
    end if
    call UnMassWeightSFG(gdipder,freq,nmodes,nexci,3)
  end if

!***************************************
! Determine the Abcissas for Integration
!***************************************

! Define the start time, end time, and number of steps for the 
! numerical integration.  Set up the arrays to store the time 
! (x(t)) and weight (w(t)) integration values.  This information 
! is obtained by Gauss-Legendre quadrature, which tells you which 
! "x-values" to evaluate the integral at numerically.  This is a 
! good method for evaluating integrals due to not requiring each 
! "x-value" to be at a regular interval.

  tstart = 0.0
  allocate(x(tstep))
  allocate(w(tstep))
  call gauleg(tstart,tstop,x,w,tstep)

!*************************************
! Determine the Lifetime function g(t)
!*************************************

  allocate(g(nexci,tstep))

  call Lineshape(solvmod, g, x, gam, gam2, kappa, &
                     freqbath, deltabath, tstep, nexci, temper, &
                     labs, ltpa, lfluor, lraman, lhyperraman, lsfg, &
                     lsecondhyperraman) 

!******************************************************
! Determine the Total Vibrational Reorganization Energy
!******************************************************

  allocate(vibreorg(nexci))

  call VibrationalReorg(freq, delta, vibreorg, nmodes, nexci)

! READ ME PLEASE
! The overlap integrals calculated in this program always assume that
! every mode starts in the ground state.  For instance, a molecule
! with 3 modes has the initial state |0 0 0>, where each number 
! indicates the vibrational quantum number for a particular mode.
! No hot bands (lines resulting from initial states where modes are 
! excited thermally) are accounted for in this code at the moment.
! This is insignificant for the Franck-Condon terms, but it may be 
! important for the correct description of Herzberg-Teller terms.  
! The coupling between normal modes is stored BOTH in the transition 
! dipole moment derivatives and the overlap integrals.

!*********************************
! Collect Absorption Spectrum Data
!*********************************  

! Calculate the absorption cross section using the integration steps 
! from the gauleg subroutine. 

  if (labs) then
    call cpu_time(time1)
    call AbsPlot(lherztell, latermoff, x, w, tstep, osc, tdipder, &
                 freq, delta, lfreq, freq0, g, nexci, nmodes, & 
                 nabspts)
    call cpu_time(time2)
  end if 

!************************
! One-Photon Fluorescence
!************************

! This is implemented based on:
!
! Petrenko, T.; Krylova, O., Neese, F.; and Sokolowski, M.  New J. Phys.
! 11, 015001 (2009).

  if (lfluor) then
    call cpu_time(time1)
    if (solvmod(1)==1) then
      call StokesShift(freq, delta, stokes, stokessolv, freqcut, &
                  nmodes, nexci, fluorstate)
    else
      stokes = zero
    end if
    call FluorPlot(fluorstate, x, w, tstep, osc, freq, delta, lfreq, &
                   freq0, g, nexci, nmodes, nabspts, stokes)
    call cpu_time(time2)
  end if
                                 
!****************************
! Two-Photon Absorption (TPA)
!****************************

! Calculate the two-photon absorption cross sections.  The first part
! of this code uses the paper
!
! Macak, P.; Luo, Y.; and Agren, H.  Chem. Phys. Lett., 330, 447, 
!   (2000).
!
! for the relationship between two-photon transition dipole moments 
! and the two-photon absorption cross section.  The concept for the 
! multiplying by the df and dg factors was taken from the code for 
! linear polarized light from Dalton.

  if (ltpa) then
    call cpu_time(time1)
    call TwoAbsPlot2(lherztell, latermoff, x, w, tstep, &
                 stpm, stpmder, freq, delta, lfreq, &
                 freq0, g, nexci, nmodes, nabspts)
    call cpu_time(time2)
  end if

!************************
! Circular Dichroism (CD)
!************************

! Calculate the circular dichroism spectrum

  if (lcd) then
    call cpu_time(time1)
    call CDPlot(lherztell, latermoff, x, w, tstep, &
                osc, mtdip, freq, delta, lfreq, &
                freq0, g, nexci, nmodes, nabspts)
    call cpu_time(time2)
  end if

!*******************************
! Resonance Raman Data Collector
!*******************************
! Calculate the Raman transition polarizability

! In the code, dsum is the sum of the exponential of the line shape 
! function (see equation 20 in Neese and Petrenko).  The rsum term is
! the part of equation 19 in curly braces for km = 1.  The equation for 
! Ls is giving the integral in equation 31 (without the constants in the
! front of the equation).

! The part involving after evaluating the integral is the generation 
! of the polarizability tensor:
!
! alpha_xx alpha_xy alpha_xz
! alpha_yx alpha_yy alpha_yz
! alpha_zx alpha_zy alpha_zz
!
! Note that this piece adds in the transition dipoles to the equation 
! (these are the array terms called osc(iex,#)).  When we do the 
! numerical integration units of "time" are put into the equation.  For 
! the exponential, the unit of time must be inverse energy, meaning 
! that the unit of time must be cm in the equation if the frequency is 
! in wavenumbers.

  if (lraman) then
    call cpu_time(time1)
    ! Franck-Condon term fundamentals.
    allocate(pol(nmodes,3,3))
    allocate(pola(nmodes,3,3,3))
    allocate(polas(nmodes,3,3,3))
    allocate(polg(nmodes,3,3))
    allocate(polgs(nmodes,3,3))
    allocate(polC(nmodes,6,6))
    allocate(polD(nmodes,6,3))
    allocate(polDs(nmodes,3,6))

    call RRSFCFund(latermoff, x, w, freq, delta, osc, pol, freq0, &
                   lfreq, g, nexci, nmodes, tstep, pola, tquad,   &
                   polg, mtdip, polC, polD, lprintAtensor, lprintGtensor)

    ! Set script-type tensors to be identical to Roman-type tensors in
    ! the FC-approximation
    polas(:,:,:,:) = pola(:,:,:,:)
    do i = 1, nmodes
      polgs(i,1,1) = -polg(i,1,1)
      polgs(i,1,2) = -polg(i,2,1)
      polgs(i,1,3) = -polg(i,3,1)
      polgs(i,2,1) = -polg(i,1,2)
      polgs(i,2,2) = -polg(i,2,2)
      polgs(i,2,3) = -polg(i,3,2)
      polgs(i,3,1) = -polg(i,1,3)
      polgs(i,3,2) = -polg(i,3,2)
      polgs(i,3,3) = -polg(i,3,3)

      polDs(i,1,1) = -polD(i,1,1)
      polDs(i,1,2) = -polD(i,2,1)
      polDs(i,1,3) = -polD(i,3,1)
      polDs(i,1,4) = -polD(i,4,1)
      polDs(i,1,5) = -polD(i,5,1)
      polDs(i,1,6) = -polD(i,6,1)
      polDs(i,2,1) = -polD(i,1,1)
      polDs(i,2,2) = -polD(i,2,1)
      polDs(i,2,3) = -polD(i,3,1)
      polDs(i,2,4) = -polD(i,4,1)
      polDs(i,2,5) = -polD(i,5,1)
      polDs(i,2,6) = -polD(i,6,1)
      polDs(i,3,1) = -polD(i,1,1)
      polDs(i,3,2) = -polD(i,2,1)
      polDs(i,3,3) = -polD(i,3,1)
      polDs(i,3,4) = -polD(i,4,1)
      polDs(i,3,5) = -polD(i,5,1)
      polDs(i,3,6) = -polD(i,6,1)
    end do

    ! Stuff for Franck-Condon term overtones and combination bands.
    if (excnum > 1) then
      allocate(polot(excnum,nmodes,3,3))
      allocate(freqcbtmp(cbguess))
      allocate(excnumstmp(cbguess,nmodes))
      ! Call the subroutine FCfactor, which determines the excitation
      ! numbers for the modes involved in combination bands.  Here, i
      ! and j are dummy indices relating to the number of normal modes.
  
      ! For the array storing the excitation levels (excnums), k is an
      ! arbitrary length, that is simply used to store as many 
      ! combination bands as required.  As far as I'm aware, there is
      ! no simple systematic way to predict this value before counting
      ! the number of combination bands that satisfy the conditions
      ! of having a FC factor larger than fc_screen and a vibrational
      ! frequency lower than freq_screen.
      call FCfactor(freq, delta, excnum, fc_screen, &
                    freq_screen, freqcbtmp, excnumstmp, nexci, &
                    nmodes, numcb, cbguess)
      allocate(freqcb(numcb))              
      allocate(excnums(numcb,nmodes))              
      do i=1,numcb
        freqcb(i) = freqcbtmp(i)
        do j=1,nmodes
          excnums(i,j) = excnumstmp(i,j)
        end do
      end do
      call CombBands(freq, freqcb, excnums, nmodes, numcb)
      allocate(polcb(numcb,3,3))
      ! Franck-Condon term overtones and combination bands.
      call RRSFCOverCB(latermoff, x, w, freq, delta, osc, polot, &
                     polcb, freq0, lfreq, g, nexci, nmodes, &
                     tstep, excnum, numcb, excnums, factorial)
    end if
    ! Herzberg-Teller terms.
    if (lherztell) then
      allocate(Ls_htb(nexci,nmodes))
      allocate(Ls_htb2(nexci,nmodes))
      allocate(Ls_htb3(nexci,nmodes))

      if (lhtints) then
        allocate(htlfreqdiff(nexci))
        allocate(ht_pol(nmodes,nmodes,3,3))
        allocate(ht_rcrs(nmodes,nmodes))
        allocate(ht_rcrs_sum(nmodes))
      end if
      ! Herzberg-Teller term fundamentals.
      call RRSHTFund(x, w, freq, delta, osc, tdipder, pol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     Ls_htb, Ls_htb2, Ls_htb3, htlfreqdiff, &
                     ht_pol, ht_rcrs, ht_rcrs_sum, lhtints, &
                     pola, polg, polC, polD, tquad, tquadder, &
                     mtdip, mtdipder, polas, polgs, polDs, &
                     lprintAtensor, lprintGtensor)
      ! Herzberg-Teller term overtones and combination bands.
      if (excnum > 1) then
        allocate(sumr(nmodes))
        allocate(sumr1(nmodes))
        allocate(sumr2(nmodes))
        call RRSHTOverCB(x, w, freq, delta, osc, tdipder, polot, polcb, &
                       freq0, lfreq, g, nexci, nmodes, tstep, &
                       excnum, numcb, excnums, factorial, rootn, num, &
                       Ls_htb, Ls_htb2, Ls_htb3, sumr, sumr1, &
                       sumr2)
      end if
    end if
    ! Print polarizability derivatives if needed
    if (lprintPol) then
      call PrintPol(freq, lfreq, pol, nmodes, pola, polg, polC, polD, &
                    polas, polgs, polDs, lprintAtensor, lprintGtensor)
    end if
    ! Orientational averaging of the polarizability from fundamentals
    allocate(rcrs(nmodes))
    call RRSFundAvg(pol, rcrs, nmodes)
    ! Orientational averaging of the polarizability from overtones and
    ! combination bands.
    if (excnum > 1) then
      allocate(rcrsot(excnum,nmodes))
      allocate(rcrscb(numcb))
      call RRSOverCBAvg(polot, polcb, rcrsot, rcrscb, nmodes, excnum, &
                        numcb)
    end if
    ! Vibrational overlap integral products
    if (lvoverlap) then
      allocate(voverlap(nexci,nmodes,16))
      call VibOverlap(x, w, freq, delta, voverlap, freq0, nexci, &
                      nmodes, tstep)
    end if
  end if

!*******************************************
! Anti-Stokes Resonance Raman Data Collector
!*******************************************
! Calculate the anti-Stokes Raman transition polarizability.  No
! overtones/combination bands need to be calculated.  It is clearly
! implied in this portion that one of the modes is in a
! vibrational excited state.

  if (lasraman) then
    call cpu_time(time1)
    ! Franck-Condon term fundamentals
    allocate(aspol(nmodes,3,3))
    call ASRRSFCFund(latermoff, x, w, freq, delta, osc, aspol, freq0, &
                   lfreq, g, nexci, nmodes, tstep)
    ! Herzberg-Teller terms.
    if (lherztell) then
      allocate(Ls_htb(nexci,nmodes))
      allocate(Ls_htb2(nexci,nmodes))
      allocate(Ls_htb3(nexci,nmodes))
      allocate(expsum(nmodes))
      ! Herzberg-Teller term fundamentals.
      call ASRRSHTFund(x, w, freq, delta, osc, tdipder, aspol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb, Ls_htb2, Ls_htb3) 
    end if
    ! Orientational averaging of the polarizability from fundamentals
    allocate(rcrs(nmodes))
    call ASRRSFundAvg(aspol, rcrs, nmodes)
  end if

!*********************************************
! Resonance Vibrational Raman Optical Activity
!*********************************************
! Calculate the resonance VROA transition tensors

  if (lrvroa) then
    call cpu_time(time1)
    ! Franck-Condon term fundamentals.
    allocate(pol(nmodes,3,3))
    allocate(polg(nmodes,3,3))
    allocate(polgs(nmodes,3,3))
    allocate(pola(nmodes,3,3,3))
    allocate(polas(nmodes,3,3,3))
    call RVROAFCFund(latermoff, x, w, freq, delta, osc, mtdip, tquad, &
                   pol, polg, polgs, pola, polas, freq0, lfreq, g, &
                   nexci, nmodes, tstep)
  end if

!*************************************
! Resonance Hyper-Raman Data Collector
!*************************************
! Calculate the hyper-Raman transition hyperpolarizability

  if (lhyperraman) then
    call cpu_time(time1)
    ! Franck-Condon term fundamentals
    allocate(hpol(nmodes,3,3,3))
    call RHRSFCFund(latermoff, x, w, freq, delta, osc, stpm, hpol, &
                   freq0, lfreq, g, nexci, nmodes, tstep)
    ! Franck-Condon term overtones and combination bands
    if (excnum > 1) then
      allocate(hpolot(excnum,nmodes,3,3,3))
      allocate(freqcbtmp(cbguess))
      allocate(excnumstmp(cbguess,nmodes))
      ! Call the subroutine FCfactor, which determines the excitation
      ! numbers for the modes involved in combination bands.  Here, i
      ! and j are dummy indices relating to the number of normal modes.

      ! For the array storing the excitation levels (excnums), k is an
      ! arbitrary length, that is simply used to store as many 
      ! combination bands as required.  As far as I'm aware, there is
      ! no simple systematic way to predict this value before counting
      ! the number of combination bands that satisfy the conditions
      ! of having a FC factor larger than fc_screen and a vibrational
      ! frequency lower than freq_screen.
      call FCfactor(freq, delta, excnum, fc_screen, &
                    freq_screen, freqcbtmp, excnumstmp, nexci, &
                    nmodes, numcb, cbguess)
      allocate(freqcb(numcb))              
      allocate(excnums(numcb,nmodes))              
      do i=1,numcb
        freqcb(i) = freqcbtmp(i)
        do j=1,nmodes
          excnums(i,j) = excnumstmp(i,j)
        end do
      end do
      call CombBands(freq, freqcb, excnums, nmodes, numcb)
      allocate(hpolcb(numcb,3,3,3))
      call RHRSFCOverCB(latermoff, x, w, freq, delta, osc, stpm, &
                     hpolot, hpolcb, freq0, lfreq, g, nexci, &
                     nmodes, tstep, excnum, numcb, excnums, factorial)
    end if
    ! Herzberg-Teller terms
    if (lherztell.or.lfullherztell) then
      allocate(Ls_htb2(nexci,nmodes))
      allocate(Ls_htb3(nexci,nmodes))
      allocate(expsum(nmodes))

      if (lhtints) then
        allocate(htlfreqdiff(nexci))
        allocate(ht_hpol(nmodes,nmodes,3,3,3))
        allocate(ht_rcrs(nmodes,nmodes))
        allocate(ht_rcrs_sum(nmodes))
      end if
      ! Herzberg-Teller term fundamentals (B1 term)
      if (lrhrsb2.eqv..false.) then 
        call RHRSB1Fund(x, w, freq, delta, stpm, tdipder, hpol, &
                        freq0, lfreq, g, nexci, nmodes, tstep, &
                        expsum, Ls_htb2, Ls_htb3, htlfreqdiff, &
                        ht_hpol, ht_rcrs, ht_rcrs_sum, lhtints)
      end if
      ! Herzberg-Teller term fundamentals (B2 term)
      if (lfullherztell) then
        allocate(Ls_htb(nexci,nmodes))
        call RHRSB2Fund(x, w, freq, delta, osc, stpmder, hpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb, htlfreqdiff, &
                     ht_hpol, ht_rcrs, ht_rcrs_sum, lhtints)
      end if
      ! Herzberg-Teller term overtones and combination bands
      if (excnum > 1) then
        allocate(sumr(nmodes))
        allocate(sumr1(nmodes))
        allocate(sumr2(nmodes))
        ! B1 term overtones and combination bands
        if (lrhrsb2.eqv..false.) then
          call RHRSB1OverCB(x, w, freq, delta, stpm, tdipder, hpolot, &
                         hpolcb, freq0, lfreq, g, nexci, nmodes, &
                         tstep, excnum, numcb, excnums, factorial, &
                         rootn, expsum, Ls_htb2, Ls_htb3, sumr, &
                         sumr1)
        end if
        if (lfullherztell) then
          ! B2 term overtones and combination bands
          call RHRSB2OverCB(x, w, freq, delta, osc, stpmder, hpolot, &
                         hpolcb, freq0, lfreq, g, nexci, &
                         nmodes, tstep, excnum, numcb, excnums, &
                         factorial, num, expsum, Ls_htb, sumr2)
        end if
      end if
    end if  
    ! IR coupling model for RHRS
    if (lircoupling) then
      allocate(Ls_htir(nexci,nmodes))
      allocate(expsum(nmodes))
      call RHRSIRCoupling(x, w, freq, delta, osc, tdipder, hpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htir)
    end if
    ! -----------------------------------------------
    ! Zhongwei
    ! Print hyperpolarizability derivatives if needed
    if (lprintHpol) then
      call PrintHpol(freq, lfreq, hpol, nmodes)
    end if
    ! -----------------------------------------------
    ! Orientational averaging for RHRS fundamentals
    allocate(rcrs(nmodes))
    call RHRSFundAvg(hpol, rcrs, nmodes)
    ! Orientational averaging for RHRS overtones and combination bands
    if (excnum > 1) then
      allocate(rcrsot(excnum,nmodes))
      allocate(rcrscb(numcb))
      call RHRSOverCBAvg(hpolot, hpolcb, rcrsot, rcrscb, nmodes, &
                         excnum, numcb)
    end if
    ! Vibrational overlap integral products
    if (lvoverlap) then
      allocate(voverlap(nexci,nmodes,16))
      call VibOverlap(x, w, freq, delta, voverlap, freq0, nexci, &
                      nmodes, tstep)
    end if
    ! Hyperpolarizability Cartesian orientation averages
    if (lhpolavg) then
      call HPolAvg(hpol, freq, nmodes, lfreq, lherztell, &
                   lfullherztell, latermoff, lrhrsb2)
    end if
  end if

!*************************************************
! Anti-Stokes Resonance Hyper-Raman Data Collector
!*************************************************
! Calculate the anti-Stokes hyper-Raman transition hyperpolarizability.  
! No overtones/combination bands need to be calculated.  It is clearly
! implied here that one mode starts in an excited vibrational state.

  if (lashyperraman) then
    call cpu_time(time1)
    ! Franck-Condon term fundamentals
    allocate(ashpol(nmodes,3,3,3))
    call ASRHRSFCFund(latermoff, x, w, freq, delta, osc, stpm, ashpol, &
                   freq0, lfreq, g, nexci, nmodes, tstep)
    ! Herzberg-Teller terms
    if (lherztell.or.lfullherztell) then
      allocate(Ls_htb2(nexci,nmodes))
      allocate(Ls_htb3(nexci,nmodes))
      allocate(expsum(nmodes))
      ! Herzberg-Teller term fundamentals (B1 term)
      if (lrhrsb2.eqv..false.) then 
        allocate(Ls_htb(nexci,nmodes))
        call ASRHRSB1Fund(x, w, freq, delta, stpm, tdipder, ashpol, &
                        freq0, lfreq, g, nexci, nmodes, tstep, &
                        expsum, Ls_htb)
      end if
      ! Herzberg-Teller term fundamentals (B2 term)
      if (lfullherztell) then
        call ASRHRSB2Fund(x, w, freq, delta, osc, stpmder, ashpol, &
                     freq0, lfreq, g, nexci, nmodes, tstep, &
                     expsum, Ls_htb2, Ls_htb3)
      end if
    end if
    ! Orientational averaging for RHRS fundamentals
    allocate(rcrs(nmodes))
    call ASRHRSFundAvg(ashpol, rcrs, nmodes)
  end if

!********************************************
! Resonance Second Hyper-Raman Data Collector
!********************************************
! Calculate the Stokes second hyper-Raman transition second hyperpolarizability.
! Presently this is only implemented for A term scattering.

  if (lsecondhyperraman) then
    call cpu_time(time1)
    ! Franck-Condon term fundamentals
    allocate(sechpol(nmodes,3,3,3,3))
    call R2NDHRSFCFund(latermoff, x, w, freq, delta, osc, ttpm, sechpol, &
                       freq0, lfreq, g, nexci, nmodes, tstep)
    ! Orientational averaging for RSHRS fundamentals
    allocate(rcrs(nmodes))
    call R2NDHRSFundAvg(sechpol, rcrs, nmodes)
  end if

!*********************************************************
! Double-Resonance Sum-Frequency Generation Data Collector
!*********************************************************
! Calculate the sum-frequency generation transition hyperpolarizability

  if (lsfg) then
    call cpu_time(time1)
    ! Franck-Condon term fundamentals
    allocate(hpol(nmodes,3,3,3))
    call SFGFCFund(latermoff, x, w, freq, delta, osc, gdipder, hpol, &
                   freq0, lfreq, g, nexci, nmodes, tstep, IRUV)

    ! Herzberg-Teller terms
    if (lherztell) then
       allocate(Ls_htb(nexci,nmodes))
       allocate(Ls_htb2(nexci,nmodes))
       allocate(Ls_htb3(nexci,nmodes))
       allocate(expsum(nmodes))

       call SFGHTFund(x, w, freq, delta, osc, tdipder, gdipder, hpol, &
                      freq0, lfreq, g, nexci, nmodes, tstep, &
                      expsum, Ls_htb, Ls_htb2, Ls_htb3, IRUV)
    end if

    allocate(fhpol(nmodes,3,3,3))
    fhpol = hpol

    ! Add in the non-resonant SFG signal
    if (lnonRes) then
        call SFGnonRes(nmodes, hpol,width,shpol)
    end if


    ! Does the Full DR-SFG calculation
    if (lSFGFull) then
      ! The number of wavenumbers shifted from the max and min frequencies  to
      ! account for the contributions to peaks along the edge of the desired scan
      ! range
      widthrng = 100 
      ! Readjust freqend to the number of points that need to be calculated
      freqend = (freqend - freqstrt + 2*widthrng)/res + 1
      ! Change freqstrt to start of the actual calculation frequency
      freqstrt = freqstrt - widthrng - res
      allocate (nhpol(freqend, 3,3,3))
      call SFGFullCalc(freq, fhpol, nmodes, width, freqstrt, freqend, res,&
                       nhpol, lnonRes, shpol) 
      allocate (rcrs_fsfg(freqend))
      call SFGFullAvg(nhpol, psi, theta, orientation, rcrs_fsfg ,freqend)
    end if
    ! Orientational averaging for DR-SFG fundamentals
    allocate(rcrs(nmodes))
    call SFGFundAvg(hpol, width, psi, theta, orientation, rcrs, nmodes)
  end if

!***************************
! Unit sphere representation
!***************************
! This part evaluates the unit sphere for representation of a
! polarizabilty tensor.  For more information, see:
!
! Tuer, A.; Krouglov, S.; Cisek, R.; Tokarz, D.; and Barzda, V.  
!   J. Comput. Chem. 32, 1128 (2011).
! Tuer, A. E.; Krouglov, S.; Prent, N.; Cisek, R.; Sandkuijl, D.; 
!   Yasufuku, K.; Wilson, B. C.; and Barzda, V.  J. Phys. Chem. B
!   115, 12759 (2011).

  if (lraman.and.lunitsphere) then
    ! Unit sphere representation of the polarizability (for Stokes
    ! scattering)
    call UnitSpherePol(pol, freq, nmodes, lfreq, lherztell, &
                    latermoff, sphereradius) 
  else if (lasraman.and.lunitsphere) then
    ! Unit sphere representation of the polarizability (for anti-Stokes
    ! scattering)
    call UnitSpherePol(aspol, freq, nmodes, lfreq, lherztell, &
                    latermoff, sphereradius) 
  end if

  if ((lhyperraman.and.lunitsphere).or.(lsfg.and.lunitsphere)) then
    ! Unit sphere representation of the hyperpolarizability (for Stokes
    ! scattering)
    call UnitSphereHpol(hpol, freq, nmodes, lfreq, lherztell, &
                    lfullherztell, latermoff, lrhrsb2, sphereradius)
  else if (lashyperraman.and.lunitsphere) then
    ! Unit sphere representation of the hyperpolarizability (for 
    ! anti-Stokes scattering)
    call UnitSphereHpol(ashpol, freq, nmodes, lfreq, lherztell, &
                    lfullherztell, latermoff, lrhrsb2, sphereradius)
  end if

! Free up some memory.

  if (lraman.or.lhyperraman.or.lsfg.or.lasraman.or.lashyperraman.or. &
      lrvroa.or.lsecondhyperraman) then
    deallocate(delta)
    deallocate(freq0)
    deallocate(osc)
    deallocate(mtdip)
    deallocate(tquad)
    deallocate(stpm)
    deallocate(ttpm)
    deallocate(gam)
    deallocate(gdipder)
    deallocate(tdipder)
    deallocate(stpmder)
    deallocate(x)
    deallocate(w)
    if (lvoverlap) then
      deallocate(voverlap)
    end if
  end if

!**********************************************
! Evaluate the Boltzmann factors for every mode
!**********************************************

! Boltzmann factors are evaluated based on the general expression
! the population of state I.

! For RRS and RHRS, we assume the initial state is the ground 
! vibrational state for all modes.  This is a good representation
! of the overall population for most modes at room temperature, where 
! the value of k_B*T ~ 210 cm^{-1}

  if (lraman.or.lhyperraman.or.lsecondhyperraman) then
    allocate(boltzpop(nmodes))
    call Boltzmann(freq, boltzpop, nmodes, temper)
    if (excnum > 1) then
      allocate(boltzpopot(excnum, nmodes))
      allocate(boltzpopcb(numcb))
      call BoltzmannOvrt(boltzpopot, boltzpop, nmodes, excnum)
      call BoltzmannComb(freq, boltzpopcb, excnums, nmodes, numcb, temper)
    end if
  end if

! For anti-Stokes versions of RRS and RHRS, we assume the initial
! state is the first vibrational excited state of the molecule.

  if (lasraman.or.lashyperraman) then
    allocate(boltzpopas(nmodes))
    call BoltzmannAS(freq, boltzpopas, nmodes, temper)
  end if

!******************************************************
! Calculate differential Raman scattering cross section
!******************************************************
  if (lraman) then
    ! Fundamentals
    call RamanCrossSection(nmodes, boltzpop, lfreq, freq, rcrs)
    if (excnum > 1) then
      ! Overtone and combination band contributions
      call RamanCrossSectionOverCB(excnum, numcb, nmodes, freq, freqcb, &
                     lfreq, boltzpopot, rcrsot, boltzpopcb, rcrscb)
    end if
  end if

!**********************************************************************
! Calculate the anti-Stokes differential Raman scattering cross section
!**********************************************************************
  if (lasraman) then
    call ASRamanCrossSection(nmodes, boltzpopas, lfreq, freq, rcrs)
  end if

!****************************************************************
! Calculate the differential hyper-Raman scattering cross section
!****************************************************************
  if (lhyperraman) then
    call HypRamanCrossSection(nmodes, boltzpop, lfreq, freq, rcrs)
    if (excnum > 1) then
      call HypRamanCrossSectionOverCB(excnum, numcb, nmodes, freq, freqcb, &
                     lfreq, boltzpopot, rcrsot, boltzpopcb, rcrscb)
    end if
  end if

!****************************************************************************
! Calculate the anti-Stokes differential hyper-Raman scattering cross section
!****************************************************************************
  if (lashyperraman) then
    call ASHypRamanCrossSection(nmodes, boltzpopas, lfreq, freq, rcrs)
  end if

!***********************************************************************
! Calculate the differential second hyper-Raman scattering cross section
!***********************************************************************
  if (lsecondhyperraman) then
    call SecHypRamanCrossSection(nmodes, boltzpop, lfreq, freq, rcrs)
  end if

  ! Build an array containing the fundamentals, overtones, and 
  ! combination bands (FOC).
  
  if (lraman.or.lhyperraman.or.lsfg.or.lsecondhyperraman) then

    allocate(freq_foc(excnum*nmodes+numcb))
    allocate(rcrs_foc(excnum*nmodes+numcb))
  
    do i=1,nmodes
      freq_foc(i) = freq(i)
      rcrs_foc(i) = rcrs(i)
    end do
 
    do k=2,excnum
      r = k-1
      do j=r*nmodes+1,k*nmodes
        freq_foc(j) = k*freq(j-r*nmodes)
        rcrs_foc(j) = rcrsot(k,j-r*nmodes)
      end do
    end do

    do n=excnum*nmodes+1,excnum*nmodes+numcb
      freq_foc(n) = freqcb(n-excnum*nmodes)
      rcrs_foc(n) = rcrscb(n-excnum*nmodes)
    end do
 
    maxfrq = MAXVAL(freq_foc)

    ! Build an array containing the frequencies for the Full DR-SFG spectrum 
    if (lSFGFull) then
      allocate(freq_fsfg(freqend))
      do i=1, freqend
        freq_fsfg(i)= freqstrt + res*i
      end do
    end if 
  end if

  ! Build an array of fundamentals for anti-Stokes spectra.

  if (lasraman.or.lashyperraman) then
    allocate(freq_foc(nmodes))
    allocate(rcrs_foc(nmodes))

    ! Replace the vibrational frequencies with their negatives.

    do i=1,nmodes
      freq_foc(i) = -freq(i)
      rcrs_foc(i) = rcrs(i)
    end do

    minfrq = MINVAL(freq_foc)
    maxfrq = MAXVAL(freq_foc)

  end if

  !**********************************************
  ! Convolute the Plot with a Lorentzian Function
  !**********************************************
  if (.not.lnormalize) then
    ! Stokes spectra.
    if (lraman.or.lhyperraman.or.lsfg.or.lsecondhyperraman) then
      call StokesLorentzianConv(npoints, nmodes, excnum, numcb, &
                           freq_foc, rcrs_foc, width, maxfrq)
    end if
    ! Anti-Stokes spectra.
    if (lasraman.or.lashyperraman) then
      call AntiStokesLorentzianConv(npoints, nmodes, freq_foc, &
                               rcrs_foc, width, minfrq, maxfrq)
    end if

  !********************
  ! Normalize the plots
  !********************
  else if (lnormalize) then
    allocate(norm_crs(npoints*excnum+1))
    allocate(norm_frq(npoints*excnum+1))
    if (lraman.or.lhyperraman.or.lsfg.or.lsecondhyperraman) then
      call NormStokesLorentzianConv(npoints, nmodes, excnum, numcb, &
                               freq_foc, rcrs_foc, width, maxfrq, &
                               norm_frq, norm_crs, rmax)
    end if
  end if

  !*********************************************
  ! Code for Printing the Unbroadened Raman Data
  !*********************************************  

  call StickSpectra(lintensity, lraman, lhyperraman, lsfg, lsecondhyperraman, &
                latermoff, lherztell, lfullherztell, lrhrsb2, &
                lnormalize, lfreq, width, freq_foc, rcrs_foc, &
                nmodes, excnum, numcb, rmax, theta, psi, orientation)

  ! Determine if the calculation includes the FC, HT, or both terms.

  if ((lherztell.and.lraman.and.latermoff.eqv..false.).or. &
      (lfullherztell.and.lhyperraman.and.latermoff.eqv..false..and. &
       lrhrsb2.eqv..false.).or. &
      (lherztell.and.lsfg.and.latermoff.eqv..false.)) then
    charfcht = 'AB'
  else if (lherztell.and.lhyperraman.and.latermoff.eqv..false.) then
    charfcht = 'AB1'
  else if (lherztell.and.lhyperraman.and.latermoff) then
    charfcht = 'B1'
  else if (lfullherztell.and.lrhrsb2.and.lhyperraman.and. &
           latermoff.eqv..false.) then
    charfcht = 'AB2'
  else if (lrhrsb2.and.lhyperraman.and.latermoff) then
    charfcht = 'B2'
  else if ((lherztell.and.lraman.and.latermoff).or. &
           (lfullherztell.and.lhyperraman.and.latermoff).or. &
           (lherztell.and.lsfg.and.latermoff)) then
    charfcht = 'B'
  else
    charfcht = 'A'
  end if

  charfcht = adjustl(charfcht)

  ! Full DR-SFG Spectrum
  ! More or less the same way the stick spectra is implemented. Will create an
  ! additional file so you may have up to 3 output files from 1 TDSPEC
  ! calculation assuming you a full and stick calc for SFG
  if (lSFGFull) then
    ! Convert the laser frequency and theta/psi to characters.  Also
    ! convert the polarization to a character.
    intfreq = nint(1.0E+7_kindr/lfreq) ! Laser frequency
    write(charfreq,*) intfreq
    charfreq = adjustl(charfreq)
    inttheta = nint(theta*180.0_kindr/pi) ! Angle theta
    write(chartheta,*) inttheta
    chartheta = adjustl(chartheta)
    intpsi = nint(psi*180.0_kindr/pi) ! Angle psi
    write(charpsi,*) intpsi
    charpsi = adjustl(charpsi)
    write(charpol,*) orientation ! Polarization
    charpol = adjustl(charpol)
    ! Write the filename (concatenate strings)
    fname = 'DRSFG_' // trim(charfreq) // 'exc_p' // trim(charpsi) // &
            '_t' // trim(chartheta) // '_pol' // trim(charpol) // &
            '_' // trim(charfcht) // 'term_full.out'
    open(93,file=fname,form='formatted')
    ! gfortran throws out a warning due to this line.  Expressions defining
    ! the limits of loops is apparently a deprecated feature.  
    ! do i=1, freqend-(2*widthrng/res) ! Adjust the range to the specified input
    iend = freqend-(2*widthrng/res) ! Adjust the range to the specified input
    do i=1, iend
      ! Adjust the index so the starting wavenumber is the actually the one the
      ! user wanted
      j = i + widthrng/res 
      if (lnormalize.eqv..true.) then
        write(93,*) freq_fsfg(j),rcrs_fsfg(j)*two/(pi*width*rmax)
      else
        write(93,*) freq_fsfg(j),rcrs_fsfg(j)*two/(pi*width)
      end if
    end do
    close(93)
  end if

  ! Finish timing.  The final time was already determined above for 
  ! absorbance-type spectra.
  if ((labs.eqv..false.).and.(lfluor.eqv..false.).and. &
    (ltpa.eqv..false.).and.(lcd.eqv..false.)) then
    call cpu_time(time2)
  end if
  call Timing(labs, lraman, lasraman, ltpa, lhyperraman, &
             lashyperraman, lsfg, lfluor, lcd, lrvroa, lsecondhyperraman, &
             time1, time2)
  
end program TDSPEC
