module constants

! ========================================
! purpose: Define constants used in TDSPEC
! ========================================

  implicit none

!----------------------------------------
! Double precision
!----------------------------------------
  integer, parameter :: kindr = kind(0d0)
!----------------------------------------

  real(kindr), parameter :: zero  = 0.0_kindr
  real(kindr), parameter :: one   = 1.0_kindr
  real(kindr), parameter :: two   = 2.0_kindr
  real(kindr), parameter :: three = 3.0_kindr
  real(kindr), parameter :: four  = 4.0_kindr
  real(kindr), parameter :: five  = 5.0_kindr
  real(kindr), parameter :: six   = 6.0_kindr
  real(kindr), parameter :: seven = 7.0_kindr
  real(kindr), parameter :: eight = 8.0_kindr
  real(kindr), parameter :: nine  = 9.0_kindr
  real(kindr), parameter :: ten   = 10.0_kindr
  real(kindr), parameter :: eleven = 11.0_kindr
  real(kindr), parameter :: thrtn = 13.0_kindr
  real(kindr), parameter :: sixtn = 16.0_kindr
  real(kindr), parameter :: twentyone = 21.0_kindr
  real(kindr), parameter :: thrty  = 30.0_kindr
  real(kindr), parameter :: thrtfv = 35.0_kindr
  real(kindr), parameter :: fortyfv = 45.0_kindr
  real(kindr), parameter :: sixty   = 60.0_kindr
  real(kindr), parameter :: sixtythree = 63.0_kindr
  real(kindr), parameter :: sixtyfour = 64.0_kindr
  real(kindr), parameter :: seventy = 70.0_kindr
  real(kindr), parameter :: hndr   = 100.0_kindr
  real(kindr), parameter :: hndrfv = 105.0_kindr
  real(kindr), parameter :: hndrtwnty = 120.0_kindr
  real(kindr), parameter :: hndreghty = 180.0_kindr
  real(kindr), parameter :: twhndrten = 210.0_kindr
  real(kindr), parameter :: thrhndrftn = 315.0_kindr
  real(kindr), parameter :: f12 = one / two
  real(kindr), parameter :: f13 = one / three
  real(kindr), parameter :: f15 = one / five
  real(kindr), parameter :: f32 = three / two
  real(kindr), parameter :: f52 = five / two
  real(kindr), parameter :: f72 = seven / two
  real(kindr), parameter :: f92 = nine / two
  real(kindr), parameter :: f112 = eleven / two
  real(kindr), parameter :: root2 = dsqrt(two)
  real(kindr), parameter :: root8 = dsqrt(eight)

  real(kindr), parameter :: hartree2cm = 219474.63068_kindr
  real(kindr), parameter :: pi = 3.14159265358979323846_kindr
  real(kindr), parameter :: planck = 6.6260755E-34_kindr
  real(kindr), parameter :: planckbar = planck/(two*pi)
  real(kindr), parameter :: boltz = 1.3806580E-23_kindr
  real(kindr), parameter :: speed = 2.99792458E+8_kindr
  real(kindr), parameter :: epsilon0 = 8.8541878E-12_kindr
  real(kindr), parameter :: avogadro = 6.02214199E+23_kindr
  real(kindr), parameter :: amu = (one/avogadro)*1E-03_kindr
  real(kindr), parameter :: exparg = planck*speed*100.0_kindr/boltz
  real(kindr), parameter :: alphaf = one/137.03599_kindr
  real(kindr), parameter :: abscrfac = (four*pi*alphaf)/(three)
  real(kindr), parameter :: bohr2ang = 0.5291772_kindr
  real(kindr), parameter :: wvn2joules = 1.9864456E-23_kindr
  real(kindr), parameter :: joule2ev = 6.24152862E+18_kindr
  real(kindr), parameter :: dipconv = one/1.1794744E+29_kindr
  real(kindr), parameter :: fundcharge = 1.6021765E-19_kindr
  real(kindr), parameter :: bohr2cm = 5.291772E-9_kindr
  real(kindr), parameter :: bohr2m = bohr2cm/100_kindr
  real(kindr), parameter :: m2bohr = 18897261328.85643_kindr
  real(kindr), parameter :: wvn2hz2au = (100_kindr*speed*2.418884E-17_kindr)
  real(kindr), parameter :: au2gm = bohr2cm**four*2.418884E-17_kindr
  real(kindr), parameter :: top = dsqrt(planckbar/(four*pi*speed*100_kindr))
  real(kindr), parameter :: bottom = (dsqrt(amu)*bohr2m)
  real(kindr), parameter :: htmassweight = top/bottom
  real(kindr), parameter :: hartree2joules = 4.359743357E-18_kindr
  real(kindr), parameter :: joules2hartree = one/hartree2joules
  real(kindr), parameter :: temp = 300_kindr
  real(kindr), parameter :: kbT = boltz*temp/(planck*speed*hndr)
  real(kindr), parameter :: tentomeight = 1.0E-8_kindr
  real(kindr), parameter :: cm2hartree = one / hartree2cm

  !********************************
  ! Conversions associated with OPA
  !********************************

  real(kindr), parameter :: absconv = abscrfac*bohr2ang**2

  !*****************************************
  ! Conversions associated with Fluorescence
  !*****************************************

  real(kindr), parameter :: fl1 = four*speed**3*hndr**3
  real(kindr), parameter :: fl2 = three*pi*planckbar**2*speed**3*hndr**3
  real(kindr), parameter :: fl3 = bohr2ang**2*(1E-8_kindr)**2*wvn2joules 
  real(kindr), parameter :: fl4 = speed*hndr*joule2ev 
  real(kindr), parameter :: fluorconv = (fl1*fl3)/(fl2*fl4) 

  !********************************
  ! Conversions associated with RRS
  !********************************

  real(kindr), parameter :: pol2si = dipconv**2/wvn2joules
  real(kindr), parameter :: rrstop = pi**2*1E+12_kindr
  real(kindr), parameter :: rrsbottom = (fortyfv*epsilon0**2)
  real(kindr), parameter :: rrscrsfac = rrstop/rrsbottom

  !********************************
  ! Conversions associated with TPA
  !********************************

  real(kindr), parameter :: tpa1 = four*pi**3*alphaf**2*hartree2cm*au2gm
  real(kindr), parameter :: tpa2 = four*pi**2*wvn2hz2au**2
  real(kindr), parameter :: tpaconv = tpa1*tpa2

  !*********************************
  ! Conversions associated with RHRS
  !*********************************

  real(kindr), parameter :: hpol2si = dipconv**3*hartree2cm/hartree2joules**2
  real(kindr), parameter :: rhrstop = sixtn*pi**2*planck**3*alphaf**3*1E+8_kindr 
  real(kindr), parameter :: rhrsbottom = speed**2*fundcharge**6 
  real(kindr), parameter :: rhrscrstmp = rhrstop/rhrsbottom 
  real(kindr), parameter :: photonnorm = two*planck*speed*hndr
  real(kindr), parameter :: rhrscrsfac = rhrscrstmp*photonnorm 

  !***********************************
  ! Conversions associated with DR-SFG
  !***********************************

  ! To convert hyperpolarizabilities in a.u. to SI units
  ! 1 e*a0 = 8.478352730673933e-30 C*m
  ! 1 Eh = 4.359743357E-18 J
  ! au2si = [1 (e*a0)^3/Eh^2]*(8.47835e-30 C*m/1 e*a0)^3*
  !         (1 Eh/4.35974e-18 J)^2 = 3.206362e-53 C^3 m^3/J^2

  real(kindr), parameter :: au2si = 3.2063622989061867E-53_kindr

  ! To convert a.u. to so-called alternative SI units
  ! epsilon0 = 8.8541878E-12 C^2 m^{-1} J^{-1}
  ! au2altsi = 3.206362e-53 C^3 m^3/J^2 / 8.8541878E-12 C^2 m^{-1} J^{-1}
  !          = 3.621294e-42 m^4/V

  real(kindr), parameter :: au2altsi = 3.621294658902748E-42_kindr

  !*******************
  ! Unit sphere radius
  !*******************

  real(kindr), parameter :: unitradius = one 

  !**********************************
  ! Conversions associated with RSHRS
  !**********************************

  real(kindr), parameter :: sechpol2si = dipconv**4*hartree2cm/hartree2joules**3
  real(kindr), parameter :: r2ndhrstop = sixtyfour*pi**2*planck**4*alphaf**4*1E+12_kindr 
  real(kindr), parameter :: r2ndhrsbottom = speed**2*fundcharge**8 
  real(kindr), parameter :: r2ndhrscrstmp = r2ndhrstop/r2ndhrsbottom 
  real(kindr), parameter :: photonnorm2 = three*planck*speed*hndr
  real(kindr), parameter :: r2ndhrscrsfac = r2ndhrscrstmp*photonnorm2*photonnorm2 

end module constants
