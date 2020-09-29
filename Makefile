# Start of Makefile
# Defining variables
#=============
# For gfortran
#=============
f90 = gfortran
flags = -Ofast -fopenmp

#=============
# For ifort
#=============
#f90 = ifort
#flags = -O2 -xHost -openmp -heap-arrays

#=============
# For f2py
#=============
f2py = f2py
f2pyflags = -c --f90flags='${flags}' -lgomp

#=======================
# Flags for bug checking
#=======================
#flags = -g -traceback

#====================================
# Flags for more verbose bug checking
#====================================
# For ifort
#flags = -g -traceback -check all -warn all
# For gfortran
#flags = -g -fbacktrace -fdump-core -ffpe-trap=invalid,zero,overflow -fbounds-check -Wall -Wextra

# Makefile
tdpsec: fact.o constants.o error.o output.o nexstate.o \
        position.o fcfactor.o gauleg.o lifetime.o abscalc.o stokes.o \
        fluorcalc.o twoabscalc.o cdcalc.o rrscalc.o asrrscalc.o rvroacalc.o \
        rrspref.o rrsaverage.o rhrscalc.o asrhrscalc.o rhrsircoupling.o \
	rhrspref.o rhrsaverage.o sfgcalc.o sfgaverage.o htints.o \
	r2ndhrscalc.o r2ndhrsaverage.o sechyperramancrs.o \
	viboverlap.o unmassweight.o combbands.o lifetime.o \
        vibreorg.o boltzmann.o hpolaverage.o unitsphere.o \
        ramancrs.o hyperramancrs.o lconvolute.o stickspec.o timing.o tdspec.o
	$(f90) -o tdspec $(flags) fact.o constants.o error.o output.o nexstate.o \
                  position.o fcfactor.o gauleg.o lifetime.o abscalc.o \
                  stokes.o fluorcalc.o twoabscalc.o cdcalc.o rrscalc.o \
                  asrrscalc.o rvroacalc.o rrspref.o rrsaverage.o rhrscalc.o \
                  asrhrscalc.o rhrsircoupling.o rhrspref.o \
                  rhrsaverage.o sfgcalc.o sfgaverage.o htints.o \
		  r2ndhrscalc.o r2ndhrsaverage.o sechyperramancrs.o \
                  viboverlap.o unmassweight.o combbands.o \
                  vibreorg.o boltzmann.o hpolaverage.o unitsphere.o \
                  ramancrs.o hyperramancrs.o lconvolute.o stickspec.o timing.o tdspec.o
constants.mod: constants.o constants.f90
	$(f90) $(flags) -c constants.f90
constants.o: constants.f90
	$(f90) $(flags) -c constants.f90
fact.mod: fact.o fact.f90
	$(f90) $(flags) -c fact.f90
fact.o: fact.f90
	$(f90) $(flags) -c fact.f90
error.o: error.f90
	$(f90) $(flags) -c error.f90
output.o: output.f90
	$(f90) $(flags) -c output.f90
nexstate.o: nexstate.f90
	$(f90) $(flags) -c nexstate.f90
position.o: position.f90
	$(f90) $(flags) -c position.f90
fcfactor.o: fcfactor.f90
	$(f90) $(flags) -c fcfactor.f90
gauleg.o: gauleg.f90
	$(f90) $(flags) -c gauleg.f90
lifetime.o: lifetime.f90
	$(f90) $(flags) -c lifetime.f90
vibreorg.o: vibreorg.f90
	$(f90) $(flags) -c vibreorg.f90
abscalc.o: abscalc.f90
	$(f90) $(flags) -c abscalc.f90
stokes.o: stokes.f90
	$(f90) $(flags) -c stokes.f90
fluorcalc.o: fluorcalc.f90
	$(f90) $(flags) -c fluorcalc.f90
twoabscalc.o: twoabscalc.f90
	$(f90) $(flags) -c twoabscalc.f90
cdcalc.o: cdcalc.f90
	$(f90) $(flags) -c cdcalc.f90
rrscalc.o: rrscalc.f90
	$(f90) $(flags) -c rrscalc.f90
asrrscalc.o: asrrscalc.f90
	$(f90) $(flags) -c asrrscalc.f90
rvroacalc.o: rvroacalc.f90
	$(f90) $(flags) -c rvroacalc.f90
rrspref.o: rrspref.f90
	$(f90) $(flags) -c rrspref.f90
rrsaverage.o: rrsaverage.f90
	$(f90) $(flags) -c rrsaverage.f90
rhrscalc.o: rhrscalc.f90
	$(f90) $(flags) -c rhrscalc.f90
asrhrscalc.o: asrhrscalc.f90
	$(f90) $(flags) -c asrhrscalc.f90
rhrsircoupling.o: rhrsircoupling.f90
	$(f90) $(flags) -c rhrsircoupling.f90
rhrspref.o: rhrspref.f90
	$(f90) $(flags) -c rhrspref.f90
rhrsaverage.o: rhrsaverage.f90
	$(f90) $(flags) -c rhrsaverage.f90
sfgcalc.o: sfgcalc.f90
	$(f90) $(flags) -c sfgcalc.f90
sfgaverage.o: sfgaverage.f90
	$(f90) $(flags) -c sfgaverage.f90
htints.o: htints.f90
	$(f90) $(flags) -c htints.f90
r2ndhrscalc.o: r2ndhrscalc.f90
	$(f90) $(flags) -c r2ndhrscalc.f90
r2ndhrsaverage.o: r2ndhrsaverage.f90
	$(f90) $(flags) -c r2ndhrsaverage.f90
sechyperramancrs.o: sechyperramancrs.f90
	$(f90) $(flags) -c sechyperramancrs.f90
viboverlap.o: viboverlap.f90
	$(f90) $(flags) -c viboverlap.f90
unmassweight.o: unmassweight.f90
	$(f90) $(flags) -c unmassweight.f90
combbands.o: combbands.f90
	$(f90) $(flags) -c combbands.f90
boltzmann.o: boltzmann.f90
	$(f90) $(flags) -c boltzmann.f90
hpolaverage.o: hpolaverage.f90
	$(f90) $(flags) -c hpolaverage.f90
unitsphere.o: unitsphere.f90
	$(f90) $(flags) -c unitsphere.f90
ramancrs.o: ramancrs.f90
	$(f90) $(flags) -c ramancrs.f90
hyperramancrs.o: hyperramancrs.f90
	$(f90) $(flags) -c hyperramancrs.f90
lconvolute.o: lconvolute.f90
	$(f90) $(flags) -c lconvolute.f90
stickspec.o: stickspec.f90
	$(f90) $(flags) -c stickspec.f90
timing.o: timing.f90
	$(f90) $(flags) -c timing.f90
tdspec.o: tdspec.f90
	$(f90) $(flags) -c tdspec.f90

clean:
	rm fact.mod fact.o constants.mod constants.o error.o output.o nexstate.o \
           position.o fcfactor.o gauleg.o lifetime.o abscalc.o \
           stokes.o fluorcalc.o twoabscalc.o cdcalc.o rrscalc.o \
           asrrscalc.o rvroacalc.o rrspref.o rrsaverage.o rhrscalc.o \
           asrhrscalc.o rhrsircoupling.o rhrspref.o \
           rhrsaverage.o sfgcalc.o sfgaverage.o htints.o \
	   r2ndhrscalc.o r2ndhrsaverage.o sechyperramancrs.o \
           viboverlap.o unmassweight.o combbands.o vibreorg.o \
           boltzmann.o hpolaverage.o unitsphere.o ramancrs.o \
           hyperramancrs.o lconvolute.o stickspec.o timing.o tdspec.o tdspec_pymodule.so
cleanall:
	rm fact.mod fact.o constants.mod constants.o error.o output.o nexstate.o \
           position.o fcfactor.o gauleg.o lifetime.o abscalc.o \
           stokes.o fluorcalc.o twoabscalc.o cdcalc.o rrscalc.o \
           asrrscalc.o rvroacalc.o rrspref.o rrsaverage.o \
           rhrscalc.o asrhrscalc.o rhrsircoupling.o rhrspref.o \
           rhrsaverage.o sfgcalc.o sfgaverage.o htints.o \
	   r2ndhrscalc.o r2ndhrsaverage.o sechyperramancrs.o \
           viboverlap.o unmassweight.o combbands.o vibreorg.o \
           boltzmann.o tdspec.o hpolaverage.o unitsphere.o \
           ramancrs.o hyperramancrs.o lconvolute.o stickspec.o \
           timing.o tdspec

# End of Makefile
