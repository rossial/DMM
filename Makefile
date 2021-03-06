FC = gfortran
FCFLAGS = -cpp -fcray-pointer -ffixed-line-length-none -fdefault-integer-8
FFFLAGS = $(FCFLAGS) -std=legacy

MOD_OBJS = \
	iso_c_utilities.o \
	dlfcn.o

OBJS = \
	amh2.o \
	gammln.o \
	int2seq2.o \
	lyapunov.o \
	opg.o \
	simdata.o \
	amh.o \
	gck2.o \
	int2seq.o \
	main.o \
	opgh.o \
	simprior.o \
	checkdesign.o \
	gck.o \
	invfbis.o \
	markovp.o \
	opgkim.o \
	simstate2.o \
	chi2inv.o \
	harmonic2.o \
	invf.o \
	mengwong2.o \
	pprod.o \
	simstate.o \
	cumnorm.o \
	harmonic.o \
	invnormcdf.o \
	mengwong.o \
	priordir.o \
	slice2.o \
	designz.o \
	hf.o \
	kf2.o \
	missing.o \
	prior.o \
	slice.o \
	drawpsi.o \
	ikf2.o \
	kf.o \
	mvncdf.o \
	ptheta2.o \
	syminv.o \
	drawtheta2.o \
	ikf.o \
	kim.o \
	mvnpdf.o \
	ptheta.o \
	tnormi.o \
	drawtheta.o \
	initrand.o \
	ks2.o \
	neweywestcov2.o \
	var.o \
	ergodic.o \
	innov2.o \
	ks.o \
	neweywestcov.o \
	recest.o \
	findinput.o \
	innov.o \
	lemma4.o \
	ols.o \
	recpr.o \
	forecast.o \
	input.o \
	logmvnpdf.o \
	openfiles.o \
	schollu.o \
	logical2integer.o

RANDLIB_OBJS = \
	advnst.o \
	genmn.o \
	genunf.o \
	ignlgi.o \
	inrgcm.o \
	phrtsd.o \
	sdot.o \
	setgmn.o \
	sgamma.o \
	genbet.o \
	gennor.o \
	getcgn.o \
	ignuin.o \
	lennob.o \
	qrgnin.o \
	setall.o \
	setsd.o \
	snorm.o \
	gengam.o \
	genprm.o \
	getsd.o \
	initgn.o \
	mltmod.o \
	ranf.o \
	setant.o \
	sexpo.o \
	spofa.o

EXEC = dmm

VPATH := $(VPATH) randlib

ifdef DLL
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
MATLAB_LIBS = -L$(MATLABROOT)/bin/glnxa64 -Wl,-rpath-link,$(MATLABROOT)/bin/glnxa64 -Wl,-rpath,$(MATLABROOT)/bin/glnxa64
endif
ifeq ($(UNAME_S), Darwin)
MATLAB_LIBS = -L$(MATLABROOT)/bin/maci64
endif
LIBS = $(MATLAB_LIBS) -llapack -leng -lmx
DLL_OBJS += design.o setfilem.o geterrstr.o
INCLUDE = -I$(MATLABROOT)/extern/include
DEFINE = -DORIGDLL
else
LIBS = -llapack
endif

all: $(MOD_OBJS) $(OBJS) $(RANDLIB_OBJS) $(DLL_OBJS)
	$(FC) $(FCFLAGS) $^ $(LIBS) -o $(EXEC)

%.o: %.f90 %.mod
	$(FC) $(FCFLAGS) $(INCLUDE) $(DEFINE) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) $(INCLUDE) $(DEFINE) -c $<

%.mod: %.f90 %.o
	@true

%.o : %.f
	$(FC) $(FFFLAGS) $(INCLUDE) $(DEFINE) -c $<

%.o : %.for
	$(FC) $(FFFLAGS) $(INCLUDE) $(DEFINE) -c $<

clean:
	rm -f *.o $(EXEC) *.mod *.PRI *.DIS *.FST *.INN *.ML *.PAR *.UNB
