
.KEEP_STATE:

.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libisospin.a

include FlagsForMake

.c.o:
	$(CC) -c $(CFLAGS) $<
.c.a:
	$(CC) -c $(CFLAGS) $<
	$(AR) r $@ $*.o; rm $*.o

all: libisospin.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac
	   
clean:
	rm -f *.a isajet.h softsusy.h spheno.h suspect.h 2hdmc.h higgsbounds.h nmssmtools.h

distclean: 
	rm -f *.a *.o *.x isajet.h softsusy.h spheno.h suspect.h 2hdmc.h higgsbounds.h nmssmtools.h
	
libisospin.a: libisospin.a(softsusy.o) libisospin.a(leshouches.o)\
	  libisospin.a(isospin.o) libisospin.a(isajet.o) libisospin.a(wilson.o)\
	  libisospin.a(bsgamma.o) libisospin.a(alphas.o) libisospin.a(general.o)\
	  libisospin.a(excluded_masses.o) libisospin.a(quark_masses.o) libisospin.a(gmuon.o)\
	  libisospin.a(btaunu.o) libisospin.a(bdtaunu.o) libisospin.a(kmunu.o)\
	  libisospin.a(bsmumu.o) libisospin.a(dslnu.o) libisospin.a(dmunu.o)\
	  libisospin.a(spheno.o) libisospin.a(suspect.o) libisospin.a(2hdmc.o)\
	  libisospin.a(nmssmtools.o) libisospin.a(higgsbounds.o) libisospin.a(flha.o)\
	  libisospin.a(bsll.o) libisospin.a(bsll_extra.o) libisospin.a(bkll.o)\
	  $(RANL)

libisospin.a(softsusy.o): softsusy.c include.h softsusy.h
libisospin.a(leshouches.o): leshouches.c include.h
libisospin.a(isospin.o): isospin.c include.h
libisospin.a(isajet.o): isajet.c include.h isajet.h
libisospin.a(wilson.o): wilson.c include.h
libisospin.a(bsgamma.o): bsgamma.c include.h
libisospin.a(alphas.o): alphas.c include.h
libisospin.a(quark_masses.o): quark_masses.c include.h
libisospin.a(general.o): general.c include.h
libisospin.a(excluded_masses.o): excluded_masses.c include.h
libisospin.a(gmuon.o): gmuon.c include.h
libisospin.a(btaunu.o): btaunu.c include.h
libisospin.a(bdtaunu.o): bdtaunu.c include.h
libisospin.a(kmunu.o): kmunu.c include.h
libisospin.a(bsmumu.o): bsmumu.c include.h
libisospin.a(dslnu.o): dslnu.c include.h
libisospin.a(dmunu.o): dmunu.c include.h
libisospin.a(spheno.o): spheno.c include.h spheno.h
libisospin.a(suspect.o): suspect.c include.h suspect.h
libisospin.a(nmssmtools.o): nmssmtools.c include.h nmssmtools.h
libisospin.a(2hdmc.o): 2hdmc.c include.h 2hdmc.h
libisospin.a(higgsbounds.o): higgsbounds.c higgsbounds.h include.h
libisospin.a(flha.o): flha.c include.h
libisospin.a(bsll.o): bsll.c include.h
libisospin.a(bsll_extra.o): bsll_extra.c include.h
libisospin.a(bkll.o): bkll.c include.h

softsusy.h:
	@echo \#define SOFTSUSY \"$(SOFTSUSY)\" > softsusy.h;

isajet.h:
	@echo \#define ISAJET \"$(ISAJET)\" > isajet.h;

spheno.h:
	@echo \#define SPHENO \"$(SPHENO)\" > spheno.h;

suspect.h:
	@echo \#define SUSPECT \"$(SUSPECT)\" > suspect.h;
	
nmssmtools.h:
	@echo \#define NMSSMTools \"$(NMSSMTools)\" > nmssmtools.h;

2hdmc.h:
	@echo \#define THDMC \"$(THDMC)\" > 2hdmc.h;

higgsbounds.h:
	@echo \#define HBwithFH \"$(HBwithFH)\" > higgsbounds.h;
