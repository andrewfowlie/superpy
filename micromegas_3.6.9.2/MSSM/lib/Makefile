
SPHENO   = $(CURDIR)/../../Packages/SPheno3.1.3/bin
SOFTSUSY = $(CURDIR)/../../Packages/softsusy-3.1.6
ISAJET   = $(CURDIR)/../../Packages/Isajet
CERNLIB  = $(HOME)/Packages/cern/pro/lib
#======================
SUSPECT= $(CURDIR)/../../Packages/SuSpect_2.41

#read compiler options ..

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif 

.PHONY: all clean

VPATH = SuSpect:../../sources

SuSpectSRC = $(wildcard SuSpect/*.f)

special = isajetLesH.c

aLibOBJ = $(subst .c,.o,$(filter-out $(special),$(wildcard *.c))) 

ifneq ($(FC),)
   all:  suspect.exe  aLib.a 
else
   all: aLib.a
endif

suspect.exe: $(SuSpectSRC) 
	$(MAKE) -C $(SUSPECT)
		
aLib.a(suspect_call.o):suspect_path.h
aLib.a(isasugra_call.o):isajet_path.h
aLib.a(ssusy_call.o): ssusy_path.h
aLib.a(spheno_call.o):spheno_path.h

ssusy_path.h: Makefile
	echo \#define SOFTSUSY \"$(SOFTSUSY)\"  > ssusy_path.h
spheno_path.h: Makefile
	echo \#define SPHENO \"$(SPHENO)\"  > spheno_path.h
isajet_path.h:Makefile
	echo \#define ISAJET \"`pwd`\"  > isajet_path.h
suspect_path.h:Makefile
	echo \#define SUSPECT \"$(SUSPECT)\"  > suspect_path.h
aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

# for Isajet before 7.79 add  flag  -Dlt779 for FC
isajet_slha: isasugra_slha.F isajetLesH.c
	$(CC) $(CFLAGS)   -c  isajetLesH.c
	$(FC) $(FFLAGS) -o isajet_slha  isasugra_slha.F isajetLesH.o \
   -L$(ISAJET) -L$(CERNLIB) -lisajet -lmathlib 

clean: 
	rm -f *.o *.a  *.so *.dll *.exe *_path.h isajet_slha so_locations 
	$(MAKE) -C nngg     clean