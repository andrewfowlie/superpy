# --- variables defined by configure ---

ifneq ($(MAKECMDGOALS),clean)
include ../../../CalcHEP_src/FlagsForMake
endif


all:lGamma.exe


FFLAGS +=   -ffixed-line-length-none

.PHONY: force squared_me_gg squared_me_gz  renconst util

export VPATH := $(CURDIR):$(CURDIR)/renconst

INCLUDE := $(patsubst %,-I%,$(subst :, ,$(VPATH)))

FFLAGS += $(INCLUDE)

ifdef DEBUG
FFLAGS += -g -DDEBUG
endif

export FC
export FFLAGS

LIBTOOLS=$(CURDIR)/../../../Packages/LoopTools-2.1


$(LIBTOOLS)/libooptools.a:
	$(MAKE) -C $(LIBTOOLS)
#subdir = squared_me_gg squared_me_gz renconst  util 
subdir = util
#OBJ= $(notdir $(subst .c,.o, $(filter-out $(special),$(wildcard *.c diagram/*.c out/*.c))))

nnggOBJ = $(subst .F,.o,$(wildcard *.F squared_me_gg/*.F squared_me_gz/*.F renconst/*.F))
nngg.a : $(subdir) nngg.a(run_gg.o run_gz.o mdl_ini.o interface.o interfacegz.o)

nngg.a :  nngg.a($(nnggOBJ))

nngg.a(run.o): 2to2.h renconst.h looptools.h model.h

nngg.a(mdl_ini.o):model.h

lGamma.exe : main.f  interfacegz.o interface.o   nngg.a $(LIBTOOLS)/libooptools.a 
	$(FC) $(FFLAGS) -o lGamma.exe  $^ ../../../CalcHEP_src/lib/libSLHAplus.a   ../../../CalcHEP_src/lib/faux.o


$(subdir) : 
	$(MAKE) -C  $@

clean:	
	rm -f nngg.a lGamma.exe *.o */*.o
