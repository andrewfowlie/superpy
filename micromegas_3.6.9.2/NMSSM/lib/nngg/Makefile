# --- variables defined by configure ---

.PHONY: clean all

CHEP   =  ../../../CalcHEP_src
LIBTOOLS = ../../../Packages/LoopTools-2.1

ifneq ($(MAKECMDGOALS),clean)
include   $(CHEP)/FlagsForMake
endif 

FFLAGS += -g  -ffixed-line-length-none 


.PHONY: force squared_me_gg squared_me_gz  renconst util

export VPATH := $(CURDIR):$(CURDIR)/renconst

INCLUDE := $(patsubst %,-I%,$(subst :, ,$(VPATH)))

FFLAGS += $(INCLUDE)

ifdef DEBUG
FFLAGS += -g -DDEBUG
endif

export FC
export FFLAGS

all:lGamma.exe



$(LIBTOOLS)/libooptools.a:
	$(MAKE) -C $(LIBTOOLS)
#subdir = squared_me_gg squared_me_gz renconst  util 
subdir = util
#OBJ= $(notdir $(subst .c,.o, $(filter-out $(special),$(wildcard *.c diagram/*.c out/*.c))))

nnggOBJ = $(subst .F,.o,$(wildcard *.F squared_me_gg/*.F squared_me_gz/*.F renconst/*.F)) interfacegz.o interface.o nmssmf.o
nngg.a : $(subdir) nngg.a(run_gg.o run_gz.o mdl_ini.o interface.o interfacegz.o)

nngg.a : $(subdir)  nngg.a($(nnggOBJ))

nngg.a(run.o): 2to2.h renconst.h looptools.h model.h

nngg.a(mdl_ini.o):model.h

lGamma.exe : main.f   nngg.a $(LIBTOOLS)/libooptools.a 
	$(FC) $(FFLAGS) -o lGamma.exe  $^ $(CHEP)/lib/libSLHAplus.a  $(CHEP)/lib/faux.o


$(subdir) : 
	$(MAKE) -C  $@

clean:	
	rm -f nngg.a lGamma.exe *.o */*.o
