
include ../../FlagsForMake 

override CFLAGS += -Iinclude -I../SLHAplus/include  -I../service2/include  -I../../include -DMAKE=\"$(MAKE)\" -Wmissing-prototypes

special =

#OBJ= $(notdir $(subst .c,.o, $(filter-out $(special),$(wildcard *.c diagram/*.c out/*.c))))

OBJ = dynamic_cs.o kin4.o vp.o VVdecay.o
lib = ../../lib

all: $(lib)/dynamic_me.a $(lib)/dynamic_vp.a

$(lib)/dynamic_me.a:$(lib)/dynamic_me.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)

$(lib)/dynamic_vp.a:$(lib)/dynamic_vp.a(vp_dynam.o)
	$(RANLIB) $@
	@rm -f $(OBJ)



