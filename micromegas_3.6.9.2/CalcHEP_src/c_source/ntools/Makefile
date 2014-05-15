.PHONY: all

include ../../FlagsForMake

Idir = .. chep_crt plot  service2 getmem dynamicME SLHAplus num

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir))

special = 

OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c strfun/*.c))))



lib = $(CALCHEP)/lib


all: $(lib)/ntools.a 

$(lib)/ntools.a:$(lib)/ntools.a($(OBJ))
	$(RANLIB) $(lib)/ntools.a
	@rm -f $(OBJ)

