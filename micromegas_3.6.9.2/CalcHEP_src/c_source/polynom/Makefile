
include ../../FlagsForMake

Idir = chep_crt getmem service2 symb

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir)) $(SNUM)


OBJ = $(subst .c,.o,$(wildcard *.c))
lib = ../../lib

$(lib)/symb.a:$(lib)/symb.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)
