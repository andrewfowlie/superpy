.PHONY: all

include ../../FlagsForMake

Idir = .. chep_crt getmem

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir))

special = main_calc.c

OBJ = $(subst .c,.o,$(filter-out $(special),$(wildcard *.c)))

lib = ../../lib
bin = ../../bin

all: $(lib)/serv.a  $(bin)/calc

$(lib)/serv.a:$(lib)/serv.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)

$(bin)/calc: main_calc.c  $(lib)/serv.a
	$(CC) $(CFLAGS)  -o $@    $^ -lm
