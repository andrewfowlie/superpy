.PHONY: all

include ../../FlagsForMake

override CFLAGS += -Iinclude -I../chep_crt/include

OBJ= $(subst .c,.o,$(wildcard *.c))

lib = ../../lib

all:$(lib)/serv.a

$(lib)/serv.a:$(lib)/serv.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)
