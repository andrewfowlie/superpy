.PHONY: all

include ../../FlagsForMake

Idir= .. num chep_crt service2 getmem

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir))

special = view_tab.c

bin = ../../bin
lib = ../../lib

OBJ = $(subst .c,.o,$(filter-out $(special),$(wildcard *.c)))
 
$(bin)/plot_view:  view_tab.c  $(lib)/serv.a
	$(CC) $(CFLAGS) -o $@ $^  $(LX11) -lm

$(lib)/serv.a :$(lib)/serv.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)
