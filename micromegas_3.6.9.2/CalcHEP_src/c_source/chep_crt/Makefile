.PHONY: all

lib = ../../lib
bin = ../../bin

all:  $(lib)/serv.a  $(lib)/servNoX11.a

include ../../FlagsForMake

Idir= getmem



override CFLAGS += $(LOCK) -Iinclude $(patsubst %,-I../%/include ,$(Idir)) $(HX11)

OBJ= $(subst .c,.o,$(wildcard *.c)) X11_crt0.o

ifeq (,$(LX11))
 vpath X11_crt0.c noxwin
else 
 vpath X11_crt0.c xwin
endif


$(lib)/serv.a:$(lib)/serv.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)
        
$(lib)/servNoX11.a:$(lib)/serv.a noxwin/X11_crt0.c
	cp $(lib)/serv.a $(lib)/servNoX11.a   	
	$(CC) -c $(CFLAGS) -o noxwin/X11_crt0.o noxwin/X11_crt0.c
	ar d $(lib)/servNoX11.a X11_crt0.o
	ar r $(lib)/servNoX11.a noxwin/X11_crt0.o
	$(RANLIB) $@
	rm noxwin/X11_crt0.o 
	 