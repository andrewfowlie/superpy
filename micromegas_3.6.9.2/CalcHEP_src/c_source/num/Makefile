.PHONY: all

include ../../FlagsForMake

Idir = .. chep_crt plot  service2 getmem dynamicME SLHAplus ntools

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir))

special = n_calchep.c usrfun_dummy.c pdf_dummy.c usrFF_dummy.c userInt.c ferror.c sqme_aux.c

OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c strfun/*.c))))

VPATH = strfun

lib = $(CALCHEP)/lib
bin = ../../bin

all: $(lib)/num_c.a $(lib)/dummy.a $(lib)/n_calchep.o  $(lib)/sqme_aux.$(SO)   $(bin)/Int

$(lib)/num_c.a:$(lib)/num_c.a($(OBJ))
	$(RANLIB) $(lib)/num_c.a
	@rm -f $(OBJ)

$(lib)/dummy.a:$(lib)/dummy.a(usrfun_dummy.o  usrFF_dummy.o pdf_dummy.o)
	$(RANLIB) $@
	@rm -f usrfun_dummy.o pdf_dummy.o usrFF_dummy.o

 $(lib)/n_calchep.o: n_calchep.c
	$(CC) $(CFLAGS) -c -o $@  $<

$(bin)/Int:  userInt.c $(lib)/serv.a $(lib)/num_c.a $(lib)/ntools.a
	$(CC) -o $@  $(CFLAGS) $^  -lm

ifneq ($(SONAME),)  
  soname_sqme_aux := $(SONAME)$(CALCHEP)/lib/sqme_aux.$(SO)
else
   soname_sqme_aux :=  
endif


$(lib)/sqme_aux.$(SO): sqme_aux.c ferror.c
	$(CC) $(CFLAGS) -DforMICROMEGAS $(SHARED) -o $@  $(soname_sqme_aux)   $^ -lm

