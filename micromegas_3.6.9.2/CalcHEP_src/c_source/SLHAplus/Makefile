.PHONY: all

include ../../FlagsForMake

ifeq ($(FC),)
all: $(CALCHEP)/lib/libSLHAplus.a 
else
all: $(CALCHEP)/lib/libSLHAplus.a  $(CALCHEP)/lib/faux.o
endif

OBJ = alpha_s.o aux.o chDiag.o ferror.o fixArg.o hgg.o jacobi.o polint.o SLHAreader.o fortran.o

override CFLAGS += -I./include

$(CALCHEP)/lib/libSLHAplus.a:$(CALCHEP)/lib/libSLHAplus.a($(OBJ))
	$(RANLIB) $(CALCHEP)/lib/libSLHAplus.a
	rm -f *.o
	 

 $(CALCHEP)/lib/faux.o:faux.f
	$(FC) $(FFLAGS) -c -o $(CALCHEP)/lib/faux.o faux.f
	