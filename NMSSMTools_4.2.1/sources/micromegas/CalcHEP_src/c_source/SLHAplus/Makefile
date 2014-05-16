.PHONY: all

include ../../FlagsForMake

all: $(CALCHEP)/lib/libSLHAplus.a  $(CALCHEP)/lib/model_aux.a

ifeq ($(FC),)
OBJ = alpha_s.o aux.o chDiag.o ferror.o fixArg.o hgg.o jacobi.o polint.o SLHAreader.o
else 
OBJ = alpha_s.o aux.o chDiag.o ferror.o fixArg.o hgg.o jacobi.o polint.o SLHAreader.o fortran.o faux.o
endif

override CFLAGS += -I./include

$(CALCHEP)/lib/libSLHAplus.a:$(CALCHEP)/lib/libSLHAplus.a($(OBJ))
	$(RANLIB) $(CALCHEP)/lib/libSLHAplus.a
	rm -f *.o
 
$(CALCHEP)/lib/model_aux.a:$(CALCHEP)/lib/libSLHAplus.a
	ln -s $(CALCHEP)/lib/libSLHAplus.a $(CALCHEP)/lib/model_aux.a 
