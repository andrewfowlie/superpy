ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif 

.PHONY: all clean

suspectOBJ = bsg.o suspect2.o  twoloophiggs.o

all:  suspect.exe 
suspect.exe:  suspect2_call.f suspect.a
	$(FC) -o suspect.exe $(FFLAGS) $^

suspect.a:suspect.a($(suspectOBJ))
	$(RANLIB) suspect.a
	rm -f $(suspectOBJ)

clean:
	rm -f suspect.a suspect.exe
