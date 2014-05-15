.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include  ../../CalcHEP_src/FlagsForMake
endif

SRC := $(wildcard *.f)

 
cpsuperh.exe: $(SRC)
	$(FC) $(FFLAGS) -o cpsuperh.exe $(SRC)

clean:
	 rm -f cpsuperh.exe 
