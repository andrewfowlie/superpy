.PHONY: all clean


ifneq ($(MAKECMDGOALS),clean)
include ../../../CalcHEP_src/FlagsForMake
endif 



all: f.out c.out f_.out

f.out: lambdas.F
	$(FC) -o f.out lambdas.F  ../../../CalcHEP_src/lib/libSLHAplus.a

f_.out: lambdas_.F
	$(FC) -o f_.out lambdas_.F  ../../../CalcHEP_src/lib/libSLHAplus.a

c.out: lambdas.c
	$(CC) -o c.out lambdas.c  ../../../CalcHEP_src/lib/libSLHAplus.a -lm
        	
clean: 
	rm -f f.out c.out 