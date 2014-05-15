
.PHONY:  clean all

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

all: work_aux.a 

workOBJ =  VandP.o  path.o
override CFLAGS += -I../../sources/micromegas/CalcHEP_src/include

work_aux.a: work_aux.a($(workOBJ))
	$(RANLIB) work_aux.a 
	rm -f $(workOBJ)

VandP.c: models/vars1.mdl models/func1.mdl models/prtcls1.mdl
	../../CalcHEP_src/bin/make_VandP models 1 

path.c:
	echo char \* WORK=\"`pwd`\"\; >  path.c


clean:
	rm -f *.a *.o  lock_ LOCK tmp/*  so_generated/* lock_
	rm -f  path.c VandP.c  autoprot.h
	rm -rf _*_*_  results/*
