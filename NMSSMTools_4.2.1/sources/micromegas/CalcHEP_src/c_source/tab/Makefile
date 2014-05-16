.PHONY: all

include ../../FlagsForMake

Idir = .. num chep_crt service2 getmem mix_events dynamicME

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir))

lib = ../../lib
bin = ../../bin

#LIBS = $(lib)/num_c.a $(lib)/serv.a $(lib)/dummy.a  $(lib)/sqme_aux.so

#LIBS = $(lib)/serv.a  
LIBS =

all: $(bin)/events2tab $(bin)/show_distr $(bin)/sum_distr $(bin)/lhe2tab
#$(lib)/num_c.a  $(lib)/serv.a $(lib)/dummy.a
$(bin)/events2tab:events2tab.c  dummy_output.c  $(lib)/num_c.a    $(lib)/serv.a ../../include/VandPgate.c   $(lib)/dummy.a $(lib)/ntools.a 
	$(CC) $(CFLAGS)  -o $@ $^ $(LX11) $(lDL)  -lm 

$(bin)/show_distr:show_distr.c  dummy_output.c   $(lib)/num_c.a  ../../include/VandPgate.c    $(lib)/serv.a $(lib)/dummy.a $(lib)/ntools.a
	$(CC) $(CFLAGS)  -o $@ $^ $(LX11) $(lDL)  -lm

$(bin)/sum_distr:sum_distr.c  dummy_output.c  $(lib)/num_c.a  $(lib)/serv.a ../../include/VandPgate.c   $(lib)/dummy.a $(lib)/ntools.a
	$(CC) $(CFLAGS)  -o $@ $^ $(LX11) $(lDL) -lm
	
$(bin)/lhe2tab:lhe2tab.c  ../mix_events/readLHE.c  dummy_output.c   $(lib)/num_c.a  $(lib)/serv.a ../../include/VandPgate.c   $(lib)/dummy.a $(lib)/ntools.a
	$(CC) $(CFLAGS)  -o $@ $^  $(LX11) $(lDL) -lm
