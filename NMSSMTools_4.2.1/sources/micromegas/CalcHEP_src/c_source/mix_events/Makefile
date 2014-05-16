.PHONY: all

include ../../FlagsForMake

lib = ../../lib
bin = ../../bin

#CERN=/usr/lib64/cernlib/2006-g77/lib/
#CERN=/usr/lib64/cernlib/2006/lib/
ifeq (,$(CERN))
   nt_maker =
else
ifeq (,$(FC))
   nt_maker =
else    
   nt_maker = $(bin)/nt_maker
endif 
endif

all: $(bin)/event_mixer  $(nt_maker) $(bin)/event2lhe 

INC = -Iinclude -I../SLHAplus/include  -I../symb/include  -I../polynom/include 

$(bin)/event_mixer:event_mixer.c event2pyth.c ch_events.c   $(lib)/libSLHAplus.a
	$(CC) $(CFLAGS) $(INC) -o $@ $^ -lm 

$(bin)/event2lhe:event2lhe.c   ch_events.c 
	$(CC) $(CFLAGS) $(INC) -o $@ $^ -lm 


$(bin)/nt_maker:nt_maker_02.F readLHE.c 
	- $(FC) -static  -Iinclude   -o $@ $^  -L$(CERN) -lpacklib -lkernlib
	


