
ROOT=
ifeq (,$(ROOT))
$(error Define path to Root, please)
endif

include ../../FlagsForMake

ifeq (,$(PATH))
  PATH=$(ROOT)/bin 
else    
  PATH:= $(PATH):$(ROOT)/bin
endif 
 

ifeq (,$(LD_LIBRARY_PATH))
   LD_LIBRARY_PATH=$(ROOT)/lib
else    
  LD_LIBRARY_PATH:=$(LD_LIBRARY_PATH):$(ROOT)/lib
endif

ifeq (,$(DYLD_LIBRARY_PATH))
   DYLD_LIBRARY_PATH = $(ROOT)/lib
else     
   DYLD_LIBRARY_PATH:=$(DYLD_LIBRARY_PATH):$(ROOT)/lib
endif    

export PATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH ROOT


lib=$(CALCHEP)/lib
INCL= -I$(ROOT)/include -I$(CALCHEP)/include -I$(CALCHEP)/c_source/dynamicME/include

$(lib)/ch_root.so:ch_dict.cc ch_class.h ch_class.cc VandPgate.o $(lib)/dynamic_vp.a  $(lib)/dynamic_me.a
	$(CXX)  $(SHARED) -o $(lib)/ch_root.so  $(CXXFLAGS)  $(INCL) ch_class.cc ch_dict.cc  $(lib)/dynamic_vp.a  $(lib)/dynamic_me.a \
	$(lib)/num_c.a $(lib)/serv.a $(lib)/sqme_aux.so VandPgate.o  $(lib)/libSLHAplus.a  $(lib)/dummy.a $(lib)/ntools.a \
	-L$(ROOT)/lib -lCore -lCint  $(lDL) -lm $(LX11) 
	rm VandPgate.o

VandPgate.o:$(CALCHEP)/include/VandPgate.c  
	$(CC) $(CFLAGS) -c $(CALCHEP)/include/VandPgate.c
		
ch_dict.cc:ch_class.h
	rootcint -f  ch_dict.cc -c $(INCL)  ch_class.h
	