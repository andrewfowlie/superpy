
.PHONY: all  COMPILE clean

ifeq ($(MAKECMDGOALS),clean)
clean :
	@if(test -r FlagsForSh) then echo FlagsForSh - compiler flags; rm -i  FlagsForSh; fi
	./sbin/setPath " "
	rm -f  include/rootDir.h
	rm -f  FlagsForMake  so_locations CMessage
	cd lib; rm -rf *.*
	cd bin; rm -f Int calc events2tab s_calchep show_distr sum_distr plot_view make_VandP event_info event_mixer lhe2tab nt_maker showHelp *.exe
	rm -f c_source/*/*.o c_source/*/so_location
	rm -fr work/results/* work/tmp/*  work/Processes work/Events work/html
	rm -f c_source/Root/ch_dict.cc c_source/Root/ch_dict.h c_source/Root/*.o
	chmod 644 mkWORKdir
endif

all:FlagsForMake COMPILE


flags:
	./getFlags

FlagsForMake: flags

COMPILE:FlagsForMake
	./sbin/setPath $(CURDIR)
	$(MAKE) -C c_source/getmem 
	$(MAKE) -C c_source/chep_crt
	$(MAKE) -C c_source/service2
	$(MAKE) -C c_source/polynom  
	$(MAKE) -C c_source/ntools
	$(MAKE) -C c_source/dynamicME
	$(MAKE) -C c_source/plot 
	$(MAKE) -C c_source/num 
	$(MAKE) -C c_source/tab
	$(MAKE) -C c_source/SLHAplus
	$(MAKE) -C c_source/mix_events
	$(MAKE) -C c_source/symb
	@if( test ! -d work/bin) then ln -s  `pwd`/bin  `pwd`/work/bin; fi
	chmod 755 mkWORKdir
	#--------------------------------------------------------
	@if(test -z "`grep lX11 FlagsForMake`") then  echo "# Warning:";echo "#    CalcHEP has compiled without X11."; echo "#    This is not fatal, but it means that CalcHEP can only"; echo "#    be used in blind mode.  Please see the CalcHEP manual"; echo "#    for further information."; fi
	# CalcHEP has compiled successfuly and can be started.
	# The manual can be found on the CalcHEP website:
	#      http://theory.sinp.msu.ru/~pukhov/calchep.html
	# The next step is typically to run 
	#      ./mkWORKdir  <new_dir>
	# where <new_dir> is the new directory where you will do
	# your calculations.  After creating this directory, you
	# should cd into it and run calchep or calchep_batch.
	# Please see the manual for further details.
	#---------------------------------------------------------
