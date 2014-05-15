
.PHONY: all clean flags

all:sources/microPath.h
	$(MAKE) -C CalcHEP_src
	$(MAKE) -C sources

sources/microPath.h:
	echo \#define micrO \"$(CURDIR)\"  > sources/microPath.h
        
   
clean:  
	rm -f sources/microPath.h
	./clean

flags: 
	$(MAKE) -C CalcHEP_src flags