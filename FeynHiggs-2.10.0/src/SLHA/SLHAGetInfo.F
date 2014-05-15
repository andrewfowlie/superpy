* SLHAGetInfo.F
* retrieve informational messages from block SP/DCINFO
* last modified 30 Nov 11 th

#include "SLHAIntern.h"

* the following two definitions must match the ones in SLHA2.m
#define MAXMSG 15
#define MSGLEN 5

* codes are, according to the Accord,
*   1 = program name
*   2 = program version
*   3 = warning messages
*   4 = error messages


	subroutine SLHAGetInfo(slhablock, i, text)
	implicit none
	ComplexType slhablock(*)
	integer i
	character*(*) text

	integer n

	n = int(slhablock(1))
	if( n .eq. invalid .or. i .le. 0 .or. i .gt. n ) then
	  print *, "Bad message number"
	  return
	endif

	call SLHACopyString(
     &    SPInfo_Len, slhablock((i-1)*MSGLEN + MAXMSG + 3),
     &    len(text), text )
	end

