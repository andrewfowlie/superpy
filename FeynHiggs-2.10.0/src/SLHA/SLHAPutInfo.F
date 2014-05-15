* SLHAPutInfo.F
* record informational messages for write-out to block SP/DCINFO
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


	subroutine SLHAPutInfo(slhablock, code, text)
	implicit none
	ComplexType slhablock(*)
	integer code
	character*(*) text

	integer n

	n = int(slhablock(1))
	if( n .eq. invalid ) then
	  n = 0
	else if( n .ge. MAXMSG ) then
	  print *, "No space left for messages"
	  return
	endif

	slhablock(1) = n + 1
	slhablock(n + 3) = code
	slhablock(2) = max(int(slhablock(2)), code)
	call SLHACopyString(
     &    len(text), text,
     &    SPInfo_Len, slhablock(n*MSGLEN + MAXMSG + 3) )
	end

