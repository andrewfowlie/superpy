* SLHAAddDecay.F
* adds a channel to a parent node in the SLHA data structure
* last modified 30 Nov 11 th

#include "SLHAIntern.h"


	subroutine SLHAAddDecay(sd, width, parent,
     &    nchildren, child1_id, child2_id, child3_id, child4_id)
	implicit none
	RealType sd(*), width
	integer parent
	integer nchildren, child1_id, child2_id, child3_id, child4_id

	integer child, id(4), i, j, p

	if( parent .le. 0 .or. width .eq. 0 ) return

	child = parent + Decay_Next(parent)

	if( child + nchildren + 1 .ge. LengthDecay ) then
	  print *, "Out of memory in SLHAAddDecay"
	  parent = 0
	  return
	endif

	if( nchildren .ge. 1 ) id(1) = child1_id
	if( nchildren .ge. 2 ) id(2) = child2_id
	if( nchildren .ge. 3 ) id(3) = child3_id
	if( nchildren .ge. 4 ) id(4) = child4_id

	do i = 1, nchildren
	  p = i
	  do j = i + 1, nchildren
	    if( id(j) .lt. id(p) ) p = j
	  enddo
	  Decay(child) = Decay_Entry(id(p), nchildren + 1)
	  child = child + 1
	  id(p) = id(i)
	enddo

	Decay(child) = width
	Decay(child+1) = invalid
	Decay(parent) = Decay_Entry(Decay_Id(parent), child - parent + 1)
	end

