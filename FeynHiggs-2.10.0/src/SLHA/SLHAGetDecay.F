* SLHAGetDecay.F
* searches a particular decay in the SLHA data structure
* last modified 30 Nov 11 th

#include "SLHAIntern.h"


	RealType function SLHAGetDecay(sd, parent_id,
     &    nchildren, child1_id, child2_id, child3_id, child4_id)
	implicit none
	RealType sd(*)
	integer parent_id
	integer nchildren, child1_id, child2_id, child3_id, child4_id

	integer decay, id(4), i, j, p, nextparent, nextchild

	SLHAGetDecay = 0

	decay = 1
	do while( Decay_Id(decay) .ne. parent_id )
	  if( Decay(decay) .eq. invalid ) return
	  decay = decay + Decay_Next(decay)
	enddo

	if( nchildren .eq. 0 ) goto 1
	if( nchildren .ge. 1 ) id(1) = child1_id
	if( nchildren .ge. 2 ) id(2) = child2_id
	if( nchildren .ge. 3 ) id(3) = child3_id
	if( nchildren .ge. 4 ) id(4) = child4_id

	do i = 1, nchildren
	  p = i
	  do j = i + 1, nchildren
	    if( id(j) .lt. id(p) ) p = j
	  enddo
	  if( p .ne. i ) then
	    j = id(i)
	    id(i) = id(p)
	    id(p) = j
	  endif
	enddo

	nextparent = decay + Decay_Next(decay)
	decay = decay + 2
	do while( decay .lt. nextparent )
	  nextchild = decay + Decay_Next(decay)
	  if( decay + nchildren + 1 .eq. nextchild ) then
	    if( Decay_Id(decay) .eq. id(1) ) then
	      if( decay + 2 .eq. nextchild ) goto 1
	      decay = decay + 1
	      if( Decay_Id(decay) .eq. id(2) ) then
	        if( decay + 2 .eq. nextchild ) goto 1
	        decay = decay + 1
	        if( Decay_Id(decay) .eq. id(3) ) then
	          if( decay + 2 .eq. nextchild ) goto 1
	          decay = decay + 1
	          if( Decay_Id(decay) .eq. id(4) ) then
	            if( decay + 2 .eq. nextchild ) goto 1
	          endif
	        endif
	      endif
	    endif
	  endif
	  decay = nextchild
	enddo

	return

1	SLHAGetDecay = Decay(decay + 1)
	end

