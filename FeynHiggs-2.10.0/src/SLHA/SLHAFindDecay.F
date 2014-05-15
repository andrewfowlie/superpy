* SLHAFindDecay.F
* locate an existing decay for the purpose of adding channels to it
* based on J. Guasch's SLHAOldDecay
* last modified 30 Nov 11 th

#include "SLHAIntern.h"


	integer function SLHAFindDecay(sd, parent_id)
	implicit none
	RealType sd(*)
	integer parent_id

	integer decay, next, i
	RealType tmp(LengthDecay)

	decay = 1
	do while( Decay(decay) .ne. invalid )
	  next = Decay_Next(decay)

* move a previous instance of this decay to the end
	  if( Decay_Id(decay) .eq. parent_id ) then
	    if( Decay(decay + next) .ne. invalid ) then
	      do i = 1, next
	        tmp(i) = Decay(decay + i - 1)
	      enddo
	      do while( Decay(decay + next) .ne. invalid )
	        Decay(decay) = Decay(decay + next)
	        decay = decay + 1
	      enddo
	      do i = 1, next
	        Decay(decay + i - 1) = tmp(i)
	      enddo
	    endif
	    SLHAFindDecay = decay
	    return
	  endif

	  decay = decay + next
	enddo

	SLHAFindDecay = 0
	end

