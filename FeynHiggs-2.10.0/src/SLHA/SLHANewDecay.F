* SLHANewDecay.F
* initiates setting of a particular decay in the SLHA data structure
* last modified 30 Nov 11 th

#include "SLHAIntern.h"


	integer function SLHANewDecay(sd, width, parent_id)
	implicit none
	RealType sd(*), width
	integer parent_id

	integer par_id, decay, next

	par_id = abs(parent_id)
	decay = 1
	do while( Decay(decay) .ne. invalid )
	  next = Decay_Next(decay)

* remove a previous instance of this decay, if any
	  if( abs(Decay_Id(decay)) .eq. par_id ) then
	    do while( Decay(decay + next) .ne. invalid )
	      Decay(decay) = Decay(decay + next)
	      decay = decay + 1
	    enddo
	    goto 1
	  endif

	  decay = decay + next
	enddo

1	if( decay + 2 .ge. LengthDecay ) then
	  print *, "Out of memory in SLHANewDecay"
	  SLHANewDecay = 0
	  return
	endif

	Decay(decay) = Decay_Entry(parent_id, 2)
	Decay(decay + 1) = width
	Decay(decay + 2) = invalid
	SLHANewDecay = decay
	end

