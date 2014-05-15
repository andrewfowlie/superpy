* SLHADecayTable.F
* fills an array with the decays of a given particle
* last modified 30 Nov 11 th

#include "SLHAIntern.h"


	integer function SLHADecayTable(sd, parent_id,
     &    width, id, maxparticles, maxchannels)
	implicit none
	integer parent_id, maxparticles, maxchannels
	RealType sd(*), width(maxchannels)
	integer id(0:maxparticles,maxchannels)

	integer decay, parent, nextparent, nextchild
	integer nparticles, channel, particle

	channel = 0
	decay = 1

* parent_id = 0 is a wildcard for all decays
	if( parent_id .ne. 0 ) then
	  do while( Decay_Id(decay) .ne. parent_id )
	    if( Decay(decay) .eq. invalid ) goto 1
	    decay = decay + Decay_Next(decay)
	  enddo
	endif

	do while( Decay(decay) .ne. invalid )
	  if( channel .ge. maxchannels ) then
	    print *, "maxchannels too small in SLHADecayTable."
	    goto 1
	  endif
	  channel = channel + 1

	  parent = Decay_Id(decay)
	  nextparent = decay + Decay_Next(decay)
	  width(channel) = Decay(decay + 1)
	  decay = decay + 2

	  if( maxparticles .eq. 0 ) then
	    id(0,channel) = parent
	    goto 2
	  endif

	  id(0,channel) = 1
	  id(1,channel) = parent
	  do particle = 2, maxparticles
	    id(particle,channel) = 0
	  enddo

	  do while( decay .lt. nextparent )
	    nextchild = decay + Decay_Next(decay)
	    nparticles = nextchild - decay

	    if( nparticles .ge. maxparticles ) then
	      print *, "maxparticles too small in SLHADecayTable."
	    else
	      if( channel .ge. maxchannels ) then
	        print *, "maxchannels too small in SLHADecayTable."
	        goto 1
	      endif
	      channel = channel + 1
	      id(0,channel) = nparticles
	      id(1,channel) = parent
	      do particle = 2, nparticles
	        id(particle,channel) = Decay_Id(decay)
	        decay = decay + 1
	      enddo
	      do particle = nparticles + 1, maxparticles
	        id(particle,channel) = 0
	      enddo
	      width(channel) = Decay(decay)
	    endif

	    decay = nextchild
	  enddo

2	  decay = nextparent
	enddo

1	SLHADecayTable = channel
	end

