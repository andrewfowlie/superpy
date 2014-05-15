* Delta.f
* the Kronecker delta
* this file is part of FormCalc
* last modified 15 Jun 04 th


	integer function Delta(i, j)
	implicit none
	integer i, j

	if( i .eq. j ) then
	  Delta = 1
	else
	  Delta = 0
	endif
	end

