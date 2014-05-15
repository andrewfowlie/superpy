* fortranflush.F
* C-callable function to flush the Fortran stdout (unit 6)
* this file is part of FeynHiggs
* last modified 20 Mar 09 th


	subroutine fortranflush()
	implicit none

	call flush(6)
	end

