* GetPara.F
* get the derived SUSY parameters from the internal arrays
* this file is part of FeynHiggs
* last modified 31 Jul 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHGetPara(error,
     &    fv_, MSf_, USf_, MASf_, UASf_,
     &    MCha_, UCha_, VCha_, MNeu_, ZNeu_,
     &    Deltab_, MGl_, MHtree_, SAtree_)
	implicit none
	integer error
	integer fv_
	RealType MSf_(2,5,3), MASf_(6,5)
	RealType MCha_(2), MNeu_(4), MGl_
	ComplexType USf_(2,2,5,3), UASf_(6*6,5)
	ComplexType UCha_(2*2), VCha_(2*2)
	ComplexType ZNeu_(4*4), Deltab_
	RealType MHtree_(4), SAtree_

#include "FH.h"
#define __SUBROUTINE__ "FHGetPara"

	integer g, i, t, ty

	error = 0
	CheckSf(error)
	if( error .gt. 0 ) return

	fv_ = fv

	do t = 1, 5
	  ty = t + isB(t - 1)*(bBR - 5)
	  do g = 1, 3
	    MSf_(1,t,g) = MSf(1,ty,g)
	    MSf_(2,t,g) = MSf(2,ty,g)
	    USf_(1,1,t,g) = USf(1,1,ty,g)
	    USf_(2,1,t,g) = USf(2,1,ty,g)
	    USf_(1,2,t,g) = USf(1,2,ty,g)
	    USf_(2,2,t,g) = USf(2,2,ty,g)
	  enddo
	  do i = 1, 6
	    MASf_(i,t) = MASf(i,ty)
	  enddo
	  do i = 1, 6*6
	    UASf_(i,t) = UASf_flat(i,ty)
	  enddo
	enddo

	do i = 1, 2
	  MCha_(i) = MCha(i)
	enddo
	do i = 1, 2*2
	  UCha_(i) = UCha_flat(i)
	enddo
	do i = 1, 2*2
	  VCha_(i) = VCha_flat(i)
	enddo

	do i = 1, 4
	  MNeu_(i) = MNeu(i)
	enddo
	do i = 1, 4*4
	  ZNeu_(i) = ZNeu_flat(i)
	enddo

	Deltab_ = Deltaf(bTR,3)
	MGl_ = MGl

	MHtree_(1) = MHtree(1)
	MHtree_(2) = MHtree(2)
	MHtree_(3) = MHtree(3)
	MHtree_(4) = MHtree(4)
	SAtree_ = SA
	end

