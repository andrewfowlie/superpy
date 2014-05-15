* inline.h
* inline versions of the util functions
* this file is part of FormCalc
* last modified 9 Mar 13 th


#ifndef INLINE_H
#define INLINE_H
#else

	RealType SInvariant, TInvariant
	ComplexType Pair, Eps_, Eps
	ComplexType SxS
	ComplexType SxV1, SxV2, SxB1, SxB2
	ComplexType VxS1, VxS2, BxS1, BxS2
	ComplexType SxVxB1, SxVxB2, SxBxV1, SxBxV2
	ComplexType BxVxS1, BxVxS2, VxBxS1, VxBxS2
	ComplexType SxVxBxV1, SxVxBxV2, SxBxVxB1, SxBxVxB2
	ComplexType VxBxVxS1, VxBxVxS2, BxVxBxS1, BxVxBxS2
	ComplexType ChainV0, ChainB0
	ComplexType ChainV1, ChainB1
	ComplexType ChainV2, ChainB2
	ComplexType ChainV3, ChainB3
	ComplexType ChainV4, ChainB4
	ComplexType ChainV5, ChainB5
	ComplexType ChainV6, ChainB6
	integer IndexSign, IndexEps
	RealType SqDiff, ThreeMom

	ComplexType l1_, l2_, r1_, r2_
	integer a_, b_, c_, d_, e_, f_
	integer iL_, iR_, eL_, eR_
	RealType sqrtS_, ma_, mb_

	SInvariant(a_, b_) =
     &    (Re(vec(1,1,a_)) + Re(vec(1,1,b_)))*
     &    (Re(vec(2,2,a_)) + Re(vec(2,2,b_))) -
     &    Sq(vec(1,2,a_) + vec(1,2,b_))

	TInvariant(a_, b_) =
     &    (Re(vec(1,1,a_)) - Re(vec(1,1,b_)))*
     &    (Re(vec(2,2,a_)) - Re(vec(2,2,b_))) -
     &    Sq(vec(1,2,a_) - vec(1,2,b_))

	Pair(a_, b_) = .5D0*(
     &    vec(1,1,a_)*vec(2,2,b_) + vec(2,2,a_)*vec(1,1,b_) -
     &    vec(1,2,a_)*vec(2,1,b_) - vec(2,1,a_)*vec(1,2,b_) )

	Eps_(a_, b_, c_, d_) =
     &    (vec(1,1,a_)*vec(2,2,b_) - vec(2,2,a_)*vec(1,1,b_))*
     &    (vec(2,1,c_)*vec(1,2,d_) - vec(1,2,c_)*vec(2,1,d_))
	Eps(a_, b_, c_, d_) = .25D0*(
     &    Eps_(a_, b_, c_, d_) + Eps_(c_, d_, a_, b_) -
     &    Eps_(a_, c_, b_, d_) - Eps_(b_, d_, a_, c_) +
     &    Eps_(a_, d_, b_, c_) + Eps_(b_, c_, a_, d_) )

	SxS(l1_,l2_, r1_,r2_) = l1_*r1_ + l2_*r2_

	SxV1(l1_,l2_, a_) = l1_*vec(1,1,a_) + l2_*vec(2,1,a_)
	SxV2(l1_,l2_, a_) = l2_*vec(2,2,a_) + l1_*vec(1,2,a_)

	SxB1(l1_,l2_, a_) = l1_*vec(2,2,a_) - l2_*vec(2,1,a_)
	SxB2(l1_,l2_, a_) = l2_*vec(1,1,a_) - l1_*vec(1,2,a_)

	VxS1(a_, r1_,r2_) = vec(1,1,a_)*r1_ + vec(1,2,a_)*r2_
	VxS2(a_, r1_,r2_) = vec(2,1,a_)*r1_ + vec(2,2,a_)*r2_

	BxS1(a_, r1_,r2_) = vec(2,2,a_)*r1_ - vec(1,2,a_)*r2_
	BxS2(a_, r1_,r2_) = vec(1,1,a_)*r2_ - vec(2,1,a_)*r1_

	SxVxB1(l1_,l2_, a_, b_) =
     &    SxB1(SxV1(l1_,l2_, a_),SxV2(l1_,l2_, a_), b_)
	SxVxB2(l1_,l2_, a_, b_) =
     &    SxB2(SxV1(l1_,l2_, a_),SxV2(l1_,l2_, a_), b_)

	SxBxV1(l1_,l2_, a_, b_) =
     &    SxV1(SxB1(l1_,l2_, a_),SxB2(l1_,l2_, a_), b_)
	SxBxV2(l1_,l2_, a_, b_) =
     &    SxV2(SxB1(l1_,l2_, a_),SxB2(l1_,l2_, a_), b_)

	BxVxS1(b_, a_, r1_,r2_) =
     &    BxS1(b_, VxS1(a_, r1_,r2_),VxS2(a_, r1_,r2_))
	BxVxS2(b_, a_, r1_,r2_) =
     &    BxS2(b_, VxS1(a_, r1_,r2_),VxS2(a_, r1_,r2_))

	VxBxS1(b_, a_, r1_,r2_) =
     &    VxS1(b_, BxS1(a_, r1_,r2_),BxS2(a_, r1_,r2_))
	VxBxS2(b_, a_, r1_,r2_) =
     &    VxS2(b_, BxS1(a_, r1_,r2_),BxS2(a_, r1_,r2_))

	SxVxBxV1(l1_,l2_, a_, b_, c_) =
     &    SxBxV1(SxV1(l1_,l2_, a_),SxV2(l1_,l2_, a_), b_, c_)
	SxVxBxV2(l1_,l2_, a_, b_, c_) =
     &    SxBxV2(SxV1(l1_,l2_, a_),SxV2(l1_,l2_, a_), b_, c_)

	SxBxVxB1(l1_,l2_, a_, b_, c_) =
     &    SxVxB1(SxB1(l1_,l2_, a_),SxB2(l1_,l2_, a_), b_, c_)
	SxBxVxB2(l1_,l2_, a_, b_, c_) =
     &    SxVxB2(SxB1(l1_,l2_, a_),SxB2(l1_,l2_, a_), b_, c_)

	VxBxVxS1(c_, b_, a_, r1_,r2_) =
     &    VxBxS1(c_, b_, VxS1(a_, r1_,r2_),VxS2(a_, r1_,r2_))
	VxBxVxS2(c_, b_, a_, r1_,r2_) =
     &    VxBxS2(c_, b_, VxS1(a_, r1_,r2_),VxS2(a_, r1_,r2_))

	BxVxBxS1(c_, b_, a_, r1_,r2_) =
     &    BxVxS1(c_, b_, BxS1(a_, r1_,r2_),BxS2(a_, r1_,r2_))
	BxVxBxS2(c_, b_, a_, r1_,r2_) =
     &    BxVxS2(c_, b_, BxS1(a_, r1_,r2_),BxS2(a_, r1_,r2_))

#ifndef SpiLV
#define SpiLV(iL,eL) (1-2*eL)*vec(1+eL,1+eL,iL), vec(2-eL,1+eL,iL)
#define SpiLB(iL,eL) (1-2*eL)*vec(1+eL,2-eL,iL), vec(2-eL,2-eL,iL)
#define SpiRV(eR,iR) vec(1+eR,1+eR,iR), (1-2*eR)*vec(2-eR,1+eR,iR)
#define SpiRB(eR,iR) vec(1+eR,2-eR,iR), (1-2*eR)*vec(2-eR,2-eR,iR)
#endif

	ChainV0(iL_,eL_, eR_,iR_) = SxS(
     &    SpiLB(iL_,eL_),
     &    SpiRV(eR_,iR_) )
	ChainB0(iL_,eL_, eR_,iR_) = SxS(
     &    SpiLV(iL_,eL_),
     &    SpiRB(eR_,iR_) )

	ChainV1(iL_,eL_, a_, eR_,iR_) = SxS(
     &    SxV1(SpiLB(iL_,eL_), a_),
     &    SxV2(SpiLB(iL_,eL_), a_),
     &    SpiRB(eR_,iR_) )
	ChainB1(iL_,eL_, a_, eR_,iR_) = SxS(
     &    SxB1(SpiLV(iL_,eL_), a_),
     &    SxB2(SpiLV(iL_,eL_), a_),
     &    SpiRV(eR_,iR_) )

	ChainV2(iL_,eL_, a_, b_, eR_,iR_) = SxS(
     &    SxV1(SpiLB(iL_,eL_), a_),
     &    SxV2(SpiLB(iL_,eL_), a_),
     &    BxS1(b_, SpiRV(eR_,iR_)),
     &    BxS2(b_, SpiRV(eR_,iR_)) )
	ChainB2(iL_,eL_, a_, b_, eR_,iR_) = SxS(
     &    SxB1(SpiLV(iL_,eL_), a_),
     &    SxB2(SpiLV(iL_,eL_), a_),
     &    VxS1(b_, SpiRB(eR_,iR_)),
     &    VxS2(b_, SpiRB(eR_,iR_)) )

	ChainV3(iL_,eL_, a_, b_, c_, eR_,iR_) = SxS(
     &    SxVxB1(SpiLB(iL_,eL_), a_, b_),
     &    SxVxB2(SpiLB(iL_,eL_), a_, b_),
     &    VxS1(c_, SpiRB(eR_,iR_)),
     &    VxS2(c_, SpiRB(eR_,iR_)) )
	ChainB3(iL_,eL_, a_, b_, c_, eR_,iR_) = SxS(
     &    SxBxV1(SpiLV(iL_,eL_), a_, b_),
     &    SxBxV2(SpiLV(iL_,eL_), a_, b_),
     &    BxS1(c_, SpiRV(eR_,iR_)),
     &    BxS2(c_, SpiRV(eR_,iR_)) )

	ChainV4(iL_,eL_, a_, b_, c_, d_, eR_,iR_) = SxS(
     &    SxVxB1(SpiLB(iL_,eL_), a_, b_),
     &    SxVxB2(SpiLB(iL_,eL_), a_, b_),
     &    VxBxS1(c_, d_, SpiRV(eR_,iR_)),
     &    VxBxS2(c_, d_, SpiRV(eR_,iR_)) )
	ChainB4(iL_,eL_, a_, b_, c_, d_, eR_,iR_) = SxS(
     &    SxBxV1(SpiLV(iL_,eL_), a_, b_),
     &    SxBxV2(SpiLV(iL_,eL_), a_, b_),
     &    BxVxS1(c_, d_, SpiRB(eR_,iR_)),
     &    BxVxS2(c_, d_, SpiRB(eR_,iR_)) )

	ChainV5(iL_,eL_, a_, b_, c_, d_, e_, eR_,iR_) = SxS(
     &    SxVxBxV1(SpiLB(iL_,eL_), a_, b_, c_),
     &    SxVxBxV2(SpiLB(iL_,eL_), a_, b_, c_),
     &    BxVxS1(d_, e_, SpiRB(eR_,iR_)),
     &    BxVxS2(d_, e_, SpiRB(eR_,iR_)) )
	ChainB5(iL_,eL_, a_, b_, c_, d_, e_, eR_,iR_) = SxS(
     &    SxBxVxB1(SpiLV(iL_,eL_), a_, b_, c_),
     &    SxBxVxB2(SpiLV(iL_,eL_), a_, b_, c_),
     &    VxBxS1(d_, e_, SpiRV(eR_,iR_)),
     &    VxBxS2(d_, e_, SpiRV(eR_,iR_)) )

	ChainV6(iL_,eL_, a_, b_, c_, d_, e_, f_, eR_,iR_) = SxS(
     &    SxVxBxV1(SpiLB(iL_,eL_), a_, b_, c_),
     &    SxVxBxV2(SpiLB(iL_,eL_), a_, b_, c_),
     &    BxVxBxS1(d_, e_, f_, SpiRV(eR_,iR_)),
     &    BxVxBxS2(d_, e_, f_, SpiRV(eR_,iR_)) )
	ChainB6(iL_,eL_, a_, b_, c_, d_, e_, f_, eR_,iR_) = SxS(
     &    SxBxVxB1(SpiLV(iL_,eL_), a_, b_, c_),
     &    SxBxVxB2(SpiLV(iL_,eL_), a_, b_, c_),
     &    VxBxVxS1(d_, e_, f_, SpiRB(eR_,iR_)),
     &    VxBxVxS2(d_, e_, f_, SpiRB(eR_,iR_)) )

	IndexSign(a_) = signbit(ior(a_, -a_)) - 2*signbit(a_)
	IndexEps(a_, b_, c_) =
     &    IndexSign(a_ - b_)*IndexSign(c_ - b_)*IndexSign(a_ - c_)

	SqDiff(ma_, mb_) = (ma_ - mb_)*(ma_ + mb_)
	ThreeMom(sqrtS_, ma_, mb_) = sqrt(SqDiff(
     &    .5D0*(sqrtS_ - SqDiff(ma_, mb_)/sqrtS_), mb_ ))

#endif

