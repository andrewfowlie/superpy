* auxCD.F
* auxillary functions used by the three- and four-point integrals
* these functions are adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 26 Nov 13 th

#include "externals.h"
#include "types.h"


	ComplexType function spence(i_in, z_in, s)
	implicit none
	integer i_in
	ComplexType z_in
	RealType s

#include "lt.h"

	ComplexType Li2series
	external Li2series

	ComplexType z(0:1)
	RealType az1

	z(i_in) = z_in
	z(1-i_in) = 1 - z_in

#ifdef WARNINGS
	if( s .eq. 0 .and. 
     &    Im(z) .eq. 0 .and. abs(Re(z1)) .lt. acc )
     &    print *, "spence: argument on cut"
#endif

	if( Re(z(0)) .lt. .5D0 ) then
	  if( abs(z(0)) .lt. 1 ) then
	    spence = Li2series(z(1), s)
	  else
	    spence = -zeta2 -
     &        .5D0*ln(-z(0), -s)**2 - Li2series(-z(1)/z(0), -s)
	  endif
	else
	  az1 = abs(z(1))
	  if( az1 .lt. 1D-15 ) then
	    spence = zeta2
	  else if( az1 .lt. 1 ) then
	    spence = zeta2 -
     &        ln(z(0), s)*ln(z(1), -s) - Li2series(z(0), -s)
	  else
	    spence = 2*zeta2 +
     &        .5D0*ln(-z(1), -s)**2 - ln(z(0), s)*ln(z(1), -s) +
     &        Li2series(-z(0)/z(1), s)
	  endif
	endif
	end

************************************************************************

	ComplexType function Li2series(z1, s)
	implicit none
	ComplexType z1
	RealType s

	ComplexType xm, x2, new
	integer j

#include "lt.h"

* these are the even-n Bernoulli numbers, already divided by (n + 1)!
* as in Table[BernoulliB[n]/(n + 1)!, {n, 2, 50, 2}]
	RealType b(25)
	data b /
     &    0.02777777777777777777777777777777777777777778774D0, 
     &    -0.000277777777777777777777777777777777777777777778D0, 
     &    4.72411186696900982615268329554043839758125472D-6, 
     &    -9.18577307466196355085243974132863021751910641D-8, 
     &    1.89788699889709990720091730192740293750394761D-9, 
     &    -4.06476164514422552680590938629196667454705711D-11, 
     &    8.92169102045645255521798731675274885151428361D-13, 
     &    -1.993929586072107568723644347793789705630694749D-14, 
     &    4.51898002961991819165047655285559322839681901D-16, 
     &    -1.035651761218124701448341154221865666596091238D-17, 
     &    2.39521862102618674574028374300098038167894899D-19, 
     &    -5.58178587432500933628307450562541990556705462D-21, 
     &    1.309150755418321285812307399186592301749849833D-22, 
     &    -3.087419802426740293242279764866462431595565203D-24, 
     &    7.31597565270220342035790560925214859103339899D-26, 
     &    -1.740845657234000740989055147759702545340841422D-27, 
     &    4.15763564461389971961789962077522667348825413D-29, 
     &    -9.96214848828462210319400670245583884985485196D-31, 
     &    2.394034424896165300521167987893749562934279156D-32, 
     &    -5.76834735536739008429179316187765424407233225D-34, 
     &    1.393179479647007977827886603911548331732410612D-35, 
     &    -3.372121965485089470468473635254930958979742891D-37, 
     &    8.17820877756210262176477721487283426787618937D-39, 
     &    -1.987010831152385925564820669234786567541858996D-40, 
     &    4.83577851804055089628705937311537820769430091D-42 /

	xm = -ln(z1, -s)
	x2 = xm**2
	Li2series = xm - x2/4D0
	do j = 1, 25
	  xm = xm*x2
	  new = Li2series + xm*b(j)
	  if( new .eq. Li2series ) return
	  Li2series = new
	enddo
#ifdef WARNINGS
	print *, "Li2series: bad convergence"
#endif
	end

************************************************************************

	integer function eta(c1, c2, im1x, im2x, im12x)
	implicit none
	ComplexType c1, c2
	RealType im1x, im2x, im12x

	RealType im1, im2, im12

	im1 = Im(c1)
	if( im1 .eq. 0 ) im1 = im1x
	im2 = Im(c2)
	if( im2 .eq. 0 ) im2 = im2x
	im12 = Im(c1*c2)
	if( im12 .eq. 0 ) im12 = im12x

	if( im1 .lt. 0 .and. im2 .lt. 0 .and. im12 .gt. 0 ) then
	  eta = 1
	else
     &  if( im1 .gt. 0 .and. im2 .gt. 0 .and. im12 .lt. 0 ) then
	  eta = -1
	else
	  eta = 0
#ifdef WARNINGS
	  if( .not. (im2 .eq. 0 .and. Re(c2) .gt. 0 .or.
     &               im1 .eq. 0 .and. Re(c1) .gt. 0) .and.
     &        (im1 .eq. 0 .and. Re(c1) .lt. 0 .or.
     &         im2 .eq. 0 .and. Re(c2) .lt. 0 .or.
     &         im12 .eq. 0 .and. Re(c1*c2) .lt. 0) )
     &      print *, "eta not defined"
#endif
	endif
	end

