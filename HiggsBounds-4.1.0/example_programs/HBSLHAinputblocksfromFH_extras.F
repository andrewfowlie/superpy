c-------------------------------------------------------------------
c for use with FeynHiggs 2.8.0
c-------------------------------------------------------------------
	subroutine createSLHAfilewithFHwithoutHBinputblocks(
     &          inputfile,outputfile, 
     &          g2hjbb_s,g2hjbb_p,g2hjtoptop_s,g2hjtoptop_p,       
     &          g2hjtautau_s,g2hjtautau_p,                         
     &          g2hjWW,g2hjZZ,                                     
     &          g2hjgg,g2hjggZ,g2hjhiZ_nHbynH                      )

        implicit none

	integer error

c used by FHCouplings:
#include "FHCouplings.h"
	double complex couplings(ncouplings), couplingsms(ncouplingsms)
	double precision gammas(ngammas), gammasms(ngammasms)

#include "SLHA.h"

#include "PDG.h"
        integer parent

c used by FHHiggsCorr
	double precision MHiggs(4)
	double complex SAeff, DeltaSAeff, UHiggs(3,3), ZHiggs(3,3)

c used by FHSetSLHA,FHOutputSLHA
	integer key
	double complex slhadata(nslhadata)
	character*(*) inputfile,outputfile

c used by FHRetrieveSMPara:
       double precision invAlfa, AlfasMZ, GF
       double precision ME, MU, MD, MM, MC, MS, ML, MB
       double precision MW, MZ
       double precision CKMlambda, CKMA, CKMrhobar, CKMetabar

        double precision
     &          g2hjWW(3),g2hjZZ(3),                
     &          g2hjgg(3),g2hjggZ(3),g2hjhiZ_nHbynH(3,3)
        double precision g2hjbb_s(3),g2hjbb_p(3)
        double precision g2hjtautau_s(3),g2hjtautau_p(3)
        double precision g2hjtoptop_s(3),g2hjtoptop_p(3)

c        call FHSetFlagsString(error, "400302113")
        call FHSetFlagsString(error, "400202110")
	if( error .ne. 0 ) return

	call SLHARead(error, slhadata, inputfile, 1)
	if( error .ne. 0 ) return

	call FHSetSLHA(error, slhadata)
	if( error .ne. 0 ) return

        call FHSetDebug(0)

	call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)
	if( error .ne. 0 ) return
	
    	call FHSelectUZ(error, 2, 2, 1)
    	
	if( error .ne. 0 ) return

	call FHCouplings(error,
     &    couplings, couplingsms, gammas, gammasms, 1)
	if( error .ne. 0 ) return

c turn all of them on i.e. 255=1+2+4+8+16+32+64+128
        key=255 

	call FHOutputSLHA(error, slhadata, key)
	if( error .ne. 0 ) return

	call SLHAWrite(error, slhadata, outputfile)
	if( error .ne. 0 ) return

        call couplingsfromFH
	if( error .ne. 0 ) return

c---------------------------------------------------------
c---------------------------------------------------------
        contains
c-------------------------------------------------------------------
c-------------------------------------------------------------------
	subroutine couplingsfromFH

c used by FHHiggsProd:
        double precision sqrts, prodxs(nprodxs)

c misc:
        integer i,j
        double precision norm,CW2,Pi

        Pi = 3.1415926535897932384626433832795029D0

c We would like FH to calculate LHC cross sections
        sqrts=14.0D0
        call FHHiggsProd(error, sqrts, prodxs)
	if( error .ne. 0 ) return

        call FHRetrieveSMPara(error,
     &   invAlfa, AlfasMZ, GF,
     &   ME, MU, MD, MM, MC, MS, ML, MB,
     &   MW, MZ,
     &   CKMlambda, CKMA, CKMrhobar, CKMetabar)
	if( error .ne. 0 ) return
 
        do i=1,3
   
         g2hjbb_s(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &                   /RCouplingSM(H0FF(i,4,3,3))+ 
     &                    LCoupling(H0FF(i,4,3,3))
     &                   /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0
         g2hjbb_p(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &                   /RCouplingSM(H0FF(i,4,3,3))- 
     &                    LCoupling(H0FF(i,4,3,3))
     &                   /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0

         g2hjtautau_s(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &                   /RCouplingSM(H0FF(i,2,3,3))+ 
     &                    LCoupling(H0FF(i,2,3,3))
     &                   /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0
         g2hjtautau_p(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &                   /RCouplingSM(H0FF(i,2,3,3))- 
     &                    LCoupling(H0FF(i,2,3,3))
     &                   /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0

         g2hjtoptop_s(i)=(abs(RCoupling(H0FF(i,3,3,3))
     &                   /RCouplingSM(H0FF(i,3,3,3))+ 
     &                    LCoupling(H0FF(i,3,3,3))
     &                   /LCouplingSM(H0FF(i,3,3,3)))/2.0D0)**2.0D0
         g2hjtoptop_p(i)=(abs(RCoupling(H0FF(i,3,3,3))
     &                   /RCouplingSM(H0FF(i,3,3,3))- 
     &                    LCoupling(H0FF(i,3,3,3))
     &                   /LCouplingSM(H0FF(i,3,3,3)))/2.0D0)**2.0D0

         g2hjWW(i)= dble(  Coupling(H0VV(i,4)) 
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,4)) 
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjWW(i)= WhTev(i)/WhTevSM(i
c	  g2hjWW(i)= qqhTev(i)/qqhTevSM(i)
c as long as MH>80 GeV and uzint=uzext

         g2hjZZ(i)= dble(  Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjZZ(i)= ZhTev(i)/ZhTevSM(i)
c as long as MH>80 GeV and uzint=uzext
c It is also equivalent to g2hjWW(i)

         g2hjggZ(i)=0.0D0

c tevatron gluon fusion XS is not calculated in FH is MH<=90 geV
         if(MHiggs(i).gt.90.0001D0)then      
           if(gghSM(i).gt.0.0D0)then
            g2hjgg(i) = ggh(i)/gghSM(i)
           else
            g2hjgg(i) = 0.0D0
           endif 
         else
           if(GammaSM(H0VV(i,5)).le.0.0D0)then
            g2hjgg(i)=0.0D0
           else
            g2hjgg(i)= Gamma(H0VV(i,5))/GammaSM(H0VV(i,5))
          endif
         endif

        enddo
  
	if( invAlfa .gt. 0 ) then
          CW2 = (MW/MZ)**2.0D0
          norm = 1/invAlfa*Pi/CW2/(1.0D0-CW2)
	else
          norm=GF*sqrt(2.0D0)*MZ**2.0D0
	endif

        do j=1,3
          do i=1,3   
            g2hjhiZ_nHbynH(j,i)= (
     &          dble(  Coupling(H0HV(j,i)) )**2.0D0
     &        + dimag( Coupling(H0HV(j,i)) )**2.0D0
     &          )
     &       /norm
            
          enddo
        enddo

	end subroutine

c-------------------------------------------------------------------
c-------------------------------------------------------------------
c-------------------------------------------------------------------
	end subroutine
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c-------------------------------------------------------------------

