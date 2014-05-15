!--------------------------------------------------------------------
	program HBandHSwithFH
!--------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 06/03/2013).
!
! The program runs FeynHiggs (FH version > 2.9.4 required) to evaluate
! the partonic input of HiggsBounds/HiggsSignals for a example MSSM point.
! Then, we run HiggsBounds and HiggsSignals.
!
! Note: You have to set the correct FeynHiggs path in the configure file.
!       Then, this example program can be compiled by running
!
!       make HBandHSwithFH
!
!       in the HiggsSignals main directory.
!--------------------------------------------------------------------
	implicit none

		integer error

c used by FHHiggsCorr
		double precision MHiggs(4)
		double complex SAeff, UHiggs(3,3), ZHiggs(3,3)

c used by FHSelectUZ:
		integer uzint, uzext, mfeff

c used by FHCouplings:
#include "FHCouplings.h"
		double complex couplings(ncouplings), couplingsms(ncouplingsms)
		double precision gammas(ngammas), gammasms(ngammasms)
		integer fast

c used by FHHiggsProd:
	    double precision sqrts, prodxs(nprodxs)

c used by FHGetPara:
		integer nmfv
		double precision MSf(2,4,3),MASf(6,4), MCha(2), MNeu(4)
		double complex USf(2,2,4,3),UASf(6,6,4)
		double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
		double complex DeltaMB
		double precision MGl
		double precision MHtree(4), SAtree

c used by FHRetrieveSMPara:	
		double precision invAlfa, AlfasMZ, GF
		double precision ME, MU, MD, MM, MC, MS, ML, MB
		double precision MW, MZ
		double precision CKMlambda, CKMA, CKMrhobar, CKMetabar

c used by initialize_HiggsBounds
        integer nHiggsneut,nHiggsplus
        parameter (nHiggsneut = 3)
        parameter (nHiggsplus = 1)
       	character(LEN=5) whichanalyses

c used by HiggsBounds_neutral_input_part
        double precision Mh(3),GammaTotal_hj(3)
        integer CP_value(3)
        double precision  CS_lep_hjZ_ratio(3),          
     &   CS_lep_bbhj_ratio(3),CS_lep_tautauhj_ratio(3),
     &   CS_lep_hjhi_ratio_nHbynH(3,3),               
     &   CS_gg_hj_ratio(3),CS_bb_hj_ratio(3),  
     &   CS_bg_hjb_ratio(3),                       
     &   CS_ud_hjWp_ratio(3),CS_cs_hjWp_ratio(3),
     &   CS_ud_hjWm_ratio(3),CS_cs_hjWm_ratio(3), 
     &   CS_gg_hjZ_ratio(3),
     &   CS_dd_hjZ_ratio(3),CS_uu_hjZ_ratio(3),
     &   CS_ss_hjZ_ratio(3),CS_cc_hjZ_ratio(3), 
     &   CS_bb_hjZ_ratio(3),                        
     &   CS_tev_vbf_ratio(3),CS_tev_tthj_ratio(3),
     &   CS_lhc7_vbf_ratio(3),CS_lhc7_tthj_ratio(3),
     &   CS_lhc8_vbf_ratio(3),CS_lhc8_tthj_ratio(3),
     &   BR_hjss(3),BR_hjcc(3),                         
     &   BR_hjbb(3),BR_hjmumu(3),BR_hjtautau(3),                     
     &   BR_hjWW(3),BR_hjZZ(3),BR_hjZga(3),                     
     &   BR_hjgaga(3),BR_hjgg(3),
     &   BR_hjinvisible(3),BR_hjhihi_nHbynH(3,3)

c used by HiggsBounds_charged_input
        double precision Mhplus(1),GammaTotal_Hpj(1),
     &   CS_lep_HpjHmj_ratio(1),                    
     &   BR_tWpb,BR_tHpjb(1),                     
     &   BR_Hpjcs(1),BR_Hpjcb(1),BR_Hpjtaunu(1) 

c used in HiggsBounds
        double precision dmhneut_hb(nHiggsneut)
        double precision dmhch_hb(nHiggsplus)
c used in HiggsSignals (can be different)
        double precision dmhneut_hs(nHiggsneut)
        
c return values of run_HiggsBounds
        integer HBresult,chan,ncombined
        double precision obsratio

c used by set_rate_uncertainties in HiggsSignals
        double precision dCS(5),dBR(5)
        double precision dggh, dbbh

c return values of run_HiggsSignals
        double precision Chisq_mu, Chisq_mh, Chisq, Pvalue
        integer nobs
        
c run options in HiggsSignals
        integer pdf, output_level, mass_centered_method,runmode
             
c misc:
        integer i,j,as,t
        double precision norm,CW2,Pi
        double precision
     &   g2hjbb(3),g2hjWW(3),g2hjZZ(3),                
     &   g2hjgg(3),g2hjhiZ_nHbynH(3,3)
        double precision g2hjbb_s(3),g2hjbb_p(3)
        double precision g2hjtautau_s(3),g2hjtautau_p(3)
        integer sneutrino_lspcandidate_number
        logical invisible_lsp
        double precision lspcandidate_mass 

        Pi = 3.1415926535897932384626433832795029D0

!--------------------------------------- Run Options ------------------------------------!
c The string 'whichanalyses' determines which subset of experimental 
c results are used. In this example, we've used the option 'LandH',
c which instructs HiggsBounds to use tables of results
c from LEP, Tevatron and LHC (i.e. the full set of
c results supplied with HiggsBounds).
        whichanalyses='LandH'
c Choose a Higgs mass parametrization for HiggsSignals 
c (1: box, 2:gaussian, 3:box+gaussian)
        pdf = 2 
c Choose the output level of HiggsSignals
c (0: silent, 1: screen output, 2: even more output,...)
        output_level = 0
c Set the HiggsSignals runmode (1: peak-centered, 2: mass-centered, 3: both) Chi^2 method(s)
		runmode = 1
c If the mass-centered chi^2 method is used, can specify the dm_theory treatment
c (1: mass variation, 2: smearing of mu-plot with mass pdf)        
        mass_centered_method = 1        
!---------------------------- HiggsBounds and HiggsSignals ------------------------------!
c
c If you would like to perform scans over variables, the subroutines
c initialize_HiggsBounds, initialize_HiggsSignals (and finish_HiggsBounds,
c finish HiggsSignals) should be called outside the do-loops in order to save time.
c
!---- Initialize HiggsBounds and specify the dataset it should use 				     ----!
        call initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses)
!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
		call initialize_HiggsSignals(nHiggsneut,nHiggsplus,
     &                               "latestresults")  
!------------------------------ HiggsSignals options ------------------------------------!
!---- Set the output level                                                      	 ----!
		call setup_output_level(output_level)
!---- Set the Higgs mass parametrization                                     	 	 ----!
		call setup_pdf(pdf) 
!---- If the mass-centered chi^2 method is used, can specify the dm_theory treatment ----!
!     1: mass variation, 2: smearing of mu-plot with mass pdf  
        call setup_mcmethod_dm_theory(mass_centered_method)
!---------------------------------- Run FeynHiggs ---------------------------------------!
c the subroutines setFlags, setPara, setSLHA
c are also contained in this file
	call setFlags

* either use setPara to set the parameters directly
* or use setSLHA to read them from an SLHA file
	call setPara
c	call setSLHA("tmp.slha")

c       Use this for FH < 2.9.5
c	call FHGetPara(error, nmfv, MASf, UASf,
	call FHGetPara(error, nmfv, MSf, USf,MASf, UASf,
     &    MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
     &    MHtree, SAtree)
	if( error .ne. 0 ) stop

	call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)
	if( error .ne. 0 ) stop

c NOTE: we are setting uzint=uzext
        mfeff=1
        uzint=2
        uzext=2
        call FHSelectUZ(error, uzint, uzext, mfeff)
	    if( error .ne. 0 ) stop

        fast=1
	call FHCouplings(error,
     &    couplings, couplingsms, gammas, gammasms, fast)
	if( error .ne. 0 ) stop

c We would like FH to calculate LHC cross sections
        sqrts=8.0D0
        call FHHiggsProd(error, sqrts, prodxs)
	if( error .ne. 0 ) stop 

        call FHRetrieveSMPara(error,
     &    invAlfa, AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MB,
     &    MW, MZ,
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)

!--------------- Turn FeynHiggs results into HiggsBounds/Signals input ------------------!

c Set variables needed by HiggsBounds (using results from FeynHiggs).
c See HiggsBounds documentation for definition of variables used
c as arguments to HiggsBounds_neutral_input_part and run_HiggsBounds
c and FeynHiggs documentation for all other variables.

c Note: It is slightly more accurate to use the subroutine HiggsBounds_neutral_input_part
c rather than the subroutine HiggsBounds_neutral_input_effC because the SM branching ratios
c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
c ratios used in FeynHiggs

        do i=1,3
         Mh(i)=MHiggs(i) 
         GammaTotal_hj(i) = GammaTot(i)

         BR_hjss(i)       = BR(H0FF(i,4,2,2))
         BR_hjcc(i)       = BR(H0FF(i,3,2,2))
         BR_hjbb(i)       = BR(H0FF(i,4,3,3))
         BR_hjmumu(i)     = BR(H0FF(i,2,2,2))
         BR_hjtautau(i)   = BR(H0FF(i,2,3,3)) 

         BR_hjWW(i)     = BR(H0VV(i,4))
         BR_hjgaga(i)   = BR(H0VV(i,1))
         BR_hjZga(i)    = BR(H0VV(i,2))
         BR_hjZZ(i)     = BR(H0VV(i,3))
         BR_hjgg(i)     = BR(H0VV(i,5))

         if(GammaSM(H0FF(i,4,3,3)).le.0.0D0)then
          g2hjbb(i)=0.0D0
         else
          g2hjbb(i)=Gamma(H0FF(i,4,3,3)) 
     &      /GammaSM(H0FF(i,4,3,3))
         endif

c Note that this is currently equivalent to
c         g2hjbb(i)= bbh(i)/bbhSM(i)
c         g2hjbb(i)= btagbh(i)/btagbhSM(i)
c as long as MH>80 GeV

         CS_bg_hjb_ratio(i) = g2hjbb(i)
         CS_bb_hj_ratio(i)  = g2hjbb(i)
   
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

         if(       g2hjbb_p(i).lt.1.0D-10)then
            CP_value(i) = 1
         elseif(   g2hjbb_s(i).lt.1.0D-10)then
            CP_value(i) = -1
         else
            CP_value(i) = 0
         endif

         CS_lep_bbhj_ratio(i)     = g2hjbb_s(i)+g2hjbb_p(i)
         CS_lep_tautauhj_ratio(i) = g2hjtautau_s(i)+g2hjtautau_p(i)

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

         CS_lep_hjZ_ratio(i)        = g2hjZZ(i)

         CS_gg_hjZ_ratio(i)     = 0.0D0
         CS_dd_hjZ_ratio(i)     = g2hjZZ(i)
         CS_uu_hjZ_ratio(i)     = g2hjZZ(i)
         CS_ss_hjZ_ratio(i)     = g2hjZZ(i)
         CS_cc_hjZ_ratio(i)     = g2hjZZ(i)
         CS_bb_hjZ_ratio(i)     = g2hjZZ(i)

         CS_ud_hjWp_ratio(i)    = g2hjZZ(i)
         CS_cs_hjWp_ratio(i)    = g2hjZZ(i)
         CS_ud_hjWm_ratio(i)    = g2hjZZ(i)
         CS_cs_hjWm_ratio(i)    = g2hjZZ(i)

         CS_tev_vbf_ratio(i)     = g2hjZZ(i)
         CS_lhc7_vbf_ratio(i)    = g2hjZZ(i)
         CS_lhc8_vbf_ratio(i)    = g2hjZZ(i)


         if(tthSM(i).gt.0.0D0)then
           CS_tev_tthj_ratio(i)    = tth(i)/tthSM(i)  
         else
           CS_tev_tthj_ratio(i)    = 0.0D0
         endif

         CS_lhc7_tthj_ratio(i) = CS_tev_tthj_ratio(i)
         CS_lhc8_tthj_ratio(i) = CS_tev_tthj_ratio(i)

c tevatron gluon fusion XS is not calculated in FH is MH<90 geV
         if(Mh(i).gt.90.0001D0)then      
           if(gghSM(i).gt.0.0D0)then
            CS_gg_hj_ratio(i) = ggh(i)/gghSM(i)
           else
            CS_gg_hj_ratio(i) = 0.0D0
           endif 
         else
           if(GammaSM(H0VV(i,5)).le.0.0D0)then
            CS_gg_hj_ratio(i)=0.0D0
           else
            CS_gg_hj_ratio(i)= Gamma(H0VV(i,5))/GammaSM(H0VV(i,5))
          endif
         endif

        enddo
  
        norm=GF*sqrt(2.0D0)*MZ**2.0D0

        do j=1,3
          do i=1,3   
            g2hjhiZ_nHbynH(j,i)= (
     &          dble(  Coupling(H0HV(j,i)) )**2.0D0
     &        + dimag( Coupling(H0HV(j,i)) )**2.0D0
     &          )
     &       /norm
            
            CS_lep_hjhi_ratio_nHbynH(j,i) = g2hjhiZ_nHbynH(j,i)

            BR_hjhihi_nHbynH(j,i)=BR(H0HH(j,i,i))
          enddo
        enddo

c higgs->neutralino1 neutralino1 contributes the invisible Higgs decay width
c when neutralino1 or sneutrino is the LSP

        do i=1,3  
         sneutrino_lspcandidate_number=0
         invisible_lsp=.True.

c first determine whether lightest sneutrino is lighter than the lightest neutralino
c
c sneutrino_lspcandidate_number=0 indicates that lightest neutralino is 
c lighter than all the sneutrinos
         lspcandidate_mass=MNeu(1)
         do as=1,3
            if( MASf(as,1) .lt. lspcandidate_mass )then
               lspcandidate_mass=MASf(as,1)
               sneutrino_lspcandidate_number=as
            endif
         enddo

         if(     MCha(1) .lt. lspcandidate_mass )then 
           invisible_lsp=.False.
         elseif( MGl     .lt. lspcandidate_mass )then
           invisible_lsp=.False.
         else
          do as=1,6
           do t=2,4
            if( MASf(as,t) .lt. lspcandidate_mass )then
               invisible_lsp=.False.
            endif
           enddo
          enddo      
         endif

         if(invisible_lsp)then
          if(sneutrino_lspcandidate_number.eq.0)then
            BR_hjinvisible(i) = BR(H0NeuNeu(i,1,1))
          else
            BR_hjinvisible(i) = BR(H0SfSf(i,1,1,1,as))
          endif
         else
          BR_hjinvisible(i) = 0.0D0
         endif
        enddo

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c	    Charged Higgs input

        Mhplus(1)              = MHiggs(4) 
        GammaTotal_Hpj(1)      = GammaTot(4)  
        CS_lep_HpjHmj_ratio(1) = 1.0D0
        BR_tWpb                = BR( tBF(1) ) 
        BR_tHpjb(1)            = BR( tBF(2) )
        BR_Hpjcs(1)            = BR( HpFF(2,2,2) ) 
        BR_Hpjcb(1)            = BR( HpFF(2,2,3) )
        BR_Hpjtaunu(1)         = BR( HpFF(1,3,3) )

		print*,' '
		print*,'The Higgs boson masses are (h, H, A, H+):'
		print*,Mh,Mhplus
!--------------------------------- Run HiggsBounds --------------------------------------!

c calls to HiggsBounds_neutral_input_part,HiggsBounds_charged_input,set_mass_uncertainties
c which give input to HiggsBounds

c        print*,'calling HiggsBounds_neutral_input_part'  
        call HiggsBounds_neutral_input_part(Mh,GammaTotal_hj,CP_value, 
     &          CS_lep_hjZ_ratio,                            
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,     
     &          CS_lep_hjhi_ratio_nHbynH,                    
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,       
     &          CS_bg_hjb_ratio,                         
     &          CS_ud_hjWp_ratio,CS_cs_hjWp_ratio,    
     &          CS_ud_hjWm_ratio,CS_cs_hjWm_ratio,    
     &          CS_gg_hjZ_ratio,     
     &          CS_dd_hjZ_ratio,CS_uu_hjZ_ratio,     
     &          CS_ss_hjZ_ratio,CS_cc_hjZ_ratio,     
     &          CS_bb_hjZ_ratio,                         
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,    
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,    
     &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,    
     &          BR_hjss,BR_hjcc,                             
     &          BR_hjbb,BR_hjmumu,BR_hjtautau,               
     &          BR_hjWW,BR_hjZZ,BR_hjZga, BR_hjgaga,BR_hjgg,  
     &          BR_hjinvisible,BR_hjhihi_nHbynH              )
     
c        print*,'calling HiggsBounds_charged_input'
        call HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj, 
     &          CS_lep_HpjHmj_ratio,                        
     &          BR_tWpb,BR_tHpjb,                           
     &          BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu)

c Set the theory mass uncertainties. HiggsBounds will vary the masses accordingly.
		dmhneut_hb = (/2.D0, 0.D0, 0.D0/)
		dmhch_hb   = (/0.D0/)
        call HiggsBounds_set_mass_uncertainties(dmhneut_hb,dmhch_hb)

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c call to run_HiggsBounds
        call run_HiggsBounds( HBresult,chan,obsratio,ncombined)

        print*,' '      
        print*,'*************    HiggsBounds Results  **************'
        print*,' '    
        print*,'Is this parameter point excluded by LEP, Tevatron'
        print*,'or LHC data?'         
        print*, HBresult, ',  where'
        print*,'               0 = yes, it is excluded'
        print*,'               1 = no, it has not been excluded'
        print*,'              -1 = invalid parameter set'    
        print*,' '    
        print*,'The process with the highest statistical sensitivity'
        print*,'is'
        print*, chan,'(see Key.dat)'
        print*,'This process has a theoretical rate vs. limit of'
        print*, obsratio
        print*,' '    
        print*,'The number of Higgs which have contributed to the'
        print*,'theoretical rate of this process was'
        print*, ncombined
        print*,' '    
        print*,'See HiggsBounds documentation for more information.'
        print*,'****************************************************'
        print*,' '    

!--------------------------------- Run HiggsSignals -------------------------------------!

c Set the theory mass uncertainties. This will affect the Chi^2 contribution
c from the mass measurements as well as the Higgs-to-peaks assignment in
c HiggsSignals .
		dmhneut_hs = (/2.D0, 0.D0, 0.D0/)
		call HiggsSignals_neutral_input_MassUncertainty(dmhneut_hs)
c Set the MSSM rate uncertainties:
!---dCS and dBR hold the model's rate uncertainties. Can be changed by user
!---with subroutine setup_rate_uncertainties. First, set to SM values.
		dCS = (/ 0.147D0, 0.028D0, 0.037D0, 0.051D0, 0.12D0 /)  
		dBR = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
c Now, we set relative uncertainties for the processes gg->H and bb->H
c and interpolate from the effective couplings/CS-ratios (of the lightest Higgs)
c the uncertainty of single Higgs production, dCS(1).
		dggh = 0.147D0
		dbbh = 0.200D0		
		call get_singleH_uncertainty(dCS(1),dggh,dbbh,CS_gg_hj_ratio(1),
     &                                 g2hjbb_s(1)+g2hjbb_p(1),Mh(1))
        print*,"Interpolated singleH uncertainty (in %): ",dCS(1)
		call setup_rate_uncertainties(dCS, dBR)
!---- Run HiggsSignals																 ----!
	    call run_HiggsSignals(runmode,Chisq_mu,Chisq_mh,
     &                        Chisq,nobs,Pvalue)
		
		print*,' '      
        print*,'*************  HiggsSignals Results  **************'
		print*,' '
		print*,'The chi^2 contr. from signal strength observables is'
		print*,Chisq_mu
		print*,'The chi^2 contr. from Higgs mass observables is'
		print*,Chisq_mh
		print*,'The total chi^2 is'
		print*,Chisq
		print*,'The number of observables is'
		print*,nobs
		print*,' '
		print*,'In the case that the number of observables is equal'
		print*,'to the number of degrees of freedom, this corresponds'
		print*,'to a probability of'
		print*,Pvalue
        print*,' '    
        print*,'See HiggsSignals documentation for more information.'
        print*,'****************************************************'
        print*,' '    

!---------------------- Finish HiggsBounds and HiggsSignals -----------------------------!

c deallocates arrays used by HiggsBounds:

        call finish_HiggsBounds
        call finish_HiggsSignals

	end


************************************************************************

	subroutine setFlags
	implicit none

	integer mssmpart, fieldren, tanbren, higgsmix, p2approx
	integer looplevel, runningMT, botResum, tlCplxApprox

c   Using default (recommended) values of all FH flags
	parameter (mssmpart = 4)
	parameter (fieldren = 0)
	parameter (tanbren = 0)
	parameter (higgsmix = 2)
	parameter (p2approx = 0)
	parameter (looplevel = 2)
	parameter (runningMT = 1)
	parameter (botResum = 1)
	parameter (tlCplxApprox = 0)

	integer error

	call FHSetFlags(error, mssmpart, fieldren, tanbren,
     &    higgsmix, p2approx, looplevel,
     &    runningMT, botResum, tlCplxApprox)
	if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setPara
	implicit none

	double precision invAlfa, AlfasMZ, GF
	double precision ME, MU, MD, MM, MC, MS, ML, MB, MZ, MW
	double precision CKMlambda, CKMA, CKMrhobar, CKMetabar
	parameter (invAlfa = -1)
	parameter (AlfasMZ = -1)
	parameter (GF = -1)
	parameter (ME = -1)
	parameter (MU = -1)
	parameter (MD = -1)
	parameter (MM = -1)
	parameter (MC = -1)
	parameter (MS = -1)
	parameter (ML = -1)
	parameter (MB = -1)
	parameter (MW = -1)
	parameter (MZ = -1)
	parameter (CKMlambda = -1)
	parameter (CKMA = -1)
	parameter (CKMrhobar = -1)
	parameter (CKMetabar = -1)

	double precision MT, TB, MA0, MHp
	parameter (MT = 173.2)
	parameter (TB = 8.)
	parameter (MA0 = 300)
	parameter (MHp = -1)

	double precision MSusy
	double precision M3SL, M2SL, M1SL
	double precision M3SE, M2SE, M1SE
	double precision M3SQ, M2SQ, M1SQ
	double precision M3SU, M2SU, M1SU
	double precision M3SD, M2SD, M1SD
	parameter (MSusy = 1000)
	parameter (M3SL = MSusy)
	parameter (M2SL = M3SL)
	parameter (M1SL = M2SL)
	parameter (M3SE = MSusy)
	parameter (M2SE = M3SE)
	parameter (M1SE = M2SE)
	parameter (M3SQ = MSusy)
	parameter (M2SQ = M3SQ)
	parameter (M1SQ = M2SQ)
	parameter (M3SU = MSusy)
	parameter (M2SU = M3SU)
	parameter (M1SU = M2SU)
	parameter (M3SD = MSusy)
	parameter (M2SD = M3SD)
	parameter (M1SD = M2SD)

	double complex Atau, At, Ab
	double complex Amu, Ac, As
	double complex Ae, Au, Ad
	parameter (At = 2000)
	parameter (Ab = At)
	parameter (Atau = At)
	parameter (Ac = At)
	parameter (As = Ab)
	parameter (Amu = Atau)
	parameter (Au = Ac)
	parameter (Ad = As)
	parameter (Ae = Amu)

	double complex MUE, M_1, M_2, M_3
	parameter (MUE = 200)
	parameter (M_1 = 100)
	parameter (M_2 = 200)
	parameter (M_3 = 800)

	double precision Qtau, Qt, Qb
	parameter (Qtau = 0)
	parameter (Qt = 0)
	parameter (Qb = 0)

	double precision scalefactor
	parameter (scalefactor = 1)

	integer error

	call FHSetSMPara(error,
     &    invAlfa, AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MB,
     &    MW, MZ,
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)
	if( error .ne. 0 ) stop

	call FHSetPara(error, scalefactor,
     &    MT, TB, MA0, MHp,
     &    M3SL, M3SE, M3SQ, M3SU, M3SD,
     &    M2SL, M2SE, M2SQ, M2SU, M2SD,
     &    M1SL, M1SE, M1SQ, M1SU, M1SD,
     &    MUE,
     &    Atau, At, Ab,
     &    Amu, Ac, As,
     &    Ae, Au, Ad,
     &    M_1, M_2, M_3,
     &    Qtau, Qt, Qb)
        if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setSLHA(filename)
	implicit none
	character*(*) filename

#include "SLHA.h"

	integer error
	double complex slhadata(nslhadata)

	call SLHARead(error, slhadata, filename, 1)
	if( error .ne. 0 ) stop

	call FHSetSLHA(error, slhadata)
	if( error .ne. 0 ) stop
	end


!************************************************************** 
	subroutine get_singleH_uncertainty(dCS,dggh,dbbh,g2hgg,g2hbb,mh)
	use theory_colliderSfunctions

	implicit none	
	double precision dCS
 	double precision dggh, dbbh, g2hgg, g2hbb, mh
 	double precision vsmall
 	
 	vsmall=1.0D-10
 	
	if(g2hgg.le.vsmall.and.g2hbb.le.vsmall) then
		dCS = 0.0D0
	else 
		dCS=(g2hgg*LHC8_rH_gg(mh)*dggh+
     & 		g2hbb*LHC8_rH_bb(mh)*dbbh)/
     &      (g2hgg*LHC8_rH_gg(mh)+g2hbb*LHC8_rH_bb(mh))
	endif
 
	end
!************************************************************** 

