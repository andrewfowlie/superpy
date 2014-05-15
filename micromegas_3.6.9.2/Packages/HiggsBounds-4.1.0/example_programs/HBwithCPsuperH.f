      PROGRAM CPsuperH2
************************************************************************
* This is modified version of the cpsuperh2.f file which is supplied with 
* CPsuperH2.2 
* (downloaded 07 July 2010
* from http://www.hep.man.ac.uk/u/jslee/CPsuperH.html)
* This file is part of the HiggsBounds distribution.
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
*-----------------------------------------------------------------------
*+CDE HC_ COMMON BLOCKS:
      COMMON /HC_SMPARA/ AEM_H,ASMZ_H,MZ_H,SW_H,ME_H,MMU_H,MTAU_H,MDMT_H
     .                  ,MSMT_H,MBMT_H,MUMT_H,MCMT_H,MTPOLE_H,GAMW_H
     .                  ,GAMZ_H,EEM_H,ASMT_H,CW_H,TW_H,MW_H,GW_H,GP_H
     .                  ,V_H,GF_H,MTMT_H
*
      COMMON /HC_RSUSYPARA/ TB_H,CB_H,SB_H,MQ3_H,MU3_H,MD3_H,ML3_H,ME3_H
*
      COMPLEX*16 MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
      COMMON /HC_CSUSYPARA/ MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX_H/999/
*-----------------------------------------------------------------------
*ARRAYS:
      REAL*8 SMPARA_H(19),SSPARA_H(38)
      DATA NSMIN/19/
      DATA NSSIN/38/
*
      INTEGER*8 IFLAG_H(100)
      DATA NFLAG/100/
*
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(100,3)  ! 100 = NCMAX
      REAL*8     SHC_H(100)
      COMPLEX*16 CHC_H(100)
      DATA NCMAX/100/
*
      REAL*8 GAMBRN(101,3,3)   ! 101 = IFLAG_H(20)+IFLAG_H(21)+1 = NMNH
*                                      ISMN       =ISUSYN        = 50
      REAL*8 GAMBRC(51,3)      !  51 = IFLAG_H(22)+IFLAG_H(23)+1 = NMCH
*                                      ISMC       =ISUSYC        = 25
      DATA NMNH/101/
      DATA NMCH/51/

*-----------------------------------------------------------------------
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *   *    
* used by initialize_HiggsBounds and run_HiggsBounds_part
* HB input:
        integer nHiggsneut,nHiggsplus        
        character*5 whichanalyses

        double precision  Mh(3),GammaTotal_hj(3)
        integer CP_value(3)
        double precision  CS_lep_hjZ_ratio(3),          
     &          CS_lep_bbhj_ratio(3),CS_lep_tautauhj_ratio(3),
     &          CS_lep_hjhi_ratio_nHbynH(3,3),               
     &          CS_gg_hj_ratio(3),CS_bb_hj_ratio(3),  
     &          CS_bg_hjb_ratio(3),                       
     &          CS_ud_hjWp_ratio(3),CS_cs_hjWp_ratio(3),
     &          CS_ud_hjWm_ratio(3),CS_cs_hjWm_ratio(3), 
     &          CS_gg_hjZ_ratio(3),
     &          CS_dd_hjZ_ratio(3),CS_uu_hjZ_ratio(3),
     &          CS_ss_hjZ_ratio(3),CS_cc_hjZ_ratio(3), 
     &          CS_bb_hjZ_ratio(3),                        
     &          CS_tev_vbf_ratio(3),CS_tev_tthj_ratio(3),
     &          CS_lhc7_vbf_ratio(3),CS_lhc7_tthj_ratio(3),
     &          BR_hjss(3),BR_hjcc(3),                         
     &          BR_hjbb(3),BR_hjmumu(3),BR_hjtautau(3),                     
     &          BR_hjWW(3),BR_hjZZ(3),BR_hjZga(3),                     
     &          BR_hjgaga(3),BR_hjgg(3),
     &          BR_hjinvisible(3),BR_hjhihi_nHbynH(3,3)

        double precision Mhplus(1),GammaTotal_Hpj(1), 
     &          CS_lep_HpjHmj_ratio(1),                   
     &          BR_tWpb,BR_tHpjb(1),                      
     &          BR_Hpjcs(1),BR_Hpjcb(1),BR_Hpjtaunu(1) 

* HB output:
        integer HBresult,chan,ncombined  
        double precision obsratio
* misc:
        integer i,j,n
        double precision betasq
        double precision
     &          g2hjVV(3),g2hjbb(3),               
     &          g2hjhiZ_nHbynH(3,3),               
     &          max_hjff_s,max_hjff_p
        integer sneutrino_lspcandidate_number
        logical invisible_lsp
        double precision lspcandidate_mass 

c Set the number of Higgs bosons in the MSSM:
        nHiggsneut=3
        nHiggsplus=1

c The string 'whichanalyses' determines which subset of experimental 
c results are used.
c In this example, we've used the option 'onlyL',
c which instructs HiggsBounds to use tables of results
c from LEP only (i.e. no Tevatron or LHC results).
        whichanalyses='onlyL'

c The subroutine initialize_HiggsBounds reads in all necessary
c tables etc.
c It must be called before calling the run_HiggsBounds_part subroutine.

        call initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses) 

c If you would like to perform scans over variables, the subroutine
c initialize_HiggsBounds (and finish_HiggsBounds) should be called
c outside the do-loops in order to save time.
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*-----------------------------------------------------------------------

*=======================================================================
      CALL FILLINIT2(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,NCMAX,NHC_H,SHC_H,CHC_H,NMNH,GAMBRN,NMCH,GAMBRC)
*=======================================================================
*
* To use other values for the input parameters than those in the "run"
* file, one can specifiy them here. To scan the phase of, for example,
* the gluino mass parameter M_3, one can do
*
*      DO IVAR=0,72
*       SSPARA_H(10)=5.D0*DBLE(IVAR) ! Phi_3
*       print*,'Phi_3 = ',SSPARA_H(10)
*
* Don't forget commenting in "ENDDO ! IVAR" at the end of this block.
*
*-----------------------------------------------------------------------
* For the \sqrt{s}-dependent propagators and the Higgs couplings to the
* gluons and photons, the following should be specified. If not, the
* value in FILLINIT2 is to be used:
*      RAUX_H(101)= ... !  \sqrt{s} for the subroutine FILLDHPG
*-----------------------------------------------------------------------
* One may skip the time-consuming EDM calculations by commenting in the
* follwing line:
*      ISKIP_EDM=1
*-----------------------------------------------------------------------
*
      CALL FILLCPsuperH2(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC)
*
*Error messages:
*--a stop or sbottom squared mass is negative
       IF(IFLAG_H(50).EQ.1) THEN
        print*,'ERROR! IFLAG_H(50) = ',IFLAG_H(50)
        IFLAG_H(50)=0
        GOTO 99
       ENDIF
*
*--the Higgs--boson mass matrix contains a complex or negative eigenvalue
       IF(IFLAG_H(51).EQ.1) THEN
        print*,'ERROR! IFLAG_H(51) = ',IFLAG_H(51)
        IFLAG_H(51)=0
        GOTO 99
       ENDIF
*
*--the diagonalization of the Higgs mass matrix is not successful
       IF(IFLAG_H(52).EQ.1) THEN
        print*,'ERROR! IFLAG_H(52) = ',IFLAG_H(52)
        IFLAG_H(52)=0
        GOTO 99
       ENDIF
*
*--the iteration resumming the threshold corrections is not convergent
       IF(IFLAG_H(54).EQ.1) THEN
        print*,'ERROR! IFLAG_H(54) = ',IFLAG_H(54)
        IFLAG_H(54)=0
        GOTO 99
       ENDIF
*
*--Yukawa coupling has a non--perturbative value: |h_{t,b}| > 2
       IF(IFLAG_H(55).EQ.1) THEN
        print*,'ERROR! IFLAG_H(55) = ',IFLAG_H(55)
        IFLAG_H(55)=0
        GOTO 99
       ENDIF
*
*-- 1 = a tau sneutrino or a stau squared mass is negative
*-- 2 = tachyonic stop or sbottom
*-- 3 = tachyonic scalar strange
       IF(IFLAG_H(56).GT.0) THEN
        IF(IFLAG_H(56).EQ.1) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
        IF(IFLAG_H(56).EQ.2) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
        IF(IFLAG_H(56).EQ.3) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
        IFLAG_H(56)=0
        GOTO 99
       ENDIF
*
*--one of the magnitudes of the complex input parameters is negative
       IF(IFLAG_H(57).EQ.1) THEN
        print*,'ERROR! IFLAG_H(57) = ',IFLAG_H(57)
        IFLAG_H(57)=0
        GOTO 99
       ENDIF
*
*--the iterative method for the neutral Higgs-boson pole masses fails
       IF(IFLAG_H(60).EQ.1) THEN
        print*,'ERROR! IFLAG_H(60) = ',IFLAG_H(60)
        IFLAG_H(60)=0
        GOTO 99
       ENDIF
*
*-----------------------------------------------------------------------
* Users may use the following subroutine for further analysis:
*
*      CALL AURUN(ISKIP_EDM
*     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
*     .,MCH,HMASS_H,OMIX_H
*     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
*     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
*     .,NMNH,GAMBRN,NMCH,GAMBRC)
*-----------------------------------------------------------------------

c ------------------------------------------------------------------
c Set variables needed by HiggsBounds (using results from CPsuperH).
c See HiggsBounds documentation for definition of variables used
c as arguments to run_HiggsBounds_part and CPsuperH 
c documentation for all other variables.

c Note: It is slightly more accurate to use the subroutine run_HiggsBounds_part 
c rather than the subroutine run_HiggsBounds_effC because the SM branching ratios
c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
c ratios used in CPsuperH

        do i=1,3

         Mh(i)=HMASS_H(i)         
         GammaTotal_hj(i)=GAMBRN(IFLAG_H(20)+IFLAG_H(21)+1,1,i)

         BR_hjss(i)        = GAMBRN(5,3,i)
         BR_hjcc(i)        = GAMBRN(8,3,i)
         BR_hjbb(i)        = GAMBRN(6,3,i)
         BR_hjmumu(i)      = GAMBRN(2,3,i)
         BR_hjtautau(i)    = GAMBRN(3,3,i) 
         BR_hjWW(i)        = GAMBRN(10,3,i) 
         BR_hjZZ(i)        = GAMBRN(11,3,i) 
         BR_hjgaga(i)      = GAMBRN(17,3,i) 
         BR_hjgg(i)        = GAMBRN(18,3,i) 

         sneutrino_lspcandidate_number=0
         invisible_lsp=.True.
         lspcandidate_mass=MN_H(1)

         if( SNU3MASS_H .lt. lspcandidate_mass )then
            lspcandidate_mass=SNU3MASS_H
            sneutrino_lspcandidate_number=3
         endif

         if(     MC_H(1)       .lt. lspcandidate_mass )then !chargino
            invisible_lsp=.False.
         elseif( SSPara_H(9)   .lt. lspcandidate_mass )then !gluino
            invisible_lsp=.False.
         elseif( STMASS_H(1)   .lt. lspcandidate_mass )then !stop
            invisible_lsp=.False.
         elseif( SBMASS_H(1)   .lt. lspcandidate_mass )then !sbottom
            invisible_lsp=.False.
         elseif( STAUMASS_H(1) .lt. lspcandidate_mass )then !stau
            invisible_lsp=.False.
         endif

         if(invisible_lsp)then
           if(    sneutrino_lspcandidate_number.eq.0)then
             BR_hjinvisible(i)=GAMBRN(IFLAG_H(20)+1,3,i)
           elseif(sneutrino_lspcandidate_number.eq.3)then 
             BR_hjinvisible(i)=GAMBRN(IFLAG_H(20)+27,3,i)
           endif
         else
           BR_hjinvisible(i)=0.0D0
         endif


!this branching ratio is not calculated by CPsuperH, so we set it to zero
         BR_hjZga(i) = 0.0D0

         g2hjbb(i)=     
     &      abs(NHC_H(17,i))**2.0D0   
     &    + abs(NHC_H(18,i))**2.0D0

         CS_bg_hjb_ratio(i)       = g2hjbb(i)
         CS_bb_hj_ratio(i)        = g2hjbb(i)
         CS_lep_bbhj_ratio(i)     = g2hjbb(i)
         CS_lep_tautauhj_ratio(i) =      
     &      abs(NHC_H(8,i))**2.0D0   
     &    + abs(NHC_H(9,i))**2.0D0

         g2hjVV(i)= abs(NHC_H(70,i))**2.0D0

         CS_lep_hjZ_ratio(i)    = g2hjVV(i)
         CS_dd_hjZ_ratio(i)     = g2hjVV(i)
         CS_uu_hjZ_ratio(i)     = g2hjVV(i)
         CS_ss_hjZ_ratio(i)     = g2hjVV(i)
         CS_cc_hjZ_ratio(i)     = g2hjVV(i)
         CS_bb_hjZ_ratio(i)     = g2hjVV(i)
         CS_ud_hjWp_ratio(i)    = g2hjVV(i)
         CS_cs_hjWp_ratio(i)    = g2hjVV(i)
         CS_ud_hjWm_ratio(i)    = g2hjVV(i)
         CS_cs_hjWm_ratio(i)    = g2hjVV(i)
         CS_tev_vbf_ratio(i)    = g2hjVV(i)
         CS_lhc7_vbf_ratio(i)   = g2hjVV(i)

         CS_gg_hjZ_ratio(i)     = 0.0D0

         CS_tev_tthj_ratio(i) =     
     &      abs(NHC_H(26,i))**2.0D0   
     &    + abs(NHC_H(27,i))**2.0D0
         CS_lhc7_tthj_ratio(i) = CS_tev_tthj_ratio(i)
c ------------------------------------------------------------------
!note that this is an approximation
         CS_gg_hj_ratio(i) =  GAMBRN(18,1,i)    
     &          /( SMBR_Hgg(Mh(i))    *SMGamma_h(Mh(i)) )
         if(SMGamma_h(Mh(i)).lt.0)then
           CS_gg_hj_ratio(i) = 0.0D0
           !it's ok to set this to zero, because Mh(i) is out of range anyway
         endif 
c ------------------------------------------------------------------

	 BR_hjhihi_nHbynH(i,1)=GAMBRN(14,3,i)
	 BR_hjhihi_nHbynH(i,2)=GAMBRN(16,3,i)	 
	 BR_hjhihi_nHbynH(i,3)=0.0D0
 
         max_hjff_s=max(abs(NHC_H(14,i))**2.0D0,  
     &                  abs(NHC_H(23,i))**2.0D0,
     &                  abs(NHC_H(17,i))**2.0D0,
     &                  abs(NHC_H(26,i))**2.0D0,
     &                  abs(NHC_H( 8,i))**2.0D0    )
  
         max_hjff_p=max(abs(NHC_H(15,i))**2.0D0,  
     &                  abs(NHC_H(24,i))**2.0D0,
     &                  abs(NHC_H(18,i))**2.0D0,
     &                  abs(NHC_H(27,i))**2.0D0,
     &                  abs(NHC_H( 9,i))**2.0D0  )

         if(     max_hjff_p .lt. 1.0D-16 )then !CP even
          CP_value(i) =  1
         elseif( max_hjff_s .lt. 1.0D-16 )then !CP odd
          CP_value(i) = -1
         else                              !mixed CP
          CP_value(i) =  0
         endif

        enddo
	
	do j=1,3    
	 do i=1,3 
	  if(i.lt.j)then	  
	   g2hjhiZ_nHbynH(j,i)=g2hjVV(6-j-i)
	   g2hjhiZ_nHbynH(i,j)=g2hjhiZ_nHbynH(j,i)
	  else
	   g2hjhiZ_nHbynH(j,i)=0.0D0	   
	  endif    
	 enddo
	enddo 

	do j=1,3    
	 do i=1,3        
          CS_lep_hjhi_ratio_nHbynH(j,i) = g2hjhiZ_nHbynH(j,i)
	 enddo
	enddo 

        Mhplus(1)=SSPARA_H(2)
        GammaTotal_Hpj(1)=GAMBRC(IFLAG_H(22)+IFLAG_H(23)+1,1)
        CS_lep_HpjHmj_ratio(1)=1.0D0
        BR_Hpjcs(1)      = GAMBRC(5,3)
        BR_Hpjtaunu(1)   = GAMBRC(3,3)
! this branching ratios is not calculated by CPsuperH, so set to zero:
        BR_Hpjcb(1)      = 0.0D0
! t-quark branching ratios are not calculated by CPsuperH, so set to zero:
        BR_tWpb          = 0.0D0
        BR_tHpjb(1)      = 0.0D0

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
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
     &          BR_hjss,BR_hjcc,                             
     &          BR_hjbb,BR_hjmumu,BR_hjtautau,               
     &          BR_hjWW,BR_hjZZ,BR_hjZga, BR_hjgaga,BR_hjgg, 
     &          BR_hjinvisible,BR_hjhihi_nHbynH              )

        call HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj, 
     &          CS_lep_HpjHmj_ratio,                        
     &          BR_tWpb,BR_tHpjb,                           
     &          BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu)

        call run_HiggsBounds( HBresult,chan,                  
     &                      obsratio, ncombined              )


        write(*,*)        
        write(*,*)'*************    HiggsBounds Results  **************'
        write(*,*) 
        write(*,*)'Is this parameter point excluded at 95% CL?'
        write(*,*) HBresult, ',  where'
        write(*,*)'               0 = yes, it is excluded'
        write(*,*)'               1 = no, it has not been excluded'
        write(*,*)'              -1 = invalid parameter set'    
        write(*,*)
        write(*,*)'The process with the highest statistical sensitivity'
        write(*,*)'is'
        write(*,*) chan,'(see Key.dat)'
        write(*,*)'This process has a theoretical rate vs. limit of'
        write(*,*) obsratio
        write(*,*)
        write(*,*)'The number of Higgs bosons which have contributed to'
        write(*,*)'the theoretical rate of this process was'
        write(*,*) ncombined
        write(*,*)
        write(*,*)'See HiggsBounds documentation for more information.'
        write(*,*)'****************************************************'
        write(*,*)

 99   CONTINUE
*
*      ENDDO ! IVAR
*=======================================================================
       
c ------------------------------------------------------------------
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c deallocates arrays used by HiggsBounds:

        call finish_HiggsBounds
c ------------------------------------------------------------------

      STOP
      END
