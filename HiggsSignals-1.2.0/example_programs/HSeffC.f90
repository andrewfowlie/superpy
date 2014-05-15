!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 29/01/2013).
!--------------------------------------------------------------------------------------
program HSeffC
! In this example we use the effective coupling input to scan over the effective
! Higgs-gluon-gluon and Higgs-b-b couplings. All other couplings are as in the SM, except
! for the loop-induced Higgs-gamma-gamma coupling, which is derived from the tree-level
! couplings with the function get_g2hgaga (assuming a Higgs mass of 126 GeV, see below).
! The evaluated chi^2 is saved in the datafile "results/HSeffC.dat"
! The single-Higgs production cross section is composed of gg->H and bb->H. Both 
! processes have different uncertainties and scale with different couplings. The function
! get_singleH_uncertainty employs internal ratio functions of these production processes 
! from HiggsBounds to interpolate the uncertainty of single Higgs production from the
! uncertainties of its subprocesses.
!--------------------------------------------------------------------------------------
 use theory_colliderSfunctions
 use usefulbits, only : vsmall
! use pc_chisq, only : print_corr_mu_to_file
 implicit none

 integer :: nHzero, nHplus, ndf, i, j, k, ii, jj
 double precision :: obsratio, mass, Pvalue, Chisq, mu, Chisq_mu, Chisq_mh
 double precision :: dCS(5),dBR(5)
 double precision :: SMGammaTotal
double precision :: scale_bbh, scale_ggh, dggh, dbbh
 double precision :: Mh,GammaTotal,g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p, &
&                    g2hjbb_s,g2hjbb_p,g2hjtt_s,g2hjtt_p, &
&                    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p, &
&                    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,	&
&                    g2hjhiZ,BR_hjhihi,BR_hjinvisible
 character(len=100)::filename
 double precision :: dm
 integer		  :: pdf
!-HiggsBounds internal functions to obtain SM branching ratios 
 double precision :: SMBR_Htoptop,SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Hmumu, SMBR_Htautau,&
 &                   SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg,SMGamma_h
 double precision :: Htobb_rate, singleH_rate
 nHzero=1
 nHplus=0

!--Set the (relative!) rate uncertainties for your model:
!  dCS(1) - singleH				dBR(1) - gamma gamma
!  dCS(2) - VBF					dBR(2) - W W
!  dCS(3) - HW					dBR(3) - Z Z
!  dCS(4) - HZ					dBR(4) - tau tau
!  dCS(5) - ttH					dBR(5) - b bbar
!--Here, we first set them to the SM. Later, we determine the uncertainty of singleH
!--production, dCS(1), from the effective couplings. Taken SM values from
!--https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CrossSections
 dCS = (/ 0.147D0, 0.028D0, 0.037D0, 0.051D0, 0.12D0 /)
 dBR = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
!--Set the rate uncertainties of gluon-gluon fusion and bb->H here:
 dggh = 0.147D0
 dbbh = 0.200D0

!--Setting up the output
 filename='results/HSeffC.dat'
 open(21,file=filename)
 write(21,*) '# mh   scale_ggh   scale_bbh   d(singleH)   Chisq_mu   Chisq    ndf',&
& '    rate(singleH)    rate(singleH->bb)' 
 write(21,*) '#--------------------------------------------------------------------',&
& '-------------------------------------#'

!--Enter the Higgs mass and its theory uncertainty here: 
 Mh = 126.0D0
 dm = 0.0D0

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder ----!
 call initialize_HiggsSignals(nHzero,nHplus,"latestresults")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) ----!
 call setup_output_level(0)
 !---- Set the assignment range for the peak-centered method (optional)				 ----! 
 call setup_assignmentrange(10.0D0)
 !---- Disable the use of CS and BR uncertainties via provided covariance matrices	 ----!  
 call setup_correlated_rate_uncertainties(0)
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian) ----!
 pdf = 2
 call setup_pdf(pdf) 
!---- Pass the Higgs mass uncertainty to HiggsSignals ----!
 call HiggsSignals_neutral_input_MassUncertainty(dm)
!---- Set number of free model parameters ----!
 call setup_Nparam(2)

 do i=1,31
  do j=1,21
   scale_ggh = (i-1)*0.05D0
   scale_bbh = (j-1)*0.10D0

   SMGammaTotal=SMGamma_h(Mh)

! SMGamma_h(Mh), SMBR_Hgg(Mh), SMBR_Hgg(Mh) are set to -1 if called
! with Mh out of range [0.8 GeV, 500 GeV]. The calculation is then bypassed.
   if(.not. (SMGammaTotal .lt. 0)) then
    g2hjss_s=1.0d0
    g2hjss_p=0.0d0
    g2hjcc_s=1.0d0
    g2hjcc_p=0.0d0
    g2hjbb_s=scale_bbh**2
    g2hjbb_p=0.0d0
    g2hjtt_s=1.0d0
    g2hjtt_p=0.0d0         
    g2hjmumu_s=1.0d0
    g2hjmumu_p=0.0d0  
    g2hjtautau_s=1.0d0
    g2hjtautau_p=0.0d0
    g2hjWW=1.0d0
    g2hjZZ=1.0d0
    g2hjZga=1d0
    g2hjgg=scale_ggh**2
    g2hjggZ=1d0
    g2hjhiZ=0d0
    g2hjgaga=get_g2hgaga(scale_bbh,sqrt(g2hjtt_s), &
&                        sqrt(g2hjtautau_s),sqrt(g2hjWW),sqrt(g2hjZZ))
!--n.B.: get_g2hgaga can also take negative values as arguments. Here we are interested
!		 only in the positive sector, so we scan over the squared couplings and take the
!		 sqrt here.
    BR_hjhihi=0d0
    BR_hjinvisible=0d0
      
!----Calculate the new total decay width:
    GammaTotal = SMGammaTotal*(1 + &
	&  	            (g2hjWW - 1)*SMBR_HWW(Mh)+(g2hjZZ - 1)*SMBR_HZZ(Mh) + &
	&               (g2hjgg - 1)*SMBR_Hgg(Mh)+(g2hjtt_s - 1)*SMBR_Htoptop(Mh)+ &
	&               (g2hjbb_s - 1)*SMBR_Hbb(Mh)+(g2hjtautau_s - 1)*SMBR_Htautau(Mh)+ &
	&               (g2hjss_s - 1)*SMBR_Hss(Mh)+(g2hjcc_s - 1)*SMBR_Hcc(Mh)+ &
	&               (g2hjZga - 1)*SMBR_HZgam(Mh)+(g2hjmumu_s - 1)*SMBR_Hmumu(Mh)+ &		
	&               (g2hjgaga - 1)*SMBR_Hgamgam(Mh)	)

! Calculate theoretical uncertainties of singleH production from ggh and bbh effective couplings. 
    dCS(1) = get_singleH_uncertainty(dggh, dbbh, g2hjgg, g2hjbb_s+g2hjbb_p, mh)	

    call setup_rate_uncertainties(dCS, dBR)

    call HiggsBounds_neutral_input_effC(Mh,GammaTotal, &
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p, &
     &    g2hjtt_s,g2hjtt_p, &
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p, &
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ, &
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

    call run_HiggsSignals( 1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

! This will collect the main HiggsSignals results together into one file
	call get_rates(1,3,1,(/10/),singleH_rate)
	call get_rates(1,3,1,(/15/),Htobb_rate)

    write(21,*) mh,scale_ggh,scale_bbh,dCS(1),Chisq_mu,Chisq,ndf,singleH_rate,Htobb_rate

   endif
  enddo
 enddo
 
 close(21)

 write(*,*) "Finishing HiggsSignals..."
 call finish_HiggsSignals

contains

!************************************************************** 
 function get_g2hgaga(ghbb,ghtt,ghtautau,ghWW,ghZZ)
! Evaluates g2hgaga from other effective couplings, using partial widths informations
! at a Higgs mass of 126 GeV (calculated with HDECAY and taken from
! http://people.web.psi.ch/spira/higgscoup/ ).
!**************************************************************
 double precision, intent(in) :: ghbb,ghtt,ghtautau,ghWW,ghZZ
 double precision :: get_g2hgaga
  
 get_g2hgaga = (ghtt**2)*0.70904D-01 + (ghbb**2)*0.18760D-04 + (ghWW**2)*1.5863 + &
 & ghtt*ghbb*(-0.17319D-02) + ghtt*ghWW*(-0.67074) + &
 & ghbb*ghWW*0.82093D-02 + (ghtautau**2)*0.22663E-04 + &
 & ghtt*ghtautau*(-0.18696E-02) + ghbb*ghtautau*0.41239E-04 +&
 & ghtautau*ghWW*0.88634E-02
 
 end function get_g2hgaga
!************************************************************** 
 function get_singleH_uncertainty(dggh, dbbh, g2hgg, g2hbb, mh)
!************************************************************** 
 double precision, intent(in) :: dggh, dbbh, g2hgg, g2hbb, mh
 double precision :: get_singleH_uncertainty
 
 if(g2hgg.le.vsmall.and.g2hbb.le.vsmall) then
  get_singleH_uncertainty = 0.0D0
 else 
  get_singleH_uncertainty = ( g2hgg*LHC8_rH_gg(mh)*dggh + g2hbb*LHC8_rH_bb(mh)*dbbh )/ &
 & 	                        ( g2hgg*LHC8_rH_gg(mh)      + g2hbb*LHC8_rH_bb(mh)      )
 endif
 
 end function get_singleH_uncertainty
!************************************************************** 
 end program HSeffC
