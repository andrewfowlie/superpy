!******************************************************
! This example program is part of HiggsSignals (TS 05/03/2013).
!******************************************************
program HShadr
! This example program uses the hadronic cross section input format
! to scan over the two scale factors:
!   scale_ggf scales the SM single Higgs and ttH production cross section
!   scale_VH scales the SM VBF, HZ and HW production cross section
!******************************************************
 implicit none

 integer :: nHzero,nHplus,ndf, i, j, k, ii, jj, CP_value
 double precision :: Pvalue,Chisq,Chisq_mu,Chisq_mh,scale_ggf,scale_VH,dCS(5),dBR(5)
 double precision :: SMGammaTotal, SMGamma_h
 double precision :: SMBR_Htoptop,SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Hmumu, &
&     SMBR_Htautau, SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg
 double precision :: Mh,GammaTotal,CS_lep_hjZ_ratio,                           &
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,    &   
     &          CS_lep_hjhi_ratio_nHbynH,                   &
     &          CS_tev_hj_ratio ,CS_tev_hjb_ratio,    		&
     &          CS_tev_hjW_ratio,CS_tev_hjZ_ratio,    		&
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,   		&
     &          CS_lhc7_hj_ratio ,CS_lhc7_hjb_ratio,    	&
     &          CS_lhc7_hjW_ratio,CS_lhc7_hjZ_ratio,    	&
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,   	&
     &          CS_lhc8_hj_ratio ,CS_lhc8_hjb_ratio,    	&
     &          CS_lhc8_hjW_ratio,CS_lhc8_hjZ_ratio,    	&
     &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,   	&
     &          BR_hjss,BR_hjcc,                            &
     &          BR_hjbb,                                    &
     &          BR_hjmumu,                                  &
     &          BR_hjtautau,                                &
     &          BR_hjWW,BR_hjZZ,BR_hjZga,BR_hjgaga,         &
     &          BR_hjgg, BR_hjinvisible,                    &
     &          BR_hjhihi_nHbynH                            
 character(len=100)::filename
 double precision :: dm
 
 nHzero=1
 nHplus=0
           
 dm = 0.0D0

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
 call initialize_HiggsSignals(nHzero,nHplus,"latestresults")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) 	 ----!
 call setup_output_level(0)
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)	 	 ----!
 call setup_pdf(2) 
  !---- Set the assignment range for the peak-centered method (optional)		     ----! 
 call setup_assignmentrange(10.0D0)
!---- Pass the Higgs mass uncertainty to HiggsSignals							 	 ----!
 call HiggsSignals_neutral_input_MassUncertainty(dm)
!---- Open output text file                                                          ----!
 filename='results/HShadr.dat'
 open(21,file=filename)
 write(21,*) '# Mh	scale_ggf	scale_VH	Chisq_mu	Chisq_mh	Chisq	ndf    Pvalue'
 write(21,*) '#--------------------------------------------------------------------------'

 do i=1,21
  do j=1,21
   scale_ggf = 0.0D0 +(i-1)*0.1D0
   scale_VH = 0.0D0 +(j-1)*0.1D0

! do i=22,24
!  do j=14,15
!   scale_ggf = 0.0D0 +i*0.05D0
!   scale_VH = 0.0D0 +j*0.05D0

   
   Mh=125.7D0
   SMGammaTotal=SMGamma_h(Mh)

! SMGamma_h(Mh), SMBR_Hgg(Mh), SMBR_Hgg(Mh) are set to -1 if called
! with Mh out of range [0.8 GeV, 500 GeV]. The calculation is then bypassed.
   if(.not. (SMGammaTotal .lt. 0)) then
      GammaTotal=SMGammaTotal
! CP even
      CP_value=1
	  CS_lep_hjZ_ratio=1.0D0
	  CS_lep_bbhj_ratio=1.0D0
	  CS_lep_tautauhj_ratio=1.0D0
	  CS_lep_hjhi_ratio_nHbynH=0.0D0
	  CS_tev_hj_ratio=1.0D0*scale_ggf
	  CS_tev_hjb_ratio=1.0D0
	  CS_tev_hjW_ratio=1.0D0*scale_VH
	  CS_tev_hjZ_ratio=1.0D0*scale_VH
	  CS_tev_vbf_ratio=1.0D0*scale_VH
	  CS_tev_tthj_ratio=1.0D0*scale_ggf
 	  CS_lhc7_hj_ratio=1.0D0*scale_ggf
 	  CS_lhc7_hjb_ratio=1.0D0
      CS_lhc7_hjW_ratio=1.0D0*scale_VH
      CS_lhc7_hjZ_ratio=1.0D0*scale_VH
	  CS_lhc7_vbf_ratio=1.0D0*scale_VH
	  CS_lhc7_tthj_ratio=1.0D0*scale_ggf
	  CS_lhc8_hj_ratio=1.0D0*scale_ggf
	  CS_lhc8_hjb_ratio=1.0D0
	  CS_lhc8_hjW_ratio=1.0D0*scale_VH
	  CS_lhc8_hjZ_ratio=1.0D0*scale_VH
	  CS_lhc8_vbf_ratio=1.0D0*scale_VH
	  CS_lhc8_tthj_ratio=1.0D0*scale_ggf
	  BR_hjss=SMBR_Hss(Mh)
	  BR_hjcc=SMBR_Hcc(Mh)
	  BR_hjbb=SMBR_Hbb(Mh)  
	  BR_hjmumu=SMBR_Hmumu(Mh)
	  BR_hjtautau=SMBR_Htautau(Mh)
      BR_hjWW=SMBR_HWW(Mh)
      BR_hjZZ=SMBR_HZZ(Mh)
      BR_hjZga=SMBR_HZgam(Mh)
      BR_hjgaga=SMBR_Hgamgam(Mh)
	  BR_hjgg=SMBR_Hgg(Mh)
	  BR_hjinvisible=0.0D0
      BR_hjhihi_nHbynH=0.0D0              	  

	  call HiggsBounds_neutral_input_hadr(Mh,GammaTotal,CP_value,      &
     &          CS_lep_hjZ_ratio,                           &
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,    &   
     &          CS_lep_hjhi_ratio_nHbynH,                   &
     &          CS_tev_hj_ratio ,CS_tev_hjb_ratio,    		&
     &          CS_tev_hjW_ratio,CS_tev_hjZ_ratio,    		&
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,   		&
     &          CS_lhc7_hj_ratio ,CS_lhc7_hjb_ratio,    	&
     &          CS_lhc7_hjW_ratio,CS_lhc7_hjZ_ratio,    	&
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,   	&
     &          CS_lhc8_hj_ratio ,CS_lhc8_hjb_ratio,    	&
     &          CS_lhc8_hjW_ratio,CS_lhc8_hjZ_ratio,    	&
     &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,   	&
     &          BR_hjss,BR_hjcc,                            &
     &          BR_hjbb,                                    &
     &          BR_hjmumu,                                  &
     &          BR_hjtautau,                                &
     &          BR_hjWW,BR_hjZZ,BR_hjZga,BR_hjgaga,         &
     &          BR_hjgg, BR_hjinvisible,                    &
     &          BR_hjhihi_nHbynH                            )

     call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

	 write(21,*) Mh,scale_ggf,scale_VH,Chisq_mu,Chisq_mh,Chisq,ndf,Pvalue

!	 write(*,*) Mh,scale_ggf,scale_VH,Chisq_mu,Chisq_mh,Chisq,ndf,Pvalue		  
  endif
 enddo
enddo

close(21)

call finish_HiggsSignals

end program HShadr
