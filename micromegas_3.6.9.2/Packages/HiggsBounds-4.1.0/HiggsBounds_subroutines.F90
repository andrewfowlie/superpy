! This file is part of HiggsBounds
!  -KW
!************************************************************      
subroutine initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses_in)
! This the first Higgsbounds subroutine that should be called
! by the user.
! It calls subroutines to read in the tables of Standard Model data,
! read in the tables of LEP, Tevatron and LHC data,
! set up lists of processes which should be checked against 
! the experimental results, allocate arrays etc
! Arguments (input):
!   * nHiggs= number of neutral Higgs in the model 
!     (see subroutine check_nH_nHplus in input.f90 for more details)
!   * nHiggsplus= number of singly,positively charged Higgs in the model 
!     (see subroutine check_nH_nHplus in input.f90 for more details)
!   * whichanalyses_in= which combination of experimental results to use
!     (see subroutine check_whichanalyses in input.f90 for more details)
!************************************************************
 use usefulbits, only : np,Hneut,Hplus,Chineut,Chiplus,debug,inputmethod,  &
  &   inputsub,theo,whichanalyses,HiggsBounds_info,just_after_run,&
  &   file_id_debug1,file_id_debug2,allocate_if_stats_required,run_HB_classic
 use input, only : setup_input,check_number_of_particles,check_whichanalyses
 use S95tables, only : setup_S95tables,S95_t2
 use theory_BRfunctions, only : setup_BRSM
 use channels, only : setup_channels
 use output, only : setup_output
#ifdef enableCHISQ
 use S95tables_type3, only : clsb_t3,fillt3needs_M2_gt_2M1
#endif

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

!#define FORFITTINO

 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
! integer,intent(in),optional :: nHiggsplus
! character(LEN=5),intent(in),optional :: whichanalyses_in
 integer,intent(in) :: nHiggsplus
 character(LEN=5),intent(in) :: whichanalyses_in
 !-----------------------------------internal
 integer :: i
 logical :: messages 
 !-------------------------------------------

! if((.not.present(nHiggsplus)).or.(.not.present(whichanalyses_in)))then 
   !Actually, this doesn't work as I wanted it to
   !because if initialize_HiggsBounds is called in the old way, the program
   !usually just crashes..... but leaving it in for now, in case
   !some compilers accept it
!   call attempting_to_use_an_old_HB_version('init')
! endif

#ifdef FORFITTINO
 write(*,*)'The arguments passed to initialize_HiggsBounds are:'
 write(*,*)'nHiggsneut=',nHiggsneut
 write(*,*)'nHiggsplus=',nHiggsplus
 write(*,*)'whichanalyses_in=','~'//trim(adjustl(whichanalyses_in))//'~'
#endif

#ifdef DEBUGGING
  debug=.True.
#else
  debug=.False.
#endif

 messages=debug.or.(inputmethod=='datfile')

! inputmethod='subrout' !('datfile' or 'website' are also possible, but not here)

 np(Hneut)=nHiggsneut
 np(Hplus)=nHiggsplus
 
 np(Chineut)=0! do not change this without contacting us first!
 np(Chiplus)=0! do not change this without contacting us first!

 whichanalyses=whichanalyses_in
 
 if(inputmethod=='subrout') then
  if(allocated(theo))then
   stop 'subroutine HiggsBounds_initialize has already been called once'
  endif

  if(messages)write(*,*)'doing other preliminary tasks...'      ; call flush(6)
  call setup_input

  allocate(inputsub( 4 )) !(1)np(Hneut)>0 (2)np(Hplus)>0 (3)np(Chineut)>0 (4)np(Chineut)>0 and np(Chiplus)>0
   !                                                                           |        np
   !                                                                           |Hneu Hcha Chineut Chiplus
   !                                                                           | ==0  ==0  ==0     ==0 
  inputsub(1)%desc='HiggsBounds_neutral_input_*'
  inputsub(1)%req=req(   0,   1,   1,   1)
  inputsub(2)%desc='HiggsBounds_charged_input'
  inputsub(2)%req=req(   1,   0,   1,   1)
  inputsub(3)%desc='SUSYBounds_neutralinoonly_input'
  inputsub(3)%req=req(   1,   1,   0,   1)
  inputsub(4)%desc='SUSYBounds_neutralinochargino_input'
  inputsub(4)%req=req(   1,   1,   0,   0)
  do i=1,ubound(inputsub,dim=1)
   inputsub(i)%stat=0
  enddo
 endif 

#ifndef WEBVERSION
 if(inputmethod.ne.'datfile') call HiggsBounds_info
 if (run_HB_classic.EQV..True.) then
  PRINT *, "run_HB_classic=True - HiggsBounds is running in classic mode"
 endif
  
#endif

 if(messages)write(*,*)'reading in Standard Model tables...'   ; call flush(6) 
 call setup_BRSM 
             
 if(messages)write(*,*)'reading in S95tables...'               ; call flush(6)
 call setup_S95tables                                    

!if(debug)write(*,*)'doing other preliminary tasks...'      ; call flush(6)
!call setup_input

 if(messages)then
  open(file_id_debug2,file='debug_predratio.txt')
  open(file_id_debug1,file='debug_channels.txt')
 endif

 if(messages)write(*,*)'sorting out processes to be checked...'; call flush(6)
 call setup_channels      
      
 if(messages)write(*,*)'preparing output arrays...'            ; call flush(6)
 call setup_output

#ifdef enableCHISQ
 if(allocated(allocate_if_stats_required))then
   call fillt3needs_M2_gt_2M1(clsb_t3,S95_t2)
 endif 
#endif

 just_after_run=.False.
  
 contains 
   !         |        np
   !         |Hneu Hcha Chineut Chiplus
   !         | ==0  ==0  ==0     ==0 
 function req(Hneu,Hcha, Chneu,  Chcha)
  integer, intent(in) ::Hneu,Hcha, Chneu,  Chcha
  integer :: req
  
  req=1
  if(np(Hneut)==0)  req= Hneu  * req
  if(np(Hplus)==0)  req= Hcha  * req
  if(np(Chineut)==0)req= Chneu * req
  if(np(Chiplus)==0)req= Chcha * req

 end function req 
end subroutine initialize_HiggsBounds
!************************************************************      


!************************************************************      
! Version of initialize_HiggsBounds which takes an integer as
! the third argument. More useful for library linking to 
! non-Fortran codes.
subroutine initialize_HiggsBounds_int(nHn,nHp,flag)
 
       implicit none

       integer nHn,nHp,flag

    interface
        subroutine initialize_HiggsBounds(nHiggsneut, nHiggsplus, whichanalyses_in)
     integer,intent(in) :: nHiggsneut
     integer,intent(in),optional :: nHiggsplus
     character(LEN=5),intent(in),optional :: whichanalyses_in
        end subroutine initialize_HiggsBounds
    end interface 

      IF (flag.EQ.1) then
        call initialize_HiggsBounds(nHn,nHp, "onlyL")
      elseif (flag.EQ.2) then
        call initialize_HiggsBounds(nHn,nHp, "onlyH")
      elseif (flag.EQ.3) then
        call initialize_HiggsBounds(nHn,nHp, "LandH")
      elseif (flag.EQ.4) then
        call initialize_HiggsBounds(nHn,nHp, "onlyP")
      else
        stop "Illegal value for flag in call to initialize_HB"
      endif
      

end subroutine
!************************************************************      



!************************************************************      
subroutine attempting_to_use_an_old_HB_version(subroutineid)
 use usefulbits, only : vers
 character(len=4),intent(in) :: subroutineid

 select case(subroutineid)
 case('init')
  write(*,*)'The subroutine initialize_HiggsBounds has been called with the'
  write(*,*)'wrong number of arguments. It should be called as:'
  write(*,*)'initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses)'
  write(*,*)
  write(*,*)'Note that in early versions of HiggsBounds (HB 1.*.*)'
  write(*,*)'this subroutine was called as:'
  write(*,*)'initialize_HiggsBounds(nHiggsneut,whichanalyses)'
  write(*,*)
 case('effC','part','hadr')
  write(*,*)'The subroutine run_HiggsBounds_'//subroutineid//' has been discontinued in this'
  write(*,*)'version of HiggsBounds.'
 case default
  stop'wrong input to subroutine attempting_to_use_an_old_HB_version'
 end select

 write(*,*)'If you have code written for use with HB 1.*.*, you have two choices:'
 write(*,*)
 write(*,*)' (1) You can edit your code, such that it works with this' 
 write(*,*)'     version of HiggsBounds (HB'//trim(adjustl(vers))//').'
 write(*,*)'     This has the advantage that you can test your model against many, many'
 write(*,*)'     more Higgs search limits , including charged Higgs search limits.'
 write(*,*)'     See the updated manual for more information.'
 write(*,*)
 write(*,*)' (2) You can download the most recent HB 1.*.* from the HiggsBounds'
 write(*,*)'     website. This contains the LEP Higgs search limits which are'
 write(*,*)'     generally the most useful when constraining new physics models.'
 write(*,*)'     We will continue to support this code.'

 stop'Incorrect call to a HiggsBounds subroutine.'

end subroutine attempting_to_use_an_old_HB_version
!************************************************************      
subroutine HiggsBounds_input_SLHA(infile)
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): SLHA filename
!************************************************************
 use usefulbits, only : whichinput,inputsub,infile1,theo,g2,just_after_run, &
   & np,Hneut,Hplus     
 use extra_bits_for_SLHA

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 character(len=100),intent(in) :: infile
 !--------------------------------------internal
 integer :: n
 !----------------------------------------------
 
 whichinput='SLHA'

 if(np(Hneut).gt.0)inputsub(Hneut)%stat=inputsub(Hneut)%stat+1
 if(np(Hplus).gt.0)inputsub(Hplus)%stat=inputsub(Hplus)%stat+1
 ! note: can't be used for charginos or neutralinos yet 

 n=1  

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif
 
 infile1=infile   

 call getSLHAdata(theo(n),g2(n),infile1)    

 just_after_run=.False.

end subroutine HiggsBounds_input_SLHA
!************************************************************     
subroutine HiggsBounds_neutral_input_effC(Mh,GammaTotal_hj,        &  
     &          g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,               &
     &          g2hjbb_s,g2hjbb_p,g2hjtoptop_s,g2hjtoptop_p,       &
     &          g2hjmumu_s,g2hjmumu_p,                             &
     &          g2hjtautau_s,g2hjtautau_p,                         &
     &          g2hjWW,g2hjZZ,g2hjZga,                             &
     &          g2hjgaga,g2hjgg,g2hjggZ,g2hjhiZ_nHbynH,            &
     &          BR_hjinvisible,BR_hjhihi_nHbynH                    )
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): theoretical predictions (see manual for definitions)
!************************************************************
 use usefulbits, only : theo,np,Hneut,g2,whichinput,inputsub,just_after_run

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: Mh( np(Hneut) ),GammaTotal_hj( np(Hneut) ),                                 &  
     &          g2hjss_s( np(Hneut) ),g2hjss_p( np(Hneut) ),g2hjcc_s( np(Hneut) ),g2hjcc_p( np(Hneut) ),        &
     &          g2hjbb_s( np(Hneut) ),g2hjbb_p( np(Hneut) ),g2hjtoptop_s( np(Hneut) ),g2hjtoptop_p( np(Hneut) ),&
     &          g2hjmumu_s( np(Hneut) ),g2hjmumu_p( np(Hneut) ),                                        &
     &          g2hjtautau_s( np(Hneut) ),g2hjtautau_p( np(Hneut) ),                                        &
     &          g2hjWW( np(Hneut) ),g2hjZZ( np(Hneut) ),g2hjZga( np(Hneut) ),                                 &
     &          g2hjgaga( np(Hneut) ),g2hjgg( np(Hneut) ),g2hjggZ( np(Hneut) ),g2hjhiZ_nHbynH(np(Hneut),np(Hneut)),&
     &          BR_hjinvisible( np(Hneut) ),BR_hjhihi_nHbynH(np(Hneut),np(Hneut)) 
 !--------------------------------------internal
 integer :: n
 integer :: subtype
 !----------------------------------------------
 
 whichinput='effC'
 subtype=1
 n=1  
 inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_effC should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_effC'
 endif

 theo(n)%particle(Hneut)%M       = Mh 
 theo(n)%particle(Hneut)%GammaTot= GammaTotal_hj

 g2(n)%hjss_s              = g2hjss_s
 g2(n)%hjss_p              = g2hjss_p
 g2(n)%hjcc_s              = g2hjcc_s
 g2(n)%hjcc_p              = g2hjcc_p
 g2(n)%hjbb_s              = g2hjbb_s
 g2(n)%hjbb_p              = g2hjbb_p
 g2(n)%hjtoptop_s              = g2hjtoptop_s
 g2(n)%hjtoptop_p              = g2hjtoptop_p
 g2(n)%hjmumu_s                = g2hjmumu_s
 g2(n)%hjmumu_p                = g2hjmumu_p
 g2(n)%hjtautau_s              = g2hjtautau_s
 g2(n)%hjtautau_p              = g2hjtautau_p
        
 g2(n)%hjWW              = g2hjWW
 g2(n)%hjZZ              = g2hjZZ  
 g2(n)%hjZga             = g2hjZga      
 g2(n)%hjgaga            = g2hjgaga
 g2(n)%hjgg              = g2hjgg
 g2(n)%hjggZ             = g2hjggZ

 g2(n)%hjhiZ = g2hjhiZ_nHbynH

 theo(n)%BR_hjinvisible  = BR_hjinvisible 
 theo(n)%BR_hjhihi       = BR_hjhihi_nHbynH  

 just_after_run=.False. 

end subroutine HiggsBounds_neutral_input_effC
!************************************************************      
subroutine HiggsBounds_neutral_input_part(Mh,GammaTotal_hj,CP_value,       &
     &          CS_lep_hjZ_ratio,                            &
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,     &
     &          CS_lep_hjhi_ratio_nHbynH,                    &
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,       &
     &          CS_bg_hjb_ratio,                         &
     &          CS_ud_hjWp_ratio,CS_cs_hjWp_ratio,   & 
     &          CS_ud_hjWm_ratio,CS_cs_hjWm_ratio,   & 
     &          CS_gg_hjZ_ratio,     &
     &          CS_dd_hjZ_ratio,CS_uu_hjZ_ratio,     &
     &          CS_ss_hjZ_ratio,CS_cc_hjZ_ratio,     &
     &          CS_bb_hjZ_ratio,                         &
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,    &
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,    &
     &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,    &
     &          BR_hjss,BR_hjcc,                             &
     &          BR_hjbb,BR_hjmumu,BR_hjtautau,               &
     &          BR_hjWW,BR_hjZZ,BR_hjZga, BR_hjgaga,BR_hjgg, & 
     &          BR_hjinvisible,BR_hjhihi_nHbynH              )
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! (see manual for full description)
!************************************************************
 use usefulbits, only : theo,np,Hneut,partR,whichinput,inputsub,just_after_run

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 !----------------------------------------input
 double precision,intent(in) :: Mh( np(Hneut) ),GammaTotal_hj( np(Hneut) )
 integer,intent(in) ::CP_value( np(Hneut) )
 double precision,intent(in) :: CS_lep_hjZ_ratio( np(Hneut) ),                         &
     &          CS_lep_bbhj_ratio( np(Hneut) ),CS_lep_tautauhj_ratio( np(Hneut) ),     &
     &          CS_lep_hjhi_ratio_nHbynH(np(Hneut),np(Hneut)),                         &
     &          CS_gg_hj_ratio( np(Hneut) ),CS_bb_hj_ratio( np(Hneut) ),       &
     &          CS_bg_hjb_ratio( np(Hneut) ),                                      &
     &          CS_ud_hjWp_ratio( np(Hneut) ),CS_cs_hjWp_ratio( np(Hneut) ),   & 
     &          CS_ud_hjWm_ratio( np(Hneut) ),CS_cs_hjWm_ratio( np(Hneut) ),   & 
     &          CS_gg_hjZ_ratio( np(Hneut) ),    &
     &          CS_dd_hjZ_ratio( np(Hneut) ),CS_uu_hjZ_ratio( np(Hneut) ),     &
     &          CS_ss_hjZ_ratio( np(Hneut) ),CS_cc_hjZ_ratio( np(Hneut) ),     &
     &          CS_bb_hjZ_ratio( np(Hneut) ),                                      &
     &          CS_tev_vbf_ratio( np(Hneut) ),CS_tev_tthj_ratio( np(Hneut) ),    &
     &          CS_lhc7_vbf_ratio( np(Hneut) ),CS_lhc7_tthj_ratio( np(Hneut) ),    &
     &          CS_lhc8_vbf_ratio( np(Hneut) ),CS_lhc8_tthj_ratio( np(Hneut) ),    &
     &          BR_hjss( np(Hneut) ),BR_hjcc( np(Hneut) ),                             &
     &          BR_hjbb( np(Hneut) ),BR_hjmumu( np(Hneut) ),BR_hjtautau( np(Hneut) ),  &
     &          BR_hjWW( np(Hneut) ),BR_hjZZ( np(Hneut) ),BR_hjZga( np(Hneut) ),       &
     &          BR_hjgaga( np(Hneut) ),BR_hjgg( np(Hneut) ),                           & 
     &          BR_hjinvisible( np(Hneut) ),BR_hjhihi_nHbynH(np(Hneut),np(Hneut)) 
 !---------------------------------------internal
 integer :: n
 integer :: subtype
 !-----------------------------------------------

 whichinput='part'
 subtype=1
 n=1 
 inputsub(subtype)%stat=inputsub(subtype)%stat+1
      
 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_part should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_part'
 endif

 theo(n)%particle(Hneut)%M = Mh       
 theo(n)%particle(Hneut)%GammaTot= GammaTotal_hj
 theo(n)%CP_value          = CP_value  
 theo(n)%lep%XS_hjZ_ratio       = CS_lep_hjZ_ratio
 theo(n)%lep%XS_bbhj_ratio      = CS_lep_bbhj_ratio
 theo(n)%lep%XS_tautauhj_ratio  = CS_lep_tautauhj_ratio
 theo(n)%lep%XS_hjhi_ratio      = CS_lep_hjhi_ratio_nHbynH 
 partR(n)%gg_hj             = CS_gg_hj_ratio 
 partR(n)%qq_hj(5,:)        = CS_bb_hj_ratio
 partR(n)%bg_hjb            = CS_bg_hjb_ratio                    
 partR(n)%qq_hjWp(1,:)      = CS_ud_hjWp_ratio
 partR(n)%qq_hjWp(2,:)      = CS_cs_hjWp_ratio         
 partR(n)%qq_hjWm(1,:)      = CS_ud_hjWm_ratio
 partR(n)%qq_hjWm(2,:)      = CS_cs_hjWm_ratio  
 partR(n)%gg_hjZ(:)         = CS_gg_hjZ_ratio      
 partR(n)%qq_hjZ(1,:)       = CS_dd_hjZ_ratio
 partR(n)%qq_hjZ(2,:)       = CS_uu_hjZ_ratio           
 partR(n)%qq_hjZ(3,:)       = CS_ss_hjZ_ratio
 partR(n)%qq_hjZ(4,:)       = CS_cc_hjZ_ratio             
 partR(n)%qq_hjZ(5,:)       = CS_bb_hjZ_ratio                         
 theo(n)%tev%XS_vbf_ratio  = CS_tev_vbf_ratio   
 theo(n)%tev%XS_tthj_ratio = CS_tev_tthj_ratio  
 theo(n)%lhc7%XS_vbf_ratio = CS_lhc7_vbf_ratio   
 theo(n)%lhc7%XS_tthj_ratio= CS_lhc7_tthj_ratio
 theo(n)%lhc8%XS_vbf_ratio = CS_lhc8_vbf_ratio   
 theo(n)%lhc8%XS_tthj_ratio= CS_lhc8_tthj_ratio
 theo(n)%BR_hjss           = BR_hjss  
 theo(n)%BR_hjcc           = BR_hjcc                
 theo(n)%BR_hjbb           = BR_hjbb
 theo(n)%BR_hjmumu         = BR_hjmumu
 theo(n)%BR_hjtautau       = BR_hjtautau              
 theo(n)%BR_hjWW           = BR_hjWW
 theo(n)%BR_hjZZ           = BR_hjZZ
 theo(n)%BR_hjZga          = BR_hjZga  
 theo(n)%BR_hjgaga         = BR_hjgaga
 theo(n)%BR_hjgg           = BR_hjgg  
 theo(n)%BR_hjinvisible    = BR_hjinvisible             
 theo(n)%BR_hjhihi         = BR_hjhihi_nHbynH  

 just_after_run=.False. 

end subroutine HiggsBounds_neutral_input_part
!************************************************************      
subroutine HiggsBounds_neutral_input_hadr(Mh,GammaTotal_hj,CP_value,      &
     &          CS_lep_hjZ_ratio,                           &
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,    &   
     &          CS_lep_hjhi_ratio_nHbynH,                   &
     &          CS_tev_hj_ratio ,CS_tev_hjb_ratio,    &
     &          CS_tev_hjW_ratio,CS_tev_hjZ_ratio,    &
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,   &
     &          CS_lhc7_hj_ratio ,CS_lhc7_hjb_ratio,    &
     &          CS_lhc7_hjW_ratio,CS_lhc7_hjZ_ratio,    &
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,   &
     &          CS_lhc8_hj_ratio ,CS_lhc8_hjb_ratio,    &
     &          CS_lhc8_hjW_ratio,CS_lhc8_hjZ_ratio,    &
     &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,   &
     &          BR_hjss,BR_hjcc,                            &
     &          BR_hjbb,                                    &
     &          BR_hjmumu,                                  &
     &          BR_hjtautau,                                &
     &          BR_hjWW,BR_hjZZ,BR_hjZga,BR_hjgaga,         &
     &          BR_hjgg, BR_hjinvisible,                    &
     &          BR_hjhihi_nHbynH                            )
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! (see manual for full description)
!************************************************************
 use usefulbits, only : theo,np,Hneut,whichinput,inputsub,just_after_run

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 !----------------------------------------input
 double precision,intent(in) :: Mh( np(Hneut) ),GammaTotal_hj( np(Hneut) )
 integer,intent(in) :: CP_value( np(Hneut) )
 double precision,intent(in) :: CS_lep_hjZ_ratio( np(Hneut) ),                        &
     &          CS_lep_bbhj_ratio( np(Hneut) ),CS_lep_tautauhj_ratio( np(Hneut) ),    &   
     &          CS_lep_hjhi_ratio_nHbynH(np(Hneut),np(Hneut)),                        &
     &          CS_tev_hj_ratio( np(Hneut)  ) ,CS_tev_hjb_ratio( np(Hneut) ),    &
     &          CS_tev_hjW_ratio( np(Hneut) ) ,CS_tev_hjZ_ratio( np(Hneut) ),    &
     &          CS_tev_vbf_ratio( np(Hneut) ) ,CS_tev_tthj_ratio( np(Hneut)),    &
     &          CS_lhc7_hj_ratio( np(Hneut)  ),CS_lhc7_hjb_ratio( np(Hneut) ),    &
     &          CS_lhc7_hjW_ratio( np(Hneut) ),CS_lhc7_hjZ_ratio( np(Hneut) ),    &
     &          CS_lhc7_vbf_ratio( np(Hneut) ),CS_lhc7_tthj_ratio( np(Hneut)),    &
     &          CS_lhc8_hj_ratio( np(Hneut)  ),CS_lhc8_hjb_ratio( np(Hneut) ),    &
     &          CS_lhc8_hjW_ratio( np(Hneut) ),CS_lhc8_hjZ_ratio( np(Hneut) ),    &
     &          CS_lhc8_vbf_ratio( np(Hneut) ),CS_lhc8_tthj_ratio( np(Hneut)),    &
     &          BR_hjss( np(Hneut) ),BR_hjcc( np(Hneut) ),                            &
     &          BR_hjbb( np(Hneut) ),                                                 &
     &          BR_hjmumu( np(Hneut) ),BR_hjtautau( np(Hneut) ),                      &
     &          BR_hjWW( np(Hneut) ),BR_hjZZ( np(Hneut) ),                            &
     &          BR_hjZga( np(Hneut) ),BR_hjgaga( np(Hneut) ),                         &
     &          BR_hjgg( np(Hneut) ), BR_hjinvisible( np(Hneut) ),                    &
     &          BR_hjhihi_nHbynH(np(Hneut),np(Hneut))
 !-------------------------------------internal
 integer :: n
 integer :: subtype
 !---------------------------------------------
   
 whichinput='hadr'
 subtype=1
 n=1 
 inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_hadr should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_hadr'
 endif

 theo(n)%particle(Hneut)%M       = Mh  
 theo(n)%particle(Hneut)%GammaTot= GammaTotal_hj
 theo(n)%CP_value                = CP_value             
 theo(n)%lep%XS_hjZ_ratio        = CS_lep_hjZ_ratio
 theo(n)%lep%XS_bbhj_ratio       = CS_lep_bbhj_ratio
 theo(n)%lep%XS_tautauhj_ratio   = CS_lep_tautauhj_ratio
 theo(n)%lep%XS_hjhi_ratio       = CS_lep_hjhi_ratio_nHbynH 
 theo(n)%tev%XS_hj_ratio         = CS_tev_hj_ratio
 theo(n)%tev%XS_hjb_ratio        = CS_tev_hjb_ratio   
 theo(n)%tev%XS_hjW_ratio        = CS_tev_hjW_ratio
 theo(n)%tev%XS_hjZ_ratio        = CS_tev_hjZ_ratio    
 theo(n)%tev%XS_vbf_ratio        = CS_tev_vbf_ratio      
 theo(n)%tev%XS_tthj_ratio       = CS_tev_tthj_ratio
 theo(n)%lhc7%XS_hj_ratio        = CS_lhc7_hj_ratio
 theo(n)%lhc7%XS_hjb_ratio       = CS_lhc7_hjb_ratio   
 theo(n)%lhc7%XS_hjW_ratio       = CS_lhc7_hjW_ratio
 theo(n)%lhc7%XS_hjZ_ratio       = CS_lhc7_hjZ_ratio    
 theo(n)%lhc7%XS_vbf_ratio       = CS_lhc7_vbf_ratio      
 theo(n)%lhc7%XS_tthj_ratio      = CS_lhc7_tthj_ratio
 theo(n)%lhc8%XS_hj_ratio        = CS_lhc8_hj_ratio
 theo(n)%lhc8%XS_hjb_ratio       = CS_lhc8_hjb_ratio   
 theo(n)%lhc8%XS_hjW_ratio       = CS_lhc8_hjW_ratio
 theo(n)%lhc8%XS_hjZ_ratio       = CS_lhc8_hjZ_ratio    
 theo(n)%lhc8%XS_vbf_ratio       = CS_lhc8_vbf_ratio      
 theo(n)%lhc8%XS_tthj_ratio      = CS_lhc8_tthj_ratio
 theo(n)%BR_hjss           = BR_hjss
 theo(n)%BR_hjcc           = BR_hjcc 
 theo(n)%BR_hjbb           = BR_hjbb
 theo(n)%BR_hjmumu         = BR_hjmumu
 theo(n)%BR_hjtautau       = BR_hjtautau                 
 theo(n)%BR_hjWW           = BR_hjWW
 theo(n)%BR_hjZZ           = BR_hjZZ
 theo(n)%BR_hjZga          = BR_hjZga 
 theo(n)%BR_hjgaga         = BR_hjgaga
 theo(n)%BR_hjgg           = BR_hjgg
 theo(n)%BR_hjinvisible    = BR_hjinvisible                  
 theo(n)%BR_hjhihi         = BR_hjhihi_nHbynH  

 just_after_run=.False. 
    
end subroutine HiggsBounds_neutral_input_hadr
!************************************************************      
subroutine HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj, &
     &          CS_lep_HpjHmj_ratio,                        &
     &          BR_tWpb,BR_tHpjb,                           &
     &          BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu)
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): theoretical predictions (see manual for definitions)
!************************************************************
 use usefulbits, only : theo,np,Hplus,inputsub,just_after_run

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: Mhplus( np(Hplus) ),GammaTotal_Hpj( np(Hplus) ), &
     &          CS_lep_HpjHmj_ratio( np(Hplus) ),                              &
     &          BR_tWpb,BR_tHpjb( np(Hplus) ),                                 &
     &          BR_Hpjcs( np(Hplus) ),BR_Hpjcb( np(Hplus) ),BR_Hpjtaunu( np(Hplus) ) 
 !--------------------------------------internal
 integer :: n
 integer :: subtype
 !----------------------------------------------
 
 n=1  
 subtype=2
 inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hplus).eq.0)then
  write(*,*)'subroutine HiggsBounds_charged_input should'
  write(*,*)'only be called if np(Hplus)>0'
  stop 'error in subroutine HiggsBounds_charged_input'
 endif

 theo(n)%particle(Hplus)%M       = Mhplus 
 theo(n)%particle(Hplus)%GammaTot= GammaTotal_Hpj

 theo(n)%lep%XS_HpjHmj_ratio = CS_lep_HpjHmj_ratio

 theo(n)%BR_tWpb  = BR_tWpb
 theo(n)%BR_tHpjb = BR_tHpjb 

 theo(n)%BR_Hpjcs     = BR_Hpjcs
 theo(n)%BR_Hpjcb     = BR_Hpjcb
 theo(n)%BR_Hpjtaunu  = BR_Hpjtaunu
  

 just_after_run=.False. 

end subroutine HiggsBounds_charged_input
!************************************************************      
subroutine HiggsBounds_set_mass_uncertainties(dMhneut, dMhch)
! Assigns the mass uncertainties in the subroutine version.
!
   use usefulbits, only : theo,np,Hneut,Hplus
   
  implicit none
  
  double precision, intent(in) :: dMhneut(np(Hneut))
  double precision, intent(in) :: dMhch(np(Hplus))
    
  theo(1)%particle(Hneut)%dMh = dMhneut
  theo(1)%particle(Hplus)%dMh = dMhch 
  
end subroutine HiggsBounds_set_mass_uncertainties
!************************************************************      
subroutine get_mass_variation_param(n)

  use usefulbits, only : theo,np,Hneut,Hplus,diffMhneut,diffMhch,ndmh,dmhsteps,small_mh
  implicit none

  integer, intent(in) :: n

  double precision :: dMhneut(np(Hneut))
  double precision :: dMhch(np(Hplus))      
  integer :: km(np(Hneut)+np(Hplus))
  integer :: dm(dmhsteps**(np(Hneut)+np(Hplus)),np(Hneut)+np(Hplus))

  integer i,j,k,kp

  if(np(Hneut).gt.0) dMhneut = theo(n)%particle(Hneut)%dMh
  if(np(Hplus).gt.0) dMhch = theo(n)%particle(Hplus)%dMh

  if (modulo(dmhsteps,2).NE.1) then
     stop 'Wrong number of steps in set_mass_uncertainty: must be odd (>=3)'
  endif  

  ndmh = 0
  do i=1,np(Hneut)
   
   IF (dMhneut(i).GT.small_mh) THEN
    ndmh = ndmh + 1
   ENDIF
   km(i)=-(dmhsteps-1)/2
  enddo
   
  do i=1,np(Hplus)
   IF (dMhch(i).GT.small_mh) ndmh = ndmh + 1
    km(i+np(Hneut))=-(dmhsteps-1)/2
  enddo
  
  IF (ndmh.EQ.0) THEN
    RETURN
  ENDIF
    
!  print *, "Number of mass uncertainties: ", ndmh

  if(allocated(diffMhneut)) deallocate(diffMhneut)
  if(allocated(diffMhch)) deallocate(diffMhch)
    
  allocate(diffMhneut(dmhsteps**(np(Hneut)+np(Hplus)),np(Hneut)))
  allocate(diffMhch(dmhsteps**(np(Hneut)+np(Hplus)),np(Hplus)))

  k = 1
  do i=1,dmhsteps**ndmh

    do j=1,ndmh
     dm(i,j) = km(j)
    enddo

    km(k) = km(k)+1

    do j=2,ndmh
     IF (modulo(i,dmhsteps**(j-1)).EQ.0) THEN
      km(j) = km(j)+1
      km(j-1) = -1
     ENDIF
    ENDDO

  enddo    
  
  
  do i=1,dmhsteps**ndmh
   k=1
 
    do j=1,np(Hneut)
     IF (dMhneut(j).GT.small_mh) THEN
      diffMhneut(i,j)=theo(n)%particle(Hneut)%M(j)+dm(i,k)*dMhneut(k)/((dmhsteps-1)/2)
      k = k +1
     ELSE
      diffMhneut(i,j)=theo(n)%particle(Hneut)%M(j)
     ENDIF
    enddo
    
    kp = k
    do j=1,np(Hplus)
     IF (dMhch(j).GT.small_mh) THEN
      diffMhch(i,j)=theo(n)%particle(Hplus)%M(j)+dm(i,k)*dMhch(k-(kp-1))/((dmhsteps-1)/2)
      k = k +1
     ELSE
      diffMhch(i,j)=theo(n)%particle(Hplus)%M(j)
     ENDIF
    enddo
    
  
!    print *, i, (diffMhneut(i,j),j=1,np(Hneut)),(diffMhch(i,j),j=1,np(Hplus))
  
  enddo
  
end subroutine get_mass_variation_param


subroutine SUSYBounds_neutralinoonly_input(MN,GammaTotal_N, &
     &          CS_NjNi,                        &
     &          BR_NjqqNi,BR_NjZNi                           &
     &          )
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): theoretical predictions (see manual for definitions)
!************************************************************
 use usefulbits, only : theo,np,Chineut,inputsub,just_after_run

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: MN( np(Chineut) ),GammaTotal_N( np(Chineut) ) , &
     &          CS_NjNi( np(Chineut),np(Chineut) ),                        &
     &          BR_NjqqNi( np(Chineut),np(Chineut) ),BR_NjZNi( np(Chineut),np(Chineut) )                        
 !--------------------------------------internal
 integer :: n
 integer :: subtype
 !----------------------------------------------
 
 n=1  
 subtype=3
 inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Chineut).eq.0)then
   write(*,*)'subroutine SUSYBounds_neutralinoonly_input should'
   write(*,*)'only be called if np(Chineut)>0'
   stop'error in SUSYBounds_neutralinoonly_input'
 endif

 theo(n)%particle(Chineut)%M       = MN
 theo(n)%particle(Chineut)%GammaTot= GammaTotal_N

 theo(n)%lep%XS_NjNi = CS_NjNi

 theo(n)%BR_NjqqNi  = BR_NjqqNi
 theo(n)%BR_NjZNi   = BR_NjZNi  

 just_after_run=.False. 

end subroutine SUSYBounds_neutralinoonly_input
!************************************************************      
subroutine SUSYBounds_neutralinochargino_input(MC,GammaTotal_C, &
     &          CS_CpjCmj,                                      &
     &          BR_CjqqNi,                                      &
     &          BR_CjlnuNi,                                     &
     &          BR_CjWNi                                        &
     &          )
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): theoretical predictions (see manual for definitions)
!************************************************************
 use usefulbits, only : theo,np,Chineut,Chiplus,inputsub,just_after_run

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: MC( np(Chiplus) ),GammaTotal_C( np(Chiplus) ),  &
     &          CS_CpjCmj(  np(Chiplus) ),                                     &
     &          BR_CjqqNi(  np(Chiplus),np(Chineut) ),                         &
     &          BR_CjlnuNi( np(Chiplus),np(Chineut) ),                         &
     &          BR_CjWNi(   np(Chiplus),np(Chineut) )                                        
 !--------------------------------------internal
 integer :: n
 integer :: subtype
 !----------------------------------------------
 
 n=1  
 subtype=4
 inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if((np(Chineut).eq.0).or.(np(Chiplus).eq.0))then
   write(*,*)'subroutine SUSYBounds_neutralinochargino_input should'
   write(*,*)'only be called if np(Chineut)>0 and np(Chiplus)>0'
   stop 'error in subroutine SUSYBounds_neutralinochargino_input'
 endif

 theo(n)%particle(Chineut)%M       = MC
 theo(n)%particle(Chineut)%GammaTot= GammaTotal_C

 theo(n)%lep%XS_CpjCmj = CS_CpjCmj

 theo(n)%BR_CjqqNi  = BR_CjqqNi
 theo(n)%BR_CjlnuNi = BR_CjlnuNi
 theo(n)%BR_CjWNi   = BR_CjWNi  

 just_after_run=.False. 

end subroutine SUSYBounds_neutralinochargino_input
!************************************************************


subroutine run_HiggsBounds(HBresult, chan, obsratio, ncombined)
! This subroutine can be called by the user after HiggsBounds_initialize has been called.
! The input routines, where required, should be called before calling run_HiggsBounds.
! It takes theoretical predictions for a particular parameter point 
! in the model and calls subroutines which compare these predictions 
! to the experimental limits
! Arguments (output): 
!   * HBresult = 1 if point is unexcluded, 0 if excluded, -1 if parameter point is invalid
!   * chan = number of channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!   * obsratio = ratio of the theoretical rate to the observed limit for this channel
!   * ncombined = number of Higgs combined in order to calculate this obsratio
!    (see manual for more precise definitions))
! (TS 30/01/2012): Note, that if many data points are tested at the same time (as for
!                  inputmethod==datfiles), this subroutine only returns the results of
!                  the last datapoint. The full results are saved in fullHBres.
 use usefulbits, only : np, Hneut, Hplus, run_HB_classic

 implicit none
 integer HBresult, chan, ncombined
 double precision obsratio
 
 integer hbres(0:np(Hneut)+np(Hplus)), hbchan(0:np(Hneut)+np(Hplus)), hbcomb(0:np(Hneut)+np(Hplus))
 double precision hbobs(0:np(Hneut)+np(Hplus))

! Check if we are using the old 'classic' method
 if (run_HB_classic.EQV..True.) then
   call run_HiggsBounds_classic(HBresult,chan,obsratio,ncombined)
   return
 endif
 
! Call the new ('full') method
 call run_HiggsBounds_full(hbres, hbchan, hbobs, hbcomb)
 
! Combined results are contained in the zero elements of result arrays
 HBresult = hbres(0)
 chan = hbchan(0)
 obsratio = hbobs(0)
 ncombined = hbcomb(0)
    
end subroutine run_HiggsBounds
!************************************************************



subroutine run_HiggsBounds_single(h, HBresult, chan, obsratio, ncombined)
! This subroutine can be used to get the exclusion results
! for a single Higgs boson (specified by the index h).
!
! To obtain individual results from more than one Higgs boson, it
! is more efficient to use run_HiggsBounds_full rather than this method.
 use usefulbits, only : np, Hneut, Hplus

 implicit none
 integer, intent(in) :: h
 integer, intent(out) :: HBresult, chan, ncombined
 double precision, intent(out) :: obsratio
 
 integer hbres(0:np(Hneut)+np(Hplus)), hbchan(0:np(Hneut)+np(Hplus)), hbcomb(0:np(Hneut)+np(Hplus))
 double precision hbobs(0:np(Hneut)+np(Hplus))
 
 IF (h.LT.0) stop "Illegal number of Higgs boson: h < 0"
 if (h.GT.np(Hneut)+np(Hplus)) stop "Illegal number of Higgs boson"
 
 call run_HiggsBounds_full(hbres, hbchan, hbobs, hbcomb)

 HBresult = hbres(h)
 chan = hbchan(h)
 obsratio = hbobs(h)
 ncombined = hbcomb(h)
end subroutine run_HiggsBounds_single
!************************************************************


subroutine run_HiggsBounds_full( HBresult,chan,                      &
     &                      obsratio, ncombined                 )
! This subroutine can be called by the user after HiggsBounds_initialize has been called.
! The input routines, where required, should be called before calling run_HiggsBounds.
! It takes theoretical predictions for a particular parameter point 
! in the model and calls subroutines which compare these predictions 
! to the experimental limits.
! 
! The results are given as (n+1)-component arrays (starting from 0), 
! where n is the total number of Higgs bosons in the model (neutral+charged).
! The zeroth component gives the combined results (equivalent to run_HiggsBounds).
!
! Arguments (output): 
!   * HBresult = 1 if point is unexcluded, 0 if excluded, -1 if parameter point is invalid
!   * chan = number of channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!   * obsratio = ratio of the theoretical rate to the observed limit for this channel
!   * ncombined = number of Higgs combined in order to calculate this obsratio
!    (see manual for more precise definitions))
 use usefulbits, only : theo,res,inputsub,just_after_run,ndmh, &
&                       np,Hneut,Hplus,dmhsteps,ndat,fullHBres,small_mh
 use channels, only : check_channels
 !use input, only : test_input
 use theo_manip, only : complete_theo, recalculate_theo_for_datapoint

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none                      
 !----------------------------------------output
 integer,intent(out)::  HBresult(0:np(Hneut)+np(Hplus))
 integer,intent(out)::  chan(0:np(Hneut)+np(Hplus))
 integer,intent(out)::  ncombined(0:np(Hneut)+np(Hplus))

 double precision,intent(out) :: obsratio(0:np(Hneut)+np(Hplus))
 
 double precision :: Mhneut(np(Hneut))
 double precision :: Mhch(np(Hplus))
 !-------------------------------------internal
 integer :: n,i,j,ind,part
 !---------------------------------------------

! print *, "Running HiggsBounds in Normal Mode (most sensitive limit considered for each Higgs boson)"
 
 if (lbound(HBresult,dim=1).NE.0) stop "run_HiggsBounds_full: Array HBresult must begin with element 0"
 if (ubound(HBresult,dim=1).NE.(np(Hneut)+np(Hplus))) then
  stop "run_HiggsBounds_full: Upper limit of array HBresult must be equal to number of Higgses"
 endif 
 if (lbound(chan,dim=1).NE.0) stop "run_HiggsBounds_full: Array chan must begin with element 0"
 if (ubound(chan,dim=1).NE.(np(Hneut)+np(Hplus))) then
  stop "run_HiggsBounds_full: Upper limit of array chan must be equal to number of Higgses"
 endif 
 if (lbound(obsratio,dim=1).NE.0) stop "run_HiggsBounds_full: Array obsratio must begin with element 0"
 if (ubound(obsratio,dim=1).NE.(np(Hneut)+np(Hplus))) then
  stop "run_HiggsBounds_full: Upper limit of array obsratio must be equal to number of Higgses"
 endif 
 if (lbound(ncombined,dim=1).NE.0) stop "run_HiggsBounds_full: Array ncombined must begin with element 0"
 if (ubound(ncombined,dim=1).NE.(np(Hneut)+np(Hplus))) then
  stop "run_HiggsBounds_full: Upper limit of array ncombined must be equal to number of Higgses"
 endif 

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 do i=1,ubound(inputsub,dim=1)
   if(  inputsub(i)%req .ne. inputsub(i)%stat  )then
       write(*,*)'subroutine '//trim(adjustl(inputsub(i)%desc))
       write(*,*)'should be called once and only once before each call to'
       write(*,*)'subroutine run_HiggsBounds.'
       stop'error in subroutine run_HiggsBounds'
   endif
   inputsub(i)%stat=0!now we have used this input, set back to zero
 enddo

 call complete_theo    
 
  do n=1,ndat
        
   theo(n)%particle(Hneut)%Mc = theo(n)%particle(Hneut)%M
   theo(n)%particle(Hplus)%Mc = theo(n)%particle(Hplus)%M

   call get_mass_variation_param(n)

    do i=0,ubound(Hbresult,dim=1)
     obsratio(i) = -999d0
     HBresult(i) = 1
     chan(i) = -999
     ncombined(i) = -999
    enddo

!  Do we have mass uncertainties to take care off
    IF (ndmh.GT.0) THEN
   
!    print *, "Running HiggsBounds with Higgs mass uncertainties"
!    write(*,*) theo(n)%particle(Hplus)%dM
      
     if(np(Hneut).ne.0) Mhneut =  theo(n)%particle(Hneut)%M 
     if(np(Hplus).ne.0) Mhch = theo(n)%particle(Hplus)%M

!   Loop over all Higgses
     do i=1,np(Hneut)+np(Hplus)
      obsratio(i) = 1.D23

      IF (i.LE.np(Hneut)) THEN
       ind = i
       part = Hneut
      ELSE
       ind = i-np(Hneut)
       part = Hplus
      ENDIF
          
!     Check for mass steps for this particular Higgs boson
      IF(theo(n)%particle(part)%dMh(ind).GT.small_mh) THEN

!       theo(n)%particle(part)%M(ind)=theo(n)%particle(part)%M(ind)  &
!       &                             -(dmhsteps-1)/2*theo(n)%particle(part)%dMh(ind)

       theo(n)%particle(part)%M(ind)=theo(n)%particle(part)%M(ind)  &
       &                             -theo(n)%particle(part)%dMh(ind)

       do j=1,dmhsteps 
        
!        print *, theo(n)%particle(Hneut)%M, theo(n)%particle(Hplus)%M

        call recalculate_theo_for_datapoint(n)       

        call check_channels(theo(n),res(n),i)
     	
        IF (res(n)%obsratio(1).LT.obsratio(i)) THEN
         HBresult(i)    = res(n)%allowed95(1)       
         chan(i)        = res(n)%chan(1)       
         obsratio(i)    = res(n)%obsratio(1)
         ncombined(i)   = res(n)%ncombined(1)        
        ENDIF
        
!        print *, i,theo(n)%particle(part)%M(ind),HBresult(i),chan(i),obsratio(i),ncombined(i)
        
        theo(n)%particle(part)%M(ind)= theo(n)%particle(part)%M(ind)  &
        &                             +theo(n)%particle(part)%dMh(ind)/(dmhsteps-1)*2
        
      enddo
     else
      call recalculate_theo_for_datapoint(n)
      call check_channels(theo(n),res(n),i)

      HBresult(i)    = res(n)%allowed95(1)       
      chan(i)        = res(n)%chan(1)       
      obsratio(i)    = res(n)%obsratio(1)
      ncombined(i)   = res(n)%ncombined(1)            
     endif 
     
!	 Logical OR between exclusions (one Higgs excluded = combined exclusion)
     HBresult(0) = HBresult(0) * HBresult(i)
    
!    Save the data for the Higgs that has the highest ratio of theory/obs 
     IF (obsratio(i).GT.obsratio(0)) THEN
      chan(0)        = chan(i)       
      obsratio(0)    = obsratio(i)
      ncombined(0)   = ncombined(i)        
     ENDIF
     
     theo(n)%particle(Hneut)%M = Mhneut
     theo(n)%particle(Hplus)%M = Mhch
      
    enddo
 !    return 
     
   ELSE
 
 
 !!    print *, "Running HiggsBounds without Higgs mass uncertainties"

    call recalculate_theo_for_datapoint(n)       

   do i=1,np(Hneut)+np(Hplus)
    call check_channels(theo(n),res(n),i)
    HBresult(i)    = res(n)%allowed95(1)       
    chan(i)        = res(n)%chan(1)       
    obsratio(i)    = res(n)%obsratio(1)
    ncombined(i)   = res(n)%ncombined(1)        
    
    HBresult(0) = HBresult(0) * res(n)%allowed95(1)
    
    IF (obsratio(i).GT.obsratio(0)) THEN
!     write(*,*) "hello: ", n, i
     chan(0)        = res(n)%chan(1)       
     obsratio(0)    = res(n)%obsratio(1)
     ncombined(0)   = res(n)%ncombined(1)        
    ENDIF
       
!    IF (i.LE.np(Hneut)) THEN
!     print *, i,theo(n)%particle(Hneut)%M(i),HBresult(i),chan(i),obsratio(i),ncombined(i),HBresult(0), obsratio(0)   
!    ELSE
!     print *, i,theo(n)%particle(Hplus)%M(i-np(Hneut)),HBresult(i),chan(i),obsratio(i),ncombined(i),HBresult(0), obsratio(0)   
!    endif
   enddo
  ENDIF

  fullHBres(n)%allowed95=HBresult(0)
  fullHBres(n)%chan=chan(0)
  fullHBres(n)%obsratio=obsratio(0)
  fullHBres(n)%ncombined=ncombined(0)
 
 enddo


 just_after_run=.True.   
 

end subroutine run_HiggsBounds_full
!************************************************************



subroutine run_HiggsBounds_classic( HBresult,chan,obsratio,ncombined)
! This subroutine can be called by the user after HiggsBounds_initialize has been called.
! The input routines, where required, should be called before calling run_HiggsBounds.
! It takes theoretical predictions for a particular parameter point 
! in the model and calls subroutines which compare these predictions 
! to the experimental limits
! Arguments (output): 
!   * HBresult = 1 if point is unexcluded, 0 if excluded, -1 if parameter point is invalid
!   * chan = number of channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!   * obsratio = ratio of the theoretical rate to the observed limit for this channel
!   * ncombined = number of Higgs combined in order to calculate this obsratio
!    (see manual for more precise definitions))
 use usefulbits, only : theo,res,debug,inputsub,just_after_run,ndmh,diffmhneut,diffmhch, &
 np,Hneut,Hplus,full_dmth_variation,dmhsteps, ndat,fullHBres
 use channels, only : check_channels
 !use input, only : test_input
 use theo_manip, only : complete_theo, recalculate_theo_for_datapoint

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none                      
 !----------------------------------------output
 integer,intent(out)::     HBresult,chan,ncombined
 double precision,intent(out) :: obsratio
 
 double precision :: Mhneut(np(Hneut))
 double precision :: Mhch(np(Hplus))
 !-------------------------------------internal
 integer :: n,i
 integer :: HBresult_tmp,chan_tmp,ncombined_tmp
 double precision :: obsratio_tmp
 !---------------------------------------------
! n=1
 
!  print *, "Running HiggsBounds in Classic Mode (globally most sensitive limit only)"
    
  
 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 do i=1,ubound(inputsub,dim=1)
   if(  inputsub(i)%req .ne. inputsub(i)%stat  )then
       write(*,*)'subroutine '//trim(adjustl(inputsub(i)%desc))
       write(*,*)'should be called once and only once before each call to'
       write(*,*)'subroutine run_HiggsBounds.'
       stop'error in subroutine run_HiggsBounds'
   endif
   inputsub(i)%stat=0!now we have used this input, set back to zero
 enddo

 call complete_theo    

 do n=1,ndat
        
  theo(n)%particle(Hneut)%Mc = theo(n)%particle(Hneut)%M

  call get_mass_variation_param(n)

  IF (ndmh.GT.0) THEN
   
    if(np(Hneut).ne.0) Mhneut =  theo(n)%particle(Hneut)%M     
    if(np(Hplus).ne.0) Mhch = theo(n)%particle(Hplus)%M
  
    obsratio_tmp = 10.0E6 ! Set to very large initial value
    do i=1,dmhsteps**ndmh 
     theo(n)%particle(Hneut)%M = diffMhneut(i,:)
     theo(n)%particle(Hplus)%M = diffMhch(i,:)
		 
     if(debug)write(*,*)'manipulating input...'                 ; call flush(6)
     call recalculate_theo_for_datapoint(n)
	
     if(debug)write(*,*)'compare each data point to the experimental bounds...' ; call flush(6)                  
     call check_channels(theo(n),res(n),0)
		
     HBresult    = res(n)%allowed95(1)
     chan        = res(n)%chan(1)       
     obsratio    = res(n)%obsratio(1)
     ncombined   = res(n)%ncombined(1)        
!     print *, HBresult, chan, obsratio, ncombined
		
     IF (.NOT.full_dmth_variation) THEN
      IF (HBresult.EQ.1) THEN
!		 	 theo(n)%particle(Hneut)%M = Mhneut
! 			 theo(n)%particle(Hplus)%M = Mhch
       just_after_run=.True.   
       exit
      ENDIF
     ELSE
      IF (obsratio.lt.obsratio_tmp) THEN
       HBresult_tmp  = HBresult
       chan_tmp      = chan
       obsratio_tmp  = obsratio
       ncombined_tmp = ncombined
      ENDIF
     ENDIF  
   enddo
   
   IF (full_dmth_variation) THEN
    HBresult  = HBresult_tmp
    chan      = chan_tmp
    obsratio  = obsratio_tmp
    ncombined = ncombined
!    theo(n)%particle(Hneut)%M = Mhneut
!    theo(n)%particle(Hplus)%M = Mhch
    just_after_run=.True.    
!	return
   ENDIF
   
   theo(n)%particle(Hneut)%M = Mhneut
   theo(n)%particle(Hplus)%M = Mhch
   call recalculate_theo_for_datapoint(n)
   call check_channels(theo(n),res(n),0)

 ELSE
 
   if(debug)write(*,*)'manipulating input...'                 ; call flush(6)
   call recalculate_theo_for_datapoint(n)

   if(debug)write(*,*)'compare each data point to the experimental bounds...' ; call flush(6)                  
   
   call check_channels(theo(n),res(n),0)
   
   HBresult    = res(n)%allowed95(1)
   chan        = res(n)%chan(1)       
   obsratio    = res(n)%obsratio(1)
   ncombined   = res(n)%ncombined(1)        

   just_after_run=.True.   
 ENDIF

 fullHBres(n)%allowed95=HBresult
 fullHBres(n)%chan=chan
 fullHBres(n)%obsratio=obsratio
 fullHBres(n)%ncombined=ncombined

 enddo 

 just_after_run=.True.    


end subroutine run_HiggsBounds_classic
!************************************************************


!************************************************************ 
subroutine HiggsBounds_SLHA_output
!**** ******************************************************** 
use usefulbits, only : whichinput,just_after_run
use output, only : do_output

  if(.not.just_after_run)then
   stop'subroutine run_HiggsBounds should be called before subroutine HiggsBounds_SLHA_output' 
  endif  

  select case(whichinput)
  case('SLHA') 
    call do_output
  case default
    stop'The subroutine HiggsBounds_SLHA_output should only be used when whichinput=SLHA'
  end select

end subroutine HiggsBounds_SLHA_output

#ifdef enableCHISQ
!************************************************************      
subroutine initialize_HiggsBounds_chisqtables
! use S95tables, only : S95_t2
 use S95tables_type3
 use usefulbits, only : allocate_if_stats_required,theo
 implicit none            
 
 if(allocated(theo))then
  stop 'subroutine initialize_HiggsBounds_chisqtables should be called before subroutine HiggsBounds_initialize'
 elseif(allocated(clsb_t3))then
  stop 'subroutine initialize_HiggsBounds_chisqtables has already been called once'
 endif 

 allocate(clsb_t3(ntable3))

 call initializetables_type3_blank(clsb_t3)
 call initializetables3(clsb_t3)

 call readclsbfiles_binary

 if(allocated(allocate_if_stats_required))then
   stop'error in subroutine initialize_HiggsBounds_chisqtables'
 else
   allocate(allocate_if_stats_required(1))
 endif


end subroutine initialize_HiggsBounds_chisqtables

!************************************************************
subroutine finish_HiggsBounds_chisqtables
!************************************************************
 use S95tables_type3
 use usefulbits, only : allocate_if_stats_required
 implicit none      
 integer :: x      

 if(.not.allocated(clsb_t3))then
  stop 'initialize_HiggsBounds_chisqtables should be called first'
 endif

 do x=lbound(clsb_t3,dim=1),ubound(clsb_t3,dim=1)
  deallocate(clsb_t3(x)%dat)
 enddo
 deallocate(filename)
 deallocate(clsb_t3)
 deallocate(allocate_if_stats_required)

end subroutine finish_HiggsBounds_chisqtables

!************************************************************
subroutine HB_calc_stats(theory_uncertainty_1s,chisq_withouttheory,chisq_withtheory,chan2)
!************************************************************
! this is in the middle of development! DO NOT USE!
 use usefulbits, only : res,theo,pr,just_after_run,vsmall
 use interpolate
 use S95tables_type1
 use S95tables_type3
 use S95tables
 use extra_bits_for_chisquared
 implicit none      

  integer,intent(out)::chan2
  integer :: x,c,z,y
  integer :: id
  double precision, intent(in) :: theory_uncertainty_1s
  double precision :: chisq_withouttheory,chisq_withtheory
  double precision :: low_chisq,sigma

  x=1
  low_chisq=1.0D-2

  if(.not.allocated(theo))then
   stop 'subroutine HiggsBounds_initialize must be called first'
  elseif(.not.allocated(clsb_t3))then
   stop 'subroutine initialize_HiggsBounds_chisqtables must be called first'
  elseif(.not.just_after_run)then
   stop 'subroutine run_HiggsBounds must be called first'
  endif

  sigma=theory_uncertainty_1s
  if(sigma.lt.vsmall)then 
   write(*,*)'Warning: will not calculate chi^2 with theory uncertainty'
  endif

  chisq_withtheory              = -2.0D0
  chisq_withouttheory           = -2.0D0

  z=2;
  c= res(x)%chan(z)
  chan2=c

  if(res(x)%allowed95(z).eq.-1)then! labels an unphysical parameter point
   chisq_withtheory             =-1.0D0
   chisq_withouttheory          =-1.0D0
  elseif( c.gt.0 )then ! labels a physical parameter point and a real channel                      
  
   id=S95_t1_or_S95_t2_idfromelementnumber(pr(c)%ttype,pr(c)%tlist)
   y=clsb_t3elementnumber_from_S95table(pr(c)%ttype,id)

   if(y.gt.0)then

    !------------------------------

    call get_chisq(sigma,res(x)%axis_i(z),res(x)%axis_j(z),res(x)%sfactor(z), &
      &  y,chisq_withouttheory,chisq_withtheory)

    !-------------------------------

   else
     write(*,*)'hello y=',y
     stop'problem here with y'
   endif

  else
   chisq_withtheory             =0.0D0
   chisq_withouttheory          =0.0D0
  endif  

end subroutine HB_calc_stats
#endif
!************************************************************
subroutine finish_HiggsBounds
! This subroutine needs to be called right at the end, to close files
! and deallocate arrays
!************************************************************
 use usefulbits, only : deallocate_usefulbits,debug,theo,debug,inputsub, &
   & file_id_debug1,file_id_debug2
 use S95tables, only : deallocate_S95tables
 use theory_BRfunctions, only : deallocate_BRSM

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif
      
 if(debug)then
  close(file_id_debug2)
  close(file_id_debug1)
 endif

 if(.not.allocated(theo))then
  stop 'HiggsBounds_initialize  should be called first'
 endif

 if(debug)write(*,*)'finishing off...'                      ; call flush(6)
 call deallocate_BRSM
 call deallocate_S95tables
 call deallocate_usefulbits

 if(debug)write(*,*)'finished'                              ; call flush(6)
 
 if(allocated(inputsub)) deallocate(inputsub)
end subroutine finish_HiggsBounds
!
!!************************************************************
!	
!subroutine run_HiggsBounds_effC(nHdummy,Mh,GammaTotal,             &   
!     &          g2hjbb,g2hjtautau,g2hjWW,g2hjZZ,                   &
!     &          g2hjgaga,g2hjgg,g2hjhiZ_nHbynH,                    &
!     &          BR_hjhihi_nHbynH,                                  &
!     &          HBresult,chan,                                     &
!     &          obsratio, ncombined                                )
!! Obsolete subroutine
!!************************************************************
!
! implicit none
!
! !----------------------------------------input
! integer,intent(in) :: nHdummy
! double precision,intent(in) :: Mh(nHdummy),GammaTotal(nHdummy),&   
!     &          g2hjbb(nHdummy),g2hjtautau(nHdummy),            &
!     &          g2hjWW(nHdummy),g2hjZZ(nHdummy),                &
!     &          g2hjgaga(nHdummy),g2hjgg(nHdummy),              &
!     &          g2hjhiZ_nHbynH(nHdummy,nHdummy),                &
!     &          BR_hjhihi_nHbynH(nHdummy,nHdummy) 
! !----------------------------------------output
! integer ::     HBresult,chan,ncombined
! double precision :: obsratio
! !----------------------------------------------
!	  
! call attempting_to_use_an_old_HB_version('effC')
!
!end subroutine run_HiggsBounds_effC
!!************************************************************	
!subroutine run_HiggsBounds_part(nHdummy,Mh,                  &
!     &          CS_lep_hjZ_ratio,                            &
!     &          CS_lep_hjhi_ratio_nHbynH,                    &
!     &          CS_tev_gg_hj_ratio,CS_tev_bb_hj_ratio,       &
!     &          CS_tev_bg_hjb_ratio,                         &
!     &          CS_tev_ud_hjWp_ratio,CS_tev_cs_hjWp_ratio,   & 
!     &          CS_tev_ud_hjWm_ratio,CS_tev_cs_hjWm_ratio,   & 
!     &          CS_tev_dd_hjZ_ratio,CS_tev_uu_hjZ_ratio,     &
!     &          CS_tev_ss_hjZ_ratio,CS_tev_cc_hjZ_ratio,     &
!     &          CS_tev_bb_hjZ_ratio,                         &
!     &          CS_tev_pp_vbf_ratio,                         &
!     &          BR_hjbb,BR_hjtautau,                         &
!     &          BR_hjWW,BR_hjgaga,                           & 
!     &          BR_hjhihi_nHbynH,                            &
!     &          HBresult,chan,                               &
!     &          obsratio, ncombined                          )
!! Obsolete subroutine
!!************************************************************
!
! implicit none
!
! !----------------------------------------input
! integer , intent(in) :: nHdummy
! double precision,intent(in) ::Mh(nHdummy),                                 &
!     &          CS_lep_hjZ_ratio(nHdummy),                                  &
!     &          CS_lep_hjhi_ratio_nHbynH(nHdummy,nHdummy),                  &
!     &          CS_tev_gg_hj_ratio(nHdummy),CS_tev_bb_hj_ratio(nHdummy),    &
!     &          CS_tev_bg_hjb_ratio(nHdummy),                               &
!     &          CS_tev_ud_hjWp_ratio(nHdummy),CS_tev_cs_hjWp_ratio(nHdummy),& 
!     &          CS_tev_ud_hjWm_ratio(nHdummy),CS_tev_cs_hjWm_ratio(nHdummy),& 
!     &          CS_tev_dd_hjZ_ratio(nHdummy),CS_tev_uu_hjZ_ratio(nHdummy),  &
!     &          CS_tev_ss_hjZ_ratio(nHdummy),CS_tev_cc_hjZ_ratio(nHdummy),  &
!     &          CS_tev_bb_hjZ_ratio(nHdummy),                               &
!     &          CS_tev_pp_vbf_ratio(nHdummy),                               &
!     &          BR_hjbb(nHdummy),BR_hjtautau(nHdummy),                      &
!     &          BR_hjWW(nHdummy),BR_hjgaga(nHdummy),                        &
!     &          BR_hjhihi_nHbynH(nHdummy,nHdummy)
! !---------------------------------------output
! integer ::     HBresult,chan,ncombined
! double precision :: obsratio
! !-----------------------------------------------
!	  
! call attempting_to_use_an_old_HB_version('part')
!
!end subroutine run_HiggsBounds_part
!!************************************************************	
!subroutine run_HiggsBounds_hadr(nHdummy,Mh,                  &
!     &          CS_lep_hjZ_ratio,CS_lep_hjhi_ratio_nHbynH,   &
!     &          CS_tev_pp_hj_ratio, CS_tev_pp_hjb_ratio,     &
!     &          CS_tev_pp_hjW_ratio,CS_tev_pp_hjZ_ratio,     &
!     &          CS_tev_pp_vbf_ratio,                         &
!     &          BR_hjbb,BR_hjtautau,                         &
!     &          BR_hjWW,BR_hjgaga,                           &
!     &          BR_hjhihi_nHbynH,                            &
!     &          HBresult,chan,                               &
!     &          obsratio, ncombined                          )
!! Obsolete subroutine
!!************************************************************
!
! implicit none
! !----------------------------------------input
! integer,intent(in) :: nHdummy
! double precision,intent(in) :: Mh(nHdummy),              &
!     &          CS_lep_hjZ_ratio(nHdummy),                &
!     &          CS_lep_hjhi_ratio_nHbynH(nHdummy,nHdummy),&
!     &          CS_tev_pp_hj_ratio(nHdummy),              &
!     &          CS_tev_pp_hjb_ratio(nHdummy),             &
!     &          CS_tev_pp_hjW_ratio(nHdummy),             &
!     &          CS_tev_pp_hjZ_ratio(nHdummy),             &
!     &          CS_tev_pp_vbf_ratio(nHdummy),             &
!     &          BR_hjbb(nHdummy),BR_hjtautau(nHdummy),    &
!     &          BR_hjWW(nHdummy),BR_hjgaga(nHdummy),      &
!     &          BR_hjhihi_nHbynH(nHdummy,nHdummy)
! !---------------------------------------output
! integer ::     HBresult,chan,ncombined
! double precision :: obsratio
! !---------------------------------------------
!
! call attempting_to_use_an_old_HB_version('hadr')
!
!
!end subroutine run_HiggsBounds_hadr
!!************************************************************	
