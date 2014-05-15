! This file is part of HiggsBounds
!  -KW
!******************************************************************
module S95tables_type1
!******************************************************************      

 implicit none
 
 !table type 1-----------------------------      
 type table1      
  integer :: id,nx,particle_x !see usefulbits.f90 for key to particle codes n.b. they're NOT pdg
  character(LEN=45) :: label
  character(LEN=100) :: desc
  character(LEN=3) :: expt 
  double precision :: lumi, energy
  double precision :: xmax,xmin,sep,deltax
  integer :: SMlike
  double precision, allocatable :: dat(:,:) !in dat(a,b), a=row, b=1,2 for obs,pred 
 end type
 !------------------------------------------ 

 integer,parameter :: file_id_1=10 !same as file_id_common in usefulbits.f90

 contains

 !************************************************************ 
 subroutine initializetables_type1_blank(tablet1)
 !***********************************************************  
 ! still leaves dat unallocated
  integer:: i
  type(table1) :: tablet1(:)

  do i=lbound(tablet1,dim=1),ubound(tablet1,dim=1)
   tablet1(i)%id         = -1
   tablet1(i)%nx         = -1
   tablet1(i)%particle_x = -1
   tablet1(i)%label      = ''
   tablet1(i)%desc       = ''
   tablet1(i)%expt       = ''
   tablet1(i)%lumi       = -1.0D0
   tablet1(i)%energy     = -1.0D0
   tablet1(i)%xmax       = -1.0D0
   tablet1(i)%xmin       = -1.0D0
   tablet1(i)%sep        = -1.0D0
   tablet1(i)%deltax     = -1.0D0
   tablet1(i)%SMlike     = 0
   enddo
 
 end subroutine initializetables_type1_blank

 !************************************************************ 
 subroutine copy_type1(tablet1_orig,tablet1_copy)
 !***********************************************************  
 ! note tablet1_1,tablet1_2 are not arrays 
 ! still leaves dat uncopied
  type(table1) :: tablet1_orig
  type(table1) :: tablet1_copy

  tablet1_copy%id         = tablet1_orig%id
  tablet1_copy%nx         = tablet1_orig%nx
  tablet1_copy%particle_x = tablet1_orig%particle_x
  tablet1_copy%label      = tablet1_orig%label
  tablet1_copy%expt       = tablet1_orig%expt
  tablet1_copy%xmax       = tablet1_orig%xmax
  tablet1_copy%xmin       = tablet1_orig%xmin
  tablet1_copy%sep        = tablet1_orig%sep
  tablet1_copy%deltax     = tablet1_orig%deltax
 
 end subroutine copy_type1
 !***********************************************************
 function t1elementnumberfromid(t1,id)
  !--------------------------------------input
  type(table1), intent(in) :: t1(:)
  integer, intent(in) :: id
  !-----------------------------------function
  integer :: t1elementnumberfromid
  !-----------------------------------internal
  integer :: n,x
  !-------------------------------------------

  n=0
  do x=lbound(t1,dim=1),ubound(t1,dim=1)
   if(t1(x)%id.eq.id)then
    n=n+1
    t1elementnumberfromid=x
   endif
  enddo

  if(n.ne.1)stop'problem in function t3elementnumberfromid 1'

 end function t1elementnumberfromid

 !************************************************************ 
 subroutine initializetables1(S95_t1)
 !***********************************************************  
 ! fills S95_t1
 !*********************************************************** 
  use store_pathname
  use usefulbits, only: Hneut,Hplus,file_id_common2
  implicit none

  !--------------------------------------input
  type(table1) :: S95_t1(:)
  !-----------------------------------internal
  logical :: newtables
  integer :: x,xbeg,xend      
  character(len=100),allocatable :: filename(:)
  character(LEN=pathname_length+150) :: fullfilename    
  integer :: col
  integer :: ios
  !-------------------------------------------  
  
  xbeg=lbound(S95_t1,dim=1)
  xend=ubound(S95_t1,dim=1)
  
  allocate(filename(xbeg:xend))    
  x=xbeg-1

  !instead, could read in the values of xmin,xmax,sep from the
  !files, but it's kinda nice having them all here to refer to

  newtables=.True.  ! i.e. use the recommended LEP single Higgs tables
  if(newtables)then   
   x=x+1
   S95_t1(x)%id=142
   S95_t1(x)%particle_x=Hneut
   S95_t1(x)%expt='LEP'                        
   S95_t1(x)%label='hep-ex/0602042, table 14b (LEP)' ! table 14b
   S95_t1(x)%energy=0.208D0
   S95_t1(x)%xmin=12.0D0
   S95_t1(x)%xmax=120.0D0
   S95_t1(x)%sep=0.5D0
   filename(x)='lep210_hbb'
  
   x=x+1
   S95_t1(x)%id=143
   S95_t1(x)%particle_x=Hneut
   S95_t1(x)%expt='LEP'            
   S95_t1(x)%label='hep-ex/0602042, table 14c (LEP)' ! table 14c
   S95_t1(x)%energy=0.208D0   
   S95_t1(x)%xmin=4.0D0
   S95_t1(x)%xmax=120.0D0
   S95_t1(x)%sep=0.5D0      
   filename(x)='lep210_htt_interpol'
  else
   write(*,*)'WARNING: using old LEP tables' 
   x=x+1
   S95_t1(x)%id=142
   S95_t1(x)%particle_x=Hneut
   S95_t1(x)%expt='LEP'            
   S95_t1(x)%label='LEP table 14b' ! table 14b
   S95_t1(x)%energy=0.208D0
   S95_t1(x)%xmin=1.0D0
   S95_t1(x)%xmax=140.0D0
   S95_t1(x)%sep=0.1D0
   filename(x)='old-s95_h2z_bbz'
  
   x=x+1
   S95_t1(x)%id=143 
   S95_t1(x)%particle_x=Hneut
   S95_t1(x)%expt='LEP'    
   S95_t1(x)%label='LEP table 14c' ! table 14c
   S95_t1(x)%energy=0.208D0   
   S95_t1(x)%xmin=1.0D0
   S95_t1(x)%xmax=140.0D0
   S95_t1(x)%sep=0.1D0
   filename(x)='old-s95_h2z_ttz'  
  endif

  x=x+1
  S95_t1(x)%id=300 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0206022 (OPAL)' 
  S95_t1(x)%energy=0.208D0  
  S95_t1(x)%xmin=1.0D0
  S95_t1(x)%xmax=100.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='lep_decaymodeindep'  

  x=x+1
  S95_t1(x)%id=400 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0107032v1 (LEP)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=118.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='LEP_h-invisible'  

  x=x+1
  S95_t1(x)%id=500 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='LHWG Note 2002-02' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=20.0D0
  S95_t1(x)%xmax=116.0D0
  S95_t1(x)%sep=2.0D0
  filename(x)='LEP_h-gammagamma'  

  x=x+1
  S95_t1(x)%id=600 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='LHWG (unpublished)'!uses hep-ex/0510022,hep-ex/0205055,hep-ex/0312042,hep-ex/0408097 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=20.0D0
  S95_t1(x)%xmax=128.6D0
  S95_t1(x)%sep=0.1D0
  filename(x)='LEP_h-2jets' 

  x=x+1
  S95_t1(x)%id=711 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0410017 (DELPHI)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=12.0D0
  S95_t1(x)%xmax=50.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='Delphi_yuk_h_bbbb'  

  x=x+1
  S95_t1(x)%id=713 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0410017 (DELPHI)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=12.0D0
  S95_t1(x)%xmax=50.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='Delphi_yuk_a_bbbb'  

  x=x+1
  S95_t1(x)%id=721 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0410017 (DELPHI)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=4.0D0
  S95_t1(x)%xmax=50.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='Delphi_yuk_h_bbtautau'  

  x=x+1
  S95_t1(x)%id=741 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0111010 (OPAL)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=4.0D0
  S95_t1(x)%xmax=12.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='OPAL_yuk_h_bbtautau' 

  x=x+1
  S95_t1(x)%id=723 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0410017 (DELPHI)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=4.0D0
  S95_t1(x)%xmax=50.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='Delphi_yuk_a_bbtautau' 

  x=x+1
  S95_t1(x)%id=743 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0111010 (OPAL)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=4.0D0
  S95_t1(x)%xmax=12.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='OPAL_yuk_a_bbtautau' 

  x=x+1
  S95_t1(x)%id=731 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0410017 (DELPHI)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=4.0D0
  S95_t1(x)%xmax=27.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='Delphi_yuk_h_tautautautau' 

  x=x+1
  S95_t1(x)%id=733 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0410017 (DELPHI)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=4.0D0
  S95_t1(x)%xmax=26.0D0
  S95_t1(x)%sep=1.0D0
  filename(x)='Delphi_yuk_a_tautautautau' 

  x=x+1
  S95_t1(x)%id=402 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0401022 (Delphi)' 
  S95_t1(x)%energy=0.208D0  
  S95_t1(x)%xmin=40.0D0
  S95_t1(x)%xmax=114.0D0
  S95_t1(x)%sep=2.0D0
  filename(x)='Delphi_h-invisible'

  x=x+1
  S95_t1(x)%id=403 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='hep-ex/0501033 (L3)' 
  S95_t1(x)%energy=0.208D0
  S95_t1(x)%xmin=50.0D0
  S95_t1(x)%xmax=110.0D0
  S95_t1(x)%sep=5.0D0
  filename(x)='L3_h-invisible' 

  x=x+1
  S95_t1(x)%id=401 
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='LEP'    
  S95_t1(x)%label='[hep-ex] arXiv:0707.0373 (OPAL)' 
  S95_t1(x)%energy=0.208D0  
  S95_t1(x)%xmin=5.0D0
  S95_t1(x)%xmax=115.0D0
  S95_t1(x)%sep=5.0D0
  filename(x)='OPAL_h-invisible'

  !x=x+1 
  !S95_t1(x)%id=803
  !S95_t1(x)%particle_x=Hplus
  !S95_t1(x)%expt='LEP'
  !S95_t1(x)%label='[hep-ex] arxiv:0812.0267 (OPAL)'
  !S95_t1(x)%xmin=50.0D0
  !S95_t1(x)%xmax=93.0D0
  !S95_t1(x)%sep=1.0D0 
  !filename(x)='OPAL_HpHm_taunutaunu'

  !x=x+1 
  !S95_t1(x)%id=802
  !S95_t1(x)%particle_x=Hplus
  !S95_t1(x)%expt='LEP'
  !S95_t1(x)%label='[hep-ex] arxiv:0812.0267 (OPAL)'
  !S95_t1(x)%xmin=50.0D0
  !S95_t1(x)%xmax=93.0D0
  !S95_t1(x)%sep=1.0D0 
  !filename(x)='OPAL_HpHm_qqtaunu'

  !x=x+1 
  !S95_t1(x)%id=801
  !S95_t1(x)%particle_x=Hplus
  !S95_t1(x)%expt='LEP'
  !S95_t1(x)%label='[hep-ex] arxiv:0812.0267 (OPAL)'
  !S95_t1(x)%xmin=50.0D0
  !S95_t1(x)%xmax=93.0D0
  !S95_t1(x)%sep=1.0D0 
  !filename(x)='OPAL_HpHm_qqqq'

  x=x+1 
  S95_t1(x)%id=821
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='LEP'
  S95_t1(x)%label='hep-ex/0107031 (LHWG)'
  S95_t1(x)%energy=0.208D0  
  S95_t1(x)%xmin=60.0D0
  S95_t1(x)%xmax=90.0D0
  S95_t1(x)%sep=1.0D0 
  filename(x)='LEP_HpHm_qqqq'

  x=x+1 
  S95_t1(x)%id=811
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='LEP'
  S95_t1(x)%label='hep-ex/0404012 (Delphi)'
  S95_t1(x)%energy=0.208D0  
  S95_t1(x)%xmin=52.0D0
  S95_t1(x)%xmax=94.0D0
  S95_t1(x)%sep=2.0D0 
  filename(x)='Delphi_HpHm_qqqq'

  x=x+1 
  S95_t1(x)%id=813
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='LEP'
  S95_t1(x)%label='hep-ex/0404012 (Delphi)'
  S95_t1(x)%energy=0.208D0  
  S95_t1(x)%xmin=52.0D0
  S95_t1(x)%xmax=94.0D0
  S95_t1(x)%sep=2.0D0 
  filename(x)='Delphi_HpHm_taunutaunu'

!----------------------- Z H -> l l b b -------------------------

 ! x=x+1 
!  S95_t1(x)%id=10235 
!  S95_t1(x)%particle_x=Hneut 
!  S95_t1(x)%expt='CDF' 
!  S95_t1(x)%label='CDF Note 10235'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=5.7D0
!  S95_t1(x)%xmin=100.0D0  
!  S95_t1(x)%xmax=150.0D0 
!  S95_t1(x)%sep=5.0D0  
!  filename(x)='CDF_ZH_llbb_5.7fb_10235' 
!
!  x=x+1 
!  S95_t1(x)%id=3047
!  S95_t1(x)%particle_x=Hneut 
!  S95_t1(x)%expt='CDF' 
!  S95_t1(x)%label='[hep-ex] arXiv:1009.3047 (CDF)' 
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=4.1D0
!  S95_t1(x)%xmin=100.0D0  
!  S95_t1(x)%xmax=150.0D0 
!  S95_t1(x)%sep=5.0D0  
!  filename(x)='CDF_ZH_llbb_4.1fb_3047' 

  x=x+1 
  S95_t1(x)%id=10799 
  S95_t1(x)%particle_x=Hneut 
  S95_t1(x)%expt='CDF' 
  S95_t1(x)%label='CDF Note 10799'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.45D0
  S95_t1(x)%xmin=90.0D0  
  S95_t1(x)%xmax=150.0D0 
  S95_t1(x)%sep=5.0D0  
  filename(x)='CDF_ZH_llbb_9.45fb_10799' 

!  x=x+1 
!  S95_t1(x)%id=6166
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6166'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.6D0
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  filename(x)='D0_ZH_llbb_8.6fb_6166' 

  x=x+1 
  S95_t1(x)%id=6296
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6296'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.7D0
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='D0_ZH_llbb_9.7fb_6296' 

  !x=x+1 
  !S95_t1(x)%id=6089
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 6089'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='D0_ZH_llbb_6.2fb_6089' 

  x=x+1 
  S95_t1(x)%id=3564
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:1008.3564 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=4.2D0
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='D0_ZH_llbb_4.2fb_3564' 

  !x=x+1 
  !S95_t1(x)%id=10212
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 10212'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_VH_Metbb_5.7fb_10212'

!  x=x+1 
!  S95_t1(x)%id=6087
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6087'
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_VH_bb_6.4fb_6087'

!  x=x+1
!  S95_t1(x)%id=2012015
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2012-015'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.7D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=130.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='2012015_Atlas_VH_bb_ll_lnu_nunu_4.7fb-1'

!----------------------- V H -> b b Etmiss -------------------------

!  x=x+1 
!  S95_t1(x)%id=10583
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CDF'
!  S95_t1(x)%label='CDF Note 10583'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=7.8D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_VH_Metbb_7.8fb_10583'

!  x=x+1 
!  S95_t1(x)%id=6223
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6223'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.4D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_VH_bb_8.4fb_6223'

!  x=x+1 
!  S95_t1(x)%id=3935
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CDF'
!  S95_t1(x)%label='[hep-ex] arXiv:0911.3935v4 (CDF)'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=2.1D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_VH_Metbb_2.1fb_3935_interpol'
!
!  x=x+1 
!  S95_t1(x)%id=5285
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='[hep-ex] arXiv:0912.5285 (D0)'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=5.2D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_VH_bb_5.2fb_5285'

  !x=x+1 
  !S95_t1(x)%id=6092
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt=' D0' 
  !S95_t1(x)%label='D0 Note 6092'!this is what the note says, but the website says 6082 
  !S95_t1(x)%xmin=100.0D0  
  !S95_t1(x)%xmax=150.0D0 
  !S95_t1(x)%sep=5.0D0  
  !filename(x)='D0_WH_lnubb_5.3fb_6092' 

!  x=x+1 
!  S95_t1(x)%id=10596
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CDF'
!  S95_t1(x)%label='CDF Note 10596'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=7.5D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_VH_lnubb_7.5fb_10596'

!  x=x+1 
!  S95_t1(x)%id=2011103
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-103'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=1.04D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=130.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2011103_Atlas_VH_Vbb_1.04fb-1'

  x=x+1 
  S95_t1(x)%id=10798
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10798'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.45D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_VH_Metbb_9.45fb_10798'


  x=x+1 
  S95_t1(x)%id=6299
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6299'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.5D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_VH_bb_9.5fb_6299'

  
  x=x+1 
  S95_t1(x)%id=2012161
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-161'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=17.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=130.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2012161_Atlas_VH-Vbb_17.7fb-1'

!  x=x+1
!  S95_t1(x)%id=11031
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-031'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.7D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=135.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11031_CMS_VH-bb_BDT_4.7fb-1'

!  x=x+1
!  S95_t1(x)%id=12044
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-12-044'
!  S95_t1(x)%energy=8.0D0
!  S95_t1(x)%lumi=17.1D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=135.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='12044_CMS_VH_Vbb_17.1fb-1'
!
  x=x+1
  S95_t1(x)%id=13012
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-012'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=24.D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=135.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='13012_CMS_VH_bb_24fb-1'

!----------------------- VBF(H), H -> b b -------------------------


  x=x+1
  S95_t1(x)%id=13011
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-011'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=19.D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=115.0D0
  S95_t1(x)%xmax=135.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='13011_CMS_VBF_bb_19fb-1'



!----------------------- W H -> b b -------------------------

!  x=x+1 
!  S95_t1(x)%id=6220
!  S95_t1(x)%particle_x=Hneut 
!  S95_t1(x)%expt=' D0' 
!  S95_t1(x)%label='D0 Note 6220'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.5D0
!  S95_t1(x)%xmin=100.0D0  
!  S95_t1(x)%xmax=150.0D0 
!  S95_t1(x)%sep=5.0D0  
!  filename(x)='D0_WH_lnubb_8.5fb_6220' 

  x=x+1 
  S95_t1(x)%id=6309
  S95_t1(x)%particle_x=Hneut 
  S95_t1(x)%expt=' D0' 
  S95_t1(x)%label='D0 Note 6309'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.6D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=100.0D0  
  S95_t1(x)%xmax=150.0D0 
  S95_t1(x)%sep=5.0D0
  S95_t1(x)%deltax=0.0D0       
  filename(x)='D0_WH_lnubb_9.7_6309' 

!  x=x+1 
!  S95_t1(x)%id=10239
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CDF'
!  S95_t1(x)%label='CDF Note 10239' 
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=5.7D0
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  filename(x)='CDF_WH_lnubb_5.7fb_10239' 

  x=x+1 
  S95_t1(x)%id=10796
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10796' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.45D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='CDF_WH_lnubb_9.45fb_10796' 

  x=x+1 
  S95_t1(x)%id=0874
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:1012.0874 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=5.3D0
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='D0_WH_lnubb_5.3fb_0874'

  x=x+1 
  S95_t1(x)%id=5613 
  S95_t1(x)%particle_x=Hneut 
  S95_t1(x)%expt='CDF' 
  S95_t1(x)%label='[hep-ex] arXiv:0906.5613 (CDF)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=2.7D0
  S95_t1(x)%xmin=100.0D0  
  S95_t1(x)%xmax=150.0D0 
  S95_t1(x)%sep=5.0D0  
  filename(x)='CDF_WH_lnubb_2.7fb_5613'

!----------------------- V H, H -> invisible -------------------------


  x=x+1
  S95_t1(x)%id=2013011
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2013-011'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=17.7D0
  S95_t1(x)%xmin=115.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2013011_Atlas_H-inv_17.7fb-1'

  x=x+1 
  S95_t1(x)%id=13018
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-018'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=24.7D0
  S95_t1(x)%xmin=105.0D0
  S95_t1(x)%xmax=145.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='13018_CMS_ZH-inv_24.7fb-1'


!----------------------- VBF, H -> invisible -------------------------


  x=x+1 
  S95_t1(x)%id=13013
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-013'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=19.6D0
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=400.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='13013_CMS_VBF-inv_19.6fb-1'

!----------------------- H -> W W -------------------------

  x=x+1 
  S95_t1(x)%id=5757
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 5757' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=3.0D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=115.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_ppH_WW_ll_3.0fb_5757' 

  x=x+1 
  S95_t1(x)%id=3930
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='[hep-ex] arXiv:0809.3930 (CDF)' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=3.0D0
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='CDF_ggH_WW_3.0fb_3930'
 
!  x=x+1 
!  S95_t1(x)%id=3216
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:1005.3216 (TEVNPHWG)'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=5.4D0
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=300.0D0
!  S95_t1(x)%sep=5.0D0 
!  filename(x)='CDF_D0_combined_gg-H-WW_4.8-5.4fb_3216_bayesian_interpol'

  !x=x+1 
  !S95_t1(x)%id=10102
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt='CDF' 
  !S95_t1(x)%label='CDF Note 10102' 
  !S95_t1(x)%xmin=110.0D0  
  !S95_t1(x)%xmax=200.0D0 
  !S95_t1(x)%sep=5.0D0  
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_H-WW_5.3fb_10102' 

!  x=x+1 
!  S95_t1(x)%id=6221
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6219' !this note has two results in it, both can not have the id 6219
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.1D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=300.0D0
!  S95_t1(x)%sep=5.0D0 
!  filename(x)='D0_H-WW_8.1fb_6221_interpol'

  x=x+1 
  S95_t1(x)%id=6276
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6276'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_H-VV_9.7fb_6276'
  
  x=x+1 
  S95_t1(x)%id=6301
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6301'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=115.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='D0_VH_VWW_9.7fb_6301'
    
  x=x+1
  S95_t1(x)%id=10600
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10599'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=8.2D0
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='CDF_ggH-WW_8.2fb_10600_interpol'
  
  x=x+1 
  S95_t1(x)%id=10599
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10599'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=8.2D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_H-WW_8.2fb_10599'

  x=x+1 
  S95_t1(x)%id=4468 
  S95_t1(x)%particle_x=Hneut 
  S95_t1(x)%expt='CDF' 
  S95_t1(x)%label='[hep-ex] arXiv:1001.4468 (CDF)' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=4.8D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0  
  S95_t1(x)%xmax=200.0D0 
  S95_t1(x)%sep=5.0D0  
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_H-WW_4.8fb_4468_interpol' 

!  x=x+1 
!  S95_t1(x)%id=5871
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 5871'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=4.2D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_H-WW_llnunu_4.2fb_5871'

  !x=x+1 
  !S95_t1(x)%id=6082
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt=' D0' 
  !S95_t1(x)%label='D0 Note 6082'
  !S95_t1(x)%xmin=115.0D0  
  !S95_t1(x)%xmax=200.0D0 
  !S95_t1(x)%sep=5.0D0  
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='D0_H-VV_6.7fb_6082' 


!  x=x+1 
!  S95_t1(x)%id=6219
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6219'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.1D0
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_H-VV_8.1fb_6219'


!  x=x+1 
!  S95_t1(x)%id=6179
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6179'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=7.3D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='D0_H-WW-mutau_7.3fb_6179' 

  x=x+1 
  S95_t1(x)%id=6302
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6302'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=115.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='D0_H-WW_9.7fb_6302' 

  x=x+1
  S95_t1(x)%id=6183
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6183'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=8.2D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=130.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0
  S95_t1(x)%deltax=0.0D0  
  filename(x)='D0_SM_combined_6183'

  x=x+1 
  S95_t1(x)%id=4481
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:1001.4481 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=5.4D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=115.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_H-WW_5.4fb_4481' 

 ! x=x+1 
!  S95_t1(x)%id=4162
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:1001.4162 (TEVNPHWG)'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=5.4D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=130.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_D0_SM_combined_H-WW_4.8-5.4fb_4162'


  x=x+1
  S95_t1(x)%id=3331
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='TCB'
  S95_t1(x)%label='[hep-ex] arXiv:1108.3331 (TEVNPHWG)'!CDF note 10608, D0 Note 6230
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=8.2D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_D0_combined_gg-H-WW_8.2fb_3331_bayesian_interpol'

!  x=x+1
!  S95_t1(x)%id=2011134
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-134'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=1.7D0  
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=300.0D0
!  S95_t1(x)%sep=1.0D0 
!  filename(x)='2011134_Atlas_H-WW-lnulnu_1.7fb-1_interpol'
  
  x=x+1
  S95_t1(x)%id=2012012
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-012'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='2012012_Atlas_H-WW-lnulnu_4.7fb-1'


!  x=x+1
!  S95_t1(x)%id=2012158
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2012-158'
!  S95_t1(x)%energy=8.0D0
!  S95_t1(x)%lumi=13.0D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='2012158_Atlas_H-WW-enumunu_13fb-1'

  x=x+1
  S95_t1(x)%id=2013030
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2013-030'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.0D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='2013030_Atlas_H-WW-lnulnu_25fb-1'


!  x=x+1 
!  S95_t1(x)%id=5429
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='[hep-ex] arXiv: 1102.5429(CMS)'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=0.036D0
!  S95_t1(x)%xmin=130.0D0
!  S95_t1(x)%xmax=400.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='5429_CMS_H-WW_36pb-1_rat_interpol'

  x=x+1
  S95_t1(x)%id=2577
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='[hep-ex] arxiv:1112.2577'
  S95_t1(x)%desc='pp->h + X->W W* + X ->l l nu nu'
  S95_t1(x)%energy=7.D0
  S95_t1(x)%lumi=2.05D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='2577_Atlas_H-WW-lnulnu_2.05fb-1'
    
!  x=x+1 
!  S95_t1(x)%id=1489
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='[hep-ex] arxiv:1202.1489 (CMS)'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.6D0
!  S95_t1(x)%SMlike=1  
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11024_CMS_H-WW-lnulnu_4.6fb-1_interpol'  

!  x=x+1 
!  S95_t1(x)%id=12042
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-12-042'
!  S95_t1(x)%energy=8.0D0
!  S95_t1(x)%lumi=16.0D0
!  S95_t1(x)%SMlike=1  
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='12042_CMS_H-WW-lnulnu_16.0fb-1'  

  x=x+1 
  S95_t1(x)%id=13003
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-003'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.0D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=5D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='13003_CMS_H-WW-lnulnu_25fb-1'  

  x=x+1 
  S95_t1(x)%id=13022
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-022'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.4D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=5D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='13022_CMS_VBF-WW_25.4fb-1'  



!--------------- H -> VV -> l nu l nu ----------------

  x=x+1
  S95_t1(x)%id=3357
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='[hep-ex] arXiv:1109.3357 (ATLAS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=1.04D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=200.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=20.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='3357_Atlas_H-ZZ-llnunu_1.04fb-1'

!  x=x+1
!  S95_t1(x)%id=2011148
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-148'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=2.05D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=200.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=20.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='2011148_Atlas_H-ZZ-lnulnu_2.05fb-1'

  x=x+1
  S95_t1(x)%id=2012016
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-016'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=200.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='2012016_Atlas_H-ZZ-llnunu_4.7fb-1'

  x=x+1
  S95_t1(x)%id=3478
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='[hep-ex] arxiv:1202.3478 (CMS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.6D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=250.0D0
  S95_t1(x)%xmax=590.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='11026_CMS_H-ZZ-llnunu_4.6fb-1'

!----------------------- H -> W W -----------------------------

!  x=x+1 
!  S95_t1(x)%id=2011052
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-052'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=0.035D0
!  S95_t1(x)%xmin=220.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=20.0D0 
!  filename(x)='2011052_Atlas_H-WW_35pb-1'

  x=x+1 
  S95_t1(x)%id=3615
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='[hep-ex] arXiv:1109.3615 (ATLAS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=0.035D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=240.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=20.0D0 
  S95_t1(x)%deltax=0.0D0     
  filename(x)='3615_Atlas_H-WW-lnuqq_1.04fb-1'  

  x=x+1 
  S95_t1(x)%id=2012018
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-018'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.7D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=300.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0     
  filename(x)='2012018_Atlas_H-WW-lnuqq_4.7fb-1'  

  x=x+1 
  S95_t1(x)%id=12046
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-12-046'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=17D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=170.0D0
  S95_t1(x)%xmax=580.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0     
  filename(x)='12046_CMS_H-WW-lnuqq_17fb-1'  


!------------------------- H -> Z Z --------------------------

  x=x+1
  S95_t1(x)%id=5064
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='[hep-ex] arXiv:1108.5064 (ATLAS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=1.04D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=200.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='5064_Atlas_H-ZZ-llqq_1.04fb-1'

  x=x+1
  S95_t1(x)%id=2012017
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-017'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.7D0
  S95_t1(x)%SMlike=1  
  S95_t1(x)%xmin=200.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='2012017_Atlas_H-ZZ-llqq_4.7fb-1'

  x=x+1
  S95_t1(x)%id=14161
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='[hep-ex] arXiv:1202.1416 (CMS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.6D0
  S95_t1(x)%SMlike=1 
  S95_t1(x)%xmin=130.0D0
  S95_t1(x)%xmax=164.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='110271_CMS_H-ZZ-llqq_4.6fb-1'

  x=x+1
  S95_t1(x)%id=14162
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='[hep-ex] arXiv:1202.1416 (CMS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.6D0
  S95_t1(x)%SMlike=1 
  S95_t1(x)%xmin=200.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='110272_CMS_H-ZZ-llqq_4.6fb-1' 

  x=x+1 
  S95_t1(x)%id=1415
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='[hep-ex] arXiv:1202.1415 (ATLAS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.8D0
  S95_t1(x)%SMlike=1 
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='1415_Atlas_H-ZZ-4l_4.8fb-1'

  x=x+1 
  S95_t1(x)%id=2012092
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-092'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=10.6D0
  S95_t1(x)%SMlike=1 
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2012092_Atlas_H-ZZ-4l_10.6fb-1'

  x=x+1 
  S95_t1(x)%id=20130131
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2013-013'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.0D0
  S95_t1(x)%SMlike=1 
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=180.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2013013-1_Atlas_H-ZZ-4l_incl_25fb-1'

  x=x+1 
  S95_t1(x)%id=20130132
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2013-013'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=21.0D0
  S95_t1(x)%xmin=200.0D0
  S95_t1(x)%xmax=1000.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='2013013-2_Atlas_H-ZZ-4l_ggF_21fb-1'

  x=x+1 
  S95_t1(x)%id=20130133
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2013-013'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=21.0D0
  S95_t1(x)%xmin=200.0D0
  S95_t1(x)%xmax=1000.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='2013013-3_Atlas_H-ZZ-4l_VBFVH_21fb-1'


  x=x+1
  S95_t1(x)%id=1997
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='[hep-ex] arxiv:1202.1997 (CMS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.7D0
  S95_t1(x)%SMlike=1 
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='11025_CMS_H-ZZ-4l_4.7fb-1'

!  x=x+1
!  S95_t1(x)%id=12041
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-12-041'
!  S95_t1(x)%energy=8.0D0
!  S95_t1(x)%lumi=17.3D0
!  S95_t1(x)%SMlike=1 
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='12041_CMS_H-ZZ-4l_17.3fb-1'

  x=x+1
  S95_t1(x)%id=130021
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-002'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.0D0
  S95_t1(x)%SMlike=1 
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='13002-1_CMS_H-ZZ-4l_lowm_25fb-1'

  x=x+1
  S95_t1(x)%id=130022
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-002'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.0D0
  S95_t1(x)%SMlike=1 
  S95_t1(x)%xmin=150.0D0
  S95_t1(x)%xmax=1000.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='13002-2_CMS_H-ZZ-4l_highm_25fb-1'


!------------------------- SM combined -----------------------

!  x=x+1 
!  S95_t1(x)%id=8961
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:0712.2383v1(TEVNPHWG)' !CDF Note 8961, D0 Note 5536
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=1.9D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_D0_SM_combined_1.0-1.9fb_8961-5536_interpol'

  x=x+1
  S95_t1(x)%id=10439
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10439'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=6.0D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_H-tautau_6.0fb_10439'

  !x=x+1 
  !S95_t1(x)%id=10133
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt='CDF' 
  !S95_t1(x)%label='CDF Note 10133' 
  !S95_t1(x)%xmin=100.0D0  
  !S95_t1(x)%xmax=150.0D0 
  !S95_t1(x)%sep=5.0D0
  !S95_t1(x)%deltax=0.0D0   
  !filename(x)='CDF_H-tautau_2.3fb_10133'

 ! x=x+1 
!  S95_t1(x)%id=4800
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='[hep-ex] arXiv:0903.4800 (D0)'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=1.0D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=105.0D0
!  S95_t1(x)%xmax=145.0D0
!  S95_t1(x)%sep=10.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_SM_tautau_1.0fb_4800'
!  
  x=x+1 
  S95_t1(x)%id=5845
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 5845'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=4.9D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=105.0D0
  S95_t1(x)%xmax=145.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_tautauqq_4.9fb_5845'

 ! x=x+1 
!  S95_t1(x)%id=9290
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:0804.3423(TEVNPHWG)' !CDF Note 9290, D0 Note 5645
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=2.4D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_D0_SM_combined_1.0-2.4fb_9290-5645_interpol'

 ! x=x+1 
!  S95_t1(x)%id=9465
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:0808.0534(TEVNPHWG)' !CDF Note 9465, D0 Note 5754
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=3.0D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=155.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  !filename(x)='CDF_D0_SM_combined_3.0fb_9465-5754_CLs'
!  filename(x)='CDF_D0_SM_combined_3.0fb_9465-5754_bayesian'

!  x=x+1 
!  S95_t1(x)%id=0598
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='[hep-ex] arXiv:0712.0598v1 (D0)'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=0.44D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_SM_combined_0598_interpol'

  x=x+1 
  S95_t1(x)%id=9999
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 9999'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=4.8D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_SM_combined_2.0-4.8fb_9999_interpol'

!  x=x+1 
!  S95_t1(x)%id=6008
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6008'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=5.4D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_SM_combined_6008'

  x=x+1 
  S95_t1(x)%id=6305
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6305'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=7.3D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=105.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_H-tautau_7.3fb_6305'

!  x=x+1 
!  S95_t1(x)%id=9998
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:0911.3930 (TEVNPHWG)'!CDF Note 9998, D0 Note 5983
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=5.4D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_D0_SM_combined_2.1-5.4fb_9998-5983_bayesian'

!  x=x+1 
!  S95_t1(x)%id=6096
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:1007.4587 (TEVNPHWG)' !CDF Note 10241, D0 Note 6096
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=6.7D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_D0_SM_combined_2.1-6.7fb_10241-6096_bayesian'

  x=x+1 
  S95_t1(x)%id=10010
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10010' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=4D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_hadr_4fb_10010'

  x=x+1
  S95_t1(x)%id=6171
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6171'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=4.3D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=105.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0                        
  filename(x)='D0_H-tautaujj-4.3fb_6171'

!  x=x+1 
!  S95_t1(x)%id=6229
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6229'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.6D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='D0_SM_combined_6229'
 
  x=x+1 
  S95_t1(x)%id=6304
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6304'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_SM_combined_6304'
  
  x=x+1 
  S95_t1(x)%id=10500
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10500'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=6.2D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_VH_tautau_6.2fb_10500'

  x=x+1 
  S95_t1(x)%id=10573
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10573'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=8.2D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=120.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_H-WW_8.2fb_10573'

  x=x+1 
  S95_t1(x)%id=6286
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6286'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=7.0D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_tautaumu_7.0fb_6286'
  
 ! x=x+1 
!  S95_t1(x)%id=3233
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:1103.3233v2 (TEVNPHWG)' !CDF Note 10432, D0 Note 6183
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.2D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=130.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  !filename(x)='CDF_D0_SM_combined_4.3-8.2fb_3233_CLs'
!  filename(x)='CDF_D0_SM_combined_4.3-8.2fb_3233_bayesian'
!
!  x=x+1 
!  S95_t1(x)%id=10606
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:1107.5518 (TEVNPHWG)' !CDF Note 10606, D0 Note 6226
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.6D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_D0_SM_combined_8.6fb_10606'

  x=x+1 
  S95_t1(x)%id=10884
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='TCB'
  S95_t1(x)%label='[hep-ex] arXiv:1207.0449 (TEVNPHWG)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=10.D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='10884_CDF_D0_SM_combined_10fb-1'

!  x=x+1 
!  S95_t1(x)%id=10607
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='TCB'
!  S95_t1(x)%label='[hep-ex] arXiv:1107.5518 (TEVNPHWG)' !CDF Note 10606, D0 Note 6226
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.6D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='CDF_D0_combined_h-bb_8.6fb_10607'

  x=x+1 
  S95_t1(x)%id=6436
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='TCB'
  S95_t1(x)%label='[hep-ex] arXiv:1207.6436 (TEVNPHWG)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='1207.6436_CDF_D0_combined_h-bb_9.7fb-1'


!  x=x+1 
!  S95_t1(x)%id=2011112
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-112'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=1.21D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2011112_Atlas_SMcombined_1.04-1.21fb-1_interpol'
!  
!  x=x+1
!  S95_t1(x)%id=2011135
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-135'
!    S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=2.3D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2011135_Atlas_SMcombined_1-2.3fb-1_interpol' 
!      
!  x=x+1 
!  S95_t1(x)%id=2011163
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-163'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.9D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=0.5D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2011163_Atlas_SMcombined_4.9fb-1'  
!  
  x=x+1
  S95_t1(x)%id=1408
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='(hep-ex) arxiv:1202.1408 (ATLAS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.9D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='1408_Atlas_SMcombined_4.9fb-1'   

  x=x+1
  S95_t1(x)%id=2012019
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-019'
  S95_t1(x)%desc='(p p)->h+...' 
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.9D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='2012019_Atlas_SMcombined_4.9fb-1'   

  x=x+1
  S95_t1(x)%id=7214
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='(hep-ex) arXiv:1207.7214 (ATLAS)'
  S95_t1(x)%desc='(p p)->h+...' 
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=10.5D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=0.1D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='7214_Atlas_SMcombined_10.5fb-1' 

  x=x+1
  S95_t1(x)%id=2011157
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2011-157, CMS-PAS-HIG-11-023'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=2.3D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='2011157_Atlas_CMS_SMcombined_2.7fb-1'
  
!  x=x+1
!  S95_t1(x)%id=11011
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-011'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=1.1D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=10.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11011_CMS_SMcomb_1.1fb-1_interpol'

  x=x+1
  S95_t1(x)%id=1488
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='[hep-ex] arxiv:1202.1488 (CMS)'
  S95_t1(x)%energy=7.D0
  S95_t1(x)%lumi=4.8D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='1488_CMS_SMcombined_4.8fb-1'
  
!  x=x+1
!  S95_t1(x)%id=12008
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-12-008'
!  S95_t1(x)%desc='(p p)->h+...' 
!  S95_t1(x)%energy=7.D0
!  S95_t1(x)%lumi=4.8D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='12008_CMS_SMcombined_4.8fb-1'

  x=x+1
  S95_t1(x)%id=12045
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-12-045'
  S95_t1(x)%desc='(p p)->h+...' 
  S95_t1(x)%energy=8.D0
  S95_t1(x)%lumi=17.3D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=600.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='12045_CMS_SMcombined_17.3fb-1'

!----------------- b H -> 3 b jets ------------------- 
      
!  x=x+1 
!  S95_t1(x)%id=10105
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CDF'
!  S95_t1(x)%label='CDF Note 10105' 
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=2.2D0
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=210.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='CDF_bbH_bb_2.2fb_10105'

!  x=x+1 
!  S95_t1(x)%id=5726
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 5726' 
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=2.6D0
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=220.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='D0_bbH_bb_2.6fb_5726_sq_interpol'

  x=x+1 
  S95_t1(x)%id=1931
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:1011.1931 (D0)' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=5.2D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='D0_bbH_bb_5.2fb_1931'

  x=x+1
  S95_t1(x)%id=4782
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='arXiv:1106.4782 (CDF)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=2.6D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=350.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='CDF_bbH_bb_2.6fb_4782'

!  x=x+1 
!  S95_t1(x)%id=12033
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-12-033'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.8D0
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=350.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='12033_CMS_bbH-bbbb_4.8fb-1'


!-------------------- H -> tau tau -------------------------

  x=x+1 
  S95_t1(x)%id=4555
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:1106.4555 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=5.4D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='D0_H-tautau_5.4fb_4555_interpol' 

!  x=x+1 
!  S95_t1(x)%id=2491
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='[hep-ex] arXiv:0805.2491 (D0)' 
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=300.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='D0_ppH_tautau_1fb_2491_interpol'

  x=x+1 
  S95_t1(x)%id=1014
  S95_t1(x)%particle_x=Hneut 
  S95_t1(x)%expt='CDF' 
  S95_t1(x)%label='[hep-ex] arXiv:0906.1014 (CDF)' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=1.8D0
  S95_t1(x)%xmin=90.0D0  
  S95_t1(x)%xmax=250.0D0 
  S95_t1(x)%sep=10.0D0  
  filename(x)='CDF_H-tautau_1.8fb_1014_interpol'

!  x=x+1 
!  S95_t1(x)%id=5740
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 5740'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=2.2D0
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=300.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='D0_H-tautau_2.2fb_5740_interpol' 

  x=x+1 
  S95_t1(x)%id=3363
  S95_t1(x)%particle_x=Hneut 
  S95_t1(x)%expt='TCB' 
  S95_t1(x)%label='[hep-ex] arXiv:1003.3363 (TEVNPHWG)' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=2.2D0
  S95_t1(x)%xmin=90.0D0  
  S95_t1(x)%xmax=200.0D0 
  S95_t1(x)%sep=10.0D0  
  filename(x)='CDF_D0_combinedH-tautau_2.2fb_3363_CLs'

!  x=x+1
!  S95_t1(x)%id=2011133
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-133'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=1.06D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=140.0D0
!  S95_t1(x)%sep=5.0D0
!  S95_t1(x)%deltax=0.0D0    
!  filename(x)='2011133_Atlas_H-tautau-ll4nu_1.06fb-1'

!  x=x+1
!  S95_t1(x)%id=2012014
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2012-014'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.7D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0
!  S95_t1(x)%deltax=0.0D0    
!  filename(x)='2012014_Atlas_H-tautau_4.7fb-1'

  x=x+1
  S95_t1(x)%id=2012160
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-160'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=17.6D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0
  S95_t1(x)%deltax=0.0D0    
  filename(x)='2012160_Atlas_H-tautau_17.6fb-1'

  x=x+1 
  S95_t1(x)%id=5003
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='[hep-ex] arXiv:1107.5003 (ATLAS)'!CERN-PH-EP-2011-104
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=0.036D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='5003_Atlas_ggh_h-tautau_36pb-1'
  
!  x=x+1
!  S95_t1(x)%id=2011132
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-132'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=1.06D0
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='2011132_Atlas_conserv_h-tautau_1.06fb-1_interpol'

  x=x+1
  S95_t1(x)%id=2012094
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-094'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.7D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=500.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='2012094_Atlas_H-tautau_4.7fb-1'
  
  
!  x=x+1
!  S95_t1(x)%id=110291
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-029'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.6D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=145.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11029_CMS_H-tautau_4.6fb-1_SM'

  x=x+1
  S95_t1(x)%id=12043
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-12-043'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=17.D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=145.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='12043_CMS_H_tautau_17fb-1'

  x=x+1 
  S95_t1(x)%id=10002
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='[hep-ex] arXiv: 1104:1619(CMS)'!CMS-PAS-HIG-10-002
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=0.036D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=500.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='10002_CMS_H-tautau_36pb-1_interpol'

  x=x+1 
  S95_t1(x)%id=110292
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='[hep-ex] arXiv: 1202.4083 (CMS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.6D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=500.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='11029_CMS_H-tautau_4.6fb-1_MSSM_interpol'

  x=x+1 
  S95_t1(x)%id=12050
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-12-050'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=17.D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=800.0D0
  S95_t1(x)%sep=1.0D0 
  filename(x)='12050_CMS_H-tautau_17fb-1'
  
!-------------------- H -> mu mu -------------------------
  
  x=x+1
  S95_t1(x)%id=2013010
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2013-010'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=21D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0     
  filename(x)='2013010_Atlas_H-mumu_21fb-1'
    
  
  
!------------------------ H W -> W W W --------------------------

  x=x+1 
  S95_t1(x)%id=5873
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 5873'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=3.6D0
  S95_t1(x)%xmin=120.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=20.0D0 
  filename(x)='D0_WH_WWW_llnunu_3.6fb_5873'

  x=x+1 
  S95_t1(x)%id=7307
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 7307 vs 3' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=2.7D0
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='CDF_WH_WWW_2.7fb_7307vs3.0'  
  !filename(x)='CDF_WH_WWW_1.9fb_7307' 

  x=x+1 
  S95_t1(x)%id=1268
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:1107.1268 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=5.3D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=115.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_VH_ll_5.3fb_1268'

!   x=x+1 
!  S95_t1(x)%id=11034
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-034'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.6D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11034_CMS_WH-WWW_4.6fb-1_interpol'

!   x=x+1 
!  S95_t1(x)%id=12039
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-12-039'
!  S95_t1(x)%energy=8.0D0
!  S95_t1(x)%lumi=10.D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='12039_CMS_WH-WWW_10fb-1'

   x=x+1 
  S95_t1(x)%id=13009
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-009'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='13009_CMS_WH-WWW_25fb-1'



   x=x+1 
  S95_t1(x)%id=12006
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-12-006'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=140.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='12006_CMS_WH-Wtautau_4.7fb-1'

   x=x+1 
  S95_t1(x)%id=12051
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-12-051'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=17.D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=145.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='12051_CMS_VH_Vtautau_17fb-1'


   x=x+1 
  S95_t1(x)%id=2012078
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-078'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=300.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2012078_Atlas_WH_WWW_4.7fb-1'


!-------------------- H -> gamma gamma --------------------

!  x=x+1 
!  S95_t1(x)%id=6177
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6177'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=8.2D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=2.5D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='D0_gaga_8.2fb_6177' 
  
  x=x+1 
  S95_t1(x)%id=6295
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6295'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=9.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=2.5D0 
  S95_t1(x)%deltax=0.0D0   
  filename(x)='D0_gaga_9.7fb_6295'   

  x=x+1 
  S95_t1(x)%id=1887
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:0901.1887 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=2.7D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_Hgaga_2.7fb_1887'

  x=x+1 
  S95_t1(x)%id=10485
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10485' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=7.0D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_SM_Hgaga_7.0fb_10485'

  x=x+1
  S95_t1(x)%id=4960
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='TCB'
  S95_t1(x)%label='[hep-ex] arXiv:1107.4960 (TEVNPHWG)'!CDF Note 10510, D0 Note 6203
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=8.2D0
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_D0_combined_SM_Hgaga_8.2fb_4960'

  x=x+1 
  S95_t1(x)%id=1414
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='[hep-ex] arXiv:1202.1414 (ATLAS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.9D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=1.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='1414_Atlas_H-gaga_4.9fb-1'

!  x=x+1 
!  S95_t1(x)%id=2012091
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2012-091'
!  S95_t1(x)%energy=8.0D0
!  S95_t1(x)%lumi=10.8D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=0.5D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2012091_Atlas_H-gaga_10.8fb-1'

  x=x+1 
  S95_t1(x)%id=2012168
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-168'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=17.8D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=0.5D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2012168_Atlas_H-gaga_17.8fb-1'


!  x=x+1 
!  S95_t1(x)%id=1487
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='[hep-ex] arXiv:1202.1487 (CMS)'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.8D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=0.5D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='1487_CMS_H-gaga_4.8fb-1'

!   x=x+1 
!  S95_t1(x)%id=12001
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-12-001'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=4.76D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=0.5D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='12001_CMS_H-gaga_BDT_4.8fb-1'
  
!   x=x+1 
!  S95_t1(x)%id=12015
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-12-015'
!  S95_t1(x)%energy=8.0D0
!  S95_t1(x)%lumi=9.9D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=0.5D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='12015_CMS_H-gaga_9.9fb-1'  
! 
   x=x+1 
  S95_t1(x)%id=13001
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-13-001'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=0.5D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='13001_CMS_H-gaga_MVA_25fb-1'  
 
  
!-------------------- H -> gamma Z --------------------
  
!  x=x+1 
!  S95_t1(x)%id=13006
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-13-006'
!  S95_t1(x)%energy=8.0D0
!  S95_t1(x)%lumi=25.D0
!  S95_t1(x)%SMlike=1
!  S95_t1(x)%xmin=120.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=1.D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='13006_CMS_H-gaZ_25fb-1'  

  x=x+1 
  S95_t1(x)%id=13075515
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='[hep-ex] arXiv:1307.5515 (CMS)'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=24.6D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=120.0D0
  S95_t1(x)%xmax=160.0D0
  S95_t1(x)%sep=0.1D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='13075515_CMS_H-Zgamma_24.6fb-1'  
    
  x=x+1 
  S95_t1(x)%id=2013009
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2013-009'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=25.D0
  S95_t1(x)%SMlike=1
  S95_t1(x)%xmin=120.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=1.D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2013009_Atlas_H-gaZ_25fb-1'
  
  
!--------------------- b H -> b tau tau ---------------------------  
  
  x=x+1 
  S95_t1(x)%id=4885
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:1106.4885 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=7.3D0  
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=320.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='D0_Hb_tautaub_7.3fb_4885'

!  x=x+1 
!  S95_t1(x)%id=5985
!  S95_t1(x)%particle_x=Hneut 
!  S95_t1(x)%expt=' D0' 
!  S95_t1(x)%label='D0 Note 5985' 
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=2.7D0
!  S95_t1(x)%xmin=90.0D0  
!  S95_t1(x)%xmax=320.0D0 
!  S95_t1(x)%sep=10.0D0  
!  filename(x)='D0_bH_btautau_2.7fb_5985' 

!  x=x+1 
!  S95_t1(x)%id=5974
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 5974'
!  S95_t1(x)%energy=1.96D0
!  S95_t1(x)%lumi=3.7D0
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=300.0D0
!  S95_t1(x)%sep=10.0D0  
!  filename(x)='D0_Hb_tautaub_3.7fb_5974' 
  
  x=x+1 
  S95_t1(x)%id=6083
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 6083' 
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=4.3D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=200.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='D0_Hb_tautaub_4.3fb_6083'
  
!-------------------- t t H -> t t b b --------------------

  x=x+1 
  S95_t1(x)%id=5739
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='D0 Note 5739'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=2.1D0
  S95_t1(x)%xmin=105.0D0
  S95_t1(x)%xmax=155.0D0
  S95_t1(x)%sep=10.0D0 
  filename(x)='D0_ttH_ttbb_2.1fb_5739'  
  
  x=x+1 
  S95_t1(x)%id=10574
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 10574'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=7.5D0  
  S95_t1(x)%xmin=100.0D0
  S95_t1(x)%xmax=170.0D0
  S95_t1(x)%sep=5.0D0 
  filename(x)='CDF_ttH_ttbb_7.5fb_10574_interpol'  

  x=x+1 
  S95_t1(x)%id=2012135
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2012-135'
  S95_t1(x)%energy=7D0
  S95_t1(x)%lumi=4.7D0  
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=140.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2012135_Atlas_ttH_Hbb_4.7fb-1'  

  x=x+1 
  S95_t1(x)%id=12025
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-12-025'
  S95_t1(x)%energy=7D0
  S95_t1(x)%lumi=5.D0  
  S95_t1(x)%xmin=110.0D0
  S95_t1(x)%xmax=140.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='12025_CMS_ttH-ttbb_5fb-1'  

!-------------------- H -> Z gamma --------------------------

  x=x+1 
  S95_t1(x)%id=0611
  S95_t1(x)%particle_x=Hneut
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:0806.0611 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=1.1D0
  S95_t1(x)%xmin=120.0D0
  S95_t1(x)%xmax=320.0D0
  S95_t1(x)%sep=20.0D0 
  filename(x)='D0_H-Zgamma_1.0-1.1fb_0611'

  !x=x+1 
  !S95_t1(x)%id=6091
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 6091'
  !S95_t1(x)%xmin=115.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='D0_VH_ll_5.4fb_6091'

  !x=x+1 
  !S95_t1(x)%id=5858
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5858'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=2.5D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='D0_gaga_4.2fb_5858'

  !x=x+1 
  !S95_t1(x)%id=10065
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 10065' 
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=10.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_SM_Hgaga_5.4fb_10065' 

!  x=x+1 
!  S95_t1(x)%id=0968
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='[hep-ex] arXiv:0912.0968 (D0)'
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=320.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='D0_Hb_tautaub_2.7fb_0968'

!---------------- charged Higgs ------------------

  x=x+1 
  S95_t1(x)%id=1811
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:0908.1811 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=1.0D0
  S95_t1(x)%xmin=80.0D0
  S95_t1(x)%xmax=155.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_Hp_qq_1.0fb_1811_interpol'

  x=x+1 
  S95_t1(x)%id=1269
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='[hep-ex] arXiv:0907.1269 (CDF) lower mass'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=2.2D0
  S95_t1(x)%xmin=60.0D0
  S95_t1(x)%xmax=70.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_t-Hplusb_1269_lowmass'

  x=x+1 
  S95_t1(x)%id=1270
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='[hep-ex] arXiv:0907.1269 (CDF) higher mass'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=2.2D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=150.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_t-Hplusb_1269_highmass'

  x=x+1 
  S95_t1(x)%id=7712
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 7712'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=0.192D0
  S95_t1(x)%xmin=80.0D0
  S95_t1(x)%xmax=160.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='CDF_t-Hplusb_tauonicHM_192pb_7712'

  x=x+1 
  S95_t1(x)%id=8353
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='CDF'
  S95_t1(x)%label='CDF Note 8353'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=0.335D0
  S95_t1(x)%xmin=80.0D0
  S95_t1(x)%xmax=120.0D0
  S95_t1(x)%sep=10.0D0
  S95_t1(x)%deltax=0.0D0    
  filename(x)='CDF_t-Hplusb_8353'

  x=x+1 
  S95_t1(x)%id=1812
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt=' D0'
  S95_t1(x)%label='[hep-ex] arXiv:0908.1811 (D0)'
  S95_t1(x)%energy=1.96D0
  S95_t1(x)%lumi=1.0D0
  S95_t1(x)%xmin=80.0D0
  S95_t1(x)%xmax=155.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='D0_Hp_taunu_1.0fb_1811_interpol'

  x=x+1
  S95_t1(x)%id=2011094
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2011-094'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=0.035D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=130.0D0
  S95_t1(x)%sep=20.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2011094_Atlas_Hplus-cs_35pb-1'

!  x=x+1 
!  S95_t1(x)%id=2011138
!  S95_t1(x)%particle_x=Hplus
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-138'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=1.03D0
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=160.0D0
!  S95_t1(x)%sep=10.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2011138_Atlas_Hp-taunu_1.03fb-1'

!  x=x+1 
!  S95_t1(x)%id=2011151
!  S95_t1(x)%particle_x=Hplus
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-151'
!  S95_t1(x)%energy=7.0D0
!  S95_t1(x)%lumi=1.03D0
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=160.0D0
!  S95_t1(x)%sep=10.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2011151_Atlas_chargedH-taunu_1.03fb-1'

  x=x+1 
  S95_t1(x)%id=2760
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='[hep-ex] arXiv:1204.2760 (ATLAS)'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=4.6D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=160.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2760_Atlas_Hp_taunu_4.6fb-1'

  x=x+1 
  S95_t1(x)%id=2013090
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='ATL'
  S95_t1(x)%label='ATLAS-CONF-2013-090'
  S95_t1(x)%energy=8.0D0
  S95_t1(x)%lumi=19.5D0
  S95_t1(x)%xmin=90.0D0
  S95_t1(x)%xmax=160.0D0
  S95_t1(x)%sep=10.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='2013090_Atlas_Hp_taunu_19.6fb-1'


  x=x+1
  S95_t1(x)%id=11008
  S95_t1(x)%particle_x=Hplus
  S95_t1(x)%expt='CMS'
  S95_t1(x)%label='CMS-PAS-HIG-11-008'
  S95_t1(x)%energy=7.0D0
  S95_t1(x)%lumi=1.1D0
  S95_t1(x)%xmin=80.0D0
  S95_t1(x)%xmax=160.0D0
  S95_t1(x)%sep=5.0D0 
  S95_t1(x)%deltax=0.0D0 
  filename(x)='11008_CMS_Hplus-taunu_1.1fb-1_interpol'







!----------------------------------------------

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!   superseded

  !x=x+1     !superseded by arXiv:1011.1931
  !S95_t1(x)%id=3556
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='[hep-ex] arXiv:0805.3556 (D0)' 
  !S95_t1(x)%xmin=90.0D0
  !S95_t1(x)%xmax=220.0D0
  !S95_t1(x)%sep=10.0D0 
  !filename(x)='D0_bbH_bb_1.0fb_3556_interpol'

  !x=x+1     !superseded by arXiv:1009.3047
  !S95_t1(x)%id=3534
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt='CDF' 
  !S95_t1(x)%label='[hep-ex] arXiv:0908.3534 (CDF)' 
  !S95_t1(x)%xmin=100.0D0  
  !S95_t1(x)%xmax=150.0D0 
  !S95_t1(x)%sep=5.0D0  
  !filename(x)='CDF_ZH_llbb_2.7fb_3534' 

  !x=x+1    !superseded by CDF 10239
  !S95_t1(x)%id=10217
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 10217' 
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_WH_lnubb_5.6fb_10217' 

  !x=x+1   !superseded by CDF Note 10133
  !S95_t1(x)%id=9248
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9248' 
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_SMcomb_tautau_2.0fb_9248_interpol'

  !x=x+1  !superseded by arXiv:1012.0874
  !S95_t1(x)%id=1970
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='[hep-ex] arXiv:0808.1970 (D0)'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='D0_WH_lnubb_1.05fb_1970'

  !x=x+1 
  !S95_t1(x)%id=6095   
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 6095'
  !S95_t1(x)%xmin=115.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='D0_H-WW_lljj_5.4fb_6095'

  !x=x+1 
  !S95_t1(x)%id=5972 
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt=' D0' 
  !S95_t1(x)%label='D0 Note 5972' 
  !S95_t1(x)%xmin=100.0D0  
  !S95_t1(x)%xmax=150.0D0 
  !S95_t1(x)%sep=5.0D0  
  !filename(x)='D0_WH_lnubb_5.0fb_5972'

  !x=x+1 
  !S95_t1(x)%id=10068
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 10068' 
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_WH_lnubb_4.8fb_10068'

  !x=x+1
  !S95_t1(x)%id=5586
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5586' 
  !S95_t1(x)%xmin=105.0D0
  !S95_t1(x)%xmax=145.0D0
  !S95_t1(x)%sep=10.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='D0_ZH_nunubb_2.1fb_5586'

  !x=x+1 
  !S95_t1(x)%id=9891
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9891'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_VH_Metbb_3.6fb_9891'

  !x=x+1 
  !S95_t1(x)%id=5876
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5876'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='D0_ZH_llbb_4.2fb_5876'

  !x=x+1 
  !S95_t1(x)%id=9889 
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt='CDF' 
  !S95_t1(x)%label='CDF Note 9889' 
  !S95_t1(x)%xmin=100.0D0  
  !S95_t1(x)%xmax=150.0D0 
  !S95_t1(x)%sep=5.0D0  
  !filename(x)='CDF_ZH_llbb_4.1fb_9889' 

  !x=x+1 
  !S95_t1(x)%id=9284
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9284' 
  !S95_t1(x)%xmin=90.0D0
  !S95_t1(x)%xmax=210.0D0
  !S95_t1(x)%sep=10.0D0 
  !filename(x)='CDF_bbH_bb_1.9fb_9284'

  !x=x+1 
  !S95_t1(x)%id=9868
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9868' 
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_WH_lnubb_4.3fb_9868' !http://www-cdf.fnal.gov/physics/new/hdg/results/whlnubb_090814/WH4.3fb.html

  !x=x+1 
  !S95_t1(x)%id=9887 !superseded by arXiv:1001.4468
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt='CDF' 
  !S95_t1(x)%label='CDF Note 9887' 
  !S95_t1(x)%xmin=110.0D0  
  !S95_t1(x)%xmax=200.0D0 
  !S95_t1(x)%sep=5.0D0  
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_H-WW_4.8fb_9887_interpol' 

  !x=x+1 
  !S95_t1(x)%id=0024
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='[hep-ex] arXiv:0811.0024v1 (D0)'
  !S95_t1(x)%xmin=90.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=10.0D0 
  !filename(x)='D0_bH_btautau_328pb_0024'

  !x=x+1 
  !S95_t1(x)%id=1266
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='[hep-ex] arXiv:0808.1266 (D0)'
  !S95_t1(x)%xmin=105.0D0
  !S95_t1(x)%xmax=135.0D0
  !S95_t1(x)%sep=10.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='D0_VH_bb_0.93fb_1266'

  !x=x+1 
  !S95_t1(x)%id=0432
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='[hep-ex] arXiv:0802.0432 (CDF)'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=140.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_VH_Metbb_1fb_0432'

  !x=x+1 
  !S95_t1(x)%id=9674
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9674'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_SM_combined_2.0-3.0fb_9674_interpol'

  !x=x+1 
  !S95_t1(x)%id=5980
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt='TCB' 
  !S95_t1(x)%label='D0 Note 5980, CDF Note 9888' 
  !S95_t1(x)%xmin=100.0D0  
  !S95_t1(x)%xmax=200.0D0 
  !S95_t1(x)%sep=10.0D0  
  !filename(x)='CDF_D0_combinedH-tautau_2.2fb_5980' 

  !x=x+1 
  !S95_t1(x)%id=9714
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='TCB'
  !S95_t1(x)%label='arXiv:0903.4001 (TEVNPHWG), MH>=155GeV'!CDF Note 9713, D0 Note 5889
  !S95_t1(x)%xmin=155.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_D0_SM_combined_4.2fb_9713-5889_bayesian_highmassonly'

  !x=x+1 
  !S95_t1(x)%id=9897
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9897'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_SM_combined_2.0-4.8fb_9897_interpol'

  !x=x+1 
  !S95_t1(x)%id=8742
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'                        
  !S95_t1(x)%label='CDF Note 8742' ! table from CDF note 8742
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0
  !filename(x)='CDF_ZH_llbb_1fb_8742_interpol' 
  
!   x=x+1
!   S95_t1(x)%id=7081
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='CDF'
!   S95_t1(x)%label='arxiv 070810 (CDF)'
!   S95_t1(x)%xmin=110.0D0
!   S95_t1(x)%xmax=150.0D0
!   S95_t1(x)%sep=5.0D0
!     S95_t1(x)%deltax=0.0D0 
!   filename(x)='CDF_ZH_nunubb_1.7fb_070810_interpol'
  
! !   x=x+1
! !   S95_t1(x)%id=8958 
! !   S95_t1(x)%particle_x=Hneut
! !   S95_t1(x)%expt='CDF' 
! !   S95_t1(x)%label='CDF Note 8958' 
! !   S95_t1(x)%xmin=110.0D0
! !   S95_t1(x)%xmax=200.0D0
! !   S95_t1(x)%sep=10.0D0
! !   filename(x)='CDF_ggH_WW_1.9fb_8958' 
   
!   x=x+1
!   S95_t1(x)%id=8957
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='CDF'
!   S95_t1(x)%label='CDF Note 8957' 
!   S95_t1(x)%xmin=110.0D0
!   S95_t1(x)%xmax=150.0D0
!   S95_t1(x)%sep=5.0D0
!   filename(x)='CDF_WH_lnubb_1.7fb_8957_interpol'

!   x=x+1
!   S95_t1(x)%id=5489
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 5489' 
!   S95_t1(x)%xmin=120.0D0
!   S95_t1(x)%xmax=200.0D0
!   S95_t1(x)%sep=20.0D0
!   filename(x)='D0_ppH_WW_emu_0.6fb_5489'
  
  !x=x+1
  !S95_t1(x)%id=5537
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5537' 
  !S95_t1(x)%xmin=120.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=20.0D0
  !filename(x)='D0_ppH_WW_ll_1.7fb_5537'
  
!   x=x+1
!   S95_t1(x)%id=5482
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 5482' 
!   S95_t1(x)%xmin=105.0D0
!   S95_t1(x)%xmax=145.0D0
!   S95_t1(x)%sep=10.0D0 
!   filename(x)='D0_ZH_llbb_1.1fb_5482' 
 
!   x=x+1 
!   S95_t1(x)%id=5624
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 5624' 
!   S95_t1(x)%xmin=115.0D0
!   S95_t1(x)%xmax=200.0D0
!   S95_t1(x)%sep=5.0D0 
!   filename(x)='D0_ppH_WW_ll_2.3fb_5624'

  !x=x+1 
  !S95_t1(x)%id=5472
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5472' 
  !S95_t1(x)%xmin=105.0D0
  !S95_t1(x)%xmax=145.0D0
  !S95_t1(x)%sep=10.0D0 
  !!filename(x)='D0_WH_lnubb_1.7fb_5472' 
  !filename(x)='D0_WH_lnubb_1.7fb_5472_Fig.6a'
  
  !x=x+1 
  !S95_t1(x)%id=5502
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5502' 
  !S95_t1(x)%xmin=120.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=20.0D0 
  !filename(x)='D0_ppH_WW_ee_0.6fb_5502'  

  !x=x+1 
  !S95_t1(x)%id=5332
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5332' 
  !S95_t1(x)%xmin=120.0D0
  !S95_t1(x)%xmax=180.0D0
  !S95_t1(x)%sep=20.0D0 
  !filename(x)='D0_ppH_WW_mutau_1fb_5332' 
  
  !x=x+1 
  !S95_t1(x)%id=5485
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5485' 
  !S95_t1(x)%xmin=120.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=20.0D0 
  !filename(x)='D0_WH_WWW_llnunu_1fb_5485'

  !x=x+1 
  !S95_t1(x)%id=9071
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9071' 
  !S95_t1(x)%xmin=90.0D0
  !S95_t1(x)%xmax=250.0D0
  !S95_t1(x)%sep=10.0D0 
  !filename(x)='CDF_ppH_tautau_1.8fb_9071_interpol'

  !x=x+1 
  !S95_t1(x)%id=5331
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5331' 
  !S95_t1(x)%xmin=90.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=10.0D0 
  !filename(x)='D0_ppH_tautau_1fb_5331_interpol'

!   x=x+1 
!   S95_t1(x)%id=5503
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 5503' 
!   S95_t1(x)%xmin=100.0D0
!   S95_t1(x)%xmax=170.0D0
!   S95_t1(x)%sep=10.0D0 
!   filename(x)='D0_bbH_bb_0.9fb_5503_interpol'

  !x=x+1 
  !S95_t1(x)%id=8954
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 8954' 
  !S95_t1(x)%xmin=90.0D0
  !S95_t1(x)%xmax=210.0D0
  !S95_t1(x)%sep=10.0D0 
  !filename(x)='CDF_bbH_bb_1fb_8954'

!   x=x+1 
!   S95_t1(x)%id=9236
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='CDF'
!   S95_t1(x)%label='CDF Note 9236' 
!   S95_t1(x)%xmin=110.0D0
!   S95_t1(x)%xmax=200.0D0
!   S95_t1(x)%sep=10.0D0 
!   filename(x)='CDF_ggH_WW_2.4fb_9236'

!   x=x+1 
!   S95_t1(x)%id=9219
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='CDF'
!   S95_t1(x)%label='CDF Note 9219' 
!   S95_t1(x)%xmin=110.0D0
!   S95_t1(x)%xmax=150.0D0
!   S95_t1(x)%sep=5.0D0 
!   filename(x)='CDF_WH_lnubb_1.9fb_9219'    

!   x=x+1 
!   S95_t1(x)%id=5601
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 5601' 
!   S95_t1(x)%xmin=100.0D0
!   S95_t1(x)%xmax=150.0D0
!   S95_t1(x)%sep=10.0D0 
!   filename(x)='D0_Hgaga_2.27fb_5601'

  !taken out until we sort out a fermiophobic-ness test
  !x=x+1 
  !S95_t1(x)%id=1514
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='[hep-ex] arXiv:0803.1514v1 (D0)' 
  !S95_t1(x)%xmin=70.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=10.0D0 
  !filename(x)='D0_gaga_1.1fb_1514'
 
  
!   x=x+1 
!   S95_t1(x)%id=9166
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='CDF'
!   S95_t1(x)%label='CDF Note 9166' 
!   S95_t1(x)%xmin=110.0D0
!   S95_t1(x)%xmax=150.0D0
!   S95_t1(x)%sep=5.0D0 
!    S95_t1(x)%deltax=0.0D0 
!   filename(x)='CDF_VH_Metbb_1.97fb_9166_interpol'  

  !x=x+1 
  !S95_t1(x)%id=9463
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9463' 
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_WH_lnubb_2.7fb_9463'

  !x=x+1 
  !S95_t1(x)%id=5570
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5570' 
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='D0_ZH_llbb_2.3fb_5570'

  !x=x+1 
  !S95_t1(x)%id=4493
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='[hep-ex] arXiv:0807.4493 (CDF)' 
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_ZH_llbb_1fb_4493_interpol'

  !x=x+1 
  !S95_t1(x)%id=9475
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9475' 
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_ZH_llbb_2.4fb_9475'

  !x=x+1 
  !S95_t1(x)%id=5737
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 5737'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='D0_gaga_2.68fb_5737'

  !x=x+1
  !S95_t1(x)%id=9483
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF note 9483'
  !S95_t1(x)%xmin=105.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_VH_Metbb_2.1fb_9483'

  !x=x+1 
  !S95_t1(x)%id=9500
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9500'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_ggH_WW_3.0fb_9500_interpol'


  !x=x+1 
  !S95_t1(x)%id=9713
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='TCB'
  !S95_t1(x)%label='CDF Note 9713, D0 Note 5889'
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_D0_SM_combined_4.2fb_9713-5889_bayesian'

!   x=x+1 
!   S95_t1(x)%id=1024
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='CDF'
!   S95_t1(x)%label='CDF Note XXXX, zhllbb_081024'
!   S95_t1(x)%xmin=100.0D0
!   S95_t1(x)%xmax=150.0D0
!   S95_t1(x)%sep=5.0D0 
!   filename(x)='CDF_ZH_llbb_2.7fb_1024'

  !x=x+1 
  !S95_t1(x)%id=9596
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9596' 
  !S95_t1(x)%xmin=100.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_WH_lnubb_2.7fb_9596' 

!   x=x+1 
!   S95_t1(x)%id=5828
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 5828'
!   S95_t1(x)%xmin=100.0D0
!   S95_t1(x)%xmax=150.0D0
!   S95_t1(x)%sep=5.0D0 
!   filename(x)='D0_WH_lnubb_2.7fb_5828'

  !x=x+1 
  !S95_t1(x)%id=9642
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9642'
  !S95_t1(x)%xmin=105.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_VH_Metbb_2.1fb_9642'

!   x=x+1 
!   S95_t1(x)%id=9022
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='CDF'
!   S95_t1(x)%label='CDF Note 9764,"old" ggH XS'
!   S95_t1(x)%xmin=110.0D0
!   S95_t1(x)%xmax=200.0D0
!   S95_t1(x)%sep=5.0D0 
!   S95_t1(x)%deltax=0.0D0 
!   filename(x)='CDF_H-WW_3.6fb_9022_interpol'

  !NEW! not yet in OB's code
  !x=x+1 
  !S95_t1(x)%id=9023
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 9764,"new" ggH XS'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
!   S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_H-WW_3.6fb_9023_interpol'

!   x=x+1 
!   S95_t1(x)%id=3493
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='CDF'
!   S95_t1(x)%label='[hep-ex] arXiv:0803.3493v1 (CDF)'
!   S95_t1(x)%xmin=110.0D0
!   S95_t1(x)%xmax=150.0D0
!   S95_t1(x)%sep=5.0D0 
!   filename(x)='CDF_WH_lnubb_1fb_3493_interpol'

  !x=x+1 
  !S95_t1(x)%id=0710 
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt='CDF' 
  !S95_t1(x)%label='CDF Note XXXX,hwwmenn_090710' 
  !S95_t1(x)%xmin=110.0D0  
  !S95_t1(x)%xmax=200.0D0 
  !S95_t1(x)%sep=5.0D0  
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_H-WW_4.8fb_0710_interpol' 

  !x=x+1 
  !S95_t1(x)%id=5984
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt=' D0' 
  !S95_t1(x)%label='D0 Note 5984,high mass' 
  !S95_t1(x)%xmin=155.0D0  
  !S95_t1(x)%xmax=200.0D0 
  !S95_t1(x)%sep=5.0D0  
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='D0_SM_combined_5984_highmassonly' 

  !x=x+1 
  !S95_t1(x)%id=3155 
  !S95_t1(x)%particle_x=Hneut 
  !S95_t1(x)%expt='CDF' 
  !S95_t1(x)%label='[hep-ex] arXiv:0905.3155 (CDF)' 
  !S95_t1(x)%xmin=110.0D0  
  !S95_t1(x)%xmax=150.0D0 
  !S95_t1(x)%sep=5.0D0  
  !filename(x)='CDF_WH_lnubb_1.9fb_3155_interpol' 

!   x=x+1 
!   S95_t1(x)%id=9997
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt='TCB'
!   S95_t1(x)%label='arXiv:0911.3930 (TEVNPHWG), MH>=160GeV'!CDF Note 9998, D0 Note 5983
!   S95_t1(x)%xmin=160.0D0
!   S95_t1(x)%xmax=200.0D0
!   S95_t1(x)%sep=5.0D0 
!   S95_t1(x)%deltax=0.0D0 
!   filename(x)='CDF_D0_SM_combined_2.1-5.4fb_9998-5983_bayesian_highmassonly'

  !x=x+1 
  !S95_t1(x)%id=6039
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='TCB'
  !S95_t1(x)%label='CDF Note 10101, D0 Note 6039'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=260.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_D0_combined_gg-H-WW_4.8-5.4fb_6039_bayesian_interpol'

  !x=x+1 
  !S95_t1(x)%id=6006
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt=' D0'
  !S95_t1(x)%label='D0 Note 6006'
  !S95_t1(x)%xmin=115.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='D0_H-WW_ll_5.4fb_6006'


  !x=x+1 
  !S95_t1(x)%id=2011026
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-026'
  !S95_t1(x)%xmin=200.0D0
  !S95_t1(x)%xmax=600.0D0
  !S95_t1(x)%sep=20.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='2011026_Atlas_H-VV_35pb-1'

  !x=x+1     !superseded by 2011085
  !S95_t1(x)%id=2011025
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-025'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=140.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='2011025_Atlas_H-gaga_37.6pb-1_interpol'

  !x=x+1 
  !S95_t1(x)%id=2011005
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-005'
  !S95_t1(x)%xmin=120.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='2011005_Atlas_H-WW_35pb-1_interpol'

  !x=x+1
  !S95_t1(x)%id=10433
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 10432'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=300.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='CDF_ggH-WW_7.1fb_10433_SMnormalized'
  
  !x=x+1 
  !S95_t1(x)%id=10432
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='CDF'
  !S95_t1(x)%label='CDF Note 10432'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=200.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='CDF_H-WW_7.1fb_10432'  
  
!  x=x+1 
!  S95_t1(x)%id=6170
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6170'
!  S95_t1(x)%xmin=100.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='D0_ZH_bb_6.2fb_6170' 

!  x=x+1 
!  S95_t1(x)%id=6182
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt=' D0'
!  S95_t1(x)%label='D0 Note 6182'
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=200.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='D0_H-VV_8.1fb_6182' 

  !x=x+1 superseded by 2011131
  !S95_t1(x)%id=2011048
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-048'
  !S95_t1(x)%xmin=120.0D0
  !S95_t1(x)%xmax=600.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='2011048_Atlas_H-ZZ-4l_40pb-1_interpol'

  !x=x+1 
  !S95_t1(x)%id=2011085
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-085'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=140.0D0
  !S95_t1(x)%sep=1.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='2011085_Atlas_H-gaga_209pb-1'

!  x=x+1 
!  S95_t1(x)%id=2748
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='[hep-ex] arXiv:1106.2748 (ATLAS)'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2748_Atlas_SMcombined_40pb-1'

  !x=x+1 
  !S95_t1(x)%id=11002
  !S95_t1(x)%particle_x=Hplus
  !S95_t1(x)%expt='CMS'
  !S95_t1(x)%label='CMS-PAS-HIG-11-002'
  !S95_t1(x)%xmin=80.0D0
  !S95_t1(x)%xmax=160.0D0
  !S95_t1(x)%sep=5.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='11002_CMS_t-Hplusb_36pb-1'

  !x=x+1 
  !S95_t1(x)%id=2011020
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-020, lower masses'
  !S95_t1(x)%xmin=6.1D0
  !S95_t1(x)%xmax=9.0D0
  !S95_t1(x)%sep=0.1D0 
  !filename(x)='2011020_Atlas_H-mumu_lowerMa_35pb-1'

  !x=x+1 
  !S95_t1(x)%id=2011021
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-020, higher masses'
  !S95_t1(x)%xmin=11.0D0
  !S95_t1(x)%xmax=11.9D0
  !S95_t1(x)%sep=0.1D0 
  !filename(x)='2011021_Atlas_H-mumu_higherMa_35pb-1'

!  x=x+1
!  S95_t1(x)%id=11003
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-003'
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11003_CMS_H-WW_1.1fb-1_interpol'

!  x=x+1
!  S95_t1(x)%id=11004
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-004'
!  S95_t1(x)%xmin=120.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=10.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11004_CMS_H-ZZ-llll_1.13fb-1'

!  x=x+1
!  S95_t1(x)%id=11005
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-005'
!  S95_t1(x)%xmin=250.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=10.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11005_CMS_H-ZZ-llnunu_1.1fb-1'

!  x=x+1
!  S95_t1(x)%id=11006
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-006'
!  S95_t1(x)%xmin=230.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11006_CMS_H-ZZ-llqq_1fb-1'

!  x=x+1
!  S95_t1(x)%id=11009
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-009'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=140.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11009_CMS_SM_H-tautau_1.1fb-1'

!  x=x+1
!  S95_t1(x)%id=11010
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-010'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=140.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11010_CMS_H-gammagamma_1.09fb-1'

  !x=x+1
  !S95_t1(x)%id=2011111
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-111'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=240.0D0
  !S95_t1(x)%sep=5.0D0 
  !filename(x)='2011111_Atlas_H-WW_1.04fb-1_interpol'

  !x=x+1
  !S95_t1(x)%id=2011131
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-131'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=600.0D0
  !S95_t1(x)%sep=2.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='2011131_Atlas_H-ZZ_1.1fb-1_interpol'

!  x=x+1
!  S95_t1(x)%id=11012
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-012'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=135.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11012_CMS_H-bb_1.1fb-1'

!  x=x+1
!  S95_t1(x)%id=11013
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-013'
!  S95_t1(x)%xmin=180.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=10.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11013_CMS_H-ZZ-lltautau_1.1fb-1'

!  x=x+1
!  S95_t1(x)%id=11014
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-014'
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11014_CMS_H-WW_1.55fb-1_interpol'

!  x=x+1
!  S95_t1(x)%id=11015
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-015'
!  S95_t1(x)%xmin=115.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11015_CMS_H-ZZ-4l_1.66fb-1_interpol'

!  x=x+1
!  S95_t1(x)%id=11016
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-016'
!  S95_t1(x)%xmin=250.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11016_CMS_H-ZZ-llnunu_1.6fb-1_interpol'

!  x=x+1
!  S95_t1(x)%id=11017
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-017'
!  S95_t1(x)%xmin=225.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11017_CMS_H-ZZ-llqq_1.6fb-1_interpol'

!  x=x+1
!  S95_t1(x)%id=11020
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-020'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=145.0D0
!  S95_t1(x)%sep=5.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11020_CMS_H-tautau_1.6fb-1_SM'

!  x=x+1
!  S95_t1(x)%id=11021
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-021'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11021_CMS_SM_H-gaga_1.66fb-1'

!n.b.: Need fermiophobic-ness check for this analysis
!  x=x+1
!  S95_t1(x)%id=110212
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-021'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='110212_CMS_Fermiophob_H-gaga_1.66fb-1'

!  x=x+1
!  S95_t1(x)%id=11022
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-022'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11022_CMS_SMcombined_1.1-1.7fb-1'

!  x=x+1 
!  S95_t1(x)%id=5895
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='[hep-ex] arXiv:1108.5895 (ATLAS)'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='5895_Atlas_H-gaga_1.08fb-1'

!  x=x+1 
!  S95_t1(x)%id=110201
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-020'
!  S95_t1(x)%xmin=90.0D0
!  S95_t1(x)%xmax=500.0D0
!  S95_t1(x)%sep=10.0D0 
!  filename(x)='11020_CMS_H-tautau_1.6fb-1_MSSM_interpol'

  !x=x+1 
  !S95_t1(x)%id=2011161
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-161'
  !S95_t1(x)%xmin=110.0D0
  !S95_t1(x)%xmax=150.0D0
  !S95_t1(x)%sep=1.0D0 
  !S95_t1(x)%deltax=0.0D0 
  !filename(x)='2011161_Atlas_H-gaga_4.9fb-1'  
  
!  x=x+1
!  S95_t1(x)%id=2011162
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='ATL'
!  S95_t1(x)%label='ATLAS-CONF-2011-162'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='2011162_Atlas_H-ZZ-4l_4.8fb-1'
  
!  x=x+1 
!  S95_t1(x)%id=11024
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-024'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=5D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11024_CMS_H-WW-lnulnu_4.6fb-1_interpol'  
  
!  x=x+1
!  S95_t1(x)%id=11025
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-025'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0 
!  filename(x)='11025_CMS_H-ZZ-4l_4.7fb-1'

!  x=x+1
!  S95_t1(x)%id=11026
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-026'
!  S95_t1(x)%xmin=250.0D0
!  S95_t1(x)%xmax=590.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11026_CMS_H-ZZ-llnunu_4.6fb-1'
  
!  x=x+1
!  S95_t1(x)%id=110271
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-027'
!  S95_t1(x)%xmin=130.0D0
!  S95_t1(x)%xmax=164.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='110271_CMS_H-ZZ-llqq_4.6fb-1'

!  x=x+1
!  S95_t1(x)%id=110272
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-027'
!  S95_t1(x)%xmin=200.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='110272_CMS_H-ZZ-llqq_4.6fb-1'  
  
!  x=x+1
!  S95_t1(x)%id=11028
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-028'
!  S95_t1(x)%xmin=190.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11028_CMS_H-ZZ-lltautau_4.7fb-1'
  
!  x=x+1
!  S95_t1(x)%id=11030
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-030'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=150.0D0
!  S95_t1(x)%sep=0.5D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11030_CMS_H-gaga_4.8fb-1'
   
!  x=x+1
!  S95_t1(x)%id=11032
!  S95_t1(x)%particle_x=Hneut
!  S95_t1(x)%expt='CMS'
!  S95_t1(x)%label='CMS-PAS-HIG-11-032'
!  S95_t1(x)%xmin=110.0D0
!  S95_t1(x)%xmax=600.0D0
!  S95_t1(x)%sep=1.0D0 
!  S95_t1(x)%deltax=0.0D0   
!  filename(x)='11032_CMS_SMcombined_4.7fb-1'

  !x=x+1
  !S95_t1(x)%id=2011150
  !S95_t1(x)%particle_x=Hneut
  !S95_t1(x)%expt='ATL'
  !S95_t1(x)%label='ATLAS-CONF-2011-150'
  !S95_t1(x)%xmin=200.0D0
  !S95_t1(x)%xmax=600.0D0
  !S95_t1(x)%sep=20.0D0 
  !S95_t1(x)%deltax=0.0D0   
  !filename(x)='2011150_Atlas_H-ZZ-llqq_2.05fb-1'
                
!   x=x+1 
!   S95_t1(x)%id=6224
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 6227, Br(h->tautau)>6%'
!   S95_t1(x)%xmin=90.0D0
!   S95_t1(x)%xmax=300.0D0
!   S95_t1(x)%sep=10.0D0 
!   filename(x)='D0_h-bb_h-tautau_comb_Br0.06_5.2-7.3fb_6224' 
! 
!   x=x+1 
!   S95_t1(x)%id=6225
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 6227, Br(h->tautau)>10%'
!   S95_t1(x)%xmin=90.0D0
!   S95_t1(x)%xmax=300.0D0
!   S95_t1(x)%sep=10.0D0 
!   filename(x)='D0_h-bb_h-tautau_comb_Br0.10_5.2-7.3fb_6225'
! 
!   x=x+1 
!   S95_t1(x)%id=6226
!   S95_t1(x)%particle_x=Hneut
!   S95_t1(x)%expt=' D0'
!   S95_t1(x)%label='D0 Note 6227, Br(h->tautau)>14%'
!   S95_t1(x)%xmin=90.0D0
!   S95_t1(x)%xmax=300.0D0
!   S95_t1(x)%sep=10.0D0 
!   filename(x)='D0_h-bb_h-tautau_comb_Br0.14_5.2-7.3fb_6226'

  ! checks we've filled the whole array
  if(x.ne.xend)then
   stop'error in initializetables1 (a)'
  endif 
 
  ! do loop to read in S95 tables 
  col=3
  do x=xbeg,xend
   S95_t1(x)%nx=nint((S95_t1(x)%xmax-S95_t1(x)%xmin)/S95_t1(x)%sep)+1
   allocate(S95_t1(x)%dat(S95_t1(x)%nx,col-1))
  enddo

  open(file_id_common2,file = trim(adjustl(pathname))//'Expt_tables/' // &
      &  'S95_t1.binary',form='unformatted')

  read(file_id_common2,iostat=ios)S95_t1(xbeg)%dat

  if(ios.eq.0)then

    do x=xbeg+1,xend
     read(file_id_common2)S95_t1(x)%dat
    enddo

  else           
    rewind(file_id_common2)
    do x=xbeg,xend
     fullfilename=trim(adjustl(pathname))//'Expt_tables/' &
              &   //trim(adjustl(S95_t1(x)%expt))//'tables/' &
            &   //trim(filename(x))//'.txt'
      

     call read_tabletype1(S95_t1(x),5,col,fullfilename) 
#ifndef WEBVERSION     
     write(file_id_common2)S95_t1(x)%dat
#endif     
    enddo
  endif   
       
  close(file_id_common2)

  deallocate(filename)
             
 end subroutine initializetables1
 
 !************************************************************            
 subroutine read_tabletype1(t1,skip,col,fullfilename)
 !************************************************************       
 !fills t1%dat
  !--------------------------------------input
  type(table1) :: t1  
  integer :: skip,col
  character(LEN=*) :: fullfilename
  !-----------------------------------internal
  integer :: i,n      
  double precision :: xdummy,xdummy_store
  !-------------------------------------------

  t1%dat=0.0D0

  open(file_id_1, file=(trim(fullfilename)))
   
  do i=1,skip
   read(file_id_1,*) !skip lines
  enddo 

  xdummy_store = t1%xmin-t1%sep
  do i=1,t1%nx
   read(file_id_1,*)xdummy,(t1%dat(i,n),n=1,col-1) 

   ! checks that x are evenly spaced as expected
   if((abs(xdummy-xdummy_store-t1%sep).gt.1.0D-7) &
     &  .or.(abs(xdummy-(t1%xmin+dble(i-1)*t1%sep)).gt.1.0D-7))then
    write(*,*)i,t1%id,xdummy,t1%xmin+dble(i-1)*t1%sep
    stop 'error in read_tabletype1 (a1)'
   endif           

   xdummy_store=xdummy       
       
  enddo  
   
  if(abs(xdummy-t1%xmax).gt.1.0D-7)stop 'error in read_tabletype1 (a2)'

  close(file_id_1) 
  
 end subroutine read_tabletype1
 !************************************************************
end module S95tables_type1
!************************************************************
