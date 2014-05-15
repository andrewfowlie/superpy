! This file is part of HiggsBounds
!  -KW
!******************************************************************
module usefulbits
!******************************************************************            
 implicit none 
 logical :: debug
 logical :: full_dmth_variation = .True.
 integer :: dmhsteps = 3
! Mass uncertainties smaller than 0.1 GeV are not considered
 double precision :: small_mh = 0.1D0
 logical :: run_HB_classic = .False.

! For the LEP chisq extension:
 logical :: chisqcut_at_mumax = .False.
 
 character(LEN=5) :: whichanalyses
 character(LEN=4) :: whichinput 
 character(LEN=7) :: inputmethod = 'subrout'
 character(LEN=9),parameter :: vers='4.1.0'
 integer :: n_additional                       
 character(len=100) :: infile1,infile2
 integer,parameter :: file_id_common=10 
 integer,parameter :: file_id_common2=12
 integer,parameter :: file_id_common3=133
 integer,parameter :: file_id_debug1=444
 integer,parameter :: file_id_debug2=45

 integer, allocatable :: analysislist(:)

 !read from http://pdg.lbl.gov/ 22.10.2009
 double precision,parameter :: ms=0.105D0
 double precision,parameter :: mc=1.27D0
 double precision,parameter :: mbmb=4.20D0
 double precision,parameter :: mmu=105.7D-3
 double precision,parameter :: mtau=1.777D0

 double precision,parameter :: MZ=91.1876D0 !PDG 2009
 double precision,parameter :: MW=80.398D0 !PDG 2009
 double precision,parameter :: GF=1.16637D-5
 double precision,parameter :: pi=3.14159265358979323846264338328D0
 double precision,parameter :: alphas=0.118D0

 double precision,parameter :: small=1.0D-6
 double precision,parameter :: vsmall=1.0D-16
 double precision,parameter :: vvsmall=1.0D-100

 type particledescriptions
  character(LEN=10) :: short
  character(LEN=30) :: long
 end type
 
! particle codes: (n.b. these are NOT pdg)
 integer,parameter :: not_a_particle   = 0
 integer,parameter :: Hneut   = 1 !either Mhi, Mh2 or Mh3 (says nothing about CP properties)
 integer,parameter :: Hplus   = 2 !single charged Higgs
 integer,parameter :: Chineut = 3 !either neutralino1, neutralino2, neutralino3 or neutralino4
 integer,parameter :: Chiplus = 4 !either chargino1 or chargino2
 integer :: np(0:4)=1 !e.g np(Hneut) holds number of neutral Higgs considered
 type(particledescriptions),allocatable :: pdesc(:)



 !for subroutine version--------------------
 type inputsubroutineinfo
  integer :: stat
  character(LEN=40) :: desc
  integer :: req
 end type
 type(inputsubroutineinfo),allocatable :: inputsub(:)

 logical :: just_after_run  

 !associated with 'channels'----------------
 integer :: ntot                        
            
 type listprocesses            
  integer :: tlist,ttype
  integer :: findi,findj
  integer :: corresponding_clsb_table_element             
 end type
            
 type(listprocesses), allocatable :: pr(:)
 type(listprocesses), allocatable :: prsep(:,:)
 !-------------------------------------------
 
 !associated with 'input'--------------------

 type particlemasses            
  double precision, allocatable  :: M(:)
! Central value for mass with uncertainties
  double precision, allocatable  :: Mc(:)
  double precision, allocatable  :: GammaTot(:)
! Mass uncertainties used for HiggsSignals
  double precision, allocatable  :: dM(:)
! Mass uncertainties (variation) used in HB
  double precision, allocatable  :: dMh(:)   
 end type

  double precision, allocatable :: diffMhneut(:,:)
  double precision, allocatable :: diffMhch(:,:)
  double precision, allocatable :: dmn(:)
  double precision, allocatable :: dmch(:)
  
  integer ndmh

 integer ndat

 type lepdataset
  double precision, allocatable :: XS_hjZ_ratio(:)
  double precision, allocatable :: XS_bbhj_ratio(:)
  double precision, allocatable :: XS_tautauhj_ratio(:)
  double precision, allocatable :: XS_hjhi_ratio(:,:)  
  double precision, allocatable :: XS_HpjHmj_ratio(:)  
  double precision, allocatable :: XS_CpjCmj(:)   
  double precision, allocatable :: XS_NjNi(:,:)   
 end type 

 type hadroncolliderdataset
  double precision, allocatable :: XS_hj_ratio(:) 
  double precision, allocatable :: XS_hjZ_ratio(:)
  double precision, allocatable :: XS_hjW_ratio(:) 
  double precision, allocatable :: XS_hjb_ratio(:) 
  double precision, allocatable :: XS_tthj_ratio(:) 
  double precision, allocatable :: XS_vbf_ratio(:)  
  double precision, allocatable :: XS_HZ_SM(:)
  double precision, allocatable :: XS_HW_SM(:)
  double precision, allocatable :: XS_H_SM(:)
  !double precision, allocatable :: XS_H_SM_9713(:),XS_H_SM_9674(:)
  double precision, allocatable :: XS_ttH_SM(:)
  double precision, allocatable :: XS_vbf_SM(:)
  ! Higgs produced in association with b, where b is tagged, comes uncut and with various cuts
  ! see subroutines in theory_XS_SM_functions.f90 for details
  double precision, allocatable :: XS_Hb_SM(:),XS_Hb_c1_SM(:),XS_Hb_c2_SM(:)
  double precision, allocatable :: XS_Hb_c3_SM(:),XS_Hb_c4_SM(:)
 end type

 type dataset        
  logical :: gooddataset
  integer, allocatable :: CP_value(:)
  double precision, allocatable :: additional(:)
  type(particlemasses), allocatable :: particle(:)
  double precision, allocatable :: BR_hjss(:),BR_hjcc(:)
  double precision, allocatable :: BR_hjbb(:),BR_hjmumu(:),BR_hjtautau(:)   
  double precision, allocatable :: BR_hjinvisible(:)
  double precision, allocatable :: BR_hjhihi(:,:)  
  type(lepdataset) :: lep 
  !-------------------------------------------
  double precision, allocatable :: BR_hjWW(:),BR_hjgaga(:)
  double precision, allocatable :: BR_hjZga(:)
  double precision, allocatable :: BR_hjZZ(:),BR_hjgg(:) 

  double precision :: BR_tWpb
  double precision, allocatable :: BR_tHpjb(:)

  double precision, allocatable :: BR_Hpjcs(:)
  double precision, allocatable :: BR_Hpjcb(:)
  double precision, allocatable :: BR_Hpjtaunu(:)
 
  double precision, allocatable :: BR_CjqqNi(:,:)
  double precision, allocatable :: BR_CjlnuNi(:,:)
  double precision, allocatable :: BR_CjWNi(:,:)
  double precision, allocatable :: BR_NjqqNi(:,:)
  double precision, allocatable :: BR_NjZNi(:,:)
 
  type(hadroncolliderdataset) :: tev 
  type(hadroncolliderdataset) :: lhc7 
  type(hadroncolliderdataset) :: lhc8
  
  double precision, allocatable :: BR_Hbb_SM(:),BR_Hcc_SM(:), BR_Hss_SM(:),BR_Hmumu_SM(:),BR_Htautau_SM(:)  
  double precision, allocatable :: BR_HWW_SM(:),BR_HZZ_SM(:),BR_HZga_SM(:),BR_Hgaga_SM(:),BR_Hgg_SM(:)
  double precision, allocatable :: BR_Hjets_SM(:)
  double precision, allocatable :: GammaTot_SM(:) 
  !------------------------------------------- 
 end type            
            
 type(dataset), allocatable :: theo(:) 

 type sqcouplratio 
  double precision, allocatable :: hjss_s(:),hjss_p(:)
  double precision, allocatable :: hjcc_s(:),hjcc_p(:)
  double precision, allocatable :: hjbb_s(:),hjbb_p(:)
  double precision, allocatable :: hjtoptop_s(:),hjtoptop_p(:)
  double precision, allocatable :: hjmumu_s(:),hjmumu_p(:)
  double precision, allocatable :: hjtautau_s(:),hjtautau_p(:)

  double precision, allocatable :: hjWW(:),hjZZ(:)
  double precision, allocatable :: hjZga(:)
  double precision, allocatable :: hjgaga(:),hjgg(:),hjggZ(:)
  double precision, allocatable :: hjhiZ(:,:)  
 end type
 
 type(sqcouplratio), allocatable :: g2(:)  
 
 type hadroncolliderextras  
  !nq_hjWp,nq_hjWm,nq_hj,nq_hjZ are set in allocate_hadroncolliderextras_parts below 
  double precision, allocatable :: qq_hjWp(:,:)
  integer :: nq_hjWp!=2 i.e. (u dbar), (c sbar)  e.g. allocate(tR%qq_hjWp(tR%nq_hjWp,np(Hneut))) 
  double precision, allocatable :: qq_hjWm(:,:)
  integer :: nq_hjWm!=2 i.e. (ubar d), (cbar s) 
   
  double precision, allocatable :: gg_hj(:)
  double precision, allocatable :: qq_hj(:,:)
  integer :: nq_hj!=5 i.e.(d dbar), (u ubar), (s sbar), (c cbar), (b bbar) 
     
  double precision, allocatable :: gg_hjZ(:)
  double precision, allocatable :: qq_hjZ(:,:)
  integer :: nq_hjZ!=5 i.e.(d dbar), (u ubar), (s sbar), (c cbar), (b bbar)
     
  double precision, allocatable :: bg_hjb(:)
 end type
 
 type(hadroncolliderextras), allocatable :: partR(:)
 !-------------------------------------------  
  
 !associated with 'output'--------------------
 integer rep
 
 type results                  
  integer, allocatable :: chan(:)            
  double precision, allocatable :: obsratio(:)
  double precision, allocatable :: sfactor(:)
  double precision, allocatable :: axis_i(:)
  double precision, allocatable :: axis_j(:)
  integer, allocatable :: allowed95(:)
  integer, allocatable :: ncombined(:) 
  character(LEN=4), allocatable :: channelselection(:)                                                  
 end type            
                        
 type(results), allocatable :: res(:)

 !--new in HB-4:
 type fullresults
  integer :: chan = 0
  integer :: ncombined = 0
  integer :: allowed95 = 1
  double precision :: obsratio = 0.0D0
 end type
 
 type(fullresults), allocatable :: fullHBres(:)

 integer, allocatable :: allocate_if_stats_required(:)
 !-------------------------------------------
             
 contains

 subroutine HiggsBounds_info
  implicit none

  write(*,*)
  write(*,*)"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(*,*)"~                                                        ~"
  write(*,*)"~                   HiggsBounds "//adjustl(vers)//"                ~"
  write(*,*)"~                                                        ~"
  write(*,*)"~     Philip Bechtle, Oliver Brein, Sven Heinemeyer,     ~"
  write(*,*)"~       Oscar Stål, Tim Stefaniak, Georg Weiglein,       ~"
  write(*,*)"~                   Karina E. Williams                   ~"
  write(*,*)"~                                                        ~"
  write(*,*)"~            arXiv:0811.4169, arXiv:1102.1898,           ~"
  write(*,*)"~            arXiv:1301.2345, arXiv:1311.0055            ~"  
  write(*,*)"~            http://higgsbounds.hepforge.org             ~"
  write(*,*)"~                                                        ~"
  write(*,*)"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(*,*)
  write(*,*)"HiggsBounds collects together results from "
  write(*,*)
  write(*,*)"    * the LEP collaborations and LEP Higgs Working Group"
  write(*,*)"    * the CDF and D0 Collaborations"
  write(*,*)"    * the ATLAS and CMS Collaborations"
  write(*,*)"    * the program HDECAY (arXiv:hep-ph/9704448)"
  write(*,*)"    * TeV4LHC Higgs Working Group report"
  write(*,*)"      (see arXiv:hep-ph/0612172 and ref. therein)"
  write(*,*)"    * LHC Higgs Cross Section Working Group"
  write(*,*)"      (arXiv:1101.0593, arXiv:1201.3084 and ref. therein)"
  write(*,*)
  
  
 end subroutine HiggsBounds_info

 !**********************************************************            
 function div(a,b,divlimit,div0res)
 !**********************************************************
 ! be careful about using this - not a mathematical limit
  double precision :: div  
  !--------------------------------------input
  double precision :: a,b,divlimit,div0res
  !-----------------------------------internal
  double precision :: small1,small2
  !-------------------------------------------  
  small1 = 1.0D-28 
  small2 = 1.0D-20  
     
  if(abs(b).gt.small1)then
   div=a/b
  elseif(abs(a).lt.small2)then  
   div=divlimit  
   if(div.lt.0)stop 'error type divA (see function div in module usefulbits)'   
  else  
   div=div0res
   if(div.lt.0)stop 'error type divB (see function div in module usefulbits)'    
  endif
    
 end function

!--TESTING
 !**********************************************************            
 subroutine iselementofarray(value, array, output)
 !**********************************************************
   implicit none
  !-------------------------------------input and output   
   double precision, intent(in) :: value
   double precision, allocatable, dimension(:), intent(in) :: array
   integer, intent(out) :: output
  !---------------------------------------------internal      
   integer :: i
   double precision :: small
  !-----------------------------------------------------
   small = 1.0D-20
   output = -1

  if(allocated(array)) then
    do i=lbound(array,dim=1),ubound(array,dim=1)
      if(abs(value-array(i)).le.small) output = 1
    enddo
  else 
    stop'error: Passing an unallocated array to subroutine iselementofarray!'
  endif

 end subroutine iselementofarray
!----
 !**********************************************************            
 subroutine fill_pdesc
 !**********************************************************
   integer :: x

   if(ubound(np,dim=1).ne.4)stop'error: have made a mistake in subroutine fill_pdesc (1)'

   x=0
   allocate(  pdesc(  ubound(np,dim=1) )  )

   x=x+1
   pdesc(x)%short='h'
   pdesc(x)%long ='neutral Higgs boson'

   x=x+1
   pdesc(x)%short='hplus'
   pdesc(x)%long ='charged Higgs boson'

   x=x+1
   pdesc(x)%short='N'
   pdesc(x)%long ='neutralino'

   x=x+1
   pdesc(x)%short='C'
   pdesc(x)%long ='chargino'

   if(x.ne.ubound(np,dim=1))stop'error: have made a mistake in subroutine fill_pdesc (2)'

 end subroutine fill_pdesc
 !**********************************************************            
 subroutine allocate_dataset_parts(d,n_addit)
 !**********************************************************
  implicit none
  !------------------------------------------- 
  type(dataset) :: d(:)
  !--------------------------------------input 
  integer, intent(in) :: n_addit
  !-----------------------------------internal
  integer :: n_add,x,y
  integer, allocatable :: np_t(:) 
  !-------------------------------------------      
  allocate(np_t(lbound(np,dim=1):ubound(np,dim=1)))

  np_t=np
  do x=lbound(np_t,dim=1),ubound(np_t,dim=1)
   if(np(x)>0)then
    np_t(x)=np(x)
   elseif(np(x).eq.0)then
    np_t(x)=1
   else
    write(*,*)'np=',np
    stop'error in subroutine allocate_dataset_parts (1)'
   endif
  enddo

  if(n_addit>0)then
    n_add=n_addit
  elseif(n_addit.eq.0)then
    n_add=1
  else
    stop'error in subroutine allocate_dataset_parts (2)'
  endif

  do x=lbound(d,dim=1),ubound(d,dim=1)
  
   allocate(d(x)%additional(n_add))   
   
   allocate(d(x)%particle( ubound(np_t,dim=1)  ))
   do y= 1,ubound(np_t,dim=1)
    allocate(d(x)%particle(y)%M(         np_t(y)        ))
    allocate(d(x)%particle(y)%Mc(         np_t(y)        ))
    allocate(d(x)%particle(y)%GammaTot(  np_t(y)        ))   
    allocate(d(x)%particle(y)%dM(        np_t(y)        ))
    allocate(d(x)%particle(y)%dMh(        np_t(y)        ))
   enddo

   allocate(d(x)%lep%XS_hjZ_ratio(       np_t(Hneut)    ))
   allocate(d(x)%lep%XS_bbhj_ratio(      np_t(Hneut)    ))
   allocate(d(x)%lep%XS_tautauhj_ratio(  np_t(Hneut)    ))
   allocate(d(x)%lep%XS_hjhi_ratio(      np_t(Hneut),np_t(Hneut)    ))
   allocate(d(x)%lep%XS_HpjHmj_ratio(    np_t(Hplus)    ))
   allocate(d(x)%lep%XS_CpjCmj(          np_t(Chiplus)  ))
   allocate(d(x)%lep%XS_NjNi(            np_t(Chineut),np_t(Chineut)  ))

   allocate(d(x)%BR_hjss(                np_t(Hneut)    ))
   allocate(d(x)%BR_hjcc(                np_t(Hneut)    ))
   allocate(d(x)%BR_hjbb(                np_t(Hneut)    ))        
   allocate(d(x)%BR_hjmumu(              np_t(Hneut)    ))      
   allocate(d(x)%BR_hjtautau(            np_t(Hneut)    ))            
 
   allocate(d(x)%BR_hjhihi(              np_t(Hneut),np_t(Hneut)    ))
   allocate(d(x)%BR_hjWW(                np_t(Hneut)    ))
   allocate(d(x)%BR_hjZZ(                np_t(Hneut)    ))
   allocate(d(x)%BR_hjZga(               np_t(Hneut)    ))
   allocate(d(x)%BR_hjgaga(              np_t(Hneut)    ))
   allocate(d(x)%BR_hjgg(                np_t(Hneut)    ))
   allocate(d(x)%BR_hjinvisible(         np_t(Hneut)    ))

   allocate(d(x)%BR_tHpjb(               np_t(Hplus)    ))
   allocate(d(x)%BR_Hpjcs(               np_t(Hplus)    ))
   allocate(d(x)%BR_Hpjcb(               np_t(Hplus)    ))
   allocate(d(x)%BR_Hpjtaunu(            np_t(Hplus)    ))

   allocate(d(x)%BR_CjqqNi(              np_t(Chiplus),np_t(Chineut)    ))
   allocate(d(x)%BR_CjlnuNi(             np_t(Chiplus),np_t(Chineut)    ))
   allocate(d(x)%BR_CjWNi(               np_t(Chiplus),np_t(Chineut)    ))
   allocate(d(x)%BR_NjqqNi(              np_t(Chineut),np_t(Chineut)    ))
   allocate(d(x)%BR_NjZNi(               np_t(Chineut),np_t(Chineut)    ))  

   allocate(d(x)%tev%XS_hjb_ratio(       np_t(Hneut)    ))
   allocate(d(x)%tev%XS_tthj_ratio(      np_t(Hneut)    ))
   allocate(d(x)%tev%XS_vbf_ratio(       np_t(Hneut)    ))
    
   allocate(d(x)%tev%XS_hjZ_ratio(       np_t(Hneut)    ))   
   allocate(d(x)%tev%XS_hjW_ratio(       np_t(Hneut)    ))  
   allocate(d(x)%tev%XS_hj_ratio(        np_t(Hneut)    ))

   allocate(d(x)%lhc7%XS_hjb_ratio(       np_t(Hneut)    ))
   allocate(d(x)%lhc7%XS_tthj_ratio(      np_t(Hneut)    ))
   allocate(d(x)%lhc7%XS_vbf_ratio(       np_t(Hneut)    ))
    
   allocate(d(x)%lhc7%XS_hjZ_ratio(       np_t(Hneut)    ))   
   allocate(d(x)%lhc7%XS_hjW_ratio(       np_t(Hneut)    ))  
   allocate(d(x)%lhc7%XS_hj_ratio(        np_t(Hneut)    ))

   allocate(d(x)%lhc8%XS_hjb_ratio(       np_t(Hneut)    ))
   allocate(d(x)%lhc8%XS_tthj_ratio(      np_t(Hneut)    ))
   allocate(d(x)%lhc8%XS_vbf_ratio(       np_t(Hneut)    ))
    
   allocate(d(x)%lhc8%XS_hjZ_ratio(       np_t(Hneut)    ))   
   allocate(d(x)%lhc8%XS_hjW_ratio(       np_t(Hneut)    ))  
   allocate(d(x)%lhc8%XS_hj_ratio(        np_t(Hneut)    ))


   allocate(d(x)%CP_value(               np_t(Hneut)    ))

   do y= 1,ubound(np_t,dim=1)
    d(x)%particle(y)%M        =-1.0D0
    d(x)%particle(y)%Mc       =-1.0D0
    d(x)%particle(y)%GammaTot =0.0D0
    d(x)%particle(y)%dM       =0.0D0
    d(x)%particle(y)%dMh      =0.0D0
   enddo

   d(x)%lep%XS_hjZ_ratio      =0.0D0
   d(x)%lep%XS_bbhj_ratio     =0.0D0
   d(x)%lep%XS_tautauhj_ratio =0.0D0
   d(x)%lep%XS_hjhi_ratio     =0.0D0 
   d(x)%lep%XS_HpjHmj_ratio   =0.0D0 
   d(x)%lep%XS_CpjCmj         =0.0D0 
   d(x)%lep%XS_NjNi           =0.0D0 

   d(x)%BR_hjss               =0.0D0
   d(x)%BR_hjcc               =0.0D0 
   d(x)%BR_hjbb               =0.0D0
   d(x)%BR_hjmumu             =0.0D0
   d(x)%BR_hjtautau           =0.0D0
   d(x)%BR_hjhihi             =0.0D0
   d(x)%BR_hjWW               =0.0D0
   d(x)%BR_hjZZ               =0.0D0
   d(x)%BR_hjZga              =0.0D0
   d(x)%BR_hjgaga             =0.0D0
   d(x)%BR_hjgg               =0.0D0
   d(x)%BR_hjinvisible        =0.0D0

   d(x)%BR_tWpb               =0.0D0
   d(x)%BR_tHpjb              =0.0D0
   d(x)%BR_Hpjcs              =0.0D0
   d(x)%BR_Hpjcb              =0.0D0
   d(x)%BR_Hpjtaunu           =0.0D0

   d(x)%BR_CjqqNi             =0.0D0
   d(x)%BR_CjlnuNi            =0.0D0
   d(x)%BR_CjWNi              =0.0D0
   d(x)%BR_NjqqNi             =0.0D0
   d(x)%BR_NjZNi              =0.0D0

   d(x)%tev%XS_hjb_ratio      =0.0D0
   d(x)%tev%XS_tthj_ratio     =0.0D0
   d(x)%tev%XS_vbf_ratio      =0.0D0

   d(x)%tev%XS_hj_ratio       =0.0D0
   d(x)%tev%XS_hjW_ratio      =0.0D0
   d(x)%tev%XS_hjZ_ratio      =0.0D0

   d(x)%lhc7%XS_hjb_ratio      =0.0D0
   d(x)%lhc7%XS_tthj_ratio     =0.0D0
   d(x)%lhc7%XS_vbf_ratio      =0.0D0

   d(x)%lhc7%XS_hj_ratio       =0.0D0
   d(x)%lhc7%XS_hjW_ratio      =0.0D0
   d(x)%lhc7%XS_hjZ_ratio      =0.0D0

   d(x)%lhc8%XS_hjb_ratio      =0.0D0
   d(x)%lhc8%XS_tthj_ratio     =0.0D0
   d(x)%lhc8%XS_vbf_ratio      =0.0D0

   d(x)%lhc8%XS_hj_ratio       =0.0D0
   d(x)%lhc8%XS_hjW_ratio      =0.0D0
   d(x)%lhc8%XS_hjZ_ratio      =0.0D0


   d(x)%additional            =0.0D0

   d(x)%CP_value=0  
  enddo

  select case(whichanalyses)
  case('onlyH','LandH','onlyP','list ')
   do x=lbound(d,dim=1),ubound(d,dim=1)
    allocate(d(x)%tev%XS_HZ_SM(         np_t(Hneut)    ))
    allocate(d(x)%tev%XS_HW_SM(         np_t(Hneut)    ))
    allocate(d(x)%tev%XS_H_SM(          np_t(Hneut)    ))  
    allocate(d(x)%tev%XS_ttH_SM(        np_t(Hneut)    ))  
    allocate(d(x)%tev%XS_vbf_SM(        np_t(Hneut)    ))
     
    allocate(d(x)%tev%XS_Hb_SM(         np_t(Hneut)    )) 
    allocate(d(x)%tev%XS_Hb_c1_SM(      np_t(Hneut)    )) 
    allocate(d(x)%tev%XS_Hb_c2_SM(      np_t(Hneut)    )) 
    allocate(d(x)%tev%XS_Hb_c3_SM(      np_t(Hneut)    ))
    allocate(d(x)%tev%XS_Hb_c4_SM(      np_t(Hneut)    ))
 
    allocate(d(x)%lhc7%XS_HZ_SM(         np_t(Hneut)    ))
    allocate(d(x)%lhc7%XS_HW_SM(         np_t(Hneut)    ))
    allocate(d(x)%lhc7%XS_H_SM(          np_t(Hneut)    ))  
    allocate(d(x)%lhc7%XS_ttH_SM(        np_t(Hneut)    ))  
    allocate(d(x)%lhc7%XS_vbf_SM(        np_t(Hneut)    ))
     
    allocate(d(x)%lhc7%XS_Hb_SM(         np_t(Hneut)    )) 
    allocate(d(x)%lhc7%XS_Hb_c1_SM(      np_t(Hneut)    )) 
    allocate(d(x)%lhc7%XS_Hb_c2_SM(      np_t(Hneut)    )) 
    allocate(d(x)%lhc7%XS_Hb_c3_SM(      np_t(Hneut)    )) 

    allocate(d(x)%lhc8%XS_HZ_SM(         np_t(Hneut)    ))
    allocate(d(x)%lhc8%XS_HW_SM(         np_t(Hneut)    ))
    allocate(d(x)%lhc8%XS_H_SM(          np_t(Hneut)    ))  
    allocate(d(x)%lhc8%XS_ttH_SM(        np_t(Hneut)    ))  
    allocate(d(x)%lhc8%XS_vbf_SM(        np_t(Hneut)    ))
     
    allocate(d(x)%lhc8%XS_Hb_SM(         np_t(Hneut)    )) 
    allocate(d(x)%lhc8%XS_Hb_c1_SM(      np_t(Hneut)    )) 
    allocate(d(x)%lhc8%XS_Hb_c2_SM(      np_t(Hneut)    )) 
    allocate(d(x)%lhc8%XS_Hb_c3_SM(      np_t(Hneut)    )) 


    allocate(d(x)%BR_Hbb_SM(            np_t(Hneut)    ))   
    allocate(d(x)%BR_Hcc_SM(            np_t(Hneut)    ))   
    allocate(d(x)%BR_Hss_SM(            np_t(Hneut)    ))       
    allocate(d(x)%BR_Hmumu_SM(          np_t(Hneut)    ))
    allocate(d(x)%BR_Htautau_SM(        np_t(Hneut)    ))  
    allocate(d(x)%BR_HWW_SM(            np_t(Hneut)    ))  
    allocate(d(x)%BR_HZZ_SM(            np_t(Hneut)    ))  
    allocate(d(x)%BR_HZga_SM(           np_t(Hneut)    )) 
    allocate(d(x)%BR_Hgaga_SM(          np_t(Hneut)    ))  
    allocate(d(x)%BR_Hgg_SM(          np_t(Hneut)    ))  
    allocate(d(x)%BR_Hjets_SM(          np_t(Hneut)    )) 
    allocate(d(x)%GammaTot_SM(          np_t(Hneut)    ))
   enddo
  case('onlyL')
  case default
   stop'error in allocate_dataset_parts (3)'
  end select
  deallocate(np_t)
 end subroutine allocate_dataset_parts

 !**********************************************************            
 subroutine allocate_sqcouplratio_parts(gsq)
 ! to use this, gsq must be an array
 !********************************************************** 
  implicit none
  !------------------------------------------- 
  type(sqcouplratio) :: gsq(:)
  !-----------------------------------internal
  integer :: x
  integer :: nHiggsneut
  !-------------------------------------------       
  
  if(np(Hneut)>0)then
    nHiggsneut=np(Hneut)
  elseif(np(Hneut).eq.0)then
    nHiggsneut=1
  else
    stop'error in subroutine allocate_sqcouplratio_parts (1)'
  endif

  do x=lbound(gsq,dim=1),ubound(gsq,dim=1) 
   allocate(gsq(x)%hjss_s(nHiggsneut)    ,gsq(x)%hjss_p(nHiggsneut))
   allocate(gsq(x)%hjcc_s(nHiggsneut)    ,gsq(x)%hjcc_p(nHiggsneut))
   allocate(gsq(x)%hjbb_s(nHiggsneut)    ,gsq(x)%hjbb_p(nHiggsneut))
   allocate(gsq(x)%hjtoptop_s(nHiggsneut),gsq(x)%hjtoptop_p(nHiggsneut))
   allocate(gsq(x)%hjmumu_s(nHiggsneut)  ,gsq(x)%hjmumu_p(nHiggsneut))
   allocate(gsq(x)%hjtautau_s(nHiggsneut),gsq(x)%hjtautau_p(nHiggsneut))

   allocate(gsq(x)%hjWW(nHiggsneut)    ,gsq(x)%hjZZ(nHiggsneut)    )
   allocate(gsq(x)%hjZga(nHiggsneut)                               )
   allocate(gsq(x)%hjgaga(nHiggsneut)  ,gsq(x)%hjgg(nHiggsneut)    )
   allocate(gsq(x)%hjggZ(nHiggsneut)    )
   allocate(gsq(x)%hjhiZ(nHiggsneut,nHiggsneut)                 )  

   gsq(x)%hjss_s           =0.0D0
   gsq(x)%hjss_p           =0.0D0
   gsq(x)%hjcc_s           =0.0D0
   gsq(x)%hjcc_p           =0.0D0
   gsq(x)%hjbb_s           =0.0D0
   gsq(x)%hjbb_p           =0.0D0
   gsq(x)%hjtoptop_s       =0.0D0
   gsq(x)%hjtoptop_p       =0.0D0
   gsq(x)%hjmumu_s         =0.0D0
   gsq(x)%hjmumu_p         =0.0D0
   gsq(x)%hjtautau_s       =0.0D0
   gsq(x)%hjtautau_p       =0.0D0

   gsq(x)%hjWW           =0.0D0
   gsq(x)%hjZZ           =0.0D0
   gsq(x)%hjZga          =0.0D0
   gsq(x)%hjgaga         =0.0D0
   gsq(x)%hjgg           =0.0D0
   gsq(x)%hjggZ          =0.0D0
   gsq(x)%hjhiZ          =0.0D0
  enddo 
 
 end subroutine allocate_sqcouplratio_parts

!  !**********************************************************            
!  subroutine deallocate_sqcouplratio_parts(gsq)
!  !********************************************************** 
!   implicit none
!   !--------------------------------------input  
!   type(sqcouplratio) :: gsq(:)
!   !-----------------------------------internal
!   integer :: x
!   !-------------------------------------------       
!   
!   do x=lbound(gsq,dim=1),ubound(gsq,dim=1) 
!    deallocate(gsq(x)%hjbb     )
!    deallocate(gsq(x)%hjtautau )
!    deallocate(gsq(x)%hjWW     )
!    deallocate(gsq(x)%hjZZ     )
!    deallocate(gsq(x)%hjgaga   )
!    deallocate(gsq(x)%hjgg     )
!    deallocate(gsq(x)%hjggZ    )
!    deallocate(gsq(x)%hjhiZ    )  
!   enddo 
!  
!  end subroutine deallocate_sqcouplratio_parts
! 
 !**********************************************************            
 subroutine allocate_hadroncolliderextras_parts(tR)
 !********************************************************** 
  implicit none
  !-------------------------------------------  
  type(hadroncolliderextras) :: tR(:)
  !-----------------------------------internal
  integer :: x
  integer :: nHiggsneut
  !-------------------------------------------       
    
  if(np(Hneut)>0)then
    nHiggsneut=np(Hneut)
  elseif(np(Hneut).eq.0)then
    nHiggsneut=1
  else
    stop'error in subroutine allocate_hadroncolliderextras_parts (1)'
  endif

  tR%nq_hjWp=2 ! (u dbar), (c sbar)  e.g
  tR%nq_hjWm=2 ! (ubar d), (cbar s) 
  tR%nq_hj=5   !(d dbar), (u ubar), (s sbar), (c cbar), (b bbar) 
  tR%nq_hjZ=5  !(d dbar), (u ubar), (s sbar), (c cbar), (b bbar)
  
  do x=lbound(tR,dim=1),ubound(tR,dim=1) 
   allocate(tR(x)%qq_hjWp(tR(x)%nq_hjWp,nHiggsneut))
   allocate(tR(x)%qq_hjWm(tR(x)%nq_hjWm,nHiggsneut))
   allocate(tR(x)%gg_hj(nHiggsneut))
   allocate(tR(x)%qq_hj(tR(x)%nq_hj,nHiggsneut)) 
   allocate(tR(x)%gg_hjZ(nHiggsneut))
   allocate(tR(x)%qq_hjZ(tR(x)%nq_hjZ,nHiggsneut)) 
   allocate(tR(x)%bg_hjb(nHiggsneut)) 

   tR(x)%qq_hjWp   =0.0D0
   tR(x)%qq_hjWm   =0.0D0
   tR(x)%gg_hj     =0.0D0
   tR(x)%qq_hj     =0.0D0
   tR(x)%gg_hjZ    =0.0D0
   tR(x)%qq_hjZ    =0.0D0 
   tR(x)%bg_hjb    =0.0D0 
  enddo 
 
 end subroutine allocate_hadroncolliderextras_parts
 
 !**********************************************************            
 subroutine deallocate_hadroncolliderextras_parts(tR)
 !********************************************************** 
  implicit none
  !--------------------------------------input  
  type(hadroncolliderextras) :: tR(:)
  !-----------------------------------internal
  integer :: x
  !-------------------------------------------       
  
  do x=lbound(tR,dim=1),ubound(tR,dim=1) 
   deallocate(tR(x)%qq_hjWp)
   deallocate(tR(x)%qq_hjWm)
   deallocate(tR(x)%gg_hj)
   deallocate(tR(x)%qq_hj) 
   deallocate(tR(x)%gg_hjZ)
   deallocate(tR(x)%qq_hjZ)   
   deallocate(tR(x)%bg_hjb)  
  enddo 
 
 end subroutine deallocate_hadroncolliderextras_parts

 !**********************************************************            
 subroutine deallocate_usefulbits
 !**********************************************************
 ! deallocates theo,res (and everything inside)
 ! deallocates c,predratio,fact 
 !************************************************************
  implicit none
  !-----------------------------------internal
  integer x,y
  !-------------------------------------------
  deallocate(pdesc)!allocated in fill_pdesc

  !these are allocated in subroutine do_input
  do x=lbound(theo,dim=1),ubound(theo,dim=1)
   deallocate(theo(x)%additional)  

   do y= 1,ubound(np,dim=1)
    deallocate(theo(x)%particle(y)%M)
    deallocate(theo(x)%particle(y)%GammaTot)
    deallocate(theo(x)%particle(y)%dM)
    deallocate(theo(x)%particle(y)%dMh)
   enddo   
   deallocate(theo(x)%particle)
   
   deallocate(theo(x)%lep%XS_hjZ_ratio) 
   deallocate(theo(x)%lep%XS_bbhj_ratio) 
   deallocate(theo(x)%lep%XS_tautauhj_ratio) 
   deallocate(theo(x)%lep%XS_hjhi_ratio) 
   deallocate(theo(x)%lep%XS_HpjHmj_ratio) 
   deallocate(theo(x)%lep%XS_CpjCmj)
   deallocate(theo(x)%lep%XS_NjNi)
   
   deallocate(theo(x)%BR_hjss)
   deallocate(theo(x)%BR_hjcc)
   deallocate(theo(x)%BR_hjbb)
   deallocate(theo(x)%BR_hjmumu)
   deallocate(theo(x)%BR_hjtautau)
   deallocate(theo(x)%BR_hjhihi)
   deallocate(theo(x)%BR_hjWW)   
   deallocate(theo(x)%BR_hjZZ)   
   deallocate(theo(x)%BR_hjZga)  
   deallocate(theo(x)%BR_hjgaga)    
   deallocate(theo(x)%BR_hjgg)  
   deallocate(theo(x)%BR_hjinvisible)  

   deallocate(theo(x)%BR_tHpjb) 
   deallocate(theo(x)%BR_Hpjcs) 
   deallocate(theo(x)%BR_Hpjcb)
   deallocate(theo(x)%BR_Hpjtaunu)

   deallocate(theo(x)%BR_CjqqNi)
   deallocate(theo(x)%BR_CjlnuNi)
   deallocate(theo(x)%BR_CjWNi)
   deallocate(theo(x)%BR_NjqqNi)
   deallocate(theo(x)%BR_NjZNi)  
 
   deallocate(theo(x)%tev%XS_hjb_ratio)  
   deallocate(theo(x)%tev%XS_tthj_ratio)   
   deallocate(theo(x)%tev%XS_vbf_ratio)       
                 
   deallocate(theo(x)%tev%XS_hjZ_ratio)
   deallocate(theo(x)%tev%XS_hjW_ratio)
   deallocate(theo(x)%tev%XS_hj_ratio)      

   deallocate(theo(x)%lhc7%XS_hjb_ratio)  
   deallocate(theo(x)%lhc7%XS_tthj_ratio)   
   deallocate(theo(x)%lhc7%XS_vbf_ratio)       
                 
   deallocate(theo(x)%lhc7%XS_hjZ_ratio)
   deallocate(theo(x)%lhc7%XS_hjW_ratio)
   deallocate(theo(x)%lhc7%XS_hj_ratio) 

   deallocate(theo(x)%lhc8%XS_hjb_ratio)  
   deallocate(theo(x)%lhc8%XS_tthj_ratio)   
   deallocate(theo(x)%lhc8%XS_vbf_ratio)       
                 
   deallocate(theo(x)%lhc8%XS_hjZ_ratio)
   deallocate(theo(x)%lhc8%XS_hjW_ratio)
   deallocate(theo(x)%lhc8%XS_hj_ratio) 


   !deallocate(theo(x)%inLEPrange_Hpj) 
   !deallocate(theo(x)%inTEVrange_Hpj)

   deallocate(theo(x)%CP_value)                     
  enddo 
    
  select case(whichanalyses)
  case('onlyH','LandH','onlyP','list ')            
   do x=lbound(theo,dim=1),ubound(theo,dim=1)        
    deallocate(theo(x)%BR_Hbb_SM)       
    deallocate(theo(x)%BR_Hss_SM)           
    deallocate(theo(x)%BR_Hcc_SM)       
    deallocate(theo(x)%BR_Hmumu_SM)   
    deallocate(theo(x)%BR_Htautau_SM)     
    deallocate(theo(x)%BR_HWW_SM)    
    deallocate(theo(x)%BR_HZZ_SM)    
    deallocate(theo(x)%BR_HZga_SM)    
    deallocate(theo(x)%BR_Hgaga_SM)   
    deallocate(theo(x)%BR_Hgg_SM)       
    deallocate(theo(x)%BR_Hjets_SM)  
    deallocate(theo(x)%GammaTot_SM)
                    
    deallocate(theo(x)%tev%XS_HZ_SM)
    deallocate(theo(x)%tev%XS_HW_SM)
    deallocate(theo(x)%tev%XS_H_SM) 
    deallocate(theo(x)%tev%XS_ttH_SM)
    deallocate(theo(x)%tev%XS_vbf_SM) 
    !deallocate(theo(x)%tev%XS_H_SM_9713) 
    !deallocate(theo(x)%tev%XS_H_SM_9674) 
    
    deallocate(theo(x)%tev%XS_Hb_SM) 
    deallocate(theo(x)%tev%XS_Hb_c1_SM) 
    deallocate(theo(x)%tev%XS_Hb_c2_SM)  
    deallocate(theo(x)%tev%XS_Hb_c3_SM)    
    deallocate(theo(x)%tev%XS_Hb_c4_SM)

    deallocate(theo(x)%lhc7%XS_HZ_SM)
    deallocate(theo(x)%lhc7%XS_HW_SM)
    deallocate(theo(x)%lhc7%XS_H_SM) 
    deallocate(theo(x)%lhc7%XS_ttH_SM)
    deallocate(theo(x)%lhc7%XS_vbf_SM) 
    
    deallocate(theo(x)%lhc7%XS_Hb_SM) 
    deallocate(theo(x)%lhc7%XS_Hb_c1_SM) 
    deallocate(theo(x)%lhc7%XS_Hb_c2_SM)  
    deallocate(theo(x)%lhc7%XS_Hb_c3_SM)   

    deallocate(theo(x)%lhc8%XS_HZ_SM)
    deallocate(theo(x)%lhc8%XS_HW_SM)
    deallocate(theo(x)%lhc8%XS_H_SM) 
    deallocate(theo(x)%lhc8%XS_ttH_SM)
    deallocate(theo(x)%lhc8%XS_vbf_SM) 
    
    deallocate(theo(x)%lhc8%XS_Hb_SM) 
    deallocate(theo(x)%lhc8%XS_Hb_c1_SM) 
    deallocate(theo(x)%lhc8%XS_Hb_c2_SM)  
    deallocate(theo(x)%lhc8%XS_Hb_c3_SM)   

                
   enddo 
  case('onlyL')
  case default 
   stop'error in deallocate_usefulbits'
  end select 
  
  deallocate(theo) !allocated in subroutine do_input  

  !allocated in subroutine setup_output
  if(allocated(res)) then
   do x=lbound(res,dim=1),ubound(res,dim=1)  
    deallocate(res(x)%chan)   
    deallocate(res(x)%obsratio)
    deallocate(res(x)%axis_i)
    deallocate(res(x)%axis_j)
    deallocate(res(x)%sfactor) 
    deallocate(res(x)%allowed95) 
    deallocate(res(x)%ncombined) 
   enddo     
   deallocate(res) !allocated in subroutine setup_output                 
  endif

  if (allocated(fullHBres)) then
   deallocate(fullHBres)
  endif
    
!  call deallocate_sqcouplratio_parts(g2)
  do x=lbound(g2,dim=1),ubound(g2,dim=1) 

   deallocate(g2(x)%hjss_s)
   deallocate(g2(x)%hjss_p)
   deallocate(g2(x)%hjcc_s)
   deallocate(g2(x)%hjcc_p)
   deallocate(g2(x)%hjbb_s)
   deallocate(g2(x)%hjbb_p)
   deallocate(g2(x)%hjtoptop_s)
   deallocate(g2(x)%hjtoptop_p)
   deallocate(g2(x)%hjmumu_s)
   deallocate(g2(x)%hjmumu_p)
   deallocate(g2(x)%hjtautau_s)
   deallocate(g2(x)%hjtautau_p)

   deallocate(g2(x)%hjWW)  
   deallocate(g2(x)%hjZZ)
   deallocate(g2(x)%hjZga)
   deallocate(g2(x)%hjgaga)
   deallocate(g2(x)%hjgg)
   deallocate(g2(x)%hjggZ)
   deallocate(g2(x)%hjhiZ) 
  enddo 
  deallocate(g2) 
  
  !these are allocated in subroutine do_input 
  call deallocate_hadroncolliderextras_parts(partR)          
  deallocate(partR) !allocated in subroutine do_input      

  if(allocated(pr)) deallocate(pr)  !allocated in subroutine fill_pr or fill_pr_select
  if(allocated(prsep)) deallocate(prsep)  !allocated in subroutine fill_pr or fill_pr_select
  
  if(allocated(diffMhneut)) deallocate(diffMhneut)
  if(allocated(diffMhch)) deallocate(diffMhch)
  
  if(allocated(dmn)) deallocate(dmn)
  if(allocated(dmch)) deallocate(dmch)
  
 end subroutine deallocate_usefulbits
 !**********************************************************            

end module usefulbits                        
!******************************************************************
