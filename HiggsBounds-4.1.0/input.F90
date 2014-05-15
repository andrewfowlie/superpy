! This file is part of HiggsBounds
!  -KW
!******************************************************************
module input
!******************************************************************      
#ifdef NAGf90Fortran
 use F90_UNIX_ENV, only : iargc,getarg
 use F90_UNIX_IO, only : flush
#endif

 use usefulbits, only : Hneut,Hplus,Chineut,Chiplus
 implicit none

 integer,parameter :: nHmax=9 
  !nHneut>9 not allowed in this version
  !(if you really need nHneut>9,
  !change nHmax in subroutine setup_input,          
  !and adapt subroutine outputproc_t1/t2
  !and do_output, by changing the size of
  !the character variables nHchar,i,j
  !since LEN=1 will no longer be sufficient)
 integer :: f_orig=15 ! file id
                      ! n.b. debug_predratio.txt is 14, debug_channels.txt is 12 
                      ! (see HiggsBounds.f90 and HiggsBounds_subroutines.F90) 
 character(LEN=50),allocatable :: stem_array(:)
 logical,allocatable :: required(:)
 logical,allocatable :: isoptional(:) 
 integer :: nargs_datfile
 logical,parameter :: official=.True.!do not change this without contacting us first!
 !logical,parameter :: official=.False. 

 contains


 !************************************************************      
 subroutine setup_input
 !************************************************************
 ! * if inputmethod='datfile' or 'website', finds whichanalyses,whichinput,np
 ! * if inputmethod='datfile', finds infile1 (prefix for input/output filenames)
 ! * sets ndat (number of parameter points considered)
 ! * sets n_additional (number of additional data values for each parameter point)
 ! * allocates theo, g2, partR
 !************************************************************
  use usefulbits, only : np,ndat,n_additional,theo,g2,partR, &
         &               inputmethod,whichinput, &
         &               debug,infile1,infile2, &
         &               allocate_hadroncolliderextras_parts,allocate_dataset_parts, &
         &               allocate_sqcouplratio_parts,fill_pdesc

  use theory_BRfunctions      
  implicit none
  !-----------------------------------internal      
  integer :: f,ios,cc,g
  character(LEN=50) :: stem
  !-------------------------------------------      

  call fill_pdesc

  select case(inputmethod)
  case('datfile')
   call getbasiccommandline!get whichanalyses,whichinput,numbers of particles
   call getshortcommandline(infile1,infile2)! get infile1 
  case('website')
   call getbasiccommandline!get whichanalyses,whichinput,numbers of particles
   !don't need infile1
  case('subrout') 
   !np(Hneut),np(Hplus)  are already set
   call check_number_of_particles
   call check_whichanalyses
   !haven't yet set whichanalyses,whichinput
   infile1=''
  case default
   stop'incorrect value for inputmethod'
  end select

  !getting ndat and n_additional...
  select case(inputmethod)
  case('datfile')   
   f=f_orig
   n_additional=0
   select case(whichinput)
   case('part','effC','hadr') 
 
    call fill_stem_array

    g=1
    stem=stem_array(g)
    f=f_orig+g

    open(f,file=trim(infile1)//trim(stem)//'.dat',status='old',action='read',iostat=ios)
    if(ios.ne.0)then
     call file_name_msg(f)
     write(*,*)'Check that <prefix> was specified correctly'
     write(*,*)'and that file exists.'
     call flush(6)
     stop'Problem opening file: see standard output for more info'
    endif

    ndat=getfilelength(f) ! number of data sets (i.e. lines) in input 
 
    if((ndat.le.0))stop 'error getting ndat'         
    close(f)

    g=ubound(stem_array,dim=1)
    stem=stem_array(g)
    f=f_orig+g
    !the last element in stem_array should be 'additional'
    if(trim(stem_array(g)).ne.'additional')stop'Error in subroutine setup_input (a)'

    required(g)=.False.

    open(f,file=trim(infile1)//trim(stem)//'.dat',status='old',action='read',iostat=ios)  

    if(ios.eq.0)then    

     cc=count_columns(f)

     if(cc.gt.1)then
      if(getfilelength(f).eq.ndat)then
       n_additional=cc-1  
       required(g)=.True.
      endif
     endif
    endif
    close(f)

   case('SLHA')
    ndat=1  ;  n_additional=0 
   case default
    stop'error in subroutine do_input (1a)'
   end select
  case('website','subrout')
   ndat=1   ;  n_additional=0
  case default
   stop'error in subroutine do_input (1b)'
  end select
  !...finished getting ndat and n_additional

  if(debug)then
   write(*,*)'np(Hneut)=',np(Hneut)
   write(*,*)'ndat=',ndat
   write(*,*)'n_additional=',n_additional
  endif

  allocate(theo(ndat))

  call allocate_dataset_parts(theo,n_additional)                
  
  allocate(g2(ndat)) 
  call allocate_sqcouplratio_parts(g2)
  
  allocate(partR(ndat)) 
  call allocate_hadroncolliderextras_parts(partR)  

 end subroutine setup_input

!************************************************************      
 subroutine do_input
 ! only used for inputmethod='datfile' or 'website'
 ! determines what input is needed and calls the appropriate
 ! subroutines to get this input
 !************************************************************
  use usefulbits, only : ndat,inputmethod,whichinput, &
         &               theo,g2,partR,infile1
  use extra_bits_for_web, only : getlongcommandline2web  
  use extra_bits_for_SLHA, only : getSLHAdata
  use theory_BRfunctions      
  implicit none
  !-----------------------------------internal      
  integer :: n
  logical :: webdebugmode
  !-------------------------------------------      

  select case(inputmethod)
  case('datfile')   
   select case(whichinput)
   case('part','effC','hadr')

    do n=1,ubound(stem_array,dim=1)
      call readthefile(n) 
    enddo
    deallocate(stem_array)
    deallocate(required)
    deallocate(isoptional)

   case('SLHA') 
     n=1      
     if(ndat.ne.1)then
       stop'error in subroutine do_input (4): need to specify infile1 for each SLHA file somehow'
     endif          
     call getSLHAdata(theo(n),g2(n),infile1)
     !call test_input(n)
   case default
    stop'error in subroutine do_input (1)'
   end select

  case('website')   
   n=1
   call getlongcommandline2web(theo(n),g2(n),partR(n),webdebugmode)

   if(webdebugmode) call test_input(n)

  case default
   stop'error in subroutine do_input (3)'
  end select    

 end subroutine do_input
 !************************************************************
 subroutine fill_stem_array
 ! lists all possible files 
 ! note that stem_array(ubound(stem_array,dim=1)) should be 'additional'
 !************************************************************
 use usefulbits, only : np,whichinput,whichanalyses
 implicit none
 integer :: n,nt
 
   n=1
   nt=33

   allocate(stem_array(nt))
   allocate(required(nt))
   allocate(isoptional(nt))   

   ! the following table sets input filenames and says which input files are relevant to each option
   !                                      |               |                         |        np
   !                                      |   whichinput  |       whichanalyses     |Hneu Hcha Chineut Chiplus
   !         stem                         |part hadr effC |LandH onlyL onlyH  onlyP | ==0  ==0  ==0     ==0 
   call fill('MH_GammaTot'                ,   1,   1,   1,     1,    1,    1,    1,    0, 1,     1,    1) 
   call fill('MHall_uncertainties'        ,  -1,  -1,  -1,     1,    1,    1,    1,    1, 1,     1,    1)    
   call fill('MHplus_GammaTot'            ,   1,   1,   1,     1,    1,    1,    1,    1, 0,     1,    1) 
   call fill('MC_GammaTot'                ,   1,   1,   1,     1,    1,    1,    1,    1, 1,     0,    0) 
   call fill('MN_GammaTot'                ,   1,   1,   1,     1,    1,    1,    1,    1, 1,     0,    1) 
   call fill('effC'                       ,   0,   0,   1,     1,    1,    1,    1,    0, 1,     1,    1)
 
   call fill('BR_H_OP'                    ,   1,   1,   0,     1,    1,    1,    1,    0, 1,     1,    1) 
   call fill('BR_H_NP'                    ,   1,   1,   1,     1,    1,    1,    1,    0, 1,     1,    1) 
   call fill('BR_t'                       ,   1,   1,   1,     1,    0,    1,    1,    1, 0,     1,    1) 
   call fill('BR_Hplus'                   ,   1,   1,   1,     1,    1,    1,    1,    1, 0,     1,    1) 
   call fill('BR_C'                       ,   1,   1,   1,     1,    1,    1,    1,    1, 1,     0,    0) 
   call fill('BR_N'                       ,   1,   1,   1,     1,    1,    1,    1,    1, 1,     0,    1)

   call fill('LEP_HZ_CS_ratios'           ,   1,   1,   0,     1,    1,    0,    1,    0, 1,     1,    1) 
   call fill('LEP_H_ff_CS_ratios'         ,   1,   1,   0,     1,    1,    0,    1,    0, 1,     1,    1) 
   call fill('LEP_2H_CS_ratios'           ,   1,   1,   0,     1,    1,    0,    1,    0, 1,     1,    1) 
   call fill('LEP_HpHm_CS_ratios'         ,   1,   1,   1,     1,    1,    0,    1,    1, 0,     1,    1) 
   call fill('LEP_CpCm_CS'                ,   1,   1,   1,     1,    1,    0,    1,    1, 1,     0,    0) 
   call fill('LEP_2N_CS'                  ,   1,   1,   1,     1,    1,    0,    1,    1, 1,     0,    1) 

   call fill('TEVLHC_H_0jet_partCS_ratios',   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('TEVLHC_H_1jet_partCS_ratios',   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('TEVLHC_HW_partCS_ratios'    ,   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('TEVLHC_HZ_partCS_ratios'    ,   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 

   call fill('TEV_H_vbf_hadCS_ratios'     ,   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('TEV_H_tt_hadCS_ratios'      ,   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('TEV_1H_hadCS_ratios'        ,   0,   1,   0,     1,    0,    1,    1,    0, 1,     1,    1) 

   call fill('LHC7_H_vbf_hadCS_ratios'    ,   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('LHC7_H_tt_hadCS_ratios'     ,   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('LHC7_1H_hadCS_ratios'       ,   0,   1,   0,     1,    0,    1,    1,    0, 1,     1,    1)

   call fill('LHC8_H_vbf_hadCS_ratios'    ,   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('LHC8_H_tt_hadCS_ratios'     ,   1,   0,   0,     1,    0,    1,    1,    0, 1,     1,    1) 
   call fill('LHC8_1H_hadCS_ratios'       ,   0,   1,   0,     1,    0,    1,    1,    0, 1,     1,    1)

   call fill('CP_values'                  ,   1,   1,   0,     1,    1,    1,    1,    0, 1,     1,    1) 
   call fill('additional'                 ,   1,   1,   1,     1,    1,    1,    1,    1, 1,     1,    1) 

   if(n.ne.nt+1)stop'Error in subroutine fill_stem_array A'
 
!write(*,*)'hello whichinput',whichinput
!write(*,*)'hello whichanalyses',whichanalyses
!write(*,*)'hello np(Hneut)',np(Hneut)
!write(*,*)'hello np(Hplus)',np(Hplus)
!do n=1,nt
! write(*,*)'hello ',stem_array(n),required(n)
!enddo 
   !stop'hello ending here for now'


   contains
   !                                   |               |                        |        np
   !                                   |   whichinput  |       whichanalyses    |Hneu Hcha Chineut Chiplus
   !         stem                      |part hadr effC |LandH onlyL onlyH onlyP | ==0  ==0  ==0     ==0 

   subroutine fill(stem ,part,hadr,effC, LandH,onlyL,onlyH,onlyP, Hneu,Hcha, Chneu,  Chcha)
   !nb required(i) for 'additional' is not set until later
    implicit none

    character(LEN=*), intent(in):: stem
    integer, intent(in) :: part,hadr,effC,LandH,onlyL,onlyH,onlyP,Hneu,Hcha,Chneu,Chcha
    integer :: req

    stem_array(n)=stem 

    req=1
    select case(whichinput)
    case('part')
     req= part * req
    case('hadr')
     req= hadr * req
    case('effC')
     req= effC * req
    case default
     stop'error in subroutine fill(whichinput)'
    end select

    select case(whichanalyses)
    case('LandH')
     req= LandH * req
    case('onlyL')
     req= onlyL * req
    case('onlyH')
     req= onlyH * req
    case('onlyP')
     req= onlyP * req
    case default
     stop'error in subroutine fill(whichanalyses)'
    end select

    if(np(Hneut)==0)  req= Hneu  * req
    if(np(Hplus)==0)  req= Hcha  * req
    if(np(Chineut)==0)req= Chneu * req
    if(np(Chiplus)==0)req= Chcha * req

    select case(req)
    case(0)
     required(n)=.False.
     isoptional(n)=.False.
    case(1)
     required(n)=.True.
     isoptional(n)=.False.
    case(-1)
     required(n)=.False.
     isoptional(n)=.True. 
    case default
     stop'error in subroutine fill(req)' 
    end select
   
    n=n+1

   end subroutine fill

 end subroutine fill_stem_array

 !************************************************************      
 subroutine readthefile(n_stem)
 !************************************************************
 ! Opens and reads the input data from files
 ! checks the number of lines and columns in the file first
 ! First column is always line number. Last line is checked to make
 ! sure line number is correct
 !************************************************************ 
  use usefulbits, only : np,theo,ndat,infile1,g2,partR, &
                   &     n_additional
  implicit none  
  !--------------------------------------input
  integer,intent(in) :: n_stem
  !-----------------------------------internal
  character(LEN=50) :: stem
  integer :: jj,i,j,f,k,q,ios,x,y
  double precision :: nc
  logical :: needed, opt
  double precision, allocatable :: BR_hjhihi_in(:),BR_NjqqNi_in(:),BR_NjZNi_in(:)
  character(LEN=500) :: line
  integer :: numskip
  double precision, allocatable :: numbersfromline(:)
  logical :: readHneut, readHplus
  !-------------------------------------------      
  stem=stem_array(n_stem)  
  needed=required(n_stem)
  opt=isoptional(n_stem)
  
  if(needed.or.opt)then
   f=f_orig+n_stem
   open(f,file=trim(infile1)//trim(stem)//'.dat',status='old',action='read',iostat=ios)      

    if(ios.ne.0)then
     if(needed) then
     call file_name_msg(f)
     stop'problem opening file: see standard output for more info'
     else 
      write(*,*) 'Optional datafile '//trim(infile1)//trim(stem)//'.dat'//' not found.'//&
      &          ' Using default values.'
      return
     endif 
    endif

   if((ndat.ne.getfilelength(f)))then
     write(*,*)'wrong no. lines in input file'
     write(*,*)'It had',getfilelength(f),'lines'
     write(*,*)'but should have been',ndat
     call file_name_msg(f)
     stop 'error in input file length (see standard output for filename).'
   endif

!--NEW TO ALLOW DIFFERENT FILE READINGS OF DMHALL_UNCERTAINTIES.dat
   numskip=0
   readHneut=.True.
   readHplus=.True.
   if(isoptional(n_stem)) then
    if(get_ncol(stem).ne.count_columns(f))then
     if(get_ncol(stem).eq.np(Hneut)+1) then
      readHneut=.True.
      readHplus=.False.
      numskip=0
     else if(get_ncol(stem).eq.np(Hplus)+1) then
      readHneut=.False.
      readHplus=.True.
      numskip=count_columns(f)-1-np(Hplus)      
     else     
      readHneut=.False.
      readHplus=.False.
      numskip=0
	 endif
	endif 
   else
    if(get_ncol(stem).ne.count_columns(f))then
      write(*,*)'wrong no. columns in input file'
      write(*,*)'It had',count_columns(f),'columns'
      write(*,*)'but should have been',get_ncol(stem)
      write(*,*)'including line ID number'
      call file_name_msg(f)
      stop 'error in input file format (see standard output for filename).'
    endif
   endif 

   select case(trim(stem))
   case('MH_GammaTot') 
     x=Hneut
     do jj=1,ndat     
       read(f,*)  nc,  (theo(jj)%particle(x)%M(i)       ,i=1,np(x)), &
         &           (theo(jj)%particle(x)%GammaTot(i),i=1,np(x))
     enddo
   case('MHall_uncertainties') 
     x=Hneut
     y=Hplus
     do jj=1,ndat
      read(f,'(A)') line
!!      write(*,*) "Line: ", line
      call extractnumbersfromline(line, nc, numbersfromline)
!------Copy the HB mass uncertainties to the HS mass uncertainties for datfiles input.
      if(readHneut) then
       do i=1,np(x)
!        theo(jj)%particle(x)%dM(i) = theo(jj)%particle(x)%dMh(i)
        theo(jj)%particle(x)%dM(i)  = numbersfromline(i)
        theo(jj)%particle(x)%dMh(i) = numbersfromline(i)
       enddo 
       numskip=np(x)
      endif
      if(readHplus) then 
       do i=1,np(y)
!        theo(jj)%particle(y)%dM(i) = theo(jj)%particle(y)%dMh(i)
        theo(jj)%particle(y)%dM(i)  = numbersfromline(numskip+i)
        theo(jj)%particle(y)%dMh(i) = numbersfromline(numskip+i)
       enddo        
      endif 
      deallocate(numbersfromline)
!!      write(*,*) nc, theo(jj)%particle(Hneut)%dM, theo(jj)%particle(Hplus)%dM
     enddo
   case('MHplus_GammaTot')
     x=Hplus
     do jj=1,ndat     
       read(f,*)  nc,  (theo(jj)%particle(x)%M(i)       ,i=1,np(x)), &
         &           (theo(jj)%particle(x)%GammaTot(i),i=1,np(x))
     enddo
   case('MC_GammaTot')
     x=Chiplus
     do jj=1,ndat     
       read(f,*)  nc,  (theo(jj)%particle(x)%M(i)       ,i=1,np(x)), &
         &           (theo(jj)%particle(x)%GammaTot(i),i=1,np(x))
     enddo
   case('MN_GammaTot')
     x=Chineut
     do jj=1,ndat     
       read(f,*)  nc,  (theo(jj)%particle(x)%M(i)       ,i=1,np(x)), &
         &           (theo(jj)%particle(x)%GammaTot(i),i=1,np(x))
     enddo
   case('effC')        
     do jj=1,ndat                                    
      read(f,*)   nc,  (g2(jj)%hjss_s(i)                   ,i=1,np(Hneut)), &
           &           (g2(jj)%hjss_p(i)                   ,i=1,np(Hneut)), &
           &           (g2(jj)%hjcc_s(i)                   ,i=1,np(Hneut)), &
           &           (g2(jj)%hjcc_p(i)                   ,i=1,np(Hneut)), &
           &           (g2(jj)%hjbb_s(i)                   ,i=1,np(Hneut)), &
           &           (g2(jj)%hjbb_p(i)                   ,i=1,np(Hneut)), &
           &           (g2(jj)%hjtoptop_s(i)               ,i=1,np(Hneut)), &
           &           (g2(jj)%hjtoptop_p(i)               ,i=1,np(Hneut)), &
           &           (g2(jj)%hjmumu_s(i)                 ,i=1,np(Hneut)), &
           &           (g2(jj)%hjmumu_p(i)                 ,i=1,np(Hneut)), &
           &           (g2(jj)%hjtautau_s(i)               ,i=1,np(Hneut)), &
           &           (g2(jj)%hjtautau_p(i)               ,i=1,np(Hneut)), &
           &           (g2(jj)%hjWW(i)                     ,i=1,np(Hneut)), &
           &           (g2(jj)%hjZZ(i)                     ,i=1,np(Hneut)), &
           &           (g2(jj)%hjZga(i)                    ,i=1,np(Hneut)), &
           &           (g2(jj)%hjgaga(i)                   ,i=1,np(Hneut)), &
           &           (g2(jj)%hjgg(i)                     ,i=1,np(Hneut)), &
           &           (g2(jj)%hjggZ(i)                    ,i=1,np(Hneut)), &
           &          ((g2(jj)%hjhiZ(j,i)                  ,i=1,j),j=1,np(Hneut))   
    enddo     
   case('LEP_HZ_CS_ratios')                 
     do jj=1,ndat    
      read(f,*)  nc, (theo(jj)%lep%XS_hjZ_ratio(i)        ,i=1,np(Hneut)) 
     enddo      
   case('LEP_H_ff_CS_ratios')                 
     do jj=1,ndat    
      read(f,*)  nc, (theo(jj)%lep%XS_bbhj_ratio(i)       ,i=1,np(Hneut)), &
           &         (theo(jj)%lep%XS_tautauhj_ratio(i)   ,i=1,np(Hneut))
     enddo
   case('LEP_2H_CS_ratios')
     do jj=1,ndat                                   
      read(f,*)  nc, ((theo(jj)%lep%XS_hjhi_ratio(j,i),i=1,j),j=1,np(Hneut))   
     enddo 
   case('LEP_HpHm_CS_ratios')
     do jj=1,ndat                                   
      read(f,*)  nc, (theo(jj)%lep%XS_HpjHmj_ratio(i)     ,i=1,np(Hplus))   
     enddo 
   case('LEP_CpCm_CS')
     do jj=1,ndat                                   
      read(f,*)  nc, (theo(jj)%lep%XS_CpjCmj(i)     ,i=1,np(Chiplus))   
     enddo
   case('LEP_2N_CS')
     do jj=1,ndat                                   
      read(f,*)  nc, ((theo(jj)%lep%XS_NjNi(j,i),i=1,j),j=1,np(Chineut))   
     enddo
   case('TEVLHC_H_0jet_partCS_ratios')
    do jj=1,ndat     
     read(f,*)  nc, (partR(jj)%gg_hj(i)           ,i=1,np(Hneut)), &
         &          (partR(jj)%qq_hj(5,i)         ,i=1,np(Hneut))    
    enddo   
   case('TEVLHC_H_1jet_partCS_ratios')           
    do jj=1,ndat     
     read(f,*)  nc, (partR(jj)%bg_hjb(i),i=1,np(Hneut))   
    enddo    
   case('TEVLHC_HW_partCS_ratios')
    do jj=1,ndat     
     read(f,*)  nc, ((partR(jj)%qq_hjWp(q,i)       ,i=1,np(Hneut)), q=1,partR(jj)%nq_hjWp), &
         &          ((partR(jj)%qq_hjWm(q,i)       ,i=1,np(Hneut)), q=1,partR(jj)%nq_hjWm)  
    enddo      
   case('TEVLHC_HZ_partCS_ratios')
    do jj=1,ndat     
     read(f,*)  nc,(partR(jj)%gg_hjZ(i)          ,i=1,np(Hneut)), &
         &       (( partR(jj)%qq_hjZ(q,i)        ,i=1,np(Hneut)), q=1,partR(jj)%nq_hjZ ) 
    enddo   
   case('TEV_H_vbf_hadCS_ratios')    
    do jj=1,ndat     
     read(f,*)  nc, (theo(jj)%tev%XS_vbf_ratio(i),i=1,np(Hneut)) 
    enddo    
   case('TEV_H_tt_hadCS_ratios')    
    do jj=1,ndat     
     read(f,*)  nc, (theo(jj)%tev%XS_tthj_ratio(i),i=1,np(Hneut)) 
    enddo  
   case('TEV_1H_hadCS_ratios')  
     do jj=1,ndat    
      read(f,*)  nc, (theo(jj)%tev%XS_hj_ratio(i)   ,i=1,np(Hneut)) , &
                  &  (theo(jj)%tev%XS_hjb_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%tev%XS_hjW_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%tev%XS_hjZ_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%tev%XS_vbf_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%tev%XS_tthj_ratio(i) ,i=1,np(Hneut))   
     enddo     
   case('LHC7_H_vbf_hadCS_ratios')    
    do jj=1,ndat     
     read(f,*)  nc, (theo(jj)%lhc7%XS_vbf_ratio(i),i=1,np(Hneut)) 
    enddo    
   case('LHC7_H_tt_hadCS_ratios')    
    do jj=1,ndat     
     read(f,*)  nc, (theo(jj)%lhc7%XS_tthj_ratio(i),i=1,np(Hneut)) 
    enddo  
   case('LHC7_1H_hadCS_ratios')  
     do jj=1,ndat    
      read(f,*)  nc, (theo(jj)%lhc7%XS_hj_ratio(i)   ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc7%XS_hjb_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc7%XS_hjW_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc7%XS_hjZ_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc7%XS_vbf_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc7%XS_tthj_ratio(i) ,i=1,np(Hneut))   
     enddo            
   case('LHC8_H_vbf_hadCS_ratios')    
    do jj=1,ndat     
     read(f,*)  nc, (theo(jj)%lhc8%XS_vbf_ratio(i),i=1,np(Hneut)) 
    enddo    
   case('LHC8_H_tt_hadCS_ratios')    
    do jj=1,ndat     
     read(f,*)  nc, (theo(jj)%lhc8%XS_tthj_ratio(i),i=1,np(Hneut)) 
    enddo  
   case('LHC8_1H_hadCS_ratios')  
     do jj=1,ndat    
      read(f,*)  nc, (theo(jj)%lhc8%XS_hj_ratio(i)   ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc8%XS_hjb_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc8%XS_hjW_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc8%XS_hjZ_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc8%XS_vbf_ratio(i)  ,i=1,np(Hneut)) , &
                  &  (theo(jj)%lhc8%XS_tthj_ratio(i) ,i=1,np(Hneut))   
     enddo            
   case('BR_H_OP')       
     do jj=1,ndat          
       read(f,*)    nc, (theo(jj)%BR_hjss(i)       ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjcc(i)         ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjbb(i)         ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjmumu(i)       ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjtautau(i)     ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjWW(i)         ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjZZ(i)         ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjZga(i)        ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjgaga(i)       ,i=1,np(Hneut)), &
           &          (theo(jj)%BR_hjgg(i)         ,i=1,np(Hneut))
     enddo  
   case('BR_H_NP')
     if(np(Hneut)>1)then !because if np(Hneut)=1, then matrix BR_hjhihi has no off diagonal entries
       allocate(BR_hjhihi_in(np(Hneut)**2-np(Hneut)))   
       BR_hjhihi_in  =0.0D0  

       do jj=1,ndat      
        read(f,*)  nc,  (theo(jj)%BR_hjinvisible(i)    ,i=1,np(Hneut)), &
             &        (         BR_hjhihi_in(i)      ,i=1,int(np(Hneut)**2-np(Hneut)))  

        k=0
        do j=1,np(Hneut)
          do i=1,np(Hneut)
           if(i.ne.j)then
          k=k+1
          theo(jj)%BR_hjhihi(j,i) = BR_hjhihi_in(k)         
           endif           
          enddo
        enddo    
      
       enddo
       deallocate(BR_hjhihi_in)

     else
       do jj=1,ndat      
        read(f,*)  nc,  (theo(jj)%BR_hjinvisible(i)    ,i=1,np(Hneut))
       enddo
     endif
  case('BR_t')
    do jj=1,ndat     
     read(f,*)  nc,  theo(jj)%BR_tWpb,          &
         &          (theo(jj)%BR_tHpjb(i)       ,i=1,np(Hplus))  
    enddo 
  case('BR_Hplus')
    do jj=1,ndat     
     read(f,*)  nc, (theo(jj)%BR_Hpjcs(i)       ,i=1,np(Hplus)), &
         &          (theo(jj)%BR_Hpjcb(i)       ,i=1,np(Hplus)), &
         &          (theo(jj)%BR_Hpjtaunu(i)    ,i=1,np(Hplus))  
    enddo
  case('BR_C')
    do jj=1,ndat     
     read(f,*)  nc, ((theo(jj)%BR_CjqqNi(j,i)       ,i=1,np(Chineut)),j=1,np(Chiplus)), &
         &          ((theo(jj)%BR_CjlnuNi(j,i)      ,i=1,np(Chineut)),j=1,np(Chiplus)), &
         &          ((theo(jj)%BR_CjWNi(j,i)        ,i=1,np(Chineut)),j=1,np(Chiplus))  
    enddo
  case('BR_N')

     if(np(Chineut)>1)then !because if np(Chineut)=1, then BR matrix has no off diagonal entries
       allocate( BR_NjqqNi_in( np(Chineut)**2-np(Chineut) ) )
       allocate(  BR_NjZNi_in( np(Chineut)**2-np(Chineut) ) )
   
       BR_NjqqNi_in  =0.0D0
       BR_NjZNi_in    =0.0D0  

       do jj=1,ndat      
        read(f,*)  nc,  (BR_NjqqNi_in(i)    ,i=1,int(np(Chineut)**2-np(Chineut)))  , &
             &        (BR_NjZNi_in(i)      ,i=1,int(np(Chineut)**2-np(Chineut)))  

        k=0
        do j=1,np(Chineut)
          do i=1,np(Chineut)
           if(i.ne.j)then
          k=k+1
          theo(jj)%BR_NjqqNi(j,i) = BR_NjqqNi_in(k) 
          theo(jj)%BR_NjZNi(j,i)   = BR_NjZNi_in(k)        
           endif           
          enddo
        enddo    
      
       enddo

       deallocate(BR_NjqqNi_in)
       deallocate(BR_NjZNi_in)
     endif
  case('CP_values')
    do jj=1,ndat     
     read(f,*)  nc, (theo(jj)%CP_value(i)    ,i=1,np(Hneut)) 
    enddo
  case('additional')          
    do jj=1,ndat 
     read(f,*)  nc, (theo(jj)%additional(i),i=1,n_additional)                 
    enddo 
  case default
    stop'problem in subroutine readthefile (2)'
  end select 
  
  close(f)

  if(ndat.ne.nint(nc))then 
     write(*,*)'last line read in was not labled ndat'
     call file_name_msg(f)
     stop 'error in input file (see standard output for filename and more details).'    
  endif 
  
  endif
 end subroutine readthefile
 !************************************************************      
 subroutine extractnumbersfromline(line, nc, numbers)
 ! Scans a line and extracts numbers separated by whitespaces
 !************************************************************
  implicit none
  character(LEN=500), intent(in) :: line
  double precision, allocatable, intent(out) :: numbers(:)
  double precision, intent(out) :: nc
!  integer, intent(out) :: ii

  integer :: i, indx, prev, beginning, N
  double precision :: dbltmp
  
  prev      = 0 
  beginning = 1 
  N         = 0
 
  do i=1,len(line)
   indx = index('0123456789.EeDd-+', line(i:i))
   if (indx.eq.0 .and. prev.gt.0) then
    read(line(beginning:i-1), *) dbltmp
    N=N+1
   else if (indx.gt.0 .and. prev.eq.0) then
    beginning = i 
   end if
   prev = indx
  end do

  allocate(numbers(N-1))
  N=0
  do i=1,len(line)
   indx = index('0123456789.EeDd-+', line(i:i))
   if (indx.eq.0 .and. prev.gt.0) then
    N=N+1
    if(N.eq.1) then 
     read(line(beginning:i-1), *) nc
    else
     read(line(beginning:i-1), *) numbers(N-1)
    endif    
   else if (indx.gt.0 .and. prev.eq.0) then
    beginning = i 
   end if
   prev = indx
  end do
 
 end subroutine extractnumbersfromline
 !************************************************************      
 subroutine getbasiccommandline
 !************************************************************
 ! finds whichanalyses,whichinput,np from command line 
 !************************************************************ 
  !nb iargc and getarg are non-standard
  use usefulbits, only : np,pdesc,whichanalyses,whichinput,inputmethod      
  implicit none  
  !-----------------------------------internal
#ifndef NAGf90Fortran
  integer :: iargc
#endif
  character(LEN=100) :: temp
  character(LEN=5) :: nHtemp
  integer :: i,x,xmax
  integer :: number_args 
  logical :: wrong_args 
  logical :: additionalSLHAoutput = .False.
  !-------------------------------------------      
  number_args = IARGC() 

  if(inputmethod .eq. 'datfile')then

   if(official)then !whichanalyses, whichinput, np(Hneut), np(Hplus) , prefix
    nargs_datfile =  2+         2        +1
    np(Chineut)=0
    np(Chiplus)=0
   else             !whichanalyses, whichinput, all elements of np, prefix
    nargs_datfile =  2+ ubound(np,dim=1) +1
   endif

   wrong_args= (number_args .ne. nargs_datfile)
   additionalSLHAoutput = (number_args.eq.(nargs_datfile+1))
  elseif(inputmethod .eq. 'website')then 
   wrong_args= (number_args .lt. (2+ ubound(np,dim=1)))
  else
   write(*,*)'inputmethod=',inputmethod
   stop 'error in getbasiccommandline'
  endif

  if(wrong_args.and.(.not.additionalSLHAoutput))then
   write(*,*) "Incorrect number of parameters given on command line"
   call command_line_how2
   stop "Error: command line entered incorrectly (see standard output for more info)"
  endif
      
  ! Read arguments into text strings.
  i=1
  
  temp=""
  call GETARG(i,temp)
  whichanalyses = ""
  whichanalyses = trim(temp)
  i=i+1

  temp=""
  call GETARG(i,temp)
  whichinput = ""
  whichinput = trim(temp)
  i=i+1

  if(whichinput.ne.'SLHA'.and.wrong_args.and.additionalSLHAoutput)then
   write(*,*) "Incorrect number of parameters given on command line"
   call command_line_how2
   stop "Error: command line entered incorrectly (see standard output for more info)"
  endif
  
  if((inputmethod .eq. 'website').or.(.not.official))then
    xmax=ubound(np,dim=1)
  else !datfile, official
    xmax=2 !nHneut,nHplus
  endif

  do x=1,xmax
   temp=""
   call GETARG(i,temp)
   nHtemp = ""
   nHtemp = trim(temp)
   i=i+1

   if(verify(nHtemp," 1234567890").gt.0)then ! checks that the string nHtemp just contains the characters " 1234567890"
    ! the function verify is standard in fortran 95
    write(*,*)'Incorrect n'//trim(adjustl(pdesc(x)%short))//': not a number.'
    write(*,*)'(you entered n'//trim(adjustl(pdesc(x)%short))//'="'//trim(adjustl(nHtemp))//'")'
    write(*,*)'n'//trim(adjustl(pdesc(x)%short))//' is the number of '//trim(adjustl(pdesc(x)%long))//'s'
    call command_line_how2
    stop "Error: command line entered incorrectly (see standard output for more info)"
   endif

   read(nHtemp,*) np(x)
  enddo

  call check_number_of_particles
  call check_whichinput 
  call check_whichanalyses

 end subroutine getbasiccommandline
 !************************************************************      
 subroutine command_line_how2
 !************************************************************      
  use usefulbits, only : np,inputmethod

   write(*,*)'The correct syntax for the command line is:' 

   if(inputmethod.eq.'website')then
    if(ubound(np,dim=1).eq.4)then
     write(*,*)' ./HiggsBounds whichanalyses whichinput nHneut nHplus nNeutralino nChargino debug ...'
     write(*,*)'e.g.'
     write(*,*)' ./HiggsBounds LandH part 3 1 4 2 T ...'
    else
     stop'error in subroutine command_line_how2'
    endif
   elseif(.not.official)then
    if(ubound(np,dim=1).eq.4)then
     write(*,*)' ./HiggsBounds whichanalyses whichinput nHneut nHplus nNeutralino nChargino prefix'
     write(*,*)'e.g.'
     write(*,*)' ./HiggsBounds LandH part 3 1 4 2 mhmax'
    else
     stop'error in subroutine command_line_how2'
    endif
   else !official, datfile
    write(*,*)' ./HiggsBounds whichanalyses whichinput nHneut nHplus prefix'
    write(*,*)'e.g.'
    write(*,*)' ./HiggsBounds LandH part 3 1 mhmax'
   endif

   write(*,*)'See HiggsBounds manual for more details.' 
   call flush(6)
 end subroutine command_line_how2
 !************************************************************      
 subroutine check_whichinput
 !************************************************************      
  use usefulbits, only : whichinput

  select case(whichinput)
  case('part','effC','hadr','SLHA')
  case default
   call command_line_how2
   write(*,*)'Error: This value of "whichinput" is not allowed.'
   write(*,*)'(you entered whichinput="'//trim(adjustl(whichinput))//'")'
   write(*,*)'Allowed values for "whichinput" are (see manual):'
   write(*,*)'part        masses, branching ratios, decay widths, LEP cross sections,'
   write(*,*)'            Tevatron and LHC partonic cross sections'
   write(*,*)'effC        effective coupling approx' 
   write(*,*)'hadr        masses, branching ratios, decay widths, LEP cross sections,'
   write(*,*)'            Tevatron and LHC hadronic cross sections'
   write(*,*)'SLHA        SUSY Les Houches Accord files'     
   call flush(6)                  
   stop 'error: input type selected incorrectly (see standard output for more info)'   
  end select 

 end subroutine check_whichinput 

 !************************************************************      
 subroutine check_number_of_particles
 !************************************************************      
  use usefulbits, only : np,pdesc,inputmethod,whichinput
  integer :: x

  do x=1,ubound(np,dim=1)
   if(np(x).lt.0)then     
    write(*,*) 'number of '//trim(adjustl(pdesc(x)%long))//'s must be greater than zero'         
    if(inputmethod.eq.'datfile')call command_line_how2        
    stop 'error in subroutine check_number_of_particles (a) (see standard output for more info)'
   elseif(np(x)>nHmax)then
    write(*,*) 'number of '//trim(adjustl(pdesc(x)%long))//'s must be less than',nHmax
    write(*,*) ' (if you need more than this, please contact us)'         
    if(inputmethod.eq.'datfile')call command_line_how2        
    stop 'error in subroutine check_number_of_particles (b) (see standard output for more info)'
   endif
  enddo

  if(sum(np).eq.0)then
    stop'There should be a non-zero number of particles'
  endif

  if((inputmethod.eq.'datfile').and.(official))then
    do x=3,ubound(np,dim=1) !i.e. particles other than Hneut,Hplus
     if(np(x).gt.0)then
      write(*,*)'In this version of HiggsBounds, the number of ' &
         &//trim(adjustl(pdesc(x)%long))//'s must be zero'   
      write(*,*)'and you have entered the number',np(x)  
      write(*,*)'Please contact us if you would like more information.'
      stop'error in subroutine check_number_of_particles (c) (see standard output for more info)'
     endif
    enddo
  endif

  if(whichinput.eq.'SLHA')then
    if((np(Hneut).lt.0).or.(np(Hneut).gt.5))then
      write(*,*)'If a SLHA file is used as input,'
      write(*,*)'number of neutral Higgs must be in the range 0:5'
      stop'error in subroutine check_number_of_particles (d) (see standard output for more info)'
    endif
    if((np(Hplus).lt.0).or.(np(Hplus).gt.1))then
      write(*,*)'If a SLHA file is used as input,'
      write(*,*)'number of charged Higgs must be in the range 0:1'
      stop'error in subroutine check_number_of_particles (e) (see standard output for more info)'
    endif
  endif

 end subroutine check_number_of_particles 
 !************************************************************
 subroutine check_whichanalyses
 !************************************************************      
  use usefulbits, only : whichanalyses,inputmethod
 
  select case(whichanalyses)
  case('onlyL','onlyH','LandH','onlyP','list')
  case default
   if(inputmethod.eq.'datfile')call command_line_how2
   write(*,*)'Error: This value of "whichanalyses" is not allowed' 
   write(*,*)'(you entered whichanalyses="'//trim(adjustl(whichanalyses))//'")'
   write(*,*)'Allowed values for "whichanalyses" are:'
   write(*,*)'onlyL         only LEP results used'
   write(*,*)'onlyH         only hadronic colliders i.e. only Tevatron and LHC results used'   
   write(*,*)'LandH         LEP, Tevatron and LHC results used'
   write(*,*)'onlyP         use all analyses with an arXiv number'    
   call flush(6)                  
   stop 'error: experiment selected incorrectly (see standard output for more info)'  
  end select

 end subroutine check_whichanalyses 

 !************************************************************      
 subroutine getshortcommandline(inf1,inf2)
 !************************************************************
 ! used if inputmethod='datfile'
 ! finds infile1 from command line 
 !************************************************************ 
  !nb iargc and getarg are non-standard      
  implicit none
  !--------------------------------------input
  character(LEN=*) :: inf1
  character(LEN=*),optional :: inf2   
  !-----------------------------------internal
#ifndef NAGf90Fortran
  integer :: iargc
#endif
  character(LEN=100) :: temp
  integer :: i
  integer :: number_args
  !-------------------------------------------      
   
  number_args = IARGC() 
  
  if(number_args.ne.nargs_datfile)then
   if(number_args.eq.(nargs_datfile+1)) then
    i=nargs_datfile
    temp=""
    call GETARG(i,temp)
    inf1 = ""
    inf1 = trim(temp)
    i=i+1
    temp=""
    call GETARG(i,temp)
    if(present(inf2)) then
     inf2 = ""
     inf2 = trim(temp)  
    endif 
   else
    stop "Incorrect number of parameters given (getshortcommandline)"
   endif 
  else
   ! Read last argument into text string.
   i=nargs_datfile            
   temp=""
   call GETARG(i,temp)
   inf1 = ""
   inf1 = trim(temp)
   i=i+1
   
   if(present(inf2)) then   
    inf2=inf1
   endif 
  endif
  
  write(*,*)"------"
  write(*,*)"From command line:"
  write(*,*)"prefix:",trim(inf1) 
  write(*,*)"------"
  call flush(6)

 end subroutine getshortcommandline            
            
 !****************************************************
 function getfilelength(fileid)
 !****************************************************
 ! calculates file length and checks for errors
 ! nb. files must end in a 'newline' character
 !**************************************************** 
  implicit none
  !--------------------------------------input 
  integer fileid   
  !-----------------------------------function 
  integer :: getfilelength
  !-----------------------------------internal      
  integer :: n,ios,m
  character(LEN=5) :: filechar
  character(LEN=20) :: sample
  !-------------------------------------------      
  
  write(filechar,'(I5)')fileid      

  !this will count the number of lines in the file, including the last one,
  ! even if it doesn't end in a newline character
  n = 0                        
  do
   read(fileid,'(a)',iostat=ios) sample

   if(ios.lt.0)then
    exit
   elseif(ios.gt.0) then
    call file_name_msg(fileid)
    stop 'error in input file: see standard output'
   elseif(trim(adjustl(sample)).eq.'')then
    write(*,*)'No blank lines allowed in input file.' 
    call file_name_msg(fileid)
    stop 'no blank lines allowed in input file: see standard output'
   endif
   n = n + 1               
  enddo            

  if(n.eq.0)stop'File is empty'
      
  getfilelength=n
  rewind(fileid)

  m = 0 ;  
  do     !this will count the number of lines which end in a newline character      
   read(fileid,*,iostat=ios)
   if(ios.lt.0) exit 
   m = m + 1
  enddo
  rewind(fileid)

  if(m.ne.n)then !checking that every line end in a newline character
    call file_name_msg(fileid)
    stop'Error: file needs to end with a newline character (see standard output for filename)'
  endif

 end function getfilelength

 !****************************************************
 subroutine file_name_msg(fileid)
 !****************************************************
 use usefulbits,only : infile1
 implicit none
 integer, intent(in) :: fileid
    write(*,*)'The problematic input file is called:' 
    write(*,*)'  '//trim(adjustl(infile1))//trim(adjustl(stem_array(fileid-f_orig)))//'.dat'
    call flush(6)
 end subroutine file_name_msg

 !****************************************************
 function count_columns(fileid)
 !****************************************************
 ! calculateshow many columns of numbers there are in file 
 ! with the id 'fileid'
 ! max line length= 1000 characters
 ! assumes objects in columns are separated by spaces 
 ! and that first line has the same
 ! number of columns as the rest of the file
 !****************************************************
  implicit none
  !--------------------------------------input
  integer, intent(in) :: fileid
  !-----------------------------------function
  integer :: count_columns
  !-----------------------------------internal
  integer :: n,x,a
  character(LEN=5000) :: teststr
  !-------------------------------------------

  n=0
  read(fileid,'(a)')teststr ! reads first line into a string 'teststr'

  teststr=adjustl(teststr) ! shifts characters left so that first char is not a space
  
  a=1
  do while (a.gt.0)  !replace all tab characters with spaces

   a=INDEX(teststr,'	')
   teststr=teststr(:a-1)//'   '//teststr(a+1:)

  enddo

  teststr=adjustl(teststr) ! shifts characters left so that first char is not a space

  do while (trim(adjustl(teststr)).ne.'') !while teststr is not made of spaces   
   n=n+1   
   x=INDEX(teststr,' ') !finds position of first space from left
   teststr=adjustl(teststr(x:)) !chops of everything left of the first space, then adjustl
  enddo  
  
  count_columns=n 

  rewind(fileid) 
 
 end function count_columns 
 !********************************************************
 subroutine test_input(n)
 ! prints out input (useful for testing)
 !********************************************************
  use usefulbits, only : np, theo, g2, partR
  implicit none
  !--------------------------------------input
  integer :: n
  !-----------------------------------internal
  integer :: j, i
  !-------------------------------------------

  write(*,*)'**************INPUT**************' 

#include "write_out_input.txt" 

  write(*,*)'*********************************' 
 end subroutine test_input   
 !*****************************************************
 function get_ncol(stem)
 ! calculates how many columns  HiggsBounds should expect 
 ! to find in file stem//'.dat'
  use usefulbits, only : np,n_additional

  implicit none
  !--------------------------------------input
   character(LEN=50), intent(in) :: stem
  !-----------------------------------function
   integer :: get_ncol
  !-----------------------------------internal
   integer :: x
  !-------------------------------------------

   select case(trim(stem))
   case('MH_GammaTot')
    x=2*np(Hneut)
   case('MHall_uncertainties')
    x=np(Hneut)+np(Hplus) 
   case('MHplus_GammaTot')
    x=2*np(Hplus)
   case('MC_GammaTot')
    x=2*np(Chiplus)
   case('MN_GammaTot')
    x=2*np(Chineut)
   case('BR_H_NP')
    x=1*np(Hneut)+np(Hneut)*(np(Hneut)-1)
   case('BR_H_OP')
    x=10*np(Hneut)
   case('LEP_HZ_CS_ratios') 
    x=1*np(Hneut)
   case('LEP_H_ff_CS_ratios') 
    x=2*np(Hneut)
   case('LEP_2H_CS_ratios')
    x= (np(Hneut)*(np(Hneut)+1))/2        
   case('LEP_HpHm_CS_ratios') 
    x=1*np(Hplus)
   case('LEP_CpCm_CS') 
    x=1*np(Chiplus)
   case('LEP_2N_CS') 
    x=1*(np(Chineut)*(np(Chineut)+1))/2  
   case('effC')
    x=18*np(Hneut)  + (np(Hneut)*(np(Hneut)+1))/2
   case('TEVLHC_H_0jet_partCS_ratios')
    x=2*np(Hneut)
   case('TEVLHC_H_1jet_partCS_ratios')
    x=1*np(Hneut)
   case('TEVLHC_HW_partCS_ratios')
    x=4*np(Hneut)
   case('TEVLHC_HZ_partCS_ratios')
    x=6*np(Hneut)
   case('TEV_H_vbf_hadCS_ratios')
    x=1*np(Hneut)
   case('TEV_H_tt_hadCS_ratios')
    x=1*np(Hneut)  
   case('TEV_1H_hadCS_ratios')  
    x=6*np(Hneut)
   case('LHC7_H_vbf_hadCS_ratios')
    x=1*np(Hneut)
   case('LHC7_H_tt_hadCS_ratios')
    x=1*np(Hneut)  
   case('LHC7_1H_hadCS_ratios')  
    x=6*np(Hneut)
   case('LHC8_H_vbf_hadCS_ratios')
    x=1*np(Hneut)
   case('LHC8_H_tt_hadCS_ratios')
    x=1*np(Hneut)  
   case('LHC8_1H_hadCS_ratios')  
    x=6*np(Hneut)
   case('BR_t')
    x=1+np(Hplus)
   case('BR_Hplus')
    x=3*np(Hplus)
   case('BR_C')
    x=3*np(Chiplus)*np(Chineut)
   case('BR_N')
    x=2*np(Chineut)*(np(Chineut)-1)
   case('CP_values')
    x=1*np(Hneut)
   case('additional')   
    x=n_additional 
   case default 
    write(*,*)'stem=',stem
    call flush(6)
    stop 'error in function get_ncol'
   end select

   ! each file has the line number as the first column
   get_ncol= x + 1

 end function get_ncol
 !************************************************************

end module input
!****************************************************
