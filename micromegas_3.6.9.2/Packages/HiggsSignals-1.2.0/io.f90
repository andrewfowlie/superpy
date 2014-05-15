!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 03/03/2013)
!--------------------------------------------------------------------
module io
!--------------------------------------------------------------------

 use usefulbits, only : Hneut,Hplus,Chineut,Chiplus
 use input, only : count_columns, getfilelength, stem_array, required, f_orig
 implicit none

 contains

!--------------------------------------------------------------------
 subroutine setup_input_for_hs
!--------------------------------------------------------------------
! * if inputmethod='datfile' finds np,Exptdata
! * if inputmethod='datfile', finds infile1 (prefix for input/output filenames)
! * sets ndat (number of parameter points considered)
! * sets n_additional (number of additional data values for each parameter point)
! * allocates theo, g2, partR
!--------------------------------------------------------------------
  use usefulbits, only : np,ndat,n_additional,theo,g2,partR, &
         &               inputmethod,whichinput, &
         &               debug,infile1,infile2, &
         &               allocate_hadroncolliderextras_parts,allocate_dataset_parts, &
         &               allocate_sqcouplratio_parts,fill_pdesc
  use input, only : getshortcommandline,check_number_of_particles, fill_stem_array,  &
         &          file_name_msg
  use theory_BRfunctions      
  implicit none
  !-----------------------------------internal      
  integer :: f,ios,cc,g
  character(LEN=50) :: stem
  !-------------------------------------------      

  call fill_pdesc

  select case(inputmethod)
  case('datfile')
   call getbasiccommandline_for_hs!get whichanalyses,whichinput,numbers of particles
   call getshortcommandline(infile1,infile2)! get infile1 
  case('subrout') 
   !np(Hneut),np(Hplus)  are already set
   call check_number_of_particles
!   call check_whichanalyses
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

 end subroutine setup_input_for_hs
!--------------------------------------------------------------------
 subroutine getbasiccommandline_for_hs
!--------------------------------------------------------------------
 ! finds whichanalyses,whichinput,np from command line 
!--------------------------------------------------------------------
  !nb iargc and getarg are non-standard
  use usefulbits, only : np,pdesc,whichinput,inputmethod
  use usefulbits_hs, only : Exptdir, runmode, pdf
  use input, only : check_number_of_particles, check_whichinput, nargs_datfile
  implicit none  
  !-----------------------------------internal
  character(LEN=100) :: temp
  character(LEN=5) :: nHtemp, pdftemp
  integer :: i,x
  integer :: number_args 
  logical :: wrong_args 
  !-------------------------------------------      
  number_args = IARGC() 

 !Exptdir, runmode, pdf, whichinput, np(Hneut), np(Hplus) , prefix
  nargs_datfile =  4+         2        +1
  np(Chineut)=0
  np(Chiplus)=0

  wrong_args= (number_args .ne. nargs_datfile)

  if( wrong_args )then
   write(*,*) "Incorrect number of parameters given on command line"
   call command_line_how2_for_hs
   stop "Error: command line entered incorrectly (see standard output for more info)"
  endif
      
  ! Read arguments into text strings.
  i=1
  
  temp=""
  call GETARG(i,temp)
  Exptdir = ""
  Exptdir = trim(temp)
  i=i+1
    
  temp=""
  call GETARG(i,temp)
  runmode = ""
  runmode = trim(temp)
  i=i+1
  
  temp=""
  call GETARG(i,temp)
  pdftemp = ""
  pdftemp = trim(temp)
  i=i+1  


  if(verify(pdftemp," 1234567890").gt.0)then 
    ! checks that the string nHtemp just contains the characters " 1234567890"
    ! the function verify is standard in fortran 95
    write(*,*)'Incorrect pdf: not a number.'
    write(*,*)'(you entered pdf="'//trim(adjustl(pdftemp))//'")'
    call command_line_how2_for_hs
    stop "Error: command line entered incorrectly (see standard output for more info)"
   endif

   read(pdftemp,*) pdf

  temp=""
  call GETARG(i,temp)
  whichinput = ""
  whichinput = trim(temp)
  i=i+1
  
  do x=1,2
   temp=""
   call GETARG(i,temp)
   nHtemp = ""
   nHtemp = trim(temp)
   i=i+1

   if(verify(nHtemp," 1234567890").gt.0)then 
    ! checks that the string nHtemp just contains the characters " 1234567890"
    ! the function verify is standard in fortran 95
    write(*,*)'Incorrect n'//trim(adjustl(pdesc(x)%short))//': not a number.'
    write(*,*)'(you entered n'//trim(adjustl(pdesc(x)%short))//&
&             '="'//trim(adjustl(nHtemp))//'")'
    write(*,*)'n'//trim(adjustl(pdesc(x)%short))//' is the number of '//&
&             trim(adjustl(pdesc(x)%long))//'s'
    call command_line_how2_for_hs
    stop "Error: command line entered incorrectly (see standard output for more info)"
   endif

   read(nHtemp,*) np(x)
  enddo

  call check_number_of_particles
  call check_whichinput 

 end subroutine getbasiccommandline_for_hs
!--------------------------------------------------------------------
 subroutine command_line_how2_for_hs
!--------------------------------------------------------------------
  use usefulbits, only : np,inputmethod

   write(*,*)'The correct syntax for the command line is:' 

   write(*,*)' ./HiggsSignals <expdata> <mode> <pdf> ',&
&            '<whichinput> <nHzero> <nHplus> <prefix>'
   write(*,*)'e.g.'
   write(*,*)' ./HiggsSignals latestresults peak 2 part 3 1 example_data/random/HB_randomtest50points_'

   write(*,*)'See HiggsSignals manual for more details.' 
   call flush(6)
 end subroutine command_line_how2_for_hs
!--------------------------------------------------------------------
 subroutine do_output_for_hs
 ! Writes output to file or screen, depending on whether 
 ! inputmethod='datfile' or inputmethod='website'
!--------------------------------------------------------------------
  use usefulbits, only : theo,ndat,inputmethod,np,fullHBres, &
          &              n_additional, file_id_common, file_id_common2,  & !input   
        &              whichanalyses,whichinput,pr,listprocesses,infile1  
  use usefulbits_hs, only : analyses, StrCompress, pdf, runmode,  HSres, HSvers, Exptdir,&
&  Nparam

!  use extra_bits_for_web
!  use S95tables  
!  use extra_bits_for_SLHA    

  implicit none
  !-----------------------------------internal
  integer :: i,j,ii,jj,x,n,n_sum,kk
  character(LEN=1) :: nHchar
  character(LEN=3) :: addchar,tempchar
  character(LEN=8) :: createdate
  character(LEN=10) :: createtime
  character(LEN=19) :: createdateandtime
  character(LEN=100) :: format43  
  character(LEN=500) :: columndescription
  type(listprocesses) :: proc  
  character(LEN=200):: descrip 
  double precision, allocatable :: Mhall(:,:)
  character(LEN=100) :: formatspec
  character(LEN=13) :: pdf_desc(3) = (/'box         ','gaussian    ','box+gaussian'/)
  !-------------------------------------------

 select case(inputmethod)     
 case('datfile')      
  select case(whichinput)
  case('SLHA')
   call HiggsSignals_outputSLHAdata(infile1)
  case('effC','part','hadr')

   if((np(Hneut)>0).or.(np(Hplus)>0))then
     allocate(Mhall(ndat,np(Hneut)+np(Hplus)))
     do jj=1,ndat
      ii=0
      if(np(Hneut)>0)then
       do i=1,np(Hneut)
        Mhall(jj,i)=theo(jj)%particle(Hneut)%M(i)
       enddo
       ii=np(Hneut)
      endif
      if(np(Hplus)>0)then
       do i=1,np(Hplus)
        Mhall(jj,ii+i)=theo(jj)%particle(Hplus)%M(i)
       enddo
      endif
     enddo
   else
    stop'error in subroutine do_output(1)'
   endif
 
   format43='(1I14,'
   write(nHchar,'(I1)')np(Hneut)
   if(np(Hneut)>0)format43=trim(adjustl(format43))//nHchar//'G16.6,'

   write(nHchar,'(I1)')np(Hplus)
   if(np(Hplus)>0)format43=trim(adjustl(format43))//nHchar//'G16.6,'

   format43=trim(adjustl(format43))//'3G16.6,3I6,1G16.6'

   write(addchar,'(I3)')n_additional
   if(n_additional>0)format43=trim(adjustl(format43))//','//trim(adjustl(addchar))//'G16.6'

   format43=trim(adjustl(format43))//')' 

   open(file_id_common,file=trim(infile1)//"HiggsSignals_results.dat") 
   if(runmode.eq.'peak'.or.runmode.eq.'both') then   
    open(file_id_common2,file=trim(infile1)//"peakobservables.dat") 
   endif
   
   call date_and_time(createdate,createtime)

   createdateandtime=createdate(7:8)//'.' &
         &  //createdate(5:6)//'.' &
         &  //createdate(1:4)//' at '&
         &  //createtime(1:2)//':' & 
         &  //createtime(3:4)  
   
   write(file_id_common,*)'# generated with HiggsSignals version '//&
&                         trim(adjustl(HSvers))//' on '//createdateandtime
   write(file_id_common,*)'# settings: '//adjustl(trim(Exptdir))//', '//whichinput//', '&
&                         //runmode//', '//adjustl(trim(pdf_desc(pdf)))
   write(file_id_common,*)'#'
   write(file_id_common,*)'# column abbreviations' 
   write(file_id_common,*)'#   n          : line id of input' 
   if(np(Hneut)>0)write(file_id_common,*)'#   Mh(i)      : ',&
&                                        'Neutral Higgs boson masses in GeV'
   if(np(Hplus)>0)write(file_id_common,*)'#   Mhplus(i)  : ',&
&                                        'Charged Higgs boson masses in GeV'
   write(file_id_common,*)'#   csq(mu)    : Chi^2 from the signal strengths observables'      
   write(file_id_common,*)'#   csq(mh)    : Chi^2 from the Higgs mass observables'          
   write(file_id_common,*)'#   csq(tot)   : total Chi^2' 
   write(file_id_common,*)'#   nobs(mu)   : number of signal strength observables'    
   write(file_id_common,*)'#   nobs(mh)   : number of Higgs mass observables'
   write(file_id_common,*)'#   nobs(tot)  : total number of observables'   
   write(file_id_common,"(A,I3)")'#   Pvalue     : Probability, given csq(tot) and ndf=nobs(tot)-',Nparam
   if(n_additional>0)then
    write(file_id_common,*)'#   additional : ',&
&       'optional additional data stored in <prefix>additional.dat (e.g. tan beta)'     
   endif
   write(file_id_common,*)'#'   
   
   if(runmode.eq.'peak'.or.runmode.eq.'both') then
    formatspec='(I3,8X,I10,3X,F6.2,1X,4F9.4,7X,A3,6X,F6.2,5X,F6.2,3X,I3,3X,A40,5X,A)'
    write(file_id_common2,*) "# List of peak observables used from " &
&                            //trim(adjustl(Exptdir))//":"
    write(file_id_common2,*)'#'   
    write(file_id_common2,*)'# Number        : HS internal numbering of analysis'
    write(file_id_common2,*)'# Analysis-ID   : identity number of analysis'
    write(file_id_common2,*)'# mh_obs        : Mass position of peak observable',&
&                         ' (hypothetical Higgs mass where signal strength was measured)'
    write(file_id_common2,*)'# mu_obs        : Observed (best-fit) signal strength value'
    write(file_id_common2,*)'# dmu_low       : Lower 1s (cyan band) signal strength value'
    write(file_id_common2,*)'# dmu_high      : Upper 1s (cyan band) signal strength value'
    write(file_id_common2,*)'# dmh_exp       : Experimental mass uncertainty of analysis'
    write(file_id_common2,*)'# collaboration : ',&
&                             'Exp. collaboration which published the data'
    write(file_id_common2,*)'# energy        : ',&
&                             'center-of-mass energy of the experiment (TeV)'
    write(file_id_common2,*)'# luminosity    : ',&
&                             'integrated luminosity used for the measurement (fb-1)'    
    write(file_id_common2,*)'# mass-chi2     : ',&
&                             'Mass measurement entering the Chi^2 (1/0=yes/no)'
    write(file_id_common2,*)'# description   : Higgs process of analysis'                   
    write(file_id_common2,*)'# reference     : ',&
&                             'Experimental publication of the data/analysis'
    write(file_id_common2,*)'#'       
    write(file_id_common2,*)'# Number Analysis-ID   mh_obs    mu_obs  dmu_low dmu_high ',&
& 				" dmh_exp  collaboration  energy  luminosity  mass-chi2               ", &
&               "          description     reference"
    write(file_id_common2,*) "#" 
    kk=0
    do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
     do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
      kk=kk+1
      write(file_id_common2,formatspec)kk,analyses(i)%id,analyses(i)%peaks(j)%mpeak,&
    & analyses(i)%peaks(j)%mu, analyses(i)%peaks(j)%dmulow,analyses(i)%peaks(j)%dmuup,	&
    & analyses(i)%table%deltam,analyses(i)%table%collaboration, analyses(i)%table%energy,&
    & analyses(i)%table%lumi, analyses(i)%table%mhchisq, &
    & trim(strcompress(analyses(i)%table%desc)), analyses(i)%table%label
     enddo
    enddo
   close(file_id_common2)
   endif 
      
   columndescription='#cols:  n'

   if(np(Hneut)>0)then
    do i=1,np(Hneut)
     write(tempchar,'(I3)')i
     columndescription=trim(columndescription)//'           Mh('//&
    &                  trim(adjustl(tempchar))//')  '
    enddo                                        
   endif

   if(np(Hplus)>0)then
    do i=1,np(Hplus)
     write(tempchar,'(I3)')i
     columndescription=trim(columndescription)//'       Mhplus('// &
   &                   trim(adjustl(tempchar))//')  '
    enddo
   endif

   columndescription=trim(columndescription)//'         csq(mu)         csq(mh)      '//&
&   '   csq(tot) nobs(mu) nobs(mh) nobs(tot)    Pvalue'

   if(n_additional>0)then
    do i=1,n_additional
     write(tempchar,'(I3)')i
     columndescription=trim(columndescription)//'   additional('// &
     &                 trim(adjustl(tempchar))//')   '
    enddo     
   endif

   write(file_id_common,*)trim(columndescription) 
   write(file_id_common,*)'#'  

   if(n_additional>0)then
    do jj=1,ndat                                                             
     write(file_id_common,fmt=trim(format43))jj,(Mhall(jj,i),i=1,np(Hneut)+np(Hplus)), &
          &      HSres(jj)%Chisq_mu,HSres(jj)%Chisq_mh,HSres(jj)%Chisq,              &
          &      HSres(jj)%nobs_mpred+HSres(jj)%nobs_peak_mu,HSres(jj)%nobs_peak_mh, &
          &      HSres(jj)%nobs,HSres(jj)%Pvalue,(theo(jj)%additional(i),i=1,n_additional)
    enddo  
   else
    do jj=1,ndat                  
   write(file_id_common,fmt=trim(format43))jj,(Mhall(jj,i),i=1,np(Hneut)+np(Hplus)), &
          &      HSres(jj)%Chisq_mu,HSres(jj)%Chisq_mh,HSres(jj)%Chisq,              &
          &      HSres(jj)%nobs_mpred+HSres(jj)%nobs_peak_mu,HSres(jj)%nobs_peak_mh, &
          &      HSres(jj)%nobs,HSres(jj)%Pvalue
    enddo  
   endif

  close(file_id_common)
  deallocate(Mhall) 
  case default
   stop'error in subroutine do_output (*2)'    
  end select
! case('website')               
!    do jj=1,ndat 
!-----NOTE: Need to change weboutput in extra_bits_for_web for the new all-Higgses mode!!!    
!      call weboutput(res(jj))
!    enddo 
 case('subrout') 
    select case(whichinput)
    case('SLHA')
       call HiggsSignals_outputSLHAdata(infile1)
    case default
       stop'error in subroutine do_output(c)'
    end select
 case default
    stop'error in subroutine do_output (*1)'                                              
 end select  
                                               
 end subroutine do_output_for_hs
!--------------------------------------------------------------------
subroutine setup_output_for_hs
!-Allocates the HBresults array.
!--------------------------------------------------------------------
  use usefulbits, only : ndat !input
  use usefulbits_hs, only : HSres !input
  
  implicit none

  allocate(HSres(ndat))

end subroutine setup_output_for_hs
!--------------------------------------------------------------------
subroutine get_ID_of_peakobservable(ii, ID)
 use usefulbits_hs, only : HSres
 
 implicit none
 integer, intent(in) :: ii
 integer, intent(out) :: ID

 if(ii.gt.size(HSres(1)%obsID,dim=1).or.ii.le.0) then
  write(*,*) 'WARNING in get_peakinfo_from_HSresults: Peak observable ',		&
&			 ii,' unkwown! Returning zero.'
  ID = 0
 else
  ID = HSres(1)%obsID(ii)
 endif 
 
end subroutine get_ID_of_peakobservable
!--------------------------------------------------------------------
subroutine get_peakinfo(obsID, muobs, dmuup, dmulow, mpeak, dm)
!--------------------------------------------------------------------
 use usefulbits_hs, only : obs
 implicit none
 integer, intent(in) :: obsID
 double precision, intent(out) :: muobs, dmuup, dmulow, mpeak, dm
 integer :: pos, ii
 
 pos = -1
 do ii=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(ii)%id.eq.obsID) then
   pos = ii
   exit   	
  endif
 enddo
  
 if(pos.ne.-1) then    
  mpeak = obs(pos)%peak%mpeak
  dm = obs(pos)%peak%dm  
  muobs = obs(pos)%peak%mu
  dmuup = obs(pos)%peak%dmuup
  dmulow = obs(pos)%peak%dmulow
 else
  write(*,*) "WARNING in get_peakinfo: ID unknown."
 endif
end subroutine get_peakinfo

!--------------------------------------------------------------------
subroutine get_peak_channels(obsID, Nc, IDs, efficiencies)
!--------------------------------------------------------------------
 use usefulbits_hs, only : obs
 implicit none
 integer, intent(in) :: obsID
 integer, intent(out) :: Nc
 integer, allocatable, intent(out) :: IDs(:)
 double precision, allocatable, intent(out) :: efficiencies(:)

 integer :: pos, ii
 
 pos = -1
 do ii=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(ii)%id.eq.obsID) then
   pos = ii
   exit   	
  endif
 enddo
  
 if(pos.ne.-1) then
  Nc = obs(pos)%table%Nc
  allocate(IDs(Nc),efficiencies(Nc))
  IDs = obs(pos)%table%channel_id
  efficiencies = obs(pos)%table%channel_eff
 else
  write(*,*) "WARNING in get_peakinfo: ID unknown."
 endif
end subroutine get_peak_channels
!--------------------------------------------------------------------
subroutine get_peakinfo_from_HSresults(obsID, mupred, domH, nHcomb)
!--------------------------------------------------------------------
 use usefulbits_hs, only : HSres
 
 implicit none
 integer, intent(in) :: obsID
 double precision, intent(out) :: mupred
 integer, intent(out) :: domH, nHcomb
 integer :: pos, ii
 
 pos = -1
 do ii=lbound(HSres(1)%nH,dim=1),ubound(HSres(1)%nH,dim=1)
  if(obsID.eq.HSres(1)%obsID(ii)) then
   pos = ii
   exit   	
  endif
 enddo
 
 if(pos.ne.-1) then  
  mupred = HSres(1)%mupred(ii)
  domH = HSres(1)%domH(ii)
  nHcomb = HSres(1)%nH(ii)
 else
  write(*,*) 'WARNING in get_peakinfo_from_HSresults: Peak observable ',		&
&			 obsID,' unkwown! Returning zeros.'
  mupred = 0.0D0
  domH = 0
  nHcomb = 0 
 endif
 
end subroutine get_peakinfo_from_HSresults
!--------------------------------------------------------------------
subroutine get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses)
!--------------------------------------------------------------------
 use usefulbits_hs, only : HSres
 
 implicit none
 integer, intent(out) :: ntotal, npeakmu, npeakmh, nmpred, nanalyses
  
 ntotal = HSres(1)%nobs
 npeakmu = HSres(1)%nobs_peak_mu
 npeakmh = HSres(1)%nobs_peak_mh
 nmpred = HSres(1)%nobs_mpred
 nanalyses = HSres(1)%nanalysis
 
end subroutine get_number_of_observables
!--------------------------------------------------------------------
subroutine HiggsSignals_SLHA_output(detailed)
!--------------------------------------------------------------------
use usefulbits, only : whichinput,just_after_run
use usefulbits_hs, only : SLHAdetailed
integer, intent(in) :: detailed

if(detailed.eq.0) SLHAdetailed = .False.

!use output, only : do_output
!!!integer, intent(in) :: i

  if(.not.just_after_run)then
   stop'subroutine run_HiggsSignals should be called'//&
   	   ' before subroutine HiggsSignals_SLHA_output' 
  endif  

  select case(whichinput)
  case('SLHA')
    call do_HiggsSignals_output
  case default
    stop'The subroutine HiggsSignals_SLHA_output should'//&
    ' only be used if an input SLHA file is given.'
  end select

end subroutine HiggsSignals_SLHA_output
!--------------------------------------------------------------------
 subroutine do_HiggsSignals_output
!--------------------------------------------------------------------
  use usefulbits, only : inputmethod,whichanalyses,whichinput,infile1  

  implicit none

 select case(inputmethod)     
 case('datfile')      
  select case(whichinput)
  case('SLHA')
   call HiggsSignals_outputSLHAdata(infile1)
   case default
    stop'error in subroutine do_output (*2)'    
   end select
 case('subrout') 
    select case(whichinput)
    case('SLHA')
       call HiggsSignals_outputSLHAdata(infile1)
    case default
       stop'error in subroutine do_output(c)'
    end select
 case default
    stop'error in subroutine do_output (*1)'                                              
 end select  
                                               
 end subroutine do_HiggsSignals_output
!--------------------------------------------------------------------
subroutine HiggsSignals_create_SLHA_output_default(detailed)
! Wrapper subroutine of HiggsSignals_create_SLHA_output with a fixed
! name of output file ("HS-output.slha").
!--------------------------------------------------------------------
 implicit none
 integer, intent(in) :: detailed

 call HiggsSignals_create_SLHA_output("HS-output.slha", detailed)

end subroutine HiggsSignals_create_SLHA_output_default
!--------------------------------------------------------------------
subroutine HiggsSignals_create_SLHA_output(infile, detailed)
!--------------------------------------------------------------------
 use usefulbits, only : just_after_run, file_id_common
 use usefulbits_hs, only : SLHAdetailed, newSLHAfile
 implicit none

  character(len=*),intent(in) :: infile
  integer, intent(in) :: detailed
  integer :: ios
  character(len=80) :: infile_80

  if(detailed.eq.0) SLHAdetailed = .False.

  write(infile_80,*) infile
  if(trim(adjustl(infile_80)).ne.trim(adjustl(infile))) then
   write(*,*) "WARNING: SLHA output file name is too long (> 80 elements)!"
  endif
  
  if(.not.just_after_run)then
   stop'subroutine run_HiggsSignals should be called'//&
   	   ' before subroutine HiggsSignals_SLHA_output' 
  endif  

  open(unit=file_id_common,file=trim(adjustl(infile_80)),status='new',iostat=ios)
  if(ios.ne.0) then
   write(*,*) 'Problem creating the SLHA file: $'//trim(adjustl(infile))//'$'
   write(*,*) 'This file may already exist. SLHA output writing is skipped.'
  else
   close(file_id_common)
   newSLHAfile=.True.
   call HiggsSignals_outputSLHAdata(infile_80)
   newSLHAfile=.False.
  endif

end subroutine HiggsSignals_create_SLHA_output
!************************************************************
 subroutine HiggsSignals_outputSLHAdata(infile) 
 !************************************************************ 
  use usefulbits, only : whichanalyses,file_id_common
  use extra_bits_for_SLHA, only : h
  use usefulbits_hs, only : HSvers, HSres, runmode, pdf, analyses, Exptdir, SLHAdetailed,&
 &                          newSLHAfile 
  use SLHA_manip		! n.b.: needed for readSLHAfile, writeSLHAfile_except, finishwithSLHA
  
  implicit none
  !--------------------------------------input
  character(len=80),intent(in) :: infile
  !-----------------------------------internal
  integer :: i,j,ios,a
  integer :: k_out, Nh
  character(LEN=200):: descrip, formatstring
  !-------------------------------------------              

  open(file_id_common,file=trim(adjustl(infile)),status='old',iostat=ios) 
  if(ios.ne.0)then 
    write(*,*)'problem opening the SLHA file: $'//trim(adjustl(infile))//'$'
  else

    k_out=file_id_common

    if(.not.newSLHAfile) then
     call readSLHAfile(file_id_common)
     rewind(file_id_common)
     call writeSLHAfile_except(k_out,'HiggsSignalsResults','HiggsSignalsPeakObservables',&
& 	  'HiggsSignalsMassCenteredObservables')
    endif 

    write(k_out,'(a)')'BLOCK HiggsSignalsResults'
    write(k_out,'(3X,I2,17X,A2,A,A2,15X,A)') 0,'||',trim(adjustl(HSvers)),'||',		&
&    '# HiggsSignals version'
    write(k_out,'(3X,I2,9X,A2,A,A2,15X,A)')1,'||',trim(adjustl(Exptdir)),'||',	&
&    '# experimental data set'
	select case(runmode)
	 case('peak')
    write(k_out,'(3X,I2,25X,I1,15X,A)') 2,1,								&
&   '# Chi-squared method ("peak"(1) or "mass"(2)-centered or "both"(3))'
	 case('mass')
    write(k_out,'(3X,I2,25X,I1,15X,A)') 2,2,								&
&   '# Chi-squared method ("peak"(1) or "mass"(2)-centered or "both"(3))'
	 case('both')
    write(k_out,'(3X,I2,25X,I1,15X,A)') 2,3,								&
&   '# Chi-squared method ("peak"(1) or "mass"(2)-centered or "both"(3))'
	end select
	write(k_out,'(3X,I2,11X,I15,15X,A)') 3,pdf,											&
&	'# Parametrization of Higgs mass uncertainty (1:box, 2:gaussian, 3:box+gaussian)'
    write(k_out,'(3X,I2,11X,I15,15X,A)') 4,HSres(1)%nobs_peak_mu,			&
&	'# Number of signal strength peak observables'
    write(k_out,'(3X,I2,11X,I15,15X,A)') 5,HSres(1)%nobs_peak_mh,			&
&   '# Number of Higgs mass peak observables'    
    write(k_out,'(3X,I2,11X,I15,15X,A)') 6,HSres(1)%nobs_mpred,				&
&   '# Number of mass-centered observables'
    write(k_out,'(3X,I2,11X,I15,15X,A)') 7,HSres(1)%nobs,					&
&   '# Number of observables (total)'
    write(k_out,'(3X,I2,11X,F15.8,15X,A)') 8,HSres(1)%Chisq_peak_mu,		&
&   '# chi^2 from signal strength peak observables'    
    write(k_out,'(3X,I2,11X,F15.8,15X,A)') 9,HSres(1)%Chisq_mh,				&
&   '# chi^2 from Higgs mass peak observables'
    write(k_out,'(3X,I2,11X,F15.8,15X,A)') 10,HSres(1)%Chisq_mpred,			&
&   '# chi^2 from mass-centered observables'   
    write(k_out,'(3X,I2,11X,F15.8,15X,A)') 11,HSres(1)%Chisq_mu,			&
&   '# chi^2 from signal strength (total)'
    write(k_out,'(3X,I2,11X,F15.8,15X,A)') 12,HSres(1)%Chisq,'# chi^2 (total)'
    write(k_out,'(3X,I2,11X,F15.8,15X,A)') 13,HSres(1)%Pvalue,				&
&   '# Probability (total chi^2, total number observables)'
    
    if(SLHAdetailed) then
     select case(runmode)
	  case('peak')	 
	   call writeSLHAblock_HiggsSignalsPeakObservables(k_out)
	  
	  case('mass')
	   call writeSLHAblock_HiggsSignalsMassCenteredObservables(k_out)
	  
	  case('both') 
	   call writeSLHAblock_HiggsSignalsPeakObservables(k_out)	  
	   call writeSLHAblock_HiggsSignalsMassCenteredObservables(k_out)

	  case default
	 end select 
    endif 
    SLHAdetailed=.True.

    close(file_id_common)
    if(k_out.ne.file_id_common)close(k_out)

    if(.not.newSLHAfile) then
     call finishwithSLHA
    endif 

  endif
 end subroutine HiggsSignals_outputSLHAdata
!--------------------------------------------------------------------
 subroutine writeSLHAblock_HiggsSignalsPeakObservables(k_out)
  use usefulbits_hs, only : analyses, StrCompress
  use extra_bits_for_SLHA, only : h  
  implicit none
  integer, intent(in) :: k_out
  integer :: i,j,a, Nh
  character(LEN=200):: formatstring
  

  write(k_out,'(a)')'BLOCK HiggsSignalsPeakObservables'
  write(k_out,'(A6,3X,A4,25X,A5,3X,A)')'#  OBS','FLAG','VALUE','# DESCRIPTION'	  
  a=0	
   do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
    do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
	a=a+1
    Nh = size(analyses(i)%peaks(j)%Higgs_comb)		
    write(k_out,'(3X,I3,4X,I2,1X,I30,3X,A)')a,1,analyses(i)%table%id,'# Analysis ID'
    write(k_out,'(3X,I3,4X,I2,5X,A2,A,A2,8X,A)')a,2,'||',trim(analyses(i)%table%label),	&
&    '||','# Reference to publication'
    write(k_out,'(3X,I3,4X,I2,5X,A2,A,A2,6X,A)')a,3,'||',								&
&    trim(strcompress(analyses(i)%table%desc)),'||','# Description (Search channel)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,4,analyses(i)%table%energy,			&
&   '# Center-of-mass energy'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,5,analyses(i)%table%lumi,			&
&   '# Luminosity'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,6,analyses(i)%table%dlumi*100.0D0,	&
&   '# Luminosity uncertainty (in %)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,7,analyses(i)%table%deltam,		&
&	'# Mass resolution (GeV)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,8,analyses(i)%peaks(j)%mpeak,	&
&	'# Mass value at peak position (in GeV)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,9,analyses(i)%peaks(j)%mu,	&
&	'# Observed signal strength modifier (mu)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,10,analyses(i)%peaks(j)%dmulow,	&
&	'# Lower 68%C.L. uncertainty on observed mu'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,11,analyses(i)%peaks(j)%dmuup,	&
&	'# Upper 68%C.L. uncertainty on observed mu'
!    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,12,analyses(i)%peaks(j)%dmulow0sq,	&
!&	'# Intrinsic squared lower 68%C.L. cyan band value on observed mu (corrected by HS)'
!    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,13,analyses(i)%peaks(j)%dmuup0sq,	&
!&	'# Intrinsic squared upper 68%C.L. cyan band value on observed mu (corrected by HS)'
!	write(formatstring,'(A13,I2,A2,I2,A8)') '(3X,I3,4X,I2,',31-1*Nh,'X,',Nh,'I1,3X,A)'
!	write(k_out,formatstring) a,12,analyses(i)%peaks(j)%Higgs_comb,				&
!&	'# Higgs combination (= indices of Higgs bosons which have been used)'
	write(formatstring,'(A13,I2,A3,I1,A1,I1,A6)') &
&	                   '(3X,I3,4X,I2,',31-1*Nh,'X,I',Nh,'.',Nh,',3X,A)'
	write(k_out,formatstring) a,12, &
&     convert_Higgscomb_to_bin(analyses(i)%peaks(j)%Higgs_comb), &
&     '# Assigned Higgs combination'
!-------TODO: List the pdg numbers here of these combinations!
     write(k_out,'(3X,I3,4X,I2,23X,I8,3X,A)')a,13,analyses(i)%peaks(j)%domH,&
&	 '# Index of dominant Higgs boson'
	if(analyses(i)%peaks(j)%domH.eq.0) then
     write(k_out,'(3X,I3,4X,I2,23X,A8,3X,A)')a,14,'NaN',						&
&    '# pdg number of dominant Higgs boson'
	 write(k_out,'(3X,I3,4X,I2,23X,A8,3X,A)')a,15,'NaN',						&
&	 '# Mass of the dominant Higgs boson'
	 write(k_out,'(3X,I3,4X,I2,23X,A8,3X,A)')a,16,'NaN',						&
&	 '# Signal strength modifier of the dominant Higgs boson'
    else
     write(k_out,'(3X,I3,4X,I2,23X,I8,3X,A)')a,14,h(analyses(i)%peaks(j)%domH),	&
&    '# pdg number of dominant Higgs boson'
	 write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,15,	&
&	 analyses(i)%peaks(j)%Higgses(analyses(i)%peaks(j)%domH)%m,		&
&	 '# Mass of dominant Higgs boson'
	 write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,16,	&
&	 analyses(i)%peaks(j)%Higgses(analyses(i)%peaks(j)%domH)%mu,	&
&	 '# Signal strength modifier of dominant Higgs boson'
	endif
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,17,analyses(i)%peaks(j)%total_mu,&
&	'# Total predicted signal strength modifier mu'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,18,analyses(i)%peaks(j)%chisq_mu,&
&	'# Chi-squared value (mu-part)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,19,analyses(i)%peaks(j)%chisq_mh,&
&	'# Chi-squared value (mh-part)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,20,analyses(i)%peaks(j)%chisq_tot,&
&   '# Chi-squared value (total)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,21,analyses(i)%peaks(j)%chisq_max,&
&	'# Chi-squared value for no predicted signal (mu=0)'
   enddo
  enddo 

end subroutine writeSLHAblock_HiggsSignalsPeakObservables
!--------------------------------------------------------------------
 subroutine writeSLHAblock_HiggsSignalsMassCenteredObservables(k_out)
!--------------------------------------------------------------------
  use usefulbits_hs, only : analyses, StrCompress
  use usefulbits, only : np, Hneut
!  use extra_bits_for_SLHA, only : h  
  implicit none
  integer, intent(in) :: k_out
  integer :: i,j,a, Nh
  character(LEN=200):: formatstring
  Nh = np(Hneut)

  write(k_out,'(a)')'BLOCK HiggsSignalsMassCenteredObservables'
  write(k_out,'(A6,3X,A4,25X,A5,3X,A)')'#  OBS','FLAG','VALUE','# DESCRIPTION'	  
  a=0	
  do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
   do j=lbound(analyses(i)%mpred%mp_Higgses,dim=1),										&
&	    ubound(analyses(i)%mpred%mp_Higgses,dim=1)
    a=a+1
    write(k_out,'(3X,I3,4X,I2,1X,I30,3X,A)')a,1,analyses(i)%table%id,'# Analysis ID'
    write(k_out,'(3X,I3,4X,I2,5X,A2,A,A2,8X,A)')a,2,'||',trim(analyses(i)%table%label),	&
&    '||','# Reference to publication'
    write(k_out,'(3X,I3,4X,I2,5X,A2,A,A2,6X,A)')a,3,'||',								&
&    trim(strcompress(analyses(i)%table%desc)),'||','# Description (Search channel)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,4,analyses(i)%table%energy,				&
&   '# Center-of-mass energy'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,5,analyses(i)%table%lumi,				&
&   '# Luminosity'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,6,analyses(i)%table%dlumi*100.0D0,		&
&   '# Luminosity uncertainty (in %)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,7,analyses(i)%table%deltam,				&
&	'# Mass resolution (GeV)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,8,analyses(i)%mpred%mp_Higgses(j)%m,	&
    '# Mass of tested Higgs boson (GeV)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,9,analyses(i)%mpred%mp_Higgses(j)%dm,	&
&   '# Mass uncertainty of tested Higgs boson (GeV)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,10,analyses(i)%mpred%mp_Higgses(j)%mu,	&
&   '# Signal strength of tested Higgs boson'
    write(k_out,'(3X,I3,4X,I2,23X,I8,3X,A)')a,11,										&
&   size(analyses(i)%mpred%mp_Higgses(j)%Higgses,dim=1),								&
&   '# Number of combined Higgs bosons'
	write(formatstring,'(A13,I2,A3,I1,A1,I1,A6)') &
&	                   '(3X,I3,4X,I2,',31-1*Nh,'X,I',Nh,'.',Nh,',3X,A)'
	write(k_out,formatstring) a,12, &
&     convert_Higgscombint_to_bin(analyses(i)%mpred%mp_Higgses(j)%Higgscomb), &
&     '# Combined Higgs boson code'
    write(k_out,'(3X,I3,4X,I2,23X,F8.2,3X,A)')a,13,										&
&   analyses(i)%mpred%mp_Higgses(j)%m_obs,'# Observed mass value (GeV)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,14,										&
&   analyses(i)%mpred%mp_Higgses(j)%mu_obs,												&
&   '# Observed signal strength  (from experiment)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,15,										&
&   analyses(i)%mpred%mp_Higgses(j)%dmu_low_obs,										&
&   '# Lower 68%C.L. uncertainty of observed signal strength'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,16,										&
&   analyses(i)%mpred%mp_Higgses(j)%dmu_up_obs,											&
&   '# Upper 68%C.L. uncertainty of observed signal strength'
!    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,16,										&
!&   analyses(i)%mpred%mp_Higgses(j)%dmu_low0_obs,										&
!&   '# Intrinsic lower 68%C.L. cyan band value of observed mu (corrected by HS)'
!    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,17,										&
!&   analyses(i)%mpred%mp_Higgses(j)%dmu_up0_obs,										&
!&   '# Intrinsic upper 68%C.L. cyan band value of observed mu (corrected by HS)'
    write(k_out,'(3X,I3,4X,I2,23X,F8.4,3X,A)')a,17,analyses(i)%mpred%mp_Higgses(j)%chisq,&
&	'# Chi-squared value'
	enddo
   enddo	

end subroutine writeSLHAblock_HiggsSignalsMassCenteredObservables
!--------------------------------------------------------------------
function convert_Higgscomb_to_bin(Higgscomb)
!--------------------------------------------------------------------
 implicit none
 integer, dimension(:), intent(in) :: Higgscomb
 integer :: convert_Higgscomb_to_bin
 integer :: i, j
 integer :: decnum
 
 decnum=0
 do i=lbound(Higgscomb,dim=1),ubound(Higgscomb,dim=1)
  if(Higgscomb(i).eq.0) decnum=decnum+0
  if(Higgscomb(i).eq.1) decnum=decnum+1
  if(Higgscomb(i).eq.2) decnum=decnum+2
  if(Higgscomb(i).eq.3) decnum=decnum+4
  if(Higgscomb(i).eq.4) decnum=decnum+8  
  if(Higgscomb(i).eq.5) decnum=decnum+16  
  if(Higgscomb(i).eq.6) decnum=decnum+32  
  if(Higgscomb(i).eq.7) decnum=decnum+64  
  if(Higgscomb(i).eq.8) decnum=decnum+128  
  if(Higgscomb(i).eq.9) decnum=decnum+256
 enddo
 
 convert_Higgscomb_to_bin = convert_dec_to_bin(decnum)
 
end function convert_Higgscomb_to_bin
!--------------------------------------------------------------------
function convert_Higgscombint_to_bin(Higgscomb)
!--------------------------------------------------------------------
 implicit none
 integer, intent(in) :: Higgscomb
 integer :: convert_Higgscombint_to_bin
 integer :: i, j
 integer :: decnum, Higgscombtmp
 
 Higgscombtmp = Higgscomb
 
 decnum=0
 do while(Higgscombtmp > 0)
  i=MOD(Higgscombtmp,10)
  if(i.eq.0) decnum=decnum+0
  if(i.eq.1) decnum=decnum+1
  if(i.eq.2) decnum=decnum+2
  if(i.eq.3) decnum=decnum+4
  if(i.eq.4) decnum=decnum+8  
  if(i.eq.5) decnum=decnum+16  
  if(i.eq.6) decnum=decnum+32  
  if(i.eq.7) decnum=decnum+64  
  if(i.eq.8) decnum=decnum+128  
  if(i.eq.9) decnum=decnum+256
  Higgscombtmp=(Higgscombtmp-i)/10
 enddo

 convert_Higgscombint_to_bin = convert_dec_to_bin(decnum)
  
end function convert_Higgscombint_to_bin
!------------------------------------------------------------
subroutine read_matrix_from_file(length, filename, m, status)
!------------------------------------------------------------
 use usefulbits, only : file_id_common3
 implicit none
! use store_pathname_hs, only : pathname_HS
! use usefulbits_hs, only : delta_rate
 integer, intent(in) :: length
 character(LEN=*), intent(in) :: filename
 double precision, dimension(:,:), intent(inout) :: m
 logical, intent(out) :: status
 integer :: n, ios
 double precision, dimension(length) :: dblline
  
 status=.False.
 if((size(m,dim=1).ne.length).or.(size(m,dim=2).ne.length)) then
  write(*,*) 'Problem in subroutine read_matrix_from_file!'
 else
  open(file_id_common3, file=filename,form='formatted') 
  n=0
  do
   read(file_id_common3,*, iostat=ios) dblline
   if(ios.ne.0) exit
   n=n+1
   m(n,:)=dblline(:)
  enddo
 
  if(n.ne.length) then
   write(*,*) 'Problem with reading in the matrix from ',filename 
  else
   status=.True.
  endif
 endif 
 
 close(file_id_common3)

end subroutine read_matrix_from_file
!--------------------------------------------------------------------
function convert_dec_to_bin(dec)
!--------------------------------------------------------------------
! implicit none
 integer, intent(in) :: dec
 integer :: convert_dec_to_bin
 integer :: bineq
 integer :: power
 integer :: temp

 bineq = 0
 power = 0
 temp = dec

 do while(temp > 0)
  bineq = bineq + (MOD(temp, 2) * (MOD(temp, 2) * (10 ** power)))
  power = power + 1
  temp = temp / 2
 enddo

 convert_dec_to_bin = bineq
 
end function convert_dec_to_bin
!--------------------------------------------------------------------
end module io
!--------------------------------------------------------------------