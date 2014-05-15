! This file is part of HiggsBounds
!  -KW
!********************************************************
module output
!********************************************************
 use usefulbits, only : Hneut,Hplus,file_id_common    
 contains      

 !********************************************************      
 subroutine setup_output
 ! Creates the key to the process numbers used
 ! (by the way, note that these numbers are dependent on the 
 ! setting of whichanalyses)
 ! Sets 'rep'. This is the number of processes that will be checked
 ! against the experimentally measured limit
 ! e.g. if rep=2, the process with the highest stat sens will be checked, 
 ! but so will the process with the 2nd highest stat.sens.
 ! (can be useful, despitethe statistical intepretation of results from
 !  the 2nd highest process being hazy)
 ! Allocates the array 'res', which holds stuff to do with output
 !********************************************************
  use usefulbits, only : ndat,inputmethod, & !input
                       & rep, & !output
                       & res,fullHBres, & !allocated
                       & allocate_if_stats_required,infile1,whichinput
  implicit none
  !-----------------------------------internal
  integer :: i                    
  character(len=100) :: start_of_filename !must be same size as infile1 or bigger
  !-------------------------------------------

  if(len(start_of_filename).lt.len(infile1))then
     stop'problem in subroutine setup_output 2'
  endif

  select case(inputmethod)
  case('datfile')
   !create Key to process numbers in file Key.dat   
   
   select case(whichinput)
   case('effC','part','hadr')
    start_of_filename=infile1
   case('SLHA')
    start_of_filename=''
   case default 
     stop'problem in subroutine setup_output 1'
   end select
   call createKey(start_of_filename)

   rep=1
  case('website') 
   rep=3          !if change this, also need to change subroutine weboutput
                  !in module extra_bits_for_web       
  case('subrout')
   !create Key to process numbers in file Key.dat      
#ifndef WEBVERSION
   call createKey(infile1)
#endif
   if(allocated(allocate_if_stats_required))then
     rep=2         
   else
     rep=1          !if change this, also need to change subroutine run_HiggsBounds_effC
   endif
  end select
      
  allocate(res(ndat))
  do i=1,ndat
   allocate(res(i)%chan(rep))
   allocate(res(i)%obsratio(rep))
   allocate(res(i)%axis_i(rep))
   allocate(res(i)%axis_j(rep))
   allocate(res(i)%sfactor(rep))
   allocate(res(i)%allowed95(rep))
   allocate(res(i)%ncombined(rep))
   allocate(res(i)%channelselection(rep))
  enddo

!--NEW in HS-4:
  allocate(fullHBres(ndat))


  if(allocated(allocate_if_stats_required))then
    do i=1,ndat
      res(i)%channelselection='full'
      res(i)%channelselection(2)='clsb'
    enddo
  else
    do i=1,ndat
      res(i)%channelselection='full'
    enddo
  endif


      
 end subroutine setup_output
 
 !********************************************************      
 subroutine do_output
 ! Writes output to file or screen, depending on whether 
 ! inputmethod='datfile' or inputmethod='website'
 !******************************************************** 
  use usefulbits, only : res,theo,ndat,inputmethod,np,vers,fullHBres, &
          &              n_additional, & !input   
        &              whichanalyses,whichinput,pr,listprocesses,infile1,infile2  
  use extra_bits_for_web
  use S95tables  
  use extra_bits_for_SLHA    

  implicit none
  !-----------------------------------internal
  integer :: i,ii,jj,x,n,n_sum
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
  !-------------------------------------------

 select case(inputmethod)     
 case('datfile')      
  select case(whichinput)
  case('SLHA')
   call outputSLHAdata(infile2)
  case('effC','part','hadr')

   n_sum=0
   do i=1,ndat
     n_sum=n_sum+fullHBres(i)%allowed95
   enddo

   if(n_sum.le.(-ndat))then ! i.e. no good data sets
    !write(*,*)'Check that M>0 and the sum of the branching ratios for each Higgs is <=1'
    write(*,*)'no valid data sets: check that particle mass>0 for at least one data set'
   endif

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

   format43=trim(adjustl(format43))//'2I6,G16.6,1I6'

   write(addchar,'(I3)')n_additional
   if(n_additional>0)format43=trim(adjustl(format43))//','//trim(adjustl(addchar))//'G16.6'

   format43=trim(adjustl(format43))//')' 

   open(file_id_common,file=trim(infile1)//"HiggsBounds_results.dat") 
   
   call date_and_time(createdate,createtime)

   createdateandtime=createdate(7:8)//'.' &
         &  //createdate(5:6)//'.' &
         &  //createdate(1:4)//' at '&
         &  //createtime(1:2)//':' & 
         &  //createtime(3:4)  
   
   write(file_id_common,*)'# generated with HiggsBounds version '//trim(adjustl(vers))//' on '//createdateandtime
   write(file_id_common,*)'# settings: '//whichanalyses//', '//whichinput   
   write(file_id_common,*)'#'
   write(file_id_common,*)'# column abbreviations' 
   write(file_id_common,*)'#   n          : line id of input' 
   if(np(Hneut)>0)write(file_id_common,*)'#   Mh(i)      : Neutral Higgs boson masses in GeV'
   if(np(Hplus)>0)write(file_id_common,*)'#   Mhplus(i)  : Charged Higgs boson masses in GeV'
   write(file_id_common,*)'#   HBresult   : scenario allowed flag (1: allowed, 0: excluded, -1: unphysical)'      
   write(file_id_common,*)'#   chan       : most sensitive channel (see below). chan=0 if no channel applies'          
   write(file_id_common,*)'#   obsratio   : ratio [sig x BR]_model/[sig x BR]_limit (<1: allowed, >1: excluded)' 
   write(file_id_common,*)'#   ncomb      : number of Higgs bosons combined in most sensitive channel'
   if(n_additional>0)then
    write(file_id_common,*)'#   additional : optional additional data stored in <prefix>additional.dat (e.g. tan beta)'     
   endif
   write(file_id_common,*)'#'   
   write(file_id_common,*)'# channel numbers used in this file'  
   
   do x=1,ubound(pr,dim=1)
    n=0
    do jj=1,ndat   
      if(x.eq.fullHBres(jj)%chan)n=n+1
    enddo
    if(n.ge.1)then
      proc=pr(x)   
      write(tempchar,'(I3)')x          
      call outputproc(proc,file_id_common,descrip,1)
      write(file_id_common,*)'#         '//trim(tempchar)//' :'//trim(descrip)       
    endif     
   enddo   
   
   write(file_id_common,*)'# (for full list of processes, see Key.dat)'   
   write(file_id_common,*)'#'     
      
   columndescription='#cols:  n'

   if(np(Hneut)>0)then
    do i=1,np(Hneut)
     write(tempchar,'(I3)')i
     columndescription=trim(columndescription)//'           Mh('//trim(adjustl(tempchar))//')  '
    enddo                                        
   endif

   if(np(Hplus)>0)then
    do i=1,np(Hplus)
     write(tempchar,'(I3)')i
     columndescription=trim(columndescription)//'       Mhplus('//trim(adjustl(tempchar))//')  '
    enddo
   endif

   columndescription=trim(columndescription)//'  HBresult  chan    obsratio     ncomb '

   if(n_additional>0)then
    do i=1,n_additional
     write(tempchar,'(I3)')i
     columndescription=trim(columndescription)//'   additional('//trim(adjustl(tempchar))//')   '
    enddo     
   endif

   write(file_id_common,*)trim(columndescription) 
   write(file_id_common,*)'#'  

   if(n_additional>0)then
    do jj=1,ndat                                                             
       write(file_id_common,fmt=trim(format43))jj,(Mhall(jj,i),i=1,np(Hneut)+np(Hplus)), &
              &      fullHBres(jj)%allowed95,fullHBres(jj)%chan,fullHBres(jj)%obsratio,  &
              &      fullHBres(jj)%ncombined,(theo(jj)%additional(i),i=1,n_additional)
    enddo  
   else
    do jj=1,ndat                                                             
       write(file_id_common,fmt=trim(format43))jj,(Mhall(jj,i),i=1,np(Hneut)+np(Hplus)), &
              &      fullHBres(jj)%allowed95,fullHBres(jj)%chan,fullHBres(jj)%obsratio,  &
              &      fullHBres(jj)%ncombined
    enddo  
   endif


   close(file_id_common)
   deallocate(Mhall) 
  case default
   stop'error in subroutine do_output (*2)'    
  end select
 case('website')               
    do jj=1,ndat 
!-----NOTE: Need to change weboutput in extra_bits_for_web for the new all-Higgses mode!!!    
      call weboutput(res(jj))
    enddo 
 case('subrout') 
    select case(whichinput)
    case('SLHA')
       call outputSLHAdata(infile1)
    case default
       stop'error in subroutine do_output(c)'
    end select
 case default
    stop'error in subroutine do_output (*1)'                                              
 end select  
                                               
 end subroutine do_output
 !********************************************************      
 subroutine createKey(file_beginning)
 ! writes the key to the process numbers
 ! (by the way, note that the process numbers are dependent on the 
 ! setting of whichanalyses)
 !********************************************************
  use S95tables
  use usefulbits, only : np,ntot,pr,vers,whichanalyses !input
  implicit none
  character(len=*),intent(in) :: file_beginning
  !-----------------------------------internal
  integer :: i,j,n   
  character(LEN=45) :: tableno
  character(LEN=200):: descrip 
  !-------------------------------------------      
      
  open(file_id_common,file=trim(file_beginning)//'Key.dat')

  132   format('***********         Key to Process Numbers       *************')
                            
  131   format('**************************************************************',/ &
              ,'      ',      A      ,/ )      
                    
  130   format('process ',1I4)              
  
  133   format(' ------------------------------------------------------------ ')
  134   format('**************************************************************')

  write(file_id_common,132)
  write(file_id_common,*)'This key has been generated with HiggsBounds version '//trim(adjustl(vers))//' '
  write(file_id_common,*)'with the setting whichanalyses='//whichanalyses
  write(file_id_common,*) 
  write(file_id_common,134)
  write(file_id_common,*) 
  write(file_id_common,130)0
  write(file_id_common,*)' no process applies'

       
  do n=1,ntot    
   i=pr(n)%findi
   j=pr(n)%findj
   select case(pr(n)%ttype)
   case(1)
    tableno=trim(S95_t1(pr(n)%tlist)%label)
   case(2)
    tableno=trim(S95_t2(pr(n)%tlist)%label)  
   end select
        
   if((j.eq.1).and.(i.eq.1))write(file_id_common,131)trim(tableno)
        
   write(file_id_common,130)n                

   if(pr(n)%ttype.eq.1)then
    call outputproc(pr(n),file_id_common,descrip,1)
    write(file_id_common,*)trim(descrip)
   elseif(pr(n)%ttype.eq.2)then
    if((i.eq.j).and.(S95_t2(pr(n)%tlist)%needs_M2_gt_2M1))then
    ! these processes involve Br(hj->hj+hj) and so never exist
      write(file_id_common,*)' ***'                 
    else
      call outputproc(pr(n),file_id_common,descrip,1)
      write(file_id_common,*)trim(descrip)
    endif 
   else 
    stop 'error in subroutine createKey'
   endif
        
   select case(pr(n)%ttype)
   case(1)
    if((i.eq.np( S95_t1(pr(n)%tlist)%particle_x  )).and.(j.ne.np( S95_t1(pr(n)%tlist)%particle_x  )))write(file_id_common,133) 
   case(2)       
    if((i.eq.np( S95_t2(pr(n)%tlist)%particle_x1 )).and.(j.ne.np( S95_t2(pr(n)%tlist)%particle_x2 )))write(file_id_common,133)
   end select         
  enddo     
      
  close(file_id_common)

 end subroutine createKey
!******************************************************************

end module output
!******************************************************************neut
