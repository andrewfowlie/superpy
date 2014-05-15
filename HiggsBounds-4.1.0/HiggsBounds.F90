! This file is part of HiggsBounds
!  -KW
!************************************************************      
program HiggsBounds
!************************************************************
! HiggsBounds uses theoretical cross section and branching ratios
! from theories containing an arbitrary number of neutral
! Higgs particles and compares to the experimental limits 
! from LEP, Tevatron and LHC
!************************************************************
 use usefulbits, only : np, Hneut, Hplus, whichanalyses,debug,inputmethod, &
&                       file_id_debug1,file_id_debug2,HiggsBounds_info
 use input, only : setup_input,do_input
 use output, only : do_output !setup_output

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif
       
 implicit none 
 !-----------------------------------internal      
 logical :: messages
 integer :: HBresult,chan,ncombined
 double precision :: obsratio
 !-------------------------------------------               

!#define DEBUGGING
#ifdef DEBUGGING
  debug=.True.
#else
  debug=.False.
#endif

#ifdef WEBVERSION
 inputmethod='website'  
#else
 inputmethod='datfile'    !(inputmethod='subrout' is also possible, but do not set that here)
#endif

 messages=debug.or.(inputmethod=='datfile')

 if(inputmethod.eq.'datfile')call HiggsBounds_info

 if(debug)then  
  open(file_id_debug2,file='debug_predratio.txt')
  open(file_id_debug1,file='debug_channels.txt')
 endif
                                    
 if(messages)write(*,*)'doing some preliminary tasks...'      ; call flush(6)
 call setup_input

 call initialize_HiggsBounds(np(Hneut),np(Hplus),whichanalyses)
      
 if(messages)write(*,*)'getting theoretical input...'          ; call flush(6)
 call do_input
             
 if(messages)write(*,*)'compare each data point to the experimental bounds...' 
                                                                 call flush(6)            
 call run_HiggsBounds(HBresult,chan,obsratio,ncombined)
      
 if(debug)then
  close(file_id_debug2)
  close(file_id_debug1)
 endif

 if(messages)write(*,*)'beginning output...'                   ; call flush(6)
 call do_output

 if(messages)write(*,*)'finishing off...'                      ; call flush(6)
 call finish_HiggsBounds
      
 if(messages)write(*,*)'finished'                              ; call flush(6)             
end program HiggsBounds


      
