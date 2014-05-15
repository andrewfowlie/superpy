
!******************************************************
program HBSLHAinputblocksfromFH
!
! Run with
!   ./HBSLHAinputblocksfromFH <path to SLHA file>
! where <path to SLHA file> is the path to the SLHA file providing the 
! input for FeynHiggs
!
! The program will create a new file, called
!    <path to SLHA file>.fh
! which will contain the blocks from the original SLHA file and 
! the results from FeynHiggs to the SLHA file 
! in the standard SLHA format, but also the blocks
!     HiggsBoundsInputHiggsCouplingsFermions
!       and
!     HiggsBoundsInputHiggsCouplingsBosons
! which are defined in the HiggsBounds manual.
!
!******************************************************
 use extra_bits_for_SLHA, only : addcouplingsblocktoSLHAfile
 use usefulbits, only : sqcouplratio,np,Hneut,Hplus,allocate_sqcouplratio_parts     
#ifdef NAGf90Fortran
 use F90_UNIX_ENV, only : iargc,getarg
#endif

 implicit none

 character(len=300) :: infile,outfile
 type(sqcouplratio) :: gsq(1)

 character(LEN=300) :: temp
 integer :: number_args
#ifndef NAGf90Fortran
 integer :: iargc
#endif

 integer :: i

 number_args = IARGC() 
  
 if( number_args .ne. 1)then
   write(*,*)'!******************************************************' 
   write(*,*)' '
   write(*,*)' Run with'
   write(*,*)'     ./HBSLHAinputblocksfromFH <path to SLHA file>'
   write(*,*)' where <path to SLHA file> is the path to the SLHA file providing the' 
   write(*,*)' input for FeynHiggs'
   write(*,*)''
   write(*,*)' The program will create a new file, called'
   write(*,*)'     <path to SLHA file>.fh'
   write(*,*)' which will contain the blocks from the original SLHA file and' 
   write(*,*)' the results from FeynHiggs to the SLHA file '
   write(*,*)' in the standard SLHA format, but also the blocks'
   write(*,*)'     HiggsBoundsInputHiggsCouplingsFermions'
   write(*,*)'       and'
   write(*,*)'     HiggsBoundsInputHiggsCouplingsBosons'
   write(*,*)' which are defined in the HiggsBounds manual.'
   write(*,*)''
   write(*,*)'!******************************************************'
  
   stop "Incorrect number of arguments given to HBSLHAinputblocksfromFH"
 endif

 ! Read arguments into text strings.
 i=1
 temp=""
 call GETARG(i,temp)
 infile = ""
 infile = trim(temp)

 outfile=trim(adjustl(infile))//'.fh'

 np=0
 np(Hneut)=3
 np(Hplus)=1
 call allocate_sqcouplratio_parts(gsq)

 call createSLHAfilewithFHwithoutHBinputblocks(infile,outfile,&
     &          gsq(1)%hjbb_s,gsq(1)%hjbb_p,gsq(1)%hjtoptop_s,gsq(1)%hjtoptop_p, &
     &          gsq(1)%hjtautau_s,gsq(1)%hjtautau_p,                             &
     &          gsq(1)%hjWW,gsq(1)%hjZZ,                                         &
     &          gsq(1)%hjgg,gsq(1)%hjggZ,gsq(1)%hjhiZ                            )

 call addcouplingsblocktoSLHAfile(outfile,gsq(1)) 

end program HBSLHAinputblocksfromFH



