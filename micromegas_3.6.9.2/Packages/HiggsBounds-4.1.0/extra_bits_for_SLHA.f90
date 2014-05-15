! This file is part of HiggsBounds
!  -KW
!******************************************************************
module extra_bits_for_SLHA
!******************************************************************
 use PDGnumbering
 implicit none 
 integer,parameter :: h(5)=(/h0,HH,A0,h03,A02/)
  !h(1)=h0 !lightest CP even Higgs in real MSSM
  !h(2)=HH !heaviest CP even Higgs in real MSSM
  !h(3)=A0 !CP odd Higgs in real MSSM
  !h(4)=h03 !NMSSM only
  !h(5)=A02 !NMSSM only
 integer,parameter :: neutralino(5)=(/neut1,neut2,neut3,neut4,neut5/)
  !neutrlino(5) is in NMSSM only
 integer,parameter :: chargino(2)=(/char1p,char2p/)

 contains
      
 !************************************************************      
 subroutine getSLHAdata(d,gsq,infile)
 ! looks at theory predictions needed for Higgs searches only
 !************************************************************ 
  use usefulbits, only : dataset,hadroncolliderextras,Hneut,Hplus,Chineut,Chiplus,np,sqcouplratio,file_id_common
  use SLHA_manip
  
  implicit none
  !--------------------------------------input
  type(dataset) :: d
  type(sqcouplratio) :: gsq
  character(len=80),intent(in) :: infile
  !-----------------------------------internal
  integer :: i,j,x,ios,n
  integer :: particlecontent,Rparityviolation,CPviolation
  !double precision,allocatable ::g2hjcc(:,:),g2hjss(:,:)
  double precision,allocatable ::g2hjbb(:,:)
  double precision,allocatable ::g2hjtoptop(:,:),g2hjtautau(:,:)
  double precision :: mass
  integer :: visible_lspcandidate_PDG(21),invisible_lspcandidate_PDG(7)
  logical :: is_valid_point

  type lightestsusyparticle
   integer :: pdg
   double precision :: mass
  end type

  type(lightestsusyparticle) :: lsp

  logical :: invisible_lsp

  double precision :: cofmenergy
  !-------------------------------------------              
 
  if((np(Hneut).lt.0).or.(np(Hneut).gt.5))then
    stop'can not use subroutine getSLHAdata if number of neutral Higgs is not in range 0:5'
  endif
  if((np(Hplus).lt.0).or.(np(Hplus).gt.1))then
    stop'can not use subroutine getSLHAdata if number of charged Higgs is not in range 0:1'
  endif
  if((np(Chineut).lt.0).or.(np(Chineut).gt.5))then
    stop'can not use subroutine getSLHAdata if number of neutralinos is not in range 0:5'
  endif
  if((np(Chiplus).lt.0).or.(np(Chiplus).gt.2))then
    stop'can not use subroutine getSLHAdata if number of charginos is not in range 0:2'
  endif

  open(file_id_common,file=trim(infile),status='old',action='read',iostat=ios) 
  if(ios.ne.0)then 
   write(*,*)'problem opening the SLHA file: $'//trim(adjustl(infile))//'$'
   !stop'problem opening SLHA input file'
  else

   call readSLHAfile(file_id_common)

   call check_validity(is_valid_point)
   if(is_valid_point)then

    particlecontent  =get_modsel(3)
    Rparityviolation =get_modsel(4)
    CPviolation      =get_modsel(5)

    if(Rparityviolation.ne.0)stop'HB can not yet use SLHA files with R parity violation'

    select case(particlecontent)
    case(0)
     if((np(Hneut).ne.0).and.(np(Hneut).ne.3))then
      write(*,*)'warning: modsel(3)=0 (MSSM) in SLHA file'
      write(*,*)' but number of neutral Higgs that HB was told to expect=',np(Hneut)
     endif
     if((np(Chineut).ne.0).and.(np(Chineut).ne.4))then
      write(*,*)'warning: modsel(3)=0 (MSSM) in SLHA file'
      write(*,*)' but number of neutralinos that HB was told to expect=',np(Chineut)
     endif
    case(1)
     if((np(Hneut).ne.0).and.(np(Hneut).ne.5))then
      write(*,*)'warning: modsel(3)=1 (NMSSM) in SLHA file'
      write(*,*)' but number of neutral Higgs that HB was told to expect=',np(Hneut)
     endif
     if((np(Chineut).ne.0).and.(np(Chineut).ne.5))then
      write(*,*)'warning: modsel(3)=1 (NMSSM) in SLHA file'
      write(*,*)' but number of neutralinos that HB was told to expect=',np(Chineut)
     endif
    end select

    !------------------------------------------------
    !------------ work out what the lsp is -----------
    invisible_lsp=.False.

    call fill_visible_lspcandidate_PDG
    call fill_invisible_lspcandidate_PDG

    lsp%mass=1.0D12 !set to a very big value

    !find the invisible LSP candidate with the lowest mass
    !and use it for BR_hjinvisible
    do x=1,ubound(invisible_lspcandidate_PDG,dim=1)
       n=invisible_lspcandidate_PDG(x)
       mass=get_mass( n )
       if(mass.gt.0.0D0.and.(mass.lt.lsp%mass))then !mass<0 means that this particle mass couldn't be found
         lsp%mass=mass
         lsp%pdg =n
         invisible_lsp=.True.
       endif
    enddo

    !however, if there's a charged SUSY particle or gluino with
    !lower mass than the invisible LSP candidate, 
    !this candidate is not the LSP, so set BR_hjinvisible=0
    do x=1,ubound(visible_lspcandidate_PDG,dim=1)
       n=visible_lspcandidate_PDG(x)
       mass=get_mass( n )
       if(mass.gt.0.0D0.and.(mass.lt.lsp%mass))then !mass<0 means that this particle mass couldn't be found
         lsp%mass=mass
         lsp%pdg =n
         invisible_lsp=.False.
       endif
    enddo
    !------------------------------------------------

    if(np(Hneut).gt.0)then

     !allocate(g2hjcc(    np(Hneut),2))
     !allocate(g2hjss(    np(Hneut),2))
     allocate(g2hjbb(    np(Hneut),2))
     allocate(g2hjtoptop(np(Hneut),2))
     allocate(g2hjtautau(np(Hneut),2))

     do i=1,np(Hneut)
      d%particle(Hneut)%M(i) = get_mass(h(i))
      d%particle(Hneut)%dMh(i) = get_mass_uncertainty(h(i))
      if(d%particle(Hneut)%dMh(i).lt.0.1D0) d%particle(Hneut)%dMh(i)=0.0D0
      d%particle(Hneut)%dM(i) = d%particle(Hneut)%dMh(i)
      d%particle(Hneut)%GammaTot(i) = get_totaldecaywidth( h(i))

      d%BR_hjss(i) = get_twobodybranchingratio( h(i),squark,sbar)
      d%BR_hjcc(i) = get_twobodybranchingratio( h(i),cquark,cbar)
      d%BR_hjbb(i) = get_twobodybranchingratio( h(i),bquark,bbar)

      d%BR_hjmumu(i)   = get_twobodybranchingratio( h(i),mup,mum)
      d%BR_hjtautau(i) = get_twobodybranchingratio( h(i),taup,taum)

      d%BR_hjWW(i)   = get_twobodybranchingratio( h(i),Wp ,Wm)
      d%BR_hjZZ(i)   = get_twobodybranchingratio( h(i),Z0 ,Z0)
      d%BR_hjZga(i)  = get_twobodybranchingratio( h(i),Z0    ,photon)
      d%BR_hjgaga(i) = get_twobodybranchingratio( h(i),photon,photon)
      d%BR_hjgg(i)   = get_twobodybranchingratio( h(i),gluon ,gluon)

      if(invisible_lsp)then
         if(lsp%pdg.eq.neut1) then
         d%BR_hjinvisible(i) = get_twobodybranchingratio( h(i), lsp%pdg, lsp%pdg )
         else
         d%BR_hjinvisible(i) = get_twobodybranchingratio( h(i), lsp%pdg,-lsp%pdg )
         endif
      else
         d%BR_hjinvisible(i) = 0.0D0
      endif

      !g2hjcc(i,:)      = get_HiggsCouplingsFermions( h(i), cquark,  cquark   )
      !g2hjss(i,:)      = get_HiggsCouplingsFermions( h(i), squark,  squark   )
      g2hjbb(i,:)      = get_HiggsCouplingsFermions( h(i), bquark,  bquark   )
      g2hjtoptop(i,:)  = get_HiggsCouplingsFermions( h(i), tquark,  tquark   )
      g2hjtautau(i,:)  = get_HiggsCouplingsFermions( h(i), taum,    taum     )
      gsq%hjWW(i)      = get_HiggsCouplingsBosons( h(i), Wp,     Wp        )
      gsq%hjZZ(i)      = get_HiggsCouplingsBosons( h(i), Z0,     Z0        )
      gsq%hjgg(i)      = get_HiggsCouplingsBosons( h(i), gluon,  gluon     )   
      gsq%hjggZ(i)     = get_HiggsCouplingsBosons( h(i), gluon,  gluon, Z0 )
     enddo

     do j=1,np(Hneut)
      do i=1,np(Hneut)
        gsq%hjhiZ(j,i)   = get_HiggsCouplingsBosons( h(j),   h(i),   Z0    )
        d%BR_hjhihi(j,i) = get_twobodybranchingratio( h(j), h(i), h(i) )
      enddo
     enddo

     gsq%hjZga          = 0.0D0 !these are not needed
     gsq%hjgaga         = 0.0D0
     gsq%hjmumu_s       = 0.0D0
     gsq%hjmumu_p       = 0.0D0
     gsq%hjcc_s         = 0.0D0
     gsq%hjcc_p         = 0.0D0
     gsq%hjss_s         = 0.0D0
     gsq%hjss_p         = 0.0D0

    !gsq%hjcc_s(:)          = g2hjcc(:,1)
    !gsq%hjcc_p(:)          = g2hjcc(:,2)

    !gsq%hjss_s(:)          = g2hjss(:,1)
    !gsq%hjss_p(:)          = g2hjss(:,2)

     gsq%hjbb_s(:)          = g2hjbb(:,1)
     gsq%hjbb_p(:)          = g2hjbb(:,2)
 
     gsq%hjtoptop_s(:)      = g2hjtoptop(:,1)
     gsq%hjtoptop_p(:)      = g2hjtoptop(:,2)
 
     gsq%hjtautau_s(:)      = g2hjtautau(:,1)
     gsq%hjtautau_p(:)      = g2hjtautau(:,2)

     !deallocate(g2hjcc)
     !deallocate(g2hjss)
     deallocate(g2hjbb)
     deallocate(g2hjtoptop)
     deallocate(g2hjtautau)

    endif

    if(np(Hplus).gt.0)then
     i=1
     d%particle(Hplus)%M(i) = get_mass( Hp )
     d%particle(Hplus)%GammaTot(i) = get_totaldecaywidth( Hp )
     d%particle(Hplus)%dMh(i) = get_mass_uncertainty( Hp )
     !For now, set the LEP cross section ratio to one. Later, read it in from a block.
     d%lep%XS_HpjHmj_ratio(i)=1.0
     d%BR_tWpb        = get_twobodybranchingratio( tquark, Wp,   bquark )
     d%BR_tHpjb(i)    = get_twobodybranchingratio( tquark, Hp,   bquark )
     d%BR_Hpjcs(i)    = get_twobodybranchingratio( Hp,     cquark,  sbar   )
     d%BR_Hpjcb(i)    = get_twobodybranchingratio( Hp,     cquark,  bbar   )
     d%BR_Hpjtaunu(i) = get_twobodybranchingratio( Hp,     taup,    nutau  )
    endif


    if(lsp%pdg.eq.neut1)then ! all these chargino and neutralino searches rely on the neutralino1
                             ! being the lsp
     if(np(Chineut).gt.0)then  
       cofmenergy=get_SPhenocrosssectionCMenergy(ep,em,0.0D0,0.0D0,1)
       if(abs(cofmenergy-207.0D0).gt.1.0D-3)then
          write(*,*)'Wrong centre of mass energy for chargino and neutralino LEP production XS.'
       else 

          do i=1,np(Chineut)
            d%particle(Chineut)%M(i)        = abs(get_mass(neutralino(i))) ! 'abs' because 
                !SLHA files can have negative neutralino masses (if mixing matrix is defined to
                !be real). Note that this means that the neutralino masses should never 
                !be set to -1 to indicate they should be ignored (which we can do for the Higgs).
            d%particle(Chineut)%GammaTot(i) = get_totaldecaywidth(neutralino(i))  
          enddo

          d%lep%XS_NjNi=0.0D0
          do j=1,np(Chineut) !we're only interested when one of the particles is neutralino1
            d%lep%XS_NjNi(j,1)=get_SPhenocrosssection(neutralino(j),neut1)
            d%lep%XS_NjNi(1,j)=d%lep%XS_NjNi(j,1) !should be symmetric
          enddo 

          d%BR_NjqqNi  = 0.0D0 !ignoring for the moment

          d%BR_NjZNi   = 0.0D0
          d%BR_NjqqNi  = 0.0D0
          do j=2,np(Chineut)!we're only interested when one of the daughter particles is neutralino1
                            !and parent is not neutralino1
           d%BR_NjZNi(j,1)   = get_twobodybranchingratio(neutralino(j),neut1,Z0)
           d%BR_NjqqNi(j,1)  = get_threebodybranchingratio(neutralino(j),neut1,uquark,ubar)  &
                          &  + get_threebodybranchingratio(neutralino(j),neut1,dquark,dbar)  &
                          &  + get_threebodybranchingratio(neutralino(j),neut1,cquark,cbar)  &
                          &  + get_threebodybranchingratio(neutralino(j),neut1,squark,sbar)  &
                          &  + get_threebodybranchingratio(neutralino(j),neut1,tquark,tbar)  &
                          &  + get_threebodybranchingratio(neutralino(j),neut1,bquark,ubar)  
         enddo

         if(np(Chiplus).gt.0)then

            do i=1,np(Chiplus)
             d%particle(Chiplus)%M(i)        = abs(get_mass(chargino(i)))
             d%particle(Chiplus)%GammaTot(i) = get_totaldecaywidth(chargino(i))
            enddo

            do j=1,np(Chiplus)
              d%lep%XS_CpjCmj(j)=get_SPhenocrosssection(chargino(j),-chargino(j))
            enddo

            d%BR_CjWNi   = 0.0D0
            d%BR_CjqqNi  = 0.0D0 
            d%BR_CjlnuNi = 0.0D0
            do j=1,np(Chiplus)
               d%BR_CjWNi(j,1)   = get_twobodybranchingratio(chargino(j),neut1,Wp)
               d%BR_CjqqNi(j,1)  = get_threebodybranchingratio(chargino(j),neut1,dbar,uquark) &
                               & + get_threebodybranchingratio(chargino(j),neut1,sbar,cquark) 
               d%BR_CjlnuNi(j,1) = get_threebodybranchingratio(chargino(j),neut1,ep,  nue)  &
                               & + get_threebodybranchingratio(chargino(j),neut1,mup, numu) &
                               & + get_threebodybranchingratio(chargino(j),neut1,taup,nutau) 
            enddo

         endif
       endif
     endif
    endif



   endif
   close(file_id_common)

   call finishwithSLHA
  endif

  contains

  subroutine fill_invisible_lspcandidate_PDG

   invisible_lspcandidate_PDG(1) = neut1

   invisible_lspcandidate_PDG(2) = s_nueL   
   invisible_lspcandidate_PDG(3) = s_numuL  
   invisible_lspcandidate_PDG(4) = s_nutauL 

   invisible_lspcandidate_PDG(5) = s_nutau1A
   invisible_lspcandidate_PDG(6) = s_nutau2A
   invisible_lspcandidate_PDG(7) = s_nutau3A

  end subroutine fill_invisible_lspcandidate_PDG

  subroutine fill_visible_lspcandidate_PDG
  
   visible_lspcandidate_PDG(1) = char1p
   visible_lspcandidate_PDG(2) = char2p
 
   visible_lspcandidate_PDG(3) = s_t1 
   visible_lspcandidate_PDG(4) = s_b1 
   visible_lspcandidate_PDG(5) = s_cL 
   visible_lspcandidate_PDG(6) = s_sL 
   visible_lspcandidate_PDG(7) = s_uL 
   visible_lspcandidate_PDG(8) = s_dL 

   visible_lspcandidate_PDG(9) = s_t2 
   visible_lspcandidate_PDG(10) = s_b2 
   visible_lspcandidate_PDG(11) = s_cR 
   visible_lspcandidate_PDG(12) = s_sR 
   visible_lspcandidate_PDG(13) = s_uR 
   visible_lspcandidate_PDG(14) = s_dR 

   visible_lspcandidate_PDG(15) = s_emL   
   visible_lspcandidate_PDG(16) = s_mumL  
   visible_lspcandidate_PDG(17) = s_taumL  

   visible_lspcandidate_PDG(18) = s_emR    
   visible_lspcandidate_PDG(19) = s_mumR  
   visible_lspcandidate_PDG(20) = s_taumR 

   visible_lspcandidate_PDG(21) = gluino

  end subroutine fill_visible_lspcandidate_PDG


 end subroutine getSLHAdata
 !************************************************************     
 subroutine outputSLHAdata(infile) 
 !************************************************************ 
  use usefulbits, only : whichanalyses,pr,vers,file_id_common,fullHBres,infile1
  use SLHA_manip
  use S95tables      
  
  implicit none
  !--------------------------------------input
  character(len=80),intent(in) :: infile
  !-----------------------------------internal
  integer :: x,y,ios,jj
  integer :: k_out
  character(LEN=200):: descrip 
  logical :: newfile = .False.
  logical :: exist
  !-------------------------------------------              

  if(infile.eq.infile1) then
   open(file_id_common,file=trim(infile),status='old',iostat=ios) 
  else
   inquire(file=trim(infile),exist=exist)
   if(exist) then
    open(file_id_common,file=trim(infile),status='replace',iostat=ios)
   else
    open(file_id_common,file=trim(infile),status='new',iostat=ios) 
   endif 
   newfile=.True.   
  endif 
  
   if(ios.ne.0)then 
     write(*,*)'problem opening the SLHA file: $'//trim(adjustl(infile))//'$'
   else
    k_out=file_id_common

    if(.not.newfile) then
     call readSLHAfile(file_id_common)
     rewind(file_id_common)
     call writeSLHAfile_except(k_out,'HiggsBoundsResults')
    endif 

    !write(k_out,'(a)')'#'
    write(k_out,'(a)')'Block HiggsBoundsResults      # results from HiggsBounds http://projects.hepforge.org/higgsbounds'
    write(k_out,'(a)')'# HBresult   : scenario allowed flag (1: allowed, 0: excluded, -1: unphysical)'      
    write(k_out,'(a)')'# chan id number: most sensitive channel (see below). chan=0 if no channel applies'          
    write(k_out,'(a)')'# obsratio   : ratio [sig x BR]_model/[sig x BR]_limit (<1: allowed, >1: excluded)' 
    write(k_out,'(a)')'# ncomb      : number of Higgs bosons combined in most sensitive channel'
    write(k_out,'(a)')'# Note that the HB channel id number varies depending on the HB version and setting "whichanalyses"'
    write(k_out,'(a)')'#'
    write(k_out,*)'   0    '//trim(adjustl(vers))//'     ||'//whichanalyses// &
     &  '||            # version of HB used to produce these results,the HB setting "whichanalyses"'
    write(k_out,'(a)')'#'

    jj=1
    y=1
!    x=res(jj)%chan(y)
    x=fullHBres(y)%chan
    call outputproc(pr(x),0,descrip,1)

    write(k_out,'(a)')'#CHANNELTYPE 1: channel with the highest statistical sensitivity'
!    write(k_out,*)'   1',1   ,    res(jj)%chan(y),     '                 # channel id number'
!    write(k_out,*)'   1',2   ,    res(jj)%allowed95(y),'                 # HBresult      '
!    write(k_out,*)'   1',3   ,    res(jj)%obsratio(y), '   # obsratio  '
!    write(k_out,*)'   1',4   ,    res(jj)%ncombined(y),'                 # ncombined'      
    write(k_out,*)'   1',1   ,    fullHBres(y)%chan,     '                 # channel id number'
    write(k_out,*)'   1',2   ,    fullHBres(y)%allowed95,'                 # HBresult      '
    write(k_out,*)'   1',3   ,    fullHBres(y)%obsratio, '   # obsratio  '
    write(k_out,*)'   1',4   ,    fullHBres(y)%ncombined,'                 # ncombined'      
    write(k_out,*)'   1',5   ,'||'//trim(adjustl(descrip))//'|| # text description of channel'
    write(k_out,'(a)')'#'

    close(file_id_common)
    if(k_out.ne.file_id_common)close(k_out)

    call finishwithSLHA

  endif
 end subroutine outputSLHAdata

 !************************************************************     
 subroutine addcouplingsblocktoSLHAfile(infile,gsq) 
 !************************************************************
  use usefulbits, only : file_id_common,sqcouplratio,np,Hneut
  use SLHA_manip 
  
  implicit none
  !--------------------------------------input
  type(sqcouplratio) :: gsq
  character(len=300),intent(in) :: infile
  !-----------------------------------internal
  integer :: i,ios,j
  integer :: k_out
  !-------------------------------------------              

  open(file_id_common,file=trim(infile),status='old',iostat=ios) 
  if(ios.ne.0)then 
    write(*,*)'problem opening the SLHA file: $'//trim(adjustl(infile))//'$'
  else

    call readSLHAfile(file_id_common) 

    rewind(file_id_common)
    k_out=file_id_common

    call writeSLHAfile_except(k_out, &
                   &  'HiggsBoundsInputHiggsCouplingsBosons', &
                   &  'HiggsBoundsInputHiggsCouplingsFermions')

    write(k_out,'(a)')'#'
    write(k_out,'(a)')'Block HiggsBoundsInputHiggsCouplingsBosons'
    write(k_out,'(a)')'# For exact definitions of NormEffCoupSq see HiggsBounds manual'
  
    do i=1,np(Hneut)
       write(k_out,'(G16.6,4I6,a)')gsq%hjWW(i),        3,h(i),Wp,Wp, &
           &  ' # higgs-W-W effective coupling^2, normalised to SM'
    enddo
    do i=1,np(Hneut) 
       write(k_out,'(G16.6,4I6,a)')gsq%hjZZ(i),        3,h(i),Z0,Z0, &
           &  ' # higgs-Z-Z effective coupling^2, normalised to SM'
    enddo
    do i=1,np(Hneut) 
       write(k_out,'(G16.6,4I6,a)')gsq%hjgg(i),        3,h(i),gluon,gluon, &
           &  ' # higgs-gluon-gluon effective coupling^2, normalised to SM'
    enddo
    do j=1,np(Hneut)
      do i=1,j   
        write(k_out,'(G16.6,4I6,a)')gsq%hjhiZ(j,i),3,h(j),h(i),Z0, &
           &  ' # higgs-higgs-Z effective coupling^2, normalised'
      enddo
    enddo
    do i=1,np(Hneut) 
       write(k_out,'(G16.6,5I6,a)')gsq%hjggZ(i),        4,h(i),gluon,gluon,Z0, &
           &  ' # higgs-gluon-gluon-Z effective coupling^2, normalised to SM'
    enddo

    write(k_out,'(a)')'#'
    write(k_out,'(a)')'Block HiggsBoundsInputHiggsCouplingsFermions'
    write(k_out,'(a)')'# For exact definitions of NormEffCoupSq see HiggsBounds manual'
    write(k_out,'(a)')'# ScalarNormEffCoupSq PseudoSNormEffCoupSq    NP    IP1      IP2      IP3' // &
           &   ' # Scalar, Pseudoscalar Normalised Effective Coupling Squared'
    do i=1,np(Hneut) 
       write(k_out,*)gsq%hjbb_s(i),gsq%hjbb_p(i),        3,h(i),bquark,bquark, &
           &  '# higgs-b-b eff. coupling^2, normalised to SM'
    enddo
    do i=1,np(Hneut) 
       write(k_out,*)gsq%hjtoptop_s(i),gsq%hjtoptop_p(i),3,h(i),tquark,tquark, &
           &  '# higgs-top-top eff. coupling^2, normalised to SM'
    enddo
    do i=1,np(Hneut) 
       write(k_out,*)gsq%hjtautau_s(i),gsq%hjtautau_p(i),3,h(i),taum,taum, &
           &   '# higgs-tau-tau eff. coupling^2, normalised to SM'
    enddo

    close(file_id_common)
    if(k_out.ne.file_id_common)close(k_out)

    call finishwithSLHA

  endif
 end subroutine addcouplingsblocktoSLHAfile
 !************************************************************ 
end module extra_bits_for_SLHA
!******************************************************************
