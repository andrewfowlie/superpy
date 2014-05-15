c ==================================================================== c
c                        slepton/sneutrino decays                      c
c ==================================================================== c
c These routines are generalisations of corresponding routines from
c     SDECAY: A Fortran code for the decays of the supersymmetric 
c         particles in the MSSM
c     by M. Muhlleitner (Karlsruhe, Inst. Technol.),
c        A. Djouadi (Orsay, LPT & CERN, Theory Division),
c        Y. Mambrini (Orsay, LPT),
c     Comput.Phys.Commun.168:46-70 (2005), hep-ph/0311167.
c    SDECAY should be cited whenever NMSDECAY is used.
c
c
      SUBROUTINE NS_SLEPTON
*
      IMPLICIT NONE
      INTEGER I,J
      DOUBLE PRECISION NS_lamb,PI,SQR2,multilim
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS      
      DOUBLE PRECISION sw,cw,tw,flagmulti,flagqcd,flagloop
      DOUBLE PRECISION ale(2,2),alto(2,2),alsne(2,2),blsne(2,2),
     .          alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2)
      DOUBLE PRECISION gcsntaur(2,2)
      DOUBLE PRECISION Hstaustaur(3,2,2),Astaustaur(2,2,2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION asne2,asntau2
      DOUBLE PRECISION sellneute(5),selrneute(5),
     .     sellcharnue(2),selrcharnue(2),snellneut(5),snellchar(2),
     .     sntauneut(5),sntauchar(2),sntau1hcstau(2),sntau1wstau(2),
     .     stau1neut(5),stau1char(2),stau2neut(5),
     .     stau2char(2),stau1hcsn(2),stau2hcsn(2),
     .     stau2H(3),stau2A(2),stau2ztau,stau1wsn(2),stau2wsn(2)
      DOUBLE PRECISION gmsntau(2),gmstau(2)
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .          brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .          brsntauneut(5),brsntauchar(2),brsntau1hcstau(2),
     .          brsntau1wstau(2),brstau1neut(5),brstau2neut(5),
     .          brstau1char(2),
     .          brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .          brstau1wsn(2),brstau2wsn(2),brstau2H(3),brstau2A(2),
     .          brstau2ztau
      DOUBLE PRECISION xintegslstaustarltau,xintegslstaultau
      DOUBLE PRECISION brselrstau,brselrstaustar,brstau2stau1star,
     .    brstau2stau1,brstau2stau1nn,
     .    brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .    brsellstau1star,brsellstau1,brsellstau1nutau,
     .    brsnestau1star,brsnestau1,brsnestau1nutau
      DOUBLE PRECISION xintegstau2stau1star,xintegstau2stau1,
     .    xintegstau2stau1nn
      DOUBLE PRECISION xintegsntaustau1star,xintegsntaustau1,
     .    xintegsntaustau1nutau
      DOUBLE PRECISION xintegsellstau1star,xintegsellstau1,
     .    xintegsellstau1nutau
       DOUBLE PRECISION xintegsnestau1star,xintegsnestau1,
     .    xintegsnestau1nutau
       DOUBLE PRECISION selltot,selltot2,selltot3,
     .selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     .sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
*     
      COMMON/SLEPTON_2GAMMA/sellneute,selrneute,
     .    sellcharnue,selrcharnue,snellneut,snellchar,
     .    sntauneut,sntauchar,sntau1hcstau,sntau1wstau,
     .    stau1neut,stau1char,stau2neut,stau2char,stau1hcsn,
     .    stau2hcsn,stau2H,stau2A,stau2ztau,stau1wsn,stau2wsn
      COMMON/SLEPTON_WIDTH/selltot,selltot2,selltot3,
     .selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     .sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
      COMMON/SLEPTON_BR_2BD/brsellneute,brselrneute,brsellcharnue,
     .          brselrcharnue,brsnellneut,brsnellchar,
     .          brsntauneut,brsntauchar,brsntau1hcstau,
     .          brsntau1wstau,brstau1neut,brstau2neut,brstau1char,
     .          brstau2char,brstau1hcsn,brstau2hcsn,
     .          brstau1wsn,brstau2wsn,brstau2H,brstau2A,brstau2ztau
      COMMON/SLEPTON_BR_3BD/brselrstau,brselrstaustar,brstau2stau1star,
     .   brstau2stau1,brstau2stau1nn,
     .   brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .    brsellstau1star,brsellstau1,brsellstau1nutau,
     .    brsnestau1star,brsnestau1,brsnestau1nutau
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_rhsneutr/asne2,asntau2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup5/ale,alto,alsne,blsne,alsnt,blsnt
      COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/NS_coup20/gwtb,gwntau
      COMMON/NS_HIGGSTAUTAU/Hstaustaur,Astaustaur
      COMMON/NS_hcstausntau/gcsntaur
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_multilim/multilim
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop
*    
      EXTERNAL NS_lamb

      selltot  = 0.d0
      selltot2 = 0.D0
      selltot3 = 0.d0
      selrtot2 = 0.D0
      selrtot = 0.d0
      selrtot3 = 0.d0
      
      do i =1,2
         sellcharnue(i) = 0.D0
         selrcharnue(i) = 0.D0
      enddo
      do j=1,5
        sellneute(j) = 0.D0
        selrneute(j) = 0.D0
      enddo

   
      snellchar(1) = 0.D0
      snellchar(2) = 0.D0
      do j=1,5
         snellneut(j) = 0.D0
      enddo

      stau1tot2 = 0.D0 
      stau2tot2 = 0.D0
      stau2tot3 = 0.d0
      stau2tot = 0.d0
      do j=1,5
         stau1neut(j) = 0.D0
         stau2neut(j) = 0.D0
      enddo
      do i=1,2
         stau1char(i) = 0.D0
         stau1hcsn(i) = 0.D0
         stau1wsn(i)  = 0.D0
         stau2char(i) = 0.D0
         stau2hcsn(i) = 0.D0
         stau2wsn(i)  = 0.D0
      end do
      do i=1,3
      stau2H(i) = 0.D0
      enddo
      do i=1,2
      stau2A(i) = 0.D0
      enddo
      stau2ztau = 0.D0
      sneltot2 =0.d0
      sneltot3 =0.d0
      sneltot =0.d0
      sntautot2=0.d0
      sntautot3=0.d0
      sntautot=0.d0
      do i=1,2
         sntauchar(i)    = 0.D0
         sntau1hcstau(i) = 0.D0
         sntau1wstau(i)  = 0.D0
      end do
      do j=1,5
         sntauneut(j) = 0.D0
      enddo
c ==================================================================== c
c                        selectron 2-body decays                       c
c ==================================================================== c

c  selectron_L --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + e-

      do i=1,5
         if(ase1.gt.amneut(i)) then
            sellneute(i)=g2s*ae(1,i)**2*
     .           (ase1**2-amneut(i)**2)*
     .           NS_lamb(0.D0,amneut(i)/ase1)/(16*pi*ase1)
         else
            sellneute(i)=0.D0
         endif
      enddo
      
C      WRITE(*,*) "sellneute(1,2)",sellneute(1),sellneute(2),
C     .sellneute(3),sellneute(4),sellneute(5)  
c -------------------------------------------------------------------- c
c  selectron_L --> chi-_1/chi-_2 + neutrino_e

      do i=1,2
         if(ase1.gt.amchar(i)) then
            sellcharnue(i)=g2s*ale(1,i)**2*(ase1**2-amchar(i)**2)*
     .                     NS_lamb(0.D0,amchar(i)/ase1)/(16*pi*ase1)
         else
            sellcharnue(i)=0.D0
         endif
      enddo
      
C      WRITE(*,*) "sellcharnue(1,2)",sellcharnue(1),sellcharnue(2)
c -------------------------------------------------------------------- c
c  selectron_R --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + e-

      do i=1,5
         if(ase2.gt.amneut(i)) then
            selrneute(i)=g2s*be(2,i)**2*
     .           (ase2**2-amneut(i)**2)*
     .           NS_lamb(0.D0,amneut(i)/ase2)/(16*pi*ase2)
         else
            selrneute(i)=0.D0
         endif
      enddo
      
C      WRITE(*,*) "selrneute(1,2)",selrneute(1),selrneute(2),
C     .selrneute(3),selrneute(4),selrneute(5)
c -------------------------------------------------------------------- c
c  selectron_R --> chi-_1/chi-_2 + neutrino_e
C  UE: vanishes since ale(2,i)=0

      do i=1,2
            selrcharnue(i)=0.D0
      enddo

c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

        do i=1,5
          selltot2=selltot2+sellneute(i)
          selrtot2=selrtot2+selrneute(i)
        enddo

        do i=1,2
          selltot2=selltot2+sellcharnue(i)
          selrtot2=selrtot2+selrcharnue(i)
        enddo
c--------------------------------------------------------------------- c
c ----- selectron_l 3-body decays and 3-body total widths ------------ c
c--------------------------------------------------------------------- c
c Only: selectron_l -> stau_1* + e + tau      (xintegsellstau1star)
c       selectron_l -> stau_1 + e + tau       (xintegsellstau1)
c       selectron_l -> stau_1 + nu_e + nu_tau (xintegsellstau1nutau)

        Call NS_xintegsell(xintegsellstau1star,xintegsellstau1,
     .    xintegsellstau1nutau)
        
        selltot3=xintegsellstau1star+xintegsellstau1
     .    +xintegsellstau1nutau

        if(selltot3.lt.multilim*selltot2)Then
           selltot3=0.d0
        endif
        selltot=selltot2+selltot3
      do i=1,5
         brsellneute(i)=sellneute(i)/selltot
      enddo
      do i=1,2
         brsellcharnue(i)=sellcharnue(i)/selltot
      enddo
      if (selltot3.ne.0.D0)Then
      brsellstau1star=xintegsellstau1star/selltot
      brsellstau1=xintegsellstau1/selltot
      brsellstau1nutau=xintegsellstau1nutau/selltot
      endif
c      write(*,*)"Selectron_l_3-body BRs to stau1star,stau1,stau1+neut:"
c      write(*,*) brsellstau1star,brsellstau1,brsellstau1nutau
c--------------------------------------------------------------------- c
c ----- RH selectron 3-body decays and 3-body total widths ----------- c
c--------------------------------------------------------------------- c
c Only:  selectron_r -> stau_1* + e + tau* (xintegslstaustarltau)
c        selectron_r -> stau_1 + e + tau   (xintegslstaultau)
      if(flagmulti.eq.1.D0) then
        Call NS_xintegselr(xintegslstaustarltau,xintegslstaultau)

        selrtot3=xintegslstaustarltau+xintegslstaultau

        if (selrtot3.lt.multilim*selrtot2)Then
           selrtot3=0.d0
        endif
        selrtot=selrtot2+selrtot3
        
      do i=1,5
         brselrneute(i) = selrneute(i)/selrtot
      enddo
      do i=1,2
         brselrcharnue(i) = selrcharnue(i)/selrtot
      enddo
      if(selrtot3.ne.0.D0) then
      brselrstaustar=xintegslstaustarltau/selrtot
      brselrstau=xintegslstaultau/selrtot
      endif
      endif
C     WRITE(*,*) " brsellneute(1,2)",  brsellneute(1), brsellneute(2)
C     WRITE(*,*) " brselrneute(1,2)",  brselrneute(1), brselrneute(2)
C     WRITE(*,*) " brsellcharnue(1,2)",brsellcharnue(1),brsellcharnue(2)
C     WRITE(*,*) " brselrcharnue(1,2)",brselrcharnue(1),brselrcharnue(2)
c     write(*,*) "Selectron_r_3-body BRs to stau1star,stau1:"
c     WRITE(*,*)  brselrstaustar,brselrstau
c ==================================================================== c
c                       sneutrino_el 2-body decays                     c
c ==================================================================== c
c  sneutrino_eL --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + neutrino_e

      do i=1,5
         if(asne1.gt.amneut(i)) then
            snellneut(i)=g2s*(anu(1,i)**2+bnu(2,i)**2)*
     .        (asne1**2-amneut(i)**2)*NS_lamb(0.D0,amneut(i)/asne1)
     .        /(16*pi*asne1)
         else
            snellneut(i)=0.D0
         endif
      enddo
C      WRITE(*,*) "snellneut(1,2)",snellneut(1),snellneut(2),
C     .snellneut(3),snellneut(4),snellneut(5)
c -------------------------------------------------------------------- c
c  sneutrino_eL --> chi+_1/chi+_2 + e-

      do i=1,2
         if(asne1.gt.amchar(i)) then
           snellchar(i)=g2s*alsne(1,i)**2*(asne1**2-amchar(i)**2)*
     .           NS_lamb(0.D0,amchar(i)/asne1)/(16*pi*asne1)
         else
            snellchar(i)=0.D0
         endif
      enddo

C      WRITE(*,*) "snellchar(1,2)",snellchar(1),snellchar(2)
c -------------------------------------------------------------------- c
      sneltot2=snellneut(1)+snellneut(2)+snellneut(3)+
     .         snellneut(4)+snellneut(5)+snellchar(1)+snellchar(2)
c--------------------------------------------------------------------- c
c -------- sneutrino_e 3-body decays and 3-body total widths --------- c
c--------------------------------------------------------------------- c
c Only: sne -> nu_e + stau1* + tau* (xintegsnestau1star)
c       sne -> nu_e + stau1 + tau   (xintegsnestau1)
c       sne -> e + stau1 + nu_tau   (xintegsnestau1nutau)
      if(flagmulti.eq.1.D0) then
        Call NS_xintegsne(xintegsnestau1star,xintegsnestau1,
     .    xintegsnestau1nutau)
        
        sneltot3=xintegsnestau1star+xintegsnestau1+
     .    xintegsnestau1nutau

        if (sneltot3.lt.multilim*sneltot2)Then
           sneltot3=0.d0
        endif

        sneltot=sneltot2+sneltot3

      do i=1,5
         brsnellneut(i) = snellneut(i)/sneltot
      enddo
      do i=1,2
         brsnellchar(i) = snellchar(i)/sneltot
      enddo
      If (sneltot3.ne.0d0)Then
      brsnestau1star=xintegsnestau1star/sneltot
      brsnestau1=xintegsnestau1/sneltot
      brsnestau1nutau=xintegsnestau1nutau/sneltot
      endif
      endif

C      WRITE(*,*) "brsnellneut(1,2)",brsnellneut(1),brsnellneut(2)
C      WRITE(*,*) "brsnellchar(1,2)",brsnellchar(1),brsnellchar(2)
c      write(*,*)"Sneutrino_e_3-body BRs to stau1star,stau1,stau1+neut:"
c      write(*,*) brsnestau1star,brsnestau1,brsnestau1nutau,sneltot,
c     .sneltot2,xintegsnestau1star,xintegsnestau1,
c     .    xintegsnestau1nutau
c ==================================================================== c
c                       stau1/2 2-body decays                          c
c ==================================================================== c

      gmsntau(1) = asntau1
      gmsntau(2) = asntau2
c -------------------------------------------------------------------- c
c  stau1 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + tau

      do i=1,5
         if(astau1.gt.(amneut(i)+mtau)) then
            stau1neut(i)=g2s*((atau(1,i)**2+btau(1,i)**2)*
     .           (astau1**2-amneut(i)**2-mtau**2)
     .           -4*atau(1,i)*btau(1,i)*mtau*xmneut(i)
     .           )*NS_lamb(mtau/astau1,amneut(i)/astau1)
     .        /(16.D0*pi*astau1)
         else
            stau1neut(i)=0.D0
         endif
      enddo
C      WRITE(*,*) "stau1neut(1,2)",stau1neut(1),stau1neut(2),
C     .stau1neut(3),stau1neut(4),stau1neut(5)
c -------------------------------------------------------------------- c
c  stau2 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + tau
      
      do i=1,5
         if(astau2.gt.(amneut(i)+mtau)) then
            stau2neut(i)=g2s*((atau(2,i)**2+btau(2,i)**2)*
     .           (astau2**2-amneut(i)**2-mtau**2)
     .           -4*atau(2,i)*btau(2,i)*mtau*xmneut(i)
     .           )*NS_lamb(mtau/astau2,amneut(i)/astau2)
     .           /(16.D0*pi*astau2)
         else
            stau2neut(i)=0.D0
         endif
      enddo
C      WRITE(*,*) "stau2neut(1,2)",stau2neut(1),stau2neut(2),
C     .stau2neut(3),stau2neut(4),stau2neut(5)
c ----------------------------------------------------------------- c
c  stau1 --> chi+_1/chi+_2 + nu_tau
      do i=1,2
         if(astau1.gt.amchar(i)) then
             stau1char(i)=g2s*alto(1,i)**2*(astau1**2-amchar(i)**2)*
     .        NS_lamb(0.D0,amchar(i)/astau1)/(16*pi*astau1)
         else
            stau1char(i)=0.D0
         endif
      enddo

C      WRITE(*,*) "stau1char(1)",stau1char(1),stau1char(2)
c -------------------------------------------------------------------- c
c  stau2 --> chi-_1/chi-_2 + nu_tau
      do i=1,2
         if(astau2.gt.amchar(i)) then
            stau2char(i)=g2s*alto(2,i)**2*(astau2**2-amchar(i)**2)*
     .        NS_lamb(0.D0,amchar(i)/astau2)/(16*pi*astau2)
         else
            stau2char(i)=0.D0
         endif
      enddo

C      WRITE(*,*) "stau2char(1)",stau2char(1),stau2char(2)
c -------------------------------------------------------------------- c
c  stau1 --> H- + sneutrino_tau1/2
      do i=1,2
         if(astau1.gt.(gmsntau(i)+cmass)) then
            stau1hcsn(i)=g2s*mw**2*gcsntaur(i,1)**2*
     .           NS_lamb(gmsntau(i)/astau1,cmass/astau1)
     .           /(16.D0*pi*astau1)
         else
            stau1hcsn(i)=0.D0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau2 --> H- + sneutrino_tau1/2
      do i=1,2
         if(astau2.gt.(gmsntau(i)+cmass)) then
            stau2hcsn(i)=g2s*mw**2*gcsntaur(i,2)**2*
     .           NS_lamb(gmsntau(i)/astau2,cmass/astau2)
     .           /(16.D0*pi*astau2)
         else
            stau2hcsn(i)=0.D0
         endif
      enddo
C      WRITE(*,*) "gcsntaur(1,2)",gcsntaur(1,2)
C      WRITE(*,*) "stau2hcsn(1,2)",stau2hcsn(1),stau2hcsn(2)
c -------------------------------------------------------------------- c
c  stau2 --> h(i) + stau1 [stau2H(3)]
      do I=1,3
      if(astau2.gt.(astau1+SMASS(I))) then
         stau2H(i)=g2s*mz**4/mw**2*Hstaustaur(I,2,1)**2*
     .         NS_lamb(astau1/astau2,SMASS(I)/astau2)/(16.D0*pi*astau2)
      else
         stau2H(I)=0.D0
      endif
      enddo

c      WRITE(*,*) "stau2H(1,2)",stau2H(1),stau2H(2)
c -------------------------------------------------------------------- c
c  stau2 --> A(i) + stau1 [stau2A(2)]
      do I=1,2
      if(astau2.gt.(astau1+PMASS(I))) then
            stau2A(I)=g2s*mz**4/mw**2*Astaustaur(I,2,1)**2*
     .         NS_lamb(astau1/astau2,PMASS(I)/astau2)/(16.D0*pi*astau2)
         else
           stau2A(I)=0.D0
         endif
      enddo

C      WRITE(*,*) PMASS(1),Astaustaur(1,2,1)
C      WRITE(*,*) "stau2A(1,2)",stau2A(1),stau2A(2)

c -------------------------------------------------------------------- c
c  stau2 --> Z + stau1

      if(astau2.gt.(astau1+mz)) then
        stau2ztau=g2s/64.D0/pi/cw**2/mz**2*astau2**3*gztautau(2,1)**2*
     .           NS_lamb(astau1/astau2,mz/astau2)**3
      else
         stau2ztau=0.D0
      endif

C      WRITE(*,*) "stau2ztau",stau2ztau
c -------------------------------------------------------------------- c
c  stau1 --> W- + sntau1/2
      do i=1,2
         if(astau1.gt.(gmsntau(i)+mw)) then
            stau1wsn(i)=g2s/32.D0/pi/mw**2*astau1**3*gwntau(i,1)**2*
     .                NS_lamb(gmsntau(i)/astau1,mw/astau1)**3
         else
            stau1wsn(i)=0.D0
         endif
      enddo

C      WRITE(*,*) "stau1wsn(1)",stau1wsn(1)
c -------------------------------------------------------------------- c
c  stau2 --> W- + sntau1/2
      do i=1,2
         if(astau2.gt.(gmsntau(i)+mw)) then
            stau2wsn(i)=g2s/32.D0/pi/mw**2*astau2**3*gwntau(i,2)**2*
     .                NS_lamb(gmsntau(i)/astau2,mw/astau2)**3
         else
            stau2wsn(i)=0.D0
         endif
      enddo

C      WRITE(*,*) "stau2wsn(1)",stau2wsn(1)
c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

      stau1tot2=stau1neut(1)+stau1neut(2)+stau1neut(3)+stau1neut(4)+
     .          stau1neut(5)+
     .          stau1char(1)+stau1char(2)+stau1hcsn(1)+stau1hcsn(2)+
     .          stau1wsn(1)+stau1wsn(2)


      stau2tot2=stau2neut(1)+stau2neut(2)+stau2neut(3)+stau2neut(4)+
     .          stau2neut(5)+
     .          stau2char(1)+stau2char(2)+stau2hcsn(1)+stau2hcsn(2)+
     .          stau2wsn(1)+stau2wsn(2)+stau2H(1)+stau2H(2)+stau2H(3)+
     .          stau2A(1)+stau2A(2)+stau2ztau

      do i=1,5
         brstau1neut(i)=stau1neut(i)/stau1tot2
      enddo
      do i=1,2
         brstau1char(i)=stau1char(i)/stau1tot2
         brstau1hcsn(i)=stau1hcsn(i)/stau1tot2
         brstau1wsn(i) =stau1wsn(i)/stau1tot2
      enddo

c--------------------------------------------------------------------- c
c ----- stau_2 3-body decays and 3-body total widths ----------------- c
c--------------------------------------------------------------------- c
c Only: stau_2 -> stau_1* + tau + tau*     (xintegstau2stau1star)
c       stau_2 -> stau_1 + tau + tau       (xintegstau2stau1)
c       stau_2 -> stau_1 + nu_tau + nu_tau (xintegstau2stau1nn)
      if(flagmulti.eq.1.D0) then
        Call NS_xintegstau2(xintegstau2stau1star,xintegstau2stau1,
     .    xintegstau2stau1nn)
        
        stau2tot3=xintegstau2stau1star+xintegstau2stau1
     .    +xintegstau2stau1nn

        if (stau2tot3.lt.multilim*stau2tot2)Then
           stau2tot3=0.d0
        endif
           
        stau2tot=stau2tot2+stau2tot3

      do i=1,5
         brstau2neut(i)=stau2neut(i)/stau2tot
      enddo
      do i=1,2
         brstau2char(i)=stau2char(i)/stau2tot
         brstau2hcsn(i)=stau2hcsn(i)/stau2tot
         brstau2wsn(i) =stau2wsn(i)/stau2tot
      enddo
      do I=1,3
         brstau2H(I) = stau2H(I)/stau2tot
      enddo
      do I=1,2
         brstau2A(I) = stau2A(I)/stau2tot
      enddo
      brstau2ztau = stau2ztau/stau2tot
      if (stau2tot3.ne.0d0)Then
      brstau2stau1star=xintegstau2stau1star/stau2tot
      brstau2stau1=xintegstau2stau1/stau2tot
      brstau2stau1nn=xintegstau2stau1nn/stau2tot
      endif
      endif
c      write(*,*) "Stau2_3-body BRs to stau1star,stau1,stau1+neut:"
c      write(*,*) brstau2stau1star,brstau2stau1,brstau2stau1nn

c ==================================================================== c
c                       sneutrino_tau 2-body decays                    c
c ==================================================================== c

      gmstau(1) = astau1
      gmstau(2) = astau2
c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + neutrino_tau
      do i=1,5
         if(asntau1.gt.amneut(i)) then
            sntauneut(i)=g2s*(antau(1,i)**2+bntau(1,i)**2)*
     .        (asntau1**2-amneut(i)**2)*
     .           NS_lamb(0.D0,amneut(i)/asntau1)
     .        /(16*pi*asntau1)
         else
            sntauneut(i)=0.D0
         endif
      enddo

C      write(*,*)'sntauneut(i)',sntauneut(1),sntauneut(2)
c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> chi+_1/chi+_2 + tau-
      do i=1,2
         if(asntau1.gt.(amchar(i)+mtau)) then
            sntauchar(i)=g2s*((alsnt(1,i)**2+blsnt(1,i)**2)*
     .           (asntau1**2-amchar(i)**2-mtau**2)
     .           -4.D0*alsnt(1,i)*blsnt(1,i)*xmchar(i)*mtau)*
     .           NS_lamb(mtau/asntau1,amchar(i)/asntau1)
     .           /(16*pi*asntau1)
         else
            sntauchar(i)=0.D0
         endif
      enddo

C       write(*,*)'sntauchar(i)',sntauchar(1),sntauchar(2)
c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> H- + stau1/2
c      call NS_hcstausntau(gcsntaur)
      do i=1,2
         if(asntau1.gt.(gmstau(i)+cmass)) then
           sntau1hcstau(i)=g2s*mw**2*gcsntaur(1,i)**2*
     .      NS_lamb(gmstau(i)/asntau1,cmass/asntau1)/(16.D0*pi*asntau1)
         else
            sntau1hcstau(i)=0.D0
         endif
      end do

C      write(*,*)'sntau1hcstau(i)',sntau1hcstau(1),sntau1hcstau(2)
c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> W- + stau1/2
      do i=1,2
         if(asntau1.gt.(gmstau(i)+mw)) then
            sntau1wstau(i)=g2s/32.D0/pi/mw**2*asntau1**3*
     .                gwntau(1,i)**2*
     .                NS_lamb(gmstau(i)/asntau1,mw/asntau1)**3
         else
            sntau1wstau(i)=0.D0
         endif
      enddo

C       write(*,*)'sntau1wstau(i)',sntau1wstau(1),sntau1wstau(2)
c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

      sntautot2=sntauneut(1)+sntauneut(2)+sntauneut(3)+
     .          sntauneut(4)+sntauneut(5)+sntauchar(1)+sntauchar(2)+
     .          sntau1hcstau(1)+sntau1hcstau(2)+
     .          sntau1wstau(1)+sntau1wstau(2)
c--------------------------------------------------------------------- c
c ----- sntau 3-body decays and 3-body total widths ------------------ c
c--------------------------------------------------------------------- c
c Only: sntau -> nu_tau + stau1* + tau* (xintegsntaustau1star)
c       sntau -> nu_tau + stau1 + tau   (xintegsntaustau1)
c       sntau -> tau + stau1 + nu_tau   (xintegsntaustau1nutau)
      if(flagmulti.eq.1.D0) then
        Call NS_xintegsntau(xintegsntaustau1star,xintegsntaustau1,
     .    xintegsntaustau1nutau)
        
        sntautot3=xintegsntaustau1star+xintegsntaustau1+
     .    xintegsntaustau1nutau

        if (sntautot3.lt.multilim*sntautot2)Then
           sntautot3=0.d0
        endif
        sntautot=sntautot2+sntautot3
      do i=1,5
         brsntauneut(i) = sntauneut(i)/sntautot
      enddo
      do i=1,2
         brsntauchar(i)    = sntauchar(i)/sntautot
         brsntau1wstau(i)  = sntau1wstau(i)/sntautot
         brsntau1hcstau(i) = sntau1hcstau(i)/sntautot
      enddo
      if (sntautot3.ne.0.d0)Then
      brsntaustau1star=xintegsntaustau1star/sntautot
      brsntaustau1=xintegsntaustau1/sntautot
      brsntaustau1nutau=xintegsntaustau1nutau/sntautot
      endif
      endif

c      write(*,*) "Sneutr_tau_3-body BRs to stau1star,stau1,stau1+neut:"
c      write(*,*) brsntaustau1star,brsntaustau1,brsntaustau1nutau
C      WRITE(*,*) "brsntauneut(1,2)",brsntauneut(1),brsntauneut(2)
C      WRITE(*,*) "brsntauchar(1,2)",brsntauchar(1),brsntauchar(2)
C      WRITE(*,*) "brsntau1wstau(12)",brsntau1wstau(1),brsntau1wstau(2)
C      WRITE(*,*)"brsntau1hcstau(i)",brsntau1hcstau(1),brsntau1hcstau(2)

      END
c ==================================================================== c
c                        selectron_R 3-body decays                     c
c ==================================================================== c

        SUBROUTINE NS_xintegselr(xintegslstaustarltau,xintegslstaultau)
*
        INTEGER nx1t,ny1t
        DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
        DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
        DOUBLE PRECISION xintegslstaustarltau,xintegslstaultau
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/NS_pi/PI,SQR2
        COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
        COMMON/NS_nx1/nx1t,ny1t
*
        EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
        EXTERNAL NS_slstaustarltau,NS_slstaultau

c -------------------------------------------------------------------- c
c                        selectron_R 3-body decay into stau_1*         c
c                                neutralino exchange                   c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=(MTAU/MLR)**2
      xmu3=(MSL1/MLR)**2

      if(MLR.gt.(MSL1+MTAU)) then
         call NS_integ2(NS_slstaustarltau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegslstaustarltau=g2s**2/32.D0/(2.D0*pi)**3*MLR*sum
      else 
         xintegslstaustarltau=0.D0
      endif

c -------------------------------------------------------------------- c
c                        selectron_R 3-body decay into stau_1          c
c                                neutralino exchange                   c
c -------------------------------------------------------------------- c

      if(MLR.gt.(MSL1+MTAU)) then
         call NS_integ2(NS_slstaultau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegslstaultau=g2s**2/32.D0/(2.D0*pi)**3*MLR*sum
      else 
         xintegslstaultau=0.D0
      endif
        end

c ==================================================================== c
c        selectron_R 3-body decay into stau_1*, neutralino exchange    c
c ==================================================================== c
       DOUBLE PRECISION FUNCTION NS_slstaustarltau(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 
*
        xtau=(MTAU/MLR)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MLR)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau=(MSL1/MLR)**2

        NS_slstaustarltau=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MLR).and.
     .     (dabs(mneu(l)).gt.MLR)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
          NS_slstaustarltau=NS_slstaustarltau+1d0/(dneut(k)*dneut(l))*
     .        be(2,k)*be(2,l)*(btau(1,k)*btau(1,l)*
     .           dabs(mneu(k)*mneu(l))/MLR**2*(x1+x2-1d0+xstau-xtau)
     .         +atau(1,k)*atau(1,l)*((1d0-x1)*(1d0-x2)-xstau+xtau)
     .         +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .         +btau(1,k)*atau(1,l)*dsqrt(xneut(k)*xtau)*x1)
           endif
        enddo
        enddo

        end
c ==================================================================== c
c        selectron_R 3-body decay into stau_1, neutralino exchange     c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_slstaultau(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 

        xtau=(MTAU/MLR)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MLR)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau=(MSL1/MLR)**2

        NS_slstaultau=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MLR).and.
     .     (dabs(mneu(l)).gt.MLR)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_slstaultau=NS_slstaultau+1d0/(dneut(k)*dneut(l))*
     .        be(2,k)*be(2,l)*(atau(1,k)*atau(1,l)*
     .           dabs(mneu(k)*mneu(l))/MLR**2*(x1+x2-1d0+xstau-xtau)
     .         +btau(1,k)*btau(1,l)*((1d0-x1)*(1d0-x2)-xstau+xtau)
     .         +btau(1,k)*atau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .         +atau(1,k)*btau(1,l)*dsqrt(xneut(k)*xtau)*x1)
           endif
        enddo
        enddo

        end

c ==================================================================== c
c                        stau_2 3-body decays                          c
c ==================================================================== c

        SUBROUTINE NS_xintegstau2(xintegstau2stau1star,xintegstau2stau1,
     .    xintegstau2stau1nn)
* 
        INTEGER nx1t,ny1t
        DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
        DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
        DOUBLE PRECISION xintegstau2stau1star,xintegstau2stau1,
     .    xintegstau2stau1nn
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/NS_pi/PI,SQR2
        COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
        COMMON/NS_nx1/nx1t,ny1t
*
        EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
        EXTERNAL NS_stau2stau1star,NS_stau2stau1,NS_stau2stau1nn

c -------------------------------------------------------------------- c
c --------------------- stau_2 -> stau_1* ---------------------------- c
c -------------------------------------------------------------------- c
        xmu1=MTAU**2/MSL2**2
        xmu2=MTAU**2/MSL2**2
        xmu3=MSL1**2/MSL2**2

      if(MSL2.gt.(MSL1+2.D0*MTAU)) then
         call NS_integ2(NS_stau2stau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegstau2stau1star=g2s**2/32.D0/(2.D0*pi)**3*MSL2*sum
      else 
         xintegstau2stau1star=0.D0
      endif

c -------------------------------------------------------------------- c
c --------------------- stau_2 -> stau_1 ----------------------------- c
c -------------------------------------------------------------------- c

      if(MSL2.gt.(MSL1+2.D0*MTAU)) then
         call NS_integ2(NS_stau2stau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegstau2stau1=g2s**2/32.D0/(2.D0*pi)**3*MSL2*sum
      else 
         xintegstau2stau1=0.D0
      endif

c -------------------------------------------------------------------- c
c --------------------- stau_2 -> stau_1, chargino exchange ---------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=0d0
      xmu3=MSL1**2/MSL2**2

      if(MSL2.gt.MSL1) then
         call NS_integ2(NS_stau2stau1nn,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegstau2stau1nn=g2s**2/32.D0/(2.D0*pi)**3*MSL2*sum
      else 
         xintegstau2stau1nn=0.D0
      endif

       end

c ==================================================================== c
c        stau_2 3-body decay into stau_1*, neutralino exchange         c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_stau2stau1star(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 
*
        xtau=(MTAU/MSL2)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MSL2)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau1=(MSL1/MSL2)**2

        NS_stau2stau1star=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MSL2-MTAU).and.
     .     (dabs(mneu(l)).gt.MSL2-MTAU)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_stau2stau1star=NS_stau2stau1star+1d0/dneut(k)/dneut(l)*(
     .          (btau(1,k)*btau(1,l)*atau(2,k)*atau(2,l)+
     .           atau(1,k)*atau(1,l)*btau(2,k)*btau(2,l))*
     .          ((1.D0-x1)*(1.D0-x2)-xstau1+xtau*(x1-x2+xstau1
     .           -2.D0*xtau)+xtau)
     .         +(btau(1,k)*btau(1,l)*btau(2,k)*btau(2,l)+
     .           atau(1,k)*atau(1,l)*atau(2,k)*atau(2,l))*
     .           MNEU(k)*MNEU(l)/MSL2**2*(x1+x2-1.D0+xstau1
     .           -2.D0*xtau)
     .         +(btau(1,k)*atau(1,l)*atau(2,k)*btau(2,l)+
     .           atau(1,k)*btau(1,l)*btau(2,k)*atau(2,l))*
     .           2.D0*xtau*(-1.D0+x1-xtau)
     .         +dsqrt(xtau)*(x1-2.D0*xtau)*(MNEU(k)/MSL2*
     .          (btau(1,k)*atau(1,l)*btau(2,k)*btau(2,l)+
     .           atau(1,k)*btau(1,l)*atau(2,k)*atau(2,l))
     .         +MNEU(l)/MSL2*
     .          (btau(1,k)*atau(1,l)*atau(2,k)*atau(2,l)+
     .           atau(1,k)*btau(1,l)*btau(2,k)*btau(2,l)) )
     .         +dsqrt(xtau)*(-1.D0+x1+xstau1-2.D0*xtau)* 
     .         (MNEU(k)/MSL2*
     .          (btau(1,k)*btau(1,l)*btau(2,k)*atau(2,l)+
     .           atau(1,k)*atau(1,l)*atau(2,k)*btau(2,l))
     .         +MNEU(l)/MSL2*
     .          (btau(1,k)*btau(1,l)*atau(2,k)*btau(2,l)+
     .           atau(1,k)*atau(1,l)*btau(2,k)*atau(2,l)) )
     .         +(btau(1,l)*atau(1,k)*atau(2,k)*btau(2,l)+
     .           atau(1,l)*btau(1,k)*btau(2,k)*atau(2,l))*
     .           (-2.D0)*xtau*MNEU(k)*MNEU(l)/MSL2**2 )
           endif
        enddo
        enddo
  
        end

c ==================================================================== c
c        stau_2 3-body decay into stau_1, neutralino exchange          c
c ==================================================================== c
        DOUBLE PRECISION FUNCTION NS_stau2stau1(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 

        xtau=(MTAU/MSL2)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MSL2)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau1=(MSL1/MSL2)**2

        NS_stau2stau1=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MSL2-MTAU).and.
     .     (dabs(mneu(l)).gt.MSL2-MTAU)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_stau2stau1=NS_stau2stau1+1d0/(dneut(k)*dneut(l))*(
     .          (atau(1,k)*atau(1,l)*atau(2,k)*atau(2,l)+
     .           btau(1,k)*btau(1,l)*btau(2,k)*btau(2,l))*
     .          ((1.D0-x1)*(1.D0-x2)-xstau1+xtau*(x1-x2+xstau1
     .           -2.D0*xtau)+xtau)
     .         +(atau(1,k)*atau(1,l)*btau(2,k)*btau(2,l)+
     .           btau(1,k)*btau(1,l)*atau(2,k)*atau(2,l))*
     .           MNEU(k)*MNEU(l)/MSL2**2*(x1+x2-1.D0+xstau1
     .           -2.D0*xtau)
     .         +(atau(1,k)*btau(1,l)*atau(2,k)*btau(2,l)+
     .           btau(1,k)*atau(1,l)*btau(2,k)*atau(2,l))*
     .           2.D0*xtau*(-1.D0+x1-xtau)
     .         +dsqrt(xtau)*(x1-2.D0*xtau)*(MNEU(k)/MSL2*
     .          (atau(1,k)*btau(1,l)*btau(2,k)*btau(2,l)+
     .           btau(1,k)*atau(1,l)*atau(2,k)*atau(2,l))
     .         +MNEU(l)/MSL2*
     .          (atau(1,k)*btau(1,l)*atau(2,k)*atau(2,l)+
     .           btau(1,k)*atau(1,l)*btau(2,k)*btau(2,l)) )
     .         +dsqrt(xtau)*(-1.D0+x1+xstau1-2.D0*xtau)* 
     .         (MNEU(k)/MSL2*
     .          (atau(1,k)*atau(1,l)*btau(2,k)*atau(2,l)+
     .           btau(1,k)*btau(1,l)*atau(2,k)*btau(2,l))
     .         +MNEU(l)/MSL2*
     .          (atau(1,k)*atau(1,l)*atau(2,k)*btau(2,l)+
     .           btau(1,k)*btau(1,l)*btau(2,k)*atau(2,l)) )
     .         +(atau(1,l)*btau(1,k)*atau(2,k)*btau(2,l)+
     .           btau(1,l)*atau(1,k)*btau(2,k)*atau(2,l))*
     .           (-2.D0)*xtau*MNEU(k)*MNEU(l)/MSL2**2 )
           endif
        enddo
        enddo

        end
c ==================================================================== c
c        stau_2 3-body decay into stau_1, chargino exchange            c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_stau2stau1nn(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
        DOUBLE PRECISION ale(2,2),altau(2,2),alsne(2,2),blsne(2,2),
     .     alsnt(2,2),blsnt(2,2)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
*
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup5/ale,altau,alsne,blsne,alsnt,blsnt
         COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        do i=1,2
          xchar(i)=(MCH(i)/MSL2)**2
          dchar(i)=1d0-x1-xchar(i)
        enddo

        xstau1=(MSL1/MSL2)**2

        NS_stau2stau1nn=0.d0

        do k=1,2
          do l=1,2
           if((dabs(mch(k)).gt.MSL2).and.(dabs(mch(l)).gt.MSL2)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_stau2stau1nn=NS_stau2stau1nn+1.d0/(dchar(k)*dchar(l))*
     .        altau(2,k)*altau(2,l)*altau(1,k)*altau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1)
           endif
        enddo
        enddo

        end

c ==================================================================== c
c                          sntau 3-body decays                         c
c ==================================================================== c

        SUBROUTINE NS_xintegsntau(xintegsntaustau1star,xintegsntaustau1,
     .    xintegsntaustau1nutau)
* 
        INTEGER nx1t,ny1t
        DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
        DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
        DOUBLE PRECISION xintegsntaustau1star,xintegsntaustau1,
     .    xintegsntaustau1nutau
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU

* 
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/NS_pi/PI,SQR2
        COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
        COMMON/NS_nx1/nx1t,ny1t

      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_sntaustau1star,NS_sntaustau1,NS_sntaustau1nutau

c -------------------------------------------------------------------- c
c --------------------- sntau -> stau_1* ----------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=MTAU**2/MSNT**2
      xmu3=MSL1**2/MSNT**2

      if(MSNT.gt.(MSL1+MTAU)) then
        call NS_integ2(NS_sntaustau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsntaustau1star=g2s**2/32.D0/(2.D0*pi)**3*MSNT*sum
      else 
         xintegsntaustau1star=0.D0
      endif

c -------------------------------------------------------------------- c
c --------------------- sntau -> stau_1 ------------------------------ c
c -------------------------------------------------------------------- c

      if(MSNT.gt.(MSL1+MTAU)) then
         call NS_integ2(NS_sntaustau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsntaustau1=g2s**2/32.D0/(2.D0*pi)**3*MSNT*sum
      else 
         xintegsntaustau1=0.D0
      endif

c -------------------------------------------------------------------- c
c --------------------- sntau -> stau_1, chargino exchange ----------- c
c -------------------------------------------------------------------- c

      xmu1=MTAU**2/MSNT**2
      xmu2=0d0
      xmu3=MSL1**2/MSNT**2

      if(MSNT.gt.MTAU+MSL1) then
         call NS_integ2(NS_sntaustau1nutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsntaustau1nutau=g2s**2/32.D0/(2.D0*pi)**3*MSNT*sum
      else 
         xintegsntaustau1nutau=0.D0
      endif
c -------------------------------------------------------------------- c

       end

c ==================================================================== c
c        sntau 3-body decay into stau_1*, neutralino exchange         c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_sntaustau1star(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 

        xtau=(MTAU/MSNT)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MSNT)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau1=(MSL1/MSNT)**2

        NS_sntaustau1star=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MSNT).and.
     .     (dabs(mneu(l)).gt.MSNT)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sntaustau1star=NS_sntaustau1star+1d0/dneut(k)/dneut(l)*
     .          antau(1,k)*antau(1,l)*(btau(1,k)*btau(1,l)
     .         *((1.D0-x1)*(1.D0-x2)-xstau1+xtau)
     .         +(atau(1,k)*atau(1,l))*
     .           (MNEU(k)*MNEU(l))/MSNT**2*(x1+x2-1.D0+xstau1-xtau)
     .         +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .         +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)    
           endif
        enddo
        enddo

        end

c ==================================================================== c
c        sntau 3-body decay into stau_1, neutralino exchange          c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_sntaustau1(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 

        xtau=(MTAU/MSNT)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MSNT)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau1=(MSL1/MSNT)**2

        NS_sntaustau1=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MSNT).and.
     .     (dabs(mneu(l)).gt.MSNT)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sntaustau1=NS_sntaustau1+1d0/(dneut(k)*dneut(l))*
     .          antau(1,k)*antau(1,l)*(atau(1,k)*atau(1,l)
     .         *((1.D0-x1)*(1.D0-x2)-xstau1+xtau)
     .         +(btau(1,k)*btau(1,l))*
     .           (MNEU(k)*MNEU(l))/MSNT**2*(x1+x2-1.D0+xstau1-xtau)
     .         +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .         +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
             endif
        enddo
        enddo

        end

c ==================================================================== c
c        sntau 3-body decay into stau_1, chargino exchange             c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_sntaustau1nutau(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dchar(2),xchar(2)
        DOUBLE PRECISION ale(2,2),altau(2,2),alsne(2,2),blsne(2,2),
     .     alsnt(2,2),blsnt(2,2)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup5/ale,altau,alsne,blsne,alsnt,blsnt

        xtau=(MTAU/MSNT)**2
        do i=1,2
          xchar(i)=(MCH(i)/MSNT)**2
          dchar(i)=1d0-x1+xtau-xchar(i)
        enddo
        xstau1=(MSL1/MSNT)**2

        NS_sntaustau1nutau=0d0

        do k=1,2
          do l=1,2
           if((dabs(mch(k)).gt.MSNT).and.(dabs(mch(l)).gt.MSNT)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
          NS_sntaustau1nutau=NS_sntaustau1nutau+1d0/(dchar(k)*dchar(l))*
     .alsnt(1,k)*alsnt(1,l)*altau(1,k)*altau(1,l)*
     .((1.D0-x1)*(1.D0-x2)-xstau1+xtau*(xstau1+x1-x2-xtau))

       endif

        enddo
        enddo

        end

c ==================================================================== c
c                   selectron_l 3-body decays                          c
c ==================================================================== c

        SUBROUTINE NS_xintegsell(xintegsellstau1star,xintegsellstau1,
     .    xintegsellstau1nutau)
* 
        INTEGER nx1t,ny1t
        DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
        DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
        DOUBLE PRECISION xintegsellstau1star,xintegsellstau1,
     .    xintegsellstau1nutau
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/NS_pi/PI,SQR2
        COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
        COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_sellstau1star,NS_sellstau1,NS_sellstau1nutau
c -------------------------------------------------------------------- c
c ---------------- selectron_l -> stau_1* ---------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=MTAU**2/MLL**2
      xmu3=MSL1**2/MLL**2

      if(MLL.gt.(MSL1+MTAU)) then
         call NS_integ2(NS_sellstau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsellstau1star=g2s**2/32.D0/(2.D0*pi)**3*MLL*sum
      else 
         xintegsellstau1star=0.D0
      endif

c -------------------------------------------------------------------- c
c ---------------- selectron_l -> stau_1 ----------------------------- c
c -------------------------------------------------------------------- c

      if(MLL.gt.(MSL1+MTAU)) then
         call NS_integ2(NS_sellstau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsellstau1=g2s**2/32.D0/(2.D0*pi)**3*MLL*sum
      else 
         xintegsellstau1=0.D0
      endif

c -------------------------------------------------------------------- c
c --------------- selectron_l -> stau_1, chargino exchange ----------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=0d0
      xmu3=MSL1**2/MLL**2

      if(MLL.gt.MSL1) then
         call NS_integ2(NS_sellstau1nutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsellstau1nutau=g2s**2/32.D0/(2.D0*pi)**3*MLL*sum
      else 
         xintegsellstau1nutau=0.D0
      endif
c -------------------------------------------------------------------- c

       end

c ==================================================================== c
c   selectron_left 3-body decay into stau_1*, neutralino exchange      c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_sellstau1star(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 
*
        xtau=(MTAU/MLL)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MLL)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau1=(MSL1/MLL)**2

        NS_sellstau1star=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MLL).and.
     .     (dabs(mneu(l)).gt.MLL)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sellstau1star=NS_sellstau1star+1d0/dneut(k)/dneut(l)*
     .        ae(1,k)*ae(1,l)*(btau(1,k)*btau(1,l)*
     .          ((1.D0-x1)*(1.D0-x2)-xstau1+xtau)
     .         +(atau(1,k)*atau(1,l))*
     .           MNEU(k)*MNEU(l)/MLL**2*(x1+x2-1.D0+xstau1-xtau)
     .         +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .         +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
        enddo

        end

c ==================================================================== c
c     selectron_l 3-body decay into stau_1, neutralino exchange        c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_sellstau1(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 

        xtau=(MTAU/MLL)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MLL)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau1=(MSL1/MLL)**2

        NS_sellstau1=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MLL).and.
     .     (dabs(mneu(l)).gt.MLL)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sellstau1=NS_sellstau1+1d0/(dneut(k)*dneut(l))*
     .         ae(1,k)*ae(1,l)*(atau(1,k)*atau(1,l)*
     .          ((1.D0-x1)*(1.D0-x2)-xstau1+xtau)
     .         +(btau(1,k)*btau(1,l))*
     .           MNEU(k)*MNEU(l)/MLL**2*(x1+x2-1.D0+xstau1-xtau)
     .         +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .         +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
        enddo

        end


c ==================================================================== c
c     selectron_l 3-body decay into stau_1, chargino exchange          c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_sellstau1nutau(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
        DOUBLE PRECISION ale(2,2),altau(2,2),alsne(2,2),blsne(2,2),
     .     alsnt(2,2),blsnt(2,2)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup5/ale,altau,alsne,blsne,alsnt,blsnt

        do i=1,2
          xchar(i)=(MCH(i)/MLL)**2
          dchar(i)=1d0-x1-xchar(i)
        enddo

        xstau1=(MSL1/MLL)**2

        NS_sellstau1nutau=0d0

        do k=1,2
          do l=1,2
           if((dabs(mch(k)).gt.MLL).and.(dabs(mch(l)).gt.MLL)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sellstau1nutau=NS_sellstau1nutau+1d0/(dchar(k)*dchar(l))*
     .        ale(1,k)*ale(1,l)*altau(1,k)*altau(1,l)*
     .        ((1.D0-x1)*(1.D0-x2)-xstau1)

           endif
        enddo
        enddo

        end

c ==================================================================== c
c                          sne 3-body decays                         c
c ==================================================================== c

        SUBROUTINE NS_xintegsne(xintegsnestau1star,xintegsnestau1,
     .    xintegsnestau1nutau)
* 
        INTEGER nx1t,ny1t
        DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
        DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
        DOUBLE PRECISION xintegsnestau1star,xintegsnestau1,
     .    xintegsnestau1nutau
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/NS_pi/PI,SQR2
        COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
        COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_snestau1star,NS_snestau1,NS_snestau1char
c -------------------------------------------------------------------- c
c --------------------- sne -> stau_1* ----------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=MTAU**2/MNL**2
      xmu3=MSL1**2/MNL**2

      if(MNL.gt.(MSL1+MTAU)) then
         call NS_integ2(NS_snestau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnestau1star=g2s**2/32.D0/(2.D0*pi)**3*MNL*sum
      else 
         xintegsnestau1star=0.D0
      endif

c -------------------------------------------------------------------- c
c --------------------- sne -> stau_1 -------------------------------- c
c -------------------------------------------------------------------- c

      if(MNL.gt.(MSL1+MTAU)) then
         call NS_integ2(NS_snestau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnestau1=g2s**2/32.D0/(2.D0*pi)**3*MNL*sum
      else 
         xintegsnestau1=0.D0
      endif

c -------------------------------------------------------------------- c
c --------------------- sne -> stau_1, chargino exchange ------------- c
c -------------------------------------------------------------------- c

      xmu1=MTAU**2/MNL**2
      xmu2=0d0
      xmu3=MSL1**2/MNL**2

      if(MNL.gt.MSL1) then
         call NS_integ2(NS_snestau1char,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnestau1nutau=g2s**2/32.D0/(2.D0*pi)**3*MNL*sum
      else 
         xintegsnestau1nutau=0.D0
      endif
c -------------------------------------------------------------------- c

       end

c ==================================================================== c
c         snu_e 3-body decay into stau_1*, neutralino exchange         c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_snestau1star(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 
*
        xtau=(MTAU/MNL)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MNL)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau1=(MSL1/MNL)**2

        NS_snestau1star=0d0
     
        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MNL).and.
     .     (dabs(mneu(l)).gt.MNL)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_snestau1star=NS_snestau1star+1d0/dneut(k)/dneut(l)*
     .          anu(1,k)*anu(1,l)*(btau(1,k)*btau(1,l)
     .         *((1.D0-x1)*(1.D0-x2)-xstau1+xtau)
     .         +(atau(1,k)*atau(1,l))*
     .          (MNEU(k)*MNEU(l))/MNL**2*(x1+x2-1.D0+xstau1-xtau)
     .         +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .         +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)   
           endif
        enddo
        enddo
       
        end

c ==================================================================== c
c        snu_e 3-body decay into stau_1, neutralino exchange           c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_snestau1(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
        DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau 

        xtau=(MTAU/MNL)**2
        do i=1,5
          xneut(i)=(MNEU(i)/MNL)**2
          dneut(i)=1d0-x1+xtau-xneut(i)
        enddo
        xstau1=(MSL1/MNL)**2

        NS_snestau1=0d0

        do k=1,5
          do l=1,5
           if((dabs(mneu(k)).gt.MNL).and.
     .     (dabs(mneu(l)).gt.MNL)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_snestau1=NS_snestau1+1d0/(dneut(k)*dneut(l))*
     .          anu(1,k)*anu(1,l)*(atau(1,k)*atau(1,l)
     .         *((1.D0-x1)*(1.D0-x2)-xstau1+xtau)
     .         +btau(1,k)*btau(1,l)*
     .          (MNEU(k)*MNEU(l))/MNL**2*(x1+x2-1.D0+xstau1-xtau)
     .         +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .         +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)  
          endif
        enddo
        enddo

        end

c ==================================================================== c
c        snu_e 3-body decay into stau_1, chargino exchange             c
c ==================================================================== c

        DOUBLE PRECISION FUNCTION NS_snestau1char(x1,x2)
*
        INTEGER i,k,l
        DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
        DOUBLE PRECISION ale(2,2),altau(2,2),alsne(2,2),blsne(2,2),
     .     alsnt(2,2),blsnt(2,2)

        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
*
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
        COMMON/NS_coup5/ale,altau,alsne,blsne,alsnt,blsnt

        do i=1,2
          xchar(i)=(MCH(i)/MNL)**2
          dchar(i)=1d0-x1-xchar(i)
        enddo
        xstau1=(MSL1/MNL)**2

        NS_snestau1char=0d0

        do k=1,2
          do l=1,2
           if((dabs(mch(k)).gt.MNL).and.(dabs(mch(l)).gt.MNL)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
           NS_snestau1char=NS_snestau1char+1d0/(dchar(k)*dchar(l))*
     .       (altau(1,k)*altau(1,l))
     .       *alsne(1,k)*alsne(1,l)*
     .       ((1.D0-x1)*(1.D0-x2)-xstau1)

           
           endif
        enddo
        enddo

        end
