c ==================================================================== c
c                bottom squark decays (2, 3 body)                      c
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
      SUBROUTINE NS_SBOTTOM
*
      IMPLICIT NONE
      INTEGER I,J,jsign,nj,k
      DOUBLE PRECISION amuv,lamv,amuvdiv
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION alp, nf
      DOUBLE PRECISION amurefer
      DOUBLE PRECISION amsq,scalmur
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS      
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION rmtc,rmbc,rmtauc
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION abot(2,5),bbot(2,5),alsbot(2,2),aksbot(2,2)
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2)
      DOUBLE PRECISION gctbr(2,2)
      DOUBLE PRECISION Hstaustaur(3,2,2),Astaustaur(2,2,2)
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION gmst(2)
      DOUBLE PRECISION delta11c,delta12c,delta13c,delta14c,delta15c,
     .          delta1H(3),delta2H(3),delta3H(3),delta4H(3),delta5H(3),
     .          delta1A(2),delta2A(2),delta3A(2),delta4A(2),delta5A(2),
     .          delta21c,delta22c,delta23c,delta24c,delta25c,
     .          del1,del2,del3,del4,del5,
     .          adel1,adel2,adel3,adel4,adel5,
     .          bdel1,bdel2,bdel3,bdel4,bdel5
      DOUBLE PRECISION NS_lamb,resum
      DOUBLE PRECISION NS_glbneut,NS_grbneut,NS_corrreali,NS_glbchar,
     .         NS_grbchar
      DOUBLE PRECISION NS_gamtop1,NS_gamtop2,NS_gamglui1,NS_gamglui2,
     .         NS_gamglui3,NS_gam11,NS_gam12,NS_gamvirt,NS_gamcfdec,
     .         NS_gamreal
      DOUBLE PRECISION NS_gvirtgl,NS_gvirtmix,NS_sbotstop1719,
     .         NS_dcounterhc,NS_realcorr
      DOUBLE PRECISION NS_gvirtmixdiv,NS_botneut1719,NS_dcounterneut
      DOUBLE PRECISION NS_gluonvertex,NS_gluinoWvertex,NS_gluinoZvertex
      DOUBLE PRECISION NS_wavefuncvertex,NS_quarkmixW,NS_quarkmixZ
      DOUBLE PRECISION NS_realgluonem
*
      DOUBLE PRECISION sb1neutt(5),sb2neutt(5),qcdsb1neut(5),
     .     qcdsb2neut(5),sb1chart(2),sb2chart(2),qcdsb1chart(2),
     .     qcdsb2chart(2),sb1hcst(2),sb2hcst(2),qcdsb1hcst(2),
     .     qcdsb2hcst(2),sb1glui,sb2glui,qcdsb1glui,qcdsb2glui,
     .     sb1wst(2),sb2wst(2),qcdsb1wst(2),qcdsb2wst(2),
     .     sb2H(3),qcdsb2H(3),sb2A(2),qcdsb2A(2),
     .     sb2zbot,qcdsb2zbot
      DOUBLE PRECISION xintegsbstau(2,2),xintegsbsntau(2,2),
     .     xintegsbsel(2,2),
     .     xintegsbtstsb(2,2),xintegsbtbstb(2,2),
     .     xintegsbtaustnu(2,2),xintegsbelstnu(2,2),
     .     xintegsbupstdow(2,2),xintegsbsnel(2),
     .     xintegsb2sb1bb,xintegsb2sb1starbb,xintegsb2sb1tt,
     .     xintegsb2sb1uu,xintegsb2sb1dd,xintegsb2sb1ee,
     .     xintegsb2sb1nunu,xintegsb2sb1tautau
      DOUBLE PRECISION sbottot(2),sbottot2(2),sbottotmulti(2)
      DOUBLE PRECISION brsb1neutt(5),brsb2neutt(5),brsb1chart(2),
     .          brsb2chart(2),brsb1hcst(2),brsb2hcst(2),
     .          brsb1glui,brsb2glui,brsb1wst(2),
     .          brsb2wst(2),brsb2H(3),brsb2A(2),brsb2zbot
      DOUBLE PRECISION  brsbstau(2,2),brsbsntau(2,2),brsbsel(2,2),
     .          brsbtstsb(2,2),brsbtbstb(2,2),brsbtaustnu(2,2),
     .          brsbelstnu(2,2),brsbupstdow(2,2),brsbsnel(2),
     .          brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .          brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .          brsb2sb1tautau
      DOUBLE PRECISION sbottot2lo(2),sbottot2nlo(2)
      DOUBLE PRECISION flagmulti,flagqcd,flagloop,multilim
*
      COMMON/SBOTTOM_WIDTH/sbottot,sbottot2,sbottotmulti
      COMMON/SBOTTOM_BR_2BD/brsb1neutt,brsb2neutt,brsb1chart,
     .          brsb2chart,brsb1hcst,brsb2hcst,
     .          brsb1glui,brsb2glui,brsb1wst,
     .          brsb2wst,brsb2H,brsb2A,brsb2zbot
      COMMON/SBOTTOM_BR_3BD/brsbstau,brsbsntau,brsbsel,
     .          brsbtstsb,brsbtbstb,brsbtaustnu,
     .          brsbelstnu,brsbupstdow,brsbsnel,
     .          brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .          brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .          brsb2sb1tautau
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_runmcalc/rmtc,rmbc,rmtauc
      COMMON/NS_refscale/amurefer
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/NS_coup20/gwtb,gwntau
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr  
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr 
      COMMON/NS_HIGGSTAUTAU/Hstaustaur,Astaustaur
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_multilim/multilim
*
      EXTERNAL NS_lamb
      EXTERNAL NS_glbneut,NS_grbneut,NS_corrreali,NS_glbchar,NS_grbchar
      EXTERNAL NS_gamtop1,NS_gamtop2,NS_gamglui1,NS_gamglui2,
     .         NS_gamglui3,NS_gam11,NS_gam12,NS_gamvirt,NS_gamcfdec,
     .         NS_gamreal
      EXTERNAL NS_gvirtgl,NS_gvirtmix,NS_sbotstop1719,NS_dcounterhc,
     .         NS_realcorr
      EXTERNAL NS_gvirtmixdiv,NS_botneut1719,NS_dcounterneut
      EXTERNAL NS_gluonvertex,NS_gluinoWvertex,NS_gluinoZvertex
      EXTERNAL NS_wavefuncvertex,NS_quarkmixW,NS_quarkmixZ
      EXTERNAL NS_realgluonem
      EXTERNAL resum
c -------------------------------------------------------------------- c
c                    sbottom1/2 2- and 3-body decays                   c
c -------------------------------------------------------------------- c
c -- initialization --

      do i=1,2
         sbottot(i)    = 0.D0
         sbottot2(i)   = 0.D0
         sbottot2lo(i) = 0.D0
         sbottot2nlo(i)= 0.D0
         sbottotmulti(i) = 0.D0

         sb1chart(i) = 0.D0
         sb1hcst(i)  = 0.D0
         sb1wst(i)   = 0.D0

         sb2chart(i) = 0.D0
         sb2hcst(i)  = 0.D0
         sb2wst(i)   = 0.D0

         qcdsb1chart(i) = 0.D0
         qcdsb1hcst(i)  = 0.D0
         qcdsb1wst(i)   = 0.D0

         qcdsb2chart(i) = 0.D0
         qcdsb2hcst(i)  = 0.D0
         qcdsb2wst(i)   = 0.D0

         xintegsbsnel(i) = 0.D0
         do j=1,2
            xintegsbstau(i,j)    = 0.D0
            xintegsbsntau(i,j)   = 0.D0
            xintegsbsel(i,j)     = 0.D0
            xintegsbtstsb(i,j)   = 0.D0
            xintegsbtbstb(i,j)   = 0.D0
            xintegsbtaustnu(i,j) = 0.D0
            xintegsbelstnu(i,j)  = 0.D0
            xintegsbupstdow(i,j) = 0.D0
         enddo
      enddo

      do j=1,5
         sb1neutt(j) = 0.D0
         sb2neutt(j) = 0.D0
         qcdsb1neut(j) = 0.D0
         qcdsb2neut(j) = 0.D0
      enddo

      sb1glui = 0.D0
      sb2glui = 0.D0
      sb2zbot = 0.D0
      qcdsb1glui = 0.D0
      qcdsb2glui = 0.D0
      qcdsb2zbot = 0.D0

      DO i=1,3
         sb2H(i)=0.D0
         qcdsb2H(i)=0.D0
      ENDDO
      DO i=1,2
         sb2A(i)=0.D0
         qcdsb2A(i)=0.D0
      ENDDO

      xintegsb2sb1bb     = 0.D0
      xintegsb2sb1tt     = 0.D0
      xintegsb2sb1uu     = 0.D0
      xintegsb2sb1dd     = 0.D0
      xintegsb2sb1ee     = 0.D0
      xintegsb2sb1nunu   = 0.D0
      xintegsb2sb1tautau = 0.D0
      xintegsb2sb1starbb = 0.D0

      DO i=1,5
         brsb1neutt(i) = 0.D0
         brsb2neutt(i) = 0.D0
      ENDDO
      DO i=1,2
         brsb1chart(i) = 0.D0
         brsb2chart(i) = 0.D0
         brsb1hcst(i) = 0.D0
         brsb2hcst(i) = 0.D0
      ENDDO
      brsb1glui = 0.D0
      brsb2glui = 0.D0
      brsb2zbot = 0.D0
      DO i=1,2
         brsb1wst(i) = 0.D0
         brsb2wst(i) = 0.D0
      ENDDO
      DO i=1,3
         brsb2H(i) = 0.D0
      ENDDO
      DO i=1,2
         brsb2A(i) = 0.D0
      ENDDO
      DO i=1,2
         DO j=1,2
            brsbstau(i,j) = 0.D0
            brsbsntau(i,j) = 0.D0
            brsbsel(i,j) = 0.D0
            brsbtstsb(i,j) = 0.D0
            brsbtbstb(i,j) = 0.D0
            brsbtaustnu(i,j) = 0.D0
            brsbelstnu(i,j) = 0.D0
            brsbupstdow(i,j) = 0.D0
         ENDDO
         brsbsnel(i) = 0.D0
      ENDDO
      brsb2sb1bb = 0.D0
      brsb2sb1starbb = 0.D0
      brsb2sb1tt = 0.D0
      brsb2sb1uu = 0.D0
      brsb2sb1dd = 0.D0
      brsb2sb1ee = 0.D0
      brsb2sb1nunu = 0.D0
      brsb2sb1tautau = 0.D0

      DO i=1,2
         gmst(i)=0.D0
      ENDDO
      amsq=0.D0
      DO i=1,3
         delta1H(i)=0.D0
         delta2H(i)=0.D0
         delta3H(i)=0.D0
         delta4H(i)=0.D0
         delta5H(i)=0.D0
      ENDDO
      DO i=1,2
         delta1A(i)=0.D0
         delta2A(i)=0.D0
         delta3A(i)=0.D0
         delta4A(i)=0.D0
         delta5A(i)=0.D0
      ENDDO
      delta11c=0.D0
      delta12c=0.D0
      delta13c=0.D0
      delta14c=0.D0
      delta15c=0.D0
      delta21c=0.D0
      delta22c=0.D0
      delta23c=0.D0
      delta24c=0.D0
      delta25c=0.D0
      del1=0.D0
      del2=0.D0
      del3=0.D0
      del4=0.D0
      del5=0.D0
      adel1=0.D0
      adel2=0.D0
      adel3=0.D0
      adel4=0.D0
      adel5=0.D0
      bdel1=0.D0
      bdel2=0.D0
      bdel3=0.D0
      bdel4=0.D0
      bdel5=0.D0   

c -- stop masses --
      gmst(1) = ast1
      gmst(2) = ast2
c -------------------------------------------------------------------- c
c For QCD corrections: the fixed scale is amurefer = Q, where
c the couplings are defined:
      amuvdiv = amuv/amurefer
c -------------------------------------------------------------------- c
c  sbottom1 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + bottom
      do i=1,5
         if(asb1.gt.(amneut(i)+amb)) then
            sb1neutt(i)=g2s*((abot(1,i)**2+bbot(1,i)**2)*(asb1**2-
     .           amb**2-amneut(i)**2)-4.D0*abot(1,i)*bbot(1,i)*
     .           amb*xmneut(i))*NS_lamb(amb/asb1,amneut(i)/asb1)
     .           /(16*pi*asb1)
         else
            sb1neutt(i)=0.D0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do j=1,5
         if(asb1.gt.(amneut(j)+amb)) then
            if(xmneut(j).le.0.D0) then
               jsign = 1
            else
               jsign = 0
            endif
            qcdsb1neut(j) = -g2s/24.D0/pi**2/asb1*gs2/(4.D0*pi)*
     .           ((bbot(1,j)*NS_glbneut(1,j,amuv,amuvdiv,lamv)
     .            +abot(1,j)*NS_grbneut(1,j,amuv,amuvdiv,lamv))*
     .           (asb1**2-amb**2-amneut(j)**2)
     .           -2.D0*(bbot(1,j)*NS_grbneut(1,j,amuv,amuvdiv,lamv)
     .                 +abot(1,j)*NS_glbneut(1,j,amuv,amuvdiv,lamv))*
     .           amb*xmneut(j))*NS_lamb(amb/asb1,amneut(j)/asb1) 
     .           +g2s/(6.D0*pi**2*asb1)*gs2/(4.D0*pi)*
     .           NS_corrreali(amb,amneut(j),asb1,lamv,1,jsign,1,j,2)
         else
            qcdsb1neut(j) = 0.D0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sbottom2 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + bottom
      do i=1,5
         if(asb2.gt.(amneut(i)+amb)) then
            sb2neutt(i)=g2s*((abot(2,i)**2+bbot(2,i)**2)*(asb2**2-
     .           amb**2-amneut(i)**2)-4.D0*abot(2,i)*bbot(2,i)*
     .           amb*xmneut(i))*NS_lamb(amb/asb2,amneut(i)/asb2)
     .           /(16*pi*asb2)
         else
            sb2neutt(i)=0.D0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do j=1,5
         if(asb2.gt.(amneut(j)+amb)) then
            if(xmneut(j).le.0.D0) then
               jsign = 1
            else
               jsign = 0
            endif

            qcdsb2neut(j) = -g2s/24.D0/pi**2/asb2*gs2/(4.D0*pi)*
     .           ((bbot(2,j)*NS_glbneut(2,j,amuv,amuvdiv,lamv)
     .            +abot(2,j)*NS_grbneut(2,j,amuv,amuvdiv,lamv))*
     .           (asb2**2-amb**2-amneut(j)**2)
     .           -2.D0*(bbot(2,j)*NS_grbneut(2,j,amuv,amuvdiv,lamv)
     .                 +abot(2,j)*NS_glbneut(2,j,amuv,amuvdiv,lamv))*
     .           amb*xmneut(j))*NS_lamb(amb/asb2,amneut(j)/asb2) 
     .           +g2s/(6.D0*pi**2*asb2)*gs2/(4.D0*pi)*
     .           NS_corrreali(amb,amneut(j),asb2,lamv,1,jsign,2,j,2)
         else
            qcdsb2neut(j) = 0.D0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sbottom1 --> chi-_1/chi-_2 + top
      do i=1,2
         if(asb1.gt.(amchar(i)+amt)) then
            sb1chart(i)=g2s*((alsbot(1,i)**2+aksbot(1,i)**2)*
     .           (asb1**2-amt**2-amchar(i)**2)
     .           -4.D0*alsbot(1,i)*aksbot(1,i)*
     .           amt*xmchar(i))*NS_lamb(amt/asb1,amchar(i)/asb1)
     .           /(16*pi*asb1)
         else
            sb1chart(i)=0.D0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do j=1,2,1
         if(asb1.gt.(amchar(j)+amt)) then
            jsign = 0

            qcdsb1chart(j) = -g2s/24.D0/pi**2/asb1*gs2/(4.D0*pi)*
     .           ((aksbot(1,j)*NS_glbchar(1,j,amuv,amuvdiv,lamv)
     .            +alsbot(1,j)*NS_grbchar(1,j,amuv,amuvdiv,lamv))*
     .           (asb1**2-amt**2-amchar(j)**2)
     .           -2.D0*(aksbot(1,j)*NS_grbchar(1,j,amuv,amuvdiv,lamv)
     .                 +alsbot(1,j)*NS_glbchar(1,j,amuv,amuvdiv,lamv))*
     .           amt*xmchar(j))*NS_lamb(amt/asb1,amchar(j)/asb1) 
     .           +g2s/(6.D0*pi**2*asb1)*gs2/(4.D0*pi)*
     .           NS_corrreali(amt,amchar(j),asb1,lamv,2,jsign,1,j,2)
         else
            qcdsb1chart(j) = 0.D0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sbottom2 --> chi-_1/chi-_2 + top
      do i=1,2
         if(asb2.gt.(amchar(i)+amt)) then
            sb2chart(i)=g2s*((alsbot(2,i)**2+aksbot(2,i)**2)*
     .           (asb2**2-amt**2-amchar(i)**2)
     .           -4.D0*alsbot(2,i)*aksbot(2,i)*
     .           amt*xmchar(i))*NS_lamb(amt/asb2,amchar(i)/asb2)
     .           /(16*pi*asb2)
         else
            sb2chart(i)=0.D0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do j=1,2
         if(asb2.gt.(amchar(j)+amt)) then
            jsign = 0

            qcdsb2chart(j) = -g2s/24.D0/pi**2/asb2*gs2/(4.D0*pi)*
     .           ((aksbot(2,j)*NS_glbchar(2,j,amuv,amuvdiv,lamv)
     .            +alsbot(2,j)*NS_grbchar(2,j,amuv,amuvdiv,lamv))*
     .           (asb2**2-amt**2-amchar(j)**2)
     .           -2.D0*(aksbot(2,j)*NS_grbchar(2,j,amuv,amuvdiv,lamv)
     .                 +alsbot(2,j)*NS_glbchar(2,j,amuv,amuvdiv,lamv))*
     .           amt*xmchar(j))*NS_lamb(amt/asb2,amchar(j)/asb2) 
     .           +g2s/(6.D0*pi**2*asb2)*gs2/(4.D0*pi)*
     .           NS_corrreali(amt,amchar(j),asb2,lamv,2,jsign,2,j,2)
         else
            qcdsb2chart(j) = 0.D0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sbottom1 --> gluino + bottom
      if(asb1.gt.(mgluino+amb)) then
         sb1glui = 8.D0*gs2*((asb1**2-amb**2-mgluino**2)+4.D0*amb*
     .        mgluino*dsin(theb)*dcos(theb))*
     .        NS_lamb(amb/asb1,mgluino/asb1)/(16.D0*pi*asb1)/3.D0  
      else
         sb1glui = 0.D0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asb1.gt.(mgluino+amb)) then
         amsq    = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         scalmur = amurefer
         alp     = gs2/(4.D0*pi)
         nf      = 6.D0
         qcdsb1glui = 8.D0*pi*alp/3.D0/amb**2*sb1glui*
     .        NS_gamtop1(asb1,asb2,amb,mgluino,theb,1,amuv,lamv) +
     .        16.D0*pi*alp**2*NS_lamb(amb/asb1,mgluino/asb1)/
     .        (9.D0*amb**2*asb1)*
     .        NS_gamtop2(asb1,asb2,amb,mgluino,theb,1,amuv) +
     .        4.D0*pi*alp/mgluino**2*sb1glui*(nf-2.D0)*
     .        NS_gamglui1(ast1,ast2,amsq,amt,mgluino,amuv) +
     .        2.D0*pi*alp/mgluino**2*sb1glui*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,
     .                    mgluino,1,amuv) +
     .        4.D0*pi*alp*3.D0/mgluino**2*sb1glui*
     .        NS_gamglui3(mgluino,amuv,lamv) +
     .        8.D0*4.D0/3.D0*pi*alp*sb1glui*
     .        NS_gam11(asb1,asb2,amb,mgluino,theb,1,amuv,lamv) +
     .        8.D0*16.D0/9.D0*pi*alp**2/asb1*
     .        NS_lamb(amb/asb1,mgluino/asb1)*
     .        NS_gam12(asb1,asb2,amb,mgluino,theb,1,amuv,lamv,scalmur) +
     .        alp**2*NS_lamb(amb/asb1,mgluino/asb1)/asb1*
     .        NS_gamvirt(asb1,asb2,amb,mgluino,theb,1,amuv,lamv) +
     .        alp**2*NS_gamreal(asb1,amb,mgluino,theb,1,lamv) +
     .        alp/(4.D0*pi)*sb1glui*
     .        NS_gamcfdec(ast1,ast2,amt,asb1,asb2,amb,mgluino,amsq,amuv,
     .        scalmur)
      else
         qcdsb1glui = 0.D0
      endif
      endif
      
c -------------------------------------------------------------------- c
c  sbottom2 --> gluino + bottom
      if(asb2.gt.(mgluino+amb)) then
         sb2glui = 8.D0*gs2*((asb2**2-amb**2-mgluino**2)-4.D0*amb*
     .        mgluino*dsin(theb)*dcos(theb))*
     .        NS_lamb(amb/asb2,mgluino/asb2)/(16.D0*pi*asb2)/3.D0
      else
         sb2glui = 0.D0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asb2.gt.(mgluino+amb)) then
         amsq    = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         scalmur = amurefer
         alp     = gs2/(4.D0*pi)
         nf      = 6.D0
         qcdsb2glui = 8.D0*pi*alp/3.D0/amb**2*sb2glui*
     .        NS_gamtop1(asb2,asb1,amb,mgluino,theb,2,amuv,lamv) +
     .        16.D0*pi*alp**2*NS_lamb(amb/asb2,mgluino/asb2)/
     .        (9.D0*amb**2*asb2)*
     .        NS_gamtop2(asb2,asb1,amb,mgluino,theb,2,amuv) +
     .        4.D0*pi*alp/mgluino**2*sb2glui*(nf-2.D0)*
     .        NS_gamglui1(ast2,ast1,amsq,amt,mgluino,amuv) +
     .        2.D0*pi*alp/mgluino**2*sb2glui*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,
     .                    mgluino,2,amuv) +
     .        4.D0*pi*alp*3.D0/mgluino**2*sb2glui*
     .        NS_gamglui3(mgluino,amuv,lamv) +
     .        8.D0*4.D0/3.D0*pi*alp*sb2glui*
     .        NS_gam11(asb2,asb1,amb,mgluino,theb,2,amuv,lamv) +
     .        8.D0*16.D0/9.D0*pi*alp**2/asb2*
     .        NS_lamb(amb/asb2,mgluino/asb2)*
     .        NS_gam12(asb2,asb1,amb,mgluino,theb,2,amuv,lamv,scalmur) +
     .        alp**2*NS_lamb(amb/asb2,mgluino/asb2)/asb2*
     .        NS_gamvirt(asb2,asb1,amb,mgluino,theb,2,amuv,lamv) +
     .        alp**2*NS_gamreal(asb2,amb,mgluino,theb,2,lamv) +
     .        alp/(4.D0*pi)*sb2glui*
     .        NS_gamcfdec(ast2,ast1,amt,asb2,asb1,amb,mgluino,amsq,amuv,
     .        scalmur)
      else
         qcdsb2glui = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c  sbottom1 --> H- + stop1/2
      do i=1,2
         if(asb1.gt.(gmst(i)+cmass)) then
            sb1hcst(i)=g2s*amw**2*gctbr(i,1)**2*
     .           NS_lamb(gmst(i)/asb1,cmass/asb1)/(16.D0*pi*asb1)
         else
            sb1hcst(i)=0.D0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do nj=1,2
      if(asb1.gt.(gmst(nj)+cmass)) then
         alp = gs2/(4.D0*pi)
         delta11c = -dsqrt(2.D0)*amw**2*gctbr(nj,1)*
     .        NS_gvirtgl(gmst(nj),cmass,asb1,lamv,amuv)
         delta12c =  -dsqrt(2.D0)*amw**2*gctbr(3-nj,1)*
     .        1.D0/(gmst(3-nj)**2-gmst(nj)**2)*
     .        NS_gvirtmix(ast1,ast2,gmst(nj),mgluino,rmtc,thet,amuv)
     .        -dsqrt(2.D0)*amw**2*gctbr(nj,2)*
     .        1.D0/(asb2**2-asb1**2)*
     .        NS_gvirtmix(asb1,asb2,asb1,mgluino,rmbc,theb,amuv)
         delta13c = NS_sbotstop1719(amuv,1,nj)
         delta14c = NS_dcounterhc(asb1,rmbc,theb,1,gmst(nj),rmtc,
     .                         thet,nj,mgluino,amuv,amuvdiv,lamv,nj,1)
         delta15c = NS_realcorr(cmass,asb1,gmst(nj),lamv,6,0,nj,1,asb1)

         qcdsb1hcst(nj) = -g2s*amw**2/(24.D0*dsqrt(2.D0)*pi*amw**2*
     .        asb1)*alp/pi*NS_lamb(gmst(nj)/asb1,cmass/asb1)*
     .        gctbr(nj,1)*(delta11c+delta12c+delta13c+delta14c+delta15c)

      else
         qcdsb1hcst(nj) = 0.D0
      endif
      enddo
      endif
c -------------------------------------------------------------------- c
c sbottom2 --> H(K) + sbottom1
      DO K=1,3
         if(asb2.gt.(asb1+SMASS(K))) then
         sb2H(K)=g2s*amz**4/amw**2*Hsbotsbotr(K,2,1)**2*
     .         NS_lamb(asb1/asb2,SMASS(K)/asb2)/(16.D0*pi*asb2)
      else
         sb2H(K)=0.D0
      endif
      ENDDO
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      DO K=1,3
         if(asb2.gt.(asb1+SMASS(K))) then
         alp = gs2/(4.D0*pi)
         delta1H(K) = -dsqrt(2.D0)*amz**2*Hsbotsbotr(K,2,1)*
     .        NS_gvirtgl(asb2,SMASS(K),asb1,lamv,amuv)
         delta2H(K) =  -dsqrt(2.D0)*amz**2*Hsbotsbotr(K,2,2)*
     .        1.D0/(asb2**2-asb1**2)*
     .        NS_gvirtmixdiv(asb1,asb2,asb1,mgluino,rmbc,theb,amuv)
     .        -dsqrt(2.D0)*amz**2*Hsbotsbotr(K,1,1)*
     .        1.D0/(asb1**2-asb2**2)*
     .        NS_gvirtmixdiv(asb1,asb2,asb2,mgluino,rmbc,theb,amuv)
         delta3H(K) = NS_botneut1719(K,amuv)
         delta4H(K) = NS_dcounterneut(asb1,asb2,rmbc,theb,mgluino,amuv,
     .                          amuvdiv,lamv,2,K)
         delta5H(K) = NS_realcorr(SMASS(K),asb2,asb1,lamv,K,2,2,1,asb2)


         qcdsb2H(K) = -g2s*amz**2/(24.D0*dsqrt(2.D0)*pi*amw**2*asb2)*
     .        alp/pi*NS_lamb(asb1/asb2,SMASS(K)/asb2)*
     .        Hsbotsbotr(K,2,1)*
     .        (delta1H(K)+delta2H(K)+delta3H(K)+
     .        delta4H(K)+delta5H(K))
      else
         qcdsb2H(K) = 0.D0
      endif
      ENDDO
      endif
c -------------------------------------------------------------------- c
c  sbottom2 --> A(K) + sbottom1
      DO K=1,2
         if(asb2.gt.(asb1+PMASS(K))) then
            sb2A(K)=g2s*amz**4/amw**2*Asbotsbotr(K,2,1)**2*
     .         NS_lamb(asb1/asb2,PMASS(K)/asb2)/(16.D0*pi*asb2)
         else
            sb2A(K)=0.D0
         endif
      ENDDO
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      DO K=1,2
        
      if(asb2.gt.(asb1+PMASS(K))) then
         alp = gs2/(4.D0*pi)
         delta1A(K) = dsqrt(2.D0)*amz**2*Asbotsbotr(K,2,1)*
     .        NS_gvirtgl(asb2,PMASS(K),asb1,lamv,amuv)
         delta2A(K) = 0.D0 
         delta3A(K) = NS_botneut1719(3+K,amuv)
         delta4A(K) =NS_dcounterneut(asb1,asb2,rmbc,theb,mgluino,amuv,
     .                          amuvdiv,lamv,2,K+3)
        delta5A(K) =NS_realcorr(PMASS(K),asb2,asb1,lamv,K+3,2,2,1,asb2)

        qcdsb2A(K) = g2s*amz**2/(24.D0*dsqrt(2.D0)*pi*amw**2*asb2)*
     .        alp/pi*NS_lamb(asb1/asb2,PMASS(K)/asb2)*Asbotsbotr(K,2,1)*
     .        (delta1A(K)+delta2A(K)+delta3A(K)+delta4A(K)+delta5A(K))
      else
         qcdsb2A(K) = 0.D0
      endif
      ENDDO
      endif
c -------------------------------------------------------------------- c
c  sbottom2 --> H- + stop1/2
      do i=1,2
         if(asb2.gt.(gmst(i)+cmass)) then
            sb2hcst(i)=g2s*amw**2*gctbr(i,2)**2*
     .           NS_lamb(gmst(i)/asb2,cmass/asb2)/(16.D0*pi*asb2)
         else
            sb2hcst(i)=0.D0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do nj=1,2
      if(asb2.gt.(gmst(nj)+cmass)) then
         alp = gs2/(4.D0*pi)
         delta21c = -dsqrt(2.D0)*amw**2*gctbr(nj,2)*
     .        NS_gvirtgl(gmst(nj),cmass,asb2,lamv,amuv)
         delta22c =  -dsqrt(2.D0)*amw**2*gctbr(3-nj,2)*
     .        1.D0/(gmst(3-nj)**2-gmst(nj)**2)*
     .        NS_gvirtmix(ast1,ast2,gmst(nj),mgluino,rmtc,thet,amuv)
     .        -dsqrt(2.D0)*amw**2*gctbr(nj,1)*
     .        1.D0/(asb1**2-asb2**2)*
     .        NS_gvirtmix(asb1,asb2,asb2,mgluino,rmbc,theb,amuv)
         delta23c = NS_sbotstop1719(amuv,2,nj)
         delta24c = NS_dcounterhc(asb2,rmbc,theb,2,gmst(nj),rmtc,
     .                         thet,nj,mgluino,amuv,amuvdiv,lamv,nj,2)
         delta25c = NS_realcorr(cmass,asb2,gmst(nj),lamv,6,0,nj,2,asb2)
         qcdsb2hcst(nj) = -g2s*amw**2/(24.D0*dsqrt(2.D0)*pi*amw**2*
     .        asb2)*alp/pi*NS_lamb(gmst(nj)/asb2,cmass/asb2)*
     .        gctbr(nj,2)*(delta21c+delta22c+delta23c+delta24c+delta25c)

      else
         qcdsb2hcst(nj) = 0.D0
      endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sbottom2 --> Z + sbottom1
      if(asb2.gt.(asb1+mz)) then
         sb2zbot=g2s/64.D0/pi/cw**2/mz**2*asb2**3*gzbb(2,1)**2*
     .           NS_lamb(asb1/asb2,mz/asb2)**3
      else
         sb2zbot=0.D0
      endif
c -- QCD corrections --
      if(flagqcd.eq.1.D0) then
      if(asb2.gt.(asb1+mz)) then
         alp = gs2/(4.D0*pi)
         del1 = -alp/3.D0/pi*gzbb(2,1)/2.D0/cw*
     .        NS_gluonvertex(asb2,asb1,mz,lamv,amuv)
         del2 = -alp/3.D0/pi/cw*NS_gluinoZvertex(asb2,asb1,mz,lamv,
     .        amuv,mgluino,amb,-1.D0/2.D0,-1.D0/3.D0,sw,theb)
         del3 = alp/pi*NS_wavefuncvertex(asb2,asb1,amb,amb,
     .     theb,theb,1.D0,2.D0,2,1,mgluino,lamv,amuv)
         del4 = alp/pi*NS_quarkmixZ(asb2,theb,-1.D0/2.D0,-1.D0/3.D0,
     .     asb1,asb2,amb,mgluino,amuv)
         del5 = NS_realgluonem(asb2,asb1,mz,lamv)
         qcdsb2zbot =  g2s/16.D0/pi/mz**2*asb2**3*(gzbb(2,1)/2.D0/cw)*
     .        NS_lamb(asb1/asb2,mz/asb2)**3*(2.D0*del1+2.D0*del2
     .        +2.D0*del3+2.D0*del4) + 
     .        g2s/3.D0/pi**2/asb2*alp*(gzbb(2,1)/(2.D0*cw))**2*del5
      else
         qcdsb2zbot = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c  sbottom1 --> W- + stop1/2
      do i=1,2
         if(asb1.gt.(gmst(i)+mw)) then
            sb1wst(i)=g2s/32.D0/pi/mw**2*asb1**3*gwtb(i,1)**2*
     .                NS_lamb(gmst(i)/asb1,mw/asb1)**3
         else
            sb1wst(i)=0.D0
         endif
      enddo
c -- QCD corrections --
      if(flagqcd.eq.1.D0) then
      do i=1,2
         if(asb1.gt.(gmst(i)+mw)) then
            alp = gs2/(4.D0*pi)
            del1 = -alp/3.D0/pi*gwtb(i,1)/dsqrt(2.D0)*
     .           NS_gluonvertex(asb1,gmst(i),mw,lamv,amuv)
            del2 = -dsqrt(2.D0)/3.D0*alp/pi*
     .           NS_gluinoWvertex(asb1,gmst(i),mw,lamv,amuv,mgluino,
     .           amb,amt,theb,thet,1,i)
            del3 = alp/pi*NS_wavefuncvertex(asb1,gmst(i),amb,amt,
     .           theb,thet,2.D0,2.D0,1,i,mgluino,lamv,amuv)
            del4 = alp/pi*NS_quarkmixW(asb1,theb,-1.D0/2.D0,-1.D0/3.D0,
     .             asb1,asb2,amb,gmst(i),thet,1.D0/2.D0,2.D0/3.D0,
     .             ast1,ast2,amt,1,i,mgluino,amuv)
            del5 = NS_realgluonem(asb1,gmst(i),mw,lamv)
            qcdsb1wst(i) = g2s/16.D0/pi/mw**2*asb1**3*
     .           (gwtb(i,1)/dsqrt(2.D0))*
     .           NS_lamb(gmst(i)/asb1,mw/asb1)**3*(2.D0*del1+2.D0*del2
     .           +2.D0*del3+2.D0*del4) + 
     .           g2s/3.D0/pi**2/asb1*alp*(gwtb(i,1)/dsqrt(2.D0))**2*del5
         else
            qcdsb1wst(i) = 0.D0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sbottom2 --> W- + stop1/2
      do i=1,2
         if(asb2.gt.(gmst(i)+mw)) then
            sb2wst(i)=g2s/32.D0/pi/mw**2*asb2**3*gwtb(i,2)**2*
     .                NS_lamb(gmst(i)/asb2,mw/asb2)**3
         else
            sb2wst(i)=0.D0
         endif
      enddo
c -- QCD corrections --
      if(flagqcd.eq.1.D0) then
      do i=1,2
         if(asb2.gt.(gmst(i)+mw)) then
            alp = gs2/(4.D0*pi)
            del1 = -alp/3.D0/pi*gwtb(i,2)/dsqrt(2.D0)*
     .           NS_gluonvertex(asb2,gmst(i),mw,lamv,amuv)
            del2 = -dsqrt(2.D0)/3.D0*alp/pi*
     .           NS_gluinoWvertex(asb2,gmst(i),mw,lamv,amuv,mgluino,
     .           amb,amt,theb,thet,2,i)
            del3 = alp/pi*NS_wavefuncvertex(asb2,gmst(i),amb,amt,
     .           theb,thet,2.D0,2.D0,2,i,mgluino,lamv,amuv)
            del4 = alp/pi*NS_quarkmixW(asb2,theb,-1.D0/2.D0,-1.D0/3.D0,
     .             asb1,asb2,amb,gmst(i),thet,1.D0/2.D0,2.D0/3.D0,
     .             ast1,ast2,amt,2,i,mgluino,amuv)
            del5 = NS_realgluonem(asb2,gmst(i),mw,lamv)
            qcdsb2wst(i) = g2s/16.D0/pi/mw**2*asb2**3*
     .           (gwtb(i,2)/dsqrt(2.D0))*
     .           NS_lamb(gmst(i)/asb2,mw/asb2)**3*(2.D0*del1+2.D0*del2
     .           +2.D0*del3+2.D0*del4) + 
     .           g2s/3.D0/pi**2/asb2*alp*(gwtb(i,2)/dsqrt(2.D0))**2*del5
         else
            qcdsb2wst(i) = 0.D0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c ----------------- 2-body decays and 2-body total widths ------------ c
c -------------------------------------------------------------------- c

      sbottot2lo(1)=sb1neutt(1)+sb1neutt(2)+sb1neutt(3)+sb1neutt(4)+
     .              sb1neutt(5)+
     .              sb1chart(1)+sb1chart(2)+sb1glui+sb1hcst(1)+
     .              sb1hcst(2)+sb1wst(1)+sb1wst(2)

      sbottot2lo(2)=sb2neutt(1)+sb2neutt(2)+sb2neutt(3)+sb2neutt(4)+
     .              sb2neutt(5)+
     .              sb2chart(1)+sb2chart(2)+sb2glui+sb2hcst(1)+
     .              sb2hcst(2)+sb2wst(1)+sb2wst(2)+
     .              sb2H(1)+sb2H(2)+sb2H(3)+sb2A(1)+sb2A(2)+
     .              sb2zbot

c UE: Resum if qcdcorr < -tree:
      do i=1,5
      If(qcdsb1neut(i).lt.-sb1neutt(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb1->b+chi0_",i
        qcdsb1neut(i)=resum(sb1neutt(i),qcdsb1neut(i))
      If(qcdsb2neut(i).lt.-sb2neutt(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb2->b+chi0_",i
        qcdsb2neut(i)=resum(sb2neutt(i),qcdsb2neut(i))
      enddo
      do i=1,2
      If(qcdsb1chart(i).lt.-sb1chart(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb1->t+char_",i
        qcdsb1chart(i)=resum(sb1chart(i),qcdsb1chart(i))
      If(qcdsb2chart(i).lt.-sb2chart(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb2->t+char_",i
        qcdsb2chart(i)=resum(sb2chart(i),qcdsb2chart(i))
      If(qcdsb1hcst(i).lt.-sb1hcst(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb1->hc+st_",i
        qcdsb1hcst(i)=resum(sb1hcst(i),qcdsb1hcst(i))
      If(qcdsb2hcst(i).lt.-sb2hcst(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb2->hc+st_",i
        qcdsb2hcst(i)=resum(sb2hcst(i),qcdsb2hcst(i))
      If(qcdsb1wst(i).lt.-sb1wst(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb1->W+st_",i
        qcdsb1wst(i)=resum(sb1wst(i),qcdsb1wst(i))
      If(qcdsb2wst(i).lt.-sb2wst(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb2->W+st_",i
        qcdsb2wst(i)=resum(sb2wst(i),qcdsb2wst(i))
      If(qcdsb2A(i).lt.-sb2A(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb2->sb1+A_",i
        qcdsb2A(i)=resum(sb2A(i),qcdsb2A(i))
      enddo
      If(qcdsb1glui.lt.-sb1glui) 
     .write(*,*)"Warning: large negative rad. corrs. to sb1->glui+b"
      qcdsb1glui=resum(sb1glui,qcdsb1glui)
      If(qcdsb2glui.lt.-sb2glui) 
     .write(*,*)"Warning: large negative rad. corrs. to sb2->glui+b"
      qcdsb2glui=resum(sb2glui,qcdsb2glui)
      do i=1,3
         If(qcdsb2H(i).lt.-sb2H(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sb2->sb1+H_",i
        qcdsb2H(i)=resum(sb2H(i),qcdsb2H(i))
      enddo
         If(qcdsb2zbot.lt.-sb2zbot) 
     .write(*,*)"Warning: large negative rad. corrs. to sb2->sb1+Z"
      qcdsb2zbot=resum(sb2zbot,qcdsb2zbot)
23     format(A,I1)
c UE: End resummation

      sbottot2nlo(1)=sbottot2lo(1) + qcdsb1neut(1)+qcdsb1neut(2)+
     .               qcdsb1neut(3) + qcdsb1neut(4)+
     .               qcdsb1neut(5) + qcdsb1chart(1)+
     .               qcdsb1chart(2)+ qcdsb1glui+qcdsb1hcst(1)+
     .               qcdsb1hcst(2)+qcdsb1wst(1)+qcdsb1wst(2)

      sbottot2nlo(2)=sbottot2lo(2) + qcdsb2neut(1)+qcdsb2neut(2)+
     .               qcdsb2neut(3)+qcdsb2neut(4)+qcdsb2chart(1)+
     .               qcdsb2chart(2)+qcdsb2glui+qcdsb2hcst(1)+
     .               qcdsb2hcst(2)+
     .               qcdsb2H(1)+qcdsb2H(2)+qcdsb2H(3)+
     .               qcdsb2A(2)+qcdsb2A(1)+
     .               qcdsb2zbot+qcdsb2wst(1)+qcdsb2wst(2)

      if(flagqcd.eq.0.D0) then
         do i=1,2
            sbottot2(i) = sbottot2lo(i)
         enddo
      elseif(flagqcd.eq.1.D0) then
         do i=1,2
            sbottot2(i) = sbottot2nlo(i)
         enddo
      endif

c -------------------------------------------------------------------- c
c --------------- 3-body decays and the 3-body widths ---------- c
c -------------------------------------------------------------------- c
      if(flagmulti.eq.1.D0) then

      call NS_xintegsbottom(xintegsb2sb1bb,xintegsb2sb1tt,
     .       xintegsb2sb1uu,xintegsb2sb1dd,xintegsb2sb1ee,
     .       xintegsb2sb1nunu,xintegsb2sb1tautau,xintegsb2sb1starbb,
     .       xintegsbtstsb,xintegsbtbstb,xintegsbtaustnu,xintegsbelstnu,
     .       xintegsbupstdow,xintegsbstau,xintegsbsntau,xintegsbsel,
     .       xintegsbsnel)

      do i=1,2
            sbottotmulti(i) = 0.D0
            do j=1,2
               sbottotmulti(i) = sbottotmulti(i)+xintegsbstau(i,j)+
     .              xintegsbsntau(i,j)+2.D0*xintegsbsel(i,j)+
     .              xintegsbtstsb(i,j)+xintegsbtbstb(i,j)+
     .              xintegsbtaustnu(i,j)+2.D0*xintegsbelstnu(i,j)+
     .              2.D0*xintegsbupstdow(i,j)
            enddo
            sbottotmulti(i) = sbottotmulti(i)+2.D0*xintegsbsnel(i)
            if(i.eq.2) then
               sbottotmulti(i) = sbottotmulti(i)+xintegsb2sb1bb+
     .              xintegsb2sb1tt+2.D0*xintegsb2sb1uu+
     .              2.D0*xintegsb2sb1dd+2.D0*xintegsb2sb1ee+
     .              3.D0*xintegsb2sb1nunu+xintegsb2sb1tautau+
     .              xintegsb2sb1starbb
            endif

      enddo
      endif
c Check if numerically relevant:
      do i =1,2
      if (sbottotmulti(i).lt.multilim*sbottot2(i))Then
         sbottotmulti(i)=0.d0
      endif
      enddo

c ------------------------  total widths -------------------------- c
      do i =1,2
      sbottot(i)=sbottot2(i)+sbottotmulti(i)
      enddo
c ---------------------  sbottom branching ratios ----------------- c
      if(flagqcd.eq.1.D0) then
         do i=1,5
            sb1neutt(i) = sb1neutt(i)+qcdsb1neut(i)
            sb2neutt(i) = sb2neutt(i)+qcdsb2neut(i)
         enddo
         do i=1,2
            sb1chart(i) = sb1chart(i) + qcdsb1chart(i)
            sb2chart(i) = sb2chart(i) + qcdsb2chart(i)
            sb1hcst(i)  = sb1hcst(i) + qcdsb1hcst(i)
            sb2hcst(i)  = sb2hcst(i) + qcdsb2hcst(i)
            sb1wst(i)   = sb1wst(i) + qcdsb1wst(i)
            sb2wst(i)   = sb2wst(i) + qcdsb2wst(i)
         enddo
         sb1glui = sb1glui + qcdsb1glui
         sb2glui = sb2glui + qcdsb2glui
         sb2zbot = sb2zbot + qcdsb2zbot
         DO i=1,3
            sb2H(i) = sb2H(i)+qcdsb2H(i)
         ENDDO
         DO i=1,2
            sb2A(i) = sb2A(i)+qcdsb2A(i)
         ENDDO
      endif

      do i=1,5
         brsb1neutt(i)=sb1neutt(i)/sbottot(1)
         brsb2neutt(i)=sb2neutt(i)/sbottot(2)
      enddo
      do i=1,2
         brsb1chart(i) = sb1chart(i)/sbottot(1)
         brsb2chart(i) = sb2chart(i)/sbottot(2)
         brsb1hcst(i)  = sb1hcst(i)/sbottot(1)
         brsb2hcst(i)  = sb2hcst(i)/sbottot(2)
         brsb1wst(i)   = sb1wst(i)/sbottot(1)
         brsb2wst(i)   = sb2wst(i)/sbottot(2)
      enddo
      brsb1glui = sb1glui/sbottot(1)
      brsb2glui = sb2glui/sbottot(2)
         DO i=1,3
            brsb2H(i) = sb2H(i)/sbottot(2)
         ENDDO
         DO i=1,2
            brsb2A(i) = sb2A(i)/sbottot(2)
         ENDDO
      brsb2zbot = sb2zbot/sbottot(2)

c -- 3-body decays --

      do i=1,2
         if (sbottotmulti(i).ne.0d0)Then
            do j=1,2
               brsbstau(i,j)    = xintegsbstau(i,j)/sbottot(i)
               brsbsntau(i,j)   = xintegsbsntau(i,j)/sbottot(i)
               brsbsel(i,j)     = xintegsbsel(i,j)/sbottot(i)
               brsbtstsb(i,j)   = xintegsbtstsb(i,j)/sbottot(i)
               brsbtbstb(i,j)   = xintegsbtbstb(i,j)/sbottot(i)
               brsbtaustnu(i,j) = xintegsbtaustnu(i,j)/sbottot(i)
               brsbelstnu(i,j)  = xintegsbelstnu(i,j)/sbottot(i)
               brsbupstdow(i,j) = xintegsbupstdow(i,j)/sbottot(i)
            enddo
            brsbsnel(i) = xintegsbsnel(i)/sbottot(i)
         endif
      enddo

      if (sbottotmulti(2).ne.0d0)Then
         brsb2sb1bb     = xintegsb2sb1bb/sbottot(2)
         brsb2sb1starbb = xintegsb2sb1starbb/sbottot(2)
         brsb2sb1tt     = xintegsb2sb1tt/sbottot(2)
         brsb2sb1uu     = xintegsb2sb1uu/sbottot(2)
         brsb2sb1dd     = xintegsb2sb1dd/sbottot(2)
         brsb2sb1ee     = xintegsb2sb1ee/sbottot(2)
         brsb2sb1nunu   = xintegsb2sb1nunu/sbottot(2)
         brsb2sb1tautau = xintegsb2sb1tautau/sbottot(2)
c         write(*,*)'sbotbr',xintegsb2sb1bb,xintegsb2sb1starbb,
c     .xintegsb2sb1tt, xintegsb2sb1uu,xintegsb2sb1dd,xintegsb2sb1ee
c     .,xintegsb2sb1nunu,xintegsb2sb1tautau,sbottot(2),sbottotmulti(2)
      endif

      END
c ==================================================================== c
c                       The sbottom 3-body decays                      c
c ==================================================================== c
      SUBROUTINE NS_xintegsbottom(xintegsb2sb1bb,xintegsb2sb1tt,
     .     xintegsb2sb1uu,xintegsb2sb1dd,xintegsb2sb1ee,
     .     xintegsb2sb1nunu,xintegsb2sb1tautau,xintegsb2sb1starbb,
     .     xintegsbtstsb,xintegsbtbstb,xintegsbtaustnu,xintegsbelstnu,
     .     xintegsbupstdow,xintegsbstau,xintegsbsntau,xintegsbsel,
     .     xintegsbsnel)
*
      IMPLICIT NONE
      INTEGER nx1t,ny1t,ni,nj
      DOUBLE PRECISION xintegsbtstsb(2,2),xintegsbtbstb(2,2),
     .          xintegsbtaustnu(2,2),xintegsbelstnu(2,2),
     .          xintegsbupstdow(2,2),xintegsbstau(2,2),
     .          xintegsbsntau(2,2),xintegsbsel(2,2),xintegsbsnel(2)
      DOUBLE PRECISION xintegsb2sb1bb,xintegsb2sb1tt,
     .     xintegsb2sb1uu,xintegsb2sb1dd,xintegsb2sb1ee,
     .     xintegsb2sb1nunu,xintegsb2sb1tautau,xintegsb2sb1starbb
      DOUBLE PRECISION 
     .gmsb(2),gmst(2),gmstau(2),gmsel(2),gmsne(2),gmsnt(2)
      DOUBLE PRECISION NS_ay,NS_by,NS_ax,NS_bx
      DOUBLE PRECISION NS_sb2sb1bb,NS_sb2sb1tt,NS_sb2sb1uu,NS_sb2sb1dd,
     .         NS_sb2sb1ee,NS_sb2sb1nunu,NS_sb2sb1tautau,
     .         NS_sb2sb1starbb,NS_sbtststarb,NS_sbtbstb,NS_sbtaustnu,
     .         NS_sbelstnu,NS_sbtnustau,NS_sbtsnutau,NS_sbtnusel,
     .         NS_sbtsnuel
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION asne2,asntau2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION XMU1,XMU2,XMU3,sum 
*      
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_rhsneutr/asne2,asntau2
      COMMON/NS_nx1/nx1t,ny1t
      COMMON/NS_indices/ni,nj
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_sb2sb1bb,NS_sb2sb1tt,NS_sb2sb1uu,NS_sb2sb1dd,
     .         NS_sb2sb1ee,NS_sb2sb1nunu,NS_sb2sb1tautau,
     .         NS_sb2sb1starbb,NS_sbtststarb,NS_sbtbstb,NS_sbtaustnu,
     .         NS_sbelstnu,NS_sbtnustau,NS_sbtsnutau,NS_sbtnusel,
     .         NS_sbtsnuel

      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2

      gmstau(1) = astau1
      gmstau(2) = astau2
      gmsel(1)  = ase1
      gmsel(2)  = ase2
      gmsne(1)  = asne1
      gmsne(2)  = asne2
      gmsnt(1)  = asntau1
      gmsnt(2)  = asntau2
c -------------------------------------------------------------------- c
c ------------------------- t stau neutrino_tau ---------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2,1
         do nj=1,2,1
            xmu1=amt**2/gmsb(ni)**2
            xmu2=0.D0
            xmu3=gmstau(nj)**2/gmsb(ni)**2

            if(gmsb(ni).gt.(amt+gmstau(nj))) then
               call NS_integ2(NS_sbtnustau,NS_ax,NS_bx,NS_ay,NS_by,
     .              xmu1,xmu2,xmu3,nx1t,ny1t,sum)
               xintegsbstau(ni,nj)=g2s**2/32.D0/(2.D0*pi)**3*gmsb(ni)*
     .              sum
            else
               xintegsbstau(ni,nj)=0.D0
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c ------------------------- t sneutrino_tau tau ---------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2,1
         do nj=1,2,1
            xmu1=amt**2/gmsb(ni)**2
            xmu2=amtau**2/gmsb(ni)**2
            xmu3=gmsnt(nj)**2/gmsb(ni)**2

            if(gmsb(ni).gt.(gmsnt(nj)+amt+amtau)) then
               call NS_integ2(NS_sbtsnutau,NS_ax,NS_bx,NS_ay,NS_by,
     .              xmu1,xmu2,xmu3,nx1t,ny1t,sum)
               xintegsbsntau(ni,nj)=g2s**2/32.D0/(2.D0*pi)**3*
     .              gmsb(ni)*sum
            else
               xintegsbsntau(ni,nj)=0.D0
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c ----------------------- t selectron neutrino_e --------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2,1
         do nj=1,2,1
            xmu1=amt**2/gmsb(ni)**2
            xmu2=0.D0
            xmu3=gmsel(nj)**2/gmsb(ni)**2

            if(gmsb(ni).gt.(gmsel(nj)+amt)) then
               call NS_integ2(NS_sbtnusel,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .              xmu2,xmu3,nx1t,ny1t,sum)
               xintegsbsel(ni,nj)=g2s**2/32.D0/(2.D0*pi)**3*gmsb(ni)*
     .              sum
            else
               xintegsbsel(ni,nj)=0.D0
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c ----------------------- t sneutrino_e electron --------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2,1
         xmu1=amt**2/gmsb(ni)**2
         xmu2=0.D0
         xmu3=gmsne(1)**2/gmsb(ni)**2

         if(gmsb(ni).gt.(gmsne(1)+amt)) then
            call NS_integ2(NS_sbtsnuel,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .           xmu2,xmu3,nx1t,ny1t,sum)
            xintegsbsnel(ni)=g2s**2/32.D0/(2.D0*pi)**3*gmsb(ni)*sum
         else
            xintegsbsnel(ni)=0.D0
         endif
      end do
c -------------------------------------------------------------------- c
c ------------------------- stop_1/2* bottom top --------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2,1
         do nj=1,2,1
            xmu1=amt**2/gmsb(ni)**2
            xmu2=amb**2/gmsb(ni)**2
            xmu3=gmst(nj)**2/gmsb(ni)**2

            if(gmsb(ni).gt.(amt+amb+gmst(nj))) then
               call NS_integ2(NS_sbtststarb,NS_ax,NS_bx,NS_ay,NS_by,
     .              xmu1,xmu2,xmu3,nx1t,ny1t,sum)
               xintegsbtstsb(ni,nj)=1.D0/32.D0/(2.D0*pi)**3*gmsb(ni)*sum
            else
               xintegsbtstsb(ni,nj)=0.D0
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c ----------------------- stop_1/2 bottom topbar --------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2,1
         do nj=1,2,1
            xmu1=amt**2/gmsb(ni)**2
            xmu2=amb**2/gmsb(ni)**2
            xmu3=gmst(nj)**2/gmsb(ni)**2

            if(gmsb(ni).gt.(amt+amb+gmst(nj))) then
               call NS_integ2(NS_sbtbstb,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .              xmu2,xmu3,nx1t,ny1t,sum)
               xintegsbtbstb(ni,nj)=1.D0/32.D0/(2.D0*pi)**3*gmsb(ni)*
     .              sum
            else
               xintegsbtbstb(ni,nj)=0.D0
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c ---------------------- stop_1/2 tau- nu_taubar --------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2,1
         do nj=1,2,1
            xmu1=amtau**2/gmsb(ni)**2
            xmu2=0.D0
            xmu3=gmst(nj)**2/gmsb(ni)**2

            if(gmsb(ni).gt.(amtau+gmst(nj))) then
               call NS_integ2(NS_sbtaustnu,NS_ax,NS_bx,NS_ay,NS_by,
     .              xmu1,xmu2,xmu3,nx1t,ny1t,sum)
               xintegsbtaustnu(ni,nj)=g2s**2/32.D0/(2.D0*pi)**3*
     .              gmsb(ni)*sum
            else
               xintegsbtaustnu(ni,nj)=0.D0
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c --------------------- stop_1/2 electron- nu_ebar ------------------- c
c -------------------------------------------------------------------- c

      do ni=1,2,1
         do nj=1,2,1
            xmu1=0.D0
            xmu2=0.D0
            xmu3=gmst(nj)**2/gmsb(ni)**2
            
            if(gmsb(ni).gt.gmst(nj)) then
               call NS_integ2(NS_sbelstnu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .              xmu2,xmu3,nx1t,ny1t,sum)
               xintegsbelstnu(ni,nj)=g2s**2/32.D0/(2.D0*pi)**3*
     .              gmsb(ni)*sum
            else
               xintegsbelstnu(ni,nj)=0.D0
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c ------------------------- stop_1/2 upbar down ---------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2,1
         do nj=1,2,1
            if(gmsb(ni).gt.gmst(nj)) then
               xintegsbupstdow(ni,nj)=3.D0*xintegsbelstnu(ni,nj)
            else
               xintegsbupstdow(ni,nj)=0.D0
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c --------------------------- sbottom1* b b -------------------------- c
c -------------------------------------------------------------------- c
     
      xmu1=amb**2/gmsb(2)**2
      xmu2=amb**2/gmsb(2)**2
      xmu3=gmsb(1)**2/gmsb(2)**2
      

      if(gmsb(2).gt.(gmsb(1)+2.D0*amb)) then
         call NS_integ2(NS_sb2sb1starbb,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsb2sb1starbb=1.D0/32.D0/(2.D0*pi)**3*gmsb(2)*sum
      else 
         xintegsb2sb1starbb=0.D0
      endif
      
c -------------------------------------------------------------------- c
c --------------------------- sbottom1 b bbar ------------------------ c
c -------------------------------------------------------------------- c
      xmu1=amb**2/gmsb(2)**2
      xmu2=amb**2/gmsb(2)**2
      xmu3=gmsb(1)**2/gmsb(2)**2

      if(gmsb(2).gt.(gmsb(1)+2.D0*amb)) then
         call NS_integ2(NS_sb2sb1bb,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .        xmu3,nx1t,ny1t,sum)
         xintegsb2sb1bb=1.D0/32.D0/(2.D0*pi)**3*gmsb(2)*sum
      else 
         xintegsb2sb1bb=0.D0
      endif
      
c -------------------------------------------------------------------- c
c --------------------------- sbottom1 t tbar ------------------------ c
c -------------------------------------------------------------------- c
      xmu1=amt**2/gmsb(2)**2
      xmu2=amt**2/gmsb(2)**2
      xmu3=gmsb(1)**2/gmsb(2)**2

      if((gmsb(1)+2.D0*amt).lt.gmsb(2)) then
         call NS_integ2(NS_sb2sb1tt,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .        xmu3,nx1t,ny1t,sum)
         xintegsb2sb1tt=g2s**2/32.D0/(2.D0*pi)**3*gmsb(2)*sum*3.D0
      else 
         xintegsb2sb1tt=0.D0
      endif
      
c -------------------------------------------------------------------- c
c --------------------------- sbottom1 up upbar ---------------------- c
c -------------------------------------------------------------------- c
      xmu1=0.D0
      xmu2=0.D0
      xmu3=gmsb(1)**2/gmsb(2)**2

      if(gmsb(2).gt.gmsb(1)) then
         call NS_integ2(NS_sb2sb1uu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .        xmu3,nx1t,ny1t,sum)
         xintegsb2sb1uu=g2s**2/32.D0/(2.D0*pi)**3*gmsb(2)*sum*3.D0
      else 
         xintegsb2sb1uu=0.D0
      endif
c -------------------------------------------------------------------- c
c ------------------------- sbottom1 down downbar -------------------- c
c -------------------------------------------------------------------- c
      xmu1=0.D0
      xmu2=0.D0
      xmu3=gmsb(1)**2/gmsb(2)**2

      if(gmsb(2).gt.gmsb(1)) then
         call NS_integ2(NS_sb2sb1dd,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .        xmu3,nx1t,ny1t,sum)
         xintegsb2sb1dd=g2s**2/32.D0/(2.D0*pi)**3*gmsb(2)*sum*3.D0
      else 
         xintegsb2sb1dd=0.D0
      endif
c -------------------------------------------------------------------- c
c ----------------------------- sbottom1 e+ e- ----------------------- c
c -------------------------------------------------------------------- c
      xmu1=0.D0
      xmu2=0.D0
      xmu3=gmsb(1)**2/gmsb(2)**2

      if(gmsb(2).gt.gmsb(1)) then
         call NS_integ2(NS_sb2sb1ee,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .        xmu3,nx1t,ny1t,sum)
         xintegsb2sb1ee=g2s**2/32.D0/(2.D0*pi)**3*gmsb(2)*sum*1.D0
      else 
         xintegsb2sb1ee=0.D0
      endif
c -------------------------------------------------------------------- c
c -------------------------- sbottom1 nu nubar ----------------------- c
c -------------------------------------------------------------------- c
      xmu1=0.D0
      xmu2=0.D0
      xmu3=gmsb(1)**2/gmsb(2)**2

      if(gmsb(2).gt.gmsb(1)) then
         call NS_integ2(NS_sb2sb1nunu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsb2sb1nunu=g2s**2/32.D0/(2.D0*pi)**3*gmsb(2)*sum*1.D0
      else 
         xintegsb2sb1nunu=0.D0
      endif
c -------------------------------------------------------------------- c
c ------------------------ sbottom1 tau+ tau- ------------------------ c
c -------------------------------------------------------------------- c
      xmu1=amtau**2/gmsb(2)**2
      xmu2=amtau**2/gmsb(2)**2
      xmu3=gmsb(1)**2/gmsb(2)**2

      if(gmsb(2).gt.(gmsb(1)+2.D0*amtau)) then
         call NS_integ2(NS_sb2sb1tautau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsb2sb1tautau=g2s**2/32.D0/(2.D0*pi)**3*gmsb(2)*sum*1.D0
      else 
         xintegsb2sb1tautau=0.D0
      endif
      
      end
c ==================================================================== c
c ========================= t stau neutrino_tau ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sbtnustau(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k,l
      DOUBLE PRECISION gmsb(2),xmuchar(2),gmstau(2),
     .          dchi(2),xmustau(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2)
      DOUBLE PRECISION ale(2,2),alto(2,2),alsne(2,2),blsne(2,2),
     .         alsnt(2,2),blsnt(2,2),blto(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION xmut,x1,x2
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_indices/ni,nj
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/NS_coup5/ale,alto,alsne,blsne,alsnt,blsnt
      COMMON/NS_charsbottop/alsbot,aksbot
*
      gmsb(1)=asb1
      gmsb(2)=asb2
      gmstau(1)=astau1
      gmstau(2)=astau2

      xmuchar(1) = amchar(1)**2/gmsb(ni)**2
      xmuchar(2) = amchar(2)**2/gmsb(ni)**2
      xmustau(1) = gmstau(1)**2/gmsb(ni)**2
      xmustau(2) = gmstau(2)**2/gmsb(ni)**2
      xmut       = amt**2/gmsb(ni)**2

      dchi(1)=1.D0-x1-xmuchar(1)+xmut
      dchi(2)=1.D0-x1-xmuchar(2)+xmut

      do i=1,2,1
         blto(1,i)=0.D0
         blto(2,i)=0.D0
      end do

      NS_sbtnustau=0.D0

      do k=1,2,1
         do l=1,2,1
        if ((amchar(k)+amt).gt.gmsb(ni).and.
     .(amchar(l)+amt).gt.gmsb(ni))Then
            NS_sbtnustau=NS_sbtnustau+1.D0/dchi(k)/dchi(l)*(
     .        (alsbot(ni,k)*alsbot(ni,l)*blto(nj,k)*blto(nj,l)
     .        +aksbot(ni,k)*aksbot(ni,l)*alto(nj,k)*alto(nj,l))*
     .        xmchar(k)*xmchar(l)/gmsb(ni)**2*(x1+x2-1.D0+xmustau(nj)
     .         -xmut)
     .       +(alsbot(ni,k)*alsbot(ni,l)*alto(nj,k)*alto(nj,l)
     .        +aksbot(ni,k)*aksbot(ni,l)*blto(nj,k)*blto(nj,l))*
     .        ((1.D0-x1)*(1.D0-x2)-xmustau(nj)+xmut*(xmustau(nj)+x1-x2
     .         -xmut)) 
     .       +(alsbot(ni,k)*aksbot(ni,l)*alto(nj,k)*alto(nj,l)
     .        +aksbot(ni,k)*alsbot(ni,l)*blto(nj,k)*blto(nj,l))*
     .        dsqrt(xmut)*xmchar(l)/gmsb(ni)*(-1.D0-xmut+xmustau(nj)+x1)
     .       +(alsbot(ni,k)*aksbot(ni,l)*blto(nj,k)*blto(nj,l)
     .        +aksbot(ni,k)*alsbot(ni,l)*alto(nj,k)*alto(nj,l))*
     .       dsqrt(xmut)*xmchar(k)/gmsb(ni)*(-1.D0-xmut+xmustau(nj)+x1))
            endif
         end do
      end do

      end
c ==================================================================== c
c ========================= t sneutrino_tau tau ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sbtsnutau(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,k,l
      DOUBLE PRECISION gmsb(2),xmuchar(2),dchi(2),
     .          xmusn(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2)
      DOUBLE PRECISION ale(2,2),alto(2,2),alsne(2,2),blsne(2,2)
     .          ,alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION asne2,asntau2
      DOUBLE PRECISION xmut,xmutau,x1,x2
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_indices/ni,nj
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_rhsneutr/asne2,asntau2
      COMMON/NS_coup5/ale,alto,alsne,blsne,alsnt,blsnt
      COMMON/NS_charsbottop/alsbot,aksbot
*
      gmsb(1)=asb1
      gmsb(2)=asb2
      xmuchar(1)= amchar(1)**2/gmsb(ni)**2
      xmuchar(2)= amchar(2)**2/gmsb(ni)**2
      xmusn(1)  = asntau1**2/gmsb(ni)**2
      xmusn(2)  = asntau2**2/gmsb(ni)**2
      xmut      = amt**2/gmsb(ni)**2
      xmutau    = amtau**2/gmsb(ni)**2

      dchi(1)=1.D0-x1-xmuchar(1)+xmut
      dchi(2)=1.D0-x1-xmuchar(2)+xmut

      NS_sbtsnutau=0.D0

      do k=1,2,1
         do l=1,2,1
       if ((amchar(k)+amt).gt.gmsb(ni).and.
     .(amchar(l)+amt).gt.gmsb(ni))Then
            NS_sbtsnutau=NS_sbtsnutau+1.D0/dchi(k)/dchi(l)*(
     .        (alsbot(ni,k)*alsbot(ni,l)*alsnt(nj,k)*alsnt(nj,l)+
     .         aksbot(ni,k)*aksbot(ni,l)*blsnt(nj,k)*blsnt(nj,l))*
     .      xmchar(k)*xmchar(l)/gmsb(ni)**2*(x1+x2-1.D0+xmusn(nj)
     .        -xmut-xmutau)
     .       +(alsbot(ni,k)*alsbot(ni,l)*blsnt(nj,k)*blsnt(nj,l)+
     .         aksbot(ni,k)*aksbot(ni,l)*alsnt(nj,k)*alsnt(nj,l))*
     .        ((1.D0-x1)*(1.D0-x2)-xmusn(nj)+xmut*(xmusn(nj)+x1-x2
     .        -xmut-xmutau)+xmutau)
     .       +(alsbot(ni,k)*alsbot(ni,l)*alsnt(nj,l)*blsnt(nj,k)+
     .         aksbot(ni,k)*aksbot(ni,l)*blsnt(nj,l)*alsnt(nj,k))*
     .        dsqrt(xmutau)*xmchar(l)/gmsb(ni)*(-2.D0*xmut+x1)
     .       +(alsbot(ni,k)*aksbot(ni,l)*blsnt(nj,k)*blsnt(nj,l)+
     .         aksbot(ni,k)*alsbot(ni,l)*alsnt(nj,k)*alsnt(nj,l))*
     .        dsqrt(xmut)*xmchar(l)/gmsb(ni)*
     .       (-1.D0-xmut-xmutau+xmusn(nj)+x1)
     .       +(alsbot(ni,k)*alsbot(ni,l)*alsnt(nj,k)*blsnt(nj,l)+
     .         aksbot(ni,k)*aksbot(ni,l)*blsnt(nj,k)*alsnt(nj,l))*
     .        dsqrt(xmutau)*xmchar(k)/gmsb(ni)*(-2.D0*xmut+x1)
     .       +(alsbot(ni,k)*aksbot(ni,l)*alsnt(nj,k)*alsnt(nj,l)+
     .         aksbot(ni,k)*alsbot(ni,l)*blsnt(nj,k)*blsnt(nj,l))*
     .        dsqrt(xmut)*xmchar(k)/gmsb(ni)
     .       *(-1.D0-xmut-xmutau+xmusn(nj)+x1)
     .       +(alsbot(ni,k)*aksbot(ni,l)*blsnt(nj,k)*alsnt(nj,l)+
     .         aksbot(ni,k)*alsbot(ni,l)*alsnt(nj,k)*blsnt(nj,l))*
     .        dsqrt(xmut*xmutau)*(-2.D0-2.D0*xmut+2.D0*x1)
     .       +(alsbot(ni,k)*aksbot(ni,l)*alsnt(nj,k)*blsnt(nj,l)+
     .         aksbot(ni,k)*alsbot(ni,l)*blsnt(nj,k)*alsnt(nj,l))*
     .        dsqrt(xmut*xmutau)*
     .       *xmchar(k)*xmchar(l)/gmsb(ni)**2*(-2.D0) )
            endif
         end do
      end do

      end
c ==================================================================== c
c ======================== t selectron neutrino_e ==================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sbtnusel(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k,l
      DOUBLE PRECISION gmsb(2),xmuchar(2),gmsel(2),
     .          dchi(2),xmusel(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2)
      DOUBLE PRECISION ale(2,2),altau(2,2),alsne(2,2),blsne(2,2),
     .          alsnt(2,2),blsnt(2,2),ble(2,2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION xmut,x1,x2
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_indices/ni,nj
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_coup5/ale,altau,alsne,blsne,alsnt,blsnt
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/NS_charsbottop/alsbot,aksbot
*
      gmsb(1)=asb1
      gmsb(2)=asb2
      gmsel(1)=ase1
      gmsel(2)=ase2

      xmuchar(1) = amchar(1)**2/gmsb(ni)**2
      xmuchar(2) = amchar(2)**2/gmsb(ni)**2
      xmusel(1)  = gmsel(1)**2/gmsb(ni)**2
      xmusel(2)  = gmsel(2)**2/gmsb(ni)**2
      xmut       = amt**2/gmsb(ni)**2

      dchi(1)=1.D0-x1-xmuchar(1)+xmut
      dchi(2)=1.D0-x1-xmuchar(2)+xmut

      do i=1,2,1
         ble(1,i)=0.D0
         ble(2,i)=0.D0
      end do

      NS_sbtnusel=0.D0

      do k=1,2,1
         do l=1,2,1
      if ((amchar(k)+amt).gt.gmsb(ni).and.
     .(amchar(l)+amt).gt.gmsb(ni))Then
            NS_sbtnusel=NS_sbtnusel+1.D0/dchi(k)/dchi(l)*(
     .        (alsbot(ni,k)*alsbot(ni,l)*ble(nj,k)*ble(nj,l)
     .        +aksbot(ni,k)*aksbot(ni,l)*ale(nj,k)*ale(nj,l))*
     .       xmchar(k)*xmchar(l)/gmsb(ni)**2*(x1+x2-1.D0+xmusel(nj)
     .         -xmut)
     .       +(alsbot(ni,k)*alsbot(ni,l)*ale(nj,k)*ale(nj,l)
     .        +aksbot(ni,k)*aksbot(ni,l)*ble(nj,k)*ble(nj,l))*
     .        ((1.D0-x1)*(1.D0-x2)-xmusel(nj)+xmut*(xmusel(nj)+x1-x2
     .         -xmut)) 
     .       +(alsbot(ni,k)*aksbot(ni,l)*ale(nj,k)*ale(nj,l)
     .        +aksbot(ni,k)*alsbot(ni,l)*ble(nj,k)*ble(nj,l))*
     .        dsqrt(xmut)*xmchar(l)/gmsb(ni)*(-1.D0-xmut+xmusel(nj)+x1)
     .       +(alsbot(ni,k)*aksbot(ni,l)*ble(nj,k)*ble(nj,l)
     .        +aksbot(ni,k)*alsbot(ni,l)*ale(nj,k)*ale(nj,l))*
     .       dsqrt(xmut)*xmchar(k)/gmsb(ni)*(-1.D0-xmut+xmusel(nj)+x1) )
            endif
         end do
      end do

      end
c ==================================================================== c
c ======================= t sneutrino_e electron ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sbtsnuel(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,L,K
      DOUBLE PRECISION gmsb(2),xmuchar(2),gmsn(4),
     .          dchi(2),xmusn(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2)
      DOUBLE PRECISION ale(2,2),altau(2,2),alsne(2,2),blsne(2,2)
     .          ,alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION asne2,asntau2
      DOUBLE PRECISION xmut,x1,x2
*          
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_indices/ni,nj
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_coup5/ale,altau,alsne,blsne,alsnt,blsnt
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_rhsneutr/asne2,asntau2
      COMMON/NS_charsbottop/alsbot,aksbot
*
      gmsb(1)=asb1
      gmsb(2)=asb2
      gmsn(1)=asne1
      gmsn(2)=asne2

      xmuchar(1)= amchar(1)**2/gmsb(ni)**2
      xmuchar(2)= amchar(2)**2/gmsb(ni)**2
      xmusn(1)  = gmsn(1)**2/gmsb(ni)**2
      xmusn(2)  = gmsn(2)**2/gmsb(ni)**2
      xmut      = amt**2/gmsb(ni)**2

      dchi(1)=1.D0-x1-xmuchar(1)+xmut
      dchi(2)=1.D0-x1-xmuchar(2)+xmut

      NS_sbtsnuel=0.D0

      do k=1,2,1
         do l=1,2,1
      if ((amchar(k)+amt).gt.gmsb(ni).and.
     .(amchar(l)+amt).gt.gmsb(ni))Then
            NS_sbtsnuel=NS_sbtsnuel+1.D0/dchi(k)/dchi(l)*(
     .        (alsbot(ni,k)*alsbot(ni,l)*alsne(1,k)*alsne(1,l)+
     .         aksbot(ni,k)*aksbot(ni,l)*blsne(1,k)*blsne(1,l))*
     .        xmchar(k)*xmchar(l)/gmsb(ni)**2*(x1+x2-1.D0+xmusn(1)
     .        -xmut)
     .       +(alsbot(ni,k)*alsbot(ni,l)*blsne(1,k)*blsne(1,l)+
     .         aksbot(ni,k)*aksbot(ni,l)*alsne(1,k)*alsne(1,l))*
     .        ((1.D0-x1)*(1.D0-x2)-xmusn(1)+xmut*(xmusn(1)+x1-x2
     .        -xmut))
     .       +(alsbot(ni,k)*aksbot(ni,l)*blsne(1,k)*blsne(1,l)+
     .         aksbot(ni,k)*alsbot(ni,l)*alsne(1,k)*alsne(1,l))*
     .        dsqrt(xmut)*xmchar(l)/gmsb(ni)*(-1.D0-xmut+xmusn(1)+x1)
     .       +(alsbot(ni,k)*aksbot(ni,l)*alsne(1,k)*alsne(1,l)+
     .         aksbot(ni,k)*alsbot(ni,l)*blsne(1,k)*blsne(1,l))*
     .        dsqrt(xmut)*xmchar(k)/gmsb(ni)*(-1.D0-xmut+xmusn(1)+x1) )
            endif
         end do
      end do

      end
c ==================================================================== c
c ========================= bottom stop_1/2* top ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sbtststarb(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,l,k
      DOUBLE PRECISION neutneut
      DOUBLE PRECISION gmst(2),xmuchar(2),xmuneut(5),gmsb(2),
     .          dchi(2),xmust(2),dneut(5)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2),
     .          abot(2,5),bbot(2,5),atopr(2,5),btopr(2,5)
      DOUBLE PRECISION
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      DOUBLE PRECISION xmut,XMUB,xmugl,x1,x2,dgl,charchar,gluiglui
     .,charneut
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/NS_indices/ni,nj
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
*
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_neutstoptop/atopr,btopr

      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2
      do k=1,5,1 
         xmuneut(k) = amneut(k)**2/gmsb(ni)**2
      end do
      xmuchar(1) = amchar(1)**2/gmsb(ni)**2
      xmuchar(2) = amchar(2)**2/gmsb(ni)**2
      xmust(1)   = gmst(1)**2/gmsb(ni)**2
      xmust(2)   = gmst(2)**2/gmsb(ni)**2
      xmut       = amt**2/gmsb(ni)**2
      xmub       = amb**2/gmsb(ni)**2
      xmugl      = mgluino**2/gmsb(ni)**2

      dchi(1)=1.D0-x1-xmuchar(1)+xmut
      dchi(2)=1.D0-x1-xmuchar(2)+xmut

      do k=1,5,1
         dneut(k)=1.D0-x2-xmuneut(k)+xmub
      end do

      dgl = 1.D0-x2-xmugl+xmub
c -------------------------------------------------------------------- c
c                           chargino exchange
c -------------------------------------------------------------------- c
      charchar=0.D0

      do k=1,2,1
         do l=1,2,1
      if ((amchar(k)+amt).gt.gmsb(ni).and.
     .(amchar(l)+amt).gt.gmsb(ni))Then
            charchar=charchar+3.D0*g2s**2/dchi(k)/dchi(l)*(
     .           xmchar(k)*xmchar(l)/gmsb(ni)**2*
     .           (-1.D0+x1+x2-xmut-xmub+xmust(nj))*
     .           (alsbot(ni,k)*alsbot(ni,l)*alstor(nj,k)*alstor(nj,l)+
     .            aksbot(ni,k)*aksbot(ni,l)*akstor(nj,k)*akstor(nj,l))
     .           +(xmut*(-x2+xmust(nj)-xmut+x1-xmub)-x2+x1*x2-x1+1.D0
     .             -xmust(nj)+xmub)*
     .           (alsbot(ni,k)*alsbot(ni,l)*akstor(nj,k)*akstor(nj,l)+
     .            aksbot(ni,k)*aksbot(ni,l)*alstor(nj,k)*alstor(nj,l))
     .           +dsqrt(xmub)*xmchar(k)/gmsb(ni)*(x1-2.D0*xmut)*
     .           (alsbot(ni,k)*alsbot(ni,l)*alstor(nj,k)*akstor(nj,l)+
     .            aksbot(ni,k)*aksbot(ni,l)*alstor(nj,l)*akstor(nj,k))
     .           +dsqrt(xmub)*xmchar(l)/gmsb(ni)*(x1-2.D0*xmut)*
     .           (alsbot(ni,k)*alsbot(ni,l)*alstor(nj,l)*akstor(nj,k)+
     .            aksbot(ni,k)*aksbot(ni,l)*alstor(nj,k)*akstor(nj,l))
     .           +dsqrt(xmut*xmub)*(-2.D0*xmut+2.D0*x1-2.D0)*
     .           (alsbot(ni,k)*aksbot(ni,l)*alstor(nj,l)*akstor(nj,k)+
     .            alsbot(ni,l)*aksbot(ni,k)*alstor(nj,k)*akstor(nj,l))
     .           +dsqrt(xmut)*xmchar(l)/gmsb(ni)
     .           *(x1-xmut+xmust(nj)-xmub-1.D0)*
     .           (alsbot(ni,k)*aksbot(ni,l)*akstor(nj,k)*akstor(nj,l)+
     .            alsbot(ni,l)*aksbot(ni,k)*alstor(nj,l)*alstor(nj,k))
     .           +dsqrt(xmut)*xmchar(k)/gmsb(ni)
     .          *(x1-xmut+xmust(nj)-xmub-1.D0)*
     .           (alsbot(ni,k)*aksbot(ni,l)*alstor(nj,k)*alstor(nj,l)+
     .            alsbot(ni,l)*aksbot(ni,k)*akstor(nj,k)*akstor(nj,l))
     .           +dsqrt(xmut*xmub)*xmchar(k)*xmchar(l)/gmsb(ni)**2
     .          *(-2.D0)*
     .           (alsbot(ni,k)*aksbot(ni,l)*alstor(nj,k)*akstor(nj,l)+
     .            alsbot(ni,l)*aksbot(ni,k)*alstor(nj,l)*akstor(nj,k)) )
            endif
         end do
      end do
     
c -------------------------------------------------------------------- c
c                          neutralino exchange
c -------------------------------------------------------------------- c
      neutneut=0.D0

      do k=1,5,1
         do l=1,5,1
       if ((amneut(k)+mb).gt.gmsb(ni).and.(amneut(l)+mb).gt.gmsb(ni))
     .Then
            neutneut=neutneut+3.D0*g2s**2/dneut(k)/dneut(l)*(
     .           xmneut(k)*xmneut(l)/gmsb(ni)**2*
     .           (-1.D0+x1+x2-xmut-xmub+xmust(nj))*
     .           (abot(ni,k)*abot(ni,l)*atopr(nj,k)*atopr(nj,l)+
     .            bbot(ni,k)*bbot(ni,l)*btopr(nj,k)*btopr(nj,l))
     .           +(xmub*(x2+xmust(nj)-xmut-x1-xmub)-x2+x1*x2-x1+1.D0
     .             -xmust(nj)+xmut)*
     .           (abot(ni,k)*abot(ni,l)*btopr(nj,k)*btopr(nj,l)+
     .            bbot(ni,k)*bbot(ni,l)*atopr(nj,k)*atopr(nj,l))
     .           +xmneut(k)/gmsb(ni)*dsqrt(xmut)*(x2-2.D0*xmub)*
     .           (abot(ni,k)*abot(ni,l)*atopr(nj,k)*btopr(nj,l)+
     .            bbot(ni,k)*bbot(ni,l)*atopr(nj,l)*btopr(nj,k))
     .           +xmneut(l)/gmsb(ni)*dsqrt(xmut)*(x2-2.D0*xmub)*
     .           (abot(ni,k)*abot(ni,l)*atopr(nj,l)*btopr(nj,k)+
     .            bbot(ni,k)*bbot(ni,l)*atopr(nj,k)*btopr(nj,l))
     .           +dsqrt(xmut*xmub)*(-2.D0*xmub+2.D0*x2-2.D0)*
     .           (abot(ni,k)*bbot(ni,l)*atopr(nj,l)*btopr(nj,k)+
     .            abot(ni,l)*bbot(ni,k)*atopr(nj,k)*btopr(nj,l))
     .           +dsqrt(xmub)*xmneut(l)/gmsb(ni)*
     .           (x2-xmut+xmust(nj)-xmub-1.D0)*
     .           (abot(ni,k)*bbot(ni,l)*btopr(nj,k)*btopr(nj,l)+
     .            abot(ni,l)*bbot(ni,k)*atopr(nj,l)*atopr(nj,k))
     .           +dsqrt(xmub)*xmneut(k)/gmsb(ni)*
     .           (x2-xmut+xmust(nj)-xmub-1.D0)*
     .           (abot(ni,k)*bbot(ni,l)*atopr(nj,k)*atopr(nj,l)+
     .            abot(ni,l)*bbot(ni,k)*btopr(nj,k)*btopr(nj,l))
     .           +dsqrt(xmut*xmub)*xmneut(k)*xmneut(l)/gmsb(ni)**2*
     .           (-2.D0)*
     .           (abot(ni,k)*bbot(ni,l)*atopr(nj,k)*btopr(nj,l)+
     .            abot(ni,l)*bbot(ni,k)*atopr(nj,l)*btopr(nj,k)) )
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c                            gluino exchange
c -------------------------------------------------------------------- c
      gluiglui=0.d0
      if ((mgluino+mb).gt.gmsb(ni))Then
      gluiglui= 2.D0/3.D0*gs2**2*4.D0/dgl**2*(
     .          xmugl*dsqrt(xmut*xmub)*(-4.D0)*
     .          gbr(ni)*gbl(ni)*gtr(nj)*gtl(nj)
     .         +dsqrt(xmugl*xmut)*2.D0*(x2-2.D0*xmub)*
     .          gtr(nj)*gtl(nj)*(gbr(ni)**2+gbl(ni)**2)
     .         +dsqrt(xmut*xmub)*4.D0*(x2-xmub-1.D0)*
     .          gtr(nj)*gtl(nj)*gbr(ni)*gbl(ni)
     .         +dsqrt(xmub*xmugl)*(-2.D0)*(1.D0-x2+xmut+xmub-xmust(nj))*
     .          gbr(ni)*gbl(ni)*(gtr(nj)**2+gtl(nj)**2)
     .         +xmugl*(x1+x2-xmub-xmut-1.D0+xmust(nj))*
     .          (gtr(nj)**2*gbr(ni)**2+gtl(nj)**2*gbl(ni)**2)
     .         +(x2*xmub-xmub*xmut+x1*x2-xmub**2-x1*xmub+xmub*xmust(nj)
     .           -x1-xmust(nj)+1.D0+xmut-x2)*
     .          (gtr(nj)**2*gbl(ni)**2+gtl(nj)**2*gbr(ni)**2) )
      endif
c -------------------------------------------------------------------- c
c                    chargino neutralino interference
c -------------------------------------------------------------------- c
      charneut=0.D0

      do k=1,2,1
         do l=1,5,1
      if ((amchar(k)+amt).gt.gmsb(ni).and.(amneut(l)+mb).gt.gmsb(ni)
     .)Then
            charneut=charneut+2.D0*3.D0*g2s**2/dchi(k)/dneut(l)*(
     .           xmchar(k)*xmneut(l)/gmsb(ni)**2*
     .           (-1.D0+x1+x2-xmut-xmub+xmust(nj))*
     .           (abot(ni,l)*atopr(nj,l)*alsbot(ni,k)*alstor(nj,k)+
     .            bbot(ni,l)*btopr(nj,l)*aksbot(ni,k)*akstor(nj,k))
     .           +(x2-x1*x2+x1-1.D0-xmub-xmut-2.D0*xmut*xmub+x1*xmub
     .            +x2*xmut+xmust(nj))*
     .           (abot(ni,l)*btopr(nj,l)*alstor(nj,k)*aksbot(ni,k)+
     .            bbot(ni,l)*atopr(nj,l)*alsbot(ni,k)*akstor(nj,k))
     .           +xmchar(k)/gmsb(ni)*dsqrt(xmut)*(x2-2.D0*xmub)*
     .           (abot(ni,l)*btopr(nj,l)*alsbot(ni,k)*alstor(nj,k)+
     .            bbot(ni,l)*atopr(nj,l)*aksbot(ni,k)*akstor(nj,k))
     .           +xmneut(l)/gmsb(ni)*dsqrt(xmub)*(x1-2.D0*xmut)*
     .           (abot(ni,l)*atopr(nj,l)*alsbot(ni,k)*akstor(nj,k)+
     .            bbot(ni,l)*btopr(nj,l)*aksbot(ni,k)*alstor(nj,k))
     .           +dsqrt(xmut*xmub)*(1.D0-xmut-xmub+xmust(nj))*
     .           (abot(ni,l)*btopr(nj,l)*alsbot(ni,k)*akstor(nj,k)+
     .            bbot(ni,l)*atopr(nj,l)*aksbot(ni,k)*alstor(nj,k))
     .           +dsqrt(xmub)*xmchar(k)/gmsb(ni)*
     .           (x2-xmut+xmust(nj)-xmub-1.D0)*
     .           (abot(ni,l)*btopr(nj,l)*aksbot(ni,k)*akstor(nj,k)+
     .            bbot(ni,l)*atopr(nj,l)*alsbot(ni,k)*alstor(nj,k))
     .           +dsqrt(xmut)*xmneut(l)/gmsb(ni)*
     .           (x1-xmut+xmust(nj)-xmub-1.D0)*
     .           (abot(ni,l)*atopr(nj,l)*aksbot(ni,k)*alstor(nj,k)+
     .            bbot(ni,l)*btopr(nj,l)*alsbot(ni,k)*akstor(nj,k))
     .           +dsqrt(xmut*xmub)*xmchar(k)*xmneut(l)/gmsb(ni)**2*
     .           (-2.D0)*
     .           (abot(ni,l)*atopr(nj,l)*aksbot(ni,k)*akstor(nj,k)+
     .            bbot(ni,l)*btopr(nj,l)*alsbot(ni,k)*alstor(nj,k)) )
            endif
         end do
      end do
c -------------------------------------------------------------------- c

      NS_sbtststarb=charchar+neutneut+gluiglui+charneut

      end
c ==================================================================== c
c ======================== stop_1/2 tbar bottom ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sbtbstb(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,l,k,i
      DOUBLE PRECISION gmst(2),gmsb(2),xmust(2),dw(2),dch(2),gctbr(2,2)
      DOUBLE PRECISION
     .gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),gdr(2),gdl(2)
      DOUBLE PRECISION dneut(5),xmuneut(5)
      DOUBLE PRECISION atopr(2,5),btopr(2,5),abot(2,5),bbot(2,5)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),amch    
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      DOUBLE PRECISION chtbrunr,chtbrunl,xmut,xmub,xmuw,
     .xmuch,xmugl,x2,x1,x3,dgl,stsbotww,stsbothh,stsbotgl,stsbotneut,
     .stsbothw,stsbotwneut,stsbothcneut
*      
      COMMON/NS_indices/ni,nj
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,amch
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar     
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_coup20/gwtb,gwntau
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
*
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_neutsbotbot/abot,bbot
*
      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2
      xmust(1)   = gmst(1)**2/gmsb(ni)**2
      xmust(2)   = gmst(2)**2/gmsb(ni)**2
      xmut       = amt**2/gmsb(ni)**2
      xmub       = amb**2/gmsb(ni)**2
      xmuw       = mw**2/gmsb(ni)**2
      xmuch      = amch**2/gmsb(ni)**2
      xmugl      = mgluino**2/gmsb(ni)**2
      do i=1,5,1
         xmuneut(i) = amneut(i)**2/gmsb(ni)**2
      end do

      x3 = 2.D0-x1-x2

      dw(1)   = 1.D0-x3+xmust(1)-xmuw
      dw(2)   = 1.D0-x3+xmust(2)-xmuw
      dch(1)  = 1.D0-x3+xmust(1)-xmuch
      dch(2)  = 1.D0-x3+xmust(2)-xmuch
      dgl     = 1.D0-x2+xmub-xmugl
      do i=1,5,1
         dneut(i) = 1.D0-x2+xmub-xmuneut(i)
      end do
c -------------------------------------------------------------------- c
c                              W exchange
c -------------------------------------------------------------------- c
      stsbotww=0.d0
      if ((amw+gmst(nj)).gt.gmsb(ni))Then 
      stsbotww=3.D0*g2s**2/4.D0/dw(nj)**2*gwtb(nj,ni)**2*( 
     .     4.D0*(1.D0+x1*x2-x1-x2-xmust(nj)) 
     .   + xmub*(3.D0-3.D0*x1+x2+xmust(nj)-xmub+2.D0*xmut) 
     .   + xmut*(3.D0-3.D0*x2+x1+xmust(nj)-xmut) 
     .   + xmub/xmuw**2*(1.D0-xmust(nj))**2*
     .     (-1.D0+x1+x2+xmust(nj)-xmub+2.D0*xmut)
     .   + xmut/xmuw**2*(1.D0-xmust(nj))**2*
     .     (-1.D0+x1+x2+xmust(nj)-xmut) 
     .   + 2.D0*xmub/xmuw*(1.D0-xmust(nj))*
     .     (-x1+x2-xmub-1.D0+xmust(nj)+2.D0*xmut) 
     .   + 2.D0*xmut/xmuw*(1.D0-xmust(nj))*
     .     (-x2+x1-xmut-1.D0+xmust(nj)) ) 
      endif
c -------------------------------------------------------------------- c
c                             H+ exchange
c -------------------------------------------------------------------- c
      stsbothh=0.d0
      if ((amch+gmst(nj)).gt.gmsb(ni))Then
      stsbothh=3.D0*g2s**2*gctbr(nj,ni)**2*amw**2/gmsb(ni)**2
     .     /dch(nj)**2*(
     .     (chtbrunr**2+chtbrunl**2)*(x1+x2+xmust(nj)-1.D0-xmub-xmut)
     .     -4.D0*chtbrunr*chtbrunl*dsqrt(xmut*xmub) )
      endif
c -------------------------------------------------------------------- c
c                           gluino exchange
c -------------------------------------------------------------------- c
      stsbotgl=0.d0
      if ((mgluino+amb).gt.gmsb(ni))Then
      stsbotgl=gs2*gs2*2.D0/3.D0/dgl**2*4.D0*(
     .     gbr(ni)*gbl(ni)*gtr(nj)*gtl(nj)*dsqrt(xmut*xmub)*(
     .     -4.D0*xmugl-4.D0*(1.D0+xmub-x2) ) +
     .     gbr(ni)*gbl(ni)*(gtr(nj)**2+gtl(nj)**2)*dsqrt(xmub*xmugl)*
     .     2.D0*(-xmub-xmut+xmust(nj)+x2-1.D0) +
     .     gtr(nj)*gtl(nj)*(gbr(ni)**2+gbl(ni)**2)*dsqrt(xmut*xmugl)*
     .     2.D0*(x2-2.D0*xmub) +
     .     (gbr(ni)**2*gtl(nj)**2+gbl(ni)**2*gtr(nj)**2)*xmugl*
     .     (-1.D0-xmub-xmut+xmust(nj)+x1+x2) +
     .     (gbr(ni)**2*gtr(nj)**2+gbl(ni)**2*gtl(nj)**2)*
     .     (1.D0+x1*x2-x1-x2-xmust(nj)+xmut+xmub*(-xmut+xmust(nj)+x2
     .      -x1-xmub)) )
      endif
c -------------------------------------------------------------------- c
c                        neutralino exchange
c -------------------------------------------------------------------- c
      stsbotneut = 0.D0

      do k=1,5,1
         do l=1,5,1
       If
     .((amneut(k)+amb).gt.gmsb(ni).and.(amneut(l)+amb).gt.gmsb(ni))Then
            stsbotneut=stsbotneut+3.D0*g2s**2/dneut(k)/dneut(l)*(
     .          (atopr(nj,k)*atopr(nj,l)*abot(ni,k)*abot(ni,l)+
     .           btopr(nj,k)*btopr(nj,l)*bbot(ni,k)*bbot(ni,l))*
     .          ((1.D0-x1)*(1.D0-x2)-xmust(nj)+xmub*(-x1+x2+xmust(nj)
     .           -xmub-xmut)+xmut)
     .         +(atopr(nj,k)*atopr(nj,l)*bbot(ni,k)*bbot(ni,l)+
     .           btopr(nj,k)*btopr(nj,l)*abot(ni,k)*abot(ni,l))*
     .           xmneut(k)*xmneut(l)/gmsb(ni)**2*(x1+x2-1.D0+xmust(nj)
     .           -xmub-xmut)
     .         +(atopr(nj,k)*btopr(nj,l)*abot(ni,k)*bbot(ni,l)+
     .           btopr(nj,k)*atopr(nj,l)*bbot(ni,k)*abot(ni,l))*
     .           2.D0*dsqrt(xmut*xmub)*(-1.D0+x2-xmub)
     .         +dsqrt(xmut)*(x2-2.D0*xmub)*( xmneut(k)/gmsb(ni)*
     .          (atopr(nj,k)*btopr(nj,l)*bbot(ni,k)*bbot(ni,l)+
     .           btopr(nj,k)*atopr(nj,l)*abot(ni,k)*abot(ni,l))
     .         +xmneut(l)/gmsb(ni)*
     .          (atopr(nj,k)*btopr(nj,l)*abot(ni,k)*abot(ni,l)+
     .           btopr(nj,k)*atopr(nj,l)*bbot(ni,k)*bbot(ni,l)) )
     .         +dsqrt(xmub)*(-1.D0+x2+xmust(nj)-xmut-xmub)* 
     .         (xmneut(k)/gmsb(ni)*
     .          (atopr(nj,k)*atopr(nj,l)*bbot(ni,k)*abot(ni,l)+
     .           btopr(nj,k)*btopr(nj,l)*abot(ni,k)*bbot(ni,l))
     .         +xmneut(l)/gmsb(ni)*
     .          (atopr(nj,k)*atopr(nj,l)*abot(ni,k)*bbot(ni,l)+
     .           btopr(nj,k)*btopr(nj,l)*bbot(ni,k)*abot(ni,l)) )
     .         +(atopr(nj,l)*btopr(nj,k)*abot(ni,k)*bbot(ni,l)+
     .           btopr(nj,l)*atopr(nj,k)*bbot(ni,k)*abot(ni,l))*
     .           (-2.D0)*dsqrt(xmut*xmub)*
     .             xmneut(k)*xmneut(l)/gmsb(ni)**2 )
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c                           H+ W interference
c -------------------------------------------------------------------- c
      stsbothw=0.d0
      if((amch+gmst(nj)).gt.gmsb(ni).and.(amw+gmst(nj)).gt.gmsb(ni))
     .Then
      stsbothw=-3.D0*g2s**2*gctbr(nj,ni)*(-gwtb(nj,ni))*amw/gmsb(ni)
     .     /dw(nj)/dch(nj)*
     .     (dsqrt(xmub)*chtbrunr*(xmust(nj)+xmut-xmub+x2-x1-1.D0)
     .     +dsqrt(xmut)*chtbrunl*(-xmust(nj)+xmut-xmub+x2-x1+1.D0)
     .     +dsqrt(xmub)/xmuw*chtbrunr*(xmust(nj)-1.D0)*(xmub-xmut-x1-x2
     .      +1.D0-xmust(nj))
     .     +dsqrt(xmut)/xmuw*chtbrunl*(xmust(nj)-1.D0)*(xmub-xmut+x1+x2
     .      -1.D0+xmust(nj)) )
      endif
c -------------------------------------------------------------------- c
c                       neutralino W interference
c -------------------------------------------------------------------- c
      stsbotwneut = 0.D0

      do l=1,5,1
       
      if((amneut(l)+amb).gt.gmsb(ni).and.(amw+gmst(nj)).gt.gmsb(ni))Then
         stsbotwneut=stsbotwneut
     .     +3.D0*g2s**2/dw(nj)/dneut(l)*(-gwtb(nj,ni))*(
     .      btopr(nj,l)*bbot(ni,l)*dsqrt(xmub*xmut)*(1.D0/xmuw*
     .      (xmust(nj)*(xmust(nj)+xmub-xmut-2.D0)-xmub+xmut+1.D0) +
     .      2.D0*x2-xmust(nj)-xmub+xmut-3.D0) +
     .      atopr(nj,l)*abot(ni,l)*(1.D0/xmuw*(1.D0-xmust(nj))*(xmub*
     .      (xmust(nj)+x2-xmub+xmut-1.D0)-x2*xmut) + xmub*(1.D0+x2+
     .      xmust(nj)+xmut-2.D0*x1-xmub)-2.D0*xmust(nj)+2.D0*x1*x2+2.D0
     .      -2.D0*x1-2.D0*x2+xmut*(2.D0-x2)) +
     .      btopr(nj,l)*abot(ni,l)*dsqrt(xmut)*xmneut(l)/gmsb(ni)*
     .      (1.D0/xmuw*
     .      (xmust(nj)-1.D0)*(xmust(nj)+x1+x2-1.D0-xmut+xmub)+x2-x1
     .      -xmust(nj)-xmub+xmut+1.D0) +
     .      atopr(nj,l)*bbot(ni,l)*dsqrt(xmub)*xmneut(l)/gmsb(ni)*
     .      (1.D0/xmuw*
     .      (1.D0-xmust(nj))*(xmust(nj)+x1+x2-1.D0+xmut-xmub)+x2-x1
     .      +xmust(nj)-xmub+xmut-1.D0) )
         endif

      end do
c -------------------------------------------------------------------- c
c                    neutralino Higgs interference
c -------------------------------------------------------------------- c
      stsbothcneut = 0.D0
      
      do l=1,5,1
        
      if
     .((amneut(l)+amb).gt.gmsb(ni).and.(amch+gmst(nj)).gt.gmsb(ni))Then
         stsbothcneut=stsbothcneut
     .     +3.D0*g2s**2*2.D0/dch(nj)/dneut(l)*
     .      (-gctbr(nj,ni))*amw/gmsb(ni)*(
     .      (chtbrunl*bbot(ni,l)*atopr(nj,l)+chtbrunr*abot(ni,l)*
     .       btopr(nj,l))*dsqrt(xmub*xmut)*xmneut(l)/gmsb(ni)*(-2.D0) +
     .      (chtbrunl*bbot(ni,l)*btopr(nj,l)+chtbrunr*abot(ni,l)*
     .       atopr(nj,l))*dsqrt(xmub)*(-1.D0-xmub-xmut+xmust(nj)+x2) +
     .      (chtbrunl*abot(ni,l)*btopr(nj,l)+chtbrunr*bbot(ni,l)*
     .       atopr(nj,l))*xmneut(l)/gmsb(ni)*(-1.D0-xmub-xmut+xmust(nj)
     .       +x1+x2) +
     .      (chtbrunl*abot(ni,l)*atopr(nj,l)+chtbrunr*bbot(ni,l)*
     .       btopr(nj,l))*dsqrt(xmut)*(-2.D0*xmub+x2) )
         endif

      end do

c -------------------------------------------------------------------- c

      NS_sbtbstb = stsbotww+stsbothh+stsbothw+stsbotgl+stsbotneut+
     .          stsbotwneut+stsbothcneut

      end
c ==================================================================== c
c ====================== stop_1/2 tau- nu_taubar ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sbtaustnu(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj
      DOUBLE PRECISION gmst(2),gmsb(2),xmust(2),dw(2),dch(2),gctbr(2,2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2)
      DOUBLE PRECISION 
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),amch
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION chctaunur,chctaunul
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION xmutau,xmuw,xmuch,x1,x2,x3,ststauww,ststauhh
     .,ststauhw
*      
      COMMON/NS_indices/ni,nj
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,amch
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/NS_coup14/chctaunur,chctaunul
      COMMON/NS_coup20/gwtb,gwntau
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_hcsbotstop/gctbr

      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2
*
      xmust(nj) = gmst(nj)**2/gmsb(ni)**2
      xmutau   = amtau**2/gmsb(ni)**2
      xmuw     = mw**2/gmsb(ni)**2
      xmuch    = amch**2/gmsb(ni)**2

      x3 = 2.D0-x1-x2

      dw(nj)  = 1.D0-x3+xmust(nj)-xmuw
      dch(nj) = 1.D0-x3+xmust(nj)-xmuch
 
c -------------------------------------------------------------------- c
c                             W exchange
c -------------------------------------------------------------------- c
      ststauww=0.d0

      if ((mw+gmst(nj)).gt.gmsb(ni))Then
      ststauww=1.D0/4.D0/dw(nj)**2*(-gwtb(nj,ni))**2*( 
     .     4.D0*(1.D0+x1*x2-x1-x2-xmust(nj)) + 
     .     xmutau*(3.D0-3.D0*x1+x2+xmust(nj)-xmutau) 
     .     +xmutau/xmuw**2*(1.D0-xmust(nj))**2*(-1.D0+xmust(nj)+x1+x2
     .      -xmutau)
     .     +2.D0*xmutau/xmuw*(1.D0-xmust(nj))*(-x1+x2-xmutau-1.D0+
     .      xmust(nj)) )
      endif

c -------------------------------------------------------------------- c
c                             H+ exchange
c -------------------------------------------------------------------- c
      ststauhh=0.d0
 
      if ((amch+gmst(nj)).gt.gmsb(ni))Then
      ststauhh=gctbr(nj,ni)**2*amw**2/gmsb(ni)**2/dch(nj)**2*
     .     chctaunur**2*(x1+x2+xmust(nj)-1.D0-xmutau)
      endif

c -------------------------------------------------------------------- c
c                          H+ W interference
c -------------------------------------------------------------------- c
      ststauhw=0.d0

      if(((amch+gmst(nj)).gt.gmsb(ni)).and.(mw+gmst(nj)).gt.gmsb(ni))
     .Then
      ststauhw=-gctbr(nj,ni)*(-gwtb(nj,ni))*amw/gmsb(ni)/dw(nj)/dch(nj)*
     .      dsqrt(xmutau)*chctaunur*(
     .      -1.D0-x1+x2-xmutau+xmust(nj)+
     .      1.D0/xmuw*(xmust(nj)-1.D0)*(xmutau-x1-x2+1.D0-xmust(nj))) 
      endif

c -------------------------------------------------------------------- c
      NS_sbtaustnu = ststauww+ststauhh+ststauhw

      end
c ==================================================================== c
c ======================== stop_1/2 e- nu_ebar ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sbelstnu(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj
      DOUBLE PRECISION gmst(2),gmsb(2),xmust(2),dw(2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION XMUW,X1,X2,X3,stselww
*      
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup20/gwtb,gwntau
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
*
      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2

      xmust(nj) = gmst(nj)**2/gmsb(ni)**2
      xmuw     = mw**2/gmsb(ni)**2
      x3 = 2.D0-x1-x2
      dw(nj) = 1.D0-x3+xmust(nj)-xmuw
      
c -------------------------------------------------------------------- c
c                             W exchange
c -------------------------------------------------------------------- c
      stselww=0.d0
      If((mw+gmst(nj)).gt.gmsb(ni))Then
      stselww=1.D0/4.D0/dw(nj)**2*(-gwtb(nj,ni))**2*
     .        4.D0*(1.D0+x1*x2-x1-x2-xmust(nj))

      endif

c -------------------------------------------------------------------- c
      NS_sbelstnu = stselww
     
      end
c ==================================================================== c
c ======================= sbottom1* bottom bottom ==================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sb2sb1starbb(x1,x2)
*
      IMPLICIT NONE
      INTEGER l,k,i
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION dneut(5),xmuneut(5)
      DOUBLE PRECISION
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      DOUBLE PRECISION xmub,xmugl,xmusb1,x1,x2,dgl,sb2sb1neut
     .,sb2sb1gg
*
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
c --- several definitions ---
     
      xmub   = amb**2/asb2**2
      do i=1,5,1
         xmuneut(i) = amneut(i)**2/asb2**2
      end do
      xmugl  = mgluino**2/asb2**2
      xmusb1 = asb1**2/asb2**2

      do i=1,5,1
         dneut(i) = 1.D0-x1+xmub-xmuneut(i)
      end do
      dgl = 1.D0-x1+xmub-xmugl

      NS_sb2sb1starbb=0.d0
c -------------------------------------------------------------------- c
c                          neutralino exchange
c -------------------------------------------------------------------- c
      sb2sb1neut=0.D0
      do k=1,5,1
         do l=1,5,1
      if((amneut(k)+mb).gt.asb2.and.(amneut(l)+mb).gt.asb2)Then
            sb2sb1neut=sb2sb1neut+1.D0/dneut(k)/dneut(l)*(
     .          (bbot(1,k)*bbot(1,l)*abot(2,k)*abot(2,l)+
     .           abot(1,k)*abot(1,l)*bbot(2,k)*bbot(2,l))*
     .          ((1.D0-x1)*(1.D0-x2)-xmusb1+xmub*(x1-x2+xmusb1
     .           -2.D0*xmub)+xmub)
     .         +(bbot(1,k)*bbot(1,l)*bbot(2,k)*bbot(2,l)+
     .           abot(1,k)*abot(1,l)*abot(2,k)*abot(2,l))*
     .           xmneut(k)*xmneut(l)/asb2**2*(x1+x2-1.D0+xmusb1
     .           -2.D0*xmub)
     .         +(bbot(1,k)*abot(1,l)*abot(2,k)*bbot(2,l)+
     .           abot(1,k)*bbot(1,l)*bbot(2,k)*abot(2,l))*
     .           2.D0*xmub*(-1.D0+x1-xmub)
     .         +dsqrt(xmub)*(x1-2.D0*xmub)*(xmneut(k)/asb2*
     .          (bbot(1,k)*abot(1,l)*bbot(2,k)*bbot(2,l)+
     .           abot(1,k)*bbot(1,l)*abot(2,k)*abot(2,l))
     .         +xmneut(l)/asb2*
     .          (bbot(1,k)*abot(1,l)*abot(2,k)*abot(2,l)+
     .           abot(1,k)*bbot(1,l)*bbot(2,k)*bbot(2,l)) )
     .         +dsqrt(xmub)*(-1.D0+x1+xmusb1-2.D0*xmub)* 
     .         (xmneut(k)/asb2*
     .          (bbot(1,k)*bbot(1,l)*bbot(2,k)*abot(2,l)+
     .           abot(1,k)*abot(1,l)*abot(2,k)*bbot(2,l))
     .         +xmneut(l)/asb2*
     .          (bbot(1,k)*bbot(1,l)*abot(2,k)*bbot(2,l)+
     .           abot(1,k)*abot(1,l)*bbot(2,k)*abot(2,l)) )
     .         +(bbot(1,l)*abot(1,k)*abot(2,k)*bbot(2,l)+
     .           abot(1,l)*bbot(1,k)*bbot(2,k)*abot(2,l))*
     .           (-2.D0)*xmub*xmneut(k)*xmneut(l)/asb2**2 )
            endif

         end do
      end do
c -------------------------------------------------------------------- c
c                             gluino exchange
c -------------------------------------------------------------------- c
      sb2sb1gg=0.d0
      if (xmugl.gt.1.d0)Then
        sb2sb1gg=1.D0/dgl**2*4.D0*( -4.D0*dsqrt(xmub*xmugl)*xmub*
     .     (gbl(2)*gbl(1)+gbr(1)*gbr(2))*(gbr(1)*gbl(2)+gbl(1)*gbr(2))
     .    +2.D0*dsqrt(xmub*xmugl)*(
     .     (gbr(1)*gbl(1)*(gbl(2)**2+gbr(2)**2)+gbl(2)*gbr(2)*
     .     (gbr(1)**2+gbl(1)**2))*x1 +
     .     gbl(2)*gbr(2)*(gbr(1)**2+gbl(1)**2)*(xmusb1-1.D0))
     .    +(-2.D0)*xmub*xmugl*(gbr(1)*gbr(2)+gbl(2)*gbl(1))**2
     .    +xmub*((gbr(1)**2*gbl(2)**2+gbl(1)**2*gbr(2)**2)*(1.D0+x1
     .     -x2+xmusb1)+4.D0*gbr(1)*gbl(2)*gbl(1)*gbr(2)*(x1-1.D0))
     .    +xmub**2*(-2.D0)*(gbr(1)*gbl(2)+gbl(1)*gbr(2))**2
     .    +xmugl*(gbr(1)**2*gbr(2)**2+gbl(2)**2*gbl(1)**2)*(x1+x2+
     .     xmusb1-1.D0)
     .    +(gbr(1)**2*gbl(2)**2+gbl(1)**2*gbr(2)**2)*(x1*x2-x1-x2
     .     -xmusb1+1.D0) )

        endif

              NS_sb2sb1starbb=3.D0*g2s**2*sb2sb1neut+
     .             2.D0/3.D0*gs2**2*sb2sb1gg

c      write(*,*)'sb2sb1starbb',NS_sb2sb1starbb,sb2sb1neut,
c     .sb2sb1gg,mgluino
c     .,mgluino,asb2,asb1,gbl(1),gbl(2),gbr(1),gbr(2)


      end
c ==================================================================== c
c ===================== sbottom1 bottom bottombar ==================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sb2sb1bb(x1,x2)
*
      IMPLICIT NONE
      INTEGER l,k,i,j
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION dneut(5),xmuneut(5)
      DOUBLE PRECISION gztautau(2,2),gztt(2,2),gzbb(2,2)
      DOUBLE PRECISION
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot      
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION Hbbr(3),Abbr(2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION P(2,3)
      DOUBLE PRECISION xmuh(3),xmua(2),dh(3),da(2)
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION azbb11,azbb12,azbb21,azbb22,xmub,xmuz,
     .xmugl,xmusb1,x3,x1,x2,dgl,dz,vzz,azz,sb2sb1zz,sb2sb1hk,
     .sb2sb1aa,sb2sb1higgs,sb2sb1neut,sb2sb1gg,sb2sb1neutz,
     .sb2sb1hneut,sb2sb1hz
*     
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_CPodd_MIX/P
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr 
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_phibotbot/Hbbr,Abbr

      azbb11=gzbb(1,1)
      azbb12=gzbb(1,2)
      azbb21=gzbb(2,1)
      azbb22=gzbb(2,2)
c --- several definitions ---

      xmub   = amb**2/asb2**2
      do i=1,5,1
         xmuneut(i) = amneut(i)**2/asb2**2
      end do
      xmugl  = mgluino**2/asb2**2
      xmuz   = mz**2/asb2**2
      xmusb1 = asb1**2/asb2**2
      DO i =1,3
      xmuh(i) = SMASS(i)**2/asb2**2
      enddo
      DO i =1,2
      xmua(i) = PMASS(i)**2/asb2**2
      enddo
      x3  = 2.D0-x1-x2

      do i=1,5,1
         dneut(i) = 1.D0-x1+xmub-xmuneut(i)
      end do
      dgl = 1.D0-x1+xmub-xmugl
      dz  = 1.D0-x3+xmusb1-xmuz

      DO i =1,3
      dh(i) = 1.D0-x3+xmusb1-xmuh(i)
      enddo
      DO i =1,2
      da(i) = 1.D0-x3+xmusb1-xmua(i)
      enddo
      vzz = vzzbotbot
      azz = azzbotbot
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      sb2sb1zz=0.d0
      
      If((mz+asb1).gt.asb2)Then
      sb2sb1zz = 1.D0/4.D0/cw**2/dz**2*azbb12**2*(
     .      8.D0*xmub*azz**2*(-x1-x2+xmusb1+3.D0)
     .     -8.D0*(azz**2+vzz**2)*(x1+x2-1.D0+xmusb1-x1*x2)
     .     -16.D0*azz**2*(1.D0-xmusb1)**2*xmub/xmuz
     .     +8.D0*azz**2*xmub/xmuz**2*(1.D0-xmusb1)**2*
     .     (x1+x2+xmusb1-1.D0))
      endif

c -------------------------------------------------------------------- c
c                              Higgs exchange
c -------------------------------------------------------------------- c
      sb2sb1hk=0.0d0
      sb2sb1aa=0.0d0

             DO I = 1,3
                DO J = 1,3
      if((SMASS(i)+asb1).gt.asb2.and.(SMASS(j)+asb1).gt.asb2)Then
           sb2sb1hk = sb2sb1hk+(2.D0/asb2**2/dh(I)/dh(J)*
     .(Hbbr(I)/dsqrt(2.D0))*(Hbbr(J)/dsqrt(2.D0))*
     .    amz**4/amw**2*Hsbotsbotr(I,1,2)*Hsbotsbotr(J,1,2))*
     .    (-1.D0+x1+x2+xmusb1-4.D0*xmub)
      endif
                ENDDO
             ENDDO
              DO I = 1,2
                DO J = 1,2
      if((PMASS(i)+asb1).gt.asb2.and.(PMASS(j)+asb1).gt.asb2)Then
            sb2sb1aa = sb2sb1aa+2.D0/asb2**2/da(I)/da(J)*
     .(Abbr(I)/dsqrt(2.D0))*(Abbr(J)/dsqrt(2.D0))*
     .     amz**4/amw**2*Asbotsbotr(I,1,2)*Asbotsbotr(J,1,2)
     .*(-1.D0+x1+x2+xmusb1)
      endif
                ENDDO
             ENDDO   
      sb2sb1higgs=sb2sb1hk+sb2sb1aa
c -------------------------------------------------------------------- c
c                          neutralino exchange
c -------------------------------------------------------------------- c
      sb2sb1neut=0.D0

      do k=1,5,1
         do l=1,5,1
      if((amneut(k)+mb).gt.asb2.and.(amneut(l)+mb).gt.asb2)Then
            sb2sb1neut=sb2sb1neut+1.D0/dneut(k)/dneut(l)*(
     .          (abot(1,k)*abot(1,l)*abot(2,k)*abot(2,l)+
     .           bbot(1,k)*bbot(1,l)*bbot(2,k)*bbot(2,l))*
     .          ((1.D0-x1)*(1.D0-x2)-xmusb1+xmub*(x1-x2+xmusb1
     .           -2.D0*xmub)+xmub)
     .         +(abot(1,k)*abot(1,l)*bbot(2,k)*bbot(2,l)+
     .           bbot(1,k)*bbot(1,l)*abot(2,k)*abot(2,l))*
     .           xmneut(k)*xmneut(l)/asb2**2*(x1+x2-1.D0+xmusb1
     .           -2.D0*xmub)
     .         +(abot(1,k)*bbot(1,l)*abot(2,k)*bbot(2,l)+
     .           bbot(1,k)*abot(1,l)*bbot(2,k)*abot(2,l))*
     .           2.D0*xmub*(-1.D0+x1-xmub)
     .         +dsqrt(xmub)*(x1-2.D0*xmub)*(xmneut(k)/asb2*
     .          (abot(1,k)*bbot(1,l)*bbot(2,k)*bbot(2,l)+
     .           bbot(1,k)*abot(1,l)*abot(2,k)*abot(2,l))
     .         +xmneut(l)/asb2*
     .          (abot(1,k)*bbot(1,l)*abot(2,k)*abot(2,l)+
     .           bbot(1,k)*abot(1,l)*bbot(2,k)*bbot(2,l)) )
     .         +dsqrt(xmub)*(-1.D0+x1+xmusb1-2.D0*xmub)* 
     .         (xmneut(k)/asb2*
     .          (abot(1,k)*abot(1,l)*bbot(2,k)*abot(2,l)+
     .           bbot(1,k)*bbot(1,l)*abot(2,k)*bbot(2,l))
     .         +xmneut(l)/asb2*
     .          (abot(1,k)*abot(1,l)*abot(2,k)*bbot(2,l)+
     .           bbot(1,k)*bbot(1,l)*bbot(2,k)*abot(2,l)) )
     .         +(abot(1,l)*bbot(1,k)*abot(2,k)*bbot(2,l)+
     .           bbot(1,l)*abot(1,k)*bbot(2,k)*abot(2,l))*
     .           (-2.D0)*xmub*xmneut(k)*xmneut(l)/asb2**2 )
            endif
         end do
      end do
c -------------------------------------------------------------------- c
c                             gluino exchange
c -------------------------------------------------------------------- c
      sb2sb1gg=0.d0
      if(xmugl.gt.1.0d0)Then
      sb2sb1gg=1.D0/dgl**2*4.D0*( -4.D0*dsqrt(xmub*xmugl)*xmub*
     .     (gbr(2)*gbl(1)+gbr(1)*gbl(2))*(gbr(1)*gbr(2)+gbl(1)*gbl(2))
     .    +2.D0*dsqrt(xmub*xmugl)*(
     .     (gbr(1)*gbl(1)*(gbr(2)**2+gbl(2)**2)+gbr(2)*gbl(2)*
     .     (gbr(1)**2+gbl(1)**2))*x1 +
     .     gbr(2)*gbl(2)*(gbr(1)**2+gbl(1)**2)*(xmusb1-1.D0))
     .    +(-2.D0)*xmub*xmugl*(gbr(1)*gbl(2)+gbr(2)*gbl(1))**2
     .    +xmub*((gbr(1)**2*gbr(2)**2+gbl(1)**2*gbl(2)**2)*(1.D0+x1
     .     -x2+xmusb1)+4.D0*gbr(1)*gbr(2)*gbl(1)*gbl(2)*(x1-1.D0))
     .    +xmub**2*(-2.D0)*(gbr(1)*gbr(2)+gbl(1)*gbl(2))**2
     .    +xmugl*(gbr(1)**2*gbl(2)**2+gbr(2)**2*gbl(1)**2)*(x1+x2+
     .     xmusb1-1.D0)
     .    +(gbr(1)**2*gbr(2)**2+gbl(1)**2*gbl(2)**2)*(x1*x2-x1-x2
     .     -xmusb1+1.D0) )
      endif
      
c -------------------------------------------------------------------- c
c                         neutralino Z interference
c -------------------------------------------------------------------- c
      sb2sb1neutz=0.D0

      do l=1,5,1
         if ((amneut(l)+mb).gt.asb2.and.(mz+asb1).gt.asb2)Then
         sb2sb1neutz=sb2sb1neutz+azbb12/cw/dz/dneut(l)*(
     .      xmub*(abot(1,l)*abot(2,l)*(vzz-azz)+
     .            bbot(1,l)*bbot(2,l)*(vzz+azz))*(
     .      1.D0/xmuz*(xmusb1*(xmusb1-2.D0)+1.D0)+
     .      2.D0*x1-xmusb1-3.D0)
     .     +(abot(1,l)*abot(2,l)*(vzz+azz)+
     .       bbot(1,l)*bbot(2,l)*(vzz-azz))*(
     .      1.D0/xmuz*(1.D0-xmusb1)*(xmub*(xmusb1+x1-1.D0)-x1*xmub)
     .      +xmub*(1.D0+x1+xmusb1-2.D0*x2)-2.D0*xmusb1+2.D0*x1*x2+2.D0
     .      -2.D0*(x1+x2)+xmub*(-x1+2.D0) )
     .     +dsqrt(xmub)*xmneut(l)/asb2*
     .      (abot(1,l)*bbot(2,l)*(vzz-azz)+
     .       bbot(1,l)*abot(2,l)*(vzz+azz))*(
     .      1.D0/xmuz*(xmusb1-1.D0)*(xmusb1+x1+x2-1.D0)+x1-x2-xmusb1
     .      +1.D0 )
     .     +dsqrt(xmub)*xmneut(l)/asb2*
     .      (abot(1,l)*bbot(2,l)*(vzz+azz)+
     .       bbot(1,l)*abot(2,l)*(vzz-azz))*(
     .      1.D0/xmuz*(1.D0-xmusb1)*(xmusb1+x1+x2-1.D0)+x1-x2+xmusb1
     .      -1.D0 ) )
          endif
      end do
c -------------------------------------------------------------------- c
c                       neutralino Higgs interference
c -------------------------------------------------------------------- c
      sb2sb1hneut=0.D0

      do l=1,5,1
           do I = 1,3
       if ((amneut(l)+mb).gt.asb2.and.(smass(I)+asb1).gt.asb2)Then
         sb2sb1hneut=sb2sb1hneut-2.D0*(Hbbr(I)/dsqrt(2.D0)/dneut(l)
     .        /dh(I)/asb2*amz**2/amw*Hsbotsbotr(I,1,2))
     .        *((abot(1,l)*bbot(2,l)+abot(2,l)*bbot(1,l))*
     .        xmneut(l)/asb2*(x1+x2+xmusb1-1.D0-4.D0*xmub) 
     .        +(abot(1,l)*abot(2,l)+bbot(2,l)*bbot(1,l))*
     .        dsqrt(xmub)*(2.D0*x1-1.D0+xmusb1-4.D0*xmub) )
       endif
          enddo
          ENDDO
      do l=1,5,1
           do I = 1,2
        if ((amneut(l)+mb).gt.asb2.and.(pmass(I)+asb1).gt.asb2)Then
              sb2sb1hneut=sb2sb1hneut+
     .        2.D0*Abbr(I)/dsqrt(2.D0)/dneut(l)/da(I)/asb2*amz**2/amw*
     .        (-Asbotsbotr(I,1,2))*(
     .        (abot(1,l)*abot(2,l)-bbot(2,l)*bbot(1,l))*
     .        dsqrt(xmub)*(1.D0-xmusb1) +
     .        (abot(1,l)*bbot(2,l)-abot(2,l)*bbot(1,l))*
     .        xmneut(l)/asb2*(1.D0-x1-x2-xmusb1) )
        endif
              enddo
          ENDDO
          
c -------------------------------------------------------------------- c
c                          Higgs Z interference
c -------------------------------------------------------------------- c
          sb2sb1hz=0.0d0

          do I =1,3
       if ((smass(I)+asb1).gt.asb2.and.(mz+asb1).gt.asb2)Then
      sb2sb1hz=sb2sb1hz-2.D0/2.D0/cw*
     .   (azbb12*Hbbr(I)/dsqrt(2.D0)*amz**2/amw*Hsbotsbotr(I,1,2)
     .   /asb2/dz/dh(I))*2.D0*dsqrt(xmub)*vzz*2.D0*(x1-x2)
       endif
          enddo
          
          do I =1,2
       if ((pmass(I)+asb1).gt.asb2.and.(mz+asb1).gt.asb2)Then
      sb2sb1hz=sb2sb1hz+2.D0/2.D0/cw*
     .    azbb12*Abbr(I)/dsqrt(2.D0)*amz**2/amw*
     .(-Asbotsbotr(I,1,2))/asb2/dz/da(I)*
     .   (2.D0*dsqrt(xmub)*azz*(2.D0/xmuz*(1.D0+(xmusb1-1.D0)*(x1
     .   +x2)+xmusb1**2-2.D0*xmusb1)+2.D0-2.D0*xmusb1) )
       endif
         ENDDO
         
c -------------------------------------------------------------------- c

      NS_sb2sb1bb=3.D0*g2s**2*(sb2sb1zz+sb2sb1higgs+sb2sb1neut+
     .         sb2sb1neutz+sb2sb1hneut+sb2sb1hz) +
     .         gs2**2*2.D0/3.D0*sb2sb1gg
c      write(*,*)'ji',sb2sb1zz,sb2sb1higgs,sb2sb1neut,sb2sb1gg
c     .         sb2sb1neutz,sb2sb1hneut,sb2sb1hz,sb2sb1gg,asb1,asb2

      end
c ==================================================================== c
c ========================= sbottom1 top topbar ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sb2sb1tt(x1,x2)
*
      IMPLICIT NONE
      INTEGER i,j,k,l
      DOUBLE PRECISION dchi(2),xmuchar(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alsbot1(2,2),
     .          blsbot1(2,2),alsbot2(2,2),blsbot2(2,2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      DOUBLE PRECISION Httr(3),Attr(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION P(2,3)
      DOUBLE PRECISION xmuh(3),xmua(2),dh(3),da(2)
      DOUBLE PRECISION azbb11,azbb12,azbb21,azbb22,xmut,xmuz,
     .xmusb1,x3,x1,x2,dz,vzz,azz,sb2sb1zz,sb2sb1hk,
     .sb2sb1aa,sb2sb1higgs,sb2sb1chi,sb2sb1chiz,
     .sb2sb1hz,sb2sb1hchi
*      
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr 
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar

      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_CPodd_MIX/P
      COMMON/NS_MZMWscaleQ/AMZ,AMW
*
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_phitoptop/Httr,Attr
*
      azbb11=gzbb(1,1)
      azbb12=gzbb(1,2)
      azbb21=gzbb(2,1)
      azbb22=gzbb(2,2)
      do i=1,2,1
         do j=1,2,1
            alsbot1(i,j) = alsbot(i,j)
            blsbot1(i,j) = aksbot(i,j)
            alsbot2(i,j) = alsbot(i,j)
            blsbot2(i,j) = aksbot(i,j)
         end do
      end do
     
c --- several definitions ---

      xmut       = amt**2/asb2**2
      xmuchar(1) = amchar(1)**2/asb2**2
      xmuchar(2) = amchar(2)**2/asb2**2
      xmuz       = mz**2/asb2**2
      xmusb1     = asb1**2/asb2**2

      do i =1,3
      xmuh(i) = SMASS(i)**2/asb2**2
      enddo
      do i =1,2
      xmua(i) = PMASS(i)**2/asb2**2
      enddo
      x3  = 2.D0-x1-x2

      dchi(1) = 1.D0-x1+xmut-xmuchar(1)
      dchi(2) = 1.D0-x1+xmut-xmuchar(2)
      dz      = 1.D0-x3+xmusb1-xmuz

      do i =1,3
      dh(i) = 1.D0-x3+xmusb1-xmuh(i)
      enddo
      do i =1,2
      da(i) = 1.D0-x3+xmusb1-xmua(i)
      enddo
      vzz = vzztoptop
      azz = azztoptop
c -------------------------------------------------------------------- c
c                               Z exchange
c -------------------------------------------------------------------- c
      sb2sb1zz=0.d0
      if((mz+asb1).gt.asb2)Then
      sb2sb1zz = 1.D0/4.D0/cw**2/dz**2*azbb12**2*(
     .      8.D0*xmut*azz**2*(-x1-x2+xmusb1+3.D0)
     .     -8.D0*(azz**2+vzz**2)*(x1+x2-1.D0+xmusb1-x1*x2)
     .     -16.D0*azz**2*(1.D0-xmusb1)**2*xmut/xmuz
     .     +8.D0*azz**2*xmut/xmuz**2*(1.D0-xmusb1)**2*
     .     (x1+x2+xmusb1-1.D0))
      endif
c -------------------------------------------------------------------- c
c                              Higgs exchange
c -------------------------------------------------------------------- c
      sb2sb1hk=0.0d0
      sb2sb1aa=0.0d0

       DO I = 1,3 
               DO J = 1,3
       if((smass(I)+asb1).gt.asb2.and.(smass(J)+asb1).gt.asb2)then
           sb2sb1hk = sb2sb1hk+(2.D0/asb2**2/dh(I)/dh(J)*
     .(Httr(I)/dsqrt(2.D0))*(Httr(J)/dsqrt(2.D0))*
     .    amz**4/amw**2*Hsbotsbotr(I,1,2)*Hsbotsbotr(J,1,2))*
     .    (-1.D0+x1+x2+xmusb1-4.D0*xmut)
       endif
               ENDDO
       ENDDO
       DO I = 1,2
                DO J = 1,2
      if((pmass(I)+asb1).gt.asb2.and.(pmass(J)+asb1).gt.asb2)then  
            sb2sb1aa = sb2sb1aa+2.D0/asb2**2/da(I)/da(J)*
     .(Attr(I)/dsqrt(2.D0))*(Attr(J)/dsqrt(2.D0))*
     .     amz**4/amw**2*Asbotsbotr(I,1,2)*Asbotsbotr(J,1,2)
     .*(-1.D0+x1+x2+xmusb1)
       endif
               ENDDO
       ENDDO   
      sb2sb1higgs = sb2sb1hk+sb2sb1aa
c -------------------------------------------------------------------- c
c                          chargino exchange
c -------------------------------------------------------------------- c
      sb2sb1chi=0.D0

      do k=1,2,1
         do l=1,2,1
      if((amchar(k)+amt).gt.asb2.and.(amchar(l)+amt).gt.asb2)then
            sb2sb1chi=sb2sb1chi+1.D0/dchi(k)/dchi(l)*(
     .          (alsbot2(1,k)*alsbot2(1,l)*alsbot1(2,k)*alsbot1(2,l)+
     .           blsbot2(1,k)*blsbot2(1,l)*blsbot1(2,k)*blsbot1(2,l))*
     .          ((1.D0-x1)*(1.D0-x2)-xmusb1+xmut*(x1-x2+xmusb1
     .           -2.D0*xmut)+xmut)
     .         +(alsbot2(1,k)*alsbot2(1,l)*blsbot1(2,k)*blsbot1(2,l)+
     .           blsbot2(1,k)*blsbot2(1,l)*alsbot1(2,k)*alsbot1(2,l))*
     .           xmchar(k)*xmchar(l)/asb2**2*(x1+x2-1.D0+xmusb1
     .           -2.D0*xmut)
     .         +(alsbot2(1,k)*blsbot2(1,l)*alsbot1(2,k)*blsbot1(2,l)+
     .           blsbot2(1,k)*alsbot2(1,l)*blsbot1(2,k)*alsbot1(2,l))*
     .           2.D0*xmut*(-1.D0+x1-xmut)
     .         +dsqrt(xmut)*(x1-2.D0*xmut)*(xmchar(k)/asb2*
     .          (alsbot2(1,k)*blsbot2(1,l)*blsbot1(2,k)*blsbot1(2,l)+
     .           blsbot2(1,k)*alsbot2(1,l)*alsbot1(2,k)*alsbot1(2,l))
     .         +xmchar(l)/asb2**
     .          (alsbot2(1,k)*blsbot2(1,l)*alsbot1(2,k)*alsbot1(2,l)+
     .           blsbot2(1,k)*alsbot2(1,l)*blsbot1(2,k)*blsbot1(2,l)) )
     .         +dsqrt(xmut)*(-1.D0+x1+xmusb1-2.D0*xmut)* 
     .         (xmchar(k)/asb2*
     .          (alsbot2(1,k)*alsbot2(1,l)*blsbot1(2,k)*alsbot1(2,l)+
     .           blsbot2(1,k)*blsbot2(1,l)*alsbot1(2,k)*blsbot1(2,l))
     .         +xmchar(l)/asb2*
     .          (alsbot2(1,k)*alsbot2(1,l)*alsbot1(2,k)*blsbot1(2,l)+
     .           blsbot2(1,k)*blsbot2(1,l)*blsbot1(2,k)*alsbot1(2,l)) )
     .         +(alsbot2(1,l)*blsbot2(1,k)*alsbot1(2,k)*blsbot1(2,l)+
     .           blsbot2(1,l)*alsbot2(1,k)*blsbot1(2,k)*alsbot1(2,l))*
     .           (-2.D0)*xmut*xmchar(k)/asb2*xmchar(l)/asb2 )
      endif
         end do
      end do
      
c -------------------------------------------------------------------- c
c                         chargino Z interference
c -------------------------------------------------------------------- c
      sb2sb1chiz=0.D0

      do l=1,2,1
      if((amchar(l)+amt).gt.asb2.and.(mz+asb1).gt.asb2)then
         sb2sb1chiz=sb2sb1chiz+azbb12/cw/dz/dchi(l)*(
     .      xmut*(alsbot2(1,l)*alsbot1(2,l)*(vzz-azz)+
     .            blsbot2(1,l)*blsbot1(2,l)*(vzz+azz))*(
     .      1.D0/xmuz*(xmusb1*(xmusb1-2.D0)+1.D0)+
     .      2.D0*x1-xmusb1-3.D0)
     .     +(alsbot2(1,l)*alsbot1(2,l)*(vzz+azz)+
     .       blsbot2(1,l)*blsbot1(2,l)*(vzz-azz))*(
     .      1.D0/xmuz*(1.D0-xmusb1)*(xmut*(xmusb1+x1-1.D0)-x1*xmut)
     .      +xmut*(1.D0+x1+xmusb1-2.D0*x2)-2.D0*xmusb1+2.D0*x1*x2+2.D0
     .      -2.D0*(x1+x2)+xmut*(-x1+2.D0) )
     .     +dsqrt(xmut)*xmchar(l)/asb2*
     .      (alsbot2(1,l)*blsbot1(2,l)*(vzz-azz)+
     .       blsbot2(1,l)*alsbot1(2,l)*(vzz+azz))*(
     .      1.D0/xmuz*(xmusb1-1.D0)*(xmusb1+x1+x2-1.D0)+x1-x2-xmusb1
     .      +1.D0 )
     .     +dsqrt(xmut)*xmchar(l)/asb2*
     .      (alsbot2(1,l)*blsbot1(2,l)*(vzz+azz)+
     .       blsbot2(1,l)*alsbot1(2,l)*(vzz-azz))*(
     .      1.D0/xmuz*(1.D0-xmusb1)*(xmusb1+x1+x2-1.D0)+x1-x2+xmusb1
     .      -1.D0 ) )
      endif
      end do
c -------------------------------------------------------------------- c
c                       chargino Higgs interference
c -------------------------------------------------------------------- c
      sb2sb1hchi=0.D0


      do l=1,2,1
          do I = 1,3
      if((amchar(l)+amt).gt.asb2.and.(smass(i)+asb1).gt.asb2)then
         sb2sb1hchi=sb2sb1hchi-2.D0*(Httr(I)/dsqrt(2.D0)/dchi(l)
     .        /dh(I)/asb2*amz**2/amw*Hsbotsbotr(I,1,2))*(
     .        (alsbot2(1,l)*blsbot1(2,l)+alsbot1(2,l)*blsbot2(1,l))*
     .        xmchar(l)/asb2*(x1+x2+xmusb1-1.D0-4.D0*xmut) 
     .        +(alsbot2(1,l)*alsbot1(2,l)+blsbot1(2,l)*blsbot2(1,l))*
     .        dsqrt(xmut)*(2.D0*x1-1.D0+xmusb1-4.D0*xmut) )
            endif
         enddo
         do J = 1,2
      if((amchar(l)+amt).gt.asb2.and.(pmass(J)+asb1).gt.asb2)then
            sb2sb1hchi=sb2sb1hchi
     .        +2.D0*Attr(J)/dsqrt(2.D0)/dchi(l)/da(J)/asb2*amz**2/amw*
     .        (-Asbotsbotr(I,1,2))*(
     .         (alsbot2(1,l)*alsbot1(2,l)-blsbot1(2,l)*blsbot2(1,l))*
     .        dsqrt(xmut)*(1.D0-xmusb1) +
     .        (alsbot2(1,l)*blsbot1(2,l)-alsbot1(2,l)*blsbot2(1,l))*
     .        xmchar(l)/asb2*(1.D0-x1-x2-xmusb1) )
       endif
         end do
       enddo
c -------------------------------------------------------------------- c
c                         Z Higgs interference
c -------------------------------------------------------------------- c
       sb2sb1hz=0.0d0

      DO I =1,3
         if((smass(i)+asb1).gt.asb2.and.(mz+asb1).gt.asb2)then
      sb2sb1hz=sb2sb1hz-2.D0/2.D0/cw*
     .   (azbb12*Httr(I)/dsqrt(2.D0)*amz**2/amw*Hsbotsbotr(I,1,2)
     .   /asb2/dz/dh(I))*2.D0*dsqrt(xmut)*vzz*2.D0*(x1-x2)
         endif
      ENDDO

       DO I =1,2
         if((pmass(i)+asb1).gt.asb2.and.(mz+asb1).gt.asb2)then
          sb2sb1hz=sb2sb1hz
     .   +2.D0/2.D0/cw*
     .    azbb12*Attr(I)/dsqrt(2.D0)*amz**2/amw*
     .   (-Asbotsbotr(I,1,2))/asb2/dz/da(I)*
     .   (2.D0*dsqrt(xmut)*azz*(2.D0/xmuz*(1.D0+(xmusb1-1.D0)*(x1
     .   +x2)+xmusb1**2-2.D0*xmusb1)+2.D0-2.D0*xmusb1) )
          endif
       ENDDO
c -------------------------------------------------------------------- c

      NS_sb2sb1tt=sb2sb1zz+sb2sb1higgs+sb2sb1chi+sb2sb1chiz+sb2sb1hchi+
     .         sb2sb1hz

      end
c ==================================================================== c
c ========================= sbottom1 up upbar ======================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sb2sb1uu(x1,x2)
*
      IMPLICIT NONE
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION azbb11,azbb12,azbb21,azbb22,xmuz,xmusb1,x3,x1,x2
     .,dz,azz,vzz,sb2sb1zz
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW

      azbb11=gzbb(1,1)
      azbb12=gzbb(1,2)
      azbb21=gzbb(2,1)
      azbb22=gzbb(2,2)

      xmuz   = mz**2/asb2**2
      xmusb1 = asb1**2/asb2**2

      x3  = 2.D0-x1-x2
      dz  = 1.D0-x3+xmusb1-xmuz

      vzz = vzztoptop
      azz = azztoptop
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      sb2sb1zz=0.0d0
      if((mz+asb1).gt.asb2)Then
      sb2sb1zz = 1.D0/4.D0/cw**2/dz**2*azbb12**2*
     .     2.D0*(vzz**2+azz**2)*4.D0*(1.D0+x1*x2-x1-x2-xmusb1) 
      endif
c ----------------------------------------------------------------- c
      NS_sb2sb1uu=sb2sb1zz

      end
c ==================================================================== c
c ======================= sbottom1 down downbar ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sb2sb1dd(x1,x2)
*
      IMPLICIT NONE
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION azbb11,azbb12,azbb21,azbb22,xmuz,xmusb1,x3,x1,x2
     .,dz,azz,vzz,sb2sb1zz
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW

      azbb11=gzbb(1,1)
      azbb12=gzbb(1,2)
      azbb21=gzbb(2,1)
      azbb22=gzbb(2,2)

      xmuz   = mz**2/asb2**2
      xmusb1 = asb1**2/asb2**2

      x3  = 2.D0-x1-x2
      dz  = 1.D0-x3+xmusb1-xmuz

      vzz = vzzbotbot
      azz = azzbotbot
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      sb2sb1zz=0.0d0
      if((mz+asb1).gt.asb2)Then
      sb2sb1zz = 1.D0/4.D0/cw**2/dz**2*azbb12**2*
     .     2.D0*(vzz**2+azz**2)*4.D0*(1.D0+x1*x2-x1-x2-xmusb1)
      endif
c -------------------------------------------------------------------- c

      NS_sb2sb1dd=sb2sb1zz

      end
c ==================================================================== c
c =========================== sbottom1 e+ e- ========================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sb2sb1ee(x1,x2)
*
      IMPLICIT NONE
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION azbb11,azbb12,azbb21,azbb22,xmuz,xmusb1,x3,x1,x2
     .,dz,azz,vzz,sb2sb1zz
*     
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW

      azbb11=gzbb(1,1)
      azbb12=gzbb(1,2)
      azbb21=gzbb(2,1)
      azbb22=gzbb(2,2)

      xmuz   = mz**2/asb2**2
      xmusb1 = asb1**2/asb2**2

      x3  = 2.D0-x1-x2
      dz  = 1.D0-x3+xmusb1-xmuz

      vzz = vzztautau
      azz = azztautau
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      sb2sb1zz=0.d0
      if((mz+asb1).gt.asb2)Then
      sb2sb1zz = 1.D0/4.D0/cw**2/dz**2*azbb12**2*
     .     2.D0*(vzz**2+azz**2)*4.D0*(1.D0+x1*x2-x1-x2-xmusb1)
      endif
c -------------------------------------------------------------------- c
      NS_sb2sb1ee=sb2sb1zz

      end
c ==================================================================== c
c ========================= sbottom1 nu nubar ======================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sb2sb1nunu(x1,x2)
*
      IMPLICIT NONE
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION azbb11,azbb12,azbb21,azbb22,xmuz,xmusb1,x3,x1,x2
     .,dz,azz,vzz,sb2sb1zz
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW

      azbb11=gzbb(1,1)
      azbb12=gzbb(1,2)
      azbb21=gzbb(2,1)
      azbb22=gzbb(2,2)

      xmuz   = mz**2/asb2**2
      xmusb1 = asb1**2/asb2**2

      x3  = 2.D0-x1-x2
      dz  = 1.D0-x3+xmusb1-xmuz

      vzz = vzzneutneut
      azz = azzneutneut
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      sb2sb1zz=0.d0
      if((mz+asb1).gt.asb2)Then
      sb2sb1zz = 1.D0/4.D0/cw**2/dz**2*azbb12**2*
     .     2.D0*(vzz**2+azz**2)*4.D0*(1.D0+x1*x2-x1-x2-xmusb1) 
      endif
c -------------------------------------------------------------------- c

      NS_sb2sb1nunu=sb2sb1zz

      end
c ==================================================================== c
c ========================= sbottom1 tau+ tau- ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sb2sb1tautau(x1,x2)
*
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS 
      DOUBLE PRECISION P(2,3)
      DOUBLE PRECISION xmuh(3),xmua(2),dh(3),da(2)
      DOUBLE PRECISION azbb11,azbb12,azbb21,azbb22,xmuz,xmutau,x3,x1,x2
     .,dz,azz,vzz,sb2sb1zz,xmusb1,sb2sb1hk,sb2sb1aa,sb2sb1higgs,
     .sb2sb1hz
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
*
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr 
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_CPodd_MIX/P
      COMMON/NS_MZMWscaleQ/AMZ,AMW

      azbb11=gzbb(1,1)
      azbb12=gzbb(1,2)
      azbb21=gzbb(2,1)
      azbb22=gzbb(2,2)

      xmuz   = mz**2/asb2**2
      xmutau = amtau**2/asb2**2
      xmusb1 = asb1**2/asb2**2

      DO i =1,3
      xmuh(i) = SMASS(i)**2/asb2**2
      enddo
      DO i =1,2
      xmua(i) = PMASS(i)**2/asb2**2
      enddo

      x3  = 2.D0-x1-x2
      dz  = 1.D0-x3+xmusb1-xmuz

      DO i =1,3
      dh(i) = 1.D0-x3+xmusb1-xmuh(i)
      enddo
      DO i =1,2
      da(i) = 1.D0-x3+xmusb1-xmua(i)
      enddo
      vzz = vzztautau
      azz = azztautau
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      sb2sb1zz=0.d0
      if((mz+asb1).gt.asb2)Then
      sb2sb1zz = 1.D0/4.D0/cw**2/dz**2*azbb12**2*(
     .   1.D0/xmuz**2*(4.D0*(vzz**2-azz**2)*xmutau*(1.D0-(x1+x2)*(1.D0  
     .   -xmusb1)**2+xmusb1*(-xmusb1**2+3.D0*xmusb1-3.D0))+
     .   2.D0*(vzz**2+azz**2)*(1.D0-xmusb1)**2*(2.D0*xmutau*(xmusb1
     .   +x1+x2-1.D0)) )
     .   +1.D0/xmuz*(8.D0*(vzz**2-azz**2)*xmutau*(1.D0-xmusb1)**2+
     .   4.D0*(vzz**2+azz**2)*(1.D0-xmusb1)*2.D0*xmutau*(xmusb1-1.D0))
     .   +4.D0*(vzz**2-azz**2)*xmutau*(x1+x2-xmusb1-3.D0)
     .   +2.D0*(vzz**2+azz**2)*(4.D0*(1.D0+x1*x2-x1-x2-xmusb1)+
     .   xmutau*(-2.D0*x1-2.D0*x2+2.D0*xmusb1+6.D0)) )
      endif
c -------------------------------------------------------------------- c
c                             Higgs exchange
c -------------------------------------------------------------------- c
      sb2sb1hk=0.0d0
      sb2sb1aa=0.0d0


      do I =1,3
         do J =1,3
      if((smass(i)+asb1).gt.asb2.and.(smass(j)+asb1).gt.asb2)Then
      sb2sb1hk = sb2sb1hk+(2.D0/asb2**2/dh(I)/dh(J)*
     .(-scaltau/dsqrt(2.D0)*(-S(i,2)))*(-scaltau/dsqrt(2.D0)*(-S(J,2)))*
     .    amz**4/amw**2*Hsbotsbotr(I,1,2)*Hsbotsbotr(J,1,2))*
     .    (-1.D0+x1+x2+xmusb1-4.D0*xmutau)
      endif
          ENDDO
      ENDDO
      do I =1,2
         do J =1,2
      if((pmass(i)+asb1).gt.asb2.and.(pmass(j)+asb1).gt.asb2)Then
      sb2sb1aa = sb2sb1aa+
     .     2.D0/asb2**2/da(I)/da(J)*
     . (scaltau/dsqrt(2.D0)*P(I,2))*(scaltau/dsqrt(2.D0)*P(J,2))*
     .     amz**4/amw**2*Asbotsbotr(I,1,2)*Asbotsbotr(J,1,2)*
     . (-1.D0+x1+x2+xmusb1)
      endif
         enddo
      enddo  

      sb2sb1higgs = sb2sb1hk+sb2sb1aa
c -------------------------------------------------------------------- c
c                           Higgs Z interference
c -------------------------------------------------------------------- c
      sb2sb1hz=0.0d0

       DO I = 1,3
          if((smass(i)+asb1).gt.asb2.and.(mz+asb1).gt.asb2)Then
       sb2sb1hz= sb2sb1hz-1.D0/cw*
     .   (azbb12*(-scaltau/dsqrt(2.D0)*(-S(i,2)))*
     . amz**2/amw*Hsbotsbotr(I,1,2)
     .   /asb2/dz/dh(I))*
     . 2.D0*dsqrt(xmutau)*vzz*2.D0*(x1-x2)
         endif
      enddo

      DO J = 1,2
        if((pmass(j)+asb1).gt.asb2.and.(mz+asb1).gt.asb2)Then 
         sb2sb1hz= sb2sb1hz
     .   + 1.D0/cw*
     .    azbb12*(-scaltau/dsqrt(2.D0)*(P(J,2)))*amz**2/amw*
     .   (-Asbotsbotr(J,1,2))/asb2/dz/da(j)*
     .   (2.D0*dsqrt(xmutau)*azz*(2.D0/xmuz*(1.D0+(xmusb1-1.D0)*(x1
     .   +x2)+xmusb1**2-2.D0*xmusb1)+2.D0-2.D0*xmusb1) )
        endif
      ENDDO
c -------------------------------------------------------------------- c

      NS_sb2sb1tautau=sb2sb1zz+sb2sb1higgs+sb2sb1hz

      end
