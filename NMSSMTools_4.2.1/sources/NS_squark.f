c ==================================================================== c
c                        up/down squark 2-body decays                         c
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
      SUBROUTINE NS_SQUARKS
*
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION NS_lamb
      DOUBLE PRECISION amuv,lamv
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION flagmulti,flagqcd,flagloop
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION SCALb,SCALt,scaltau,gs2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION scala,alp,ca,cf,amsq,rval
      DOUBLE PRECISION amurefer
      DOUBLE PRECISION suplneutup(5),suprneutup(5),suplchardow(2),
     .          suprchardow(2),qcdsuplneutup(5),qcdsuprneutup(5),
     .          qcdsuplchardow(2),qcdsuprchardow(2)
      DOUBLE PRECISION sdowlneutdow(5),sdowlcharup(2),sdowrneutdow(5),
     .          sdowrcharup(2),qcdsdowlneutdow(5),qcdsdowlcharup(2),
     .          qcdsdowrneutdow(5),qcdsdowrcharup(2)
      DOUBLE PRECISION supltot2,suprtot2,sdowltot2,sdowrtot2
      DOUBLE PRECISION brsuplnup(5),brsuplcdow(2),brsuplglui,
     .          brsuprnup(5),brsuprcdow(2),brsuprglui,
     .          brsdowlndow(5),brsdowlchup(2),brsdowlglui,
     .          brsdowrndow(5),brsdowrchup(2),brsdowrglui
      DOUBLE PRECISION supltot2lo,suprtot2lo,supltot2nlo,suprtot2nlo,
     .          sdowltot2lo,sdowrtot2lo,sdowltot2nlo,sdowrtot2nlo
      DOUBLE PRECISION suplglui,qcdsuplglui,qcdsuprglui,
     . sdowlglui,sdowrglui,qcdsdowlglui,qcdsdowrglui,suprglui
      DOUBLE PRECISION NS_ftotqcd,NS_gama,NS_gamfcap,NS_gamf,
     . NS_gamglui2,NS_gamrendec,resum
*
      COMMON/SQUARK_2GAMMA/suplneutup,suprneutup,suplchardow,
     .          suprchardow,qcdsuplneutup,qcdsuprneutup,
     .          qcdsuplchardow,qcdsuprchardow,
     .          sdowlneutdow,sdowlcharup,sdowrneutdow,
     .          sdowrcharup,qcdsdowlneutdow,qcdsdowlcharup,
     .          qcdsdowrneutdow,qcdsdowrcharup
      COMMON/SQUARK_WIDTH/supltot2,suprtot2,sdowltot2,sdowrtot2
      COMMON/SQUARK_BR_2BD/brsuplnup,brsuplcdow,brsuplglui,
     .          brsuprnup,brsuprcdow,brsuprglui,
     .          brsdowlndow,brsdowlchup,brsdowlglui,
     .          brsdowrndow,brsdowrchup,brsdowrglui
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_coup7/alup,aldo
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_refscale/amurefer
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_qcdscales/amuv,lamv
*
      EXTERNAL NS_lamb
      EXTERNAL NS_ftotqcd
      EXTERNAL NS_gama,NS_gamfcap,NS_gamf,NS_gamglui2,NS_gamrendec
      EXTERNAL resum

c -- initialization --

      supltot2 = 0.D0
      suprtot2 = 0.D0
      supltot2lo = 0.D0
      suprtot2lo = 0.D0
      supltot2nlo = 0.D0
      suprtot2nlo = 0.D0
      do i=1,5
         suplneutup(i) = 0.D0
         suprneutup(i) = 0.D0
         qcdsuplneutup(i) = 0.D0
         qcdsuprneutup(i) = 0.D0
      enddo
      do i=1,2
         suplchardow(i) = 0.D0
         suprchardow(i) = 0.D0
         qcdsuplchardow(i) = 0.D0
         qcdsuprchardow(i) = 0.D0
      enddo
      suplglui = 0.D0
      suprglui = 0.D0
      qcdsuplglui = 0.D0
      qcdsuprglui = 0.D0

      sdowltot2 = 0.D0
      sdowrtot2 = 0.D0
      sdowltot2lo = 0.D0
      sdowrtot2lo = 0.D0
      sdowltot2nlo = 0.D0
      sdowrtot2nlo = 0.D0
      do i=1,5
         sdowlneutdow(i) = 0.D0
         sdowrneutdow(i) = 0.D0
         qcdsdowlneutdow(i) = 0.D0
         qcdsdowrneutdow(i) = 0.D0
      enddo
      do i=1,2
         sdowlcharup(i) = 0.D0
         sdowrcharup(i) = 0.D0
         qcdsdowlcharup(i) = 0.D0
         qcdsdowrcharup(i) = 0.D0
      enddo
      sdowlglui = 0.D0
      sdowrglui = 0.D0
      qcdsdowlglui = 0.D0
      qcdsdowrglui = 0.D0
c -------------------------------------------------------------------- c
c  supl -> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + up
      do i=1,5
         if(amneut(i).le.asup1) then
            suplneutup(i)=g2s*(aup(1,i)**2+bup(1,i)**2)*(asup1**2
     .           -amneut(i)**2)*NS_lamb(0.d0,amneut(i)/asup1)
     .           /(16*pi*asup1)
         else
            suplneutup(i)=0.d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do i=1,5
         if(amneut(i).le.asup1) then
            qcdsuplneutup(i)=4.D0/3.D0*gs2/(4.D0*pi)/pi*
     .           suplneutup(i)*
     .           NS_ftotqcd(amneut(i)**2/asup1**2,mgluino**2/asup1**2)
         else
            qcdsuplneutup(i)=0.d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  supl -> chi+_1/chi+_2 + down
      do i=1,2
         if (amchar(i).le.asup1) then
            suplchardow(i)=g2s*alup(1,i)**2*(asup1**2-amchar(i)**2)*
     .           NS_lamb(0.d0,amchar(i)/asup1)/(16*pi*asup1)
         else
            suplchardow(i)=0.d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do i=1,2
         if(amchar(i).le.asup1) then
            qcdsuplchardow(i)=4.D0/3.D0*gs2/(4.D0*pi)/pi*
     .           suplchardow(i)*
     .           NS_ftotqcd(amchar(i)**2/asup1**2,mgluino**2/asup1**2)
         else
            qcdsuplchardow(i)=0.d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  supr -> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + up
      do i=1,5
         if(amneut(i).le.asup2) then
            suprneutup(i)=g2s*(aup(2,i)**2+bup(2,i)**2)*(asup2**2
     .           -amneut(i)**2)*NS_lamb(0.d0,amneut(i)/asup2)
     .           /(16*pi*asup2)
         else
            suprneutup(i)=0.d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do i=1,5
         if(amneut(i).le.asup2) then
            qcdsuprneutup(i)=4.D0/3.D0*gs2/(4.D0*pi)/pi*
     .           suprneutup(i)*
     .           NS_ftotqcd(amneut(i)**2/asup2**2,mgluino**2/asup2**2)
         else
            qcdsuprneutup(i)=0.d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  supr -> chi+_1/chi+_2 + down
      do i=1,2
         if (amchar(i).le.asup2) then
            suprchardow(i)=g2s*alup(2,i)**2*(asup2**2-amchar(i)**2)*
     .           NS_lamb(0.d0,amchar(i)/asup2)/(16*pi*asup2)
         else
            suprchardow(i)=0.d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do i=1,2
         if(amchar(i).le.asup2) then
            qcdsuprchardow(i)=4.D0/3.D0*gs2/(4.D0*pi)/pi*
     .           suprchardow(i)*
     .           NS_ftotqcd(amchar(i)**2/asup2**2,mgluino**2/asup2**2)
         else
            qcdsuprchardow(i)=0.d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  supl --> gluino + up
      if(asup1.gt.mgluino) then
         suplglui = 8.D0*gs2*(asup1**2-mgluino**2)*
     .              NS_lamb(0.D0,mgluino/asup1)/(16.D0*pi*asup1)/3.D0 
      else
         suplglui = 0.D0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asup1.gt.mgluino) then
         scala = amurefer
         alp   = gs2/(4.D0*pi)
         ca    = 3.D0
         cf    = 4.D0/3.D0
         amsq  = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         rval  = mgluino**2/amsq**2

         qcdsuplglui = suplglui*alp/pi*( ca*NS_gama(rval) + 
     .        cf*NS_gamfcap(rval) + 4.D0*NS_gamf(rval) + 
     .        2.D0*pi**2/mgluino**2*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,mgluino,
     .        1,amsq) + 
     .        NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,scala) )
      else
         qcdsuplglui = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c  supr --> gluino + up
      
      if(asup2.gt.mgluino) then
         suprglui = 8.D0*gs2*(asup2**2-mgluino**2)*
     .              NS_lamb(0.D0,mgluino/asup2)/(16.D0*pi*asup2)/3.D0
      else
         suprglui = 0.D0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asup2.gt.mgluino) then
         scala = amurefer
         alp   = gs2/(4.D0*pi)
         ca    = 3.D0
         cf    = 4.D0/3.D0
         amsq  = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         rval  = mgluino**2/amsq**2

         qcdsuprglui = suprglui*alp/pi*( ca*NS_gama(rval) + 
     .        cf*NS_gamfcap(rval) + 4.D0*NS_gamf(rval) + 
     .        2.D0*pi**2/mgluino**2*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,mgluino,
     .        2,amsq) + 
     .        NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,scala) )
      else
         qcdsuprglui = 0.D0
      endif
      endif
c ------------------------------------------------------------------- c
       supltot2lo = suplneutup(1)+suplneutup(2)+suplneutup(3)
     .           + suplneutup(4) + suplneutup(5)
     .           + suplchardow(1)+suplchardow(2)+suplglui

      suprtot2lo = suprneutup(1)+suprneutup(2)+suprneutup(3)
     .           + suprneutup(4)+suprneutup(5) 
     .           + suprchardow(1)+suprchardow(2) +suprglui

      supltot2nlo = supltot2lo + qcdsuplneutup(1)+qcdsuplneutup(2)
     .            + qcdsuplneutup(3)+qcdsuplneutup(4)+qcdsuplneutup(5)
     .            +qcdsuplchardow(1)+ qcdsuplchardow(2)+qcdsuplglui

      suprtot2nlo = suprtot2lo + qcdsuprneutup(1)+qcdsuprneutup(2)
     .            + qcdsuprneutup(3)+qcdsuprneutup(4)+qcdsuprneutup(5)
     .            + qcdsuprchardow(1)+ qcdsuprchardow(2)+qcdsuprglui

      if(flagqcd.eq.0.D0) then
         supltot2 = supltot2lo
         suprtot2 = suprtot2lo
      elseif(flagqcd.eq.1.D0) then
         supltot2 = supltot2nlo
         suprtot2 = suprtot2nlo
      endif

      if(flagqcd.eq.1.D0) then
         do i=1,5
            suplneutup(i) = suplneutup(i)+qcdsuplneutup(i)
            suprneutup(i) = suprneutup(i)+qcdsuprneutup(i)
         enddo
         do i=1,2
            suplchardow(i) = suplchardow(i)+qcdsuplchardow(i)
            suprchardow(i) = suprchardow(i)+qcdsuprchardow(i)
         enddo
         suplglui = suplglui+qcdsuplglui
         suprglui = suprglui+qcdsuprglui
      endif

      do i=1,5
         brsuplnup(i) = suplneutup(i)/supltot2
         brsuprnup(i) = suprneutup(i)/suprtot2
      enddo
      do i=1,2
         brsuplcdow(i) = suplchardow(i)/supltot2
         brsuprcdow(i) = suprchardow(i)/suprtot2
      enddo
      brsuplglui = suplglui/supltot2
      brsuprglui = suprglui/suprtot2
      
c -------------------------------------------------------------------- c
c  sdownl -> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + down

      do i=1,5
         if(amneut(i).le.asdown1) then
            sdowlneutdow(i)=g2s*((ado(1,i)**2+bdo(1,i)**2)*
     .           (asdown1**2-amneut(i)**2)
     .           )*NS_lamb(0.d0,amneut(i)/asdown1)
     .           /(16*pi*asdown1)
         else
            sdowlneutdow(i)=0.D0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do i=1,5
         if(amneut(i).le.asdown1) then
            qcdsdowlneutdow(i)=4.D0/3.D0*gs2/(4.D0*pi)/pi*
     .         sdowlneutdow(i)*
     .         NS_ftotqcd(amneut(i)**2/asdown1**2,mgluino**2/asdown1**2)
         else
            qcdsdowlneutdow(i)=0.d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sdownl -> chi-_1/chi-_2 up

      do i=1,2
         if(amchar(i).le.asdown1) then
            sdowlcharup(i)=g2s*aldo(1,i)**2*
     .           (asdown1**2-amchar(i)**2)*
     .           NS_lamb(0.d0,amchar(i)/asdown1)/(16*pi*asdown1)
         else
            sdowlcharup(i)=0.D0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do i=1,2
         if(amchar(i).le.asdown1) then
            qcdsdowlcharup(i)=4.D0/3.D0*gs2/(4.D0*pi)/pi*
     .         sdowlcharup(i)*
     .         NS_ftotqcd(amchar(i)**2/asdown1**2,mgluino**2/asdown1**2)
         else
            qcdsdowlcharup(i)=0.D0
         endif
      end do
      endif
c -------------------------------------------------------------------- c
c  sdownr -> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + down

      do i=1,5
         if(amneut(i).le.asdown2) then
            sdowrneutdow(i)=g2s*((ado(2,i)**2+bdo(2,i)**2)*
     .           (asdown2**2-amneut(i)**2)
     .           )*NS_lamb(0.d0,amneut(i)/asdown2)
     .           /(16*pi*asdown2)
         else
            sdowrneutdow(i)=0.D0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do i=1,5
         if(amneut(i).le.asdown2) then
            qcdsdowrneutdow(i)=4.D0/3.D0*gs2/(4.D0*pi)/pi*
     .         sdowrneutdow(i)*
     .         NS_ftotqcd(amneut(i)**2/asdown2**2,mgluino**2/asdown2**2)
         else
            qcdsdowrneutdow(i)=0.d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sdownr -> chi-_1/chi-_2 up

      do i=1,2
         if(amchar(i).le.asdown2) then
            sdowrcharup(i)=g2s*aldo(2,i)**2*
     .           (asdown2**2-amchar(i)**2)*
     .           NS_lamb(0.d0,amchar(i)/asdown2)/(16*pi*asdown2)
         else
            sdowrcharup(i)=0.D0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      do i=1,2
         if(amchar(i).le.asdown2) then
            qcdsdowrcharup(i)=4.D0/3.D0*gs2/(4.D0*pi)/pi*
     .         sdowrcharup(i)*
     .         NS_ftotqcd(amchar(i)**2/asdown2**2,mgluino**2/asdown2**2)
         else
            qcdsdowrcharup(i)=0.D0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sdownl --> gluino + down

      if(asdown1.gt.mgluino) then
         sdowlglui = 8.D0*gs2*(asdown1**2-mgluino**2)*
     .            NS_lamb(0.D0,mgluino/asdown1)/(16.D0*pi*asdown1)/3.D0 
      else
         sdowlglui = 0.D0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asdown1.gt.mgluino) then
         scala = amurefer
         alp   = gs2/(4.D0*pi)
         ca    = 3.D0
         cf    = 4.D0/3.D0
         amsq  = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         rval  = mgluino**2/amsq**2

         qcdsdowlglui = sdowlglui*alp/pi*( ca*NS_gama(rval) + 
     .        cf*NS_gamfcap(rval) + 4.D0*NS_gamf(rval) + 
     .        2.D0*pi**2/mgluino**2*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,mgluino,
     .        1,amsq) + 
     .        NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,scala) )
      else
         qcdsdowlglui = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c  sdownr --> gluino + down
      
      if(asdown2.gt.mgluino) then
         sdowrglui = 8.D0*gs2*(asdown2**2-mgluino**2)*
     .             NS_lamb(0.D0,mgluino/asdown2)/(16.D0*pi*asdown2)/3.D0
      else
         sdowrglui = 0.D0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asdown2.gt.mgluino) then
         scala = amurefer
         alp   = gs2/(4.D0*pi)
         ca    = 3.D0
         cf    = 4.D0/3.D0
         amsq  = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         rval  = mgluino**2/amsq**2

         qcdsdowrglui = sdowrglui*alp/pi*( ca*NS_gama(rval) + 
     .        cf*NS_gamfcap(rval) + 4.D0*NS_gamf(rval) + 
     .        2.D0*pi**2/mgluino**2*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,mgluino,
     .        2,amsq) + 
     .        NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,scala) )
      else
         qcdsdowrglui = 0.D0
      endif
      endif
c ----------------------------
      sdowltot2lo = sdowlneutdow(1)+sdowlneutdow(2)+sdowlneutdow(3)+
     .              sdowlneutdow(4)+sdowlneutdow(5)+
     .              sdowlcharup(1)+sdowlcharup(2)+sdowlglui
 
      sdowrtot2lo = sdowrneutdow(1)+sdowrneutdow(2)+sdowrneutdow(3)+
     .              sdowrneutdow(4)+sdowrneutdow(5)+
     .              sdowrcharup(1)+sdowrcharup(2)+sdowrglui

      sdowltot2nlo = sdowltot2lo + qcdsdowlneutdow(1)+
     .               qcdsdowlneutdow(2)+qcdsdowlneutdow(3)+
     .               qcdsdowlneutdow(4)+qcdsdowlneutdow(5)+
     .               qcdsdowlcharup(1)+qcdsdowlcharup(2)+qcdsdowlglui
 
      sdowrtot2nlo = sdowrtot2lo + qcdsdowrneutdow(1)+
     .               qcdsdowrneutdow(2)+qcdsdowrneutdow(3)+
     .               qcdsdowrneutdow(4)+qcdsdowrneutdow(5)+
     .               qcdsdowrcharup(1)+qcdsdowrcharup(2)+qcdsdowrglui

      if(flagqcd.eq.0.D0) then
         sdowltot2 = sdowltot2lo
         sdowrtot2 = sdowrtot2lo
      elseif(flagqcd.eq.1.D0) then
         sdowltot2 = sdowltot2nlo
         sdowrtot2 = sdowrtot2nlo
      endif

      if(flagqcd.eq.1.D0) then
c UE: Resum if qcdcorr < -tree:
         do i=1,5
      If(qcdsdowlneutdow(i).lt.-sdowlneutdow(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sdl->d+chi0_",i
       qcdsdowlneutdow(i)=resum(sdowlneutdow(i),qcdsdowlneutdow(i))
            sdowlneutdow(i) = sdowlneutdow(i)+qcdsdowlneutdow(i)
      If(qcdsdowrneutdow(i).lt.-sdowrneutdow(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sdr->d+chi0_",i
       qcdsdowrneutdow(i)=resum(sdowrneutdow(i),qcdsdowrneutdow(i))
            sdowrneutdow(i) = sdowrneutdow(i)+qcdsdowrneutdow(i)
         enddo
         do i=1,2
      If(qcdsdowlcharup(i).lt.-sdowlcharup(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sdl->u+char_",i
       qcdsdowlcharup(i)=resum(sdowlcharup(i),qcdsdowlcharup(i))
            sdowlcharup(i) = sdowlcharup(i)+qcdsdowlcharup(i)
      If(qcdsdowrcharup(i).lt.-sdowrcharup(i)) 
     .write(*,23)"Warning: large negative rad. corrs. to sdr->u+char_",i
       qcdsdowrcharup(i)=resum(sdowrcharup(i),qcdsdowrcharup(i))
            sdowrcharup(i) = sdowrcharup(i)+qcdsdowrcharup(i)
         enddo
      If(qcdsdowlglui.lt.-sdowlglui) 
     .write(*,23)"Warning: large negative rad. corrs. to sdl->d+gluino"
       qcdsdowlglui=resum(sdowlglui,qcdsdowlglui)
         sdowlglui = sdowlglui+qcdsdowlglui
      If(qcdsdowrglui.lt.-sdowrglui) 
     .write(*,23)"Warning: large negative rad. corrs. to sdr->d+gluino"
       qcdsdowrglui=resum(sdowrglui,qcdsdowrglui)
         sdowrglui = sdowrglui+qcdsdowrglui
      endif
23     format(A,I1)

      do i=1,5
         brsdowlndow(i) = sdowlneutdow(i)/sdowltot2
         brsdowrndow(i) = sdowrneutdow(i)/sdowrtot2
      enddo
      do i=1,2
         brsdowlchup(i) = sdowlcharup(i)/sdowltot2
         brsdowrchup(i) = sdowrcharup(i)/sdowrtot2
      enddo
      brsdowlglui = sdowlglui/sdowltot2
      brsdowrglui = sdowrglui/sdowrtot2

      END

