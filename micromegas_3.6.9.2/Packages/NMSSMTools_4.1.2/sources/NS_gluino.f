c ==================================================================== c
c                        gluino decays                                 c
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
      SUBROUTINE NS_GLUINO
*
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION amuv,lamv
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION  nf
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     .         gdr(2),gdl(2)
      DOUBLE PRECISION qscal,amuref,multilim,flagmulti,flagqcd,flagloop
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION gst1,gst2,gsb1,gsb2,gsupl,gsupr,
     .         gsdownl,gsdownr,qcdgst1,qcdgst2,qcdgsb1,qcdgsb2,
     .         qcdgsupl,qcdgsupr,qcdgsdownl,qcdgsdownr
      DOUBLE PRECISION glnjgluon(5)
      DOUBLE PRECISION xintegoup(5),xintegodn(5),xintegotp(5),
     .         xintegobt(5),xintegoud(2),xintegotb(2),xintegocc(2),
     .         xinteghcst1b,xintegwst1b      
      DOUBLE PRECISION gluitot,gluitot2,gluitotmulti,gluitotrad 
      DOUBLE PRECISION brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon(5)
      DOUBLE PRECISION brgoup(5),brgoch(5),brgodn(5),brgost(5),
     .         brgotp(5),
     .         brgobt(5),brgoud(2),brgocs(2),brgotb(2),brhcst1b,brwst1b
      DOUBLE PRECISION gluitot2lo,gluitot2nlo
      DOUBLE PRECISION amsq,alp,ca,cf,rval
      DOUBLE PRECISION NS_lamb,resum 
      DOUBLE PRECISION NS_gamtop1,NS_gamtop2,NS_gamglui1,NS_gamglui2,
     .NS_gamglui3,NS_gam11,NS_gam12,NS_gamvirt,NS_gamrealgl,
     .NS_gamcfdec
      DOUBLE PRECISION NS_gama,NS_gamfcap,NS_gamf,NS_gamrendec
      COMPLEX*16 NS_iint,NS_i2int,NS_jint,NS_kint
*
      COMMON/GLUINO_3GAMMA/xintegoup,xintegodn,xintegotp,
     .         xintegobt,xintegoud,xintegotb,xintegocc,
     .         xinteghcst1b,xintegwst1b
      COMMON/GLUINO_WIDTH/gluitot,gluitot2,gluitotmulti,gluitotrad
      COMMON/GLUINO_BR_2BD/brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon
      COMMON/GLUINO_BR_3BD/brgoup,brgoch,brgodn,brgost,brgotp,
     .         brgobt,brgoud,brgocs,brgotb,brhcst1b,brwst1b
*             
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_refscale/qscal
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_multilim/multilim
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop
*    
      EXTERNAL NS_iint,NS_i2int,NS_jint,NS_kint
      EXTERNAL NS_gamtop1,NS_gamtop2,NS_gamglui1,NS_gamglui2,
     .         NS_gamglui3,NS_gam11,NS_gam12,NS_gamvirt,NS_gamrealgl,
     .         NS_gamcfdec
      EXTERNAL NS_gama,NS_gamfcap,NS_gamf,NS_gamrendec
      EXTERNAL resum
*
c -- INITIALIZE ALL PARAMETERS USED 
      gluitot     = 0.D0
      gluitot2    = 0.D0
      gluitot2lo  = 0.D0
      gluitot2nlo = 0.D0
      gluitotmulti= 0.D0
      gluitotrad  = 0.D0
* 
      gst1    = 0.D0
      gst2    = 0.D0
      gsb1    = 0.D0
      gsb2    = 0.D0
      gsupl   = 0.D0
      gsupr   = 0.D0
      gsdownl = 0.D0
      gsdownr = 0.D0
*
      qcdgst1    = 0.D0
      qcdgst2    = 0.D0
      qcdgsb1    = 0.D0
      qcdgsb2    = 0.D0
      qcdgsupl   = 0.D0
      qcdgsupr   = 0.D0
      qcdgsdownl = 0.D0
      qcdgsdownr = 0.D0
*
      DO i=1,5
         xintegoup(i) = 0.D0
         xintegodn(i) = 0.D0
         xintegotp(i) = 0.D0
         xintegobt(i) = 0.D0
         glnjgluon(i) = 0.D0
      ENDDO
*
      DO i=1,2
         xintegoud(i) = 0.D0
         xintegotb(i) = 0.D0
      ENDDO
*
      xinteghcst1b = 0.D0
      xintegwst1b  = 0.D0
*
c -------------------------------------------------------------------- c
c For QCD corrections: the fixed scale is qscal = Q, where
c the couplings are defined:
      amuref = qscal
c -------------------------------------------------------------------- c

c -- begin 2 body decays
c
c  gluino --> stop1 + top

      if((ast1+amt).le.mgluino) then
         gst1=gs2/2.D0*((gtl(1)**2+gtr(1)**2)*
     .        (mgluino**2-ast1**2+amt**2)
     .        +4*gtl(1)*gtr(1)*mgluino*amt)*
     .        NS_lamb(amt/mgluino,ast1/mgluino)
     .         /(16*pi*mgluino)
      else
         gst1=0.D0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(mgluino.gt.(ast1+amt)) then
         qcdgst1 = 0.D0
         amsq    = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         alp     = gs2/(4.D0*pi)
         nf      = 6.D0
         qcdgst1 = 8.D0*pi*alp/3.D0/amt**2*gst1*
     .        NS_gamtop1(ast1,ast2,amt,mgluino,thet,1,amuv,lamv) -
     .        pi*alp**2*NS_lamb(amt/mgluino,ast1/mgluino)/
     .        (3.D0*amt**2*mgluino)*
     .        NS_gamtop2(ast1,ast2,amt,mgluino,thet,1,amuv) +
     .        4.D0*pi*alp/mgluino**2*gst1*(nf-2.D0)*
     .        NS_gamglui1(ast1,ast2,amsq,amt,mgluino,amuv) +
     .        2.D0*pi*alp/mgluino**2*gst1*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,
     .                    mgluino,1,amuv) +
     .        4.D0*pi*alp*3.D0/mgluino**2*gst1*
     .        NS_gamglui3(mgluino,amuv,lamv) +
     .        8.D0*4.D0/3.D0*pi*alp*gst1*
     .        NS_gam11(ast1,ast2,amt,mgluino,thet,1,amuv,lamv) -
     .        8.D0/3.D0*pi*alp**2/mgluino*
     .        NS_lamb(amt/mgluino,ast1/mgluino)*
     .        NS_gam12(ast1,ast2,amt,mgluino,thet,1,amuv,lamv,amuref) -
     .        3.D0/16.D0*alp**2*NS_lamb(amt/mgluino,ast1/mgluino)
     .        /mgluino*
     .        NS_gamvirt(ast1,ast2,amt,mgluino,thet,1,amuv,lamv) -
     .        3.D0/16.D0*ast1/mgluino*
     .        alp**2*NS_gamrealgl(ast1,amt,mgluino,thet,1,lamv) +
     .        alp/(4.D0*pi)*gst1*
     .        NS_gamcfdec(ast1,ast2,amt,asb1,asb2,amb,mgluino,amsq,amuv,
     .        amuref)

      else
         qcdgst1 = 0.D0
      endif
      endif
C      write(*,*)'qcdst1',gst1,qcdgst1,mgluino,amuref
c -------------------------------------------------------------------- c
c  ---gluino --> stop2 + top

      if((ast2+amt).le.mgluino) then
         gst2=gs2/2.d0*((gtl(2)**2+gtr(2)**2)*
     .        (mgluino**2-ast2**2+amt**2)
     .        +4*gtl(2)*gtr(2)*mgluino*amt)*
     .        NS_lamb(amt/mgluino,ast2/mgluino)
     .         /(16*pi*mgluino)
      else
         gst2=0.d0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(mgluino.gt.(ast2+amt)) then
         qcdgst2 = 0.D0
         amsq    = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         alp     = gs2/(4.D0*pi)
         nf      = 6.D0
         qcdgst2 = 8.D0*pi*alp/3.D0/amt**2*gst2*
     .        NS_gamtop1(ast2,ast1,amt,mgluino,thet,2,amuv,lamv) -
     .        pi*alp**2*NS_lamb(amt/mgluino,ast2/mgluino)/
     .        (3.D0*amt**2*mgluino)*
     .        NS_gamtop2(ast2,ast1,amt,mgluino,thet,2,amuv) +
     .        4.D0*pi*alp/mgluino**2*gst2*(nf-2.D0)*
     .        NS_gamglui1(ast2,ast1,amsq,amt,mgluino,amuv) +
     .        2.D0*pi*alp/mgluino**2*gst2*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,
     .                    mgluino,2,amuv) +
     .        4.D0*pi*alp*3.D0/mgluino**2*gst2*
     .        NS_gamglui3(mgluino,amuv,lamv) +
     .        8.D0*4.D0/3.D0*pi*alp*gst2*
     .        NS_gam11(ast2,ast1,amt,mgluino,thet,2,amuv,lamv) -
     .        8.D0/3.D0*pi*alp**2/mgluino*
     .        NS_lamb(amt/mgluino,ast2/mgluino)*
     .        NS_gam12(ast2,ast1,amt,mgluino,thet,2,amuv,lamv,amuref) -
     .        3.D0/16.D0*alp**2*NS_lamb(amt/mgluino,ast2/mgluino)
     .        /mgluino*
     .        NS_gamvirt(ast2,ast1,amt,mgluino,thet,2,amuv,lamv) -
     .        3.D0/16.D0*ast2/mgluino*
     .        alp**2*NS_gamrealgl(ast2,amt,mgluino,thet,2,lamv) +
     .        alp/(4.D0*pi)*gst2*
     .        NS_gamcfdec(ast2,ast1,amt,asb2,asb1,amb,mgluino,amsq,amuv,
     .        amuref)

      else
         qcdgst2 = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c --- gluino --> sbottom1 + bottom

      if((asb1+amb).le.mgluino) then
         gsb1=gs2/2.d0*((gbl(1)**2+gbr(1)**2)*
     .        (mgluino**2-asb1**2+amb**2)
     .        +4*gbl(1)*gbr(1)*mgluino*amb)*
     .        NS_lamb(amb/mgluino,asb1/mgluino)
     .         /(16*pi*mgluino)
      else
         gsb1=0.d0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
       if(mgluino.gt.(asb1+amb)) then
         qcdgsb1 = 0.D0
         amsq    = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         alp     = gs2/(4.D0*pi)
         nf      = 6.D0
         qcdgsb1 = 8.D0*pi*alp/3.D0/amb**2*gsb1*
     .        NS_gamtop1(asb1,asb2,amb,mgluino,theb,1,amuv,lamv) -
     .        pi*alp**2*NS_lamb(amb/mgluino,asb1/mgluino)/
     .        (3.D0*amb**2*mgluino)*
     .        NS_gamtop2(asb1,asb2,amb,mgluino,theb,1,amuv) +
     .        4.D0*pi*alp/mgluino**2*gsb1*(nf-2.D0)*
     .        NS_gamglui1(ast1,ast2,amsq,amt,mgluino,amuv) +
     .        2.D0*pi*alp/mgluino**2*gsb1*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,
     .                    mgluino,1,amuv) +
     .        4.D0*pi*alp*3.D0/mgluino**2*gsb1*
     .        NS_gamglui3(mgluino,amuv,lamv) +
     .        8.D0*4.D0/3.D0*pi*alp*gsb1*
     .        NS_gam11(asb1,asb2,amb,mgluino,theb,1,amuv,lamv) -
     .        8.D0/3.D0*pi*alp**2/mgluino*
     .        NS_lamb(amb/mgluino,asb1/mgluino)*
     .        NS_gam12(asb1,asb2,amb,mgluino,theb,1,amuv,lamv,amuref) -
     .        3.D0/16.D0*alp**2*NS_lamb(amb/mgluino,asb1/mgluino)
     .        /mgluino*
     .        NS_gamvirt(asb1,asb2,amb,mgluino,theb,1,amuv,lamv) -
     .        3.D0/16.D0*asb1/mgluino*
     .        alp**2*NS_gamrealgl(asb1,amb,mgluino,theb,1,lamv) +
     .        alp/(4.D0*pi)*gsb1*
     .        NS_gamcfdec(ast1,ast2,amt,asb1,asb2,amb,mgluino,amsq,amuv,
     .        amuref)
      else
         qcdgsb1 = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c ---  gluino --> sbottom2 + bottom

      if((asb2+amb).le.mgluino) then
         gsb2=gs2/2.d0*((gbl(2)**2+gbr(2)**2)*
     .        (mgluino**2-asb2**2+amb**2)
     .        +4*gbl(2)*gbr(2)*mgluino*amb)*
     .        NS_lamb(amb/mgluino,asb2/mgluino)
     .         /(16*pi*mgluino)
      else
         gsb2=0.d0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(mgluino.gt.(asb2+amb)) then
         qcdgsb2 = 0.D0
         amsq    = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         alp     = gs2/(4.D0*pi)
         nf      = 6.D0
         qcdgsb2 = 8.D0*pi*alp/3.D0/amb**2*gsb2*
     .        NS_gamtop1(asb2,asb1,amb,mgluino,theb,2,amuv,lamv) -
     .        pi*alp**2*NS_lamb(amb/mgluino,asb2/mgluino)/
     .        (3.D0*amb**2*mgluino)*
     .        NS_gamtop2(asb2,asb1,amb,mgluino,theb,2,amuv) +
     .        4.D0*pi*alp/mgluino**2*gsb2*(nf-2.D0)*
     .        NS_gamglui1(ast2,ast1,amsq,amt,mgluino,amuv) +
     .        2.D0*pi*alp/mgluino**2*gsb2*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,
     .                    mgluino,2,amuv) +
     .        4.D0*pi*alp*3.D0/mgluino**2*gsb2*
     .        NS_gamglui3(mgluino,amuv,lamv) +
     .        8.D0*4.D0/3.D0*pi*alp*gsb2*
     .        NS_gam11(asb2,asb1,amb,mgluino,theb,2,amuv,lamv) -
     .        8.D0/3.D0*pi*alp**2/mgluino*
     .        NS_lamb(amb/mgluino,asb2/mgluino)*
     .        NS_gam12(asb2,asb1,amb,mgluino,theb,2,amuv,lamv,amuref) -
     .        3.D0/16.D0*alp**2*NS_lamb(amb/mgluino,asb2/mgluino)
     .        /mgluino*
     .        NS_gamvirt(asb2,asb1,amb,mgluino,theb,2,amuv,lamv) -
     .        3.D0/16.D0*asb2/mgluino*
     .        alp**2*NS_gamrealgl(asb2,amb,mgluino,theb,2,lamv) +
     .        alp/(4.D0*pi)*gsb2*
     .        NS_gamcfdec(ast2,ast1,amt,asb2,asb1,amb,mgluino,amsq,amuv,
     .        amuref)
      else
         qcdgsb2 = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c --- gluino --> supl + up

      if(asup1.le.mgluino) then
         gsupl=gs2/2.d0*(gur(1)**2+gul(1)**2)*
     .        (mgluino**2-asup1**2)*
     .        NS_lamb(0.d0,asup1/mgluino)/(16*pi*mgluino)
      else
         gsupl=0.d0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asup1.le.mgluino) then
         qcdgsupl = 0.D0
         alp   = gs2/(4.D0*pi)
         ca    = 3.D0
         cf    = 4.D0/3.D0
         amsq  = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         rval  = mgluino**2/amsq**2

         qcdgsupl = gsupl*alp/pi*( ca*(NS_gama(rval)-pi**2) + 
     .        cf*(NS_gamfcap(rval)+pi**2) + 4.D0*NS_gamf(rval) + 
     .        2.D0*pi**2/mgluino**2*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,mgluino,
     .        1,amsq) + 
     .       NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,amuref))
      else
         qcdgsupl = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c ---  gluino --> supr + up

      if(asup2.le.mgluino) then
         gsupr=gs2/2.d0*(gul(2)**2+gur(2)**2)*
     .        (mgluino**2-asup2**2)*
     .        NS_lamb(0.d0,asup2/mgluino)/(16*pi*mgluino)
      else
         gsupr=0.d0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asup2.le.mgluino) then
         qcdgsupr = 0.D0
         alp   = gs2/(4.D0*pi)
         ca    = 3.D0
         cf    = 4.D0/3.D0
         amsq  = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         rval  = mgluino**2/amsq**2

         qcdgsupr = gsupr*alp/pi*( ca*(NS_gama(rval)-pi**2) + 
     .        cf*(NS_gamfcap(rval)+pi**2) + 4.D0*NS_gamf(rval) + 
     .        2.D0*pi**2/mgluino**2*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,mgluino,
     .        2,amsq) + 
     .       NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,amuref))
      else
         qcdgsupr = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c --- gluino --> sdownl + down

      if(asdown1.le.mgluino) then
         gsdownl=gs2/2.d0*(gdl(1)**2+gdr(1)**2)*
     .        (mgluino**2-asdown1**2)*
     .        NS_lamb(0.d0,asdown1/mgluino)/(16*pi*mgluino)
      else
         gsdownl=0.d0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asdown1.le.mgluino) then
         qcdgsdownl = 0.D0
         alp   = gs2/(4.D0*pi)
         ca    = 3.D0
         cf    = 4.D0/3.D0
         amsq  = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         rval  = mgluino**2/amsq**2

         qcdgsdownl = gsdownl*alp/pi*( ca*(NS_gama(rval)-pi**2) + 
     .        cf*(NS_gamfcap(rval)+pi**2) + 4.D0*NS_gamf(rval) + 
     .        2.D0*pi**2/mgluino**2*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,mgluino,
     .        1,amsq) + 
     .       NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,amuref))
      else
         qcdgsdownl = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c ---  gluino --> sdownr + down

      if(asdown2.le.mgluino) then
         gsdownr=gs2/2.d0*(gdl(2)**2+gdr(2)**2)*
     .        (mgluino**2-asdown2**2)*
     .        NS_lamb(0.d0,asdown2/mgluino)/(16*pi*mgluino)
      else
         gsdownr=0.d0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1.D0) then
      if(asdown2.le.mgluino) then
         qcdgsdownr = 0.D0
         alp   = gs2/(4.D0*pi)
         ca    = 3.D0
         cf    = 4.D0/3.D0
         amsq  = 2.D0*(asup1+asup2+asdown1+asdown2)/8.D0
         rval  = mgluino**2/amsq**2

         qcdgsdownr = gsdownr*alp/pi*( ca*(NS_gama(rval)-pi**2) + 
     .        cf*(NS_gamfcap(rval)+pi**2) + 4.D0*NS_gamf(rval) + 
     .        2.D0*pi**2/mgluino**2*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,mgluino,
     .        2,amsq) + 
     .       NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,amuref))
      else
         qcdgsdownr = 0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c     ---- TOTAL WIDTH
c -------------------------------------------------------------------- c

      gluitot2lo = 2.D0*(gst1+gst2+gsb1+gsb2+2.D0*gsupl+2.D0*gsupr+
     .           2.D0*gsdownl+2.D0*gsdownr)

C  Resum if qcdcorr < -tree:
      If(qcdgst1.lt.-gst1) 
     .  write(*,*)"Warning: large negative rad. corrs. to go->st1+t"
      qcdgst1=resum(gst1,qcdgst1)
      If(qcdgst2.lt.-gst2) 
     .  write(*,*)"Warning: large negative rad. corrs. to go->st2+t"
      qcdgst2=resum(gst2,qcdgst2)
      If(qcdgsb1.lt.-gsb1) 
     .  write(*,*)"Warning: large negative rad. corrs. to go->sb1+b"
      qcdgsb1=resum(gsb1,qcdgsb1)
      If(qcdgsb2.lt.-gsb2) 
     .  write(*,*)"Warning: large negative rad. corrs. to go->sb2+t"
      qcdgsb2=resum(gsb2,qcdgsb2)
      If(qcdgsupl.lt.-gsupl) 
     .  write(*,*)"Warning: large negative rad. corrs. to go->supl+u"
      qcdgsupl=resum(gsupl,qcdgsupl)
      If(qcdgsupr.lt.-gsupr) 
     .  write(*,*)"Warning: large negative rad. corrs. to go->supr+u"
      qcdgsupr=resum(gsupr,qcdgsupr)
      If(qcdgsdownl.lt.-gsdownl) 
     .  write(*,*)"Warning: large negative rad. corrs. to go->sdownl+d"
      qcdgsdownl=resum(gsdownl,qcdgsdownl)
      If(qcdgsdownr.lt.-gsdownr) 
     .  write(*,*)"Warning: large negative rad. corrs. to go->sdownr+d"
      qcdgsdownr=resum(gsdownr,qcdgsdownr)
C  End resummation
      gluitot2nlo = 2.D0*(qcdgst1+qcdgst2+qcdgsb1+qcdgsb2+2.D0*qcdgsupl
     .                   +2.D0*qcdgsupr+2.D0*qcdgsdownl+2.D0*qcdgsdownr)
     .            + gluitot2lo

      gluitot2 = gluitot2nlo
      
c      write(*,*) "gluitot2,gluitot2lo:",gluitot2,gluitot2lo
c -------------------------------------------------------------------- c
c -- the 3-body decays and 3-body total widths --
c -------------------------------------------------------------------- c
         gluitotmulti=0.D0
      if(flagmulti.eq.1.D0) then
         call NS_xinteggo

         do i=1,5,1
            gluitotmulti=gluitotmulti+2.D0*xintegoup(i)+
     .           2.D0*xintegodn(i)+xintegotp(i)+xintegobt(i)
         end do
         do i=1,2,1
            gluitotmulti = gluitotmulti+4.D0*xintegoud(i)+
     .           2.D0*xintegotb(i)
         end do
         gluitotmulti = gluitotmulti+2.D0*xinteghcst1b+2.D0*xintegwst1b
        endif

c ---- Consider 3-body decays only if BR > multilim -------------------c
         if (gluitotmulti.lt.multilim*gluitot2)Then
           gluitotmulti =0.d0
         endif
c -------------------------------------------------------------------- c

c -- loop decays gluino -> neutralino_i + gluon
      if(flagloop.eq.1.D0) then
       if(gluitot2.eq.0.D0) then
         call NS_gluiraddecay(glnjgluon)
         gluitotrad=0.D0
         do i=1,5,1
            gluitotrad=gluitotrad+glnjgluon(i)
         end do
       endif
       endif
        
c ------------------------ the total width --------------------------- c
          gluitot = gluitot2+gluitotmulti+gluitotrad

c -------------------- the gluino branching ratios ------------------- c
c -- 2-body decays --

         gst1    = gst1+qcdgst1
         gst2    = gst2+qcdgst2
         gsb1    = gsb1+qcdgsb1
         gsb2    = gsb2+qcdgsb2
         gsupl   = gsupl+qcdgsupl
         gsupr   = gsupr+qcdgsupr
         gsdownl = gsdownl+qcdgsdownl
         gsdownr = gsdownr+qcdgsdownr

      brgst1    = gst1/gluitot
      brgst2    = gst2/gluitot
      brgsb1    = gsb1/gluitot
      brgsb2    = gsb2/gluitot
      brgsupl   = gsupl/gluitot
      brgsupr   = gsupr/gluitot
      brgsdownl = gsdownl/gluitot
      brgsdownr = gsdownr/gluitot

c -- 3-body and loop decays --
      if(gluitotmulti.ne.0d0)Then
         do i=1,5,1
            brgoup(i)=xintegoup(i)/gluitot
            brgoch(i)=xintegoup(i)/gluitot
            brgodn(i)=xintegodn(i)/gluitot
            brgost(i)=xintegodn(i)/gluitot
            brgotp(i)=xintegotp(i)/gluitot
            brgobt(i)=xintegobt(i)/gluitot
         end do
         do i=1,2,1
            brgoud(i)=xintegoud(i)/gluitot
            brgocs(i)=xintegoud(i)/gluitot
            brgotb(i)=xintegotb(i)/gluitot
         end do
         brhcst1b = xinteghcst1b/gluitot
         brwst1b  = xintegwst1b/gluitot
      endif     

      if(gluitot2.eq.0.D0) then
         do i=1,5,1
            brglnjgluon(i)=glnjgluon(i)/gluitot
         end do
      endif

      END

c ==================================================================== c
c             Radiative decays gluino -> neutralino_j gluon            c
c ==================================================================== c

       SUBROUTINE NS_gluiraddecay(glnjgluon)
*
       IMPLICIT NONE
       INTEGER J
       DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
       DOUBLE PRECISION glnjgluon(5),eps(5)
       DOUBLE PRECISION PI,SQR2
       DOUBLE PRECISION gabcd,gabcd0,gijgluon
       DOUBLE PRECISION NS_gluicoupabcd,NS_gluicoupabcd0
       DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
       COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
       COMMON/NS_pi/PI,SQR2
       COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*   
       EXTERNAL NS_gluicoupabcd,NS_gluicoupabcd0
*
       do j=1,5,1
         if(mgluino.gt.amneut(j)) then
            if(xmneut(j).ge.0.D0) then
               eps(j) = 1.D0
            else
               eps(j) = -1.D0
            endif
*            
            gabcd  = NS_gluicoupabcd(xmneut(j),j)
            gabcd0 = NS_gluicoupabcd0(xmneut(j),j) 
            gijgluon = -dsqrt(g2s)*g3s/32.D0/pi**2*eps(j)*mgluino*
     .           (gabcd+gabcd0)
            glnjgluon(j) = gijgluon**2*(mgluino**2-amneut(j)**2)**3/
     .                     8.D0/pi/mgluino**5*1.D0/4.D0
         else
            glnjgluon(j) = 0.D0
         endif
       enddo
      end

C    ---------------------------------------------
      DOUBLE PRECISION FUNCTION NS_gluicoupabcd(xmnj,j)
*
      IMPLICIT NONE
      INTEGER K,J
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     . gdr(2),gdl(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION mfer,mbos,gl,gr,fl,fr,glfrgrfl,glflgrfr,xmnj
     .,epsj,amnj
      COMPLEX*16 NS_iint,NS_i2int,NS_jint,NS_kint
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_pi/PI,SQR2 
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_neutstoptop/atopr,btopr
*
      EXTERNAL NS_iint,NS_i2int,NS_jint,NS_kint

      if(xmnj.le.0.D0) then
         epsj = -1.D0
      else
         epsj = 1.D0
      endif

      amnj = dabs(xmnj)
*
      NS_gluicoupabcd  = 0.D0
*
      do k=1,4,1
         if(k.eq.1) then
            mfer = amt
            mbos = ast1
            gl   = -dsqrt(2.D0)*atopr(1,j)
            gr   = -dsqrt(2.D0)*btopr(1,j)
            fl   = -2.D0*gtr(1)
            fr   = -2.D0*gtl(1)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.2) then
            mfer = amt
            mbos = ast2
            gl   = -dsqrt(2.D0)*atopr(2,j)
            gr   = -dsqrt(2.D0)*btopr(2,j)
            fl   = -2.D0*gtr(2)
            fr   = -2.D0*gtl(2)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.3) then
            mfer = amb
            mbos = asb1
            gl   = -dsqrt(2.D0)*abot(1,j)
            gr   = -dsqrt(2.D0)*bbot(1,j)
            fl   = -2.D0*gbr(1)
            fr   = -2.D0*gbl(1)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.4) then
            mfer = amb
            mbos = asb2
            gl   = -dsqrt(2.D0)*abot(2,j)
            gr   = -dsqrt(2.D0)*bbot(2,j)
            fl   = -2.D0*gbr(2)
            fr   = -2.D0*gbl(2)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         endif

         NS_gluicoupabcd = NS_gluicoupabcd + ( 
     .        glfrgrfl*(mgluino*dreal((NS_i2int(mgluino,amnj,mfer,mbos)-
     .                                 NS_kint(mgluino,amnj,mfer,mbos)))
     .               -epsj*amnj*dreal(NS_kint(mgluino,amnj,mfer,mbos)))
     .        +mfer*glflgrfr*dreal(NS_iint(mgluino,amnj,mfer,mbos)) )
      end do

      return

      end

c -------------------------------------------------------------------- c
      DOUBLE PRECISION FUNCTION NS_gluicoupabcd0(xmnj,j)
*
      IMPLICIT NONE
      INTEGER j,k
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     . gdr(2),gdl(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION mfer,mbos,gl,gr,fl,fr,glfrgrfl,glflgrfr,xmnj
     .,epsj,amnj,amu,amd
      COMPLEX*16 NS_i2int0,NS_kint0,NS_iint
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2 
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
*
      EXTERNAL NS_i2int0,NS_kint0,NS_iint,NS_jint0

      if(xmnj.le.0.D0) then
         epsj = -1.D0
      else
         epsj = 1.D0
      endif

      amnj = dabs(xmnj)

      amu = 1.D-2
      amd = 1.D-2

      NS_gluicoupabcd0  = 0.D0

      do k=1,4,1
         if(k.eq.1) then
            mfer = amu
            mbos = asup1
            gl   = -dsqrt(2.D0)*aup(1,j)
            gr   = -dsqrt(2.D0)*bup(1,j)
            fl   = -2.D0*gur(1)
            fr   = -2.D0*gul(1)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.2) then
            mfer = amu
            mbos = asup2
            gl   = -dsqrt(2.D0)*aup(2,j)
            gr   = -dsqrt(2.D0)*bup(2,j)
            fl   = -2.D0*gur(2)
            fr   = -2.D0*gul(2)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.3) then
            mfer = amd
            mbos = asdown1
            gl   = -dsqrt(2.D0)*ado(1,j)
            gr   = -dsqrt(2.D0)*bdo(1,j)
            fl   = -2.D0*gdr(1)
            fr   = -2.D0*gdl(1)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.4) then
            mfer = amd
            mbos = asdown2
            gl   = -dsqrt(2.D0)*ado(2,j)
            gr   = -dsqrt(2.D0)*bdo(2,j)
            fl   = -2.D0*gdr(2)
            fr   = -2.D0*gdl(2)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         endif

         NS_gluicoupabcd0 = NS_gluicoupabcd0 + ( 
     .        glfrgrfl*(mgluino*dreal((NS_i2int0(mgluino,amnj,mbos)-
     .                                 NS_kint0(mgluino,amnj,mbos)))
     .                 -epsj*amnj*dreal(NS_kint0(mgluino,amnj,mbos)) ) )
      end do

      NS_gluicoupabcd0 = 2.D0*NS_gluicoupabcd0

      return
      end
c ==================================================================== c
c                         Gluino 3-body decays
c ==================================================================== c
      SUBROUTINE NS_xinteggo
*
      IMPLICIT NONE
      INTEGER nx1t,ny1t,ni,nj
      DOUBLE PRECISION xintegoup(5),xintegodn(5),xintegotp(5),
     .        xintegobt(5),xintegoud(2),xintegotb(2),xintegocc(2),
     .        xinteghcst1b,xintegwst1b
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),amch
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION SUM,xmu1,xmu2,xmu3
      DOUBLE PRECISION NS_goup,NS_godn,NS_gobt,NS_gotp,NS_goud,NS_gotb,
     .         NS_ghcst1b,NS_gwst1b,NS_gocc
      DOUBLE PRECISION NS_ay,NS_by,NS_ax,NS_bx
C
      COMMON/NS_nx1/nx1t,ny1t
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,amch
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_indices/ni,nj
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_pi/PI,SQR2
      COMMON/GLUINO_3GAMMA/xintegoup,xintegodn,xintegotp,
     .         xintegobt,xintegoud,xintegotb,xintegocc,
     .         xinteghcst1b,xintegwst1b
*
      EXTERNAL NS_goup,NS_godn,NS_gobt,NS_gotp,NS_goud,NS_gotb,
     .         NS_ghcst1b,NS_gwst1b,NS_gocc
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
c -------------------------------------------------------------------- c
c                    gluino into H+ + stop1 + b
c -------------------------------------------------------------------- c
      xmu1=amb**2/mgluino**2
      xmu2=ast1**2/mgluino**2
      xmu3=amch**2/mgluino**2

      if(mgluino.gt.(ast1+amch+amb)) then
         call NS_integ2(NS_ghcst1b,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum)
         xinteghcst1b=sum*mgluino/128.D0/(2.D0*pi)**3
      else
         xinteghcst1b=0.D0
      endif
*      write(*,*)'gluino into H+ + stop1 + b',xinteghcst1b

c -------------------------------------------------------------------- c
c                    gluino into W+ + stop1 + b
c -------------------------------------------------------------------- c
      xmu1=amb**2/mgluino**2
      xmu2=ast1**2/mgluino**2
      xmu3=mw**2/mgluino**2

      if(mgluino.gt.(ast1+mw+amb)) then
         call NS_integ2(NS_gwst1b,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum)
         xintegwst1b=sum*mgluino/128.D0/(2.D0*pi)**3
      else
         xintegwst1b=0.D0
      endif
*      write(*,*)'gluino into W+ + stop1 + b',xintegwst1b
c -------------------------------------------------------------------- c
c               gluino --> neutralino_j + up + upbar
c -------------------------------------------------------------------- c
      do nj=1,5,1
         xmu1=0.D0
         xmu2=0.D0
         xmu3=amneut(nj)**2/mgluino**2

         if(mgluino.gt.amneut(nj)) then
            call NS_integ2(NS_goup,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .                  xmu3,nx1t,ny1t,sum)
            xintegoup(nj)=sum*mgluino/(2*pi)**3/8.d0/64.D0
         else 
            xintegoup(nj)=0.D0
         endif
      end do
c -------------------------------------------------------------------- c
c             gluino --> neutralino_j + down + downbar
c -------------------------------------------------------------------- c
      do nj=1,5,1
         xmu1=0.D0
         xmu2=0.D0
         xmu3=amneut(nj)**2/mgluino**2

         if(mgluino.gt.amneut(nj)) then
            call NS_integ2(NS_godn,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .                  xmu3,nx1t,ny1t,sum)
            xintegodn(nj)=sum*mgluino/(2*pi)**3/8.D0/64.D0
         else
            xintegodn(nj)=0.D0
         endif
      end do
c -------------------------------------------------------------------- c
c            gluino --> neutralino_j + bottom + bottombar
c -------------------------------------------------------------------- c
      do nj=1,5,1
         xmu1=amb**2/mgluino**2
         xmu2=amb**2/mgluino**2
         xmu3=amneut(nj)**2/mgluino**2

         if(mgluino.gt.(amneut(nj)+2.D0*amb)) then
            call NS_integ2(NS_gobt,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .                  xmu3,nx1t,ny1t,sum)
            xintegobt(nj)=sum*mgluino/(2*pi)**3/8.d0/64.D0
         else
            xintegobt(nj)=0.D0
         endif
      end do
c -------------------------------------------------------------------- c
c              gluino --> neutralino_j + top topbar
c -------------------------------------------------------------------- c
      do nj=1,5,1
         xmu1=amt**2/mgluino**2
         xmu2=amt**2/mgluino**2
         xmu3=amneut(nj)**2/mgluino**2

         if(mgluino.gt.(amneut(nj)+2.D0*amt)) then
            call NS_integ2(NS_gotp,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .                  xmu3,nx1t,ny1t,sum)
            xintegotp(nj)=sum*mgluino/(2*pi)**3/8.d0/64.D0
         else
            xintegotp(nj)=0.D0
         endif
      end do
c -------------------------------------------------------------------- c
c               gluino --> chargino_j- + up + downbar
c -------------------------------------------------------------------- c
      do nj=1,2,1
         xmu1=0.D0
         xmu2=0.D0
         xmu3=amchar(nj)**2/mgluino**2

         if(mgluino.gt.amchar(nj)) then
            call NS_integ2(NS_goud,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .                  xmu3,nx1t,ny1t,sum)
            xintegoud(nj)=sum*mgluino/(2*pi)**3/8.d0/64.D0
         else
            xintegoud(nj)=0.D0
         endif
      end do
c -------------------------------------------------------------------- c
c              gluino --> chargino_j- + top + bottombar
c -------------------------------------------------------------------- c
      do nj=1,2,1
         xmu1=amt**2/mgluino**2
         xmu2=amb**2/mgluino**2
         xmu3=amchar(nj)**2/mgluino**2

         if(mgluino.gt.(amchar(nj)+amt+amb)) then
            call NS_integ2(NS_gotb,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .                  xmu3,nx1t,ny1t,sum)
            xintegotb(nj)=sum*mgluino/(2*pi)**3/8.d0/64.D0
         else
            xintegotb(nj)=0.D0
         endif
      end do
c -------------------------------------------------------------------- c
c              gluino --> chargino_j- + top + bottombar
c -------------------------------------------------------------------- c
      do nj=1,2,1
         xmu1=amb**2/mgluino**2
         xmu2=amt**2/mgluino**2
         xmu3=amchar(nj)**2/mgluino**2

         if(mgluino.gt.(amchar(nj)+amt+amb)) then
            call NS_integ2(NS_gocc,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .                  xmu3,nx1t,ny1t,sum)
            xintegocc(nj)=sum*mgluino/(2*pi)**3/8.d0/64.D0
         else
            xintegocc(nj)=0.D0
         endif
      end do
      end 
c ==================================================================== c
c                    gluino --> H+ + stop1 + bottom
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_ghcst1b(x1,x2)
      IMPLICIT NONE
      INTEGER I,K
      DOUBLE PRECISION dsb(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     .gdr(2),gdl(2)
      DOUBLE PRECISION gctbr(2,2)
      DOUBLE PRECISION achtop,vchtop,achtau,vchtau
      DOUBLE PRECISION chtbrunr,chtbrunl
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),amch
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION xmust1,xmub,xmut,xmuch,dt,ac,vc
     .,ghcsbot,ghctop,ghcsbtop,X1,X2,xmusb(2)
C
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,amch
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup15/achtop,vchtop,achtau,vchtau
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
      COMMON/NS_hcsbotstop/gctbr
*
      xmusb(1) = asb1**2/mgluino**2
      xmusb(2) = asb2**2/mgluino**2
      xmust1 = ast1**2/mgluino**2
      xmub   = amb**2/mgluino**2
      xmut   = amt**2/mgluino**2
      xmuch  = amch**2/mgluino**2
*      
      dsb(1) = 1-x1+xmub-xmusb(1)
      dsb(2) = 1-x1+xmub-xmusb(2)
      dt     = 1-x2+xmust1-xmut
*
      ac = -vchtop+achtop
      vc = -vchtop-achtop
c -------------------------------------------------------------------- c
c                             top exchange
c -------------------------------------------------------------------- c
      ghctop=0.D0
      if ((amt+ast1).gt.mgluino)Then 
      if(mgluino.gt.(amch+ast1+amb)) then
         ghctop=2.D0*g2s*g3s/dt**2*
     .        ( 2.D0*dsqrt(xmut*xmub)*
     .        (gtr(1)**2+gtl(1)**2)*ac*vc*(2.D0-x2)+4.D0*dsqrt(xmub)*
     .        gtl(1)*gtr(1)*vc*ac*(1.D0-x2+xmust1+xmut)+2.D0*
     .        dsqrt(xmut)*gtr(1)*gtl(1)*(ac**2+vc**2)*(-xmuch-x2+
     .        xmub+xmust1+1.D0)+(gtl(1)**2*ac**2+gtr(1)**2*vc**2)*(
     .        x1*x2-xmust1*x1-x1+x2**2+x2*xmuch-x2*xmub-x2*xmust1
     .        -3.D0*x2-2.D0*xmuch+2.D0*xmub+2.D0*xmust1+2.D0)+
     .        (gtl(1)**2*vc**2+gtr(1)**2*ac**2)*xmut*x1 )
      else 
         ghctop=0.D0
      endif
      endif
c -------------------------------------------------------------------- c
c                           sbottom exchange
c -------------------------------------------------------------------- c
      ghcsbot = 0.D0
       if(mgluino.gt.(amch+ast1+amb)) then
         do i=1,2
            do k=1,2
               if (xmusb(i).gt.1.d0.and.xmusb(k).gt.1.d0)Then
               ghcsbot=ghcsbot+2.D0*g3s*g2s/dsb(i)/dsb(k)*
     .              gctbr(1,i)*gctbr(1,k)*amw**2/mgluino**2*(
     .              (gbr(i)*gbr(k)+gbl(i)*gbl(k))*x1+
     .              2.D0*dsqrt(xmub)*(gbr(i)*gbl(k)+gbr(k)*gbl(i)) )
               endif
            enddo
         enddo
      else 
         ghcsbot=0.D0
      endif
c -------------------------------------------------------------------- c
c                       sbottom top interference
c -------------------------------------------------------------------- c
      ghcsbtop = 0.D0
     
      if(mgluino.gt.(amch+ast1+amb)) then
         do i=1,2
            If(xmusb(i).gt.0d0.and.(amt+ast1).gt.mgluino)Then
            ghcsbtop=ghcsbtop
     .           -2.D0*g3s*g2s*gctbr(1,i)*amw*2.D0/mgluino/dt/dsb(i)*(
     .           dsqrt(xmut)*x1*(gtr(1)*gbr(i)*ac+gtl(1)*gbl(i)*vc)
     .           +2.D0*dsqrt(xmub*xmut)*(gtr(1)*gbl(i)*ac+
     .           gbr(i)*gtl(1)*vc)+dsqrt(xmub)*(gtr(1)*gbr(i)*vc+
     .           gtl(1)*gbl(i)*ac)*(2.D0-x2)+(gtr(1)*gbl(i)*vc+
     .           gbr(i)*gtl(1)*ac)*(-xmuch-x2+xmub+xmust1+1.D0))
            endif
         enddo
      else 
         ghcsbtop=0.D0
      endif
      NS_ghcst1b=ghctop+ghcsbot+ghcsbtop
      end

c ==================================================================== c
c                    gluino --> W+ + stop1 + bottombar
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_gwst1b(x1,x2)
*
      IMPLICIT NONE
      INTEGER I,K
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION dsb(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),
     .gul(2),gdr(2),gdl(2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2)
      DOUBLE PRECISION vtb,xmusb(2), xmust1,xmub,xmut,xmuw,dt
     .,gwtop,gwsbot,gwsbtop,X1,X2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*      
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_coup20/gwtb,gwntau
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
*
      vtb = 0.99915D0
*
      xmusb(1) = asb1**2/mgluino**2
      xmusb(2) = asb2**2/mgluino**2
      xmust1 = ast1**2/mgluino**2
      xmub   = amb**2/mgluino**2
      xmut   = amt**2/mgluino**2
      xmuw   = mw**2/mgluino**2
*      
      dsb(1) = 1-x1+xmub-xmusb(1)
      dsb(2) = 1-x1+xmub-xmusb(2)
      dt     = 1-x2+xmust1-xmut
c -------------------------------------------------------------------- c
c                              top exchange
c -------------------------------------------------------------------- c
      gwtop=0.D0
      if(mgluino.gt.(mw+ast1+amb)) then
         if ((amt+ast1).gt.mgluino)Then
         gwtop=2.D0*g2s*g3s/dt**2*vtb**2/2.D0*( 
     .        2.D0*gtr(1)*gtl(1)*dsqrt(xmut)*
     .        ((1.D0-x2+xmub+xmust1-2.D0*xmuw)+1.D0/xmuw*
     .        (1.D0-x2+xmust1-xmub)**2) 
     .        +gtr(1)**2*(x1*xmut+
     .         xmut/xmuw*(1.D0-x2+xmust1-xmub-xmuw)*(2.D0-x1-x2))
     .        +gtl(1)**2*
     .        ( (1.D0-xmust1)*x1+(-2.D0+x2)*(-1.D0+x1+x2+
     .        xmuw-xmust1-xmub)+1.D0/xmuw*(1.D0-x2+xmust1-xmub-xmuw)*
     .        ((-2.D0+x2)*(1.D0-x1+xmub-xmust1-xmuw)+(1.D0-xmust1)*(
     .        2.D0-x1-x2))) )
         endif
      else 
         gwtop=0.D0
      endif
c -------------------------------------------------------------------- c
c                           sbottom exchange
c -------------------------------------------------------------------- c
      gwsbot = 0.D0
      
      if(mgluino.gt.(mw+ast1+amb)) then
         do i=1,2
            do k=1,2
               if(xmusb(i).gt.1.d0.and.xmusb(k).gt.1.d0)Then
               gwsbot=gwsbot+g3s*g2s/2.D0*gwtb(1,i)*gwtb(1,k)/dsb(i)
     .              /dsb(k)*2.D0*
     .              (2.D0*dsqrt(xmub)*(gbr(i)*gbl(k)+
     .              gbr(k)*gbl(i))+x1*(gbr(i)*gbr(k)+gbl(i)*gbl(k)))*
     .              (-(2.D0*(1.D0-x1+xmub)+2.D0*xmust1-xmuw)+
     .              1.D0/xmuw*(1.D0-x1+xmub-xmust1)**2)
               endif
            enddo
         enddo
      else 
         gwsbot=0.D0
      endif
  
c -------------------------------------------------------------------- c
c                       sbottom top interference
c -------------------------------------------------------------------- c
      gwsbtop = 0.D0
            if(mgluino.gt.(mw+ast1+amb)) then
         do i=1,2
            if(xmusb(i).gt.1.d0.and.(amt+ast1).gt.mgluino)Then
            gwsbtop=gwsbtop+g3s*g2s/2.D0*gwtb(1,i)*vtb*2.D0*2.D0*(
     .           gtr(1)*gbr(i)*dsqrt(xmub*xmut)*(2.D0*x2-1.D0/xmuw*(
     .           1.D0-x1+xmub-xmust1-xmuw)*(2.D0-x1-x2))+gbr(i)*
     .           gtl(1)*dsqrt(xmub)*(2.D0*x2-4.D0*xmust1-1.D0/xmuw*(
     .           1.D0-x2+xmust1-xmub+xmuw)*(1.D0-x1+xmub-xmust1
     .           -xmuw))+gtr(1)*gbl(i)*dsqrt(xmut)*(2.D0*x1-4.D0*xmub
     .           -1.D0/xmuw*(1.D0-x2+xmust1-xmub-xmuw)*(1.D0-x1+xmub
     .           -xmust1+xmuw))+gtl(1)*gbl(i)*(-2.D0-x1*x2-x1*xmust1
     .           +3.D0*x1+x2*xmub+2.D0*x2-4.D0*xmub-2.D0*xmust1+2.D0*
     .           xmuw+1.D0/xmuw*(1.D0-x1+xmub-xmust1)*(2.D0*xmub
     .           -x2*xmub+x1*x2-xmust1*x1-x1)) )
            endif
         enddo
      else 
         gwsbtop=0.D0
      endif
     
      NS_gwst1b=gwtop+gwsbot+gwsbtop
*      write(*,*)'NS_gwst1b',gwtop,gwsbot,gwsbtop
      end
c ==================================================================== c
c                    gluino --> neutralino up upbar
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_goup(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,I,K
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsup(2),dsupb(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     .gdr(2),gdl(2)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION xmusup(2),x3,y1,y2,y3,X1,X2,xmuneut1
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
*
      xmuneut1 = amneut(nj)**2/mgluino**2

      xmusup(1) = asup1**2/mgluino**2
      xmusup(2) = asup2**2/mgluino**2

      dsup(1)  = 1.D0-x1-xmusup(1)
      dsup(2)  = 1.D0-x1-xmusup(2)
      dsupb(1) = 1.D0-x2-xmusup(1)
      dsupb(2) = 1.D0-x2-xmusup(2)

      x3 = 2.D0-x1-x2
      y1 = 1.D0-xmuneut1-x1
      y2 = 1.D0-xmuneut1-x2
      y3 = 1.D0+xmuneut1-x3

      NS_goup=0.D0
      
      if (mgluino.gt.amneut(nj)) then
         do i=1,2
            do k=1,2
             if (xmusup(i).gt.1.d0.and.xmusup(k).gt.1.d0)Then  
               NS_goup=NS_goup+4.D0*g3s*g2s/dsup(i)/dsup(k)*2.D0*(
     .              (gur(i)*gur(k)+gul(i)*gul(k))*
     .              (aup(i,nj)*aup(k,nj)+bup(k,nj)*bup(i,nj))*
     .              x1*y1)
     .         +4.D0*g3s*g2s/dsupb(i)/dsupb(k)*2.D0*(
     .              (gur(i)*gur(k)+gul(i)*gul(k))*
     .              (aup(i,nj)*aup(k,nj)+bup(k,nj)*bup(i,nj))*
     .              x2*y2)
     .         +4.D0*g3s*g2s/dsup(i)/dsupb(k)*2.D0*(
     .         (aup(k,nj)*gur(k)*gul(i)*bup(i,nj)+gul(k)*
     .          aup(i,nj)*bup(k,nj)*gur(i))*(-x1*y1-x2*y2+x3*y3)
     .         +2.D0*(gul(k)*aup(k,nj)*gul(i)*aup(i,nj)+
     .                gur(k)*bup(k,nj)*gur(i)*bup(i,nj))*
     .          xmneut(nj)/mgluino*y3)
               endif
            enddo
         enddo
      else 
         NS_goup=0.d0
      endif
*	WRITE(*,*)'NS_goup',NS_goup
      end
c ==================================================================== c
c                   gluino --> neutralino down downbar
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_godn(x1,x2)
*	
      IMPLICIT NONE
      INTEGER ni,nj,I,K
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsdn(2),dsdnB(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     .gdr(2),gdl(2)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION xmuneut1,xmusd(2),x3,y1,y2,y3,X1,X2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/NS_indices/ni,nj
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu

      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
*     
      xmuneut1 = amneut(nj)**2/mgluino**2

      xmusd(1) = asdown1**2/mgluino**2
      xmusd(2) = asdown2**2/mgluino**2

      dsdn(1)  = 1-x1-xmusd(1)
      dsdn(2)  = 1-x1-xmusd(2)
      dsdnb(1) = 1-x2-xmusd(1)
      dsdnb(2) = 1-x2-xmusd(2)

      x3 = 2-x1-x2
      y1 = 1-xmuneut1-x1
      y2 = 1-xmuneut1-x2
      y3 = 1+xmuneut1-x3
    
      NS_godn=0.D0
      if (mgluino.gt.amneut(nj)) then
          do i=1,2
            do k=1,2
               if (xmusd(i).gt.1.d0.and.xmusd(k).gt.1.d0)Then
               NS_godn=NS_godn+4.D0*g3s*g2s/dsdn(i)/dsdn(k)*2.D0*(
     .              (gdr(i)*gdr(k)+gdl(i)*gdl(k))*
     .              (ado(i,nj)*ado(k,nj)+bdo(k,nj)*bdo(i,nj))*
     .              x1*y1)
     .         +4.D0*g3s*g2s/dsdnb(i)/dsdnb(k)*2.D0*(
     .              (gdr(i)*gdr(k)+gdl(i)*gdl(k))*
     .              (ado(i,nj)*ado(k,nj)+bdo(k,nj)*bdo(i,nj))*
     .              x2*y2)
     .         +4.D0*g3s*g2s/dsdn(i)/dsdnb(k)*2.D0*(
     .         (ado(k,nj)*gdr(k)*gdl(i)*bdo(i,nj)+gdr(i)*
     .          ado(i,nj)*bdo(k,nj)*gdl(k))*(-x1*y1-x2*y2+x3*y3)
     .         +2.D0*(gdl(k)*ado(k,nj)*gdl(i)*ado(i,nj)+
     .                gdr(k)*bdo(k,nj)*gdr(i)*bdo(i,nj))*
     .          xmneut(nj)/mgluino*y3)
               endif
            enddo
         enddo
      else 
         NS_godn=0.d0
      endif


      end
c ==================================================================== c
c                 gluino --> neutralino bottom bottombar
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_gobt(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,I,K,J
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsb(2),dsbb(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     .gdr(2),gdl(2)
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION abot1(2,5),bbot1(2,5),abot2(2,5),bbot2(2,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION xmuneut1,xmusb(2),X1,X2,uh,th,xmub,db11,
     .db12,db21,db22,db3,db4,ab11,ab12,ab13,ab14,ab21,ab22,ab23,ab24,
     .ab5,ab6,ab7,ab8,ab9,ab10
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_neutsbotbot/abot,bbot
*
      do i=1,2,1
         do j=1,5,1
            abot1(i,j)=abot(i,j)
            bbot1(i,j)=bbot(i,j)
            abot2(i,j)=abot(i,j)
            bbot2(i,j)=bbot(i,j)
         end do
      end do

      xmuneut1 = amneut(nj)**2/mgluino**2
     
      xmusb(1) = asb1**2/mgluino**2
      xmusb(2) = asb2**2/mgluino**2

      dsb(1)  = 1-x1-xmusb(1)+amb**2/mgluino**2
      dsb(2)  = 1-x1-xmusb(2)+amb**2/mgluino**2
      dsbb(1) = 1-x2-xmusb(1)+amb**2/mgluino**2
      dsbb(2) = 1-x2-xmusb(2)+amb**2/mgluino**2

      NS_gobt=0.D0

      uh = 1.D0-x1+amb**2/mgluino**2
      th = 1.D0-x2+amb**2/mgluino**2
      xmub = amb**2/mgluino**2
            if(mgluino.gt.(amneut(nj)+2.D0*amb)) then 
         do i=1,2
            do k=1,2
               if (xmusb(i).gt.1.d0.and.xmusb(k).gt.1.d0)Then
               db11 = abot1(k,nj)*bbot1(i,nj)+abot1(i,nj)*bbot1(k,nj)
               db12 = abot1(k,nj)*abot1(i,nj)+bbot1(k,nj)*bbot1(i,nj)
               db21 = abot2(k,nj)*bbot2(i,nj)+abot2(i,nj)*bbot2(k,nj)
               db22 = abot2(k,nj)*abot2(i,nj)+bbot2(k,nj)*bbot2(i,nj)
               db3  = gbl(k)*gbr(i)+gbl(i)*gbr(k)
               db4  = gbl(i)*gbl(k)+gbr(i)*gbr(k)

               ab11 = db11*db3
               ab12 = db11*db4
               ab13 = db12*db3
               ab14 = db12*db4

               ab21 = db21*db3
               ab22 = db21*db4
               ab23 = db22*db3
               ab24 = db22*db4

               ab5 = abot1(k,nj)*bbot2(i,nj)*gbl(k)*gbr(i)+
     .               abot2(i,nj)*bbot1(k,nj)*gbl(i)*gbr(k)
               ab6 = abot2(i,nj)*bbot1(k,nj)*gbl(k)*gbl(i)+
     .               bbot2(i,nj)*abot1(k,nj)*gbr(k)*gbr(i)
               ab7 = bbot1(k,nj)*bbot2(i,nj)*gbl(i)*gbr(k)+
     .               abot1(k,nj)*abot2(i,nj)*gbr(i)*gbl(k)
               ab8 = abot1(k,nj)*abot2(i,nj)*gbr(k)*gbr(i)+
     .               bbot1(k,nj)*bbot2(i,nj)*gbl(i)*gbl(k)
               ab9 = abot2(i,nj)*bbot1(k,nj)*gbr(k)*gbr(i)+
     .               abot1(k,nj)*bbot2(i,nj)*gbl(i)*gbl(k)
               ab10 = abot1(k,nj)*bbot2(i,nj)*gbl(i)*gbr(k)+
     .                abot2(i,nj)*bbot1(k,nj)*gbl(k)*gbr(i)
               ab11 = abot1(k,nj)*abot2(i,nj)*gbl(k)*gbl(i)+
     .                bbot1(k,nj)*bbot2(i,nj)*gbr(i)*gbr(k)
               ab12 = bbot1(k,nj)*bbot2(i,nj)*gbl(k)*gbr(i)+
     .                abot1(k,nj)*abot2(i,nj)*gbr(k)*gbl(i)

               NS_gobt=NS_gobt
     .           +8.D0*g3s*g2s/dsb(i)/dsb(k)*(
     .           -4.D0*ab11*xmneut(nj)/mgluino*xmub
     .           +2.D0*ab12*xmneut(nj)*amb/mgluino**2*(-xmub-1.D0+uh)
     .           +2.D0*ab13*amb/mgluino*(-xmub-xmuneut1+uh)
     .           +ab14*(-uh**2+uh*(1.D0+xmuneut1+2.D0*xmub)-
     .                 (xmuneut1+xmub)*(1.D0+xmub)) )
     .         +8.D0*g3s*g2s/dsbb(i)/dsbb(k)*(
     .           -4.D0*ab21*xmneut(nj)/mgluino*xmub
     .           +2.D0*ab22*xmneut(nj)*amb/mgluino**2*(-xmub-1.D0+th)
     .           +2.D0*ab23*amb/mgluino*(-xmub-xmuneut1+th)
     .           +ab24*(-th**2+th*(1.D0+xmuneut1+2.D0*xmub)-
     .                 (xmuneut1+xmub)*(1.D0+xmub)) )
     .         -2.D0*8.D0*g3s*g2s/dsb(k)/dsbb(i)*(
     .           ab5*xmub*(uh+th-2.D0*xmub) 
     .           +ab6*amb/mgluino*(th-xmuneut1-xmub) 
     .           +ab7*xmneut(nj)*amb/mgluino**2*(uh-xmub-1.D0) 
     .           -2.D0*ab8*xmneut(nj)/mgluino*xmub
     .           +ab9*amb/mgluino*(uh-xmuneut1-xmub)
     .           +ab10*(uh*th-xmuneut1-xmub**2)
     .           +ab11*xmneut(nj)/mgluino*(uh+th-xmuneut1-1.D0)
     .           +ab12*xmneut(nj)*amb/mgluino**2*(th-xmub-1.D0) )
               endif
            enddo
         enddo
      else 
         NS_gobt=0.D0
      endif

*      WRITE(*,*)'NS_gobt',NS_gobt
      end
c ==================================================================== c
c                   gluino --> neutralino top topbar
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_gotp(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,I,j,K
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION xmuneut1,xmut,xmust(2),X1,X2,uh,th,
     .db11,db12,db21,db22,db3,db4,ab11,ab12,ab13,ab14,ab21,ab22,ab23,
     .ab24,ab5,ab6,ab7,ab8,ab9,ab10,gmst(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dst(2),dstb(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     .gdr(2),gdl(2)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION atopr1(2,5),btopr1(2,5),atopr2(2,5),btopr2(2,5)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_neutstoptop/atopr,btopr
*
      do i=1,2,1
         do j=1,5,1
            atopr1(i,j)=atopr(i,j)
            btopr1(i,j)=btopr(i,j)
            atopr2(i,j)=atopr(i,j)
            btopr2(i,j)=btopr(i,j)
         end do
      end do
  
      xmuneut1 = amneut(nj)**2/mgluino**2
      
      gmst(1) = ast1
      gmst(2) = ast2
      xmust(1) = ast1**2/mgluino**2
      xmust(2) = ast2**2/mgluino**2

      dst(1)  = 1-x1-xmust(1)+amt**2/mgluino**2
      dst(2)  = 1-x1-xmust(2)+amt**2/mgluino**2
      dstb(1) = 1-x2-xmust(1)+amt**2/mgluino**2
      dstb(2) = 1-x2-xmust(2)+amt**2/mgluino**2

      NS_gotp=0.D0

      uh = 1.D0-x1+amt**2/mgluino**2
      th = 1.D0-x2+amt**2/mgluino**2
      xmut = amt**2/mgluino**2
            if(mgluino.gt.(amneut(nj)+2.D0*amt)) then 
         do i=1,2
            do k=1,2
       if ((gmst(i)+amt).gt.mgluino.and.(gmst(k)+amt).gt.mgluino)Then
              db11 = atopr1(k,nj)*btopr1(i,nj)+atopr1(i,nj)*btopr1(k,nj)
              db12 = atopr1(k,nj)*atopr1(i,nj)+btopr1(k,nj)*btopr1(i,nj)
              db21 = atopr2(k,nj)*btopr2(i,nj)+atopr2(i,nj)*btopr2(k,nj)
              db22 = atopr2(k,nj)*atopr2(i,nj)+btopr2(k,nj)*btopr2(i,nj)
               db3 = gtl(k)*gtr(i)+gtl(i)*gtr(k)
               db4 = gtl(i)*gtl(k)+gtr(i)*gtr(k)

               ab11 = db11*db3
               ab12 = db11*db4
               ab13 = db12*db3
               ab14 = db12*db4

               ab21 = db21*db3
               ab22 = db21*db4
               ab23 = db22*db3
               ab24 = db22*db4

               ab5 = atopr1(k,nj)*btopr2(i,nj)*gtl(k)*gtr(i)+
     .               atopr2(i,nj)*btopr1(k,nj)*gtl(i)*gtr(k)
               ab6 = atopr2(i,nj)*btopr1(k,nj)*gtl(k)*gtl(i)+
     .               btopr2(i,nj)*atopr1(k,nj)*gtr(k)*gtr(i)
               ab7 = btopr1(k,nj)*btopr2(i,nj)*gtl(i)*gtr(k)+
     .               atopr1(k,nj)*atopr2(i,nj)*gtr(i)*gtl(k)
               ab8 = atopr1(k,nj)*atopr2(i,nj)*gtr(k)*gtr(i)+
     .               btopr1(k,nj)*btopr2(i,nj)*gtl(i)*gtl(k)
               ab9 = atopr2(i,nj)*btopr1(k,nj)*gtr(k)*gtr(i)+
     .               atopr1(k,nj)*btopr2(i,nj)*gtl(i)*gtl(k)
               ab10 = atopr1(k,nj)*btopr2(i,nj)*gtl(i)*gtr(k)+
     .                atopr2(i,nj)*btopr1(k,nj)*gtl(k)*gtr(i)
               ab11 = atopr1(k,nj)*atopr2(i,nj)*gtl(k)*gtl(i)+
     .                btopr1(k,nj)*btopr2(i,nj)*gtr(i)*gtr(k)
               ab12 = btopr1(k,nj)*btopr2(i,nj)*gtl(k)*gtr(i)+
     .                atopr1(k,nj)*atopr2(i,nj)*gtr(k)*gtl(i)

               NS_gotp=NS_gotp
     .           +8.D0*g3s*g2s/dst(i)/dst(k)*(
     .           -4.D0*ab11*xmneut(nj)/mgluino*xmut
     .           +2.D0*ab12*xmneut(nj)*amt/mgluino**2*(-xmut-1.D0+uh)
     .           +2.D0*ab13*amt/mgluino*(-xmut-xmuneut1+uh)
     .           +ab14*(-uh**2+uh*(1.D0+xmuneut1+2.D0*xmut)-
     .                 (xmuneut1+xmut)*(1.D0+xmut)) )
     .         +8.D0*g3s*g2s/dstb(i)/dstb(k)*(
     .           -4.D0*ab21*xmneut(nj)/mgluino*xmut
     .           +2.D0*ab22*xmneut(nj)*amt/mgluino**2*(-xmut-1.D0+th)
     .           +2.D0*ab23*amt/mgluino*(-xmut-xmuneut1+th)
     .           +ab24*(-th**2+th*(1.D0+xmuneut1+2.D0*xmut)-
     .                 (xmuneut1+xmut)*(1.D0+xmut)) )
     .         -2.D0*8.D0*g3s*g2s/dst(k)/dstb(i)*(
     .           ab5*xmut*(uh+th-2.D0*xmut) 
     .           +ab6*amt/mgluino*(th-xmuneut1-xmut) 
     .           +ab7*xmneut(nj)*amt/mgluino**2*(uh-xmut-1.D0) 
     .           -2.D0*ab8*xmneut(nj)/mgluino*xmut
     .           +ab9*amt/mgluino*(uh-xmuneut1-xmut)
     .           +ab10*(uh*th-xmuneut1-xmut**2)
     .           +ab11*xmneut(nj)/mgluino*(uh+th-xmuneut1-1.D0)
     .           +ab12*xmneut(nj)*amt/mgluino**2*(th-xmut-1.D0) )
       endif
            enddo
         enddo
      else 
         NS_gotp=0.D0
      endif

*      write(*,*)'NS_gotp',NS_gotp
      end
c ==================================================================== c
c                   gluino --> chargino- up downbar
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_goud(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsup(2),dsdn(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     .gdr(2),gdl(2)
      DOUBLE PRECISION alup(2,2),aldo(2,2),blup(2,2),bldo(2,2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION x1,x2,
     .xmusdn(2),xmusup(2),xmuchar1,x3,y1,y2,y3
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_coup7/alup,aldo
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
*
      xmusdn(1) = asdown1**2/mgluino**2
      xmusdn(2) = asdown2**2/mgluino**2
      xmusup(1) = asup1**2/mgluino**2
      xmusup(2) = asup2**2/mgluino**2

      dsup(1)=1-x1-xmusup(1)
      dsup(2)=1-x1-xmusup(2)
      dsdn(1)=1-x2-xmusdn(1)
      dsdn(2)=1-x2-xmusdn(2)

      xmuchar1 = amchar(nj)**2/mgluino**2

      x3 = 2-x1-x2
      y1 = 1-xmuchar1-x1
      y2 = 1-xmuchar1-x2
      y3 = 1+xmuchar1-x3

      do i=1,2,1
         blup(1,i) = 0.D0
         blup(2,i) = 0.D0
         bldo(1,i) = 0.D0
         bldo(2,i) = 0.D0
      end do

      NS_goud=0.D0
      
      if (mgluino.gt.amchar(nj)) then
         do i=1,2
            do k=1,2
               if (xmusdn(i).gt.1.0d0.and.xmusdn(k).gt.1.0d0)Then
               NS_goud=NS_goud+4.D0*g3s*g2s/dsdn(i)/dsdn(k)*2.D0*(
     .              (gdr(i)*gdr(k)+gdl(i)*gdl(k))*
     .              (aldo(i,nj)*aldo(k,nj)+
     .               bldo(k,nj)*bldo(i,nj))*x2*y2)
               endif
               if (xmusup(i).gt.1.0d0.and.xmusup(k).gt.1.0d0)Then
               NS_goud=NS_goud
     .         +4.D0*g3s*g2s/dsup(i)/dsup(k)*2.D0*(
     .              (gur(i)*gur(k)+gul(i)*gul(k))*
     .              (alup(i,nj)*alup(k,nj)+
     .               blup(k,nj)*blup(i,nj))*x1*y1)
               endif
               if (xmusdn(i).gt.1.0d0.and.xmusup(k).gt.1.0d0)Then
               NS_goud=NS_goud
     .         +4.D0*g3s*g2s/dsdn(i)/dsup(k)*2.D0*(
     .         (blup(k,nj)*gdr(i)*gul(k)*aldo(i,nj)+gur(k)*
     .          bldo(i,nj)*alup(k,nj)*gdl(i))*(-x1*y1-x2*y2+x3*y3)
     .         +2.D0*(gul(k)*alup(k,nj)*gdl(i)*aldo(i,nj)+
     .                gur(k)*blup(k,nj)*gdr(i)*bldo(i,nj))*
     .          xmchar(nj)/mgluino*y3)
               endif
            enddo
         enddo
      else
         NS_goud=0.D0
      endif
      
*      write(*,*)'NS_goud',NS_goud
      end
c ==================================================================== c
c                  gluino --> chargino- top bottombar
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_gotb(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,I,K
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2) 
      DOUBLE PRECISION dsbt(2),dstp(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),
     .gdr(2),gdl(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION x1,x2,uh,th,xmub,xmut,db1,db2,db3,db4,ab1,
     .ab2,ab3,ab4,dt1,dt2,
     .dt3,dt4,at1,at2,at3,at4,ab5,ab6,ab7,ab8,ab9,ab10,ab11,ab12,
     .xmusb(2),xmust(2),xmuchar1,gmst(2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

      gmst(1) = ast1
      gmst(2) = ast2
      xmusb(1) = asb1**2/mgluino**2
      xmusb(2) = asb2**2/mgluino**2
      xmust(1) = ast1**2/mgluino**2
      xmust(2) = ast2**2/mgluino**2

      dstp(1)=1-x1-xmust(1)+amt**2/mgluino**2
      dstp(2)=1-x1-xmust(2)+amt**2/mgluino**2
      dsbt(1)=1-x2-xmusb(1)+amb**2/mgluino**2
      dsbt(2)=1-x2-xmusb(2)+amb**2/mgluino**2

      xmuchar1 = amchar(nj)**2/mgluino**2

      NS_gotb=0.D0

      uh = 1.D0-x1+amt**2/mgluino**2
      th = 1.D0-x2+amb**2/mgluino**2
      xmub = amb**2/mgluino**2
      xmut = amt**2/mgluino**2
      
      if(mgluino.gt.(amchar(nj)+amb+amt)) then 
         do i=1,2
            do k=1,2
               db1 = alsbot(k,nj)*aksbot(i,nj)+alsbot(i,nj)*aksbot(k,nj)
               db2 = alsbot(k,nj)*alsbot(i,nj)+aksbot(k,nj)*aksbot(i,nj)
               db3 = gbl(k)*gbr(i)+gbl(i)*gbr(k)
               db4 = gbl(i)*gbl(k)+gbr(i)*gbr(k)

               ab1 = db1*db3
               ab2 = db1*db4
               ab3 = db2*db3
               ab4 = db2*db4

               dt1 = alstor(k,nj)*akstor(i,nj)+alstor(i,nj)*akstor(k,nj)
               dt2 = alstor(k,nj)*alstor(i,nj)+akstor(k,nj)*akstor(i,nj)
               dt3 = gtl(k)*gtr(i)+gtl(i)*gtr(k)
               dt4 = gtl(i)*gtl(k)+gtr(i)*gtr(k)

               at1 = dt1*dt3
               at2 = dt1*dt4
               at3 = dt2*dt3
               at4 = dt2*dt4

               ab5 = alstor(k,nj)*aksbot(i,nj)*gtl(k)*gbr(i)+
     .               alsbot(i,nj)*akstor(k,nj)*gbl(i)*gtr(k)
               ab6 = alstor(k,nj)*aksbot(i,nj)*gtr(k)*gbr(i)+
     .               akstor(k,nj)*alsbot(i,nj)*gtl(k)*gbl(i)
               ab7 = alstor(k,nj)*alsbot(i,nj)*gbr(i)*gtl(k)+
     .               akstor(k,nj)*aksbot(i,nj)*gbl(i)*gtr(k)
               ab8 = alstor(k,nj)*alsbot(i,nj)*gtr(k)*gbr(i)+
     .               akstor(k,nj)*aksbot(i,nj)*gbl(i)*gtl(k)
               ab9 = alstor(k,nj)*aksbot(i,nj)*gtl(k)*gbl(i)+
     .               alsbot(i,nj)*akstor(k,nj)*gbr(i)*gtr(k)
               ab10 = alstor(k,nj)*aksbot(i,nj)*gbl(i)*gtr(k)+
     .                alsbot(i,nj)*akstor(k,nj)*gtl(k)*gbr(i)
               ab11 = alstor(k,nj)*alsbot(i,nj)*gtl(k)*gbl(i)+
     .                akstor(k,nj)*aksbot(i,nj)*gbr(i)*gtr(k)
               ab12 = alstor(k,nj)*alsbot(i,nj)*gtr(k)*gbl(i)+
     .                akstor(k,nj)*aksbot(i,nj)*gtl(k)*gbr(i)

               if (xmusb(i).gt.1.d0.and.xmusb(k).gt.1.d0)Then
               NS_gotb=NS_gotb               
     .           +8.D0*g3s*g2s/dsbt(i)/dsbt(k)*(
     .           -4.D0*ab1*xmchar(nj)/mgluino*amb*amt/mgluino**2
     .           +2.D0*ab2*xmchar(nj)*amt/mgluino**2*(-xmub-1.D0+th)
     .           +2.D0*ab3*amb/mgluino*(-xmut-xmuchar1+th)
     .           +ab4*(-th**2+th*(1.D0+xmuchar1+xmub+xmut)-
     .                 (xmuchar1+xmut)*(1.D0+xmub)) )
               endif
        if ((gmst(i)+amt).gt.mgluino.and.(gmst(k)+amt).gt.mgluino)Then
               NS_gotb=NS_gotb
     .         +8.D0*g3s*g2s/dstp(i)/dstp(k)*(
     .           -4.D0*at1*xmchar(nj)/mgluino*amb*amt/mgluino**2
     .           +2.D0*at2*xmchar(nj)*amb/mgluino**2*(-xmut-1.D0+uh)
     .           +2.D0*at3*amt/mgluino*(-xmub-xmuchar1+uh)
     .           +at4*(-uh**2+uh*(1.D0+xmuchar1+xmub+xmut)-
     .                 (xmuchar1+xmub)*(1.D0+xmut)) )
               endif
               if (xmusb(i).gt.1.d0.and.(gmst(k)+amt).gt.mgluino)Then
               NS_gotb=NS_gotb
     .         -2.D0*8.D0*g3s*g2s/dsbt(i)/dstp(k)*(
     .           ab5*amb*amt/mgluino**2*(uh+th-xmub-xmut) 
     .           +ab6*amb/mgluino*(th-xmuchar1-xmut) 
     .           +ab7*xmchar(nj)*amb/mgluino**2*(uh-xmut-1.D0) 
     .           -2.D0*ab8*xmchar(nj)/mgluino*amt*amb/mgluino**2
     .           +ab9*amt/mgluino*(uh-xmuchar1-xmub)
     .           +ab10*(uh*th-xmuchar1-xmub*xmut)
     .           +ab11*xmchar(nj)/mgluino*(uh+th-xmuchar1-1.D0)
     .           +ab12*xmchar(nj)*amt/mgluino**2*(th-xmub-1.D0) )
               endif
            enddo
         enddo
      else 
         NS_gotb=0.D0
      endif
      end
c----------------------------------------------------------------------c
      DOUBLE PRECISION FUNCTION NS_gocc(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION x1,x2,uh,th,xmub,xmut,db1,db2,db3,db4,ab1,
     .ab2,ab3,ab4,dt1,dt2,
     .dt3,dt4,at1,at2,at3,at4,ab5,ab6,ab7,ab8,ab9,ab10,ab11,ab12,
     .xmusb(2),xmust(2),xmuchar1,gmst(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsbt(2),dstp(2)
      DOUBLE PRECISION gtr(2),gtl(2),gbr(2),gbl(2),gur(2),
     .gul(2),gdr(2),gdl(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1,csmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

      gmst(1) = ast1
      gmst(2) = ast2
      xmusb(1) = asb1**2/mgluino**2
      xmusb(2) = asb2**2/mgluino**2
      xmust(1) = ast1**2/mgluino**2
      xmust(2) = ast2**2/mgluino**2

      dstp(1)=1-x2-xmust(1)+amt**2/mgluino**2
      dstp(2)=1-x2-xmust(2)+amt**2/mgluino**2
      dsbt(1)=1-x1-xmusb(1)+amb**2/mgluino**2
      dsbt(2)=1-x1-xmusb(2)+amb**2/mgluino**2

      xmuchar1 = amchar(nj)**2/mgluino**2

      NS_gocc=0.D0

      uh = 1.D0-x1+amb**2/mgluino**2
      th = 1.D0-x2+amt**2/mgluino**2
      xmub = amb**2/mgluino**2
      xmut = amt**2/mgluino**2
      
      if(mgluino.gt.(amchar(nj)+amb+amt)) then 
         do i=1,2
            do k=1,2
               db1 = alsbot(k,nj)*aksbot(i,nj)+alsbot(i,nj)*aksbot(k,nj)
               db2 = alsbot(k,nj)*alsbot(i,nj)+aksbot(k,nj)*aksbot(i,nj)
               db3 = gbl(k)*gbr(i)+gbl(i)*gbr(k)
               db4 = gbl(i)*gbl(k)+gbr(i)*gbr(k)

               ab1 = db1*db3
               ab2 = db1*db4
               ab3 = db2*db3
               ab4 = db2*db4

               dt1 = alstor(k,nj)*akstor(i,nj)+alstor(i,nj)*akstor(k,nj)
               dt2 = alstor(k,nj)*alstor(i,nj)+akstor(k,nj)*akstor(i,nj)
               dt3 = gtl(k)*gtr(i)+gtl(i)*gtr(k)
               dt4 = gtl(i)*gtl(k)+gtr(i)*gtr(k)

               at1 = dt1*dt3
               at2 = dt1*dt4
               at3 = dt2*dt3
               at4 = dt2*dt4

               ab5 = alstor(k,nj)*aksbot(i,nj)*gtl(k)*gbr(i)+
     .               alsbot(i,nj)*akstor(k,nj)*gbl(i)*gtr(k)
               ab6 = alstor(k,nj)*aksbot(i,nj)*gtl(k)*gbl(i)+
     .               akstor(k,nj)*alsbot(i,nj)*gtr(k)*gbr(i)
               ab7 = alstor(k,nj)*alsbot(i,nj)*gbl(i)*gtr(k)+
     .               akstor(k,nj)*aksbot(i,nj)*gbr(i)*gtl(k)
               ab8 = alstor(k,nj)*alsbot(i,nj)*gtr(k)*gbr(i)+
     .               akstor(k,nj)*aksbot(i,nj)*gbl(i)*gtl(k)
               ab9 = alstor(k,nj)*aksbot(i,nj)*gtr(k)*gbr(i)+
     .               alsbot(i,nj)*akstor(k,nj)*gbl(i)*gtl(k)
               ab10 = alstor(k,nj)*aksbot(i,nj)*gbl(i)*gtr(k)+
     .                alsbot(i,nj)*akstor(k,nj)*gtl(k)*gbr(i)
               ab11 = alstor(k,nj)*alsbot(i,nj)*gtl(k)*gbl(i)+
     .                akstor(k,nj)*aksbot(i,nj)*gbr(i)*gtr(k)
               ab12 = alstor(k,nj)*alsbot(i,nj)*gtl(k)*gbr(i)+
     .                akstor(k,nj)*aksbot(i,nj)*gtr(k)*gbl(i)
               if (xmusb(i).gt.1.d0.and.xmusb(k).gt.1.d0)Then
               NS_gocc=NS_gocc
     .           +8.D0*g3s*g2s/dsbt(i)/dsbt(k)*(
     .           -4.D0*ab1*xmchar(nj)/mgluino*amb*amt/mgluino**2
     .           +2.D0*ab2*xmchar(nj)*amt/mgluino**2*(-xmub-1.D0+uh)
     .           +2.D0*ab3*amb/mgluino*(-xmut-xmuchar1+uh)
     .           +ab4*(-uh**2+uh*(1.D0+xmuchar1+xmub+xmut)-
     .                 (xmuchar1+xmut)*(1.D0+xmub)) )
               endif
        if ((gmst(i)+amt).gt.mgluino.and.(gmst(k)+amt).gt.mgluino)Then
               NS_gocc=NS_gocc
     .         +8.D0*g3s*g2s/dstp(i)/dstp(k)*(
     .           -4.D0*at1*xmchar(nj)/mgluino*amb*amt/mgluino**2
     .           +2.D0*at2*xmchar(nj)*amb/mgluino**2*(-xmut-1.D0+th)
     .           +2.D0*at3*amt/mgluino*(-xmub-xmuchar1+th)
     .           +at4*(-th**2+th*(1.D0+xmuchar1+xmub+xmut)-
     .                 (xmuchar1+xmub)*(1.D0+xmut)) )
               endif
               if (xmusb(i).gt.1.d0.and.(gmst(k)+amt).gt.mgluino)Then
               NS_gocc=NS_gocc
     .         -2.D0*8.D0*g3s*g2s/dsbt(i)/dstp(k)*(
     .           ab5*amb*amt/mgluino**2*(uh+th-xmub-xmut) 
     .           +ab9*amb/mgluino*(uh-xmuchar1-xmut) 
     .           +ab12*xmchar(nj)*amb/mgluino**2*(th-xmut-1.D0) 
     .           -2.D0*ab8*xmchar(nj)/mgluino*amt*amb/mgluino**2
     .           +ab6*amt/mgluino*(th-xmuchar1-xmub)
     .           +ab10*(uh*th-xmuchar1-xmub*xmut)
     .           +ab11*xmchar(nj)/mgluino*(uh+th-xmuchar1-1.D0)
     .           +ab7*xmchar(nj)*amt/mgluino**2*(uh-xmub-1.D0) )
               endif
            enddo
         enddo
      else 
         NS_gocc=0.D0
      endif
      end
