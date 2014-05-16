c ------------------------------------------------------------------ c
c This routine serves to
c 1) translate parameters, masses, mixing angles and couplings
c    from NMSSMTools into the conventions originally used in
c     SDECAY: A Fortran code for the decays of the supersymmetric 
c         particles in the MSSM
c     by M. Muhlleitner (Karlsruhe, Inst. Technol.),
c        A. Djouadi (Orsay, LPT & CERN, Theory Division),
c        Y. Mambrini (Orsay, LPT),
c     Comput.Phys.Commun.168:46-70 (2005), hep-ph/0311167.
c    SDECAY should be cited whenever NMSDECAY is used.
c 2) call the subroutines (in separate files) for the various
c      sparticle decay widths
c 3) call the subroutine for the output
c ------------------------------------------------------------------ c
c
      SUBROUTINE NMSDECAY(PAR)
*      
      IMPLICIT NONE
*
      INTEGER I,J,K,nx1t,ny1t,A,B,D,NMSFLAG
      DOUBLE PRECISION PAR(*)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .     MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .     CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION Q2
      DOUBLE PRECISION MHUS,MHDS,MHSS
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ 
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asne2,asntau2
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION sdmhd2,sdmuq,sdmsq,sdmbr,sdmdr
      DOUBLE PRECISION SCALb,SCALt,scaltau,gs2 
      DOUBLE PRECISION flagmulti,flagqcd,flagloop,multilim
      DOUBLE PRECISION COSBETA,SINBETA
      DOUBLE PRECISION P(2,3)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION rmtc,rmbc,rmtauc
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION amuv,lamv,amuref
      DOUBLE PRECISION au,ad,atau_yu,amu
C ----COUPLINGS -----------
      DOUBLE PRECISION Hstaustaur(3,2,2),Astaustaur(2,2,2)
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION GHNEUNEU,GANEUNEU,GHCHACHA,GACHACHA
      DOUBLE PRECISION GHCNEUCHAL,GHCNEUCHAR
      DOUBLE PRECISION HRULUL,HRDLDL,HRURUR,HRDRDR,HRULUR
      DOUBLE PRECISION HRDLDR,HRLLLL,HRLRLR,HRLLLR
      DOUBLE PRECISION HIULUR,HIDLDR,HILLLR
      DOUBLE PRECISION LAMBDA,KAPPA,ALQ,AKQ,MUEFFQ,NUQ
      DOUBLE PRECISION HTAU12,HTAU22,HTAU11
      DOUBLE PRECISION ATAU12,ATAU22,ATAU11
      DOUBLE PRECISION HISB12,HISB11,HISB22
      DOUBLE PRECISION AISB12,AISB11,AISB22
      DOUBLE PRECISION HIST12,HIST11,HIST22
      DOUBLE PRECISION AIST12,AIST11,AIST22
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5)
      DOUBLE PRECISION hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .     achachaL(2,2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ale(2,2),altau(2,2),alsne(2,2),blsne(2,2),
     .     alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION ae(2,5),be(2,5),atau(2,5),btau(2,5),anu(2,5),
     .     bnu(2,5),antau(2,5),bntau(2,5)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION s11,s12,s21,s22
      DOUBLE PRECISION chctaunur,chctaunul
      DOUBLE PRECISION achtop,vchtop,achtau,vchtau
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .     azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2)
      DOUBLE PRECISION gul(2),gur(2),gdl(2),gdr(2),gtl(2),
     .     gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION dtlttr(3,2,2),dtattr(2,2,2),
     .datlttr(3,2,2),datattr(2,2,2),
     .dthlttr(3,2,2),dthattr(2,2,2)
      DOUBLE PRECISION dtHIST12,dtHIST11,dtHIST22
      DOUBLE PRECISION dtAIST12,dtAIST11,dtAIST22
      DOUBLE PRECISION datHIST12,datHIST11,datHIST22
      DOUBLE PRECISION datAIST12,datAIST11,datAIST22
      DOUBLE PRECISION dmixHIST12,dmixHIST11,dmixHIST22
      DOUBLE PRECISION dmixAIST12,dmixAIST11,dmixAIST22
      DOUBLE PRECISION dtHISB12,dtHISB11,dtHISB22
      DOUBLE PRECISION dtAISB12,dtAISB11,dtAISB22
      DOUBLE PRECISION dabHISB12,dabHISB11,dabHISB22
      DOUBLE PRECISION dabAISB12,dabAISB11,dabAISB22
      DOUBLE PRECISION dmixHISB12,dmixHISB11,dmixHISB22
      DOUBLE PRECISION dmixAISB12,dmixAISB11,dmixAISB22
      DOUBLE PRECISION dblbbr(3,2,2),dbabbr(2,2,2),
     .dablbbr(3,2,2),dababbr(2,2,2),
     .dthlbbr(3,2,2),dthabbr(2,2,2)
      DOUBLE PRECISION chtbrunr,chtbrunl
      DOUBLE PRECISION Httr(3),Attr(2)
      DOUBLE PRECISION Hbbr(3),Abbr(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2)
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION chcsntau(2,2),gcsntaur(2,2)
      DOUBLE PRECISION chctb(2,2),gctbr(2,2)
c
c
c ----- COMMON BLOCKS FROM NMSSMTOOLS ---------------- c
      COMMON/NMSFLAG/NMSFLAG
c Couplings at MZ:
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
c Couplings at the Susy scale Q2 (gauge couplings squared!):
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS 
c SM masses
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
c Higgs spectrum and mixing angles
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
c Sparticle spectrum and mixing angles
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .     MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .     CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
c Susy scale Q2 = Q^2
      COMMON/RENSCALE/Q2
c Soft Higgs mass terms at the Susy scale
      COMMON/SUSYMH/MHUS,MHDS,MHSS
c Parameters at the scale QSTSB ~ Q:
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
c ---------- END from NMSSMTools--------------------------- c
c
      COMMON/NS_pi/PI,SQR2
c Flags for 3-body decays, QCD corrections, loop decays:
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop
      COMMON/NS_multilim/multilim
c For QCD corrections:
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_refscale/amuref
c For integration:
      COMMON/NS_nx1/nx1t,ny1t
c MASSES:
c   MZ, MW at the Susy scale Q2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
C   Gaugino/neutralino/chargino masses + absolute values:
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
c   R.h. sneutrino:
      COMMON/NS_rhsneutr/asne2,asntau2
c   top, bottom, tau:
      COMMON/NS_runmcalc/rmtc,rmbc,rmtauc
c Mixing angles:
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_CPodd_MIX/P
c Parameters:
c   Tan(beta) at QSTSB ~ Susy scale Q2
      COMMON/NS_tanb/tanbeta
c   A-terms, mu_eff:
      COMMON/NS_trilin_mu/au,ad,atau_yu,amu
c Soft mass terms for rad. decays:
      COMMON/NS_hikasakob/sdmuq,sdmhd2
      COMMON/NS_hikasakob02/sdmsq,sdmbr,sdmdr
c
C   -------COUPLINGS-----------       
c Yukawas rescaled by g2, at the scale Q;
c    gs2 is the strong coupling squared:
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
c Chargino-slepton-lepton:
      COMMON/NS_coup5/ale,altau,alsne,blsne,alsnt,blsnt
c Chargino-squark-quark (up,down):
      COMMON/NS_coup7/alup,aldo
c Chargino-sbottom-top:
      COMMON/NS_charsbottop/alsbot,aksbot
c Chargino-stop-bottom:
      COMMON/NS_charstopbot/alstor,akstor
c Neutralino-slepton-lepton:
      COMMON/NS_coup8/ae,be,atau,btau,anu,bnu,antau,bntau
c Neutralino-squark-quark (up,down):
      COMMON/NS_coup10/aup,bup,ado,bdo
c Neutralino-sbottom-bottom:
      COMMON/NS_neutsbotbot/abot,bbot
c Neutralino-stop-top:
      COMMON/NS_neutstoptop/atopr,btopr
c Gluino:
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
c Neutral Higgs:  
c   charginos, neutralinos:
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .  achachaL
c Stops, sbottoms, staus:
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr  
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr 
      COMMON/NS_HIGGSTAUTAU/Hstaustaur,Astaustaur
c   top,bottom:
      COMMON/NS_phitoptop/Httr,Attr
      COMMON/NS_phibotbot/Hbbr,Abbr
c Charged Higgs:
c   neutralino-chargino:
      COMMON/NS_coup3/ql,qr,ol,or
c   tau left/right:
      COMMON/NS_coup14/chctaunur,chctaunul
c   top/bottom/tau:
      COMMON/NS_coup15/achtop,vchtop,achtau,vchtau
c   up/down quarks:
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
c   stau/sneutrino_tau:
      COMMON/NS_hcstausntau/gcsntaur
c   stop/sbottom:
      COMMON/NS_hcsbotstop/gctbr
c Z-boson:
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .  azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau
c W-boson:
      COMMON/NS_coup18/awff,vwff
      COMMON/NS_coup20/gwtb,gwntau
c
C  -- DERIVATIVE COUPLINGS 
      COMMON/NS_higgsst1st2deriv/
     .  dtlttr,dtattr,datlttr,datattr,dthlttr,dthattr
      COMMON/NS_higgssb1sb2deriv/dblbbr,dbabbr,
     .  dablbbr,dababbr,dthlbbr,dthabbr
c
c --------------------- END COMMONS -----------------------------
c
c -- FUNCTIONS FOR COUPLINGS as in decay.f, but at the Susy scale Q:
c
c    HIGGS-NEUTRALINOS/CHARGINOS:
        GHNEUNEU(A,B,D)= LAMBDA/SQR2
     .       * (S(A,1)*(zz(B,3)*zz(D,5)+zz(B,5)*zz(D,3))
     .       + S(A,2)*(zz(B,4)*zz(D,5)+zz(B,5)*zz(D,4))
     .       + S(A,3)*(zz(B,4)*zz(D,3)+zz(B,3)*zz(D,4)))
     .       - SQR2*KAPPA*S(A,3)*zz(B,5)*zz(D,5)
     .      + DSQRT(g1s)/2.D0*(-S(A,1)*(zz(B,1)*zz(D,4)+zz(B,4)*zz(D,1))
     .       + S(A,2)*(zz(B,1)*zz(D,3)+zz(B,3)*zz(D,1)))
     .       + DSQRT(g2s)/2.D0*(S(A,1)*(zz(B,2)*zz(D,4)+zz(B,4)*zz(D,2))
     .       - S(A,2)*(zz(B,2)*zz(D,3)+zz(B,3)*zz(D,2)))
        GANEUNEU(A,B,D)= LAMBDA/SQR2
     .       * (P(A,1)*(zz(B,3)*zz(D,5)+zz(B,5)*zz(D,3))
     .       + P(A,2)*(zz(B,4)*zz(D,5)+zz(B,5)*zz(D,4))
     .       + P(A,3)*(zz(B,4)*zz(D,3)+zz(B,3)*zz(D,4)))
     .       - SQR2*KAPPA*P(A,3)*zz(B,5)*zz(D,5)
     .       - DSQRT(g1s)/2.D0*(-P(A,1)*(zz(B,1)*zz(D,4)
     .       +zz(B,4)*zz(D,1))
     .       + P(A,2)*(zz(B,1)*zz(D,3)+zz(B,3)*zz(D,1)))
     .       - DSQRT(g2s)/2.D0*(P(A,1)*(zz(B,2)*zz(D,4)+zz(B,4)*zz(D,2))
     .       - P(A,2)*(zz(B,2)*zz(D,3)+zz(B,3)*zz(D,2)))
        GHCHACHA(A,B,D)= LAMBDA/SQR2*S(A,3)*U(B,2)*V(D,2)
     .       + DSQRT(g2s)/SQR2*(S(A,1)*U(B,1)*V(D,2)+
     .       S(A,2)*U(B,2)*V(D,1))
        GACHACHA(A,B,D)= LAMBDA/SQR2*P(A,3)*U(B,2)*V(D,2)
     .       - DSQRT(g2s)/SQR2*(P(A,1)*U(B,1)*V(D,2)+
     .       P(A,2)*U(B,2)*V(D,1))
        GHCNEUCHAL(A,B)= LAMBDA*COSBETA*zz(A,5)*U(B,2)
     .       - SINBETA/SQR2*(DSQRT(g1s)*zz(A,1)+
     .       DSQRT(g2s)*zz(A,2))*U(B,2)
     .       + SINBETA*DSQRT(g2s)*zz(A,3)*U(B,1)
        GHCNEUCHAR(A,B)= LAMBDA*SINBETA*zz(A,5)*V(B,2)
     .       + COSBETA/SQR2*(DSQRT(g1s)*zz(A,1)+
     .         DSQRT(g2s)*zz(A,2))*V(B,2)
     .       + COSBETA*DSQRT(g2s)*zz(A,4)*V(B,1)
c
c HIGGS A - STAU(1,2) - STAU(1,2)
c as in decay.f, but rescaled by 1/DSQRT(G2s)
        HRLLLL(A)= -(1/DSQRT(G2s))*(SQR2*(HTAUS**2*H2Q*S(A,2)
     .       + (-g1s/4.D0+g2s/4.D0)*(H1Q*S(A,1)-H2Q*S(A,2))))
        HRLRLR(A)= -(1/DSQRT(G2s))*(SQR2*(HTAUS**2*H2Q*S(A,2)
     .       + g1s/2.D0*(H1Q*S(A,1)-H2Q*S(A,2))))
        HRLLLR(A)= -(HTAUS/SQR2*(-MUEFFQ*S(A,1)+ATAU_YU*S(A,2)
     .       - LAMBDA*H1Q*S(A,3)))/DSQRT(G2s)
        HILLLR(A)= -(HTAUS/SQR2*(MUEFFQ*P(A,1)+ATAU_YU*P(A,2)
     .       + LAMBDA*H1Q*P(A,3)))/DSQRT(G2s)       
        HTAU12(A) = CL*SL*(HRLRLR(A)-HRLLLL(A))+
     .       (CL**2-SL**2)*HRLLLR(A)
        HTAU11(A) = CL**2*HRLLLL(A)+SL**2*HRLRLR(A)
     .       +2*CL*SL*HRLLLR(A)
        HTAU22(A) = SL**2*HRLLLL(A)+CL**2*HRLRLR(A)
     .       -2*CL*SL*HRLLLR(A)
        ATAU12(A) = (CL**2-SL**2)*HILLLR(A)
        ATAU11(A) = 2*CL*SL*HILLLR(A)
        ATAU22(A) = -2*CL*SL*HILLLR(A)
c
c HIGGS A - SBOTTOM(1,2) - SBOTTOM(1,2)
c as in decay.f, but rescaled by 1/DSQRT(G2s)
        HRDLDL(A)= -(1.0/DSQRT(G2s))*(SQR2*(HBOTS**2*H2Q*S(A,2)
     .       + (g1s/12.D0+g2s/4.D0)*(H1Q*S(A,1)-H2Q*S(A,2))))
        HRDRDR(A)= -(1.0/DSQRT(G2s))*(SQR2*(HBOTS**2*H2Q*S(A,2)
     .       + g1s/6.D0*(H1Q*S(A,1)-H2Q*S(A,2))))
        HRDLDR(A)= -(HBOTS/SQR2*(-MUEFFQ*S(A,1)+AD*S(A,2)
     .       - LAMBDA*H1Q*S(A,3)))/DSQRT(G2s)
        HIDLDR(A)= -(HBOTS/SQR2*(MUEFFQ*P(A,1)+AD*P(A,2)
     .       + LAMBDA*H1Q*P(A,3)))/DSQRT(G2s)
        HISB12(A) = CB*SB*(HRDRDR(A)-HRDLDL(A))+
     .       (CB**2-SB**2)*HRDLDR(A)
        HISB11(A) = CB**2*HRDLDL(A)+SB**2*HRDRDR(A)
     .       +2*CB*SB*HRDLDR(A)
        HISB22(A) = SB**2*HRDLDL(A)+CB**2*HRDRDR(A)
     .       -2*CB*SB*HRDLDR(A)
        AISB12(A)= (CB**2-SB**2)*HIDLDR(A)
        AISB11(A)= 2*CB*SB*HIDLDR(A)
        AISB22(A)= -2*CB*SB*HIDLDR(A)
c
c HIGGS A - STOP(1,2) - STOP(1,2)
c as in decay.f, but rescaled by 1/DSQRT(G2s)
        HRULUL(A)= -(1/DSQRT(G2s))*(SQR2*(HTOPS**2*H1Q*S(A,1)
     .       + (g1s/12.D0-g2s/4.D0)*(H1Q*S(A,1)-H2Q*S(A,2))))
        HRURUR(A)= -(1/DSQRT(G2s))*(SQR2*(HTOPS**2*H1Q*S(A,1)
     .       - g1s/3.D0*(H1Q*S(A,1)-H2Q*S(A,2))))
        HRULUR(A)= -(HTOPS/SQR2*(AU*S(A,1)-MUEFFQ*S(A,2)
     .       - LAMBDA*H2Q*S(A,3)))/DSQRT(G2s)
        HIULUR(A)= -(HTOPS/SQR2*(AU*P(A,1)+MUEFFQ*P(A,2)
     .       + LAMBDA*H2Q*P(A,3)))/DSQRT(G2s)
        HIST12(A)= CT*ST*(HRURUR(A)-HRULUL(A))+
     .       (CT**2-ST**2)*HRULUR(A)
        HIST11(A)= CT**2*HRULUL(A)+ ST**2*HRURUR(A)
     .       +2*CT*ST*HRULUR(A)
        HIST22(A)= ST**2*HRULUL(A)+ CT**2*HRURUR(A)
     .       -2*CT*ST*HRULUR(A) 
        AIST12(A)= (CT**2-ST**2)*HIULUR(A)
        AIST11(A)= 2*CT*ST*HIULUR(A)
        AIST22(A)= -2*CT*ST*HIULUR(A)
c
C  -- DERIVATIVES OF  COULINGS -- 
c
C  --- Stop1-stop2-Higgs --		
c ---- derivatives d/dmtop ----
        dtHIST12(A)= -(CT**2-ST**2)/H1Q/dsqrt(2.D0)/dsqrt(g2s)*
     .  (AU*S(A,1)-MUEFFQ*S(A,2) - LAMBDA*H2Q*S(A,3))
        dtHIST11(A)= dsqrt(2.D0)*scalt*S(A,1)-
     .  2*ST*CT/dsqrt(2.D0)/H1Q*(AU*S(A,1)-MUEFFQ*S(A,2)
     .  -LAMBDA*H2Q*S(A,3))
        dtHIST22(A)=dsqrt(2.D0)*scalt*S(A,1)+
     .  2*ST*CT/dsqrt(2.D0)/H1Q*(AU*S(A,1)-MUEFFQ*S(A,2)
     .  -LAMBDA*H2Q*S(A,3))
        dtAIST12(A)= -(CT**2-ST**2)/H1Q/dsqrt(2.D0)/dsqrt(g2s)*
     .  (AU*P(A,1)+MUEFFQ*P(A,2)+ LAMBDA*H2Q*P(A,3))
        dtAIST11(A)= -2*CT*ST/H1Q/dsqrt(2.D0)/dsqrt(g2s)*
     .  (AU*P(A,1)+MUEFFQ*P(A,2)+ LAMBDA*H2Q*P(A,3))
        dtAIST22(A)= 2*CT*ST/H1Q/dsqrt(2.D0)/dsqrt(g2s)*
     .  (AU*P(A,1)+MUEFFQ*P(A,2)+ LAMBDA*H2Q*P(A,3))
C ---- derivatives d/dAt ----
        datHIST12(A)= -(CT**2-ST**2)*scalt*S(A,1)/dsqrt(2.D0)
        datHIST11(A)= -2*CT*ST*scalt*S(A,1)/dsqrt(2.D0)
        datHIST22(A)=  2*CT*ST*scalt*S(A,1)/dsqrt(2.D0)
        datAIST12(A)= -(CT**2-ST**2)*scalt/dsqrt(2.D0)*P(A,1)
        datAIST11(A)= -2*CT*ST*scalt/dsqrt(2.D0)*P(A,1)
        datAIST22(A)=  2*CT*ST*scalt/dsqrt(2.D0)*P(A,1)
C ---- derivatives theta_t ----
        dmixHIST12(A)= (CT**2-ST**2)*(HRURUR(A)-HRULUL(A))
     .  -4*ST*CT*HRULUR(A)
        dmixHIST11(A)=2*ST*CT*(HRURUR(A)-HRULUL(A))
     .  +2*(CT**2-ST**2)*HRULUR(A)
        dmixHIST22(A)=-2*ST*CT*(HRURUR(A)-HRULUL(A))
     .  -2*(CT**2-ST**2)*HRULUR(A)
        dmixAIST12(A)= -4*CT*ST*HIULUR(A)
        dmixAIST11(A)= 2*(CT**2-ST**2)*HIULUR(A)
        dmixAIST22(A)= -2*(CT**2-ST**2)*HIULUR(A)
C  ------ sbottom1-sbottom2-Higgs couplings
C  ---- derivatives d/dmbottom ----
        dtHISB12(A)= (CB**2-SB**2)/H2Q/dsqrt(2.D0)/dsqrt(g2s)*
     .  (-MUEFFQ*S(A,1)+AD*S(A,2)
     .  - LAMBDA*H1Q*S(A,3))
        dtHISB11(A)= dsqrt(2.D0)*scalb*S(A,2)+
     .  2*CB*SB/dsqrt(2.D0)/H2Q*(-MUEFFQ*S(A,1)+AD*S(A,2)
     .  - LAMBDA*H1Q*S(A,3))
        dtHISB22(A)= dsqrt(2.D0)*scalb*S(A,2)-
     .  2*CB*SB/dsqrt(2.D0)/H2Q*(-MUEFFQ*S(A,1)+AD*S(A,2)
     .  -LAMBDA*H1Q*S(A,3))
        dtAISB12(A)= (CB**2-SB**2)/H2Q/dsqrt(2.D0)/dsqrt(g2s)*
     .  (MUEFFQ*P(A,1)+AD*P(A,2)+ LAMBDA*H1Q*P(A,3))
        dtAISB11(A)= 2*CB*SB/H2Q/dsqrt(2.D0)/dsqrt(g2s)*
     .  (MUEFFQ*P(A,1)+AD*P(A,2)+ LAMBDA*H1Q*P(A,3))
        dtAISB22(A)= -2*CB*SB/H2Q/dsqrt(2.D0)/dsqrt(g2s)*
     .  (MUEFFQ*P(A,1)+AD*P(A,2)+ LAMBDA*H1Q*P(A,3))
c ---- derivative d/dAb ----
        dabHISB12(A)=  (CB**2-SB**2)*scalb*S(A,2)/dsqrt(2.D0)
        dabHISB11(A)= 2*CB*SB*scalb*S(A,2)/dsqrt(2.D0)
        dabHISB22(A)= -2*CB*SB*scalb*S(A,2)/dsqrt(2.D0)
        dabAISB12(A)= (CB**2-SB**2)*scalb*P(A,2)/dsqrt(2.D0)
        dabAISB11(A)= 2*CB*SB*scalb*P(A,2)/dsqrt(2.D0)
        dabAISB22(A)= -2*CB*SB*scalb*P(A,2)/dsqrt(2.D0)
c ---- derivative theta_b ----
        dmixHISB12(A)= (CB**2-SB**2)*(HRDRDR(A)-HRDLDL(A))
     .  -4*CB*SB*HRDLDR(A)
        dmixHISB11(A)= 2*CB*SB*(HRDRDR(A)-HRDLDL(A))
     .  +2*(CB**2-SB**2)*HRDLDR(A)
        dmixHISB22(A)= -2*CB*SB*(HRDRDR(A)-HRDLDL(A))
     .  -2*(CB**2-SB**2)*HRDLDR(A)
        dmixAISB12(A)= -4*CB*SB*HIDLDR(A)
        dmixAISB11(A)= 2*(CB**2-SB**2)*HIDLDR(A)
        dmixAISB22(A)= -2*(CB**2-SB**2)*HIDLDR(A)
c ----------------------- END OF FUNCTIONS ----------------------
c
c ------------------------ BEGIN CODE ----------------------------
         PI = 4.D0*DATAN(1.D0)
         IF(NMSFLAG.EQ.0) RETURN
         PI = 4.D0*DATAN(1.D0)
         SQR2= DSQRT(2.D0)
C ---- FLAGS ------------ put "0" to invalidate 
         flagmulti=1.D0      ! 3-body decays
         flagqcd=1.D0        ! QCD corrections to 2-body decays
         flagloop= 1.D0      ! loop decays
c --- 3-body BRs are shown only if larger than multilim:
         multilim=0.01d0
c -- amuref: scale of the couplings set to the Susy scale Q
         amuref = DSQRT(Q2)
c -- lamv and amuv as in Sdecay:
         lamv = 1.D-15
         amuv = 1.D23
c -- Parameters and couplings at the scale Q,
c    used in the Higgs coupling functions above 
         LAMBDA=PAR(1)
         KAPPA=PAR(2)
         ALQ=PAR(5)
         AKQ=PAR(6)
         MUEFFQ=PAR(4)
         NUQ=PAR(2)*PAR(4)/PAR(1)
c -- weak mixing angle at the scale Q
         cw=DSQRT(g2s/(g1s+g2s))
         sw=DSQRT(g1s/(g1s+g2s))
         tw=sw/cw
c -- tanbeta at the scale Q ~ QSTSB
         tanbeta = TANBQ
         COSBETA = 1.D0/DSQRT(1.D0+tanbeta**2)
         SINBETA = tanbeta*COSBETA
c -- MW, MZ at the scale Q ~ QSTSB
         AMW = dsqrt(1.D0/2.D0*g2s*(H1Q**2 + H2Q**2))
         AMZ = dsqrt(1.D0/2.D0*(g1s+g2s)*(H1Q**2 + H2Q**2))
c -- DRbar Yukawa couplings rescaled by the weak coupl. at the scale Q 
         scalb = HBOTS/DSQRT(g2s) 
         scalt = HTOPS/DSQRT(g2s)
         scaltau = HTAUS/dsqrt(g2s)
c -- gs2 = strong coupling squared:
         gs2=g3s
c -- Running mt,mb,mtau masses from Yukawas at Q (but H1Q, H2Q at QSTSB)
         rmtc = HTOPS*H1Q
         rmbc = HBOTS*H2Q
         rmtauc = HTAUS*H2Q
c -- CP odd mixing angles--consistent with decay.f
         DO I=1,2
         P(I,1)= P2(I,1)*COSBETA 
         P(I,2)= P2(I,1)*SINBETA  
         P(I,3)= P2(I,2)
         ENDDO
c -- gluino mass
         mgluino = MGL
c -- neutralino & chargino masses and absolute values
         amneut(1)=DABS(MNEU(1))
         amneut(2)=DABS(MNEU(2))
         amneut(3)=DABS(MNEU(3))
         amneut(4)=DABS(MNEU(4))
         amneut(5)=DABS(MNEU(5))
         xmneut(1)=MNEU(1)
         xmneut(2)=MNEU(2)
         xmneut(3)=MNEU(3)
         xmneut(4)=MNEU(4)
         xmneut(5)=MNEU(5)
         xmchar(1)=MCHA(1)
         xmchar(2)=MCHA(2)
         amchar(1)=DABS(MCHA(1))
         amchar(2)=DABS(MCHA(2))
*      
         DO i=1,2
            DO j=1,2
               uu(i,j)=U(i,j)
               vv(i,j)=V(i,j)
            ENDDO
         ENDDO
         DO i=1,5
            zz(i,1)=N(i,1)
            zz(i,2)=N(i,2)
            zz(i,3)=N(i,4)
            zz(i,4)=N(i,3)
            zz(i,5)=N(i,5)
         ENDDO
         DO k=1,5
            zp(k,1) =  zz(k,1)*cw+zz(k,2)*sw
            zp(k,2) = -zz(k,1)*sw+zz(k,2)*cw
            zp(k,3) = zz(k,3)
            zp(k,4) = zz(k,4)
            zp(k,5) = zz(k,5)
         ENDDO
c -- FOR INTEGRATION         
         nx1t  = 32
         ny1t  = 32
c -- RH sneutrino masses
         asne2   = 1.D15
         asntau2 = 1.D15
c -- Sfermion mixing angles
         sdthet=DACOS(CST)
         sdtheb=DACOS(CSB)
         sdthel=DACOS(CSL)
         ct=CST 
         st=DSQRT(1.D0-CST**2)
         cb=CSB
         sb=DSQRT(1.D0-CSB**2)
         cl=CSL
         sl=DSQRT(1.D0-CSL**2)
         cu = 1.D0
         su = 0.D0
         cd = 1.D0
         sd = 0.D0
         ce = 1.D0
         se = 0.D0
         cn = 1.D0
         sn = 0.D0
c Soft masses for COMMON/NS_hikasakob, COMMON/NS_hikasakob02:
         sdmhd2 = MHDS
         sdmuq = DSQRT(PAR(15))
         sdmsq = DSQRT(PAR(7))
         sdmbr = DSQRT(PAR(9))
         sdmdr = DSQRT(PAR(17))
c Trilinears and mu_eff for COMMON/NS_trilin_mu:
         au= PAR(12)
         ad= PAR(13)
         atau_yu= PAR(14)
         amu   = PAR(4)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c --- COUPLINGS IN SDECAY-CONVENTIONS -------------------
C
c -------------------------- Chargino --------------------------- c
c ---- CHARGINO K - SELECTRON - NEUTRINO_E 
          DO k=1,2
             ale(1,k) = -U(k,1)*ce
             ale(2,k) = U(k,1)*se
          ENDDO
c ---- CHARGINO K - STAU 1,2 - NEUTRINO_TAU 
          DO k=1,2
             altau(1,k) = -cl*U(k,1)+sl*Scaltau*U(k,2)
             altau(2,k) =  sl*U(k,1)+cl*Scaltau*U(k,2)
          ENDDO
c ---- CHARGINO K - SNEUTRINO_EL 1,2 - ELECTRON
          DO k=1,2
             alsne(1,k) = -V(k,1)
             alsne(2,k) = 0.D0
             blsne(1,k) = 0.D0
             blsne(2,k) = 0.D0
          ENDDO
c ---- CHARGINO K - SNEUTRINO_TAU 1,2 - TAU
          DO k=1,2
             alsnt(1,k) = -V(k,1)
             alsnt(2,k) = 0.D0
             blsnt(1,k) = Scaltau*U(k,2)
             blsnt(2,k) = 0.D0
          ENDDO
c ---- CHARGINO K - SUP 1,2 - DOWN
          DO k=1,2
             alup(1,k) = -V(k,1)*cu
             alup(2,k) =  V(k,1)*su
          ENDDO
c ---- CHARGINO K - SDOWN 1,2 - UP
          DO k=1,2
             aldo(1,k) = -U(k,1)*cd
             aldo(2,k) =  U(k,1)*sd
          ENDDO
c ---- CHARGINO K - SBOTTOM 1,2 - TOP
      DO k=1,2
         alsbot(1,k)=-cb*U(k,1)+sb*scalb*U(k,2)
         alsbot(2,k)= sb*U(k,1)+cb*scalb*U(k,2)
         aksbot(1,k)= cb*scalt*V(k,2)
         aksbot(2,k)= -sb*scalt*V(k,2)
      ENDDO
c ---- CHARGINO K - STOP 1,2 - BOTTOM
      DO k=1,2
         alstor(1,k)=-ct*V(k,1)+st*scalt*V(k,2)
         alstor(2,k)= st*V(k,1)+ct*scalt*V(k,2)
         akstor(1,k)=ct*scalb*U(k,2)
         akstor(2,k)= -st*scalb*U(k,2)
      ENDDO
c
c -------------------- Neutralino -------------------------- c
c ---- NEUTRALINO K - SELECTRON L,R - ELECTRON
          DO k=1,5
             ae(1,k) = ce*dsqrt(2.D0)*(zp(k,1)*sw+(0.5D0-sw**2)*
     .            zp(k,2)/cw)
             ae(2,k) = -se*dsqrt(2.D0)*(zp(k,1)*sw+(0.5D0-sw**2)*
     .            zp(k,2)/cw)
             be(1,k) = se*dsqrt(2.D0)*sw*(zp(k,2)*tw-zp(k,1))
             be(2,k) = ce*dsqrt(2.D0)*sw*(zp(k,2)*tw-zp(k,1))
          ENDDO
c ---- NEUTRALINO K - STAU - TAU
          DO k=1,5
             atau(1,k) = cl*dsqrt(2.D0)*(zp(k,1)*sw+(0.5D0-sw**2)*
     .            zp(k,2)/cw)-sl*scaltau*zz(k,3)
             atau(2,k) = -sl*dsqrt(2.D0)*(zp(k,1)*sw+(0.5D0-sw**2)*
     .            zp(k,2)/cw)-cl*scaltau*zz(k,3)
             btau(1,k) = sl*dsqrt(2.D0)*sw*(zp(k,2)*tw-zp(k,1))
     .            -cl*scaltau*zz(k,3)
             btau(2,k) = cl*dsqrt(2.D0)*sw*(zp(k,2)*tw-zp(k,1))
     .            +sl*scaltau*zz(k,3)
          ENDDO
c ---- NEUTRALINO K - SNEUTRINO_E - NEUTRINO_E
          DO k=1,5
             anu(1,k) = -zp(k,2)/dsqrt(2.D0)/cw
             anu(2,k) = 0.D0
             bnu(1,k) = 0.D0
             bnu(2,k) = 0.D0
          ENDDO
c ---- NEUTRALINO K - SNEUTRINO_TAU 1,2 - NEUTRINO_TAU
          DO k=1,5
             antau(1,k) = -zp(k,2)/dsqrt(2.D0)/cw
             antau(2,k) = 0.D0
             bntau(1,k) = 0.D0
             bntau(2,k) = 0.D0
          ENDDO
c ---- NEUTRALINO K - SUP 1,2 - UP
          DO k=1,5
             aup(1,k) = -cu*dsqrt(2.D0)*(2.D0/3.D0*zp(k,1)*sw+(0.5D0
     .            -2.D0/3.D0*sw**2)*zp(k,2)/cw)
             aup(2,k) = -su*dsqrt(2.D0)*(-2.D0*zp(k,1)*sw/3.D0+
     .            (-0.5D0+2.D0/3.D0*sw**2)*zp(k,2)/cw)
             bup(1,k) = -2.D0*su*dsqrt(2.D0)*sw*(zp(k,2)*tw-
     .            zp(k,1))/3.D0
             bup(2,k) = -2.D0*cu*dsqrt(2.D0)*sw*(zp(k,2)*tw-
     .            zp(k,1))/3.D0
          ENDDO
c ---- NEUTRALINO K - SDOWN 1,2 - DOWN
          DO k=1,5
             ado(1,k) = cd*dsqrt(2.D0)*(zp(k,1)*sw/3.D0+(0.5D0
     .            -1.D0/3.D0*sw**2)*zp(k,2)/cw)
             ado(2,k) = -sd*dsqrt(2.D0)*(zp(k,1)*sw/3.D0+(0.5D0
     .            -1.d0/3.D0*sw**2)*zp(k,2)/cw)
             bdo(1,k) = sd*dsqrt(2.D0)*sw*(zp(k,2)*tw-zp(k,1))/3.D0
             bdo(2,k) = cd*dsqrt(2.D0)*sw*(zp(k,2)*tw-zp(k,1))/3.D0
          ENDDO
c ---- NEUTRALINO K - SBOTTOM 1,2 - BOTTOM
      DO k=1,5
         abot(1,k)=cb*dsqrt(2.d0)*(zp(k,1)*sw/3.d0+(0.5d0
     .        -1.d0/3.d0*sw**2)*zp(k,2)/cw)-sb*scalb*zz(k,3)
         abot(2,k)=-sb*dsqrt(2.d0)*(zp(k,1)*sw/3.d0+(0.5d0
     .        -1.d0/3.d0*sw**2)*zp(k,2)/cw)-cb*scalb*zz(k,3)
         bbot(1,k)=sb*dsqrt(2.d0)*sw*(zp(k,2)*tw-zp(k,1))/3.d0
     .        -cb*scalb*zz(k,3)
         bbot(2,k)=cb*dsqrt(2.d0)*sw*(zp(k,2)*tw-zp(k,1))/3.d0
     .        +sb*scalb*zz(k,3)
      ENDDO
c ---- NEUTRALINO K - STOP 1,2 - TOP
      DO k=1,5
         atopr(1,k)=ct*dsqrt(2.d0)*(-2.d0*zp(k,1)*sw/3.d0
     .        +(-0.5d0+2.d0/3.d0*sw**2)*zp(k,2)/cw)-st*scalt*zz(k,4)
         atopr(2,k)=-st*dsqrt(2.d0)*(-2.d0*zp(k,1)*sw/3.d0+(-0.5d0
     .        +2.d0/3.d0*sw**2)*zp(k,2)/cw)-ct*scalt*zz(k,4)
         btopr(1,k)=-2.d0*st*dsqrt(2.d0)*sw*(zp(k,2)*tw-zp(k,1))/3.d0
     .        -ct*scalt*zz(k,4)
         btopr(2,k)=-2.d0*ct*dsqrt(2.d0)*sw*(zp(k,2)*tw-zp(k,1))/3.d0
     .        +st*scalt*zz(k,4)
      ENDDO
c
c --------------------------------- Gluino ------------------------ c
c
c ---- GLUINO - UP - SUP
          gur(1) = su
          gur(2) = cu
          gul(1) = -cu
          gul(2) = su
c ---- GLUINO - DOWN - SDOWN
          gdr(1) =  sd
          gdr(2) =  cd
          gdl(1) =  -cd
          gdl(2) =  sd
c ---- GLUINO - TOP - STOP
          gtr(1) = st
          gtr(2) = ct
          gtl(1) = -ct
          gtl(2) = st
c ---- GLUINO - BOTTOM - SBOTTOM
          gbr(1) = sb
          gbr(2) = cb
          gbl(1) = -cb
          gbl(2) = sb
c
c ------------------------------------ Z boson ----------------- c
c
c ---- Z - NEUTRALINO I - NEUTRALINO J
          DO I=1,5
             DO J=1,5
                onr(I,J) = 1/2.D0/cw*(zz(I,3)*zz(J,3)-zz(I,4)*zz(J,4))
                onl(I,J) = -onr(I,J)
             ENDDO
          ENDDO
c ---- Z - CHARGINO I - CHARGINO J
          opl(1,1) = -V(1,1)*V(1,1)-0.5D0*V(1,2)*V(1,2)+sw**2
          opr(1,1) = -U(1,1)*U(1,1)-0.5D0*U(1,2)*U(1,2)+sw**2
          opl(1,2) = -V(1,1)*V(2,1)-0.5D0*V(1,2)*V(2,2)
          opr(1,2) = -U(1,1)*U(2,1)-0.5D0*U(1,2)*U(2,2)        
          opl(2,1) = opl(1,2)
          opr(2,1) = opr(1,2)  
          opl(2,2) = -V(2,1)*V(2,1)-0.5D0*V(2,2)*V(2,2)+sw**2
          opr(2,2) = -U(2,1)*U(2,1)-0.5D0*U(2,2)*U(2,2)+sw**2
c ---- Z - FERMION - ANTIFERMION
          azztoptop   = 1.D0/4.D0/cw
          vzztoptop   = (1.D0-8.D0/3.D0*sw**2)/4.D0/cw
          azzbotbot   =  -1.D0/4.D0/cw
          vzzbotbot   = (-1.D0+4.D0/3.D0*sw**2)/4.D0/cw
          azztautau   =  -1.D0/4.D0/cw
          vzztautau   = (-1.D0+4.D0*sw**2)/4.D0/cw
          azzneutneut = 1.D0/4.D0/cw
          vzzneutneut = 1.D0/4.D0/cw
c ---- Z - SBOTTOM 1,2 - SBOTTOM 1,2 
          gzbb(1,1) = -cb**2+2.D0/3.D0*sw**2
          gzbb(1,2) = sb*cb
          gzbb(2,1) = sb*cb
          gzbb(2,2) = -sb**2+2.D0/3.D0*sw**2
c ---- Z - STAU 1,2 - STAU 1,2 
          gztautau(1,1) = -cl**2+2.D0*sw**2
          gztautau(1,2) = sl*cl
          gztautau(2,1) = sl*cl
          gztautau(2,2) = -sl**2+2.D0*sw**2
c          WRITE(*,*)'gztautau',gztautau(1,1),gztautau(2,2),cl,sl,sw
c ---- Z - STOP 1,2 - STOP 1,2 
          gztt(1,1) = ct**2-4.D0/3.D0*sw**2
          gztt(1,2) = -st*ct
          gztt(2,1) = -st*ct
          gztt(2,2) = st**2-4.D0/3.D0*sw**2
c
c ------------------------------------ W boson ----------------- c
c
c ---- W^+ - NEUTRALINO I - CHARGINO J
          DO I=1,5
             DO J=1,2
                or(i,j) =  1/dsqrt(2.D0)*zz(i,3)*U(j,2)+zz(i,2)*U(j,1)
                ol(i,j) = -1/dsqrt(2.D0)*zz(i,4)*V(j,2)+zz(i,2)*V(j,1)
             ENDDO
          ENDDO
c ---- W - FERMION - FERMION'
          vwff = 1/2.D0/dsqrt(2.D0)
          awff = vwff
c ---- W - STOP 1,2 - SBOTTOM 1,2 
          gwtb(1,1) = cb*ct
          gwtb(1,2) = -(sb*ct)
          gwtb(2,1) = -(cb*st)
          gwtb(2,2) = sb*st
c ---- W - SNEUTRINO_TAU 1,2 - STAU1,2 
          gwntau(1,1) = cl*cn
          gwntau(1,2) = -(sl*cn)
          gwntau(2,1) = -(cl*sn)
          gwntau(2,2) = sl*sn
c
c ------------------ Neutral Higgs bosons ----------------------- c
c
c ---- CP EVEN/ODD HIGGS - TOP - TOPBAR
      DO A=1,3
         Httr(A) = scalt*S(A,1)
      ENDDO
      DO A=1,2
         Attr(A) = -scalt*P(A,1)  
      ENDDO
c ---- CP EVEN/ODD HIGGS - BOTTOM - BOTTOMBAR
      DO A=1,3
         Hbbr(A) = scalb*S(A,2)
      ENDDO
      DO A=1,2
         Abbr(A) = -scalb*P(A,2) 
      ENDDO
C ---- CP EVEN/ODD HIGGS K - NEUTRALINO I - NEUTRALINO J
        DO K = 1,3 
           DO I = 1,5
              DO J= 1,5
                 hchichi(k,i,j) = 1/(2.0*DSQRT(G2s))*GHNEUNEU(k,j,i)
              ENDDO
           ENDDO
        ENDDO  
        DO K = 1,2 
           DO I = 1,5
              DO J= 1,5
                 achichi(k,i,j) = -1/(2.0*DSQRT(G2s))*GANEUNEU(k,j,i)
              ENDDO
           ENDDO
        ENDDO  
c ---- CP EVEN/ODD HIGGS K - CHARGINO I - CHARGINO J
            DO K=1,3
              DO I=1,2
                DO J=1,2
                   hchachaR(K,I,J) = 1/DSQRT(G2s)*GHCHACHA(K,J,I)
                   hchachaL(K,I,J) = 1/DSQRT(G2s)*GHCHACHA(K,I,J)
                ENDDO
             ENDDO
          ENDDO
c
          DO K=1,2
             DO I=1,2
                DO J=1,2
                   achachaR(k,i,j) =  -1/DSQRT(G2s)*GACHACHA(K,J,I)
                   achachaL(k,i,j) =  1/DSQRT(G2s)*GACHACHA(K,I,J)
                ENDDO
             ENDDO
          ENDDO
c ---- CP EVEN/ODD HIGGS I - STAU 1,2 - STAU 1,2
          DO I = 1,3
             Hstaustaur(I,1,1)=HTAU11(I)*amw/amz**2
             Hstaustaur(I,1,2)=HTAU12(I)*amw/amz**2
             Hstaustaur(I,2,1)=Hstaustaur(I,1,2)
             Hstaustaur(I,2,2)=HTAU22(I)*amw/amz**2
          ENDDO   
          DO I = 1,2
             Astaustaur(I,1,1)=ATAU11(I)*amw/amz**2
             Astaustaur(I,1,2)=ATAU12(I)*amw/amz**2
             Astaustaur(I,2,1)=Astaustaur(I,1,2)
             Astaustaur(I,2,2)=ATAU22(I)*amw/amz**2
          ENDDO
c ---- CP EVEN/ODD HIGGS I - SBOTTOM 1,2 - SBOTTOM 1,2
          DO I = 1,3
             Hsbotsbotr(I,1,1) = HISB11(I)*amw/amz**2
             Hsbotsbotr(I,1,2) = HISB12(I)*amw/amz**2
             Hsbotsbotr(I,2,1) = Hsbotsbotr(I,1,2)
             Hsbotsbotr(I,2,2) = HISB22(I)*amw/amz**2
          ENDDO
          DO I = 1,2
             Asbotsbotr(I,1,1) = AISB11(I)*amw/amz**2
             Asbotsbotr(I,1,2) = AISB12(I)*amw/amz**2
             Asbotsbotr(I,2,1) = Asbotsbotr(I,1,2)
             Asbotsbotr(I,2,2) = AISB22(I)*amw/amz**2
          ENDDO
c ---- CP EVEN/ODD HIGGS I - STOP 1,2 - STOP 1,2
          DO I = 1,3
             Hstopstopr(I,1,1) = HIST11(I)*amw/amz**2
             Hstopstopr(I,1,2) = HIST12(I)*amw/amz**2
             Hstopstopr(I,2,1) = Hstopstopr(I,1,2)
             Hstopstopr(I,2,2) = HIST22(I)*amw/amz**2
          ENDDO
c          
          DO I = 1,2
             Astopstopr(I,1,1) = AIST11(I)*amw/amz**2
             Astopstopr(I,1,2) = AIST12(I)*amw/amz**2
             Astopstopr(I,2,1) = Astopstopr(I,1,2)
             Astopstopr(I,2,2) = AIST22(I)*amw/amz**2
          ENDDO
C   ----DERIVATIVE COUPLINGS ::: stop1/2-stop1/2-H/A     
       DO K =1,3
      dtlttr(K,1,1)=dtHIST11(K)*amw/amz**2
      dtlttr(K,1,2)=dtHIST12(K)*amw/amz**2
      dtlttr(K,2,1)=dtlttr(K,1,2)
      dtlttr(K,2,2)=dtHIST22(K)*amw/amz**2
       ENDDO
*
       DO K =1,2
      dtattr(K,1,1) = dtAIST11(K)*amw/amz**2
      dtattr(K,1,2) = dtAIST12(K)*amw/amz**2
      dtattr(K,2,1) = dtattr(K,1,2)
      dtattr(K,2,2) = dtAIST22(K)*amw/amz**2
      ENDDO
*
      DO K=1,3
      datlttr(K,1,1)=datHIST11(K)*amw/amz**2
      datlttr(K,1,2)=datHIST12(K)*amw/amz**2
      datlttr(K,2,1)=datlttr(K,1,2)
      datlttr(K,2,2)=datHIST22(K)*amw/amz**2
      ENDDO
*
      DO K=1,2
      datattr(K,1,1)= datAIST11(K)*amw/amz**2
      datattr(K,1,2)= datAIST12(K)*amw/amz**2
      datattr(K,2,1)= datattr(k,1,2)
      datattr(K,2,2)= datAIST22(K)*amw/amz**2
      ENDDO
*
      DO K=1,3
      dthlttr(K,1,1)=dmixHIST11(K)*amw/amz**2
      dthlttr(K,1,2)=dmixHIST12(K)*amw/amz**2
      dthlttr(K,2,1)=dthlttr(K,1,2)
      dthlttr(K,2,2)=dmixHIST22(K)*amw/amz**2
      ENDDO
*
      DO K=1,2
      dthattr(K,1,1)= dmixAIST11(K)*amw/amz**2 
      dthattr(K,1,2)= dmixAIST12(K)*amw/amz**2
      dthattr(K,2,1)= dthattr(K,1,2)
      dthattr(K,2,2)= dmixAIST22(K)*amw/amz**2
      ENDDO
c ----DERIVATIVE COUPLINGS ::: sbot1/2-sbot1/2-H/A     
      DO K =1,3
      dblbbr(K,1,1)=dtHISB11(K)*amw/amz**2
      dblbbr(K,1,2)=dtHISB12(K)*amw/amz**2
      dblbbr(K,2,1)=dblbbr(K,1,2)
      dblbbr(K,2,2)=dtHISB22(K)*amw/amz**2
      ENDDO
*
      DO K=1,2
      dbabbr(K,1,1)=dtAISB11(K)*amw/amz**2
      dbabbr(K,1,2)=dtAISB12(K)*amw/amz**2
      dbabbr(K,2,1)=dbabbr(K,1,2)
      dbabbr(K,2,2)=dtAISB22(K)*amw/amz**2
      ENDDO
*
      DO K=1,3
      dablbbr(K,1,1)=dabHISB11(k)*amw/amz**2
      dablbbr(K,1,2)=dabHISB12(k)*amw/amz**2
      dablbbr(K,2,1)=dablbbr(K,1,2)
      dablbbr(K,2,2)=dabHISB22(K)*amw/amz**2
      ENDDO
*
      DO K=1,2
      dababbr(K,1,1)=dabAISB11(K)*amw/amz**2
      dababbr(K,1,2)=dabAISB12(K)*amw/amz**2
      dababbr(K,2,1)=dababbr(K,1,2)
      dababbr(K,2,2)=dabAISB22(K)*amw/amz**2
      ENDDO
*
      DO K=1,3
      dthlbbr(K,1,1)=dmixHISB11(K)*amw/amz**2
      dthlbbr(K,1,2)=dmixHISB12(K)*amw/amz**2
      dthlbbr(K,2,1)=dthlbbr(K,1,2)
      dthlbbr(K,2,2)=dmixHISB22(K)*amw/amz**2
      ENDDO
*
      DO K=1,2
      dthabbr(K,1,1)=dmixAISB11(K)*amw/amz**2
      dthabbr(K,1,2)=dmixAISB12(K)*amw/amz**2
      dthabbr(K,2,1)=dthabbr(K,1,2)
      dthabbr(K,2,2)=dmixAISB22(K)*amw/amz**2
      ENDDO
c
c ------------------ Charged Higgs boson ----------------------- c
c
c ---- CHARGED H^+ - DBAR - UP
      chtbrunr = scalb*SINBETA
      chtbrunl = scalt*COSBETA
c
c ---- CHARGED HIGGS H^+ - Tau - neutrino_tau
          chctaunur = scaltau*SINBETA
          chctaunul = 0.D0
          vchtau = -scaltau/2.D0*SINBETA
          achtau = -vchtau
c
c ---- CHARGED H^+ - TOPBAR - BOTTOM
          vchtop = -scalb/2.D0*SINBETA-scalt/2.D0*COSBETA
          achtop = scalb/2.D0*SINBETA-scalt/2.D0*COSBETA
c
c ---- CHARGED HIGGS - NEUTRALINO A - CHARGINO B
          DO A=1,5
             DO B=1,2
                qr(A,B) = GHCNEUCHAL(A,B)/DSQRT(g2s) 
                ql(A,B) = GHCNEUCHAR(A,B)/DSQRT(g2s)
             ENDDO
          ENDDO
c
c ---- CHARGED HIGGS^+ - SNEUTRINO_TAU 1,2 - STAU 1,2
      s11 = 1.d0/dsqrt(2.d0)*amw*(scaltau**2-1.D0)*
     .  2d0*SINBETA*COSBETA
      s12 = scaltau*COSBETA*(atau_yu*tanbeta+MUEFFQ)
      s21 = 0.D0
      s22 = 0.D0      
      chcsntau(1,1)=(-cn*cl*s11-sn*sl*s22-cn*sl*s12-sn*cl*s21)
      chcsntau(1,2)=(cn*sl*s11-cn*cl*s12+sl*sn*s21-sn*cl*s22)
      chcsntau(2,1)=(sn*cl*s11+sn*sl*s12-cn*cl*s21-cn*sl*s22)
      chcsntau(2,2)=(-sn*sl*s11+sn*cl*s12+cn*sl*s21-cn*cl*s22)
      gcsntaur(1,1)=chcsntau(1,1)/amw
      gcsntaur(1,2)=chcsntau(1,2)/amw
      gcsntaur(2,1)=chcsntau(2,1)/amw
      gcsntaur(2,2)=chcsntau(2,2)/amw
c
c ---- CHARGED HIGGS^+ - STOP 1,2 - SBOTTOM 1,2
      s11 = 1.d0/dsqrt(2.d0)*amw*(scalb**2+scalt**2-1.D0)*
     .  2d0*SINBETA*COSBETA
      s12 = scalb*COSBETA*(ad*tanbeta+MUEFFQ)
      s21 = scalt*SINBETA*(au/tanbeta+MUEFFQ)
      s22 = dsqrt(2.d0)*amw*scalb*scalt
      chctb(1,1)=(-ct*cb*s11-st*sb*s22-ct*sb*s12-st*cb*s21)
      chctb(1,2)=(ct*sb*s11-ct*cb*s12+sb*st*s21-st*cb*s22)
      chctb(2,1)=(st*cb*s11+st*sb*s12-ct*cb*s21-ct*sb*s22)
      chctb(2,2)=(-st*sb*s11+st*cb*s12+ct*sb*s21-ct*cb*s22)
      gctbr(1,1)=chctb(1,1)/amw
      gctbr(1,2)=chctb(1,2)/amw
      gctbr(2,1)=chctb(2,1)/amw
      gctbr(2,2)=chctb(2,2)/amw
c
c ---------------------------------------------------------- c
c                   CALL OF DECAY ROUTINES
c -----------------------------------------------------------c
      CALL NS_GLUINO
      CALL NS_STOP
      CALL NS_SBOTTOM
      CALL NS_SQUARKS
      CALL NS_SLEPTON
      CALL NS_CHARGINO
      CALL NS_NEUTRALINO
c
      END
