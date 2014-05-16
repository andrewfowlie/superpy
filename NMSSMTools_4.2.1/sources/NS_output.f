c ==================================================================== c
c                        OUTPUT (SLHA) for DECAYS                      c
c                                                                      c
c ==================================================================== c     
      SUBROUTINE NS_output

      IMPLICIT DOUBLE PRECISION (a-h,m,o-z)
      INTEGER nout
      DOUBLE PRECISION flagmulti,flagqcd,flagloop
      DOUBLE PRECISION chartot2(2),chartot(2),chartot3(2)
      DOUBLE PRECISION brcharst1(2),brcharst2(2),brcharsb1(2),
     .     brcharsb2(2),brcharsupl(2),brcharsupr(2),
     .     brcharsdownl(2),
     .     brcharsdownr(2),brcharsnel(2),brcharsn1(2),
     .     brcharsn2(2),brcharsell(2),brcharstau1(2),brcharselr(2),
     .     brcharstau2(2),brcharhcneut(2,5),brcharwneut(2,5),
     .     brcharzchic,brcharHchic(3),brcharAchic(2)
      DOUBLE PRECISION brntaunut(2,5),brnelnue(2,5),brnmunumu(2,5),
     .     brnupdb(2,5),brnchsb(2,5),brntopbb(2,5),
     .     brglupdb(2),brglchsb(2),brgltopbb(2),
     .     brchee,brchmumu,brchtautau,brchnene,
     .     brchnmunmu,brchntauntau,brchupup,brchdodo,
     .     brchchch,brchstst,brchtoptop,brchbotbot
*
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttot3(5),neuttotrad(5)
      DOUBLE PRECISION brneutst1(5),brneutst2(5),brneutsb1(5),
     .          brneutsb2(5),
     .          brneutsupl(5),brneutsupr(5),brneutsdownl(5),
     .          brneutsdownr(5),brneutsnel(5),brneutsn1(5),
     .          brneutsn2(5),brneutsell(5),brneutselr(5),
     .          brneutstau1(5),brneutstau2(5),brneutwchar(5,2),
     .          brneuthcchar(5,2),brneutzneut(5,5),
     .          brneutHneut(5,5,3),brneutAneut(5,5,2),brnraddec(5,5)
      DOUBLE PRECISION brneutup(5,5),brneutdow(5,5),brneutch(5,5),
     .          brneutst(5,5),brneutbot(5,5),brneuttop(5,5),
     .          brneutel(5,5),brneutmu(5,5),brneuttau(5,5),
     .          brneutnue(5,5),brneutnumu(5,5),brneutnutau(5,5),
     .          brchubd(5,2),brchcbs(5,2),brchtbb(5,2),brchelne(5,2),
     .          brchmunmu(5,2),brchtauntau(5,2),brglup(5),brgldo(5),
     .          brglch(5),brglst(5),brgltop(5),brglbot(5)
*
      DOUBLE PRECISION selltot,selltot2,selltot3,
     . selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     . sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .          brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .          brsntauneut(5),brsntauchar(2),brsntau1hcstau(2),
     .          brsntau1wstau(2),brstau1neut(5),brstau2neut(5),
     .          brstau1char(2),
     .          brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .          brstau1wsn(2),brstau2wsn(2),brstau2H(3),brstau2A(2),
     .          brstau2ztau
      DOUBLE PRECISION brselrstau,brselrstaustar,brstau2stau1star,
     .    brstau2stau1,brstau2stau1nn,
     .    brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .    brsellstau1star,brsellstau1,brsellstau1nutau,
     .    brsnestau1star,brsnestau1,brsnestau1nutau
*
      DOUBLE PRECISION supltot2,suprtot2,sdowltot2,sdowrtot2
      DOUBLE PRECISION brsuplnup(5),brsuplcdow(2),brsuplglui,
     .          brsuprnup(5),brsuprcdow(2),brsuprglui,
     .          brsdowlndow(5),brsdowlchup(2),brsdowlglui,
     .          brsdowrndow(5),brsdowrchup(2),brsdowrglui
*
      DOUBLE PRECISION stoptot(2),stoptot2(2),stoptot3(2),stoptotrad(2)
      DOUBLE PRECISION brst1neutt(5),brst2neutt(5),brst1charb(2),
     .          brst2charb(2),brst1hcsb(2),brst2hcsb(2),brst1wsb(2),
     .          brst2wsb(2),brst1glui,brst2glui,brst2H(3),brst2A(2),
     .          brst2ztop,brgamma,brgammaup,brgammagluino
      DOUBLE PRECISION brstopw(2,5),brstoph(2,5),brststau(2,2),
     .          brstsntau(2,2),brstsel(2,2),brstsnel(2),
     .          brstbsbst(2,2),brstbbsbt(2,2),brsttausbnu(2,2),
     .          brstelsbnu(2,2),brstupsbdow(2,2),
     .          brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .          brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      DOUBLE PRECISION sbottot(2),sbottot2(2),sbottot3(2)
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
*
      DOUBLE PRECISION gluitot,gluitot2,gluitot3,gluitotrad
      DOUBLE PRECISION brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon(5)
      DOUBLE PRECISION brgoup(5),brgoch(5),brgodn(5),brgost(5),
     .         brgotp(5),
     .         brgobt(5),brgoud(2),brgocs(2),brgotb(2),brhcst1b,brwst1b
**********************************************************
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartot3
      COMMON/CHARGINO_BR_2BD/brcharst1,brcharst2,brcharsb1,
     .          brcharsb2,brcharsupl,brcharsupr,brcharsdownl,
     .          brcharsdownr,brcharsnel,brcharsn1,
     .          brcharsn2,brcharsell,brcharstau1,brcharselr,
     .          brcharstau2,brcharhcneut,brcharwneut,
     .          brcharzchic,brcharHchic,brcharAchic
      COMMON/CHARGINO_BR_3BD/brntaunut,brnelnue,brnmunumu,
     .          brnupdb,brnchsb,brntopbb,
     .          brglupdb,brglchsb,brgltopbb,
     .          brchee,brchmumu,brchtautau,brchnene,
     .          brchnmunmu,brchntauntau,brchupup,brchdodo,
     .          brchchch,brchstst,brchtoptop,brchbotbot
*
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttot3,neuttotrad
      COMMON/NEUTRALINO_BR_2BD/brneutst1,brneutst2,brneutsb1,brneutsb2,
     .         brneutsupl,brneutsupr,brneutsdownl,brneutsdownr,
     .         brneutsnel,brneutsn1,brneutsn2,brneutsell,brneutselr,
     .         brneutstau1,brneutstau2,brneutwchar,brneuthcchar,
     .         brneutzneut,brneutHneut,brneutAneut,brnraddec
      COMMON/NEUTRALINO_BR_3BD/brneutup,brneutdow,brneutch,brneutst,
     .         brneutbot,brneuttop,brneutel,brneutmu,brneuttau,
     .         brneutnue,brneutnumu,brneutnutau,brchubd,brchcbs, 
     .         brchtbb,brchelne,brchmunmu,brchtauntau,brglup,brgldo,
     .         brglch,brglst,brgltop,brglbot
**
      COMMON/GLUINO_WIDTH/gluitot,gluitot2,gluitot3,gluitotrad
      COMMON/GLUINO_BR_2BD/brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon
      COMMON/GLUINO_BR_3BD/brgoup,brgoch,brgodn,brgost,brgotp,
     .         brgobt,brgoud,brgocs,brgotb,brhcst1b,brwst1b
**
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
     .   brsellstau1star,brsellstau1,brsellstau1nutau,
     .   brsnestau1star,brsnestau1,brsnestau1nutau
**
      COMMON/SQUARK_WIDTH/supltot2,suprtot2,sdowltot2,sdowrtot2
      COMMON/SQUARK_BR_2BD/brsuplnup,brsuplcdow,brsuplglui,
     .          brsuprnup,brsuprcdow,brsuprglui,
     .          brsdowlndow,brsdowlchup,brsdowlglui,
     .          brsdowrndow,brsdowrchup,brsdowrglui
**
      COMMON/STOP_WIDTH/stoptot,stoptot2,stoptot3,stoptotrad
      COMMON/STOP_BR_2BD/brst1neutt,brst2neutt,brst1charb,
     .          brst2charb,brst1hcsb,brst2hcsb,brst1wsb,
     .          brst2wsb,brst1glui,brst2glui,brst2H,brst2A,
     .          brst2ztop,brgamma,brgammaup,brgammagluino
      COMMON/STOP_BR_3BD/brstopw,brstoph,brststau,
     .          brstsntau,brstsel,brstsnel,
     .          brstbsbst,brstbbsbt,brsttausbnu,
     .          brstelsbnu,brstupsbdow,
     .          brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .          brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
**
      COMMON/SBOTTOM_WIDTH/sbottot,sbottot2,sbottot3
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
        
c ==================================================================== c
c                          The output file                             c
c ==================================================================== c

c ------------------ output a la Les Houches accord ------------------ c

      id =1
      idb=-1
      iu =2
      iub=-2
      is =3
      isb=-3
      ic =4
      icb=-4
      ib =5
      ibb=-5
      it =6
      itb=-6

      ie   =11
      ine  =12
      imu  =13
      inmu =14
      itau =15
      intau=16

      ihH1=25
      ihH2=35
      ihH3=45  
      ihA1=36
      ihA2=46  
      ihc=37
      igl=21
      iga=22
      iz =23
      iwc=24

      isdl=1000001
      isdr=2000001
      isul=1000002
      isur=2000002
      issl=1000003
      issr=2000003
      iscl=1000004
      iscr=2000004
      isb1=1000005
      isb2=2000005
      ist1=1000006
      ist2=2000006

      iglo=1000021
      in1 =1000022
      in2 =1000023
      in3 =1000025
      in4 =1000035
      in5 =1000045  
      ic1 =1000024
      ic2 =1000037

      intau1=1000016 
      intau2=2000016 
      inel  =1000012
      iner  =2000012
      inmul =1000014
      inmur =2000014
      
      isell =1000011
      iselr =2000011
      ismul =1000013
      ismur =2000013
      istau1=1000015
      istau2=2000015

      igrav =1000039
      if1 = 3000001
      if2 = 4000001 
      idrbar          = 1
c -- opening the output file
       nout=18
c --------- c
c W+ decays c
c --------- c
         write(nout,99)
         write(nout,100) iwc,2.085D0,'W+ (measured)'
            write(nout,102) 0.1165D0,2,-ie,ine,  
     .'BR(W+ -> e+ nu_e)'
            write(nout,102) 0.1165D0,2,-imu,inmu,  
     .'BR(W+ -> mu+ nu_mu)'
            write(nout,102) 0.112D0,2,-itau,intau,  
     .'BR(W+ -> tau+ nu_tau)'
            write(nout,102) 0.365D0,2,iu,idb,  
     .'BR(W+ -> u db)'
            write(nout,102) 0.31D0,2,ic,isb,  
     .'BR(W+ -> c sb)'
c --------- c
c Z decays c
c --------- c
         write(nout,99)
         write(nout,100) iz,2.4952D0,'Z (measured)'
            write(nout,102) 0.200D0,2,-ine,ine,  
     .'BR(Z -> invisible)'
            write(nout,102) 0.03365D0,2,-ie,ie,  
     .'BR(Z -> e+ e-)'
            write(nout,102) 0.03365D0,2,-imu,imu,  
     .'BR(Z -> mu+ mu-)'
            write(nout,102) 0.0337D0,2,-itau,itau,  
     .'BR(Z -> tau+ tau-)'
            write(nout,102) 0.111D0,2,iu,-iu,
     .'BR(Z -> u ub)'
            write(nout,102) 0.1585D0,2,id,-id,
     .'BR(Z -> d db)'
            write(nout,102) 0.1585D0,2,is,-is,
     .'BR(Z -> s sb)'
            write(nout,102) 0.12D0,2,ic,-ic,
     .'BR(Z -> c cb)'
            write(nout,102) 0.151D0,2,ib,-ib,
     .'BR(Z -> b bb)'
c ---------------- c
c chargino1 decays c
c ---------------- c
         write(nout,99)
         write(nout,100) 1000024,chartot(1),'chargino1'
      if(chartot2(1).ne.0.D0) then
      write(nout,49) 'chargino1 2-body decays'
         write(nout,101)
      if(brcharsupl(1).ne.0.D0) then
            write(nout,102) brcharsupl(1),2,isul,idb,  
     .'BR(~chi_1+ -> ~u_L  db)'
         endif
      if(brcharsupr(1).ne.0.D0) then
      write(nout,102) brcharsupr(1),2,isur,idb,  'BR(~chi_1+ -> ~u_R  
     . db)'
      endif
      if(brcharsdownl(1).ne.0.D0) then
      write(nout,102) brcharsdownl(1),2,-isdl,iu,'BR(~chi_1+ -> ~d_L* 
     . u )'
      endif
      if(brcharsdownr(1).ne.0.D0) then
      write(nout,102) brcharsdownr(1),2,-isdr,iu,'BR(~chi_1+ -> ~d_R* 
     . u )'
      endif
      if(brcharsupl(1).ne.0.D0) then
      write(nout,102) brcharsupl(1),2,iscl,isb,  'BR(~chi_1+ -> ~c_L  
     . sb)'
      endif
      if(brcharsupr(1).ne.0.D0) then
      write(nout,102) brcharsupr(1),2,iscr,isb,  'BR(~chi_1+ -> ~c_R  
     . sb)'
      endif
      if(brcharsdownl(1).ne.0.D0) then
      write(nout,102) brcharsdownl(1),2,-issl,ic,'BR(~chi_1+ -> ~s_L* 
     . c )'
      endif
      if(brcharsdownr(1).ne.0.D0) then
      write(nout,102) brcharsdownr(1),2,-issr,ic,'BR(~chi_1+ -> ~s_R* 
     . c )'
      endif
      if(brcharst1(1).ne.0.D0) then
      write(nout,102) brcharst1(1),2,ist1,ibb,   'BR(~chi_1+ -> ~t_1  
     . bb)'
      endif
      if(brcharst2(1).ne.0.D0) then
      write(nout,102) brcharst2(1),2,ist2,ibb,   'BR(~chi_1+ -> ~t_2  
     . bb)'
      endif
      if(brcharsb1(1).ne.0.D0) then
      write(nout,102) brcharsb1(1),2,-isb1,it,   'BR(~chi_1+ -> ~b_1* 
     . t )'
      endif
      if(brcharsb2(1).ne.0.D0) then
      write(nout,102) brcharsb2(1),2,-isb2,it,   'BR(~chi_1+ -> ~b_2* 
     . t )'
      endif
      if(brcharsnel(1).ne.0.D0) then
      write(nout,102) brcharsnel(1),2,inel,-ie,  'BR(~chi_1+ -> ~nu_eL 
     . e+  )'
      write(nout,102) brcharsnel(1),2,inmul,-imu,'BR(~chi_1+ -> ~nu_muL 
     . mu+ )'
      endif
      if(brcharsn1(1).ne.0.D0) then
      write(nout,102) brcharsn1(1),2,intau1,-itau,'BR(~chi_1+ -> ~nu_tau
     .1 tau+)'
      endif
      if(brcharsell(1).ne.0.D0) then
      write(nout,102) brcharsell(1),2,-isell,ine,'BR(~chi_1+ -> ~e_L+   
     . nu_e)'
      endif
      if(brcharselr(1).ne.0.D0) then
      write(nout,102) brcharselr(1),2,-iselr,ine,'BR(~chi_1+ -> ~e_R+   
     . nu_e)'
      endif
      if(brcharsell(1).ne.0.D0) then
      write(nout,102) brcharsell(1),2,-ismul,inmu,'BR(~chi_1+ -> ~mu_L+ 
     .  nu_mu)'
      endif
      if(brcharselr(1).ne.0.D0) then
      write(nout,102) brcharselr(1),2,-ismur,inmu,'BR(~chi_1+ -> ~mu_R+ 
     .  nu_mu)'
      endif
      if(brcharstau1(1).ne.0.D0) then
      write(nout,102) brcharstau1(1),2,-istau1,intau,'BR(~chi_1+ -> ~tau
     ._1+  nu_tau)'
      endif
      if(brcharstau2(1).ne.0.D0) then
      write(nout,102) brcharstau2(1),2,-istau2,intau,'BR(~chi_1+ -> ~tau
     ._2+  nu_tau)'
      endif
      if(brcharwneut(1,1).ne.0.D0) then
      write(nout,102) brcharwneut(1,1),2,in1,iwc,    'BR(~chi_1+ -> ~chi
     ._10  W+)'
      endif
      if(brcharwneut(1,2).ne.0.D0) then
      write(nout,102) brcharwneut(1,2),2,in2,iwc,    'BR(~chi_1+ -> ~chi
     ._20  W+)'
      endif
      if(brcharwneut(1,3).ne.0.D0) then
      write(nout,102) brcharwneut(1,3),2,in3,iwc,    'BR(~chi_1+ -> ~chi
     ._30  W+)'
      endif
      if(brcharwneut(1,4).ne.0.D0) then
      write(nout,102) brcharwneut(1,4),2,in4,iwc,    'BR(~chi_1+ -> ~chi
     ._40  W+)'
      endif
      if(brcharwneut(1,5).ne.0.D0) then
      write(nout,102) brcharwneut(1,5),2,in5,iwc,    'BR(~chi_1+ -> ~chi
     ._50  W+)'
      endif
      if(brcharhcneut(1,1).ne.0.D0) then
      write(nout,102) brcharhcneut(1,1),2,in1,ihc,   'BR(~chi_1+ -> ~chi
     ._10  H+)'
      endif
      if(brcharhcneut(1,2).ne.0.D0) then
      write(nout,102) brcharhcneut(1,2),2,in2,ihc,   'BR(~chi_1+ -> ~chi
     ._20  H+)'
      endif
      if(brcharhcneut(1,3).ne.0.D0) then
      write(nout,102) brcharhcneut(1,3),2,in3,ihc,   'BR(~chi_1+ -> ~chi
     ._30  H+)'
      endif
      if(brcharhcneut(1,4).ne.0.D0) then
      write(nout,102) brcharhcneut(1,4),2,in4,ihc,   'BR(~chi_1+ -> ~chi
     ._40  H+)'
      endif
      if(brcharhcneut(1,5).ne.0.D0) then
      write(nout,102) brcharhcneut(1,5),2,in5,ihc,   'BR(~chi_1+ -> ~chi
     ._50  H+)'
      endif
      endif
C   ========================
*   ---CHARGINO1 : Three body
C   ========================

      if(chartot3(1).ne.0.D0) then 
      write(nout,49) 'chargino1 3-body decays'
      write(nout,103)
      if(brnupdb(1,1).ne.0.D0) then
      write(nout,104) brnupdb(1,1),3,in1,iu,idb,     'BR(~chi_1+ -> ~chi
     ._10 u    db)'
      endif
      if(brnupdb(1,2).ne.0.D0) then
      write(nout,104) brnupdb(1,2),3,in2,iu,idb,     'BR(~chi_1+ -> ~chi
     ._20 u    db)'
      endif
      if(brnupdb(1,3).ne.0.D0) then
      write(nout,104) brnupdb(1,3),3,in3,iu,idb,     'BR(~chi_1+ -> ~chi
     ._30 u    db)'
      endif
      if(brnupdb(1,4).ne.0.D0) then
      write(nout,104) brnupdb(1,4),3,in4,iu,idb,     'BR(~chi_1+ -> ~chi
     ._40 u    db)'
      endif
      if(brnupdb(1,5).ne.0.D0) then
      write(nout,104) brnupdb(1,5),3,in5,iu,idb,     'BR(~chi_1+ -> ~chi
     ._50 u    db)'
      endif
      if(brnupdb(1,1).ne.0.D0) then
      write(nout,104) brnupdb(1,1),3,in1,ic,isb,     'BR(~chi_1+ -> ~chi
     ._10 c    sb)'
      endif
      if(brnupdb(1,2).ne.0.D0) then
      write(nout,104) brnupdb(1,2),3,in2,ic,isb,     'BR(~chi_1+ -> ~chi
     ._20 c    sb)'
      endif
      if(brnupdb(1,3).ne.0.D0) then
      write(nout,104) brnupdb(1,3),3,in3,ic,isb,     'BR(~chi_1+ -> ~chi
     ._30 c    sb)'
      endif
      if(brnupdb(1,4).ne.0.D0) then
      write(nout,104) brnupdb(1,4),3,in4,ic,isb,     'BR(~chi_1+ -> ~chi
     ._40 c    sb)'
      endif
      if(brnupdb(1,5).ne.0.D0) then
      write(nout,104) brnupdb(1,5),3,in5,ic,isb,     'BR(~chi_1+ -> ~chi
     ._50 c    sb)'
      endif
      if(brntopbb(1,1).ne.0.D0) then
      write(nout,104) brntopbb(1,1),3,in1,it,ibb,    'BR(~chi_1+ -> ~chi
     ._10 t    bb)'
      endif
      if(brntopbb(1,2).ne.0.D0) then
      write(nout,104) brntopbb(1,2),3,in2,it,ibb,    'BR(~chi_1+ -> ~chi
     ._20 t    bb)'
      endif
      if(brntopbb(1,3).ne.0.D0) then
      write(nout,104) brntopbb(1,3),3,in3,it,ibb,    'BR(~chi_1+ -> ~chi
     ._30 t    bb)'
      endif
      if(brntopbb(1,4).ne.0.D0) then
      write(nout,104) brntopbb(1,4),3,in4,it,ibb,    'BR(~chi_1+ -> ~chi
     ._40 t    bb)'
      endif
      if(brntopbb(1,5).ne.0.D0) then
      write(nout,104) brntopbb(1,5),3,in5,it,ibb,    'BR(~chi_1+ -> ~chi
     ._50 t    bb)'
      endif
      if(brnelnue(1,1).ne.0.D0) then
      write(nout,104) brnelnue(1,1),3,in1,-ie,ine,   'BR(~chi_1+ -> ~chi
     ._10 e+   nu_e)'
      endif
      if(brnelnue(1,2).ne.0.D0) then
      write(nout,104) brnelnue(1,2),3,in2,-ie,ine,   'BR(~chi_1+ -> ~chi
     ._20 e+   nu_e)'
      endif
      if(brnelnue(1,3).ne.0.D0) then
      write(nout,104) brnelnue(1,3),3,in3,-ie,ine,   'BR(~chi_1+ -> ~chi
     ._30 e+   nu_e)'
      endif
      if(brnelnue(1,4).ne.0.D0) then
      write(nout,104) brnelnue(1,4),3,in4,-ie,ine,   'BR(~chi_1+ -> ~chi
     ._40 e+   nu_e)'
      endif
      if(brnelnue(1,5).ne.0.D0) then
      write(nout,104) brnelnue(1,5),3,in5,-ie,ine,   'BR(~chi_1+ -> ~chi
     ._50 e+   nu_e)'
      endif
      if(brnmunumu(1,1).ne.0.D0) then
      write(nout,104) brnmunumu(1,1),3,in1,-imu,inmu,'BR(~chi_1+ -> ~chi
     ._10 mu+  nu_mu)'
      endif
      if(brnmunumu(1,2).ne.0.D0) then
      write(nout,104) brnmunumu(1,2),3,in2,-imu,inmu,'BR(~chi_1+ -> ~chi
     ._20 mu+  nu_mu)'
      endif
      if(brnmunumu(1,3).ne.0.D0) then
      write(nout,104) brnmunumu(1,3),3,in3,-imu,inmu,'BR(~chi_1+ -> ~chi
     ._30 mu+  nu_mu)'
      endif
      if(brnmunumu(1,4).ne.0.D0) then
      write(nout,104) brnmunumu(1,4),3,in4,-imu,inmu,'BR(~chi_1+ -> ~chi
     ._40 mu+  nu_mu)'
      endif
      if(brnmunumu(1,5).ne.0.D0) then
      write(nout,104) brnmunumu(1,5),3,in5,-imu,inmu,'BR(~chi_1+ -> ~chi
     ._50 mu+  nu_mu)'
      endif
      if(brntaunut(1,1).ne.0.D0) then
      write(nout,104) brntaunut(1,1),3,in1,-itau,intau,'BR(~chi_1+ -> ~c
     .hi_10 tau+ nu_tau)'
      endif
      if(brntaunut(1,2).ne.0.D0) then
      write(nout,104) brntaunut(1,2),3,in2,-itau,intau,'BR(~chi_1+ -> ~c
     .hi_20 tau+ nu_tau)'
      endif
      if(brntaunut(1,3).ne.0.D0) then
      write(nout,104) brntaunut(1,3),3,in3,-itau,intau,'BR(~chi_1+ -> ~c
     .hi_30 tau+ nu_tau)'
      endif
      if(brntaunut(1,4).ne.0.D0) then
      write(nout,104) brntaunut(1,4),3,in4,-itau,intau,'BR(~chi_1+ -> ~c
     .hi_40 tau+ nu_tau)'
      endif
      if(brntaunut(1,5).ne.0.D0) then
      write(nout,104) brntaunut(1,5),3,in5,-itau,intau,'BR(~chi_1+ -> ~c
     .hi_50 tau+ nu_tau)'
      endif
      if(brglupdb(1).ne.0.D0) then
      write(nout,104) brglupdb(1),3,iglo,iu,idb,       'BR(~chi_1+ -> ~g
     .      u    db)'
      endif
      if(brglchsb(1).ne.0.D0) then
      write(nout,104) brglchsb(1),3,iglo,ic,isb,       'BR(~chi_1+ -> ~g
     .      c    sb)'
      endif
      if(brgltopbb(1).ne.0.D0) then
      write(nout,104) brgltopbb(1),3,iglo,it,ibb,      'BR(~chi_1+ -> ~g
     .      t    bb)'
      endif
      endif
 
c ---------------- c
c chargino2 decays c
c ---------------- c
      write(nout,99)
      write(nout,100) 1000037,chartot(2),'chargino2'
      if(chartot2(2).ne.0.D0) then
      write(nout,49) 'chargino2 2-body decays'
      write(nout,101)
      if(brcharsupl(2).ne.0.D0) then
      write(nout,102) brcharsupl(2),2,isul,idb,  'BR(~chi_2+ -> ~u_L  
     . db)'
      endif
      if(brcharsupr(2).ne.0.D0) then
      write(nout,102) brcharsupr(2),2,isur,idb,  'BR(~chi_2+ -> ~u_R  
     . db)'
      endif
      if(brcharsdownl(2).ne.0.D0) then
      write(nout,102) brcharsdownl(2),2,-isdl,iu,'BR(~chi_2+ -> ~d_L* 
     . u )'
      endif
      if(brcharsdownr(2).ne.0.D0) then
      write(nout,102) brcharsdownr(2),2,-isdr,iu,'BR(~chi_2+ -> ~d_R* 
     . u )'
      endif
      if(brcharsupl(2).ne.0.D0) then
      write(nout,102) brcharsupl(2),2,iscl,isb,  'BR(~chi_2+ -> ~c_L  
     . sb)'
      endif
      if(brcharsupr(2).ne.0.D0) then
      write(nout,102) brcharsupr(2),2,iscr,isb,  'BR(~chi_2+ -> ~c_R  
     . sb)'
      endif
      if(brcharsdownl(2).ne.0.D0) then
      write(nout,102) brcharsdownl(2),2,-issl,ic,'BR(~chi_2+ -> ~s_L* 
     . c )'
      endif
      if(brcharsdownr(2).ne.0.D0) then
      write(nout,102) brcharsdownr(2),2,-issr,ic,'BR(~chi_2+ -> ~s_R* 
     . c )'
      endif
      if(brcharst1(2).ne.0.D0) then
      write(nout,102) brcharst1(2),2,ist1,ibb,   'BR(~chi_2+ -> ~t_1  
     . bb)'
      endif
      if(brcharst2(2).ne.0.D0) then
      write(nout,102) brcharst2(2),2,ist2,ibb,   'BR(~chi_2+ -> ~t_2  
     . bb)'
      endif
      if(brcharsb1(2).ne.0.D0) then
      write(nout,102) brcharsb1(2),2,-isb1,it,   'BR(~chi_2+ -> ~b_1* 
     . t )'
      endif
      if(brcharsb2(2).ne.0.D0) then
      write(nout,102) brcharsb2(2),2,-isb2,it,   'BR(~chi_2+ -> ~b_2* 
     . t )'
      endif
      if(brcharsnel(2).ne.0.D0) then
      write(nout,102) brcharsnel(2),2,inel,-ie,  'BR(~chi_2+ -> ~nu_eL 
     . e+  )'
      write(nout,102) brcharsnel(2),2,inmul,-imu,'BR(~chi_2+ -> ~nu_muL 
     . mu+ )'
      endif
      if(brcharsn1(2).ne.0.D0) then
      write(nout,102) brcharsn1(2),2,intau1,-itau,'BR(~chi_2+ -> ~nu_tau
     .1 tau+)'
      endif
      if(brcharsell(2).ne.0.D0) then
      write(nout,102) brcharsell(2),2,-isell,ine,'BR(~chi_2+ -> ~e_L+   
     . nu_e)'
      endif
      if(brcharselr(2).ne.0.D0) then
      write(nout,102) brcharselr(2),2,-iselr,ine,'BR(~chi_2+ -> ~e_R+   
     . nu_e)'
      endif
      if(brcharsell(2).ne.0.D0) then
      write(nout,102) brcharsell(2),2,-ismul,inmu,'BR(~chi_2+ -> ~mu_L+ 
     .  nu_mu)'
      endif
      if(brcharselr(2).ne.0.D0) then
      write(nout,102) brcharselr(2),2,-ismur,inmu,'BR(~chi_2+ -> ~mu_R+ 
     .  nu_mu)'
      endif
      if(brcharstau1(2).ne.0.D0) then
      write(nout,102) brcharstau1(2),2,-istau1,intau,'BR(~chi_2+ -> ~tau
     ._1+  nu_tau)'
      endif
      if(brcharstau2(2).ne.0.D0) then
      write(nout,102) brcharstau2(2),2,-istau2,intau,'BR(~chi_2+ -> ~tau
     ._2+  nu_tau)'
      endif
      if(brcharzchic.ne.0.D0) then
      write(nout,102) brcharzchic,2,ic1,iz          ,'BR(~chi_2+ -> ~chi
     ._1+  Z )' 
      endif
      if(brcharwneut(2,1).ne.0.D0) then
      write(nout,102) brcharwneut(2,1),2,in1,iwc,    'BR(~chi_2+ -> ~chi
     ._10  W+)'
      endif
      if(brcharwneut(2,2).ne.0.D0) then
      write(nout,102) brcharwneut(2,2),2,in2,iwc,    'BR(~chi_2+ -> ~chi
     ._20  W+)'
      endif
      if(brcharwneut(2,3).ne.0.D0) then
      write(nout,102) brcharwneut(2,3),2,in3,iwc,    'BR(~chi_2+ -> ~chi
     ._30  W+)'
      endif
      if(brcharwneut(2,4).ne.0.D0) then
      write(nout,102) brcharwneut(2,4),2,in4,iwc,    'BR(~chi_2+ -> ~chi
     ._40  W+)'
      endif
      if(brcharwneut(2,5).ne.0.D0) then
      write(nout,102) brcharwneut(2,5),2,in5,iwc,    'BR(~chi_2+ -> ~chi
     ._50  W+)'
      endif
      if(brcharHchic(1).ne.0.D0) then
      write(nout,102) brcharHchic(1),2,ic1,ihH1,     'BR(~chi_2+ -> ~chi
     ._1+  H_1 )'
      endif
      if(brcharHchic(2).ne.0.D0) then
      write(nout,102) brcharHchic(2),2,ic1,ihH2,     'BR(~chi_2+ -> ~chi
     ._1+  H_2 )'
      endif
      if(brcharHchic(3).ne.0.D0) then
      write(nout,102) brcharHchic(3),2,ic1,ihH3,     'BR(~chi_2+ -> ~chi
     ._1+  H_3 )'
      endif
      if(brcharAchic(1).ne.0.D0) then
      write(nout,102) brcharAchic(1),2,ic1,ihA1,     'BR(~chi_2+ -> ~chi
     ._1+  A_1 )'
      endif
      if(brcharAchic(2).ne.0.D0) then
      write(nout,102) brcharAchic(2),2,ic1,ihA2,     'BR(~chi_2+ -> ~chi
     ._1+  A_2 )'
      endif

      if(brcharhcneut(2,1).ne.0.D0) then
      write(nout,102) brcharhcneut(2,1),2,in1,ihc,   'BR(~chi_2+ -> ~chi
     ._10  H+)'
      endif
      if(brcharhcneut(2,2).ne.0.D0) then
      write(nout,102) brcharhcneut(2,2),2,in2,ihc,   'BR(~chi_2+ -> ~chi
     ._20  H+)'
      endif
      if(brcharhcneut(2,3).ne.0.D0) then
      write(nout,102) brcharhcneut(2,3),2,in3,ihc,   'BR(~chi_2+ -> ~chi
     ._30  H+)'
      endif
      if(brcharhcneut(2,4).ne.0.D0) then
      write(nout,102) brcharhcneut(2,4),2,in4,ihc,   'BR(~chi_2+ -> ~chi
     ._40  H+)'
      endif
      if(brcharhcneut(2,5).ne.0.D0) then
      write(nout,102) brcharhcneut(2,5),2,in5,ihc,   'BR(~chi_2+ -> ~chi
     ._50  H+)'
      endif
      endif
C    =======================
*     CHARGINO2 : Three body
C    =======================
c
      if(chartot3(2).ne.0.D0) then
      write(nout,49) 'chargino2 3-body decays'
      write(nout,103)
      if(brnupdb(2,1).ne.0.D0) then
      write(nout,104) brnupdb(2,1),3,in1,iu,idb,     'BR(~chi_2+ -> ~chi
     ._10 u    db)'
      endif
      if(brnupdb(2,2).ne.0.D0) then
      write(nout,104) brnupdb(2,2),3,in2,iu,idb,     'BR(~chi_2+ -> ~chi
     ._20 u    db)'
      endif
      if(brnupdb(2,3).ne.0.D0) then
      write(nout,104) brnupdb(2,3),3,in3,iu,idb,     'BR(~chi_2+ -> ~chi
     ._30 u    db)'
      endif
      if(brnupdb(2,4).ne.0.D0) then
      write(nout,104) brnupdb(2,4),3,in4,iu,idb,     'BR(~chi_2+ -> ~chi
     ._40 u    db)'
      endif
      if(brnupdb(2,1).ne.0.D0) then
      write(nout,104) brnupdb(2,1),3,in1,ic,isb,     'BR(~chi_2+ -> ~chi
     ._10 c    sb)'
      endif
      if(brnupdb(2,2).ne.0.D0) then
      write(nout,104) brnupdb(2,2),3,in2,ic,isb,     'BR(~chi_2+ -> ~chi
     ._20 c    sb)'
      endif
      if(brnupdb(2,3).ne.0.D0) then
      write(nout,104) brnupdb(2,3),3,in3,ic,isb,     'BR(~chi_2+ -> ~chi
     ._30 c    sb)'
      endif
      if(brnupdb(2,4).ne.0.D0) then
      write(nout,104) brnupdb(2,4),3,in4,ic,isb,     'BR(~chi_2+ -> ~chi
     ._40 c    sb)'
      endif
      if(brntopbb(2,1).ne.0.D0) then
      write(nout,104) brntopbb(2,1),3,in1,it,ibb,    'BR(~chi_2+ -> ~chi
     ._10 t    bb)'
      endif
      if(brntopbb(2,2).ne.0.D0) then
      write(nout,104) brntopbb(2,2),3,in2,it,ibb,    'BR(~chi_2+ -> ~chi
     ._20 t    bb)'
      endif
      if(brntopbb(2,3).ne.0.D0) then
      write(nout,104) brntopbb(2,3),3,in3,it,ibb,    'BR(~chi_2+ -> ~chi
     ._30 t    bb)'
      endif
      if(brntopbb(2,4).ne.0.D0) then
      write(nout,104) brntopbb(2,4),3,in4,it,ibb,    'BR(~chi_2+ -> ~chi
     ._40 t    bb)'
      endif
      if(brnelnue(2,1).ne.0.D0) then
      write(nout,104) brnelnue(2,1),3,in1,-ie,ine,   'BR(~chi_2+ -> ~chi
     ._10 e+   nu_e)'
      endif
      if(brnelnue(2,2).ne.0.D0) then
      write(nout,104) brnelnue(2,2),3,in2,-ie,ine,   'BR(~chi_2+ -> ~chi
     ._20 e+   nu_e)'
      endif
      if(brnelnue(2,3).ne.0.D0) then
      write(nout,104) brnelnue(2,3),3,in3,-ie,ine,   'BR(~chi_2+ -> ~chi
     ._30 e+   nu_e)'
      endif
      if(brnelnue(2,4).ne.0.D0) then
      write(nout,104) brnelnue(2,4),3,in4,-ie,ine,   'BR(~chi_2+ -> ~chi
     ._40 e+   nu_e)'
      endif
      if(brnmunumu(2,1).ne.0.D0) then
      write(nout,104) brnmunumu(2,1),3,in1,-imu,inmu,'BR(~chi_2+ -> ~chi
     ._10 mu+  nu_mu)'
      endif
      if(brnmunumu(2,2).ne.0.D0) then
      write(nout,104) brnmunumu(2,2),3,in2,-imu,inmu,'BR(~chi_2+ -> ~chi
     ._20 mu+  nu_mu)'
      endif
      if(brnmunumu(2,3).ne.0.D0) then
      write(nout,104) brnmunumu(2,3),3,in3,-imu,inmu,'BR(~chi_2+ -> ~chi
     ._30 mu+  nu_mu)'
      endif
      if(brnmunumu(2,4).ne.0.D0) then
      write(nout,104) brnmunumu(2,4),3,in4,-imu,inmu,'BR(~chi_2+ -> ~chi
     ._40 mu+  nu_mu)'
      endif
      if(brntaunut(2,1).ne.0.D0) then
      write(nout,104) brntaunut(2,1),3,in1,-itau,intau,'BR(~chi_2+ -> ~c
     .hi_10 tau+ nu_tau)'
      endif
      if(brntaunut(2,2).ne.0.D0) then
      write(nout,104) brntaunut(2,2),3,in2,-itau,intau,'BR(~chi_2+ -> ~c
     .hi_20 tau+ nu_tau)'
      endif
      if(brntaunut(2,3).ne.0.D0) then
      write(nout,104) brntaunut(2,3),3,in3,-itau,intau,'BR(~chi_2+ -> ~c
     .hi_30 tau+ nu_tau)'
      endif
      if(brntaunut(2,4).ne.0.D0) then
      write(nout,104) brntaunut(2,4),3,in4,-itau,intau,'BR(~chi_2+ -> ~c
     .hi_40 tau+ nu_tau)'
      endif
      if(brchupup.ne.0.D0) then
      write(nout,104) brchupup,3,ic1,iu,iub,           'BR(~chi_2+ -> ~c
     .hi_1+ u    ub)'
      endif
      if(brchdodo.ne.0.D0) then
      write(nout,104) brchdodo,3,ic1,id,idb,           'BR(~chi_2+ -> ~c
     .hi_1+ d    db)'
      endif
      if(brchchch.ne.0.D0) then
      write(nout,104) brchchch,3,ic1,ic,icb,           'BR(~chi_2+ -> ~c
     .hi_1+ c    cb)'
      endif
      if(brchstst.ne.0.D0) then
      write(nout,104) brchstst,3,ic1,is,isb,           'BR(~chi_2+ -> ~c
     .hi_1+ s    sb)'
      endif
      if(brchtoptop.ne.0.D0) then
      write(nout,104) brchtoptop,3,ic1,it,itb,         'BR(~chi_2+ -> ~c
     .hi_1+ t    tb)'
      endif
      if(brchbotbot.ne.0.D0) then
      write(nout,104) brchbotbot,3,ic1,ib,ibb,         'BR(~chi_2+ -> ~c
     .hi_1+ b    bb)'
      endif
      if(brchee.ne.0.D0) then
      write(nout,104) brchee,3,ic1,-ie,ie,             'BR(~chi_2+ -> ~c
     .hi_1+ e+   e-)'
      endif
      if(brchmumu.ne.0.D0) then
      write(nout,104) brchmumu,3,ic1,-imu,imu,         'BR(~chi_2+ -> ~c
     .hi_1+ mu+  mu-)'
      endif
      if(brchtautau.ne.0.D0) then
      write(nout,104) brchtautau,3,ic1,-itau,itau,     'BR(~chi_2+ -> ~c
     .hi_1+ tau+ tau-)'
      endif
      if(brchnene.ne.0.D0) then
      write(nout,104) brchnene,3,ic1,ine,-ine,         'BR(~chi_2+ -> ~c
     .hi_1+ nu_e   nu_eb)'
      endif
      if(brchnmunmu.ne.0.D0) then
      write(nout,104) brchnmunmu,3,ic1,inmu,-inmu,     'BR(~chi_2+ -> ~c
     .hi_1+ nu_mu  nu_mub)'
      endif
      if(brchntauntau.ne.0.D0) then
      write(nout,104) brchntauntau,3,ic1,intau,-intau, 'BR(~chi_2+ -> ~c
     .hi_1+ nu_tau nu_taub)'
      endif
      if(brglupdb(2).ne.0.D0) then
      write(nout,104) brglupdb(2),3,iglo,iu,idb,       'BR(~chi_2+ -> ~g
     .      u    db)'
      endif
      if(brglchsb(2).ne.0.D0) then
      write(nout,104) brglchsb(2),3,iglo,ic,isb,       'BR(~chi_2+ -> ~g
     .      c    sb)'
      endif
      if(brgltopbb(2).ne.0.D0) then
      write(nout,104) brgltopbb(2),3,iglo,it,ibb,      'BR(~chi_2+ -> ~g
     .      t    bb)'
      endif
      endif
c ------------------ c
c neutralino1 decays c
c ------------------ c
         write(nout,99)
         write(nout,100) 1000022,neuttot(1),'neutralino1'
      if(neuttot2(1).ne.0.D0) then
         write(nout,49) 'neutralino1 2-body decays'
         write(nout,101)
      if(brneutwchar(1,1).ne.0.D0) then
      write(nout,102) brneutwchar(1,1),2,ic1,-iwc,     'BR(~chi_10 -> ~c
     .hi_1+   W-)'
      write(nout,102) brneutwchar(1,1),2,-ic1,iwc,     'BR(~chi_10 -> ~c
     .hi_1-   W+)'
      endif
      if(brneutwchar(1,2).ne.0.D0) then
      write(nout,102) brneutwchar(1,2),2,ic2,-iwc,     'BR(~chi_10 -> ~c
     .hi_2+   W-)'
      write(nout,102) brneutwchar(1,2),2,-ic2,iwc,     'BR(~chi_10 -> ~c
     .hi_2-   W+)'
      endif
      if(brneuthcchar(1,1).ne.0.D0) then
      write(nout,102) brneuthcchar(1,1),2,ic1,-ihc,    'BR(~chi_10 -> ~c
     .hi_1+   H-)'
      write(nout,102) brneuthcchar(1,1),2,-ic1,ihc,    'BR(~chi_10 -> ~c
     .hi_1-   H+)'
      endif
      if(brneuthcchar(1,2).ne.0.D0) then
      write(nout,102) brneuthcchar(1,2),2,ic2,-ihc,    'BR(~chi_10 -> ~c
     .hi_2+   H-)'
      write(nout,102) brneuthcchar(1,2),2,-ic2,ihc,    'BR(~chi_10 -> ~c
     .hi_2-   H+)'
      endif
      if(brneutsupl(1).ne.0.D0) then
      write(nout,102) brneutsupl(1),2,isul,iub,        'BR(~chi_10 -> ~u
     ._L      ub)'
      write(nout,102) brneutsupl(1),2,-isul,iu,        'BR(~chi_10 -> ~u
     ._L*     u )'
      endif
      if(brneutsupr(1).ne.0.D0) then
      write(nout,102) brneutsupr(1),2,isur,iub,        'BR(~chi_10 -> ~u
     ._R      ub)'
      write(nout,102) brneutsupr(1),2,-isur,iu,        'BR(~chi_10 -> ~u
     ._R*     u )'
      endif
      if(brneutsdownl(1).ne.0.D0) then
      write(nout,102) brneutsdownl(1),2,isdl,idb,      'BR(~chi_10 -> ~d
     ._L      db)'
      write(nout,102) brneutsdownl(1),2,-isdl,id,      'BR(~chi_10 -> ~d
     ._L*     d )'
      endif
      if(brneutsdownr(1).ne.0.D0) then
      write(nout,102) brneutsdownr(1),2,isdr,idb,      'BR(~chi_10 -> ~d
     ._R      db)'
      write(nout,102) brneutsdownr(1),2,-isdr,id,      'BR(~chi_10 -> ~d
     ._R*     d )'
      endif
      if(brneutsupl(1).ne.0.D0) then
      write(nout,102) brneutsupl(1),2,iscl,icb,        'BR(~chi_10 -> ~c
     ._L      cb)'
      write(nout,102) brneutsupl(1),2,-iscl,ic,        'BR(~chi_10 -> ~c
     ._L*     c )'
      endif
      if(brneutsupr(1).ne.0.D0) then
      write(nout,102) brneutsupr(1),2,iscr,icb,        'BR(~chi_10 -> ~c
     ._R      cb)'
      write(nout,102) brneutsupr(1),2,-iscr,ic,        'BR(~chi_10 -> ~c
     ._R*     c )'
      endif
      if(brneutsdownl(1).ne.0.D0) then
      write(nout,102) brneutsdownl(1),2,issl,isb,      'BR(~chi_10 -> ~s
     ._L      sb)'
      write(nout,102) brneutsdownl(1),2,-issl,is,      'BR(~chi_10 -> ~s
     ._L*     s )'
      endif
      if(brneutsdownr(1).ne.0.D0) then
      write(nout,102) brneutsdownr(1),2,issr,isb,      'BR(~chi_10 -> ~s
     ._R      sb)'
      write(nout,102) brneutsdownr(1),2,-issr,is,      'BR(~chi_10 -> ~s
     ._R*     s )'
      endif
      if(brneutst1(1).ne.0.D0) then
      write(nout,102) brneutst1(1),2,ist1,itb,         'BR(~chi_10 -> ~t
     ._1      tb)'
      write(nout,102) brneutst1(1),2,-ist1,it,         'BR(~chi_10 -> ~t
     ._1*     t )'
      endif
      if(brneutst2(1).ne.0.D0) then
      write(nout,102) brneutst2(1),2,ist2,itb,         'BR(~chi_10 -> ~t
     ._2      tb)'
      write(nout,102) brneutst2(1),2,-ist2,it,         'BR(~chi_10 -> ~t
     ._2*     t )'
      endif
      if(brneutsb1(1).ne.0.D0) then
      write(nout,102) brneutsb1(1),2,isb1,ibb,         'BR(~chi_10 -> ~b
     ._1      bb)'
      write(nout,102) brneutsb1(1),2,-isb1,ib,         'BR(~chi_10 -> ~b
     ._1*     b )'
      endif
      if(brneutsb2(1).ne.0.D0) then
      write(nout,102) brneutsb2(1),2,isb2,ibb,         'BR(~chi_10 -> ~b
     ._2      bb)'
      write(nout,102) brneutsb2(1),2,-isb2,ib,         'BR(~chi_10 -> ~b
     ._2*     b )'
      endif
      if(brneutsell(1).ne.0.D0) then
      write(nout,102) brneutsell(1),2,isell,-ie,       'BR(~chi_10 -> ~e
     ._L-     e+)'
      write(nout,102) brneutsell(1),2,-isell,ie,       'BR(~chi_10 -> ~e
     ._L+     e-)'
      endif
      if(brneutselr(1).ne.0.D0) then
      write(nout,102) brneutselr(1),2,iselr,-ie,       'BR(~chi_10 -> ~e
     ._R-     e+)'
      write(nout,102) brneutselr(1),2,-iselr,ie,       'BR(~chi_10 -> ~e
     ._R+     e-)'
      endif
      if(brneutsell(1).ne.0.D0) then
      write(nout,102) brneutsell(1),2,ismul,-imu,      'BR(~chi_10 -> ~m
     .u_L-    mu+)' 
      write(nout,102) brneutsell(1),2,-ismul,imu,      'BR(~chi_10 -> ~m
     .u_L+    mu-)' 
      endif
      if(brneutselr(1).ne.0.D0) then
      write(nout,102) brneutselr(1),2,ismur,-imu,      'BR(~chi_10 -> ~m
     .u_R-    mu+)' 
      write(nout,102) brneutselr(1),2,-ismur,imu,      'BR(~chi_10 -> ~m
     .u_R+    mu-)'
      endif
      if(brneutstau1(1).ne.0.D0) then
      write(nout,102) brneutstau1(1),2,istau1,-itau,   'BR(~chi_10 -> ~t
     .au_1-   tau+)'
      write(nout,102) brneutstau1(1),2,-istau1,itau,   'BR(~chi_10 -> ~t
     .au_1+   tau-)'
      endif
      if(brneutstau2(1).ne.0.D0) then
      write(nout,102) brneutstau2(1),2,istau2,-itau,   'BR(~chi_10 -> ~t
     .au_2-   tau+)'
      write(nout,102) brneutstau2(1),2,-istau2,itau,   'BR(~chi_10 -> ~t
     .au_2+   tau-)'
      endif
      if(brneutsnel(1).ne.0.D0) then
      write(nout,102) brneutsnel(1),2,inel,-ine,       'BR(~chi_10 -> ~n
     .u_eL    nu_eb)'
      write(nout,102) brneutsnel(1),2,-inel,ine,       'BR(~chi_10 -> ~n
     .u_eL*   nu_e )'
      write(nout,102) brneutsnel(1),2,inmul,-inmu,     'BR(~chi_10 -> ~n
     .u_muL   nu_mub)'
      write(nout,102) brneutsnel(1),2,-inmul,inmu,     'BR(~chi_10 -> ~n
     .u_muL*  nu_mu )'
      endif
      if(brneutsn1(1).ne.0.D0) then
      write(nout,102) brneutsn1(1),2,intau1,-intau,    'BR(~chi_10 -> ~n
     .u_tau1  nu_taub)'
      write(nout,102) brneutsn1(1),2,-intau1,intau,    'BR(~chi_10 -> ~n
     .u_tau1* nu_tau )'
      endif
      endif

C     ========================
*     NEUTRALINO1 : Three body
C     ========================

      if(neuttot3(1).ne.0.D0) then 
      write(nout,49) 'neutralino1 3-body decays'
      write(nout,103)
      if(brchubd(1,1).ne.0.D0) then
      write(nout,104) brchubd(1,1),3,ic1,iub,id,       'BR(~chi_10 -> ~c
     .hi_1+ ub      d)'
      write(nout,104) brchubd(1,1),3,-ic1,idb,iu,      'BR(~chi_10 -> ~c
     .hi_1- db      u)'
      endif
      if(brchubd(1,2).ne.0.D0) then      
      write(nout,104) brchubd(1,2),3,ic2,iub,id,       'BR(~chi_10 -> ~c
     .hi_2+ ub      d)'
      write(nout,104) brchubd(1,2),3,-ic2,idb,iu,      'BR(~chi_10 -> ~c
     .hi_2- db      u)'
      endif
      if(brchcbs(1,1).ne.0.D0) then
      write(nout,104) brchcbs(1,1),3,ic1,icb,is,       'BR(~chi_10 -> ~c
     .hi_1+ cb      s)'
      write(nout,104) brchcbs(1,1),3,-ic1,isb,ic,      'BR(~chi_10 -> ~c
     .hi_1- sb      c)'
      endif
      if(brchcbs(1,2).ne.0.D0) then
      write(nout,104) brchcbs(1,2),3,ic2,icb,is,       'BR(~chi_10 -> ~c
     .hi_2+ cb      s)'
      write(nout,104) brchcbs(1,2),3,-ic2,isb,ic,      'BR(~chi_10 -> ~c
     .hi_2- sb      c)'
      endif
      if(brchtbb(1,1).ne.0.D0) then
      write(nout,104) brchtbb(1,1),3,ic1,itb,ib,       'BR(~chi_10 -> ~c
     .hi_1+ tb      b)'
      write(nout,104) brchtbb(1,1),3,-ic1,ibb,it,      'BR(~chi_10 -> ~c
     .hi_1- bb      t)'
      endif
      if(brchtbb(1,2).ne.0.D0) then
      write(nout,104) brchtbb(1,2),3,ic2,itb,ib,       'BR(~chi_10 -> ~c
     .hi_2+ tb      b)'
      write(nout,104) brchtbb(1,2),3,-ic2,ibb,it,      'BR(~chi_10 -> ~c
     .hi_2- bb      t)'
      endif
      if(brchelne(1,1).ne.0.D0) then
      write(nout,104) brchelne(1,1),3,ic1,-ine,ie,     'BR(~chi_10 -> ~c
     .hi_1+ nu_eb   e-)'
      write(nout,104) brchelne(1,1),3,-ic1,ine,-ie,    'BR(~chi_10 -> ~c
     .hi_1- nu_e    e+)'
      endif
      if(brchelne(1,2).ne.0.D0) then
      write(nout,104) brchelne(1,2),3,ic2,-ine,ie,     'BR(~chi_10 -> ~c
     .hi_2+ nu_eb   e-)'
      write(nout,104) brchelne(1,2),3,-ic2,ine,-ie,    'BR(~chi_10 -> ~c
     .hi_2- nu_e    e+)'
      endif
      if(brchmunmu(1,1).ne.0.D0) then
      write(nout,104) brchmunmu(1,1),3,ic1,-inmu,imu,  'BR(~chi_10 -> ~c
     .hi_1+ nu_mub  mu-)'
      write(nout,104) brchmunmu(1,1),3,-ic1,inmu,-imu, 'BR(~chi_10 -> ~c
     .hi_1- nu_mu   mu+)'
      endif
      if(brchmunmu(1,2).ne.0.D0) then
      write(nout,104) brchmunmu(1,2),3,ic2,-inmu,imu,  'BR(~chi_10 -> ~c
     .hi_2+ nu_mub  mu-)'
      write(nout,104) brchmunmu(1,2),3,-ic2,inmu,-imu, 'BR(~chi_10 -> ~c
     .hi_2- nu_mu   mu+)'
      endif
      if(brchtauntau(1,1).ne.0.D0) then
      write(nout,104) brchtauntau(1,1),3,ic1,-intau,itau, 'BR(~chi_10 ->
     . ~chi_1+ nu_taub tau-)'
      write(nout,104) brchtauntau(1,1),3,-ic1,intau,-itau,'BR(~chi_10 ->
     . ~chi_1- nu_tau  tau+)'
      endif
      if(brchtauntau(1,2).ne.0.D0) then
      write(nout,104) brchtauntau(1,2),3,ic2,-intau,itau, 'BR(~chi_10 ->
     . ~chi_2+ nu_taub tau-)'
      write(nout,104) brchtauntau(1,2),3,-ic2,intau,-itau,'BR(~chi_10 ->
     . ~chi_2- nu_tau  tau+)'
      endif
      if(brglup(1).ne.0.D0) then
      write(nout,104) brglup(1),3,iglo,iub,iu,            'BR(~chi_10 ->
     . ~g      ub      u)'
      endif
      if(brgldo(1).ne.0.D0) then
      write(nout,104) brgldo(1),3,iglo,idb,id,            'BR(~chi_10 ->
     . ~g      db      d)'
      endif
      if(brglch(1).ne.0.D0) then
      write(nout,104) brglch(1),3,iglo,icb,ic,            'BR(~chi_10 ->
     . ~g      cb      c)'
      endif
      if(brglst(1).ne.0.D0) then
      write(nout,104) brglst(1),3,iglo,isb,is,            'BR(~chi_10 ->
     . ~g      sb      s)'
      endif
      if(brgltop(1).ne.0.D0) then
      write(nout,104) brgltop(1),3,iglo,itb,it,           'BR(~chi_10 ->
     . ~g      tb      t)'
      endif
      if(brglbot(1).ne.0.D0) then
      write(nout,104) brglbot(1),3,iglo,ibb,ib,           'BR(~chi_10 ->
     . ~g      bb      b)'
      endif
      endif
c ------------------ c
c neutralino2 decays c
c ------------------ c

      write(nout,99)
      write(nout,100) 1000023,neuttot(2),'neutralino2'
      if(neuttot2(2).ne.0.D0) then
         write(nout,49) 'neutralino2 2-body decays'
         write(nout,101)
      if(brneutzneut(2,1).ne.0.D0) then
      write(nout,102) brneutzneut(2,1),2,in1,iz,       'BR(~chi_20 -> ~c
     .hi_10   Z )'
      endif
      if(brneutwchar(2,1).ne.0.D0) then
      write(nout,102) brneutwchar(2,1),2,ic1,-iwc,     'BR(~chi_20 -> ~c
     .hi_1+   W-)'
      write(nout,102) brneutwchar(2,1),2,-ic1,iwc,     'BR(~chi_20 -> ~c
     .hi_1-   W+)'
      endif
      if(brneutwchar(2,2).ne.0.D0) then
      write(nout,102) brneutwchar(2,2),2,ic2,-iwc,     'BR(~chi_20 -> ~c
     .hi_2+   W-)'
      write(nout,102) brneutwchar(2,2),2,-ic2,iwc,     'BR(~chi_20 -> ~c
     .hi_2-   W+)'
      endif
      if(brneutHneut(2,1,1).ne.0.D0) then
      write(nout,102) brneutHneut(2,1,1),2,in1,ihH1,   'BR(~chi_20 -> ~c
     .hi_10   H_1 )'
      endif
      if(brneutHneut(2,1,2).ne.0.D0) then
      write(nout,102) brneutHneut(2,1,2),2,in1,ihH2,   'BR(~chi_20 -> ~c
     .hi_10   H_2 )'
      endif
      if(brneutHneut(2,1,3).ne.0.D0) then
      write(nout,102) brneutHneut(2,1,3),2,in1,ihH3,   'BR(~chi_20 -> ~c
     .hi_10   H_3 )'
      endif
      if(brneutAneut(2,1,1).ne.0.D0) then
      write(nout,102) brneutAneut(2,1,1),2,in1,ihA1,   'BR(~chi_20 -> ~c
     .hi_10   A_1 )'
      endif
      if(brneutAneut(2,1,2).ne.0.D0) then
      write(nout,102) brneutAneut(2,1,2),2,in1,ihA2,   'BR(~chi_20 -> ~c
     .hi_10   A_2 )'
      endif
      if(brneuthcchar(2,1).ne.0.D0) then
      write(nout,102) brneuthcchar(2,1),2,ic1,-ihc,    'BR(~chi_20 -> ~c
     .hi_1+   H-)'
      write(nout,102) brneuthcchar(2,1),2,-ic1,ihc,    'BR(~chi_20 -> ~c
     .hi_1-   H+)'
      endif
      if(brneuthcchar(2,2).ne.0.D0) then
      write(nout,102) brneuthcchar(2,2),2,ic2,-ihc,    'BR(~chi_20 -> ~c
     .hi_2+   H-)'
      write(nout,102) brneuthcchar(2,2),2,-ic2,ihc,    'BR(~chi_20 -> ~c
     .hi_2-   H+)'
      endif
      if(brneutsupl(2).ne.0.D0) then
      write(nout,102) brneutsupl(2),2,isul,iub,        'BR(~chi_20 -> ~u
     ._L      ub)'
      write(nout,102) brneutsupl(2),2,-isul,iu,        'BR(~chi_20 -> ~u
     ._L*     u )'
      endif
      if(brneutsupr(2).ne.0.D0) then
      write(nout,102) brneutsupr(2),2,isur,iub,        'BR(~chi_20 -> ~u
     ._R      ub)'
      write(nout,102) brneutsupr(2),2,-isur,iu,        'BR(~chi_20 -> ~u
     ._R*     u )'
      endif
      if(brneutsdownl(2).ne.0.D0) then
      write(nout,102) brneutsdownl(2),2,isdl,idb,      'BR(~chi_20 -> ~d
     ._L      db)'
      write(nout,102) brneutsdownl(2),2,-isdl,id,      'BR(~chi_20 -> ~d
     ._L*     d )'
      endif
      if(brneutsdownr(2).ne.0.D0) then
      write(nout,102) brneutsdownr(2),2,isdr,idb,      'BR(~chi_20 -> ~d
     ._R      db)'
      write(nout,102) brneutsdownr(2),2,-isdr,id,      'BR(~chi_20 -> ~d
     ._R*     d )'
      endif
      if(brneutsupl(2).ne.0.D0) then
      write(nout,102) brneutsupl(2),2,iscl,icb,        'BR(~chi_20 -> ~c
     ._L      cb)'
      write(nout,102) brneutsupl(2),2,-iscl,ic,        'BR(~chi_20 -> ~c
     ._L*     c )'
      endif
      if(brneutsupr(2).ne.0.D0) then
      write(nout,102) brneutsupr(2),2,iscr,icb,        'BR(~chi_20 -> ~c
     ._R      cb)'
      write(nout,102) brneutsupr(2),2,-iscr,ic,        'BR(~chi_20 -> ~c
     ._R*     c )'
      endif
      if(brneutsdownl(2).ne.0.D0) then
      write(nout,102) brneutsdownl(2),2,issl,isb,      'BR(~chi_20 -> ~s
     ._L      sb)'
      write(nout,102) brneutsdownl(2),2,-issl,is,      'BR(~chi_20 -> ~s
     ._L*     s )'
      endif
      if(brneutsdownr(2).ne.0.D0) then
      write(nout,102) brneutsdownr(2),2,issr,isb,      'BR(~chi_20 -> ~s
     ._R      sb)'
      write(nout,102) brneutsdownr(2),2,-issr,is,      'BR(~chi_20 -> ~s
     ._R*     s )'
      endif
      if(brneutst1(2).ne.0.D0) then
      write(nout,102) brneutst1(2),2,ist1,itb,         'BR(~chi_20 -> ~t
     ._1      tb)'
      write(nout,102) brneutst1(2),2,-ist1,it,         'BR(~chi_20 -> ~t
     ._1*     t )'
      endif
      if(brneutst2(2).ne.0.D0) then
      write(nout,102) brneutst2(2),2,ist2,itb,         'BR(~chi_20 -> ~t
     ._2      tb)'
      write(nout,102) brneutst2(2),2,-ist2,it,         'BR(~chi_20 -> ~t
     ._2*     t )'
      endif
      if(brneutsb1(2).ne.0.D0) then
      write(nout,102) brneutsb1(2),2,isb1,ibb,         'BR(~chi_20 -> ~b
     ._1      bb)'
      write(nout,102) brneutsb1(2),2,-isb1,ib,         'BR(~chi_20 -> ~b
     ._1*     b )'
      endif
      if(brneutsb2(2).ne.0.D0) then
      write(nout,102) brneutsb2(2),2,isb2,ibb,         'BR(~chi_20 -> ~b
     ._2      bb)'
      write(nout,102) brneutsb2(2),2,-isb2,ib,         'BR(~chi_20 -> ~b
     ._2*     b )'
      endif
      if(brneutsell(2).ne.0.D0) then
      write(nout,102) brneutsell(2),2,isell,-ie,       'BR(~chi_20 -> ~e
     ._L-     e+)'
      write(nout,102) brneutsell(2),2,-isell,ie,       'BR(~chi_20 -> ~e
     ._L+     e-)'
      endif
      if(brneutselr(2).ne.0.D0) then
      write(nout,102) brneutselr(2),2,iselr,-ie,       'BR(~chi_20 -> ~e
     ._R-     e+)'
      write(nout,102) brneutselr(2),2,-iselr,ie,       'BR(~chi_20 -> ~e
     ._R+     e-)'
      endif
      if(brneutsell(2).ne.0.D0) then
      write(nout,102) brneutsell(2),2,ismul,-imu,      'BR(~chi_20 -> ~m
     .u_L-    mu+)' 
      write(nout,102) brneutsell(2),2,-ismul,imu,      'BR(~chi_20 -> ~m
     .u_L+    mu-)' 
      endif
      if(brneutselr(2).ne.0.D0) then
      write(nout,102) brneutselr(2),2,ismur,-imu,      'BR(~chi_20 -> ~m
     .u_R-    mu+)' 
      write(nout,102) brneutselr(2),2,-ismur,imu,      'BR(~chi_20 -> ~m
     .u_R+    mu-)'
      endif
      if(brneutstau1(2).ne.0.D0) then
      write(nout,102) brneutstau1(2),2,istau1,-itau,   'BR(~chi_20 -> ~t
     .au_1-   tau+)'
      write(nout,102) brneutstau1(2),2,-istau1,itau,   'BR(~chi_20 -> ~t
     .au_1+   tau-)'
      endif
      if(brneutstau2(2).ne.0.D0) then
      write(nout,102) brneutstau2(2),2,istau2,-itau,   'BR(~chi_20 -> ~t
     .au_2-   tau+)'
      write(nout,102) brneutstau2(2),2,-istau2,itau,   'BR(~chi_20 -> ~t
     .au_2+   tau-)'
      endif
      if(brneutsnel(2).ne.0.D0) then
      write(nout,102) brneutsnel(2),2,inel,-ine,       'BR(~chi_20 -> ~n
     .u_eL    nu_eb)'
      write(nout,102) brneutsnel(2),2,-inel,ine,       'BR(~chi_20 -> ~n
     .u_eL*   nu_e )'
      write(nout,102) brneutsnel(2),2,inmul,-inmu,     'BR(~chi_20 -> ~n
     .u_muL   nu_mub)'
      write(nout,102) brneutsnel(2),2,-inmul,inmu,     'BR(~chi_20 -> ~n
     .u_muL*  nu_mu )'
      endif
      if(brneutsn1(2).ne.0.D0) then
      write(nout,102) brneutsn1(2),2,intau1,-intau,    'BR(~chi_20 -> ~n
     .u_tau1  nu_taub)'
      write(nout,102) brneutsn1(2),2,-intau1,intau,    'BR(~chi_20 -> ~n
     .u_tau1* nu_tau )'
      endif
      endif
C     =============================
*     NEUTRALINO2 : Radiative decay
C     =============================
            if(flagloop.eq.1.D0) then
       if(brnraddec(2,1).ne.0.D0) then
         write(nout,102) brnraddec(2,1),2,in1,iga,  
     .   'BR(~chi_20 -> ~chi_10 gam)'
       endif
            endif
C     ========================
*     NEUTRALINO2 : Three body
C     ========================

      if(neuttot3(2).ne.0.D0) then
      write(nout,49) 'neutralino2 3-body decays'
      write(nout,103)
      if(brneutup(2,1).ne.0.D0) then
      write(nout,104) brneutup(2,1),3,in1,iub,iu,      
     .'BR(~chi_20 -> ~chi_10 ub      u)'
      endif
      if(brneutdow(2,1).ne.0.D0) then
      write(nout,104) brneutdow(2,1),3,in1,idb,id,     'BR(~chi_20 -> ~c
     .hi_10 db      d)'
      endif
      if(brneutch(2,1).ne.0.D0) then
      write(nout,104) brneutch(2,1),3,in1,icb,ic,      'BR(~chi_20 -> ~c
     .hi_10 cb      c)'
      endif
      if(brneutst(2,1).ne.0.D0) then
      write(nout,104) brneutst(2,1),3,in1,isb,is,      'BR(~chi_20 -> ~c
     .hi_10 sb      s)'
      endif
      if(brneuttop(2,1).ne.0.D0) then
      write(nout,104) brneuttop(2,1),3,in1,itb,it,     'BR(~chi_20 -> ~c
     .hi_10 tb      t)'
      endif
      if(brneutbot(2,1).ne.0.D0) then
      write(nout,104) brneutbot(2,1),3,in1,ibb,ib,     'BR(~chi_20 -> ~c
     .hi_10 bb      b)'
      endif
      if(brneutel(2,1).ne.0.D0) then
      write(nout,104) brneutel(2,1),3,in1,-ie,ie,      'BR(~chi_20 -> ~c
     .hi_10 e+      e-)'
      endif
      if(brneutmu(2,1).ne.0.D0) then
      write(nout,104) brneutmu(2,1),3,in1,-imu,imu,    'BR(~chi_20 -> ~c
     .hi_10 mu+     mu-)'
      endif
      if(brneuttau(2,1).ne.0.D0) then
      write(nout,104) brneuttau(2,1),3,in1,-itau,itau, 'BR(~chi_20 -> ~c
     .hi_10 tau+    tau-)'
      endif
      if(brneutnue(2,1).ne.0.D0) then
      write(nout,104) brneutnue(2,1),3,in1,-ine,ine,   'BR(~chi_20 -> ~c
     .hi_10 nu_eb   nu_e)'
      endif
      if(brneutnumu(2,1).ne.0.D0) then
      write(nout,104) brneutnumu(2,1),3,in1,-inmu,inmu,'BR(~chi_20 -> ~c
     .hi_10 nu_mub  nu_mu)'
      endif
      if(brneutnutau(2,1).ne.0.D0) then
      write(nout,104) brneutnutau(2,1),3,in1,-intau,intau,'BR(~chi_20 ->
     . ~chi_10 nu_taub nu_tau)'
      endif
      if(brchubd(2,1).ne.0.D0) then
      write(nout,104) brchubd(2,1),3,ic1,iub,id,       'BR(~chi_20 -> ~c
     .hi_1+ ub      d)'
      write(nout,104) brchubd(2,1),3,-ic1,idb,iu,      'BR(~chi_20 -> ~c
     .hi_1- db      u)'
      endif
      if(brchubd(2,2).ne.0.D0) then
      write(nout,104) brchubd(2,2),3,ic2,iub,id,       'BR(~chi_20 -> ~c
     .hi_2+ ub      d)'
      write(nout,104) brchubd(2,2),3,-ic2,idb,iu,      'BR(~chi_20 -> ~c
     .hi_2- db      u)'
      endif
      if(brchcbs(2,1).ne.0.D0) then
      write(nout,104) brchcbs(2,1),3,ic1,icb,is,       'BR(~chi_20 -> ~c
     .hi_1+ cb      s)'
      write(nout,104) brchcbs(2,1),3,-ic1,isb,ic,      'BR(~chi_20 -> ~c
     .hi_1- sb      c)'
      endif
      if(brchcbs(2,2).ne.0.D0) then
      write(nout,104) brchcbs(2,2),3,ic2,icb,is,       'BR(~chi_20 -> ~c
     .hi_2+ cb      s)'
      write(nout,104) brchcbs(2,2),3,-ic2,isb,ic,      'BR(~chi_20 -> ~c
     .hi_2- sb      c)'
      endif
      if(brchtbb(2,1).ne.0.D0) then
      write(nout,104) brchtbb(2,1),3,ic1,itb,ib,       'BR(~chi_20 -> ~c
     .hi_1+ tb      b)'
      write(nout,104) brchtbb(2,1),3,-ic1,ibb,it,      'BR(~chi_20 -> ~c
     .hi_1- bb      t)'
      endif
      if(brchtbb(2,2).ne.0.D0) then
      write(nout,104) brchtbb(2,2),3,ic2,itb,ib,       'BR(~chi_20 -> ~c
     .hi_2+ tb      b)'
      write(nout,104) brchtbb(2,2),3,-ic2,ibb,it,      'BR(~chi_20 -> ~c
     .hi_2- bb      t)'
      endif
      if(brchelne(2,1).ne.0.D0) then
      write(nout,104) brchelne(2,1),3,ic1,-ine,ie,     'BR(~chi_20 -> ~c
     .hi_1+ nu_eb   e-)'
      write(nout,104) brchelne(2,1),3,-ic1,ine,-ie,    'BR(~chi_20 -> ~c
     .hi_1- nu_e    e+)'
      endif
      if(brchelne(2,2).ne.0.D0) then
      write(nout,104) brchelne(2,2),3,ic2,-ine,ie,     'BR(~chi_20 -> ~c
     .hi_2+ nu_eb   e-)'
      write(nout,104) brchelne(2,2),3,-ic2,ine,-ie,    'BR(~chi_20 -> ~c
     .hi_2- nu_e    e+)'
      endif
      if(brchmunmu(2,1).ne.0.D0) then
      write(nout,104) brchmunmu(2,1),3,ic1,-inmu,imu,  'BR(~chi_20 -> ~c
     .hi_1+ nu_mub  mu-)'
      write(nout,104) brchmunmu(2,1),3,-ic1,inmu,-imu, 'BR(~chi_20 -> ~c
     .hi_1- nu_mu   mu+)'
      endif
      if(brchmunmu(2,2).ne.0.D0) then
      write(nout,104) brchmunmu(2,2),3,ic2,-inmu,imu,  'BR(~chi_20 -> ~c
     .hi_2+ nu_mub  mu-)'
      write(nout,104) brchmunmu(2,2),3,-ic2,inmu,-imu, 'BR(~chi_20 -> ~c
     .hi_2- nu_mu   mu+)'
      endif
      if(brchtauntau(2,1).ne.0.D0) then
      write(nout,104) brchtauntau(2,1),3,ic1,-intau,itau, 'BR(~chi_20 ->
     . ~chi_1+ nu_taub tau-)'
      write(nout,104) brchtauntau(2,1),3,-ic1,intau,-itau,'BR(~chi_20 ->
     . ~chi_1- nu_tau  tau+)'
      endif
      if(brchtauntau(2,2).ne.0.D0) then
      write(nout,104) brchtauntau(2,2),3,ic2,-intau,itau, 'BR(~chi_20 ->
     . ~chi_2+ nu_taub tau-)'
      write(nout,104) brchtauntau(2,2),3,-ic2,intau,-itau,'BR(~chi_20 ->
     . ~chi_2- nu_tau  tau+)'
      endif
      if(brglup(2).ne.0.D0) then
      write(nout,104) brglup(2),3,iglo,iub,iu,            'BR(~chi_20 ->
     . ~g      ub      u)'
      endif
      if(brgldo(2).ne.0.D0) then
      write(nout,104) brgldo(2),3,iglo,idb,id,            'BR(~chi_20 ->
     . ~g      db      d)'
      endif
      if(brglch(2).ne.0.D0) then
      write(nout,104) brglch(2),3,iglo,icb,ic,            'BR(~chi_20 ->
     . ~g      cb      c)'
      endif
      if(brglst(2).ne.0.D0) then
      write(nout,104) brglst(2),3,iglo,isb,is,            'BR(~chi_20 ->
     . ~g      sb      s)'
      endif
      if(brgltop(2).ne.0.D0) then
      write(nout,104) brgltop(2),3,iglo,itb,it,           'BR(~chi_20 ->
     . ~g      tb      t)'
      endif
      if(brglbot(2).ne.0.D0) then
      write(nout,104) brglbot(2),3,iglo,ibb,ib,           'BR(~chi_20 ->
     . ~g      bb      b)'
      endif
      endif
c ------------------ c
c neutralino3 decays c
c ------------------ c
      write(nout,99)
      write(nout,100) 1000025,neuttot(3),'neutralino3'
      if(neuttot2(3).ne.0.D0) then
      write(nout,49) 'neutralino3 2-body decays'
      write(nout,101)
      if(brneutzneut(3,1).ne.0.D0) then
      write(nout,102) brneutzneut(3,1),2,in1,iz,       'BR(~chi_30 -> ~c
     .hi_10   Z )'
      endif
      if(brneutzneut(3,2).ne.0.D0) then
      write(nout,102) brneutzneut(3,2),2,in2,iz,       'BR(~chi_30 -> ~c
     .hi_20   Z )'
      endif
      if(brneutwchar(3,1).ne.0.D0) then
      write(nout,102) brneutwchar(3,1),2,ic1,-iwc,     'BR(~chi_30 -> ~c
     .hi_1+   W-)'
      write(nout,102) brneutwchar(3,1),2,-ic1,iwc,     'BR(~chi_30 -> ~c
     .hi_1-   W+)'
      endif
      if(brneutwchar(3,2).ne.0.D0) then
      write(nout,102) brneutwchar(3,2),2,ic2,-iwc,     'BR(~chi_30 -> ~c
     .hi_2+   W-)'
      write(nout,102) brneutwchar(3,2),2,-ic2,iwc,     'BR(~chi_30 -> ~c
     .hi_2-   W+)'
      endif
c
      if(brneutHneut(3,1,1).ne.0.D0) then
      write(nout,102) brneutHneut(3,1,1),2,in1,ihH1,   'BR(~chi_30 -> ~c
     .hi_10   H_1 )'
      endif
      if(brneutHneut(3,1,2).ne.0.D0) then
      write(nout,102) brneutHneut(3,1,2),2,in1,ihH2,   'BR(~chi_30 -> ~c
     .hi_10   H_2 )'
      endif
      if(brneutHneut(3,1,3).ne.0.D0) then
      write(nout,102) brneutHneut(3,1,3),2,in1,ihH3,   'BR(~chi_30 -> ~c
     .hi_10   H_3 )'
      endif
      if(brneutAneut(3,1,1).ne.0.D0) then
      write(nout,102) brneutAneut(3,1,1),2,in1,ihA1,   'BR(~chi_30 -> ~c
     .hi_10   A_1 )'
      endif
      if(brneutAneut(3,1,2).ne.0.D0) then
      write(nout,102) brneutAneut(3,1,2),2,in1,ihA2,   'BR(~chi_30 -> ~c
     .hi_10   A_2 )'
      endif
c
      if(brneutHneut(3,2,1).ne.0.D0) then
      write(nout,102) brneutHneut(3,2,1),2,in2,ihH1,   'BR(~chi_30 -> ~c
     .hi_20   H_1 )'
      endif
      if(brneutHneut(3,2,2).ne.0.D0) then
      write(nout,102) brneutHneut(3,2,2),2,in2,ihH2,   'BR(~chi_30 -> ~c
     .hi_20   H_2 )'
      endif
      if(brneutHneut(3,2,3).ne.0.D0) then
      write(nout,102) brneutHneut(3,2,3),2,in2,ihH3,   'BR(~chi_30 -> ~c
     .hi_20   H_3 )'
      endif
      if(brneutAneut(3,2,1).ne.0.D0) then
      write(nout,102) brneutAneut(3,2,1),2,in2,ihA1,   'BR(~chi_30 -> ~c
     .hi_20   A_1 )'
      endif
      if(brneutAneut(3,2,2).ne.0.D0) then
      write(nout,102) brneutAneut(3,2,2),2,in2,ihA2,   'BR(~chi_30 -> ~c
     .hi_20   A_2 )'
      endif

      if(brneuthcchar(3,1).ne.0.D0) then
      write(nout,102) brneuthcchar(3,1),2,ic1,-ihc,    'BR(~chi_30 -> ~c
     .hi_1+   H-)'
      write(nout,102) brneuthcchar(3,1),2,-ic1,ihc,    'BR(~chi_30 -> ~c
     .hi_1-   H+)'
      endif
      if(brneuthcchar(3,2).ne.0.D0) then
      write(nout,102) brneuthcchar(3,2),2,ic2,-ihc,    'BR(~chi_30 -> ~c
     .hi_2+   H-)'
      write(nout,102) brneuthcchar(3,2),2,-ic2,ihc,    'BR(~chi_30 -> ~c
     .hi_2-   H+)'
      endif
      if(brneutsupl(3).ne.0.D0) then
      write(nout,102) brneutsupl(3),2,isul,iub,        'BR(~chi_30 -> ~u
     ._L      ub)'
      write(nout,102) brneutsupl(3),2,-isul,iu,        'BR(~chi_30 -> ~u
     ._L*     u )'
      endif
      if(brneutsupr(3).ne.0.D0) then
      write(nout,102) brneutsupr(3),2,isur,iub,        'BR(~chi_30 -> ~u
     ._R      ub)'
      write(nout,102) brneutsupr(3),2,-isur,iu,        'BR(~chi_30 -> ~u
     ._R*     u )'
      endif
      if(brneutsdownl(3).ne.0.D0) then
      write(nout,102) brneutsdownl(3),2,isdl,idb,      'BR(~chi_30 -> ~d
     ._L      db)'
      write(nout,102) brneutsdownl(3),2,-isdl,id,      'BR(~chi_30 -> ~d
     ._L*     d )'
      endif
      if(brneutsdownr(3).ne.0.D0) then
      write(nout,102) brneutsdownr(3),2,isdr,idb,      'BR(~chi_30 -> ~d
     ._R      db)'
      write(nout,102) brneutsdownr(3),2,-isdr,id,      'BR(~chi_30 -> ~d
     ._R*     d )'
      endif
      if(brneutsupl(3).ne.0.D0) then
      write(nout,102) brneutsupl(3),2,iscl,icb,        'BR(~chi_30 -> ~c
     ._L      cb)'
      write(nout,102) brneutsupl(3),2,-iscl,ic,        'BR(~chi_30 -> ~c
     ._L*     c )'
      endif
      if(brneutsupr(3).ne.0.D0) then
      write(nout,102) brneutsupr(3),2,iscr,icb,        'BR(~chi_30 -> ~c
     ._R      cb)'
      write(nout,102) brneutsupr(3),2,-iscr,ic,        'BR(~chi_30 -> ~c
     ._R*     c )'
      endif
      if(brneutsdownl(3).ne.0.D0) then
      write(nout,102) brneutsdownl(3),2,issl,isb,      'BR(~chi_30 -> ~s
     ._L      sb)'
      write(nout,102) brneutsdownl(3),2,-issl,is,      'BR(~chi_30 -> ~s
     ._L*     s )'
      endif
      if(brneutsdownr(3).ne.0.D0) then
      write(nout,102) brneutsdownr(3),2,issr,isb,      'BR(~chi_30 -> ~s
     ._R      sb)'
      write(nout,102) brneutsdownr(3),2,-issr,is,      'BR(~chi_30 -> ~s
     ._R*     s )'
      endif
      if(brneutst1(3).ne.0.D0) then
      write(nout,102) brneutst1(3),2,ist1,itb,         'BR(~chi_30 -> ~t
     ._1      tb)'
      write(nout,102) brneutst1(3),2,-ist1,it,         'BR(~chi_30 -> ~t
     ._1*     t )'
      endif
      if(brneutst2(3).ne.0.D0) then
      write(nout,102) brneutst2(3),2,ist2,itb,         'BR(~chi_30 -> ~t
     ._2      tb)'
      write(nout,102) brneutst2(3),2,-ist2,it,         'BR(~chi_30 -> ~t
     ._2*     t )'
      endif
      if(brneutsb1(3).ne.0.D0) then
      write(nout,102) brneutsb1(3),2,isb1,ibb,         'BR(~chi_30 -> ~b
     ._1      bb)'
      write(nout,102) brneutsb1(3),2,-isb1,ib,         'BR(~chi_30 -> ~b
     ._1*     b )'
      endif
      if(brneutsb2(3).ne.0.D0) then
      write(nout,102) brneutsb2(3),2,isb2,ibb,         'BR(~chi_30 -> ~b
     ._2      bb)'
      write(nout,102) brneutsb2(3),2,-isb2,ib,         'BR(~chi_30 -> ~b
     ._2*     b )'
      endif
      if(brneutsell(3).ne.0.D0) then
      write(nout,102) brneutsell(3),2,isell,-ie,       'BR(~chi_30 -> ~e
     ._L-     e+)'
      write(nout,102) brneutsell(3),2,-isell,ie,       'BR(~chi_30 -> ~e
     ._L+     e-)'
      endif
      if(brneutselr(3).ne.0.D0) then
      write(nout,102) brneutselr(3),2,iselr,-ie,       'BR(~chi_30 -> ~e
     ._R-     e+)'
      write(nout,102) brneutselr(3),2,-iselr,ie,       'BR(~chi_30 -> ~e
     ._R+     e-)'
      endif
      if(brneutsell(3).ne.0.D0) then
      write(nout,102) brneutsell(3),2,ismul,-imu,      'BR(~chi_30 -> ~m
     .u_L-    mu+)' 
      write(nout,102) brneutsell(3),2,-ismul,imu,      'BR(~chi_30 -> ~m
     .u_L+    mu-)' 
      endif
      if(brneutselr(3).ne.0.D0) then
      write(nout,102) brneutselr(3),2,ismur,-imu,      'BR(~chi_30 -> ~m
     .u_R-    mu+)' 
      write(nout,102) brneutselr(3),2,-ismur,imu,      'BR(~chi_30 -> ~m
     .u_R+    mu-)'
      endif
      if(brneutstau1(3).ne.0.D0) then
      write(nout,102) brneutstau1(3),2,istau1,-itau,   'BR(~chi_30 -> ~t
     .au_1-   tau+)'
      write(nout,102) brneutstau1(3),2,-istau1,itau,   'BR(~chi_30 -> ~t
     .au_1+   tau-)'
      endif
      if(brneutstau2(3).ne.0.D0) then
      write(nout,102) brneutstau2(3),2,istau2,-itau,   'BR(~chi_30 -> ~t
     .au_2-   tau+)'
      write(nout,102) brneutstau2(3),2,-istau2,itau,   'BR(~chi_30 -> ~t
     .au_2+   tau-)'
      endif
      if(brneutsnel(3).ne.0.D0) then
      write(nout,102) brneutsnel(3),2,inel,-ine,       'BR(~chi_30 -> ~n
     .u_eL    nu_eb)'
      write(nout,102) brneutsnel(3),2,-inel,ine,       'BR(~chi_30 -> ~n
     .u_eL*   nu_e )'
      write(nout,102) brneutsnel(3),2,inmul,-inmu,     'BR(~chi_30 -> ~n
     .u_muL   nu_mub)'
      write(nout,102) brneutsnel(3),2,-inmul,inmu,     'BR(~chi_30 -> ~n
     .u_muL*  nu_mu )'
      endif
      if(brneutsn1(3).ne.0.D0) then
      write(nout,102) brneutsn1(3),2,intau1,-intau,    'BR(~chi_30 -> ~n
     .u_tau1  nu_taub)'
      write(nout,102) brneutsn1(3),2,-intau1,intau,    'BR(~chi_30 -> ~n
     .u_tau1* nu_tau )'
      endif
      endif
C     =============================
*     NEUTRALINO3 : Radiative decay
C     =============================
            if(flagloop.eq.1.D0) then
      if(brnraddec(3,1).ne.0.D0) then
      write(nout,102) brnraddec(3,1),2,in1,iga,        'BR(~chi_30 -> ~c
     .hi_10 gam)'
      endif
      if(brnraddec(3,2).ne.0.D0) then
      write(nout,102) brnraddec(3,2),2,in2,iga,        'BR(~chi_30 -> ~c
     .hi_20 gam)'
      endif
            endif
C     =============================
*     NEUTRALINO3 : Three body
C     =============================

      if(neuttot3(3).ne.0.D0) then 
      write(nout,49) 'neutralino3 3-body decays'
      write(nout,103)
      if(brneutup(3,1).ne.0.D0) then
      write(nout,104) brneutup(3,1),3,in1,iub,iu,      'BR(~chi_30 -> ~c
     .hi_10 ub      u)'
      endif
      if(brneutdow(3,1).ne.0.D0) then
      write(nout,104) brneutdow(3,1),3,in1,idb,id,     'BR(~chi_30 -> ~c
     .hi_10 db      d)'
      endif
      if(brneutch(3,1).ne.0.D0) then
      write(nout,104) brneutch(3,1),3,in1,icb,ic,      'BR(~chi_30 -> ~c
     .hi_10 cb      c)'
      endif
      if(brneutst(3,1).ne.0.D0) then
      write(nout,104) brneutst(3,1),3,in1,isb,is,      'BR(~chi_30 -> ~c
     .hi_10 sb      s)'
      endif
      if(brneuttop(3,1).ne.0.D0) then
      write(nout,104) brneuttop(3,1),3,in1,itb,it,     'BR(~chi_30 -> ~c
     .hi_10 tb      t)'
      endif
      if(brneutbot(3,1).ne.0.D0) then
      write(nout,104) brneutbot(3,1),3,in1,ibb,ib,     'BR(~chi_30 -> ~c
     .hi_10 bb      b)'
      endif
      if(brneutel(3,1).ne.0.D0) then
      write(nout,104) brneutel(3,1),3,in1,-ie,ie,      'BR(~chi_30 -> ~c
     .hi_10 e+      e-)'
      endif
      if(brneutmu(3,1).ne.0.D0) then
      write(nout,104) brneutmu(3,1),3,in1,-imu,imu,    'BR(~chi_30 -> ~c
     .hi_10 mu+     mu-)'
      endif
      if(brneuttau(3,1).ne.0.D0) then
      write(nout,104) brneuttau(3,1),3,in1,-itau,itau, 'BR(~chi_30 -> ~c
     .hi_10 tau+    tau-)'
      endif
      if(brneutnue(3,1).ne.0.D0) then
      write(nout,104) brneutnue(3,1),3,in1,-ine,ine,   'BR(~chi_30 -> ~c
     .hi_10 nu_eb   nu_e)'
      endif
      if(brneutnumu(3,1).ne.0.D0) then
      write(nout,104) brneutnumu(3,1),3,in1,-inmu,inmu,'BR(~chi_30 -> ~c
     .hi_10 nu_mub  nu_mu)'
      endif
      if(brneutnutau(3,1).ne.0.D0) then
      write(nout,104) brneutnutau(3,1),3,in1,-intau,intau,'BR(~chi_30 ->
     . ~chi_10 nu_taub nu_tau)'
      endif
      if(brneutup(3,2).ne.0.D0) then
      write(nout,104) brneutup(3,2),3,in2,iub,iu,      'BR(~chi_30 -> ~c
     .hi_20 ub      u)'
      endif
      if(brneutdow(3,2).ne.0.D0) then
      write(nout,104) brneutdow(3,2),3,in2,idb,id,     'BR(~chi_30 -> ~c
     .hi_20 db      d)'
      endif
      if(brneutch(3,2).ne.0.D0) then
      write(nout,104) brneutch(3,2),3,in2,icb,ic,      'BR(~chi_30 -> ~c
     .hi_20 cb      c)'
      endif
      if(brneutst(3,2).ne.0.D0) then
      write(nout,104) brneutst(3,2),3,in2,isb,is,      'BR(~chi_30 -> ~c
     .hi_20 sb      s)'
      endif
      if(brneuttop(3,2).ne.0.D0) then
      write(nout,104) brneuttop(3,2),3,in2,itb,it,     'BR(~chi_30 -> ~c
     .hi_20 tb      t)'
      endif
      if(brneutbot(3,2).ne.0.D0) then
      write(nout,104) brneutbot(3,2),3,in2,ibb,ib,     'BR(~chi_30 -> ~c
     .hi_20 bb      b)'
      endif
      if(brneutel(3,2).ne.0.D0) then
      write(nout,104) brneutel(3,2),3,in2,-ie,ie,      'BR(~chi_30 -> ~c
     .hi_20 e+      e-)'
      endif
      if(brneutmu(3,2).ne.0.D0) then
      write(nout,104) brneutmu(3,2),3,in2,-imu,imu,    'BR(~chi_30 -> ~c
     .hi_20 mu+     mu-)'
      endif
      if(brneuttau(3,2).ne.0.D0) then
      write(nout,104) brneuttau(3,2),3,in2,-itau,itau, 'BR(~chi_30 -> ~c
     .hi_20 tau+    tau-)'
      endif
      if(brneutnue(3,2).ne.0.D0) then
      write(nout,104) brneutnue(3,2),3,in2,-ine,ine,   'BR(~chi_30 -> ~c
     .hi_20 nu_eb   nu_e)'
      endif
      if(brneutnumu(3,2).ne.0.D0) then
      write(nout,104) brneutnumu(3,2),3,in2,-inmu,inmu,'BR(~chi_30 -> ~c
     .hi_20 nu_mub  nu_mu)'
      endif
      if(brneutnutau(3,2).ne.0.D0) then
      write(nout,104) brneutnutau(3,2),3,in2,-intau,intau,'BR(~chi_30 ->
     . ~chi_20 nu_taub nu_tau)'
      endif
      if(brchubd(3,1).ne.0.D0) then
      write(nout,104) brchubd(3,1),3,ic1,iub,id,       'BR(~chi_30 -> ~c
     .hi_1+ ub      d)'
      write(nout,104) brchubd(3,1),3,-ic1,idb,iu,      'BR(~chi_30 -> ~c
     .hi_1- db      u)'
      endif
      if(brchubd(3,2).ne.0.D0) then
      write(nout,104) brchubd(3,2),3,ic2,iub,id,       'BR(~chi_30 -> ~c
     .hi_2+ ub      d)'
      write(nout,104) brchubd(3,2),3,-ic2,idb,iu,      'BR(~chi_30 -> ~c
     .hi_2- db      u)'
      endif
      if(brchcbs(3,1).ne.0.D0) then
      write(nout,104) brchcbs(3,1),3,ic1,icb,is,       'BR(~chi_30 -> ~c
     .hi_1+ cb      s)'
      write(nout,104) brchcbs(3,1),3,-ic1,isb,ic,      'BR(~chi_30 -> ~c
     .hi_1- sb      c)'
      endif
      if(brchcbs(3,2).ne.0.D0) then
      write(nout,104) brchcbs(3,2),3,ic2,icb,is,       'BR(~chi_30 -> ~c
     .hi_2+ cb      s)'
      write(nout,104) brchcbs(3,2),3,-ic2,isb,ic,      'BR(~chi_30 -> ~c
     .hi_2- sb      c)'
      endif
      if(brchtbb(3,1).ne.0.D0) then
      write(nout,104) brchtbb(3,1),3,ic1,itb,ib,       'BR(~chi_30 -> ~c
     .hi_1+ tb      b)'
      write(nout,104) brchtbb(3,1),3,-ic1,ibb,it,      'BR(~chi_30 -> ~c
     .hi_1- bb      t)'
      endif
      if(brchtbb(3,2).ne.0.D0) then
      write(nout,104) brchtbb(3,2),3,ic2,itb,ib,       'BR(~chi_30 -> ~c
     .hi_2+ tb      b)'
      write(nout,104) brchtbb(3,2),3,-ic2,ibb,it,      'BR(~chi_30 -> ~c
     .hi_2- bb      t)'
      endif
      if(brchelne(3,1).ne.0.D0) then
      write(nout,104) brchelne(3,1),3,ic1,-ine,ie,     'BR(~chi_30 -> ~c
     .hi_1+ nu_eb   e-)'
      write(nout,104) brchelne(3,1),3,-ic1,ine,-ie,    'BR(~chi_30 -> ~c
     .hi_1- nu_e    e+)'
      endif
      if(brchelne(3,2).ne.0.D0) then
      write(nout,104) brchelne(3,2),3,ic2,-ine,ie,     'BR(~chi_30 -> ~c
     .hi_2+ nu_eb   e-)'
      write(nout,104) brchelne(3,2),3,-ic2,ine,-ie,    'BR(~chi_30 -> ~c
     .hi_2- nu_e    e+)'
      endif
      if(brchmunmu(3,1).ne.0.D0) then
      write(nout,104) brchmunmu(3,1),3,ic1,-inmu,imu,  'BR(~chi_30 -> ~c
     .hi_1+ nu_mub  mu-)'
      write(nout,104) brchmunmu(3,1),3,-ic1,inmu,-imu, 'BR(~chi_30 -> ~c
     .hi_1- nu_mu   mu+)'
      endif
      if(brchmunmu(3,2).ne.0.D0) then
      write(nout,104) brchmunmu(3,2),3,ic2,-inmu,imu,  'BR(~chi_30 -> ~c
     .hi_2+ nu_mub  mu-)'
      write(nout,104) brchmunmu(3,2),3,-ic2,inmu,-imu, 'BR(~chi_30 -> ~c
     .hi_2- nu_mu   mu+)'
      endif
      if(brchtauntau(3,1).ne.0.D0) then
      write(nout,104) brchtauntau(3,1),3,ic1,-intau,itau, 'BR(~chi_30 ->
     . ~chi_1+ nu_taub tau-)'
      write(nout,104) brchtauntau(3,1),3,-ic1,intau,-itau,'BR(~chi_30 ->
     . ~chi_1- nu_tau  tau+)'
      endif
      if(brchtauntau(3,2).ne.0.D0) then
      write(nout,104) brchtauntau(3,2),3,ic2,-intau,itau, 'BR(~chi_30 ->
     . ~chi_2+ nu_taub tau-)'
      write(nout,104) brchtauntau(3,2),3,-ic2,intau,-itau,'BR(~chi_30 ->
     . ~chi_2- nu_tau  tau+)'
      endif
      if(brglup(3).ne.0.D0) then
      write(nout,104) brglup(3),3,iglo,iub,iu,            'BR(~chi_30 ->
     . ~g      ub      u)'
      endif
      if(brgldo(3).ne.0.D0) then
      write(nout,104) brgldo(3),3,iglo,idb,id,            'BR(~chi_30 ->
     . ~g      db      d)'
      endif
      if(brglch(3).ne.0.D0) then
      write(nout,104) brglch(3),3,iglo,icb,ic,            'BR(~chi_30 ->
     . ~g      cb      c)'
      endif
      if(brglst(3).ne.0.D0) then
      write(nout,104) brglst(3),3,iglo,isb,is,            'BR(~chi_30 ->
     . ~g      sb      s)'
      endif
      if(brgltop(3).ne.0.D0) then
      write(nout,104) brgltop(3),3,iglo,itb,it,           'BR(~chi_30 ->
     . ~g      tb      t)'
      endif
      if(brglbot(3).ne.0.D0) then
      write(nout,104) brglbot(3),3,iglo,ibb,ib,           'BR(~chi_30 ->
     . ~g      bb      b)'
      endif
      endif

c ------------------ c
c neutralino4 decays c
c ------------------ c
      write(nout,99)
      write(nout,100) 1000035,neuttot(4),'neutralino4'
      if(neuttot2(4).ne.0.D0) then
      write(nout,49) 'neutralino4 2-body decays'
      write(nout,101)
      if(brneutzneut(4,1).ne.0.D0) then
      write(nout,102) brneutzneut(4,1),2,in1,iz,       'BR(~chi_40 -> ~c
     .hi_10   Z )'
      endif
      if(brneutzneut(4,2).ne.0.D0) then
      write(nout,102) brneutzneut(4,2),2,in2,iz,       'BR(~chi_40 -> ~c
     .hi_20   Z )'
      endif
      if(brneutzneut(4,3).ne.0.D0) then
      write(nout,102) brneutzneut(4,3),2,in3,iz,       'BR(~chi_40 -> ~c
     .hi_30   Z )'
      endif
      if(brneutwchar(4,1).ne.0.D0) then
      write(nout,102) brneutwchar(4,1),2,ic1,-iwc,     'BR(~chi_40 -> ~c
     .hi_1+   W-)'
      write(nout,102) brneutwchar(4,1),2,-ic1,iwc,     'BR(~chi_40 -> ~c
     .hi_1-   W+)'
      endif
      if(brneutwchar(4,2).ne.0.D0) then
      write(nout,102) brneutwchar(4,2),2,ic2,-iwc,     'BR(~chi_40 -> ~c
     .hi_2+   W-)'
      write(nout,102) brneutwchar(4,2),2,-ic2,iwc,     'BR(~chi_40 -> ~c
     .hi_2-   W+)'
      endif
c
      if(brneutHneut(4,1,1).ne.0.D0) then
      write(nout,102) brneutHneut(4,1,1),2,in1,ihH1,   'BR(~chi_40 -> ~c
     .hi_10   H_1 )'
      endif
      if(brneutHneut(4,1,2).ne.0.D0) then
      write(nout,102) brneutHneut(4,1,2),2,in1,ihH2,   'BR(~chi_40 -> ~c
     .hi_10   H_2 )'
      endif
      if(brneutHneut(4,1,3).ne.0.D0) then
      write(nout,102) brneutHneut(4,1,3),2,in1,ihH3,   'BR(~chi_40 -> ~c
     .hi_10   H_3 )'
      endif
      if(brneutAneut(4,1,1).ne.0.D0) then
      write(nout,102) brneutAneut(4,1,1),2,in1,ihA1,   'BR(~chi_40 -> ~c
     .hi_10   A_1 )'
      endif
      if(brneutAneut(4,1,2).ne.0.D0) then
      write(nout,102) brneutAneut(4,1,2),2,in1,ihA2,   'BR(~chi_40 -> ~c
     .hi_10   A_2 )'
      endif
c
      if(brneutHneut(4,2,1).ne.0.D0) then
      write(nout,102) brneutHneut(4,2,1),2,in2,ihH1,   'BR(~chi_40 -> ~c
     .hi_20   H_1 )'
      endif
      if(brneutHneut(4,2,2).ne.0.D0) then
      write(nout,102) brneutHneut(4,2,2),2,in2,ihH2,   'BR(~chi_40 -> ~c
     .hi_20   H_2 )'
      endif
      if(brneutHneut(4,2,3).ne.0.D0) then
      write(nout,102) brneutHneut(4,2,3),2,in2,ihH3,   'BR(~chi_40 -> ~c
     .hi_20   H_3 )'
      endif
      if(brneutAneut(4,2,1).ne.0.D0) then
      write(nout,102) brneutAneut(4,2,1),2,in2,ihA1,   'BR(~chi_40 -> ~c
     .hi_20   A_1 )'
      endif
      if(brneutAneut(4,2,2).ne.0.D0) then
      write(nout,102) brneutAneut(4,2,2),2,in2,ihA2,   'BR(~chi_40 -> ~c
     .hi_20   A_2 )'
      endif
c
      if(brneutHneut(4,3,1).ne.0.D0) then
      write(nout,102) brneutHneut(4,3,1),2,in3,ihH1,   'BR(~chi_40 -> ~c
     .hi_30   H_1 )'
      endif
      if(brneutHneut(4,3,2).ne.0.D0) then
      write(nout,102) brneutHneut(4,3,2),2,in3,ihH2,   'BR(~chi_40 -> ~c
     .hi_30   H_2 )'
      endif
      if(brneutHneut(4,3,3).ne.0.D0) then
      write(nout,102) brneutHneut(4,3,3),2,in3,ihH3,   'BR(~chi_40 -> ~c
     .hi_30   H_3 )'
      endif
      if(brneutAneut(4,3,1).ne.0.D0) then
      write(nout,102) brneutAneut(4,3,1),2,in3,ihA1,   'BR(~chi_40 -> ~c
     .hi_30   A_1 )'
      endif
      if(brneutAneut(4,3,2).ne.0.D0) then
      write(nout,102) brneutAneut(4,3,2),2,in3,ihA2,   'BR(~chi_40 -> ~c
     .hi_30   A_2 )'
      endif
      if(brneuthcchar(4,1).ne.0.D0) then
      write(nout,102) brneuthcchar(4,1),2,ic1,-ihc,    'BR(~chi_40 -> ~c
     .hi_1+   H-)'
      write(nout,102) brneuthcchar(4,1),2,-ic1,ihc,    'BR(~chi_40 -> ~c
     .hi_1-   H+)'
      endif
      if(brneuthcchar(4,2).ne.0.D0) then
      write(nout,102) brneuthcchar(4,2),2,ic2,-ihc,    'BR(~chi_40 -> ~c
     .hi_2+   H-)'
      write(nout,102) brneuthcchar(4,2),2,-ic2,ihc,    'BR(~chi_40 -> ~c
     .hi_2-   H+)'
      endif
      if(brneutsupl(4).ne.0.D0) then
      write(nout,102) brneutsupl(4),2,isul,iub,        'BR(~chi_40 -> ~u
     ._L      ub)'
      write(nout,102) brneutsupl(4),2,-isul,iu,        'BR(~chi_40 -> ~u
     ._L*     u )'
      endif
      if(brneutsupr(4).ne.0.D0) then
      write(nout,102) brneutsupr(4),2,isur,iub,        'BR(~chi_40 -> ~u
     ._R      ub)'
      write(nout,102) brneutsupr(4),2,-isur,iu,        'BR(~chi_40 -> ~u
     ._R*     u )'
      endif
      if(brneutsdownl(4).ne.0.D0) then
      write(nout,102) brneutsdownl(4),2,isdl,idb,      'BR(~chi_40 -> ~d
     ._L      db)'
      write(nout,102) brneutsdownl(4),2,-isdl,id,      'BR(~chi_40 -> ~d
     ._L*     d )'
      endif
      if(brneutsdownr(4).ne.0.D0) then
      write(nout,102) brneutsdownr(4),2,isdr,idb,      'BR(~chi_40 -> ~d
     ._R      db)'
      write(nout,102) brneutsdownr(4),2,-isdr,id,      'BR(~chi_40 -> ~d
     ._R*     d )'
      endif
      if(brneutsupl(4).ne.0.D0) then
      write(nout,102) brneutsupl(4),2,iscl,icb,        'BR(~chi_40 -> ~c
     ._L      cb)'
      write(nout,102) brneutsupl(4),2,-iscl,ic,        'BR(~chi_40 -> ~c
     ._L*     c )'
      endif
      if(brneutsupr(4).ne.0.D0) then
      write(nout,102) brneutsupr(4),2,iscr,icb,        'BR(~chi_40 -> ~c
     ._R      cb)'
      write(nout,102) brneutsupr(4),2,-iscr,ic,        'BR(~chi_40 -> ~c
     ._R*     c )'
      endif
      if(brneutsdownl(4).ne.0.D0) then
      write(nout,102) brneutsdownl(4),2,issl,isb,      'BR(~chi_40 -> ~s
     ._L      sb)'
      write(nout,102) brneutsdownl(4),2,-issl,is,      'BR(~chi_40 -> ~s
     ._L*     s )'
      endif
      if(brneutsdownr(4).ne.0.D0) then
      write(nout,102) brneutsdownr(4),2,issr,isb,      'BR(~chi_40 -> ~s
     ._R      sb)'
      write(nout,102) brneutsdownr(4),2,-issr,is,      'BR(~chi_40 -> ~s
     ._R*     s )'
      endif
      if(brneutst1(4).ne.0.D0) then
      write(nout,102) brneutst1(4),2,ist1,itb,         'BR(~chi_40 -> ~t
     ._1      tb)'
      write(nout,102) brneutst1(4),2,-ist1,it,         'BR(~chi_40 -> ~t
     ._1*     t )'
      endif
      if(brneutst2(4).ne.0.D0) then
      write(nout,102) brneutst2(4),2,ist2,itb,         'BR(~chi_40 -> ~t
     ._2      tb)'
      write(nout,102) brneutst2(4),2,-ist2,it,         'BR(~chi_40 -> ~t
     ._2*     t )'
      endif
      if(brneutsb1(4).ne.0.D0) then
      write(nout,102) brneutsb1(4),2,isb1,ibb,         'BR(~chi_40 -> ~b
     ._1      bb)'
      write(nout,102) brneutsb1(4),2,-isb1,ib,         'BR(~chi_40 -> ~b
     ._1*     b )'
      endif
      if(brneutsb2(4).ne.0.D0) then
      write(nout,102) brneutsb2(4),2,isb2,ibb,         'BR(~chi_40 -> ~b
     ._2      bb)'
      write(nout,102) brneutsb2(4),2,-isb2,ib,         'BR(~chi_40 -> ~b
     ._2*     b )'
      endif
      if(brneutsell(4).ne.0.D0) then
      write(nout,102) brneutsell(4),2,isell,-ie,       'BR(~chi_40 -> ~e
     ._L-     e+)'
      write(nout,102) brneutsell(4),2,-isell,ie,       'BR(~chi_40 -> ~e
     ._L+     e-)'
      endif
      if(brneutselr(4).ne.0.D0) then
      write(nout,102) brneutselr(4),2,iselr,-ie,       'BR(~chi_40 -> ~e
     ._R-     e+)'
      write(nout,102) brneutselr(4),2,-iselr,ie,       'BR(~chi_40 -> ~e
     ._R+     e-)'
      endif
      if(brneutsell(4).ne.0.D0) then
      write(nout,102) brneutsell(4),2,ismul,-imu,      'BR(~chi_40 -> ~m
     .u_L-    mu+)' 
      write(nout,102) brneutsell(4),2,-ismul,imu,      'BR(~chi_40 -> ~m
     .u_L+    mu-)' 
      endif
      if(brneutselr(4).ne.0.D0) then      
      write(nout,102) brneutselr(4),2,ismur,-imu,      'BR(~chi_40 -> ~m
     .u_R-    mu+)' 
      write(nout,102) brneutselr(4),2,-ismur,imu,      'BR(~chi_40 -> ~m
     .u_R+    mu-)'
      endif
      if(brneutstau1(4).ne.0.D0) then
      write(nout,102) brneutstau1(4),2,istau1,-itau,   'BR(~chi_40 -> ~t
     .au_1-   tau+)'
      write(nout,102) brneutstau1(4),2,-istau1,itau,   'BR(~chi_40 -> ~t
     .au_1+   tau-)'
      endif
      if(brneutstau2(4).ne.0.D0) then
      write(nout,102) brneutstau2(4),2,istau2,-itau,   'BR(~chi_40 -> ~t
     .au_2-   tau+)'
      write(nout,102) brneutstau2(4),2,-istau2,itau,   'BR(~chi_40 -> ~t
     .au_2+   tau-)'
      endif
      if(brneutsnel(4).ne.0.D0) then
      write(nout,102) brneutsnel(4),2,inel,-ine,       'BR(~chi_40 -> ~n
     .u_eL    nu_eb)'
      write(nout,102) brneutsnel(4),2,-inel,ine,       'BR(~chi_40 -> ~n
     .u_eL*   nu_e )'
      write(nout,102) brneutsnel(4),2,inmul,-inmu,     'BR(~chi_40 -> ~n
     .u_muL   nu_mub)'
      write(nout,102) brneutsnel(4),2,-inmul,inmu,     'BR(~chi_40 -> ~n
     .u_muL*  nu_mu )'
      endif
      if(brneutsn1(4).ne.0.D0) then
      write(nout,102) brneutsn1(4),2,intau1,-intau,    'BR(~chi_40 -> ~n
     .u_tau1  nu_taub)'
      write(nout,102) brneutsn1(4),2,-intau1,intau,    'BR(~chi_40 -> ~n
     .u_tau1* nu_tau )'
      endif
      endif
C     =============================
*     NEUTRALINO4 : Radiative decay
C     =============================

            if(flagloop.eq.1.D0) then
      if(brnraddec(4,1).ne.0.D0) then
      write(nout,102) brnraddec(4,1),2,in1,iga,        'BR(~chi_40 -> ~c
     .hi_10 gam)'
      endif
      if(brnraddec(4,2).ne.0.D0) then
      write(nout,102) brnraddec(4,2),2,in2,iga,        'BR(~chi_40 -> ~c
     .hi_20 gam)'
      endif
      if(brnraddec(4,3).ne.0.D0) then
      write(nout,102) brnraddec(4,3),2,in3,iga,        'BR(~chi_40 -> ~c
     .hi_30 gam)'
      endif
            endif
C     =============================
*     NEUTRALINO4 : Three body
C     =============================

      if(neuttot3(4).ne.0.D0) then 
      write(nout,49) 'neutralino4 3-body decays'
      write(nout,103)
      if(brneutup(4,1).ne.0.D0) then
      write(nout,104) brneutup(4,1),3,in1,iub,iu,      'BR(~chi_40 -> ~c
     .hi_10 ub      u)'
      endif
      if(brneutdow(4,1).ne.0.D0) then
      write(nout,104) brneutdow(4,1),3,in1,idb,id,     'BR(~chi_40 -> ~c
     .hi_10 db      d)'
      endif
      if(brneutch(4,1).ne.0.D0) then
      write(nout,104) brneutch(4,1),3,in1,icb,ic,      'BR(~chi_40 -> ~c
     .hi_10 cb      c)'
      endif
      if(brneutst(4,1).ne.0.D0) then
      write(nout,104) brneutst(4,1),3,in1,isb,is,      'BR(~chi_40 -> ~c
     .hi_10 sb      s)'
      endif
      if(brneuttop(4,1).ne.0.D0) then
      write(nout,104) brneuttop(4,1),3,in1,itb,it,     'BR(~chi_40 -> ~c
     .hi_10 tb      t)'
      endif
      if(brneutbot(4,1).ne.0.D0) then
      write(nout,104) brneutbot(4,1),3,in1,ibb,ib,     'BR(~chi_40 -> ~c
     .hi_10 bb      b)'
      endif
      if(brneutel(4,1).ne.0.D0) then
      write(nout,104) brneutel(4,1),3,in1,-ie,ie,      'BR(~chi_40 -> ~c
     .hi_10 e+      e-)'
      endif
      if(brneutmu(4,1).ne.0.D0) then
      write(nout,104) brneutmu(4,1),3,in1,-imu,imu,    'BR(~chi_40 -> ~c
     .hi_10 mu+     mu-)'
      endif
      if(brneuttau(4,1).ne.0.D0) then
      write(nout,104) brneuttau(4,1),3,in1,-itau,itau, 'BR(~chi_40 -> ~c
     .hi_10 tau+    tau-)'
      endif
      if(brneutnue(4,1).ne.0.D0) then
      write(nout,104) brneutnue(4,1),3,in1,-ine,ine,   'BR(~chi_40 -> ~c
     .hi_10 nu_eb   nu_e)'
      endif
      if(brneutnumu(4,1).ne.0.D0) then
      write(nout,104) brneutnumu(4,1),3,in1,-inmu,inmu,'BR(~chi_40 -> ~c
     .hi_10 nu_mub  nu_mu)'
      endif
      if(brneutnutau(4,1).ne.0.D0) then
      write(nout,104) brneutnutau(4,1),3,in1,-intau,intau,'BR(~chi_40 ->
     . ~chi_10 nu_taub nu_tau)'
      endif
      if(brneutup(4,2).ne.0.D0) then
      write(nout,104) brneutup(4,2),3,in2,iub,iu,      'BR(~chi_40 -> ~c
     .hi_20 ub      u)'
      endif
      if(brneutdow(4,2).ne.0.D0) then
      write(nout,104) brneutdow(4,2),3,in2,idb,id,     'BR(~chi_40 -> ~c
     .hi_20 db      d)'
      endif
      if(brneutch(4,2).ne.0.D0) then
      write(nout,104) brneutch(4,2),3,in2,icb,ic,      'BR(~chi_40 -> ~c
     .hi_20 cb      c)'
      endif
      if(brneutst(4,2).ne.0.D0) then
      write(nout,104) brneutst(4,2),3,in2,isb,is,      'BR(~chi_40 -> ~c
     .hi_20 sb      s)'
      endif
      if(brneuttop(4,2).ne.0.D0) then
      write(nout,104) brneuttop(4,2),3,in2,itb,it,     'BR(~chi_40 -> ~c
     .hi_20 tb      t)'
      endif
      if(brneutbot(4,2).ne.0.D0) then
      write(nout,104) brneutbot(4,2),3,in2,ibb,ib,     'BR(~chi_40 -> ~c
     .hi_20 bb      b)'
      endif
      if(brneutel(4,2).ne.0.D0) then
      write(nout,104) brneutel(4,2),3,in2,-ie,ie,      'BR(~chi_40 -> ~c
     .hi_20 e+      e-)'
      endif
      if(brneutmu(4,2).ne.0.D0) then
      write(nout,104) brneutmu(4,2),3,in2,-imu,imu,    'BR(~chi_40 -> ~c
     .hi_20 mu+     mu-)'
      endif
      if(brneuttau(4,2).ne.0.D0) then
      write(nout,104) brneuttau(4,2),3,in2,-itau,itau, 'BR(~chi_40 -> ~c
     .hi_20 tau+    tau-)'
      endif
      if(brneutnue(4,2).ne.0.D0) then
      write(nout,104) brneutnue(4,2),3,in2,-ine,ine,   'BR(~chi_40 -> ~c
     .hi_20 nu_eb   nu_e)'
      endif
      if(brneutnumu(4,2).ne.0.D0) then
      write(nout,104) brneutnumu(4,2),3,in2,-inmu,inmu,'BR(~chi_40 -> ~c
     .hi_20 nu_mub  nu_mu)'
      endif
      if(brneutnutau(4,2).ne.0.D0) then
      write(nout,104) brneutnutau(4,2),3,in2,-intau,intau,'BR(~chi_40 ->
     . ~chi_20 nu_taub nu_tau)'
      endif
      if(brneutup(4,3).ne.0.D0) then
      write(nout,104) brneutup(4,3),3,in3,iub,iu,      'BR(~chi_40 -> ~c
     .hi_30 ub      u)'
      endif
      if(brneutdow(4,3).ne.0.D0) then
      write(nout,104) brneutdow(4,3),3,in3,idb,id,     'BR(~chi_40 -> ~c
     .hi_30 db      d)'
      endif
      if(brneutch(4,3).ne.0.D0) then
      write(nout,104) brneutch(4,3),3,in3,icb,ic,      'BR(~chi_40 -> ~c
     .hi_30 cb      c)'
      endif
      if(brneutst(4,3).ne.0.D0) then
      write(nout,104) brneutst(4,3),3,in3,isb,is,      'BR(~chi_40 -> ~c
     .hi_30 sb      s)'
      endif
      if(brneuttop(4,3).ne.0.D0) then
      write(nout,104) brneuttop(4,3),3,in3,itb,it,     'BR(~chi_40 -> ~c
     .hi_30 tb      t)'
      endif
      if(brneutbot(4,3).ne.0.D0) then
      write(nout,104) brneutbot(4,3),3,in3,ibb,ib,     'BR(~chi_40 -> ~c
     .hi_30 bb      b)'
      endif
      if(brneutel(4,3).ne.0.D0) then
      write(nout,104) brneutel(4,3),3,in3,-ie,ie,      'BR(~chi_40 -> ~c
     .hi_30 e+      e-)'
      endif
      if(brneutmu(4,3).ne.0.D0) then
      write(nout,104) brneutmu(4,3),3,in3,-imu,imu,    'BR(~chi_40 -> ~c
     .hi_30 mu+     mu-)'
      endif
      if(brneuttau(4,3).ne.0.D0) then
      write(nout,104) brneuttau(4,3),3,in3,-itau,itau, 'BR(~chi_40 -> ~c
     .hi_30 tau+    tau-)'
      endif
      if(brneutnue(4,3).ne.0.D0) then
      write(nout,104) brneutnue(4,3),3,in3,-ine,ine,   'BR(~chi_40 -> ~c
     .hi_30 nu_eb   nu_e)'
      endif
      if(brneutnumu(4,3).ne.0.D0) then
      write(nout,104) brneutnumu(4,3),3,in3,-inmu,inmu,'BR(~chi_40 -> ~c
     .hi_30 nu_mub  nu_mu)'
      endif
      if(brneutnutau(4,3).ne.0.D0) then
      write(nout,104) brneutnutau(4,3),3,in3,-intau,intau,'BR(~chi_40 ->
     . ~chi_30 nu_taub nu_tau)'
      endif
      if(brchubd(4,1).ne.0.D0) then
      write(nout,104) brchubd(4,1),3,ic1,iub,id,       'BR(~chi_40 -> ~c
     .hi_1+ ub      d)'
      write(nout,104) brchubd(4,1),3,-ic1,idb,iu,      'BR(~chi_40 -> ~c
     .hi_1- db      u)'
      endif
      if(brchubd(4,2).ne.0.D0) then
      write(nout,104) brchubd(4,2),3,ic2,iub,id,       'BR(~chi_40 -> ~c
     .hi_2+ ub      d)'
      write(nout,104) brchubd(4,2),3,-ic2,idb,iu,      'BR(~chi_40 -> ~c
     .hi_2- db      u)'
      endif
      if(brchcbs(4,1).ne.0.D0) then
      write(nout,104) brchcbs(4,1),3,ic1,icb,is,       'BR(~chi_40 -> ~c
     .hi_1+ cb      s)'
      write(nout,104) brchcbs(4,1),3,-ic1,isb,ic,      'BR(~chi_40 -> ~c
     .hi_1- sb      c)'
      endif
      if(brchcbs(4,2).ne.0.D0) then
      write(nout,104) brchcbs(4,2),3,ic2,icb,is,       'BR(~chi_40 -> ~c
     .hi_2+ cb      s)'
      write(nout,104) brchcbs(4,2),3,-ic2,isb,ic,      'BR(~chi_40 -> ~c
     .hi_2- sb      c)'
      endif
      if(brchtbb(4,1).ne.0.D0) then
      write(nout,104) brchtbb(4,1),3,ic1,itb,ib,       'BR(~chi_40 -> ~c
     .hi_1+ tb      b)'
      write(nout,104) brchtbb(4,1),3,-ic1,ibb,it,      'BR(~chi_40 -> ~c
     .hi_1- bb      t)'
      endif
      if(brchtbb(4,2).ne.0.D0) then
      write(nout,104) brchtbb(4,2),3,ic2,itb,ib,       'BR(~chi_40 -> ~c
     .hi_2+ tb      b)'
      write(nout,104) brchtbb(4,2),3,-ic2,ibb,it,      'BR(~chi_40 -> ~c
     .hi_2- bb      t)'
      endif
      if(brchelne(4,1).ne.0.D0) then
      write(nout,104) brchelne(4,1),3,ic1,-ine,ie,     'BR(~chi_40 -> ~c
     .hi_1+ nu_eb   e-)'
      write(nout,104) brchelne(4,1),3,-ic1,ine,-ie,    'BR(~chi_40 -> ~c
     .hi_1- nu_e    e+)'
      endif
      if(brchelne(4,2).ne.0.D0) then
      write(nout,104) brchelne(4,2),3,ic2,-ine,ie,     'BR(~chi_40 -> ~c
     .hi_2+ nu_eb   e-)'
      write(nout,104) brchelne(4,2),3,-ic2,ine,-ie,    'BR(~chi_40 -> ~c
     .hi_2- nu_e    e+)'
      endif
      if(brchmunmu(4,1).ne.0.D0) then
      write(nout,104) brchmunmu(4,1),3,ic1,-inmu,imu,  'BR(~chi_40 -> ~c
     .hi_1+ nu_mub  mu-)'
      write(nout,104) brchmunmu(4,1),3,-ic1,inmu,-imu, 'BR(~chi_40 -> ~c
     .hi_1- nu_mu   mu+)'
      endif
      if(brchmunmu(4,2).ne.0.D0) then
      write(nout,104) brchmunmu(4,2),3,ic2,-inmu,imu,  'BR(~chi_40 -> ~c
     .hi_2+ nu_mub  mu-)'
      write(nout,104) brchmunmu(4,2),3,-ic2,inmu,-imu, 'BR(~chi_40 -> ~c
     .hi_2- nu_mu   mu+)'
      endif
      if(brchtauntau(4,1).ne.0.D0) then
      write(nout,104) brchtauntau(4,1),3,ic1,-intau,itau, 'BR(~chi_40 ->
     . ~chi_1+ nu_taub tau-)'
      write(nout,104) brchtauntau(4,1),3,-ic1,intau,-itau,'BR(~chi_40 ->
     . ~chi_1- nu_tau  tau+)'
      endif
      if(brchtauntau(4,2).ne.0.D0) then
      write(nout,104) brchtauntau(4,2),3,ic2,-intau,itau, 'BR(~chi_40 ->
     . ~chi_2+ nu_taub tau-)'
      write(nout,104) brchtauntau(4,2),3,-ic2,intau,-itau,'BR(~chi_40 ->
     . ~chi_2- nu_tau  tau+)'
      endif
      if(brglup(4).ne.0.D0) then
      write(nout,104) brglup(4),3,iglo,iub,iu,            'BR(~chi_40 ->
     . ~g      ub      u)'
      endif
      if(brgldo(4).ne.0.D0) then
      write(nout,104) brgldo(4),3,iglo,idb,id,            'BR(~chi_40 ->
     . ~g      db      d)'
      endif
      if(brglch(4).ne.0.D0) then
      write(nout,104) brglch(4),3,iglo,icb,ic,            'BR(~chi_40 ->
     . ~g      cb      c)'
      endif
      if(brglst(4).ne.0.D0) then
      write(nout,104) brglst(4),3,iglo,isb,is,            'BR(~chi_40 ->
     . ~g      sb      s)'
      endif
      if(brgltop(4).ne.0.D0) then
      write(nout,104) brgltop(4),3,iglo,itb,it,           'BR(~chi_40 ->
     . ~g      tb      t)'
      endif
      if(brglbot(4).ne.0.D0)  then
      write(nout,104) brglbot(4),3,iglo,ibb,ib,           'BR(~chi_40 ->
     . ~g      bb      b)'
      endif
      endif
     
c ------------------ c
c neutralino5 decays c
c ------------------ c
      write(nout,99)
      write(nout,100) 1000045,neuttot(5),'neutralino5'
      write(nout,49) 'neutralino5 2-body decays'
      if(neuttot2(5).ne.0.D0) then
      write(nout,101)
      if(brneutzneut(5,1).ne.0.D0) then
      write(nout,102) brneutzneut(5,1),2,in1,iz,       'BR(~chi_50 -> ~c
     .hi_10   Z )'
      endif
      if(brneutzneut(5,2).ne.0.D0) then
      write(nout,102) brneutzneut(5,2),2,in2,iz,       'BR(~chi_50 -> ~c
     .hi_20   Z )'
      endif
      if(brneutzneut(5,3).ne.0.D0) then
      write(nout,102) brneutzneut(5,3),2,in3,iz,       'BR(~chi_50 -> ~c
     .hi_30   Z )'
      endif
c
      if(brneutzneut(5,4).ne.0.D0) then
      write(nout,102) brneutzneut(5,4),2,in4,iz,       'BR(~chi_50 -> ~c
     .hi_40   Z )'
      endif
      if(brneutwchar(5,1).ne.0.D0) then
      write(nout,102) brneutwchar(5,1),2,ic1,-iwc,     'BR(~chi_50 -> ~c
     .hi_1+   W-)'
      write(nout,102) brneutwchar(5,1),2,-ic1,iwc,     'BR(~chi_50 -> ~c
     .hi_1-   W+)'
      endif
      if(brneutwchar(5,2).ne.0.D0) then
      write(nout,102) brneutwchar(5,2),2,ic2,-iwc,     'BR(~chi_50 -> ~c
     .hi_2+   W-)'
      write(nout,102) brneutwchar(5,2),2,-ic2,iwc,     'BR(~chi_50 -> ~c
     .hi_2-   W+)'
      endif
c
      if(brneutHneut(5,1,1).ne.0.D0) then
      write(nout,102) brneutHneut(5,1,1),2,in1,ihH1,   'BR(~chi_50 -> ~c
     .hi_10   H_1 )'
      endif
      if(brneutHneut(5,1,2).ne.0.D0) then
      write(nout,102) brneutHneut(5,1,2),2,in1,ihH2,   'BR(~chi_50 -> ~c
     .hi_10   H_2 )'
      endif
      if(brneutHneut(5,1,3).ne.0.D0) then
      write(nout,102) brneutHneut(5,1,3),2,in1,ihH3,   'BR(~chi_50 -> ~c
     .hi_10   H_3 )'
      endif
      if(brneutAneut(5,1,1).ne.0.D0) then
      write(nout,102) brneutAneut(5,1,1),2,in1,ihA1,   'BR(~chi_50 -> ~c
     .hi_10   A_1 )'
      endif
      if(brneutAneut(5,1,2).ne.0.D0) then
      write(nout,102) brneutAneut(5,1,2),2,in1,ihA2,   'BR(~chi_50 -> ~c
     .hi_10   A_2 )'
      endif
c
      if(brneutHneut(5,2,1).ne.0.D0) then
      write(nout,102) brneutHneut(5,2,1),2,in2,ihH1,   'BR(~chi_50 -> ~c
     .hi_20   H_1 )'
      endif
      if(brneutHneut(5,2,2).ne.0.D0) then
      write(nout,102) brneutHneut(5,2,2),2,in2,ihH2,   'BR(~chi_50 -> ~c
     .hi_20   H_2 )'
      endif
      if(brneutHneut(5,2,3).ne.0.D0) then
      write(nout,102) brneutHneut(5,2,3),2,in2,ihH3,   'BR(~chi_50 -> ~c
     .hi_20   H_3 )'
      endif
      if(brneutAneut(5,2,1).ne.0.D0) then
      write(nout,102) brneutAneut(5,2,1),2,in2,ihA1,   'BR(~chi_50 -> ~c
     .hi_20   A_1 )'
      endif
      if(brneutAneut(5,2,2).ne.0.D0) then
      write(nout,102) brneutAneut(5,2,2),2,in2,ihA2,   'BR(~chi_50 -> ~c
     .hi_20   A_2 )'
      endif
c
      if(brneutHneut(5,3,1).ne.0.D0) then
      write(nout,102) brneutHneut(5,3,1),2,in3,ihH1,   'BR(~chi_50 -> ~c
     .hi_30   H_1 )'
      endif
      if(brneutHneut(5,3,2).ne.0.D0) then
      write(nout,102) brneutHneut(5,3,2),2,in3,ihH2,   'BR(~chi_50 -> ~c
     .hi_30   H_2 )'
      endif
      if(brneutHneut(5,3,3).ne.0.D0) then
      write(nout,102) brneutHneut(5,3,3),2,in3,ihH3,   'BR(~chi_50 -> ~c
     .hi_30   H_3 )'
      endif
      if(brneutAneut(5,3,1).ne.0.D0) then
      write(nout,102) brneutAneut(5,3,1),2,in3,ihA1,   'BR(~chi_50 -> ~c
     .hi_30   A_1 )'
      endif
      if(brneutAneut(5,3,2).ne.0.D0) then
      write(nout,102) brneutAneut(5,3,2),2,in3,ihA2,   'BR(~chi_50 -> ~c
     .hi_30   A_2 )'
      endif
c
      if(brneutHneut(5,4,1).ne.0.D0) then
      write(nout,102) brneutHneut(5,4,1),2,in4,ihH1,   'BR(~chi_50 -> ~c
     .hi_40   H_1 )'
      endif
      if(brneutHneut(5,4,2).ne.0.D0) then
      write(nout,102) brneutHneut(5,4,2),2,in4,ihH2,   'BR(~chi_50 -> ~c
     .hi_40   H_2 )'
      endif
      if(brneutHneut(5,4,3).ne.0.D0) then
      write(nout,102) brneutHneut(5,4,3),2,in4,ihH3,   'BR(~chi_50 -> ~c
     .hi_40   H_3 )'
      endif
      if(brneutAneut(5,4,1).ne.0.D0) then
      write(nout,102) brneutAneut(5,4,1),2,in4,ihA1,   'BR(~chi_50 -> ~c
     .hi_40   A_1 )'
      endif
      if(brneutAneut(5,4,2).ne.0.D0) then
      write(nout,102) brneutAneut(5,4,2),2,in4,ihA2,   'BR(~chi_50 -> ~c
     .hi_40   A_2 )'
      endif

      if(brneuthcchar(5,1).ne.0.D0) then
      write(nout,102) brneuthcchar(5,1),2,ic1,-ihc,    'BR(~chi_50 -> ~c
     .hi_1+   H-)'
      write(nout,102) brneuthcchar(5,1),2,-ic1,ihc,    'BR(~chi_50 -> ~c
     .hi_1-   H+)'
      endif
      if(brneuthcchar(5,2).ne.0.D0) then
      write(nout,102) brneuthcchar(5,2),2,ic2,-ihc,    'BR(~chi_50 -> ~c
     .hi_2+   H-)'
      write(nout,102) brneuthcchar(5,2),2,-ic2,ihc,    'BR(~chi_50 -> ~c
     .hi_2-   H+)'
      endif
      if(brneutsupl(5).ne.0.D0) then
      write(nout,102) brneutsupl(5),2,isul,iub,        'BR(~chi_50 -> ~u
     ._L      ub)'
      write(nout,102) brneutsupl(5),2,-isul,iu,        'BR(~chi_50 -> ~u
     ._L*     u )'
      endif
      if(brneutsupr(5).ne.0.D0) then
      write(nout,102) brneutsupr(5),2,isur,iub,        'BR(~chi_50 -> ~u
     ._R      ub)'
      write(nout,102) brneutsupr(5),2,-isur,iu,        'BR(~chi_50 -> ~u
     ._R*     u )'
      endif
      if(brneutsdownl(5).ne.0.D0) then
      write(nout,102) brneutsdownl(5),2,isdl,idb,      'BR(~chi_50 -> ~d
     ._L      db)'
      write(nout,102) brneutsdownl(5),2,-isdl,id,      'BR(~chi_50 -> ~d
     ._L*     d )'
      endif
      if(brneutsdownr(5).ne.0.D0) then
      write(nout,102) brneutsdownr(5),2,isdr,idb,      'BR(~chi_50 -> ~d
     ._R      db)'
      write(nout,102) brneutsdownr(5),2,-isdr,id,      'BR(~chi_50 -> ~d
     ._R*     d )'
      endif
      if(brneutsupl(5).ne.0.D0) then
      write(nout,102) brneutsupl(5),2,iscl,icb,        'BR(~chi_50 -> ~c
     ._L      cb)'
      write(nout,102) brneutsupl(5),2,-iscl,ic,        'BR(~chi_50 -> ~c
     ._L*     c )'
      endif
      if(brneutsupr(5).ne.0.D0) then
      write(nout,102) brneutsupr(5),2,iscr,icb,        'BR(~chi_50 -> ~c
     ._R      cb)'
      write(nout,102) brneutsupr(5),2,-iscr,ic,        'BR(~chi_50 -> ~c
     ._R*     c )'
      endif
      if(brneutsdownl(5).ne.0.D0) then
      write(nout,102) brneutsdownl(5),2,issl,isb,      'BR(~chi_50 -> ~s
     ._L      sb)'
      write(nout,102) brneutsdownl(5),2,-issl,is,      'BR(~chi_50 -> ~s
     ._L*     s )'
      endif
      if(brneutsdownr(5).ne.0.D0) then
      write(nout,102) brneutsdownr(5),2,issr,isb,      'BR(~chi_50 -> ~s
     ._R      sb)'
      write(nout,102) brneutsdownr(5),2,-issr,is,      'BR(~chi_50 -> ~s
     ._R*     s )'
      endif
      if(brneutst1(5).ne.0.D0) then
      write(nout,102) brneutst1(5),2,ist1,itb,         'BR(~chi_50 -> ~t
     ._1      tb)'
      write(nout,102) brneutst1(5),2,-ist1,it,         'BR(~chi_50 -> ~t
     ._1*     t )'
      endif
      if(brneutst2(5).ne.0.D0) then
      write(nout,102) brneutst2(5),2,ist2,itb,         'BR(~chi_50 -> ~t
     ._2      tb)'
      write(nout,102) brneutst2(5),2,-ist2,it,         'BR(~chi_50 -> ~t
     ._2*     t )'
      endif
      if(brneutsb1(5).ne.0.D0) then
      write(nout,102) brneutsb1(5),2,isb1,ibb,         'BR(~chi_50 -> ~b
     ._1      bb)'
      write(nout,102) brneutsb1(5),2,-isb1,ib,         'BR(~chi_50 -> ~b
     ._1*     b )'
      endif
      if(brneutsb2(5).ne.0.D0) then
      write(nout,102) brneutsb2(5),2,isb2,ibb,         'BR(~chi_50 -> ~b
     ._2      bb)'
      write(nout,102) brneutsb2(5),2,-isb2,ib,         'BR(~chi_50 -> ~b
     ._2*     b )'
      endif
      if(brneutsell(5).ne.0.D0) then
      write(nout,102) brneutsell(5),2,isell,-ie,       'BR(~chi_50 -> ~e
     ._L-     e+)'
      write(nout,102) brneutsell(5),2,-isell,ie,       'BR(~chi_50 -> ~e
     ._L+     e-)'
      endif
      if(brneutselr(5).ne.0.D0) then
      write(nout,102) brneutselr(5),2,iselr,-ie,       'BR(~chi_50 -> ~e
     ._R-     e+)'
      write(nout,102) brneutselr(5),2,-iselr,ie,       'BR(~chi_50 -> ~e
     ._R+     e-)'
      endif
      if(brneutsell(5).ne.0.D0) then
      write(nout,102) brneutsell(5),2,ismul,-imu,      'BR(~chi_50 -> ~m
     .u_L-    mu+)' 
      write(nout,102) brneutsell(5),2,-ismul,imu,      'BR(~chi_50 -> ~m
     .u_L+    mu-)' 
      endif
      if(brneutselr(5).ne.0.D0) then      
      write(nout,102) brneutselr(5),2,ismur,-imu,      'BR(~chi_50 -> ~m
     .u_R-    mu+)' 
      write(nout,102) brneutselr(5),2,-ismur,imu,      'BR(~chi_50 -> ~m
     .u_R+    mu-)'
      endif
      if(brneutstau1(5).ne.0.D0) then
      write(nout,102) brneutstau1(5),2,istau1,-itau,   'BR(~chi_50 -> ~t
     .au_1-   tau+)'
      write(nout,102) brneutstau1(5),2,-istau1,itau,   'BR(~chi_50 -> ~t
     .au_1+   tau-)'
      endif
      if(brneutstau2(5).ne.0.D0) then
      write(nout,102) brneutstau2(5),2,istau2,-itau,   'BR(~chi_50 -> ~t
     .au_2-   tau+)'
      write(nout,102) brneutstau2(5),2,-istau2,itau,   'BR(~chi_50 -> ~t
     .au_2+   tau-)'
      endif
      if(brneutsnel(5).ne.0.D0) then
      write(nout,102) brneutsnel(5),2,inel,-ine,       'BR(~chi_50 -> ~n
     .u_eL    nu_eb)'
      write(nout,102) brneutsnel(5),2,-inel,ine,       'BR(~chi_50 -> ~n
     .u_eL*   nu_e )'
      write(nout,102) brneutsnel(5),2,inmul,-inmu,     'BR(~chi_50 -> ~n
     .u_muL   nu_mub)'
      write(nout,102) brneutsnel(5),2,-inmul,inmu,     'BR(~chi_50 -> ~n
     .u_muL*  nu_mu )'
      endif
      if(brneutsn1(5).ne.0.D0) then
      write(nout,102) brneutsn1(5),2,intau1,-intau,    'BR(~chi_50 -> ~n
     .u_tau1  nu_taub)'
      write(nout,102) brneutsn1(5),2,-intau1,intau,    'BR(~chi_50 -> ~n
     .u_tau1* nu_tau )'
      endif
      endif
C     =============================
*     NEUTRALINO5 : Radiative decay
C     =============================

            if(flagloop.eq.1.D0) then
      if(brnraddec(5,1).ne.0.D0) then
      write(nout,102) brnraddec(5,1),2,in1,iga,        'BR(~chi_50 -> ~c
     .hi_10 gam)'
      endif
      if(brnraddec(5,2).ne.0.D0) then
      write(nout,102) brnraddec(5,2),2,in2,iga,        'BR(~chi_50 -> ~c
     .hi_20 gam)'
      endif
      if(brnraddec(5,3).ne.0.D0) then
      write(nout,102) brnraddec(5,3),2,in3,iga,        'BR(~chi_50 -> ~c
     .hi_30 gam)'
      endif
      if(brnraddec(5,4).ne.0.D0) then
      write(nout,102) brnraddec(5,4),2,in4,iga,        'BR(~chi_50 -> ~c
     .hi_40 gam)'
      endif
            endif
C     =============================
*     NEUTRALINO5 : Three body
C     =============================

      if(neuttot3(5).ne.0.D0) then 
      write(nout,49) 'neutralino5 3-body decays'
      write(nout,103)
      if(brneutup(5,1).ne.0.D0) then
      write(nout,104) brneutup(5,1),3,in1,iub,iu,      'BR(~chi_50 -> ~c
     .hi_10 ub      u)'
      endif
      if(brneutdow(5,1).ne.0.D0) then
      write(nout,104) brneutdow(5,1),3,in1,idb,id,     'BR(~chi_50 -> ~c
     .hi_10 db      d)'
      endif
      if(brneutch(5,1).ne.0.D0) then
      write(nout,104) brneutch(5,1),3,in1,icb,ic,      'BR(~chi_50 -> ~c
     .hi_10 cb      c)'
      endif
      if(brneutst(5,1).ne.0.D0) then
      write(nout,104) brneutst(5,1),3,in1,isb,is,      'BR(~chi_50 -> ~c
     .hi_10 sb      s)'
      endif
      if(brneuttop(5,1).ne.0.D0) then
      write(nout,104) brneuttop(5,1),3,in1,itb,it,     'BR(~chi_50 -> ~c
     .hi_10 tb      t)'
      endif
      if(brneutbot(5,1).ne.0.D0) then
      write(nout,104) brneutbot(5,1),3,in1,ibb,ib,     'BR(~chi_50 -> ~c
     .hi_10 bb      b)'
      endif
      if(brneutel(5,1).ne.0.D0) then
      write(nout,104) brneutel(5,1),3,in1,-ie,ie,      'BR(~chi_50 -> ~c
     .hi_10 e+      e-)'
      endif
      if(brneutmu(5,1).ne.0.D0) then
      write(nout,104) brneutmu(5,1),3,in1,-imu,imu,    'BR(~chi_50 -> ~c
     .hi_10 mu+     mu-)'
      endif
      if(brneuttau(5,1).ne.0.D0) then
      write(nout,104) brneuttau(5,1),3,in1,-itau,itau, 'BR(~chi_50 -> ~c
     .hi_10 tau+    tau-)'
      endif
      if(brneutnue(5,1).ne.0.D0) then
      write(nout,104) brneutnue(5,1),3,in1,-ine,ine,   'BR(~chi_50 -> ~c
     .hi_10 nu_eb   nu_e)'
      endif
      if(brneutnumu(5,1).ne.0.D0) then
      write(nout,104) brneutnumu(5,1),3,in1,-inmu,inmu,'BR(~chi_50 -> ~c
     .hi_10 nu_mub  nu_mu)'
      endif
      if(brneutnutau(5,1).ne.0.D0) then
      write(nout,104) brneutnutau(5,1),3,in1,-intau,intau,'BR(~chi_50 ->
     . ~chi_10 nu_taub nu_tau)'
      endif
      if(brneutup(5,2).ne.0.D0) then
      write(nout,104) brneutup(5,2),3,in2,iub,iu,      'BR(~chi_50 -> ~c
     .hi_20 ub      u)'
      endif
      if(brneutdow(5,2).ne.0.D0) then
      write(nout,104) brneutdow(5,2),3,in2,idb,id,     'BR(~chi_50 -> ~c
     .hi_20 db      d)'
      endif
      if(brneutch(5,2).ne.0.D0) then
      write(nout,104) brneutch(5,2),3,in2,icb,ic,      'BR(~chi_50 -> ~c
     .hi_20 cb      c)'
      endif
      if(brneutst(5,2).ne.0.D0) then
      write(nout,104) brneutst(5,2),3,in2,isb,is,      'BR(~chi_50 -> ~c
     .hi_20 sb      s)'
      endif
      if(brneuttop(5,2).ne.0.D0) then
      write(nout,104) brneuttop(5,2),3,in2,itb,it,     'BR(~chi_50 -> ~c
     .hi_20 tb      t)'
      endif
      if(brneutbot(5,2).ne.0.D0) then
      write(nout,104) brneutbot(5,2),3,in2,ibb,ib,     'BR(~chi_50 -> ~c
     .hi_20 bb      b)'
      endif
      if(brneutel(5,2).ne.0.D0) then
      write(nout,104) brneutel(5,2),3,in2,-ie,ie,      'BR(~chi_50 -> ~c
     .hi_20 e+      e-)'
      endif
      if(brneutmu(5,2).ne.0.D0) then
      write(nout,104) brneutmu(5,2),3,in2,-imu,imu,    'BR(~chi_50 -> ~c
     .hi_20 mu+     mu-)'
      endif
      if(brneuttau(5,2).ne.0.D0) then
      write(nout,104) brneuttau(5,2),3,in2,-itau,itau, 'BR(~chi_50 -> ~c
     .hi_20 tau+    tau-)'
      endif
      if(brneutnue(5,2).ne.0.D0) then
      write(nout,104) brneutnue(5,2),3,in2,-ine,ine,   'BR(~chi_50 -> ~c
     .hi_20 nu_eb   nu_e)'
      endif
      if(brneutnumu(5,2).ne.0.D0) then
      write(nout,104) brneutnumu(5,2),3,in2,-inmu,inmu,'BR(~chi_50 -> ~c
     .hi_20 nu_mub  nu_mu)'
      endif
      if(brneutnutau(5,2).ne.0.D0) then
      write(nout,104) brneutnutau(5,2),3,in2,-intau,intau,'BR(~chi_50 ->
     . ~chi_20 nu_taub nu_tau)'
      endif
      if(brneutup(5,3).ne.0.D0) then
      write(nout,104) brneutup(5,3),3,in3,iub,iu,      'BR(~chi_50 -> ~c
     .hi_30 ub      u)'
      endif
      if(brneutdow(5,3).ne.0.D0) then
      write(nout,104) brneutdow(5,3),3,in3,idb,id,     'BR(~chi_50 -> ~c
     .hi_30 db      d)'
      endif
      if(brneutch(5,3).ne.0.D0) then
      write(nout,104) brneutch(5,3),3,in3,icb,ic,      'BR(~chi_50 -> ~c
     .hi_30 cb      c)'
      endif
      if(brneutst(5,3).ne.0.D0) then
      write(nout,104) brneutst(5,3),3,in3,isb,is,      'BR(~chi_50 -> ~c
     .hi_30 sb      s)'
      endif
      if(brneuttop(5,3).ne.0.D0) then
      write(nout,104) brneuttop(5,3),3,in3,itb,it,     'BR(~chi_50 -> ~c
     .hi_30 tb      t)'
      endif
      if(brneutbot(5,3).ne.0.D0) then
      write(nout,104) brneutbot(5,3),3,in3,ibb,ib,     'BR(~chi_50 -> ~c
     .hi_30 bb      b)'
      endif
      if(brneutel(5,3).ne.0.D0) then
      write(nout,104) brneutel(5,3),3,in3,-ie,ie,      'BR(~chi_50 -> ~c
     .hi_30 e+      e-)'
      endif
      if(brneutmu(5,3).ne.0.D0) then
      write(nout,104) brneutmu(5,3),3,in3,-imu,imu,    'BR(~chi_50 -> ~c
     .hi_30 mu+     mu-)'
      endif
      if(brneuttau(5,3).ne.0.D0) then
      write(nout,104) brneuttau(5,3),3,in3,-itau,itau, 'BR(~chi_50 -> ~c
     .hi_30 tau+    tau-)'
      endif
      if(brneutnue(5,3).ne.0.D0) then
      write(nout,104) brneutnue(5,3),3,in3,-ine,ine,   'BR(~chi_50 -> ~c
     .hi_30 nu_eb   nu_e)'
      endif
      if(brneutnumu(5,3).ne.0.D0) then
      write(nout,104) brneutnumu(5,3),3,in3,-inmu,inmu,'BR(~chi_50 -> ~c
     .hi_30 nu_mub  nu_mu)'
      endif
      if(brneutnutau(5,3).ne.0.D0) then
      write(nout,104) brneutnutau(5,3),3,in3,-intau,intau,'BR(~chi_50 ->
     . ~chi_30 nu_taub nu_tau)'
      endif

      if(brneutup(5,4).ne.0.D0) then
      write(nout,104) brneutup(5,4),3,in4,iub,iu,      'BR(~chi_50 -> ~c
     .hi_40 ub      u)'
      endif
      if(brneutdow(5,4).ne.0.D0) then
      write(nout,104) brneutdow(5,4),3,in4,idb,id,     'BR(~chi_50 -> ~c
     .hi_40 db      d)'
      endif
      if(brneutch(5,4).ne.0.D0) then
      write(nout,104) brneutch(5,4),3,in4,icb,ic,      'BR(~chi_50 -> ~c
     .hi_40 cb      c)'
      endif
      if(brneutst(5,4).ne.0.D0) then
      write(nout,104) brneutst(5,4),3,in4,isb,is,      'BR(~chi_50 -> ~c
     .hi_40 sb      s)'
      endif
      if(brneuttop(5,4).ne.0.D0) then
      write(nout,104) brneuttop(5,4),3,in4,itb,it,     'BR(~chi_50 -> ~c
     .hi_40 tb      t)'
      endif
      if(brneutbot(5,4).ne.0.D0) then
      write(nout,104) brneutbot(5,4),3,in4,ibb,ib,     'BR(~chi_50 -> ~c
     .hi_40 bb      b)'
      endif
      if(brneutel(5,4).ne.0.D0) then
      write(nout,104) brneutel(5,4),3,in4,-ie,ie,      'BR(~chi_50 -> ~c
     .hi_40 e+      e-)'
      endif
      if(brneutmu(5,4).ne.0.D0) then
      write(nout,104) brneutmu(5,4),3,in4,-imu,imu,    'BR(~chi_50 -> ~c
     .hi_40 mu+     mu-)'
      endif
      if(brneuttau(5,4).ne.0.D0) then
      write(nout,104) brneuttau(5,4),3,in4,-itau,itau, 'BR(~chi_50 -> ~c
     .hi_40 tau+    tau-)'
      endif
      if(brneutnue(5,4).ne.0.D0) then
      write(nout,104) brneutnue(5,4),3,in4,-ine,ine,   'BR(~chi_50 -> ~c
     .hi_40 nu_eb   nu_e)'
      endif
      if(brneutnumu(5,4).ne.0.D0) then
      write(nout,104) brneutnumu(5,4),3,in4,-inmu,inmu,'BR(~chi_50 -> ~c
     .hi_40 nu_mub  nu_mu)'
      endif
      if(brneutnutau(5,4).ne.0.D0) then
      write(nout,104) brneutnutau(5,4),3,in4,-intau,intau,'BR(~chi_50 ->
     . ~chi_40 nu_taub nu_tau)'
      endif
*    
      if(brchubd(5,1).ne.0.D0) then
      write(nout,104) brchubd(5,1),3,ic1,iub,id,       'BR(~chi_50 -> ~c
     .hi_1+ ub      d)'
      write(nout,104) brchubd(5,1),3,-ic1,idb,iu,      'BR(~chi_50 -> ~c
     .hi_1- db      u)'
      endif
      if(brchubd(5,2).ne.0.D0) then
      write(nout,104) brchubd(5,2),3,ic2,iub,id,       'BR(~chi_50 -> ~c
     .hi_2+ ub      d)'
      write(nout,104) brchubd(5,2),3,-ic2,idb,iu,      'BR(~chi_50 -> ~c
     .hi_2- db      u)'
      endif
      if(brchcbs(5,1).ne.0.D0) then
      write(nout,104) brchcbs(5,1),3,ic1,icb,is,       'BR(~chi_50 -> ~c
     .hi_1+ cb      s)'
      write(nout,104) brchcbs(5,1),3,-ic1,isb,ic,      'BR(~chi_50 -> ~c
     .hi_1- sb      c)'
      endif
      if(brchcbs(5,2).ne.0.D0) then
      write(nout,104) brchcbs(5,2),3,ic2,icb,is,       'BR(~chi_50 -> ~c
     .hi_2+ cb      s)'
      write(nout,104) brchcbs(5,2),3,-ic2,isb,ic,      'BR(~chi_50 -> ~c
     .hi_2- sb      c)'
      endif
      if(brchtbb(5,1).ne.0.D0) then
      write(nout,104) brchtbb(5,1),3,ic1,itb,ib,       'BR(~chi_50 -> ~c
     .hi_1+ tb      b)'
      write(nout,104) brchtbb(5,1),3,-ic1,ibb,it,      'BR(~chi_50 -> ~c
     .hi_1- bb      t)'
      endif
      if(brchtbb(5,2).ne.0.D0) then
      write(nout,104) brchtbb(5,2),3,ic2,itb,ib,       'BR(~chi_50 -> ~c
     .hi_2+ tb      b)'
      write(nout,104) brchtbb(5,2),3,-ic2,ibb,it,      'BR(~chi_50 -> ~c
     .hi_2- bb      t)'
      endif
      if(brchelne(5,1).ne.0.D0) then
      write(nout,104) brchelne(5,1),3,ic1,-ine,ie,     'BR(~chi_50 -> ~c
     .hi_1+ nu_eb   e-)'
      write(nout,104) brchelne(5,1),3,-ic1,ine,-ie,    'BR(~chi_50 -> ~c
     .hi_1- nu_e    e+)'
      endif
      if(brchelne(5,2).ne.0.D0) then
      write(nout,104) brchelne(5,2),3,ic2,-ine,ie,     'BR(~chi_50 -> ~c
     .hi_2+ nu_eb   e-)'
      write(nout,104) brchelne(5,2),3,-ic2,ine,-ie,    'BR(~chi_50 -> ~c
     .hi_2- nu_e    e+)'
      endif
      if(brchmunmu(5,1).ne.0.D0) then
      write(nout,104) brchmunmu(5,1),3,ic1,-inmu,imu,  'BR(~chi_50 -> ~c
     .hi_1+ nu_mub  mu-)'
      write(nout,104) brchmunmu(5,1),3,-ic1,inmu,-imu, 'BR(~chi_50 -> ~c
     .hi_1- nu_mu   mu+)'
      endif
      if(brchmunmu(5,2).ne.0.D0) then
      write(nout,104) brchmunmu(5,2),3,ic2,-inmu,imu,  'BR(~chi_50 -> ~c
     .hi_2+ nu_mub  mu-)'
      write(nout,104) brchmunmu(5,2),3,-ic2,inmu,-imu, 'BR(~chi_50 -> ~c
     .hi_2- nu_mu   mu+)'
      endif
      if(brchtauntau(5,1).ne.0.D0) then
      write(nout,104) brchtauntau(5,1),3,ic1,-intau,itau, 'BR(~chi_50 ->
     . ~chi_1+ nu_taub tau-)'
      write(nout,104) brchtauntau(5,1),3,-ic1,intau,-itau,'BR(~chi_50 ->
     . ~chi_1- nu_tau  tau+)'
      endif
      if(brchtauntau(5,2).ne.0.D0) then
      write(nout,104) brchtauntau(5,2),3,ic2,-intau,itau, 'BR(~chi_50 ->
     . ~chi_2+ nu_taub tau-)'
      write(nout,104) brchtauntau(5,2),3,-ic2,intau,-itau,'BR(~chi_50 ->
     . ~chi_2- nu_tau  tau+)'
      endif
      if(brglup(5).ne.0.D0) then
      write(nout,104) brglup(5),3,iglo,iub,iu,            'BR(~chi_50 ->
     . ~g      ub      u)'
      endif
      if(brgldo(5).ne.0.D0) then
      write(nout,104) brgldo(5),3,iglo,idb,id,            'BR(~chi_50 ->
     . ~g      db      d)'
      endif
      if(brglch(5).ne.0.D0) then
      write(nout,104) brglch(5),3,iglo,icb,ic,            'BR(~chi_50 ->
     . ~g      cb      c)'
      endif
      if(brglst(5).ne.0.D0) then
      write(nout,104) brglst(5),3,iglo,isb,is,            'BR(~chi_50 ->
     . ~g      sb      s)'
      endif
      if(brgltop(5).ne.0.D0) then
      write(nout,104) brgltop(5),3,iglo,itb,it,           'BR(~chi_50 ->
     . ~g      tb      t)'
      endif
      if(brglbot(5).ne.0.D0)  then
      write(nout,104) brglbot(5),3,iglo,ibb,ib,           'BR(~chi_50 ->
     . ~g      bb      b)'
      endif
      endif
C 
c ------------- c
c Gluino decays c
c ------------- c
      write(nout,99)
      write(nout,100) iglo,gluitot,'gluino'
      if(gluitot2.ne.0.D0) then
      write(nout,49) 'gluino 2-body decays'
         write(nout,101)
      if(brgsdownl.ne.0.D0) then
      write(nout,102) brgsdownl,2,isdl,idb,'BR(~g -> ~d_L  db)'
      write(nout,102) brgsdownl,2,-isdl,id,'BR(~g -> ~d_L* d )'
      endif
      if(brgsdownr.ne.0.D0) then
      write(nout,102) brgsdownr,2,isdr,idb,'BR(~g -> ~d_R  db)'
      write(nout,102) brgsdownr,2,-isdr,id,'BR(~g -> ~d_R* d )'
      endif
      if(brgsupl.ne.0.D0) then
      write(nout,102) brgsupl,2,isul,iub  ,'BR(~g -> ~u_L  ub)'
      write(nout,102) brgsupl,2,-isul,iu  ,'BR(~g -> ~u_L* u )'
      endif
      if(brgsupr.ne.0.D0) then
      write(nout,102) brgsupr,2,isur,iub  ,'BR(~g -> ~u_R  ub)'
      write(nout,102) brgsupr,2,-isur,iu  ,'BR(~g -> ~u_R* u )'
      endif
      if(brgsdownl.ne.0.D0) then
      write(nout,102) brgsdownl,2,issl,isb,'BR(~g -> ~s_L  sb)'
      write(nout,102) brgsdownl,2,-issl,is,'BR(~g -> ~s_L* s )'
      endif
      if(brgsdownr.ne.0.D0) then
      write(nout,102) brgsdownr,2,issr,isb,'BR(~g -> ~s_R  sb)'
      write(nout,102) brgsdownr,2,-issr,is,'BR(~g -> ~s_R* s )'
      endif
      if(brgsupl.ne.0.D0) then
      write(nout,102) brgsupl,2,iscl,icb  ,'BR(~g -> ~c_L  cb)'
      write(nout,102) brgsupl,2,-iscl,ic  ,'BR(~g -> ~c_L* c )'
      endif
      if(brgsupr.ne.0.D0) then
      write(nout,102) brgsupr,2,iscr,icb  ,'BR(~g -> ~c_R  cb)'
      write(nout,102) brgsupr,2,-iscr,ic  ,'BR(~g -> ~c_R* c )'
      endif
      if(brgsb1.ne.0.D0) then
      write(nout,102) brgsb1,2,isb1,ibb   ,'BR(~g -> ~b_1  bb)'
      write(nout,102) brgsb1,2,-isb1,ib   ,'BR(~g -> ~b_1* b )'
      endif
      if(brgsb2.ne.0.D0) then
      write(nout,102) brgsb2,2,isb2,ibb   ,'BR(~g -> ~b_2  bb)'
      write(nout,102) brgsb2,2,-isb2,ib   ,'BR(~g -> ~b_2* b )'
      endif
      if(brgst1.ne.0.D0) then
      write(nout,102) brgst1,2,ist1,itb   ,'BR(~g -> ~t_1  tb)'
      write(nout,102) brgst1,2,-ist1,it   ,'BR(~g -> ~t_1* t )'
      endif
      if(brgst2.ne.0.D0) then
      write(nout,102) brgst2,2,ist2,itb   ,'BR(~g -> ~t_2  tb)'
      write(nout,102) brgst2,2,-ist2,it   ,'BR(~g -> ~t_2* t )'
      endif
      endif
C     =============================
*     Gluino: Radiative decay
C     =============================
      if(gluitotrad.ne.0.d0) then
            if(flagloop.eq.1.D0) then
      if(gluitot2.eq.0.D0) then
        write(nout,49) 'gluino 2-body decays'
        write(nout,101)
      endif
      if(brglnjgluon(1).ne.0.D0) then
      write(nout,102) brglnjgluon(1),2,in1,igl,'BR(~g -> ~chi_10 g)'
      endif
      if(brglnjgluon(2).ne.0.D0) then
      write(nout,102) brglnjgluon(2),2,in2,igl,'BR(~g -> ~chi_20 g)'
      endif
      if(brglnjgluon(3).ne.0.D0) then
      write(nout,102) brglnjgluon(3),2,in3,igl,'BR(~g -> ~chi_30 g)'
      endif
      if(brglnjgluon(4).ne.0.D0) then
      write(nout,102) brglnjgluon(4),2,in4,igl,'BR(~g -> ~chi_40 g)'
      endif
      if(brglnjgluon(5).ne.0.D0) then
      write(nout,102) brglnjgluon(5),2,in5,igl,'BR(~g -> ~chi_50 g)'
      endif
            endif
          endif
C     =============================
*     GLUINO: 3-body
C     =============================
      if(gluitot3.ne.0.D0) then
         write(nout,49) 'gluino 3-body decays'
      write(nout,103)
      if(brgodn(1).ne.0.D0) then
      write(nout,104) brgodn(1),3,in1,id,idb  ,'BR(~g -> ~chi_10 d  db)'
      endif
      if(brgodn(2).ne.0.D0) then
      write(nout,104) brgodn(2),3,in2,id,idb  ,'BR(~g -> ~chi_20 d  db)'
      endif
      if(brgodn(3).ne.0.D0) then
      write(nout,104) brgodn(3),3,in3,id,idb  ,'BR(~g -> ~chi_30 d  db)'
      endif
      if(brgodn(4).ne.0.D0) then
      write(nout,104) brgodn(4),3,in4,id,idb  ,'BR(~g -> ~chi_40 d  db)'
      endif
      if(brgodn(5).ne.0.D0) then
      write(nout,104) brgodn(5),3,in5,id,idb  ,'BR(~g -> ~chi_50 d  db)'
      endif
c     ----
      if(brgoup(1).ne.0.D0) then
      write(nout,104) brgoup(1),3,in1,iu,iub  ,'BR(~g -> ~chi_10 u  ub)'
      endif
      if(brgoup(2).ne.0.D0) then
      write(nout,104) brgoup(2),3,in2,iu,iub  ,'BR(~g -> ~chi_20 u  ub)'
      endif
      if(brgoup(3).ne.0.D0) then
      write(nout,104) brgoup(3),3,in3,iu,iub  ,'BR(~g -> ~chi_30 u  ub)'
      endif
      if(brgoup(4).ne.0.D0) then
      write(nout,104) brgoup(4),3,in4,iu,iub  ,'BR(~g -> ~chi_40 u  ub)'
      endif
      if(brgoup(5).ne.0.D0) then
      write(nout,104) brgoup(5),3,in5,iu,iub  ,'BR(~g -> ~chi_50 u  ub)'
      endif
c    ---
      if(brgodn(1).ne.0.D0) then
      write(nout,104) brgodn(1),3,in1,is,isb  ,'BR(~g -> ~chi_10 s  sb)'
      endif
      if(brgodn(2).ne.0.D0) then
      write(nout,104) brgodn(2),3,in2,is,isb  ,'BR(~g -> ~chi_20 s  sb)'
      endif
      if(brgodn(3).ne.0.D0) then
      write(nout,104) brgodn(3),3,in3,is,isb  ,'BR(~g -> ~chi_30 s  sb)'
      endif
      if(brgodn(4).ne.0.D0) then
      write(nout,104) brgodn(4),3,in4,is,isb  ,'BR(~g -> ~chi_40 s  sb)'
      endif
      if(brgodn(5).ne.0.D0) then
      write(nout,104) brgodn(5),3,in5,is,isb  ,'BR(~g -> ~chi_50 s  sb)'
      endif
c   ---
      if(brgoup(1).ne.0.D0) then
      write(nout,104) brgoup(1),3,in1,ic,icb  ,'BR(~g -> ~chi_10 c  cb)'
      endif
      if(brgoup(2).ne.0.D0) then
      write(nout,104) brgoup(2),3,in2,ic,icb  ,'BR(~g -> ~chi_20 c  cb)'
      endif
      if(brgoup(3).ne.0.D0) then
      write(nout,104) brgoup(3),3,in3,ic,icb  ,'BR(~g -> ~chi_30 c  cb)'
      endif
      if(brgoup(4).ne.0.D0) then
      write(nout,104) brgoup(4),3,in4,ic,icb  ,'BR(~g -> ~chi_40 c  cb)'
      endif
      if(brgoup(5).ne.0.D0) then
      write(nout,104) brgoup(5),3,in5,ic,icb  ,'BR(~g -> ~chi_50 c  cb)'
      endif
c    ------
      if(brgobt(1).ne.0.D0) then
      write(nout,104) brgobt(1),3,in1,ib,ibb  ,'BR(~g -> ~chi_10 b  bb)'
      endif
      if(brgobt(2).ne.0.D0) then
      write(nout,104) brgobt(2),3,in2,ib,ibb  ,'BR(~g -> ~chi_20 b  bb)'
      endif
      if(brgobt(3).ne.0.D0) then
      write(nout,104) brgobt(3),3,in3,ib,ibb  ,'BR(~g -> ~chi_30 b  bb)'
      endif
      if(brgobt(4).ne.0.D0) then
      write(nout,104) brgobt(4),3,in4,ib,ibb  ,'BR(~g -> ~chi_40 b  bb)'
      endif
      if(brgobt(5).ne.0.D0) then
      write(nout,104) brgobt(5),3,in5,ib,ibb  ,'BR(~g -> ~chi_50 b  bb)'
      endif
c    -----
      if(brgotp(1).ne.0.D0) then
      write(nout,104) brgotp(1),3,in1,it,itb  ,'BR(~g -> ~chi_10 t  tb)'
      endif
      if(brgotp(2).ne.0.D0) then
      write(nout,104) brgotp(2),3,in2,it,itb  ,'BR(~g -> ~chi_20 t  tb)'
      endif
      if(brgotp(3).ne.0.D0) then
      write(nout,104) brgotp(3),3,in3,it,itb  ,'BR(~g -> ~chi_30 t  tb)'
      endif
      if(brgotp(4).ne.0.D0) then
      write(nout,104) brgotp(4),3,in4,it,itb  ,'BR(~g -> ~chi_40 t  tb)'
      endif
      if(brgotp(5).ne.0.D0) then
      write(nout,104) brgotp(5),3,in5,it,itb  ,'BR(~g -> ~chi_50 t  tb)'
      endif
c    -----
      if(brgoud(1).ne.0.D0) then
      write(nout,104) brgoud(1),3,ic1,id,iub  ,'BR(~g -> ~chi_1+ d  ub)'
      write(nout,104) brgoud(1),3,-ic1,iu,idb ,'BR(~g -> ~chi_1- u  db)'
      endif
      if(brgoud(2).ne.0.D0) then
      write(nout,104) brgoud(2),3,ic2,id,iub  ,'BR(~g -> ~chi_2+ d  ub)'
      write(nout,104) brgoud(2),3,-ic2,iu,idb ,'BR(~g -> ~chi_2- u  db)'
      endif
      if(brgoud(1).ne.0.D0) then
      write(nout,104) brgoud(1),3,ic1,is,icb  ,'BR(~g -> ~chi_1+ s  cb)'
      write(nout,104) brgoud(1),3,-ic1,ic,isb ,'BR(~g -> ~chi_1- c  sb)'
      endif
      if(brgoud(2).ne.0.D0) then
      write(nout,104) brgoud(2),3,ic2,is,icb  ,'BR(~g -> ~chi_2+ s  cb)'
      write(nout,104) brgoud(2),3,-ic2,ic,isb ,'BR(~g -> ~chi_2- c  sb)'
      endif
      if(brgotb(1).ne.0.D0) then
      write(nout,104) brgotb(1),3,ic1,ib,itb  ,'BR(~g -> ~chi_1+ b  tb)'
      write(nout,104) brgotb(1),3,-ic1,it,ibb ,'BR(~g -> ~chi_1- t  bb)'
      endif
      if(brgotb(2).ne.0.D0) then
      write(nout,104) brgotb(2),3,ic2,ib,itb  ,'BR(~g -> ~chi_2+ b  tb)'
      write(nout,104) brgotb(2),3,-ic2,it,ibb ,'BR(~g -> ~chi_2- t  bb)'
      endif
      if(brwst1b.ne.0.D0) then
      write(nout,104) brwst1b,3,ist1,ibb,-iwc ,'BR(~g -> ~t_1    bb W-)'
      write(nout,104) brwst1b,3,-ist1,ib,iwc  ,'BR(~g -> ~t_1*   b  W+)'
      endif
      if(brhcst1b.ne.0.D0) then
      write(nout,104) brhcst1b,3,ist1,ibb,-ihc,'BR(~g -> ~t_1    bb H-)'
      write(nout,104) brhcst1b,3,-ist1,ib,ihc ,'BR(~g -> ~t_1*   b  H+)'
      endif
      endif

c ------------------ c
c Selectron_L decays c
c ------------------ c
      write(nout,99)
      write(nout,100) 1000011,selltot,'selectron_L'
      if(selltot2.ne.0.D0) then
      write(nout,49) 'selectron_L 2-body decays'
          write(nout,101)
      if(brsellneute(1).ne.0.D0) then
      write(nout,102) brsellneute(1),2,in1,ie,'BR(~e_L -> ~chi_10 e-)
     .'
      endif
      if(brsellneute(2).ne.0.D0) then
      write(nout,102) brsellneute(2),2,in2,ie,'BR(~e_L -> ~chi_20 e-)
     .'
      endif
      if(brsellneute(3).ne.0.D0) then
      write(nout,102) brsellneute(3),2,in3,ie, 'BR(~e_L -> ~chi_30 e-)
     .'
      endif
      if(brsellneute(4).ne.0.D0) then
      write(nout,102) brsellneute(4),2,in4,ie, 'BR(~e_L -> ~chi_40 e-)
     .'
      endif
      if(brsellneute(5).ne.0.D0) then
      write(nout,102) brsellneute(5),2,in5,ie, 'BR(~e_L -> ~chi_50 e-)
     .'
      endif
      if(brsellcharnue(1).ne.0.D0) then
      write(nout,102) brsellcharnue(1),2,-ic1,ine,'BR(~e_L -> ~chi_1- nu
     ._e)'
      endif
      if(brsellcharnue(2).ne.0.D0) then
      write(nout,102) brsellcharnue(2),2,-ic2,ine,'BR(~e_L -> ~chi_2- nu
     ._e)'
      endif
      endif
C-----------
* THREE BODY 
C-----------
      if(selltot3.ne.0.D0) then
      write(nout,49) 'selectron_L 3-body decays'
       write(nout,103)
      if(brsellstau1star.ne.0.D0) then
      write(nout,104) brsellstau1star,3,-istau1,ie,itau,
     .'BR(~e_L -> ~tau_1* e- tau)' 
      endif
      if(brsellstau1.ne.0.D0) then
      write(nout,104) brsellstau1,3,istau1,ie,-itau,
     .'BR(~e_L -> ~tau_1 e- taub)' 
      endif
      if(brsellstau1nutau.ne.0.D0) then
      write(nout,104)brsellstau1nutau,3,istau1,ine,-intau,
     . 'BR(~e_L -> ~tau_1 nu_e nu_taub)' 
      endif
      endif

c ------------------ c
c Selectron_R decays c
c ------------------ c
      write(nout,99)
      write(nout,100) 2000011,selrtot,'selectron_R'
      if(selrtot2.ne.0.D0) then
      write(nout,49) 'selectron_R 2-body decays'
        write(nout,101)
      if(brselrneute(1).ne.0.D0) then
      write(nout,102) brselrneute(1),2,in1,ie   ,'BR(~e_R -> ~chi_10 e-)
     .'
      endif
      if(brselrneute(2).ne.0.D0) then
      write(nout,102) brselrneute(2),2,in2,ie   ,'BR(~e_R -> ~chi_20 e-)
     .'
      endif
      if(brselrneute(3).ne.0.D0) then
      write(nout,102) brselrneute(3),2,in3,ie   ,'BR(~e_R -> ~chi_30 e-)
     .'
      endif
      if(brselrneute(4).ne.0.D0) then
      write(nout,102) brselrneute(4),2,in4,ie   ,'BR(~e_R -> ~chi_40 e-)
     .'
      endif
      if(brselrneute(5).ne.0.D0) then
      write(nout,102) brselrneute(5),2,in5,ie   ,'BR(~e_R -> ~chi_50 e-)
     .'
      endif
      if(brselrcharnue(1).ne.0.D0) then
      write(nout,102) brselrcharnue(1),2,-ic1,ine,'BR(~e_R -> ~chi_1- nu
     ._e)'
      endif
      if(brselrcharnue(2).ne.0.D0) then
      write(nout,102) brselrcharnue(2),2,-ic2,ine,'BR(~e_R -> ~chi_2- nu
     ._e)'
      endif
      endif
C-----------
* THREE BODY 
C-----------

      if(selrtot3.ne.0.D0) then
      write(nout,49) 'selectron_R 3-body decays'
       write(nout,103)
      if(brselrstaustar.ne.0.D0) then
      write(nout,104) brselrstaustar,3,-istau1,ie,-itau,
     .'BR(~e_R -> ~tau_1* e- taub)' 
      endif
      if(brselrstau.ne.0.D0) then
      write(nout,104) brselrstau,3,istau1,ie,itau,
     .'BR(~e_R -> ~tau_1 e- tau-)' 
      endif
      endif
      
c -------------- c
c Smuon_L decays c
c -------------- c
      write(nout,99)
      write(nout,100) 1000013,selltot,'smuon_L'
      if(selltot2.ne.0.D0) then
      write(nout,49) 'smuon_L 2-body decays'
          write(nout,101)
      if(brsellneute(1).ne.0.D0) then
      write(nout,102) brsellneute(1),2,in1,imu,'BR(~mu_L -> ~chi_10 mu-)
     .'
      endif
      if(brsellneute(2).ne.0.D0) then
      write(nout,102) brsellneute(2),2,in2,imu,'BR(~mu_L -> ~chi_20 mu-)
     .'
      endif
      if(brsellneute(3).ne.0.D0) then
      write(nout,102) brsellneute(3),2,in3,imu,'BR(~mu_L -> ~chi_30 mu-)
     .'
      endif
      if(brsellneute(4).ne.0.D0) then
      write(nout,102) brsellneute(4),2,in4,imu,'BR(~mu_L -> ~chi_40 mu-)
     .'
      endif
      if(brsellneute(5).ne.0.D0) then
      write(nout,102) brsellneute(5),2,in5,imu,'BR(~mu_L -> ~chi_50 mu-)
     .'
      endif
      if(brsellcharnue(1).ne.0.D0) then
      write(nout,102) brsellcharnue(1),2,-ic1,inmu,'BR(~mu_L -> ~chi_1- 
     .nu_mu)'
      endif
      if(brsellcharnue(2).ne.0.D0) then
      write(nout,102) brsellcharnue(2),2,-ic2,inmu,'BR(~mu_L -> ~chi_2- 
     .nu_mu)'
      endif
      endif
C-----------
* THREE BODY 
C-----------

            if(selltot3.ne.0.D0) then
      write(nout,49) 'smuon_L 3-body decays'
       write(nout,103)
      if(brsellstau1star.ne.0.D0) then
      write(nout,104) brsellstau1star,3,-istau1,imu,itau,
     .'BR(~e_L -> ~tau_1* mu- tau)' 
      endif
      if(brsellstau1.ne.0.D0) then
      write(nout,104) brsellstau1,3,istau1,imu,-itau,
     .'BR(~e_L -> ~tau_1 mu- taub)' 
      endif
      if(brsellstau1nutau.ne.0.D0) then
      write(nout,104)brsellstau1nutau,3,istau1,inmu,-intau,
     . 'BR(~e_L -> ~tau_1 nu_mu nu_taub)' 
      endif
      endif

c ------------------ c
c Smuon_R decays c
c ------------------ c
        write(nout,99)
        write(nout,100) 2000013,selrtot,'smuon_R'
      if(selrtot2.ne.0.D0) then
        write(nout,49) 'smuon_R 2-body decays'
        write(nout,101)
      if(brselrneute(1).ne.0.D0) then
      write(nout,102) brselrneute(1),2,in1,imu,'BR(~mu_R -> ~chi_10 mu-)
     .'
      endif
      if(brselrneute(2).ne.0.D0) then
      write(nout,102) brselrneute(2),2,in2,imu,'BR(~mu_R -> ~chi_20 mu-)
     .'
      endif
      if(brselrneute(3).ne.0.D0) then
      write(nout,102) brselrneute(3),2,in3,imu,'BR(~mu_R -> ~chi_30 mu-)
     .'
      endif
      if(brselrneute(4).ne.0.D0) then
      write(nout,102) brselrneute(4),2,in4,imu,'BR(~mu_R -> ~chi_40 mu-)
     .'
      endif
      if(brselrneute(5).ne.0.D0) then
      write(nout,102) brselrneute(5),2,in5,imu,'BR(~mu_R -> ~chi_50 mu-)
     .'
      endif
      if(brselrcharnue(1).ne.0.D0) then
      write(nout,102) brselrcharnue(1),2,-ic1,inmu,'BR(~mu_R -> ~chi_1- 
     .nu_mu)'
      endif
      if(brselrcharnue(2).ne.0.D0) then
      write(nout,102) brselrcharnue(2),2,-ic2,inmu,'BR(~mu_R -> ~chi_2- 
     .nu_mu)'
      endif
      endif
C-----------
* THREE BODY 
C-----------
      if(selrtot3.ne.0.D0) then
          write(nout,103)
      write(nout,49) 'smuon_R 3-body decays'
      if(brselrstaustar.ne.0.D0) then
      write(nout,104) brselrstaustar,3,-istau1,imu,-itau,
     .'BR(~mu_R -> ~tau_1* mu- taub)' 
      endif
      if(brselrstau.ne.0.D0) then
      write(nout,104) brselrstau,3,istau1,imu,itau,
     .'BR(~mu_R -> ~tau_1 mu- tau-)' 
      endif
      endif

c ------------- c
c Stau_1 decays c
c ------------- c
      write(nout,99)
      write(nout,100) 1000015,stau1tot2,'stau_1'
      if(stau1tot2.ne.0.D0) then
        write(nout,49) 'stau1 2-body decays'
        write(nout,101)
      if(brstau1neut(1).ne.0.D0) then
      write(nout,102) brstau1neut(1),2,in1,itau,  'BR(~tau_1 -> ~chi_10 
     . tau-)'
      endif
      if(brstau1neut(2).ne.0.D0) then
      write(nout,102) brstau1neut(2),2,in2,itau,  'BR(~tau_1 -> ~chi_20 
     . tau-)'
      endif
      if(brstau1neut(3).ne.0.D0) then
      write(nout,102) brstau1neut(3),2,in3,itau,  'BR(~tau_1 -> ~chi_30 
     . tau-)'
      endif
      if(brstau1neut(4).ne.0.D0) then
      write(nout,102) brstau1neut(4),2,in4,itau,  'BR(~tau_1 -> ~chi_40 
     . tau-)'
      endif
      if(brstau1neut(5).ne.0.D0) then
      write(nout,102) brstau1neut(5),2,in5,itau,  'BR(~tau_1 -> ~chi_50 
     . tau-)'
      endif
      if(brstau1char(1).ne.0.D0) then
      write(nout,102) brstau1char(1),2,-ic1,intau,'BR(~tau_1 -> ~chi_1- 
     . nu_tau)'
      endif
      if(brstau1char(2).ne.0.D0) then
      write(nout,102) brstau1char(2),2,-ic2,intau,'BR(~tau_1 -> ~chi_2- 
     . nu_tau)'
      endif
      if(brstau1hcsn(1).ne.0.D0) then
      write(nout,102) brstau1hcsn(1),2,intau1,-ihc,'BR(~tau_1 -> ~nu_tau
     .L H-)'
      endif
      if(brstau1wsn(1).ne.0.D0) then
      write(nout,102) brstau1wsn(1),2,intau1,-iwc,'BR(~tau_1 -> ~nu_tauL
     . W-)'
      endif
      endif

c ------------- c
c Stau_2 decays c
c ------------- c
      write(nout,99)
      write(nout,100) 2000015,stau2tot,'stau_2'
      if(stau2tot2.ne.0.D0) then
        write(nout,49) 'stau_2 2-body decays'
        write(nout,101)
      if(brstau2neut(1).ne.0.D0) then
      write(nout,102) brstau2neut(1),2,in1,itau,  'BR(~tau_2 -> ~chi_10 
     . tau-)'
      endif
      if(brstau2neut(2).ne.0.D0) then
      write(nout,102) brstau2neut(2),2,in2,itau,  'BR(~tau_2 -> ~chi_20 
     . tau-)'
      endif
      if(brstau2neut(3).ne.0.D0) then
      write(nout,102) brstau2neut(3),2,in3,itau,  'BR(~tau_2 -> ~chi_30 
     . tau-)'
      endif
      if(brstau2neut(4).ne.0.D0) then
      write(nout,102) brstau2neut(4),2,in4,itau,  'BR(~tau_2 -> ~chi_40 
     . tau-)'
      endif
      if(brstau2neut(5).ne.0.D0) then
      write(nout,102) brstau2neut(5),2,in5,itau,  'BR(~tau_2 -> ~chi_50 
     . tau-)'
      endif
      if(brstau2char(1).ne.0.D0) then
      write(nout,102) brstau2char(1),2,-ic1,intau,'BR(~tau_2 -> ~chi_1- 
     . nu_tau)'
      endif
      if(brstau2char(2).ne.0.D0) then
      write(nout,102) brstau2char(2),2,-ic2,intau,'BR(~tau_2 -> ~chi_2- 
     . nu_tau)'
      endif
      if(brstau2hcsn(1).ne.0.D0) then
      write(nout,102) brstau2hcsn(1),2,intau1,-ihc,'BR(~tau_2 -> ~nu_tau
     .L H-)'
      endif
      if(brstau2wsn(1).ne.0.D0) then
      write(nout,102) brstau2wsn(1),2,intau1,-iwc,'BR(~tau_2 -> ~nu_tauL
     . W-)'
      endif
      if(brstau2H(1).ne.0.D0) then
      write(nout,102) brstau2H(1),2,istau1,ihH1,  'BR(~tau_2 -> ~tau_1
     . H_1)'
      endif
      if(brstau2H(2).ne.0.D0) then
      write(nout,102) brstau2H(2),2,istau1,ihH2,  'BR(~tau_2 -> ~tau_1
     . H_2)'
      endif
      if(brstau2H(3).ne.0.D0) then
      write(nout,102) brstau2H(3),2,istau1,ihH3,  'BR(~tau_2 -> ~tau_1
     . H_3)'
      endif
      if(brstau2A(1).ne.0.D0) then
      write(nout,102) brstau2A(1),2,istau1,ihA1,  'BR(~tau_2 -> ~tau_1
     . A_1)'
      endif
      if(brstau2A(2).ne.0.D0) then
      write(nout,102) brstau2A(2),2,istau1,ihA2,  'BR(~tau_2 -> ~tau_1
     . A_2)'
      endif
      if(brstau2ztau.ne.0.D0) then
      write(nout,102) brstau2ztau,2,istau1,iz,    'BR(~tau_2 -> ~tau_1
     . Z)'
      endif
      endif
C-----------
* THREE BODY 
C-----------
      if(stau2tot3.ne.0.D0) then
      write(nout,49) 'stau_2 3-body decays'
      write(nout,103)
      if( brstau2stau1star.ne.0.D0) then
      write(nout,104) brstau2stau1star,3,-istau1,itau,-itau,    
     .'BR(~tau_2 -> ~tau_1* tau taub)'
      endif
      if( brstau2stau1.ne.0.D0) then
      write(nout,104) brstau2stau1,3,istau1,itau,itau,    
     .'BR(~tau_2 -> ~tau_1 tau tau)'
      endif
      if(brstau2stau1nn.ne.0.D0) then
      write(nout,104) brstau2stau1nn,3,istau1,intau,intau,    
     .'BR(~tau_2 -> ~tau_1 nu_tau nu_tau)'
      endif
      endif

c -------------------- c
c Snu_electronL decays c
c -------------------- c
      write(nout,99)
      write(nout,100) 1000012,sneltot,'snu_eL'
      if(sneltot2.ne.0.D0) then
        write(nout,49) 'snu_eL 2-body decays'
         write(nout,101)
      if(brsnellneut(1).ne.0.D0) then
      write(nout,102) brsnellneut(1),2,in1,ine, 'BR(~nu_eL -> ~chi_10 nu
     ._e)'
      endif
      if(brsnellneut(2).ne.0.D0) then
      write(nout,102) brsnellneut(2),2,in2,ine, 'BR(~nu_eL -> ~chi_20 nu
     ._e)'
      endif
      if(brsnellneut(3).ne.0.D0) then
      write(nout,102) brsnellneut(3),2,in3,ine, 'BR(~nu_eL -> ~chi_30 nu
     ._e)'
      endif
      if(brsnellneut(4).ne.0.D0) then
      write(nout,102) brsnellneut(4),2,in4,ine, 'BR(~nu_eL -> ~chi_40 nu
     ._e)'
      endif
      if(brsnellneut(5).ne.0.D0) then
      write(nout,102) brsnellneut(5),2,in5,ine, 'BR(~nu_eL -> ~chi_50 nu
     ._e)'
      endif
      if(brsnellchar(1).ne.0.D0) then
      write(nout,102) brsnellchar(1),2,ic1,ie,  'BR(~nu_eL -> ~chi_1+ e-
     .)'
      endif
      if(brsnellchar(2).ne.0.D0) then
      write(nout,102) brsnellchar(2),2,ic2,ie,  'BR(~nu_eL -> ~chi_2+ e-
     .)'
      endif
      endif
C-----------
* THREE BODY 
C-----------
      if(sneltot3.ne.0.D0) then
       write(nout,49) 'snu_eL 3-body decays'
       write(nout,103)
      if(brsnestau1star.ne.0.D0) then
      write(nout,104) brsnestau1star,3,ine,-istau1,-itau,
     .'BR(~nu_eL -> nu_e ~tau_1* taub)' 
      endif
      if(brsnestau1.ne.0.D0) then
      write(nout,104) brsnestau1,3,ine,istau1,itau,
     .'BR(~nu_eL -> nu_e ~tau_1 tau)' 
      endif
      if(brsnestau1nutau.ne.0.D0) then
      write(nout,104) brsnestau1nutau,3,ie,istau1,intau,
     .'BR(~nu_eL -> e- ~tau_1 nu_tau)' 
      endif
      endif

c ---------------- c
c Snu_muonL decays c
c ---------------- c
      write(nout,99)
      write(nout,100) 1000014,sneltot,'snu_muL'
      if(sneltot2.ne.0.D0) then
        write(nout,49) 'snu_muL 2-body decays'
        write(nout,101)
      if(brsnellneut(1).ne.0.D0) then
      write(nout,102) brsnellneut(1),2,in1,inmu, 'BR(~nu_muL -> ~chi_10 
     .nu_mu)'
      endif
      if(brsnellneut(2).ne.0.D0) then
      write(nout,102) brsnellneut(2),2,in2,inmu, 'BR(~nu_muL -> ~chi_20 
     .nu_mu)'
      endif
      if(brsnellneut(3).ne.0.D0) then
      write(nout,102) brsnellneut(3),2,in3,inmu, 'BR(~nu_muL -> ~chi_30 
     .nu_mu)'
      endif
      if(brsnellneut(4).ne.0.D0) then
      write(nout,102) brsnellneut(4),2,in4,inmu, 'BR(~nu_muL -> ~chi_40 
     .nu_mu)'
      endif
      if(brsnellneut(5).ne.0.D0) then
      write(nout,102) brsnellneut(5),2,in5,inmu, 'BR(~nu_muL -> ~chi_50 
     .nu_mu)'
      endif
      if(brsnellchar(1).ne.0.D0) then
      write(nout,102) brsnellchar(1),2,ic1,imu,  'BR(~nu_muL -> ~chi_1+ 
     .mu-)'
      endif
      if(brsnellchar(2).ne.0.D0) then
      write(nout,102) brsnellchar(2),2,ic2,imu,  'BR(~nu_muL -> ~chi_2+ 
     .mu-)'
      endif
      endif
C-----------
* THREE BODY 
C-----------
      if(sneltot3.ne.0.D0) then
      write(nout,49) 'snu_muL 3-body decays'
      write(nout,103)
      if(brsnestau1star.ne.0.D0) then
      write(nout,104) brsnestau1star,3,inmu,-istau1,-itau,
     .'BR(~nu_muL -> nu_mu ~tau_1* taub)' 
      endif
      if(brsnestau1.ne.0.D0) then
      write(nout,104) brsnestau1,3,inmu,istau1,itau,
     .'BR(~nu_muL -> nu_mu ~tau_1 tau)' 
      endif
      if(brsnestau1nutau.ne.0.D0) then
      write(nout,104) brsnestau1nutau,3,imu,istau1,intau,
     .'BR(~nu_muL -> mu- ~tau_1 nu_tau)' 
      endif
      endif

c --------------- c
c Snu_tauL decays c
c --------------- c
      write(nout,99)
      write(nout,100) 1000016,sntautot,'snu_tauL'
      if(sntautot2.ne.0.D0) then
        write(nout,49) 'sbu_tauL 2-body decays'
        write(nout,101)
      if(brsntauneut(1).ne.0.D0) then
      write(nout,102) brsntauneut(1),2,in1,intau, 'BR(~nu_tauL -> ~chi_1
     .0 nu_tau)'
      endif
      if(brsntauneut(2).ne.0.D0) then
      write(nout,102) brsntauneut(2),2,in2,intau, 'BR(~nu_tauL -> ~chi_2
     .0 nu_tau)'
      endif
      if(brsntauneut(3).ne.0.D0) then
      write(nout,102) brsntauneut(3),2,in3,intau, 'BR(~nu_tauL -> ~chi_3
     .0 nu_tau)'
      endif
      if(brsntauneut(4).ne.0.D0) then
      write(nout,102) brsntauneut(4),2,in4,intau, 'BR(~nu_tauL -> ~chi_4
     .0 nu_tau)'
      endif
      if(brsntauneut(5).ne.0.D0) then
      write(nout,102) brsntauneut(5),2,in5,intau, 'BR(~nu_tauL -> ~chi_5
     .0 nu_tau)'
      endif
      if(brsntauchar(1).ne.0.D0) then
      write(nout,102) brsntauchar(1),2,ic1,itau,  'BR(~nu_tauL -> ~chi_1
     .+ tau-)'
      endif
      if(brsntauchar(2).ne.0.D0) then
      write(nout,102) brsntauchar(2),2,ic2,itau,  'BR(~nu_tauL -> ~chi_2
     .+ tau-)'
      endif
      if(brsntau1hcstau(1).ne.0.D0) then
      write(nout,102) brsntau1hcstau(1),2,-istau1,-ihc,'BR(~nu_tauL -> ~
     .tau_1+ H-)'
      endif
      if(brsntau1hcstau(2).ne.0.D0) then
      write(nout,102) brsntau1hcstau(2),2,-istau2,-ihc,'BR(~nu_tauL -> ~
     .tau_2+ H-)'
      endif
      if(brsntau1wstau(1).ne.0.D0) then
      write(nout,102) brsntau1wstau(1),2,-istau1,-iwc, 'BR(~nu_tauL -> ~
     .tau_1+ W-)'
      endif
      if(brsntau1wstau(2).ne.0.D0) then
      write(nout,102) brsntau1wstau(2),2,-istau2,-iwc, 'BR(~nu_tauL -> ~
     .tau_2+ W-)'
      endif
      endif

C-----------
* THREE BODY 
C-----------
      if(sntautot3.ne.0.D0) then
        write(nout,49) 'snu_tauL 3-body decays'
        write(nout,103)
      if(brsntaustau1star.ne.0.D0) then
      write(nout,104) brsntaustau1star,3,intau,-istau1,-itau, 
     .'BR(~nu_tauL -> nu_tau ~tau_1* taub)'
      endif
      if(brsntaustau1.ne.0.D0) then
      write(nout,104) brsntaustau1,3,intau,istau1,itau, 
     .'BR(~nu_tauL -> nu_tau ~tau_1 tau)'
      endif
      if(brsntaustau1nutau.ne.0.D0) then
      write(nout,104) brsntaustau1nutau,3,itau,istau1,intau, 
     .'BR(~nu_tauL -> tau ~tau_1 nu_tau)'
      endif
      endif
c ------------ c
c Sup_L decays c
c ------------ c
      write(nout,99)
      write(nout,100) 1000002,supltot2,'sup_L'
      if(supltot2.ne.0.D0) then
        write(nout,49) 'sup_L 2-body decays'
        write(nout,101)
      if(brsuplnup(1).ne.0.D0) then
      write(nout,102) brsuplnup(1),2,in1,iu   ,'BR(~u_L -> ~chi_10 u)'
      endif
      if(brsuplnup(2).ne.0.D0) then
      write(nout,102) brsuplnup(2),2,in2,iu   ,'BR(~u_L -> ~chi_20 u)'
      endif
      if(brsuplnup(3).ne.0.D0) then
      write(nout,102) brsuplnup(3),2,in3,iu   ,'BR(~u_L -> ~chi_30 u)'
      endif
      if(brsuplnup(4).ne.0.D0) then
      write(nout,102) brsuplnup(4),2,in4,iu   ,'BR(~u_L -> ~chi_40 u)'
      endif
      if(brsuplnup(5).ne.0.D0) then
      write(nout,102) brsuplnup(5),2,in5,iu   ,'BR(~u_L -> ~chi_50 u)'
      endif
      if(brsuplcdow(1).ne.0.D0) then
      write(nout,102) brsuplcdow(1),2,ic1,id  ,'BR(~u_L -> ~chi_1+ d)'
      endif
      if(brsuplcdow(2).ne.0.D0) then
      write(nout,102) brsuplcdow(2),2,ic2,id  ,'BR(~u_L -> ~chi_2+ d)'
      endif
      if(brsuplglui.ne.0.D0) then
      write(nout,102) brsuplglui,2,iglo,iu    ,'BR(~u_L -> ~g      u)'
      endif
      endif

c ------------ c
c Sup_R decays c
c ------------ c
      write(nout,99)
      write(nout,100) 2000002,suprtot2,'sup_R'
      if(suprtot2.ne.0.D0) then
        write(nout,49) 'sup_R 2-body decays'
       write(nout,101)
      if(brsuprnup(1).ne.0.D0) then
      write(nout,102) brsuprnup(1),2,in1,iu   ,'BR(~u_R -> ~chi_10 u)'
      endif
      if(brsuprnup(2).ne.0.D0) then
      write(nout,102) brsuprnup(2),2,in2,iu   ,'BR(~u_R -> ~chi_20 u)'
      endif
      if(brsuprnup(3).ne.0.D0) then
      write(nout,102) brsuprnup(3),2,in3,iu   ,'BR(~u_R -> ~chi_30 u)'
      endif
      if(brsuprnup(4).ne.0.D0) then
      write(nout,102) brsuprnup(4),2,in4,iu   ,'BR(~u_R -> ~chi_40 u)'
      endif
      if(brsuprnup(5).ne.0.D0) then
      write(nout,102) brsuprnup(5),2,in5,iu   ,'BR(~u_R -> ~chi_50 u)'
      endif
      if(brsuprcdow(1).ne.0.D0) then
      write(nout,102) brsuprcdow(1),2,ic1,id  ,'BR(~u_R -> ~chi_1+ d)'
      endif
      if(brsuprcdow(2).ne.0.D0) then
      write(nout,102) brsuprcdow(2),2,ic2,id  ,'BR(~u_R -> ~chi_2+ d)'
      endif
      if(brsuprglui.ne.0.D0) then
      write(nout,102) brsuprglui,2,iglo,iu    ,'BR(~u_R -> ~g      u)'
      endif
      endif

c -------------- c
c Sdown_L decays c
c -------------- c
      write(nout,99)
      write(nout,100) 1000001,sdowltot2,'sdown_L'
      if(sdowltot2.ne.0.D0) then
        write(nout,49) 'sdown_L 2-body decays'
        write(nout,101)
      if(brsdowlndow(1).ne.0.D0) then
      write(nout,102) brsdowlndow(1),2,in1,id  ,'BR(~d_L -> ~chi_10 d)'
      endif
      if(brsdowlndow(2).ne.0.D0) then
      write(nout,102) brsdowlndow(2),2,in2,id  ,'BR(~d_L -> ~chi_20 d)'
      endif
      if(brsdowlndow(3).ne.0.D0) then
      write(nout,102) brsdowlndow(3),2,in3,id  ,'BR(~d_L -> ~chi_30 d)'
      endif
      if(brsdowlndow(4).ne.0.D0) then
      write(nout,102) brsdowlndow(4),2,in4,id  ,'BR(~d_L -> ~chi_40 d)'
      endif
      if(brsdowlndow(5).ne.0.D0) then
      write(nout,102) brsdowlndow(5),2,in5,id  ,'BR(~d_L -> ~chi_50 d)'
      endif
      if(brsdowlchup(1).ne.0.D0) then
      write(nout,102) brsdowlchup(1),2,-ic1,iu ,'BR(~d_L -> ~chi_1- u)'
      endif
      if(brsdowlchup(2).ne.0.D0) then
      write(nout,102) brsdowlchup(2),2,-ic2,iu ,'BR(~d_L -> ~chi_2- u)'
      endif
      if(brsdowlglui.ne.0.D0) then
      write(nout,102) brsdowlglui,2,iglo,id    ,'BR(~d_L -> ~g      d)'
      endif
      endif

c -------------- c
c Sdown_R decays c
c -------------- c
      write(nout,99)
      write(nout,100) 2000001,sdowrtot2,'sdown_R'
      if(sdowrtot2.ne.0.D0) then
        write(nout,49) 'sdown_R 2-body decays'
        write(nout,101)
      if(brsdowrndow(1).ne.0.D0) then
      write(nout,102) brsdowrndow(1),2,in1,id  ,'BR(~d_R -> ~chi_10 d)'
      endif
      if(brsdowrndow(2).ne.0.D0) then
      write(nout,102) brsdowrndow(2),2,in2,id  ,'BR(~d_R -> ~chi_20 d)'
      endif
      if(brsdowrndow(3).ne.0.D0) then
      write(nout,102) brsdowrndow(3),2,in3,id  ,'BR(~d_R -> ~chi_30 d)'
      endif
      if(brsdowrndow(4).ne.0.D0) then
      write(nout,102) brsdowrndow(4),2,in4,id  ,'BR(~d_R -> ~chi_40 d)'
      endif
      if(brsdowrndow(5).ne.0.D0) then
      write(nout,102) brsdowrndow(5),2,in5,id  ,'BR(~d_R -> ~chi_50 d)'
      endif
      if(brsdowrchup(1).ne.0.D0) then
      write(nout,102) brsdowrchup(1),2,-ic1,iu ,'BR(~d_R -> ~chi_1- u)'
      endif
      if(brsdowrchup(2).ne.0.D0) then
      write(nout,102) brsdowrchup(2),2,-ic2,iu ,'BR(~d_R -> ~chi_2- u)'
      endif
      if(brsdowrglui.ne.0.D0) then
      write(nout,102) brsdowrglui,2,iglo,id    ,'BR(~d_R -> ~g      d)'
      endif
      endif

c --------------- c
c Scharm_L decays c
c --------------- c
      write(nout,99)
      write(nout,100) 1000004,supltot2,'scharm_L'
      if(supltot2.ne.0.D0) then
        write(nout,49) 'scharm_L 2-body decays'
        write(nout,101)
      if(brsuplnup(1).ne.0.D0) then
      write(nout,102) brsuplnup(1),2,in1,ic   ,'BR(~c_L -> ~chi_10 c)'
      endif
      if(brsuplnup(2).ne.0.D0) then
      write(nout,102) brsuplnup(2),2,in2,ic   ,'BR(~c_L -> ~chi_20 c)'
      endif
      if(brsuplnup(3).ne.0.D0) then
      write(nout,102) brsuplnup(3),2,in3,ic   ,'BR(~c_L -> ~chi_30 c)'
      endif
      if(brsuplnup(4).ne.0.D0) then
      write(nout,102) brsuplnup(4),2,in4,ic   ,'BR(~c_L -> ~chi_40 c)'
      endif
      if(brsuplnup(5).ne.0.D0) then
      write(nout,102) brsuplnup(5),2,in5,ic   ,'BR(~c_L -> ~chi_50 c)'
      endif
      if(brsuplcdow(1).ne.0.D0) then
      write(nout,102) brsuplcdow(1),2,ic1,is  ,'BR(~c_L -> ~chi_1+ s)'
      endif
      if(brsuplcdow(2).ne.0.D0) then
      write(nout,102) brsuplcdow(2),2,ic2,is  ,'BR(~c_L -> ~chi_2+ s)'
      endif
      if(brsuplglui.ne.0.D0) then
      write(nout,102) brsuplglui,2,iglo,ic    ,'BR(~c_L -> ~g      c)'
      endif
      endif

c --------------- c
c Scharm_R decays c
c --------------- c
      write(nout,99)
      write(nout,100) 2000004,suprtot2,'scharm_R'
      if(suprtot2.ne.0.D0) then
        write(nout,49) 'scharm_R 2-body decays'
        write(nout,101)
      if(brsuprnup(1).ne.0.D0) then
      write(nout,102) brsuprnup(1),2,in1,ic   ,'BR(~c_R -> ~chi_10 c)'
      endif
      if(brsuprnup(2).ne.0.D0) then
      write(nout,102) brsuprnup(2),2,in2,ic   ,'BR(~c_R -> ~chi_20 c)'
      endif
      if(brsuprnup(3).ne.0.D0) then      
      write(nout,102) brsuprnup(3),2,in3,ic   ,'BR(~c_R -> ~chi_30 c)'
      endif
      if(brsuprnup(4).ne.0.D0) then
      write(nout,102) brsuprnup(4),2,in4,ic   ,'BR(~c_R -> ~chi_40 c)'
      endif
      if(brsuprnup(5).ne.0.D0) then
      write(nout,102) brsuprnup(5),2,in5,ic   ,'BR(~c_R -> ~chi_50 c)'
      endif
      if(brsuprcdow(1).ne.0.D0) then
      write(nout,102) brsuprcdow(1),2,ic1,is  ,'BR(~c_R -> ~chi_1+ s)'
      endif
      if(brsuprcdow(2).ne.0.D0) then
      write(nout,102) brsuprcdow(2),2,ic2,is  ,'BR(~c_R -> ~chi_2+ s)'
      endif
      if(brsuprglui.ne.0.D0) then
      write(nout,102) brsuprglui,2,iglo,ic    ,'BR(~c_R -> ~g      c)'
      endif
      endif

c ----------------- c
c Sstrange_L decays c
c ----------------- c
      write(nout,99)
      write(nout,100) 1000003,sdowltot2,'sstrange_L'
      if(sdowltot2.ne.0.D0) then
        write(nout,49) 'sstrange_L 2-body decays'
        write(nout,101)
      if(brsdowlndow(1).ne.0.D0) then
      write(nout,102) brsdowlndow(1),2,in1,is  ,'BR(~s_L -> ~chi_10 s)'
      endif
      if(brsdowlndow(2).ne.0.D0) then
      write(nout,102) brsdowlndow(2),2,in2,is  ,'BR(~s_L -> ~chi_20 s)'
      endif
      if(brsdowlndow(3).ne.0.D0) then
      write(nout,102) brsdowlndow(3),2,in3,is  ,'BR(~s_L -> ~chi_30 s)'
      endif
      if(brsdowlndow(4).ne.0.D0) then
      write(nout,102) brsdowlndow(4),2,in4,is  ,'BR(~s_L -> ~chi_40 s)'
      endif
      if(brsdowlndow(5).ne.0.D0) then
      write(nout,102) brsdowlndow(5),2,in5,is  ,'BR(~s_L -> ~chi_50 s)'
      endif
      if(brsdowlchup(1).ne.0.D0) then
      write(nout,102) brsdowlchup(1),2,-ic1,ic ,'BR(~s_L -> ~chi_1- c)'
      endif
      if(brsdowlchup(2).ne.0.D0) then
      write(nout,102) brsdowlchup(2),2,-ic2,ic ,'BR(~s_L -> ~chi_2- c)'
      endif
      if(brsdowlglui.ne.0.D0) then
      write(nout,102) brsdowlglui,2,iglo,is    ,'BR(~s_L -> ~g      s)'
      endif
      endif

c ----------------- c
c Sstrange_R decays c
c ----------------- c
      write(nout,99)
      write(nout,100) 2000003,sdowrtot2,'sstrange_R'
      if(sdowrtot2.ne.0.D0) then
        write(nout,49) 'sstrange_R 2-body decays'
        write(nout,101)
      if(brsdowrndow(1).ne.0.D0) then
      write(nout,102) brsdowrndow(1),2,in1,is  ,'BR(~s_R -> ~chi_10 s)'
      endif
      if(brsdowrndow(2).ne.0.D0) then
      write(nout,102) brsdowrndow(2),2,in2,is  ,'BR(~s_R -> ~chi_20 s)'
      endif
      if(brsdowrndow(3).ne.0.D0) then
      write(nout,102) brsdowrndow(3),2,in3,is  ,'BR(~s_R -> ~chi_30 s)'
      endif
      if(brsdowrndow(4).ne.0.D0) then
      write(nout,102) brsdowrndow(4),2,in4,is  ,'BR(~s_R -> ~chi_40 s)'
      endif
      if(brsdowrndow(5).ne.0.D0) then
      write(nout,102) brsdowrndow(5),2,in5,is  ,'BR(~s_R -> ~chi_50 s)'
      endif
      if(brsdowrchup(1).ne.0.D0) then
      write(nout,102) brsdowrchup(1),2,-ic1,ic ,'BR(~s_R -> ~chi_1- c)'
      endif
      if(brsdowrchup(2).ne.0.D0) then
      write(nout,102) brsdowrchup(2),2,-ic2,ic ,'BR(~s_R -> ~chi_2- c)'
      endif
      if(brsdowrglui.ne.0.D0) then
      write(nout,102) brsdowrglui,2,iglo,is    ,'BR(~s_R -> ~g      s)'
      endif
      endif

c ------------ c
c Stop1 decays c
c ------------ c
      write(nout,99)
      write(nout,100) 1000006,stoptot(1),'stop1'
      if(stoptot2(1).ne.0.D0) then
        write(nout,49) 'stop1 2-body decays'
        write(nout,101)
      if(brst1neutt(1).ne.0.D0) then
      write(nout,102) brst1neutt(1),2,in1,it  ,'BR(~t_1 -> ~chi_10 t )'
      endif
      if(brst1neutt(2).ne.0.D0) then
      write(nout,102) brst1neutt(2),2,in2,it  ,'BR(~t_1 -> ~chi_20 t )'
      endif
      if(brst1neutt(3).ne.0.D0) then
      write(nout,102) brst1neutt(3),2,in3,it  ,'BR(~t_1 -> ~chi_30 t )'
      endif
      if(brst1neutt(4).ne.0.D0) then
      write(nout,102) brst1neutt(4),2,in4,it  ,'BR(~t_1 -> ~chi_40 t )'
      endif
      if(brst1neutt(5).ne.0.D0) then
      write(nout,102) brst1neutt(5),2,in5,it  ,'BR(~t_1 -> ~chi_50 t )'
      endif
      if(brst1charb(1).ne.0.D0) then
      write(nout,102) brst1charb(1),2,ic1,ib  ,'BR(~t_1 -> ~chi_1+ b )'
      endif
      if(brst1charb(2).ne.0.D0) then
      write(nout,102) brst1charb(2),2,ic2,ib  ,'BR(~t_1 -> ~chi_2+ b )'
      endif
      if(brst1glui.ne.0.D0) then
      write(nout,102) brst1glui    ,2,iglo,it ,'BR(~t_1 -> ~g      t )'
      endif
      if(brst1hcsb(1).ne.0.D0) then
      write(nout,102) brst1hcsb(1) ,2,isb1,ihc,'BR(~t_1 -> ~b_1    H+)'
      endif
      if(brst1hcsb(2).ne.0.D0) then
      write(nout,102) brst1hcsb(2) ,2,isb2,ihc,'BR(~t_1 -> ~b_2    H+)'
      endif
      if(brst1wsb(1).ne.0.D0) then
      write(nout,102) brst1wsb(1)  ,2,isb1,iwc,'BR(~t_1 -> ~b_1    W+)'
      endif
      if(brst1wsb(2).ne.0.D0) then
      write(nout,102) brst1wsb(2)  ,2,isb2,iwc,'BR(~t_1 -> ~b_2    W+)'
      endif
      endif
C     =============================
*     STOP1: Radiative decay
C     =============================
            if(stoptotrad(1).ne.0.D0)Then
            if(flagloop.eq.1.D0) then
      if(brgamma.ne.0.D0) then
      write(nout,102) brgamma      ,2,in1,ic  ,'BR(~t_1 -> ~chi_10 c )'
      endif
      if(brgammaup.ne.0.D0) then
      write(nout,102) brgammaup    ,2,in1,iu  ,'BR(~t_1 -> ~chi_10 u )'
      endif
      if(brgammagluino.ne.0.D0) then
      write(nout,102) brgammagluino,2,iglo,ic ,'BR(~t_1 -> ~g      c )'
      endif
            endif
            endif
C   ========================
*   STOP1: Three body
C   ========================
      if(stoptot3(1).ne.0.D0) then
         write(nout,49) 'stop1 3-body decays'   
         write(nout,103)
      if(brstopw(1,1).ne.0.D0) then
      write(nout,104) brstopw(1,1),3,in1,ib,iwc,       'BR(~t_1 -> ~chi_
     .10  b  W+)'
      endif
      if(brstopw(1,2).ne.0.D0) then
      write(nout,104) brstopw(1,2),3,in2,ib,iwc,       'BR(~t_1 -> ~chi_
     .20  b  W+)'
      endif
      if(brstopw(1,3).ne.0.D0) then
      write(nout,104) brstopw(1,3),3,in3,ib,iwc,       'BR(~t_1 -> ~chi_
     .30  b  W+)'
      endif
      if(brstopw(1,4).ne.0.D0) then
      write(nout,104) brstopw(1,4),3,in4,ib,iwc,       'BR(~t_1 -> ~chi_
     .40  b  W+)'
      endif
      if(brstopw(1,5).ne.0.D0) then
      write(nout,104) brstopw(1,5),3,in5,ib,iwc,       'BR(~t_1 -> ~chi_
     .50  b  W+)'
      endif
      if(brstoph(1,1).ne.0.D0) then
      write(nout,104) brstoph(1,1),3,in1,ib,ihc,       'BR(~t_1 -> ~chi_
     .10  b  H+)'
      endif
      if(brstoph(1,2).ne.0.D0) then
      write(nout,104) brstoph(1,2),3,in2,ib,ihc,       'BR(~t_1 -> ~chi_
     .20  b  H+)'
      endif
      if(brstoph(1,3).ne.0.D0) then
      write(nout,104) brstoph(1,3),3,in3,ib,ihc,       'BR(~t_1 -> ~chi_
     .30  b  H+)'
      endif
      if(brstoph(1,4).ne.0.D0) then
      write(nout,104) brstoph(1,4),3,in4,ib,ihc,       'BR(~t_1 -> ~chi_
     .40  b  H+)'
      endif
      if(brstoph(1,5).ne.0.D0) then
      write(nout,104) brstoph(1,5),3,in5,ib,ihc,       'BR(~t_1 -> ~chi_
     .50  b  H+)'
      endif
      if(brstsntau(1,1).ne.0.D0) then
      write(nout,104) brstsntau(1,1),3,intau1,ib,-itau,'BR(~t_1 -> ~nu_t
     .auL b  tau+)'
      endif
      if(brstsnel(1).ne.0.D0) then
      write(nout,104) brstsnel(1),3,inel,ib,-ie,       'BR(~t_1 -> ~nu_e
     .L   b  e+)'
      write(nout,104) brstsnel(1),3,inmul,ib,-imu,     'BR(~t_1 -> ~nu_m
     .uL  b  mu+)'
      endif
      if(brststau(1,1).ne.0.D0) then
      write(nout,104) brststau(1,1),3,-istau1,ib,intau,'BR(~t_1 -> ~tau_
     .1+  b  nu_tau)'
      endif
      if(brststau(1,2).ne.0.D0) then
      write(nout,104) brststau(1,2),3,-istau2,ib,intau,'BR(~t_1 -> ~tau_
     .2+  b  nu_tau)'
      endif
      if(brstsel(1,1).ne.0.D0) then
      write(nout,104) brstsel(1,1),3,-isell,ib,ine,    'BR(~t_1 -> ~e_L+
     .    b  nu_e)'
      endif
      if(brstsel(1,2).ne.0.D0) then
      write(nout,104) brstsel(1,2),3,-iselr,ib,ine,    'BR(~t_1 -> ~e_R+
     .    b  nu_e)'
      endif
      if(brstsel(1,1).ne.0.D0) then
      write(nout,104) brstsel(1,1),3,-ismul,ib,inmu,   'BR(~t_1 -> ~mu_L
     .+   b  nu_mu)'
      endif
      if(brstsel(1,2).ne.0.D0) then
      write(nout,104) brstsel(1,2),3,-ismur,ib,inmu,   'BR(~t_1 -> ~mu_R
     .+   b  nu_mu)'
      endif
      if(brstbsbst(1,1).ne.0.D0) then
      write(nout,104) brstbsbst(1,1),3,-isb1,ib,it,    'BR(~t_1 -> ~b_1*
     .    b  t)'
      endif
      if(brstbsbst(1,2).ne.0.D0) then
      write(nout,104) brstbsbst(1,2),3,-isb2,ib,it,    'BR(~t_1 -> ~b_2*
     .    b  t)'
      endif
      if(brstbbsbt(1,1).ne.0.D0) then
      write(nout,104) brstbbsbt(1,1),3,isb1,ibb,it,    'BR(~t_1 -> ~b_1
     .    bb t)'
      endif
      if(brstbbsbt(1,2).ne.0.D0) then
      write(nout,104) brstbbsbt(1,2),3,isb2,ibb,it,    'BR(~t_1 -> ~b_2
     .    bb t)'
      endif
      if(brstupsbdow(1,1).ne.0.D0) then
      write(nout,104) brstupsbdow(1,1),3,isb1,idb,iu,  'BR(~t_1 -> ~b_1
     .    db u)'
      endif
      if(brstupsbdow(1,2).ne.0.D0) then
      write(nout,104) brstupsbdow(1,2),3,isb2,idb,iu,  'BR(~t_1 -> ~b_2
     .    db u)'
      endif
      if(brstupsbdow(1,1).ne.0.D0) then
      write(nout,104) brstupsbdow(1,1),3,isb1,isb,ic,  'BR(~t_1 -> ~b_1
     .    sb c)'
      endif
      if(brstupsbdow(1,2).ne.0.D0) then
      write(nout,104) brstupsbdow(1,2),3,isb2,isb,ic,  'BR(~t_1 -> ~b_2
     .    sb c)'
      endif
      if(brsttausbnu(1,1).ne.0.D0) then
      write(nout,104) brsttausbnu(1,1),3,isb1,-itau,intau,'BR(~t_1 -> ~b
     ._1     tau+ nu_tau)'
      endif
      if(brsttausbnu(1,2).ne.0.D0) then
      write(nout,104) brsttausbnu(1,2),3,isb2,-itau,intau,'BR(~t_1 -> ~b
     ._2     tau+ nu_tau)'
      endif
      if(brstelsbnu(1,1).ne.0.D0) then
      write(nout,104) brstelsbnu(1,1),3,isb1,-ie,ine,  'BR(~t_1 -> ~b_1 
     .    e+   nu_e)'
      endif
      if(brstelsbnu(1,2).ne.0.D0) then
      write(nout,104) brstelsbnu(1,2),3,isb2,-ie,ine,  'BR(~t_1 -> ~b_2 
     .    e+   nu_e)'
      endif
      if(brstelsbnu(1,1).ne.0.D0) then
      write(nout,104) brstelsbnu(1,1),3,isb1,-imu,inmu,'BR(~t_1 -> ~b_1 
     .    mu+  nu_mu)'
      endif
      if(brstelsbnu(1,2).ne.0.D0) then
      write(nout,104) brstelsbnu(1,2),3,isb2,-imu,inmu,'BR(~t_1 -> ~b_2 
     .    mu+  nu_mu)'
      endif
      endif

c ------------ c
c Stop2 decays c
c ------------ c
      write(nout,99)
      write(nout,100) 2000006,stoptot(2),'stop2'
      if(stoptot2(2).ne.0.D0) then
        write(nout,49) 'stop2 2-body decays'
        write(nout,101)
      if(brst2neutt(1).ne.0.D0) then
      write(nout,102) brst2neutt(1),2,in1,it  ,'BR(~t_2 -> ~chi_10 t )'
      endif
      if(brst2neutt(2).ne.0.D0) then
      write(nout,102) brst2neutt(2),2,in2,it  ,'BR(~t_2 -> ~chi_20 t )'
      endif
      if(brst2neutt(2).ne.0.D0) then
      write(nout,102) brst2neutt(3),2,in3,it  ,'BR(~t_2 -> ~chi_30 t )'
      endif
      if(brst2neutt(4).ne.0.D0) then
      write(nout,102) brst2neutt(4),2,in4,it  ,'BR(~t_2 -> ~chi_40 t )'
      endif
      if(brst2neutt(5).ne.0.D0) then
      write(nout,102) brst2neutt(5),2,in5,it  ,'BR(~t_2 -> ~chi_50 t )'
      endif
      if(brst2charb(1).ne.0.D0) then
      write(nout,102) brst2charb(1),2,ic1,ib  ,'BR(~t_2 -> ~chi_1+ b )'
      endif
      if(brst2charb(2).ne.0.D0) then
      write(nout,102) brst2charb(2),2,ic2,ib  ,'BR(~t_2 -> ~chi_2+ b )'
      endif
      if(brst2glui.ne.0.D0) then
      write(nout,102) brst2glui    ,2,iglo,it ,'BR(~t_2 -> ~g      t )'
      endif
      if(brst2H(1).ne.0.D0) then
      write(nout,102) brst2H(1)  ,2,ist1,ihH1, 'BR(~t_2 -> ~t_1    H1 )'
      endif
      if(brst2H(2).ne.0.D0) then
      write(nout,102) brst2H(2)  ,2,ist1,ihH2,'BR(~t_2 -> ~t_1    H2 )'
      endif
      if(brst2H(3).ne.0.D0) then
      write(nout,102) brst2H(3)  ,2,ist1,ihH3,'BR(~t_2 -> ~t_1    H3 )'
      endif
      if(brst2A(1).ne.0.D0) then
      write(nout,102) brst2A(1)  ,2,ist1,ihA1,'BR(~t_2 -> ~t_1    A1 )'
      endif
      if(brst2A(2).ne.0.D0) then
      write(nout,102) brst2A(2)  ,2,ist1,ihA2,'BR(~t_2 -> ~t_1    A2 )'
      endif
      if(brst2hcsb(1).ne.0.D0) then
      write(nout,102) brst2hcsb(1) ,2,isb1,ihc,'BR(~t_2 -> ~b_1    H+)'
      endif
      if(brst2hcsb(2).ne.0.D0) then
      write(nout,102) brst2hcsb(2) ,2,isb2,ihc,'BR(~t_2 -> ~b_2    H+)'
      endif
      if(brst2ztop.ne.0.D0) then
      write(nout,102) brst2ztop    ,2,ist1,iz ,'BR(~t_2 -> ~t_1    Z )'
      endif
      if(brst2wsb(1).ne.0.D0) then
      write(nout,102) brst2wsb(1)  ,2,isb1,iwc,'BR(~t_2 -> ~b_1    W+)'
      endif
      if(brst2wsb(2).ne.0.D0) then
      write(nout,102) brst2wsb(2)  ,2,isb2,iwc,'BR(~t_2 -> ~b_2    W+)'
      endif
      endif
c
C   ========================
*   -STOP2: Three body
C   ========================
      if(stoptot3(2).ne.0.D0)then
      write(nout,49) 'stop2 3-body decays'
      write(nout,103)
      if(brstopw(2,1).ne.0.D0) then
      write(nout,104) brstopw(2,1),3,in1,ib,iwc,       'BR(~t_2 -> ~chi_
     .10  b  W+)'
      endif
      if(brstopw(2,2).ne.0.D0) then
      write(nout,104) brstopw(2,2),3,in2,ib,iwc,       'BR(~t_2 -> ~chi_
     .20  b  W+)'
      endif
      if(brstopw(2,3).ne.0.D0) then
      write(nout,104) brstopw(2,3),3,in3,ib,iwc,       'BR(~t_2 -> ~chi_
     .30  b  W+)'
      endif
      if(brstopw(2,4).ne.0.D0) then
      write(nout,104) brstopw(2,4),3,in4,ib,iwc,       'BR(~t_2 -> ~chi_
     .40  b  W+)'
      endif
      if(brstopw(2,5).ne.0.D0) then
      write(nout,104) brstopw(2,5),3,in5,ib,iwc,       'BR(~t_2 -> ~chi_
     .50  b  W+)'
      endif
      if(brstoph(2,1).ne.0.D0) then
      write(nout,104) brstoph(2,1),3,in1,ib,ihc,       'BR(~t_2 -> ~chi_
     .10  b  H+)'
      endif
      if(brstoph(2,2).ne.0.D0) then
      write(nout,104) brstoph(2,2),3,in2,ib,ihc,       'BR(~t_2 -> ~chi_
     .20  b  H+)'
      endif
      if(brstoph(2,3).ne.0.D0) then
      write(nout,104) brstoph(2,3),3,in3,ib,ihc,       'BR(~t_2 -> ~chi_
     .30  b  H+)'
      endif
      if(brstoph(2,4).ne.0.D0) then
      write(nout,104) brstoph(2,4),3,in4,ib,ihc,       'BR(~t_2 -> ~chi_
     .40  b  H+)'
      endif
      if(brstoph(2,5).ne.0.D0) then
      write(nout,104) brstoph(2,5),3,in5,ib,ihc,       'BR(~t_2 -> ~chi_
     .50  b  H+)'
      endif
      if(brstsntau(2,1).ne.0.D0) then
      write(nout,104) brstsntau(2,1),3,intau1,ib,-itau,'BR(~t_2 -> ~nu_t
     .auL b  tau+)'
      endif
      if(brstsnel(2).ne.0.D0) then
      write(nout,104) brstsnel(2),3,inel,ib,-ie,       'BR(~t_2 -> ~nu_e
     .L   b  e+)'
      write(nout,104) brstsnel(2),3,inmul,ib,-imu,     'BR(~t_2 -> ~nu_m
     .uL  b  mu+)'
      endif
      if(brststau(2,1).ne.0.D0) then
      write(nout,104) brststau(2,1),3,-istau1,ib,intau,'BR(~t_2 -> ~tau_
     .1+  b  nu_tau)'
      endif
      if(brststau(2,2).ne.0.D0) then
      write(nout,104) brststau(2,2),3,-istau2,ib,intau,'BR(~t_2 -> ~tau_
     .2+  b  nu_tau)'
      endif
      if(brstsel(2,1).ne.0.D0) then
      write(nout,104) brstsel(2,1),3,-isell,ib,ine,    'BR(~t_2 -> ~e_L+
     .    b  nu_e)'
      endif
      if(brstsel(2,2).ne.0.D0) then
      write(nout,104) brstsel(2,2),3,-iselr,ib,ine,    'BR(~t_2 -> ~e_R+
     .    b  nu_e)'
      endif
      if(brstsel(2,1).ne.0.D0) then      
      write(nout,104) brstsel(2,1),3,-ismul,ib,inmu,   'BR(~t_2 -> ~mu_L
     .+   b  nu_mu)'
      endif
      if(brstsel(2,2).ne.0.D0) then
      write(nout,104) brstsel(2,2),3,-ismur,ib,inmu,   'BR(~t_2 -> ~mu_R
     .+   b  nu_mu)'
      endif
      if(brstbsbst(2,1).ne.0.D0) then
      write(nout,104) brstbsbst(2,1),3,-isb1,ib,it,    'BR(~t_2 -> ~b_1*
     .    b  t)'
      endif
      if(brstbsbst(2,2).ne.0.D0) then
      write(nout,104) brstbsbst(2,2),3,-isb2,ib,it,    'BR(~t_2 -> ~b_2*
     .    b  t)'
      endif
      if(brstbbsbt(2,1).ne.0.D0) then
      write(nout,104) brstbbsbt(2,1),3,isb1,ibb,it,    'BR(~t_2 -> ~b_1
     .    bb t)'
      endif
      if(brstbbsbt(2,2).ne.0.D0) then
      write(nout,104) brstbbsbt(2,2),3,isb2,ibb,it,    'BR(~t_2 -> ~b_2
     .    bb t)'
      endif
      if(brstupsbdow(2,1).ne.0.D0) then
      write(nout,104) brstupsbdow(2,1),3,isb1,idb,iu,  'BR(~t_2 -> ~b_1
     .    db u)'
      endif
      if(brstupsbdow(2,2).ne.0.D0) then
      write(nout,104) brstupsbdow(2,2),3,isb2,idb,iu,  'BR(~t_2 -> ~b_2
     .    db u)'
      endif
      if(brstupsbdow(2,1).ne.0.D0) then
      write(nout,104) brstupsbdow(2,1),3,isb1,isb,ic,  'BR(~t_2 -> ~b_1
     .    sb c)'
      endif
      if(brstupsbdow(2,2).ne.0.D0) then
      write(nout,104) brstupsbdow(2,2),3,isb2,isb,ic,  'BR(~t_2 -> ~b_2
     .    sb c)'
      endif
      if(brsttausbnu(2,1).ne.0.D0) then
      write(nout,104) brsttausbnu(2,1),3,isb1,-itau,intau,'BR(~t_2 -> ~b
     ._1     tau+ nu_tau)'
      endif
      if(brsttausbnu(2,2).ne.0.D0) then
      write(nout,104) brsttausbnu(2,2),3,isb2,-itau,intau,'BR(~t_2 -> ~b
     ._2     tau+ nu_tau)'
      endif
      if(brstelsbnu(2,1).ne.0.D0) then
      write(nout,104) brstelsbnu(2,1),3,isb1,-ie,ine,  'BR(~t_2 -> ~b_1 
     .    e+   nu_e)'
      endif
      if(brstelsbnu(2,2).ne.0.D0) then
      write(nout,104) brstelsbnu(2,2),3,isb2,-ie,ine,  'BR(~t_2 -> ~b_2 
     .    e+   nu_e)'
      endif
      if(brstelsbnu(2,1).ne.0.D0) then
      write(nout,104) brstelsbnu(2,1),3,isb1,-imu,inmu,'BR(~t_2 -> ~b_1 
     .    mu+  nu_mu)'
      endif
      if(brstelsbnu(2,2).ne.0.D0) then
      write(nout,104) brstelsbnu(2,2),3,isb2,-imu,inmu,'BR(~t_2 -> ~b_2 
     .    mu+  nu_mu)'
      endif
      if(brst2st1tt.ne.0.D0) then
      write(nout,104) brst2st1tt,3,ist1,it,itb,        'BR(~t_2 -> ~t_1
     .    t    tb)'
      write(nout,104) brst2st1tt,3,-ist1,it,it,        'BR(~t_2 -> ~t_1*
     .    t    t )'
      endif
      if(brst2st1bb.ne.0.D0) then
      write(nout,104) brst2st1bb,3,ist1,ib,ibb,        'BR(~t_2 -> ~t_1
     .    b    bb)'
      endif
      if(brst2st1uu.ne.0.D0) then
      write(nout,104) brst2st1uu,3,ist1,iu,iub,        'BR(~t_2 -> ~t_1
     .    u    ub)'
      endif
      if(brst2st1dd.ne.0.D0) then
      write(nout,104) brst2st1dd,3,ist1,id,idb,        'BR(~t_2 -> ~t_1
     .    d    db)'
      endif
      if(brst2st1uu.ne.0.D0) then
      write(nout,104) brst2st1uu,3,ist1,ic,icb,        'BR(~t_2 -> ~t_1
     .    c    cb)'
      endif
      if(brst2st1dd.ne.0.D0) then
      write(nout,104) brst2st1dd,3,ist1,is,isb,        'BR(~t_2 -> ~t_1
     .    s    sb)'
      endif
      if(brst2st1ee.ne.0.D0) then
      write(nout,104) brst2st1ee,3,ist1,ie,-ie,        'BR(~t_2 -> ~t_1
     .    e+   e-)'
      write(nout,104) brst2st1ee,3,ist1,imu,-imu,      'BR(~t_2 -> ~t_1
     .    mu+  mu-)'
      endif
      if(brst2st1tautau.ne.0.D0) then
      write(nout,104) brst2st1tautau,3,ist1,itau,-itau,'BR(~t_2 -> ~t_1
     .    tau+ tau-)'
      endif
      if(brst2st1nunu.ne.0.D0) then
      write(nout,104) brst2st1nunu,3,ist1,ine,-ine,    'BR(~t_2 -> ~t_1
     .    nu_e   nu_eb)'
      write(nout,104) brst2st1nunu,3,ist1,inmu,-inmu,  'BR(~t_2 -> ~t_1
     .    nu_mu  nu_mub)'
      write(nout,104) brst2st1nunu,3,ist1,intau,-intau,'BR(~t_2 -> ~t_1
     .    nu_tau nu_taub)'
      endif
      endif
c --------------- c
c Sbottom1 decays c
c --------------- c
      write(nout,99)
      write(nout,100) 1000005,sbottot(1),'sbottom1'   
      if(sbottot2(1).ne.0.D0) then
        write(nout,49) 'sbottom1 2-body decays'
        write(nout,101)
      if(brsb1neutt(1).ne.0.D0) then
      write(nout,102) brsb1neutt(1),2,in1,ib  ,'BR(~b_1 -> ~chi_10 b )'
      endif
      if(brsb1neutt(2).ne.0.D0) then
      write(nout,102) brsb1neutt(2),2,in2,ib  ,'BR(~b_1 -> ~chi_20 b )'
      endif
      if(brsb1neutt(3).ne.0.D0) then
      write(nout,102) brsb1neutt(3),2,in3,ib  ,'BR(~b_1 -> ~chi_30 b )'
      endif
      if(brsb1neutt(4).ne.0.D0) then
      write(nout,102) brsb1neutt(4),2,in4,ib  ,'BR(~b_1 -> ~chi_40 b )'
      endif
      if(brsb1neutt(5).ne.0.D0) then
      write(nout,102) brsb1neutt(5),2,in5,ib  ,'BR(~b_1 -> ~chi_50 b )'
      endif
      if(brsb1chart(1).ne.0.D0) then
      write(nout,102) brsb1chart(1),2,-ic1,it ,'BR(~b_1 -> ~chi_1- t )'
      endif
      if(brsb1chart(2).ne.0.D0) then
      write(nout,102) brsb1chart(2),2,-ic2,it ,'BR(~b_1 -> ~chi_2- t )'
      endif
      if(brsb1glui.ne.0.D0) then
      write(nout,102) brsb1glui,2,iglo,ib     ,'BR(~b_1 -> ~g      b )'
      endif
      if(brsb1hcst(1).ne.0.D0) then
      write(nout,102) brsb1hcst(1),2,ist1,-ihc,'BR(~b_1 -> ~t_1    H-)'
      endif
      if(brsb1hcst(2).ne.0.D0) then
      write(nout,102) brsb1hcst(2),2,ist2,-ihc,'BR(~b_1 -> ~t_2    H-)'
      endif
      if(brsb1wst(1).ne.0.D0) then
      write(nout,102) brsb1wst(1),2,ist1,-iwc ,'BR(~b_1 -> ~t_1    W-)'
      endif
      if(brsb1wst(2).ne.0.D0) then
      write(nout,102) brsb1wst(2),2,ist2,-iwc ,'BR(~b_1 -> ~t_2    W-)'
      endif
      endif
C   ========================
*   SBOTTOM1: Three body
C   ========================
      if(sbottot3(1).ne.0.D0) then
      write(nout,49) 'sbottom1 3-body decays'
      write(nout,103)
      if(brsbsntau(1,1).ne.0.D0) then
      write(nout,104) brsbsntau(1,1),3,-intau1,it,itau,'BR(~b_1 -> ~nu_t
     .auL* t    tau-)'
      endif
      if(brsbsnel(1).ne.0.D0) then
      write(nout,104) brsbsnel(1),3,-inel,it,ie,       'BR(~b_1 -> ~nu_e
     .L*   t    e-)'
      write(nout,104) brsbsnel(1),3,-inmul,it,imu,     'BR(~b_1 -> ~nu_m
     .uL*  t    mu-)'
      endif
      if(brsbstau(1,1).ne.0.D0) then
      write(nout,104) brsbstau(1,1),3,istau1,it,-intau,'BR(~b_1 -> ~tau_
     .1-   t    nu_taub)'
      endif
      if(brsbstau(1,2).ne.0.D0) then
      write(nout,104) brsbstau(1,2),3,istau2,it,-intau,'BR(~b_1 -> ~tau_
     .2-   t    nu_taub)'
      endif
      if(brsbsel(1,1).ne.0.D0) then
      write(nout,104) brsbsel(1,1),3,isell,it,-ine,    'BR(~b_1 -> ~e_L-
     .     t    nu_eb)'
      endif
      if(brsbsel(1,2).ne.0.D0) then
      write(nout,104) brsbsel(1,2),3,iselr,it,-ine,    'BR(~b_1 -> ~e_R-
     .     t    nu_eb)'
      endif
      if(brsbsel(1,1).ne.0.D0) then
      write(nout,104) brsbsel(1,1),3,ismul,it,-inmu,   'BR(~b_1 -> ~mu_L
     .-    t    nu_mub)'
      endif
      if(brsbsel(1,2).ne.0.D0) then
      write(nout,104) brsbsel(1,2),3,ismur,it,-inmu,   'BR(~b_1 -> ~mu_R
     .-    t    nu_mub)'
      endif
      if(brsbtstsb(1,1).ne.0.D0) then
      write(nout,104) brsbtstsb(1,1),3,-ist1,it,ib,    'BR(~b_1 -> ~t_1*
     .     t    b)'
      endif
      if(brsbtstsb(1,2).ne.0.D0) then
      write(nout,104) brsbtstsb(1,2),3,-ist2,it,ib,    'BR(~b_1 -> ~t_2*
     .     t    b)'
      endif
      if(brsbtbstb(1,1).ne.0.D0) then
      write(nout,104) brsbtbstb(1,1),3,ist1,-it,ib,    'BR(~b_1 -> ~t_1
     .     tb   b)'
      endif
      if(brsbtbstb(1,2).ne.0.D0) then
      write(nout,104) brsbtbstb(1,2),3,ist2,-it,ib,    'BR(~b_1 -> ~t_2
     .     tb   b)'
      endif
      if(brsbupstdow(1,1).ne.0.D0) then
      write(nout,104) brsbupstdow(1,1),3,ist1,iub,id,  'BR(~b_1 -> ~t_1 
     .     ub   d)'
      endif
      if(brsbupstdow(1,2).ne.0.D0) then
      write(nout,104) brsbupstdow(1,2),3,ist2,iub,id,  'BR(~b_1 -> ~t_2 
     .     ub   d)'
      endif
      if(brsbupstdow(1,1).ne.0.D0) then
      write(nout,104) brsbupstdow(1,1),3,ist1,icb,is,  'BR(~b_1 -> ~t_1 
     .     cb   s)'
      endif
      if(brsbupstdow(1,2).ne.0.D0) then
      write(nout,104) brsbupstdow(1,2),3,ist2,icb,is,  'BR(~b_1 -> ~t_2 
     .     cb   s)'
      endif
      if(brsbtaustnu(1,1).ne.0.D0) then
      write(nout,104) brsbtaustnu(1,1),3,ist1,itau,-intau,'BR(~b_1 -> ~t
     ._1      tau- nu_taub)'
      endif
      if(brsbtaustnu(1,2).ne.0.D0) then
      write(nout,104) brsbtaustnu(1,2),3,ist2,itau,-intau,'BR(~b_1 -> ~t
     ._2      tau- nu_taub)'
      endif
      if(brsbelstnu(1,1).ne.0.D0) then
      write(nout,104) brsbelstnu(1,1),3,ist1,ie,-ine,  'BR(~b_1 -> ~t_1
     .     e-   nu_eb)'
      endif
      if(brsbelstnu(1,2).ne.0.D0) then
      write(nout,104) brsbelstnu(1,2),3,ist1,ie,-ine,  'BR(~b_1 -> ~t_2
     .     e-   nu_eb)'
      endif
      if(brsbelstnu(1,1).ne.0.D0) then
      write(nout,104) brsbelstnu(1,1),3,ist1,imu,-inmu,'BR(~b_1 -> ~t_1
     .     mu-  nu_mub)'
      endif
      if(brsbelstnu(1,2).ne.0.D0) then
      write(nout,104) brsbelstnu(1,2),3,ist1,imu,-inmu,'BR(~b_1 -> ~t_2
     .     mu-  nu_mub)'
      endif
      endif
c --------------- c
c Sbottom2 decays c
c --------------- c
       write(nout,99)
      write(nout,100) 2000005,sbottot(2),'sbottom2'
      if(sbottot2(2).ne.0.D0) then
        write(nout,49) 'sbottom2 2-body decays'
        write(nout,101)
      if(brsb2neutt(1).ne.0.D0) then
      write(nout,102) brsb2neutt(1),2,in1,ib  ,'BR(~b_2 -> ~chi_10 b )'
      endif
      if(brsb2neutt(2).ne.0.D0) then
      write(nout,102) brsb2neutt(2),2,in2,ib  ,'BR(~b_2 -> ~chi_20 b )'
      endif
      if(brsb2neutt(3).ne.0.D0) then
      write(nout,102) brsb2neutt(3),2,in3,ib  ,'BR(~b_2 -> ~chi_30 b )'
      endif
      if(brsb2neutt(4).ne.0.D0) then
      write(nout,102) brsb2neutt(4),2,in4,ib  ,'BR(~b_2 -> ~chi_40 b )'
      endif
      if(brsb2neutt(5).ne.0.D0) then
      write(nout,102) brsb2neutt(5),2,in5,ib  ,'BR(~b_2 -> ~chi_50 b )'
      endif
      if(brsb2chart(1).ne.0.D0) then
      write(nout,102) brsb2chart(1),2,-ic1,it ,'BR(~b_2 -> ~chi_1- t )'
      endif
      if(brsb2chart(2).ne.0.D0) then
      write(nout,102) brsb2chart(2),2,-ic2,it ,'BR(~b_2 -> ~chi_2- t )'
      endif
      if(brsb2glui.ne.0.D0) then
      write(nout,102) brsb2glui,2,iglo,ib     ,'BR(~b_2 -> ~g      b )'
      endif
      if(brsb2H(1).ne.0.D0) then
      write(nout,102) brsb2H(1),2,isb1,ihH1   ,'BR(~b_2 -> ~b_1    H1 )'
      endif
      if(brsb2H(2).ne.0.D0) then
      write(nout,102) brsb2H(2),2,isb1,ihH2   ,'BR(~b_2 -> ~b_1    H2 )'
      endif
      if(brsb2H(3).ne.0.D0) then
      write(nout,102) brsb2H(3),2,isb1,ihH3   ,'BR(~b_2 -> ~b_1    H3 )'
      endif
      if(brsb2A(1).ne.0.D0) then
      write(nout,102) brsb2A(1)  ,2,isb1,ihA1,'BR(~b_2 -> ~b_1    A1 )'
      endif
      if(brsb2A(2).ne.0.D0) then
      write(nout,102) brsb2A(2)  ,2,isb1,ihA2,'BR(~b_2 -> ~b_1    A2 )'
      endif
      if(brsb2hcst(1).ne.0.D0) then
      write(nout,102) brsb2hcst(1),2,ist1,-ihc,'BR(~b_2 -> ~t_1    H-)'
      endif
      if(brsb2hcst(2).ne.0.D0) then
      write(nout,102) brsb2hcst(2),2,ist2,-ihc,'BR(~b_2 -> ~t_2    H-)'
      endif
      if(brsb2zbot.ne.0.D0) then
      write(nout,102) brsb2zbot,2,isb1,iz     ,'BR(~b_2 -> ~b_1    Z )'
      endif
      if(brsb2wst(1).ne.0.D0) then
      write(nout,102) brsb2wst(1),2,ist1,-iwc ,'BR(~b_2 -> ~t_1    W-)'
      endif
      if(brsb2wst(2).ne.0.D0) then
      write(nout,102) brsb2wst(2),2,ist2,-iwc ,'BR(~b_2 -> ~t_2    W-)'
      endif
      endif
c
C   ========================
*   SBOTTOM2: Three body
C   ========================
      if(sbottot3(2).ne.0.D0) then
      write(nout,49) 'sbottom2 3-body decays'
      write(nout,103)
      if(brsbsntau(2,1).ne.0.D0) then
      write(nout,104) brsbsntau(2,1),3,-intau1,it,itau,'BR(~b_2 -> ~nu_t
     .auL* t      tau-)'
      endif
      if(brsbsnel(2).ne.0.D0) then
      write(nout,104) brsbsnel(2),3,-inel,it,ie,       'BR(~b_2 -> ~nu_e
     .L*   t      e-)'
      write(nout,104) brsbsnel(2),3,-inmul,it,imu,     'BR(~b_2 -> ~nu_m
     .uL*  t      mu-)'
      endif
      if(brsbstau(2,1).ne.0.D0) then
      write(nout,104) brsbstau(2,1),3,istau1,it,-intau,'BR(~b_2 -> ~tau_
     .1-   t      nu_taub)'
      endif
      if(brsbstau(2,2).ne.0.D0) then
      write(nout,104) brsbstau(2,2),3,istau2,it,-intau,'BR(~b_2 -> ~tau_
     .2-   t      nu_taub)'
      endif
      if(brsbsel(2,1).ne.0.D0) then
      write(nout,104) brsbsel(2,1),3,isell,it,-ine,    'BR(~b_2 -> ~e_L-
     .     t      nu_eb)'
      endif
      if(brsbsel(2,2).ne.0.D0) then
      write(nout,104) brsbsel(2,2),3,iselr,it,-ine,    'BR(~b_2 -> ~e_R-
     .     t      nu_eb)'
      endif
      if(brsbsel(2,1).ne.0.D0) then
      write(nout,104) brsbsel(2,1),3,ismul,it,-inmu,   'BR(~b_2 -> ~mu_L
     .-    t      nu_mub)'
      endif
      if(brsbsel(2,2).ne.0.D0) then
      write(nout,104) brsbsel(2,2),3,ismur,it,-inmu,   'BR(~b_2 -> ~mu_R
     .-    t      nu_mub)'
      endif
      if(brsbtstsb(2,1).ne.0.D0) then
      write(nout,104) brsbtstsb(2,1),3,-ist1,it,ib,    'BR(~b_2 -> ~t_1*
     .     t      b)'
      endif
      if(brsbtstsb(2,2).ne.0.D0) then
      write(nout,104) brsbtstsb(2,2),3,-ist2,it,ib,    'BR(~b_2 -> ~t_2*
     .     t      b)'
      endif
      if(brsbtbstb(2,1).ne.0.D0) then
      write(nout,104) brsbtbstb(2,1),3,ist1,-it,ib,    'BR(~b_2 -> ~t_1
     .     tb     b)'
      endif
      if(brsbtbstb(2,2).ne.0.D0) then
      write(nout,104) brsbtbstb(2,2),3,ist2,-it,ib,    'BR(~b_2 -> ~t_2
     .     tb     b)'
      endif
      if(brsbupstdow(2,1).ne.0.D0) then
      write(nout,104) brsbupstdow(2,1),3,ist1,iub,id,  'BR(~b_2 -> ~t_1 
     .     ub     d)'
      endif
      if(brsbupstdow(2,2).ne.0.D0) then
      write(nout,104) brsbupstdow(2,2),3,ist2,iub,id,  'BR(~b_2 -> ~t_2 
     .     ub     d)'
      endif
      if(brsbupstdow(2,1).ne.0.D0) then
      write(nout,104) brsbupstdow(2,1),3,ist1,icb,is,  'BR(~b_2 -> ~t_1 
     .     cb     s)'
      endif
      if(brsbupstdow(2,2).ne.0.D0) then
      write(nout,104) brsbupstdow(2,2),3,ist2,icb,is,  'BR(~b_2 -> ~t_2 
     .     cb     s)'
      endif
      if(brsbtaustnu(2,1).ne.0.D0) then
      write(nout,104) brsbtaustnu(2,1),3,ist1,itau,-intau,'BR(~b_2 -> ~t
     ._1      tau-   nu_taub)'
      endif
      if(brsbtaustnu(2,2).ne.0.D0) then
      write(nout,104) brsbtaustnu(2,2),3,ist2,itau,-intau,'BR(~b_2 -> ~t
     ._2      tau-   nu_taub)'
      endif
      if(brsbelstnu(2,1).ne.0.D0) then
      write(nout,104) brsbelstnu(2,1),3,ist1,ie,-ine,  'BR(~b_2 -> ~t_1
     .     e-     nu_eb)'
      endif
      if(brsbelstnu(2,2).ne.0.D0) then
      write(nout,104) brsbelstnu(2,2),3,ist1,ie,-ine,  'BR(~b_2 -> ~t_2
     .     e-     nu_eb)'
      endif
      if(brsbelstnu(2,1).ne.0.D0) then
      write(nout,104) brsbelstnu(2,1),3,ist1,imu,-inmu,'BR(~b_2 -> ~t_1
     .     mu-    nu_mub)'
      endif
      if(brsbelstnu(2,2).ne.0.D0) then
      write(nout,104) brsbelstnu(2,2),3,ist1,imu,-inmu,'BR(~b_2 -> ~t_2
     .     mu-    nu_mub)' 
      endif
      if(brsb2sb1bb.ne.0.D0) then
      write(nout,104) brsb2sb1bb,3,isb1,ib,ibb,        'BR(~b_2 -> ~b_1 
     .     b      bb)'
      endif
      if(brsb2sb1starbb.ne.0.D0) then
      write(nout,104) brsb2sb1starbb,3,-isb1,ib,ib,    'BR(~b_2 -> ~b_1*
     .     b      b)'
      endif
      if(brsb2sb1tt.ne.0.D0) then
      write(nout,104) brsb2sb1tt,3,isb1,it,itb,        'BR(~b_2 -> ~b_1 
     .     t      tb)'
      endif
      if(brsb2sb1uu.ne.0.D0) then
      write(nout,104) brsb2sb1uu,3,isb1,iu,iub,        'BR(~b_2 -> ~b_1 
     .     u      ub)'
      endif
      if(brsb2sb1dd.ne.0.D0) then
      write(nout,104) brsb2sb1dd,3,isb1,id,idb,        'BR(~b_2 -> ~b_1 
     .     d      db)'
      endif
      if(brsb2sb1uu.ne.0.D0) then
      write(nout,104) brsb2sb1uu,3,isb1,ic,icb,        'BR(~b_2 -> ~b_1 
     .     c      cb)'
      endif
      if(brsb2sb1dd.ne.0.D0) then
      write(nout,104) brsb2sb1dd,3,isb1,is,isb,        'BR(~b_2 -> ~b_1 
     .     s      sb)'
      endif
      if(brsb2sb1ee.ne.0.D0) then
      write(nout,104) brsb2sb1ee,3,isb1,ie,-ie,        'BR(~b_2 -> ~b_1 
     .     e-     e+)'
      write(nout,104) brsb2sb1ee,3,isb1,imu,-imu,      'BR(~b_2 -> ~b_1 
     .     mu-    mu+)'
      endif
      if(brsb2sb1tautau.ne.0.D0) then
      write(nout,104) brsb2sb1tautau,3,isb1,itau,-itau,'BR(~b_2 -> ~b_1 
     .     tau-   tau+)'
      endif
      if(brsb2sb1nunu.ne.0.D0) then
      write(nout,104) brsb2sb1nunu,3,isb1,ine,-ine,    'BR(~b_2 -> ~b_1 
     .     nu_e   nu_eb)'
      write(nout,104) brsb2sb1nunu,3,isb1,inmu,-inmu,  'BR(~b_2 -> ~b_1 
     .     nu_mu  nu_mub)'
      write(nout,104) brsb2sb1nunu,3,isb1,intau,-intau,'BR(~b_2 -> ~b_1 
     .     nu_tau nu_taub)'
      endif
      endif
c ==================================================================== c
c                       end of the output file                          
c ==================================================================== c

 49   format('#',20x,A,E16.8)
 99   format('#',9x,'PDG',12x,'Width')
 100  format('DECAY',1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 101  format('#',10x,'BR',9x,'NDA',6x,'ID1',7x,'ID2')
 102  format(3x,1P,E16.8,0P,3x,I2,3x,(I9,1x),(I9,1x),2x,'#',1x,A)
 103  format('#',11x,'BR',9x,'NDA',6x,'ID1',7x,'ID2',7x,'ID3')
 104  format(3x,1P,E16.8,0P,3x,I2,3x,(I9,1x),(I9,1x),(I9,1x),2x,'#',
     .1x,A)

      end
