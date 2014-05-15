c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                    The calling program suspect_call.f 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  VERSION 2.41
c  Last changes : August 1, 2008 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  For the users manual, updated information, changes, maintenance, see 
c  Home page: http://w3.lpta.univ-montp2.fr/~kneur/Suspect/
c
c  This program is the example routine calling the main program SuSpect.f.
c  It has to be compiled together with suspect2.f in all cases, but it is
c  particularly useful when performing e.g. a scan of the parameter space   
c  and/or to interface SuSpect with another program. In this routine you have  
c  to set the four control parameters which are the inputs arguments of the  
c  main program:
c         SUBROUTINE SuSpect2(iknowl,input,ichoice,errmess) 
c  The input are (see details in the comments beginning of SuSpect2.f): 
c  IKNOWL: sets degree of control on various parts of the algorithm:
c  =0: blind use of the program, no control on parameters and no warning.
c  =1:  warning/error messages in output file
c  INPUT: sets input (and output) control, it offers now 4 possibilities:
c  =0: model and option input parameters ONLY read in file suspect2.in.
c  (output generated in both suspect2.out and SLHA format suspect2_lha.out) 
c  =1: define yourself IN THIS FILE all relevant input and parameters.
c  (i.e. NO reading of input files): see example of  input models below.
c  Maybe more convenient e.g. for scan over the model parameter space.
c  (output generated in both suspect2.out and SLHA format suspect2_lha.out) 
c  =2: same as input =0 but read SLHA format input file: suspect2_lha.in
c (it writes also all output in the SLHA format file: suspect2_lha.out)
c  =11: same as input=1, but NO output file(s) suspect*.out generated
c
c  ICHOICE: initialises various model/accuracy options to be considred: 
c  - ICHOICE(1): Choice of the model to be considered.
c  - ICHOICE(2): For the perturbative order (1 or 2 loop) of the RGEs. 
c  - ICHOICE(3): To impose or not the GUT scale. 
c  - ICHOICE(4): For the accuracy of the RGEs.
c  - ICHOICE(5): To impose or not the radiative EWSB. 
c  - ICHOICE(6): To chose different (scalr sector) input in general MSSM.
c  - ICHOICE(7): For the radiative corrections to the (s)particles masses. 
c  - ICHOICE(8): To set the value of the EWSB scale.
c  - ICHOICE(9): For the number of (long: RGE + full spectrum) iterations: 
c  - ICHOICE(10): For the routine calculating the Higgs boson masses.
c  - ICHOICE(11): A new option for the (higher order) Higgs mass R.C. scheme 
c  ERRMESS provides a useful set of warning/error message flags in output file.
c  - ERRMESS(i)= 0: Everything is fine.
c  - ERRMESS(1)=-1: tachyon 3rd gen. sfermion from RGE
c  - ERRMESS(2)=-1: tachyon 1,2 gen. sfermion from RGE
c  - ERRMESS(3)=-1: tachyon A    (maybe temporary: see final mass) 
c  - ERRMESS(4)=-1: tachyon 3rd gen. sfermion from mixing
c  - ERRMESS(5)=-1: mu(M_GUT) guess inconsistent 
c  - ERRMESS(6)=-1: non-convergent mu from EWSB 
c  - ERRMESS(7)=-1: EWSB maybe inconsistent  (!but RG-improved only check)
c  - ERRMESS(8)=-1: V_Higgs maybe UFB or CCB (!but RG-improved only check)
c  - ERRMESS(9)=-1: Higgs boson masses are NaN 
c  - ERRMESS(10)=-1: RGE problems (non-pert and/or Landau poles)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                      The program starts here
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PROGRAM main
c 
      implicit real*8(a-h,m,o-z)
      real*8 nl,nq
      dimension ichoice(11),errmess(10)
      dimension u(2,2),vv(2,2),z(4,4),xmn(4)
c
c
c=========================  COMMONs for inputs ============================
c   These are the commons for the parameters that can be read in the file 
c   suspect2.in (together with the various ichoices). 
c   !Important note: to interface your program with SuSpect2.f, these 
c    commons (plus the output ones below) are the only ones needed.
c
c   "Standard model" INPUT parameters (couplings and fermion masses):
       COMMON/SU_SMPAR/alfinv,sw2,alphas,mt,mb,mc,mtau
c  !MODIF!   mt,mtau pole masses; while mb is NOW mb(mb)_MSbar
c   RG evolution scale parameters (EWSB scale, high and low RGE ends):
       COMMON/SU_RGSCAL/qewsb,ehigh,elow
c   MSSM scalar sector parameters:     
       COMMON/SU_MSSMHPAR/mhu2,mhd2,ma,mu
c   The U(1), SU(2), SU(3) SUSY-breaking gaugino masses
       COMMON/SU_MSSMGPAR/m1,m2,m3 
c   The soft-SUSY breaking slepton mass terms (3d and then 1/2 gen.): 
       COMMON/SU_MSSMSLEP/msl,mtaur,mel,mer
c   The soft-SUSY breaking squark mass terms (3d and then 1/2 gen.):
       COMMON/SU_MSSMSQUA/msq,mtr,mbr,muq,mur,mdr
c   The soft-SUSY breaking trilinear couplings (3d and then 1/2 gen.):
       COMMON/SU_ATRI3/atau,at,ab
       COMMON/SU_ATRI12/al,au,ad
c   mSUGRA case input parameters:
       COMMON/SU_MSUGRA/m0,mhalf,a0
       COMMON/SU_RADEWSB/sgnmu0,tgbeta
c   GMSB case input parameters:
       COMMON/SU_GMSB/mgmmess,mgmsusy,nl,nq
c   AMSB case input parameters:
       COMMON/SU_AMSB/m32,am0,cq,cu,cd,cl,ce,chu,chd
c
c========================  COMMONs for output ============================= 
c
      COMMON/SU_OUTHIGGS/ml,mh,mch,alfa
c   light, heavy, charged Higgs masses, Higgs mix angle alpha 
       COMMON/SU_OUTGINOS/mc1,mc2,mn1,mn2,mn3,mn4,gluino
c   charginos 1,2 masses, neutralinos 1-4 masses, gluino mass 
       COMMON/SU_OUTSQU/mst1,mst2,msu1,msu2
c   stop 1,2 and sup 1,2 = scharm 1,2 masses
       COMMON/SU_OUTSQD/msb1,msb2,msd1,msd2
c   sbottom 1,2 and sdown 1,2 = sstrange 1,2 masses
       COMMON/SU_OUTSLEP/msl1,msl2,mse1,mse2,msn1,msntau
c   stau 1,2 ; selectron (=smuon) 1,2; sneut_e,mu, sneut_tau masses
       COMMON/SU_OUTMIX/thet,theb,thel
c   stop, sbottom, stau mixing angles
       COMMON/SU_MATINO/u,vv,z,xmn
c   U,V chargino and Z neutralino diagonalizing matrices and masses
      COMMON/SU_YUKAEWSB/ytauewsb,ybewsb,ytewsb,alsewsb,g2ewsb,g1ewsb
c   (final) bottom, top tau masses and gauge couplings at EWSB scale 
      common/su_runmewsb/mtaurewsb,mbrewsb,mtrewsb
c  (running masses of tau,b,top in DRbar scheme at Q= EWSB scale: 
c   may be useful as output)
c low-energy contrained parameter values: rho-1, g_mu-2, Br(b->s gamma):
      COMMON/SU_lowen/crho,gmuon,brsg
c Note for soft terms: OUTPUT values are contained in the same commons as
c for their input values: su_MSSMhpar, su_MSSMgpar, su_MSSMslep,
c                         su_MSSMsqua, su_Atri3, su_Atri12
c 
c (NB if meaning ambiguous, see detailed parameter definitions and 
c  conventions in suspect2.f file)
c
c========================  Setting the running input ======================
c   Here you set your command for reading the input SuSpect2.in and/or writing 
c   (see functions in the comments above) in the output file SuSpect.out
       IKNOWL = 1
       INPUT  = 2
c  Now 4 possible choices: 
c  INPUT = 0: input read in original SuSpect format file suspect2.in
c  (Output generated in both suspect2.out and SLHA format suspect2_lha.out) 
c
c  INPUT = 1: all input parameters are defined in this calling file
c (convenient to define e.g. loops on input parameters for scanning etc:
c  in this case user should also set ICHOICE(1)-(11) and all other input 
c  parameters, see mSUGRA, GMSB, AMSB or general MSSM examples below. 
c  Output generated in both suspect2.out and SLHA format suspect2_lha.out) 
c
c  INPUT = 2: same as input =0 but SLHA format INPUT file suspect2_lha.in
c  (output generated in both suspect2.out and SLHA format suspect2_lha.out) 
c  ! in this case any parameter defined below in this file is ignored !
c  INPUT = 11: Same as INPUT = 1, but NO OUTPUT File(s) generated
c (convenient e.g. for scan on MSSM/mSUGRA parameters)
c=======  Example of choice for the model/accuracy, etc, parameters  
      if(input.eq.0.or.input.eq.2) goto 99
c control parameters input:
c ! only relevant if INPUT=1 or 11 above, i.e. if input read from this file
c (For similar model choices with SLHA conventions use suspect2_lha.in)
      ichoice(1) = 10
c ICHOICE(1): Choice of the model: ! old Suspect input conventions are:
c            Arbitrary soft-terms at low scale         : 0
c            Arbitrary soft-terms at high scale        : 1
c            SUGRA  (cMSSM)                            : 10
c            GMSB   (cMSSM)                            : 11
c            AMSB   (cMSSM)                            : 12
c	    bottom-up RGE (arbitrary MSSM input at low energy): 2
c
      ichoice(2) = 21
c (2-loop RGE for gauge, yukawas, gauginos)
      ichoice(3) = 1
c (ichoice(3)= 0: GUT scale imposed (then EHIGH = input!); 
c            = 1: gauge unif scale calculated consistently from
c              gauge couplings input (RECOMMANDED choice!)
c
      ichoice(4) = 2
c  (RG accuracy: 1: moderately accurate and fast (generally sufficient)
c                 2: very accurate but rather slow!
      ichoice(5) = 1
c (consistent EWSB)
      ichoice(6) = 1
c (M_Hu, M_Hd (= m_0 in mSUGRA) input)
      ichoice(7) = 2
c ICHOICE(7):  SUSY radiative corrections to the (s)particles masses: 
c                      No Radiative corrections     : 0 
c                      only in mb,mt,mtau +Yukawas : 1  
c                      all squarks + gaugino R.C. in addition: 2
      ichoice(8) = 1
c (ichoice(8) = 1 for default EWSB scale=(m_t_L*m_t_R)^(1/2), =0 if not)
c then IF ichoice(8)=0 EWSB scale is set by user from input file/calling
c routine by the value of Qewsb
      ichoice(9) = 2
c !MODIFS! ichoice(9) =1: 1%    accuracy in final spectrum calculations;
c                     =2  0.01% accuracy  
      ichoice(10) = 2
c ICHOICE(10): Higgs mass Rad. Corr. calculation options:
c       approximate m_h,H,A calculation (fast but not accurate)  : 0
c       Full one-loop Higgs R.C. a la PBMZ                       : 1
c       Full one-loop a la PBMZ + dominant 2-loop (P. Slavich et al):2
c 
c                      
      ichoice(11) = 0
c ICHOICE(11): ! higher order Higgs 'scheme' in rad. corr. at mZ:
c              (NB formally a higher (2-loop) order difference)
c          RUNNING DRbar Higgs masses at loop-level at mZ (preferred!): 0
c          POLE          Higgs masses at loop-level at mZ             : 1
c
c======= Then define the needed SM and SUSY input parameters (example below)
c        (these are the parameters contained in the commons for input above):
c  "SM-like" input:
      alfinv= 127.934d0
      alphas =.1172d0
      mt =175.d0
      mb = 4.25d0     ! mb(mb)_MSbar input
      mtau =1.777d0
c
c RG evolution parameters: 
c      ehigh = 1.9d16
c      elow =91.19d0
c      qewsb = 200.d0
c (!! qewsb value only relevant if ichoice(8) = 0, see above)
c   
c  minimal SUGRA case input sample (SPS4 here):
c
      m0 = 100.d0
      mhalf = 250.d0
      A0 = -100.d0
      sgnmu0 = 1.d0
      tgbeta = 10.d0
c
c   
c  Gauge Mediated Supersymmetry Breaking (GMSB) input sample (SPS8):
c (simply uncomment the CC if needed):
CC      mgmmess = 200.d3
CC      mgmsusy = 100.d3
CC      nl = 1
CC      nq = 1
CC      tgbeta = 15.d0
CC      sgnmu0 = 1.d0
c
c   
c  Anomaly Mediated Supersymmetry Breaking (AMSB) input sample (SPS9):
c (simply uncomment the CC if needed):
CC      m32 = 60.d3
CC      am0  = 450.d0
CC      tgbeta = 10.d0
CC      sgnmu0 = 1.d0
CC      cq  =1.d0
CC      cu = 1.d0
CC      cd = 1.d0
CC      cl = 1.d0
CC      ce = 1.d0
CC      chu =1.d0
CC      chd =1.d0
c
      if(ichoice(1).ge.10) goto 99
c   non-universal case input sample: (here in fact it's msugra sps1a)
c (simply uncomment the CC if needed):
        tgbeta = 10d0
        MHU2     = 1.d4
        MHD2     = 1.d4
        M1       = 250.D0
        M2       = 250.D0
        M3       = 250.D0
        MSL      = 1.D2
        MTAUR    = 1.D2
        MSQ      = 1.d2
        MTR      = 1.d2
        MBR      = 1.d2
c
        MEL      = 1.D2
        MER      = 1.D2
        MUQ      = 1.D2
        MUR      = 1.D2
        MDR      = 1.D2
        Atau     = -1.D2
        At       = -100.d0
        Ab       = -1.D2
        AL       = -100.d0
        AU       = -100.d0
        AD       = -100.d0
c case of MA, MU input (instead of MHU2, MHD2) (set ichoice(6)=0):
cc        MA       = 1000.d0
cc        MU       = 1000.d0
c      
  99  continue
c    At this stage you can call the main subroutine suspect:
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
       CALL suspect2(IKNOWL, INPUT, ICHOICE,ERRMESS)
c
c  (ALL relevant OUTPUT will be written in suspect.out file;
c   except if INPUT=11 chosen, and you may continue with output 
c   values within this program)
c
c   ......sequel of your own program continues e.g. here 
c   ....
c   ... Bonne route!
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
