#
#                              ======================
#                              | THE SUSYHIT OUTPUT |
#                              ======================
#
#
#              ------------------------------------------------------
#              |     This is the output of the SUSY-HIT package     |
#              |  created by A.Djouadi, M.Muehlleitner and M.Spira. |
#              |  In case of problems with SUSY-HIT email to        |
#              |           margarete.muehlleitner@cern.ch           |
#              |           michael.spira@psi.ch                     |
#              |           abdelhak.djouadi@cern.ch                 |
#              ------------------------------------------------------
#
#              ------------------------------------------------------
#              |  SUSY Les Houches Accord - MSSM Spectrum + Decays  |
#              |              based on the decay programs           |
#              |                                                    |
#              |                     SDECAY 1.3b                    |
#              |                                                    |
#              |  Authors: M.Muhlleitner, A.Djouadi and Y.Mambrini  |
#              |  Ref.:    Comput.Phys.Commun.168(2005)46           |
#              |           [hep-ph/0311167]                         |
#              |                                                    |
#              |                     HDECAY 3.4                     |
#              |                                                    |
#              |  By: A.Djouadi,J.Kalinowski,M.Muhlleitner,M.Spira  |
#              |  Ref.:    Comput.Phys.Commun.108(1998)56           |
#              |           [hep-ph/9704448]                         |
#              |                                                    |
#              |                                                    |
#              |  If not stated otherwise all DRbar couplings and   |
#              |  soft SUSY breaking masses are given at the scale  |
#              |  Q=  0.10000000E+01
#              |                                                    |
#              ------------------------------------------------------
#
#
BLOCK DCINFO  # Decay Program information
     1   SDECAY/HDECAY # decay calculator
     2   1.3b  /3.4    # version number
#
BLOCK SPINFO  # Spectrum calculator information
     1   SOFTSUSY    # spectrum calculator (modified by MM
     2   3.1.6         # version number
#
BLOCK MODSEL  # Model selection
     1    10   # # tgm
#
BLOCK SMINPUTS  # Standard Model inputs
         1     1.27925000E+02   # alpha_em-1(M_Z)^MSbar
         2     1.16637000E-05   # G_F [GeV-2]
         3     1.17600000E-01   # alpha_S(M_Z)^MSbar
         4     9.11876000E+01   # M_Z pole mass
         5     4.20000000E+00   # mb(mb)^MSbar
         6     1.73300000E+02   # mt pole mass
         7     1.77699000E+00   # mtau pole mass
#
BLOCK MINPAR  # Input parameters - minimal models
         1     5.00000000E+02   # m0
         2     1.00000000E+04   # m12
         3     5.00000000E+00   # tanb
         4     1.00000000E+00   # sign(mu)
         5     8.00000000E-01   # cosu
         6     9.00000000E-01   # cosd
         7     1.00000000E+00   # cd
         8     1.00000000E+00   # cs
         9     1.00000000E+00   # cb
        10     1.00000000E+00   # ce
        11     1.00000000E+00   # cmu
        12     1.00000000E+00   # ctau
        13    -8.00000000E-01   # k
#
BLOCK MASS  # Mass Spectrum
# PDG code           mass       particle
        24     8.04863425E+01   # W+
        25     1.01961581E+02   # h
        35     7.57239970E+02   # H
        36     7.56640495E+02   # A
        37     7.60997634E+02   # H+
         5     4.82976602E+00   # b-quark pole mass calculated from mb(mb)_Msbar
   1000001     6.77477549E+02   # ~d_L
   2000001     8.28587252E+02   # ~d_R
   1000002     6.73219167E+02   # ~u_L
   2000002     6.49162879E+02   # ~u_R
   1000003     6.77477549E+02   # ~s_L
   2000003     8.28587252E+02   # ~s_R
   1000004     6.73219167E+02   # ~c_L
   2000004     6.49162879E+02   # ~c_R
   1000005     6.11231771E+02   # ~b_1
   2000005     8.27456950E+02   # ~b_2
   1000006     5.04395007E+02   # ~t_1
   2000006     6.52611275E+02   # ~t_2
   1000011     7.22099642E+02   # ~e_L
   2000011     5.27369677E+02   # ~e_R
   1000012     7.17670941E+02   # ~nu_eL
   1000013     7.22099642E+02   # ~mu_L
   2000013     5.27369677E+02   # ~mu_R
   1000014     7.17670941E+02   # ~nu_muL
   1000015     5.25835388E+02   # ~tau_1
   2000015     7.21624487E+02   # ~tau_2
   1000016     7.17127029E+02   # ~nu_tauL
   1000021     4.85274639E+02   # ~g
   1000022     7.34032753E+01   # ~chi_10
   1000023     1.42236693E+02   # ~chi_20
   1000025    -5.33209842E+02   # ~chi_30
   1000035     5.42612034E+02   # ~chi_40
   1000024     1.42170758E+02   # ~chi_1+
   1000037     5.42446306E+02   # ~chi_2+
   1000039     7.50000000E+00   # ~gravitino
#
BLOCK NMIX  # Neutralino Mixing Matrix
  1  1     9.93686529E-01   # N_11
  1  2    -5.18832253E-02   # N_12
  1  3     9.43103999E-02   # N_13
  1  4    -3.16348012E-02   # N_14
  2  1     6.85309653E-02   # N_21
  2  2     9.82550023E-01   # N_22
  2  3    -1.57874334E-01   # N_23
  2  4     7.05312251E-02   # N_24
  3  1    -4.05813500E-02   # N_31
  3  2     6.45553040E-02   # N_32
  3  3     7.00954427E-01   # N_33
  3  4     7.09118226E-01   # N_34
  4  1    -7.90173529E-02   # N_41
  4  2     1.66541875E-01   # N_42
  4  3     6.89089351E-01   # N_43
  4  4    -7.00839445E-01   # N_44
#
BLOCK UMIX  # Chargino Mixing Matrix U
  1  1     9.73855900E-01   # U_11
  1  2    -2.27166650E-01   # U_12
  2  1     2.27166650E-01   # U_21
  2  2     9.73855900E-01   # U_22
#
BLOCK VMIX  # Chargino Mixing Matrix V
  1  1     9.94619798E-01   # V_11
  1  2    -1.03592743E-01   # V_12
  2  1     1.03592743E-01   # V_21
  2  2     9.94619798E-01   # V_22
#
BLOCK STOPMIX  # Stop Mixing Matrix
  1  1     4.27359702E-01   # cos(theta_t)
  1  2     9.04081681E-01   # sin(theta_t)
  2  1    -9.04081681E-01   # -sin(theta_t)
  2  2     4.27359702E-01   # cos(theta_t)
#
BLOCK SBOTMIX  # Sbottom Mixing Matrix
  1  1     9.99710333E-01   # cos(theta_b)
  1  2     2.40676311E-02   # sin(theta_b)
  2  1    -2.40676311E-02   # -sin(theta_b)
  2  2     9.99710333E-01   # cos(theta_b)
#
BLOCK STAUMIX  # Stau Mixing Matrix
  1  1     1.94849200E-02   # cos(theta_tau)
  1  2     9.99810151E-01   # sin(theta_tau)
  2  1    -9.99810151E-01   # -sin(theta_tau)
  2  2     1.94849200E-02   # cos(theta_tau)
#
BLOCK ALPHA  # Higgs mixing
          -2.08956351E-01   # Mixing angle in the neutral Higgs boson sector
#
BLOCK HMIX Q=  1.00000000E+00  # DRbar Higgs Parameters
         1     4.99439163E+02   # mu(Q)MSSM
         2     5.42132999E+00   # tan
         3     2.65697615E+02   # higgs
         4     7.96472235E+05   # mA2(Q)MSSM
#
BLOCK GAUGE Q=  1.00000000E+00  # The gauge couplings
     1     3.41675669E-01   # gprime(Q) DRbar
     2     6.34237384E-01   # g(Q) DRbar
     3     1.29538165E+00   # g3(Q) DRbar
#
BLOCK AU Q=  1.00000000E+00  # The trilinear couplings
  1  1    -6.22310114E+02   # A_u(Q) DRbar
  2  2    -6.22306847E+02   # A_c(Q) DRbar
  3  3    -4.57760140E+02   # A_t(Q) DRbar
#
BLOCK AD Q=  1.00000000E+00  # The trilinear couplings
  1  1    -7.80731412E+02   # A_d(Q) DRbar
  2  2    -7.80729617E+02   # A_s(Q) DRbar
  3  3    -7.24593856E+02   # A_b(Q) DRbar
#
BLOCK AE Q=  1.00000000E+00  # The trilinear couplings
  1  1    -9.47229972E+01   # A_e(Q) DRbar
  2  2    -9.47225105E+01   # A_mu(Q) DRbar
  3  3    -9.45801923E+01   # A_tau(Q) DRbar
#
BLOCK Yu Q=  1.00000000E+00  # The Yukawa couplings
  1  1     0.00000000E+00   # y_u(Q) DRbar
  2  2     0.00000000E+00   # y_c(Q) DRbar
  3  3     1.02749203E+00   # y_t(Q) DRbar
#
BLOCK Yd Q=  1.00000000E+00  # The Yukawa couplings
  1  1     0.00000000E+00   # y_d(Q) DRbar
  2  2     0.00000000E+00   # y_s(Q) DRbar
  3  3     9.66689866E-02   # y_b(Q) DRbar
#
BLOCK Ye Q=  1.00000000E+00  # The Yukawa couplings
  1  1     0.00000000E+00   # y_e(Q) DRbar
  2  2     0.00000000E+00   # y_mu(Q) DRbar
  3  3     5.36194067E-02   # y_tau(Q) DRbar
#
BLOCK MSOFT Q=  1.00000000E+00  # The soft SUSY breaking masses at the scale Q
         1     6.82299730E+01   # M_1(Q)
         2     1.39516728E+02   # M_2(Q)
         3     5.98081375E+02   # M_3(Q)
        21     2.75357101E+05   # mH12(Q)
        22    -3.70879043E+05   # mH22(Q)
        31     7.20877420E+02   # meL(Q)
        32     7.20875080E+02   # mmuL(Q)
        33     7.20190624E+02   # mtauL(Q)
        34     5.26402516E+02   # meR(Q)
        35     5.26396107E+02   # mmuR(Q)
        36     5.24519190E+02   # mtauR(Q)
        41     7.72433108E+02   # mqL1(Q)
        42     7.72430919E+02   # mqL2(Q)
        43     6.89610404E+02   # mqL3(Q)
        44     7.49137691E+02   # muR(Q)
        45     7.49134559E+02   # mcR(Q)
        46     5.68243358E+02   # mtR(Q)
        47     9.08151561E+02   # mdR(Q)
        48     9.08150422E+02   # msR(Q)
        49     9.06015961E+02   # mbR(Q)
#
#
#
#                             =================
#                             |The decay table|
#                             =================
#
# - The QCD corrections to the decays gluino -> squark  + quark
#                                     squark -> gaugino + quark_prime
#                                     squark -> squark_prime + Higgs
#                                     squark -> gluino  + quark
#   are included.
#
# - The multi-body decays for the inos, stops and sbottoms are included.
#
# - The loop induced decays for the gluino, neutralinos and stops
#   are included.
#
# - The SUSY decays of the top quark are included.
#
# - Possible decays of the NLSP in GMSB models are included.
#
#
#         PDG            Width
DECAY         6     1.41530907E+00   # top decays
#          BR         NDA      ID1       ID2
     1.00000000E+00    2           5        24   # BR(t ->  b    W+)
#
#         PDG            Width
DECAY   1000021     1.00044179E-02   # gluino decays
#          BR         NDA      ID1       ID2
     5.69845051E-06    2     1000022        21   # BR(~g -> ~chi_10 g)
     1.69089504E-03    2     1000023        21   # BR(~g -> ~chi_20 g)
#          BR         NDA      ID1       ID2
     5.50566814E-26    2     1000039        21   # BR(~g -> ~G      g)
#           BR         NDA      ID1       ID2       ID3
     7.85960381E-03    3     1000022         1        -1   # BR(~g -> ~chi_10 d  db)
     5.52957020E-02    3     1000023         1        -1   # BR(~g -> ~chi_20 d  db)
     5.18471631E-02    3     1000022         2        -2   # BR(~g -> ~chi_10 u  ub)
     6.03439505E-02    3     1000023         2        -2   # BR(~g -> ~chi_20 u  ub)
     7.85960381E-03    3     1000022         3        -3   # BR(~g -> ~chi_10 s  sb)
     5.52957020E-02    3     1000023         3        -3   # BR(~g -> ~chi_20 s  sb)
     5.18471631E-02    3     1000022         4        -4   # BR(~g -> ~chi_10 c  cb)
     6.03439505E-02    3     1000023         4        -4   # BR(~g -> ~chi_20 c  cb)
     1.11311548E-02    3     1000022         5        -5   # BR(~g -> ~chi_10 b  bb)
     9.94104870E-02    3     1000023         5        -5   # BR(~g -> ~chi_20 b  bb)
     2.23071982E-03    3     1000022         6        -6   # BR(~g -> ~chi_10 t  tb)
     1.15859462E-01    3     1000024         1        -2   # BR(~g -> ~chi_1+ d  ub)
     1.15859462E-01    3    -1000024         2        -1   # BR(~g -> ~chi_1- u  db)
     1.15859462E-01    3     1000024         3        -4   # BR(~g -> ~chi_1+ s  cb)
     1.15859462E-01    3    -1000024         4        -3   # BR(~g -> ~chi_1- c  sb)
     3.57001791E-02    3     1000024         5        -6   # BR(~g -> ~chi_1+ b  tb)
     3.57001791E-02    3    -1000024         6        -5   # BR(~g -> ~chi_1- t  bb)
#
#         PDG            Width
DECAY   1000006     7.29739564E-01   # stop1 decays
#          BR         NDA      ID1       ID2
     1.06317822E+00    2     1000022         6   # BR(~t_1 -> ~chi_10 t )
     1.17511527E-02    2     1000023         6   # BR(~t_1 -> ~chi_20 t )
    -7.49293735E-02    2     1000024         5   # BR(~t_1 -> ~chi_1+ b )
#
#         PDG            Width
DECAY   2000006     6.58374558E+00   # stop2 decays
#          BR         NDA      ID1       ID2
    -2.86938841E-03    2     1000022         6   # BR(~t_2 -> ~chi_10 t )
     2.79144509E-01    2     1000023         6   # BR(~t_2 -> ~chi_20 t )
     0.00000000E+00    2     1000025         6   # BR(~t_2 -> ~chi_30 t )
     6.57749897E-01    2     1000024         5   # BR(~t_2 -> ~chi_1+ b )
    -3.17807284E-02    2     1000037         5   # BR(~t_2 -> ~chi_2+ b )
     5.05252062E-02    2     1000006        25   # BR(~t_2 -> ~t_1    h )
     4.72305046E-02    2     1000006        23   # BR(~t_2 -> ~t_1    Z )
#
#         PDG            Width
DECAY   1000005     9.29698555E+00   # sbottom1 decays
#          BR         NDA      ID1       ID2
     1.25608850E-02    2     1000022         5   # BR(~b_1 -> ~chi_10 b )
     2.05736106E-01    2     1000023         5   # BR(~b_1 -> ~chi_20 b )
     1.50904334E-04    2     1000025         5   # BR(~b_1 -> ~chi_30 b )
     4.70504515E-04    2     1000035         5   # BR(~b_1 -> ~chi_40 b )
     3.44599859E-01    2    -1000024         6   # BR(~b_1 -> ~chi_1- t )
     4.09570486E-01    2     1000021         5   # BR(~b_1 -> ~g      b )
     2.69112549E-02    2     1000006       -24   # BR(~b_1 -> ~t_1    W-)
#
#         PDG            Width
DECAY   2000005     1.13547108E+01   # sbottom2 decays
#          BR         NDA      ID1       ID2
     3.42932376E-02    2     1000022         5   # BR(~b_2 -> ~chi_10 b )
     1.95073810E-04    2     1000023         5   # BR(~b_2 -> ~chi_20 b )
     6.24993273E-04    2     1000025         5   # BR(~b_2 -> ~chi_30 b )
     6.38942003E-04    2     1000035         5   # BR(~b_2 -> ~chi_40 b )
     6.89528659E-05    2    -1000024         6   # BR(~b_2 -> ~chi_1- t )
     5.11106845E-04    2    -1000037         6   # BR(~b_2 -> ~chi_2- t )
     9.61261241E-01    2     1000021         5   # BR(~b_2 -> ~g      b )
     1.01172034E-04    2     1000005        25   # BR(~b_2 -> ~b_1    h )
     8.23896779E-04    2     1000005        23   # BR(~b_2 -> ~b_1    Z )
     8.43158795E-04    2     1000006       -24   # BR(~b_2 -> ~t_1    W-)
     6.38224997E-04    2     2000006       -24   # BR(~b_2 -> ~t_2    W-)
#
#         PDG            Width
DECAY   1000002     1.28363029E+01   # sup_L decays
#          BR         NDA      ID1       ID2
     3.06791057E-03    2     1000022         2   # BR(~u_L -> ~chi_10 u)
     1.79910050E-01    2     1000023         2   # BR(~u_L -> ~chi_20 u)
     1.15684923E-04    2     1000025         2   # BR(~u_L -> ~chi_30 u)
     7.37521421E-04    2     1000035         2   # BR(~u_L -> ~chi_40 u)
     3.59674454E-01    2     1000024         1   # BR(~u_L -> ~chi_1+ d)
     6.83271578E-04    2     1000037         1   # BR(~u_L -> ~chi_2+ d)
     4.55811107E-01    2     1000021         2   # BR(~u_L -> ~g      u)
#
#         PDG            Width
DECAY   2000002     5.96899203E+00   # sup_R decays
#          BR         NDA      ID1       ID2
     2.01725288E-01    2     1000022         2   # BR(~u_R -> ~chi_10 u)
     9.06541179E-04    2     1000023         2   # BR(~u_R -> ~chi_20 u)
     4.94038544E-05    2     1000025         2   # BR(~u_R -> ~chi_30 u)
     1.64846106E-04    2     1000035         2   # BR(~u_R -> ~chi_40 u)
     7.97153921E-01    2     1000021         2   # BR(~u_R -> ~g      u)
#
#         PDG            Width
DECAY   1000001     1.29133800E+01   # sdown_L decays
#          BR         NDA      ID1       ID2
     1.01699732E-02    2     1000022         1   # BR(~d_L -> ~chi_10 d)
     1.71387948E-01    2     1000023         1   # BR(~d_L -> ~chi_20 d)
     1.88711046E-04    2     1000025         1   # BR(~d_L -> ~chi_30 d)
     1.07935567E-03    2     1000035         1   # BR(~d_L -> ~chi_40 d)
     3.45357018E-01    2    -1000024         2   # BR(~d_L -> ~chi_1- u)
     3.41679800E-03    2    -1000037         2   # BR(~d_L -> ~chi_2- u)
     4.68400196E-01    2     1000021         1   # BR(~d_L -> ~g      d)
#
#         PDG            Width
DECAY   2000001     1.38627378E+01   # sdown_R decays
#          BR         NDA      ID1       ID2
     2.82388820E-02    2     1000022         1   # BR(~d_R -> ~chi_10 d)
     1.29913250E-04    2     1000023         1   # BR(~d_R -> ~chi_20 d)
     1.87597877E-05    2     1000025         1   # BR(~d_R -> ~chi_30 d)
     6.80008589E-05    2     1000035         1   # BR(~d_R -> ~chi_40 d)
     9.71544444E-01    2     1000021         1   # BR(~d_R -> ~g      d)
#
#         PDG            Width
DECAY   1000004     1.28363029E+01   # scharm_L decays
#          BR         NDA      ID1       ID2
     3.06791057E-03    2     1000022         4   # BR(~c_L -> ~chi_10 c)
     1.79910050E-01    2     1000023         4   # BR(~c_L -> ~chi_20 c)
     1.15684923E-04    2     1000025         4   # BR(~c_L -> ~chi_30 c)
     7.37521421E-04    2     1000035         4   # BR(~c_L -> ~chi_40 c)
     3.59674454E-01    2     1000024         3   # BR(~c_L -> ~chi_1+ s)
     6.83271578E-04    2     1000037         3   # BR(~c_L -> ~chi_2+ s)
     4.55811107E-01    2     1000021         4   # BR(~c_L -> ~g      c)
#
#         PDG            Width
DECAY   2000004     5.96899203E+00   # scharm_R decays
#          BR         NDA      ID1       ID2
     2.01725288E-01    2     1000022         4   # BR(~c_R -> ~chi_10 c)
     9.06541179E-04    2     1000023         4   # BR(~c_R -> ~chi_20 c)
     4.94038544E-05    2     1000025         4   # BR(~c_R -> ~chi_30 c)
     1.64846106E-04    2     1000035         4   # BR(~c_R -> ~chi_40 c)
     7.97153921E-01    2     1000021         4   # BR(~c_R -> ~g      c)
#
#         PDG            Width
DECAY   1000003     1.29133800E+01   # sstrange_L decays
#          BR         NDA      ID1       ID2
     1.01699732E-02    2     1000022         3   # BR(~s_L -> ~chi_10 s)
     1.71387948E-01    2     1000023         3   # BR(~s_L -> ~chi_20 s)
     1.88711046E-04    2     1000025         3   # BR(~s_L -> ~chi_30 s)
     1.07935567E-03    2     1000035         3   # BR(~s_L -> ~chi_40 s)
     3.45357018E-01    2    -1000024         4   # BR(~s_L -> ~chi_1- c)
     3.41679800E-03    2    -1000037         4   # BR(~s_L -> ~chi_2- c)
     4.68400196E-01    2     1000021         3   # BR(~s_L -> ~g      s)
#
#         PDG            Width
DECAY   2000003     1.38627378E+01   # sstrange_R decays
#          BR         NDA      ID1       ID2
     2.82388820E-02    2     1000022         3   # BR(~s_R -> ~chi_10 s)
     1.29913250E-04    2     1000023         3   # BR(~s_R -> ~chi_20 s)
     1.87597877E-05    2     1000025         3   # BR(~s_R -> ~chi_30 s)
     6.80008589E-05    2     1000035         3   # BR(~s_R -> ~chi_40 s)
     9.71544444E-01    2     1000021         3   # BR(~s_R -> ~g      s)
#
#         PDG            Width
DECAY   1000011     8.56579317E+00   # selectron_L decays
#          BR         NDA      ID1       ID2
     7.72123205E-02    2     1000022        11   # BR(~e_L -> ~chi_10 e-)
     3.23898633E-01    2     1000023        11   # BR(~e_L -> ~chi_20 e-)
     1.27140591E-04    2     1000025        11   # BR(~e_L -> ~chi_30 e-)
     9.82553486E-04    2     1000035        11   # BR(~e_L -> ~chi_40 e-)
     5.91170838E-01    2    -1000024        12   # BR(~e_L -> ~chi_1- nu_e)
     6.60851368E-03    2    -1000037        12   # BR(~e_L -> ~chi_2- nu_e)
#
#         PDG            Width
DECAY   2000011     2.33589208E+00   # selectron_R decays
#          BR         NDA      ID1       ID2
     9.95765284E-01    2     1000022        11   # BR(~e_R -> ~chi_10 e-)
     4.23471600E-03    2     1000023        11   # BR(~e_R -> ~chi_20 e-)
#
#         PDG            Width
DECAY   1000013     8.56579317E+00   # smuon_L decays
#          BR         NDA      ID1       ID2
     7.72123205E-02    2     1000022        13   # BR(~mu_L -> ~chi_10 mu-)
     3.23898633E-01    2     1000023        13   # BR(~mu_L -> ~chi_20 mu-)
     1.27140591E-04    2     1000025        13   # BR(~mu_L -> ~chi_30 mu-)
     9.82553486E-04    2     1000035        13   # BR(~mu_L -> ~chi_40 mu-)
     5.91170838E-01    2    -1000024        14   # BR(~mu_L -> ~chi_1- nu_mu)
     6.60851368E-03    2    -1000037        14   # BR(~mu_L -> ~chi_2- nu_mu)
#
#         PDG            Width
DECAY   2000013     2.33589208E+00   # smuon_R decays
#          BR         NDA      ID1       ID2
     9.95765284E-01    2     1000022        13   # BR(~mu_R -> ~chi_10 mu-)
     4.23471600E-03    2     1000023        13   # BR(~mu_R -> ~chi_20 mu-)
#
#         PDG            Width
DECAY   1000015     2.33646660E+00   # stau_1 decays
#          BR         NDA      ID1       ID2
     9.92401527E-01    2     1000022        15   # BR(~tau_1 -> ~chi_10  tau-)
     5.34324676E-03    2     1000023        15   # BR(~tau_1 -> ~chi_20  tau-)
     2.25522593E-03    2    -1000024        16   # BR(~tau_1 -> ~chi_1-  nu_tau)
     3.52170807E-28    2     1000039        15   # BR(~tau_1 -> ~G       tau-)
#
#         PDG            Width
DECAY   2000015     8.56820670E+00   # stau_2 decays
#          BR         NDA      ID1       ID2
     7.72063210E-02    2     1000022        15   # BR(~tau_2 -> ~chi_10  tau-)
     3.23314865E-01    2     1000023        15   # BR(~tau_2 -> ~chi_20  tau-)
     6.29911734E-04    2     1000025        15   # BR(~tau_2 -> ~chi_30  tau-)
     1.46419917E-03    2     1000035        15   # BR(~tau_2 -> ~chi_40  tau-)
     5.89874792E-01    2    -1000024        16   # BR(~tau_2 -> ~chi_1-  nu_tau)
     6.67053739E-03    2    -1000037        16   # BR(~tau_2 -> ~chi_2-  nu_tau)
     4.72902442E-04    2     1000015        25   # BR(~tau_2 -> ~tau_1   h)
     3.66471590E-04    2     1000015        23   # BR(~tau_2 -> ~tau_1   Z)
#
#         PDG            Width
DECAY   1000012     8.62276394E+00   # snu_eL decays
#          BR         NDA      ID1       ID2
     1.12440266E-01    2     1000022        12   # BR(~nu_eL -> ~chi_10 nu_e)
     2.74865196E-01    2     1000023        12   # BR(~nu_eL -> ~chi_20 nu_e)
     4.99141288E-04    2     1000025        12   # BR(~nu_eL -> ~chi_30 nu_e)
     2.67199401E-03    2     1000035        12   # BR(~nu_eL -> ~chi_40 nu_e)
     6.08209745E-01    2     1000024        11   # BR(~nu_eL -> ~chi_1+ e-)
     1.31365748E-03    2     1000037        11   # BR(~nu_eL -> ~chi_2+ e-)
#
#         PDG            Width
DECAY   1000014     8.62276394E+00   # snu_muL decays
#          BR         NDA      ID1       ID2
     1.12440266E-01    2     1000022        14   # BR(~nu_muL -> ~chi_10 nu_mu)
     2.74865196E-01    2     1000023        14   # BR(~nu_muL -> ~chi_20 nu_mu)
     4.99141288E-04    2     1000025        14   # BR(~nu_muL -> ~chi_30 nu_mu)
     2.67199401E-03    2     1000035        14   # BR(~nu_muL -> ~chi_40 nu_mu)
     6.08209745E-01    2     1000024        13   # BR(~nu_muL -> ~chi_1+ mu-)
     1.31365748E-03    2     1000037        13   # BR(~nu_muL -> ~chi_2+ mu-)
#
#         PDG            Width
DECAY   1000016     8.63039365E+00   # snu_tauL decays
#          BR         NDA      ID1       ID2
     1.12252120E-01    2     1000022        16   # BR(~nu_tauL -> ~chi_10 nu_tau)
     2.74380018E-01    2     1000023        16   # BR(~nu_tauL -> ~chi_20 nu_tau)
     4.96460247E-04    2     1000025        16   # BR(~nu_tauL -> ~chi_30 nu_tau)
     2.65681495E-03    2     1000035        16   # BR(~nu_tauL -> ~chi_40 nu_tau)
     6.07330539E-01    2     1000024        15   # BR(~nu_tauL -> ~chi_1+ tau-)
     2.14930286E-03    2     1000037        15   # BR(~nu_tauL -> ~chi_2+ tau-)
     7.34745174E-04    2    -1000015       -24   # BR(~nu_tauL -> ~tau_1+ W-)
#
#         PDG            Width
DECAY   1000024     2.33788687E-05   # chargino1+ decays
#          BR         NDA      ID1       ID2
     8.78186524E-27    2     1000039        24   # BR(~chi_1+ -> ~G       W+)
#           BR         NDA      ID1       ID2       ID3
     3.36745774E-01    3     1000022         2        -1   # BR(~chi_1+ -> ~chi_10 u    db)
     3.36745774E-01    3     1000022         4        -3   # BR(~chi_1+ -> ~chi_10 c    sb)
     1.08951570E-01    3     1000022       -11        12   # BR(~chi_1+ -> ~chi_10 e+   nu_e)
     1.08951570E-01    3     1000022       -13        14   # BR(~chi_1+ -> ~chi_10 mu+  nu_mu)
     1.08605312E-01    3     1000022       -15        16   # BR(~chi_1+ -> ~chi_10 tau+ nu_tau)
#
#         PDG            Width
DECAY   1000037     4.00427892E+00   # chargino2+ decays
#          BR         NDA      ID1       ID2
     6.60733474E-02    2     1000006        -5   # BR(~chi_2+ -> ~t_1     bb)
     1.19598014E-05    2    -1000015        16   # BR(~chi_2+ -> ~tau_1+  nu_tau)
     2.85445081E-01    2     1000024        23   # BR(~chi_2+ -> ~chi_1+  Z )
     7.15842014E-02    2     1000022        24   # BR(~chi_2+ -> ~chi_10  W+)
     3.04926257E-01    2     1000023        24   # BR(~chi_2+ -> ~chi_20  W+)
     2.71959154E-01    2     1000024        25   # BR(~chi_2+ -> ~chi_1+  h )
     1.14330754E-28    2     1000039        24   # BR(~chi_2+ -> ~G       W+)
#
#         PDG            Width
DECAY   1000022     3.15278859E-32   # neutralino1 decays
#          BR         NDA      ID1       ID2
     1.00000000E+00    2     1000039        22   # BR(~chi_10 -> ~G        gam)
#
#         PDG            Width
DECAY   1000023     6.67892093E-07   # neutralino2 decays
#          BR         NDA      ID1       ID2
     2.85538261E-08    2     1000022        22   # BR(~chi_20 -> ~chi_10 gam)
#          BR         NDA      ID1       ID2
     4.94237916E-25    2     1000039        22   # BR(~chi_20 -> ~G        gam)
     1.49843511E-25    2     1000039        23   # BR(~chi_20 -> ~G        Z)
     6.54593174E-29    2     1000039        25   # BR(~chi_20 -> ~G        h)
#           BR         NDA      ID1       ID2       ID3
     1.39706837E-01    3     1000022        -2         2   # BR(~chi_20 -> ~chi_10 ub      u)
     1.99618576E-01    3     1000022        -1         1   # BR(~chi_20 -> ~chi_10 db      d)
     1.39706837E-01    3     1000022        -4         4   # BR(~chi_20 -> ~chi_10 cb      c)
     1.99618576E-01    3     1000022        -3         3   # BR(~chi_20 -> ~chi_10 sb      s)
     2.12995855E-01    3     1000022        -5         5   # BR(~chi_20 -> ~chi_10 bb      b)
     1.10071868E-02    3     1000022       -11        11   # BR(~chi_20 -> ~chi_10 e+      e-)
     1.10071868E-02    3     1000022       -13        13   # BR(~chi_20 -> ~chi_10 mu+     mu-)
     1.20758093E-02    3     1000022       -15        15   # BR(~chi_20 -> ~chi_10 tau+    tau-)
     2.47691021E-02    3     1000022       -12        12   # BR(~chi_20 -> ~chi_10 nu_eb   nu_e)
     2.47691021E-02    3     1000022       -14        14   # BR(~chi_20 -> ~chi_10 nu_mub  nu_mu)
     2.47249029E-02    3     1000022       -16        16   # BR(~chi_20 -> ~chi_10 nu_taub nu_tau)
     2.83805103E-12    3     1000024        -2         1   # BR(~chi_20 -> ~chi_1+ ub      d)
     2.83805103E-12    3    -1000024        -1         2   # BR(~chi_20 -> ~chi_1- db      u)
     2.83805103E-12    3     1000024        -4         3   # BR(~chi_20 -> ~chi_1+ cb      s)
     2.83805103E-12    3    -1000024        -3         4   # BR(~chi_20 -> ~chi_1- sb      c)
     9.44292154E-13    3     1000024       -12        11   # BR(~chi_20 -> ~chi_1+ nu_eb   e-)
     9.44292154E-13    3    -1000024        12       -11   # BR(~chi_20 -> ~chi_1- nu_e    e+)
     9.44292154E-13    3     1000024       -14        13   # BR(~chi_20 -> ~chi_1+ nu_mub  mu-)
     9.44292154E-13    3    -1000024        14       -13   # BR(~chi_20 -> ~chi_1- nu_mu   mu+)
#
#         PDG            Width
DECAY   1000025     3.78601318E+00   # neutralino3 decays
#          BR         NDA      ID1       ID2
     8.95209991E-02    2     1000022        23   # BR(~chi_30 -> ~chi_10   Z )
     2.42846366E-01    2     1000023        23   # BR(~chi_30 -> ~chi_20   Z )
     3.06413739E-01    2     1000024       -24   # BR(~chi_30 -> ~chi_1+   W-)
     3.06413739E-01    2    -1000024        24   # BR(~chi_30 -> ~chi_1-   W+)
     2.00056197E-02    2     1000022        25   # BR(~chi_30 -> ~chi_10   h )
     3.47942435E-02    2     1000023        25   # BR(~chi_30 -> ~chi_20   h )
     2.55665052E-07    2     2000011       -11   # BR(~chi_30 -> ~e_R-     e+)
     2.55665052E-07    2    -2000011        11   # BR(~chi_30 -> ~e_R+     e-)
     2.55665052E-07    2     2000013       -13   # BR(~chi_30 -> ~mu_R-    mu+)
     2.55665052E-07    2    -2000013        13   # BR(~chi_30 -> ~mu_R+    mu-)
     2.13582324E-06    2     1000015       -15   # BR(~chi_30 -> ~tau_1-   tau+)
     2.13582324E-06    2    -1000015        15   # BR(~chi_30 -> ~tau_1+   tau-)
     6.08404279E-33    2     1000039        22   # BR(~chi_30 -> ~G        gam)
     3.48362822E-29    2     1000039        23   # BR(~chi_30 -> ~G        Z)
     7.06735155E-29    2     1000039        25   # BR(~chi_30 -> ~G        h)
#
#         PDG            Width
DECAY   1000035     3.83040990E+00   # neutralino4 decays
#          BR         NDA      ID1       ID2
     2.29190758E-02    2     1000022        23   # BR(~chi_40 -> ~chi_10   Z )
     3.85524640E-02    2     1000023        23   # BR(~chi_40 -> ~chi_20   Z )
     3.17008636E-01    2     1000024       -24   # BR(~chi_40 -> ~chi_1+   W-)
     3.17008636E-01    2    -1000024        24   # BR(~chi_40 -> ~chi_1-   W+)
     7.60563440E-02    2     1000022        25   # BR(~chi_40 -> ~chi_10   h )
     2.28404702E-01    2     1000023        25   # BR(~chi_40 -> ~chi_20   h )
     6.30299072E-06    2     2000011       -11   # BR(~chi_40 -> ~e_R-     e+)
     6.30299072E-06    2    -2000011        11   # BR(~chi_40 -> ~e_R+     e-)
     6.30299072E-06    2     2000013       -13   # BR(~chi_40 -> ~mu_R-    mu+)
     6.30299072E-06    2    -2000013        13   # BR(~chi_40 -> ~mu_R+    mu-)
     1.24648204E-05    2     1000015       -15   # BR(~chi_40 -> ~tau_1-   tau+)
     1.24648204E-05    2    -1000015        15   # BR(~chi_40 -> ~tau_1+   tau-)
     2.23112901E-32    2     1000039        22   # BR(~chi_40 -> ~G        gam)
     8.18880623E-29    2     1000039        23   # BR(~chi_40 -> ~G        Z)
     3.20500028E-29    2     1000039        25   # BR(~chi_40 -> ~G        h)
#
#         PDG            Width
DECAY        25     3.33954751E-03   # h decays
#          BR         NDA      ID1       ID2
     8.34667112E-01    2           5        -5   # BR(h -> b       bb     )
     8.25965544E-02    2         -15        15   # BR(h -> tau+    tau-   )
     2.92542576E-04    2         -13        13   # BR(h -> mu+     mu-    )
     6.29569749E-04    2           3        -3   # BR(h -> s       sb     )
     2.08249646E-02    2           4        -4   # BR(h -> c       cb     )
     4.76829733E-02    2          21        21   # BR(h -> g       g      )
     1.31091654E-03    2          22        22   # BR(h -> gam     gam    )
     6.56524716E-05    2          22        23   # BR(h -> Z       gam    )
     1.08749951E-02    2          24       -24   # BR(h -> W+      W-     )
     1.05471899E-03    2          23        23   # BR(h -> Z       Z      )
#
#         PDG            Width
DECAY        35     4.36424999E+00   # H decays
#          BR         NDA      ID1       ID2
     7.01440631E-02    2           5        -5   # BR(H -> b       bb     )
     1.04575824E-02    2         -15        15   # BR(H -> tau+    tau-   )
     3.69729066E-05    2         -13        13   # BR(H -> mu+     mu-    )
     5.39852937E-05    2           3        -3   # BR(H -> s       sb     )
     3.69799022E-06    2           4        -4   # BR(H -> c       cb     )
     3.03867421E-01    2           6        -6   # BR(H -> t       tb     )
     6.08706801E-04    2          21        21   # BR(H -> g       g      )
     3.51939482E-06    2          22        22   # BR(H -> gam     gam    )
     2.04337170E-07    2          23        22   # BR(H -> Z       gam    )
     2.15017419E-02    2          24       -24   # BR(H -> W+      W-     )
     1.05464573E-02    2          23        23   # BR(H -> Z       Z      )
     5.21590603E-03    2          25        25   # BR(H -> h       h      )
     1.75628079E-22    2          36        36   # BR(H -> A       A      )
     1.91193787E-14    2          23        36   # BR(H -> Z       A      )
     4.70332692E-02    2     1000024  -1000024   # BR(H -> ~chi_1+ ~chi_1-)
     1.46390795E-01    2     1000024  -1000037   # BR(H -> ~chi_1+ ~chi_2-)
     1.46390795E-01    2     1000037  -1000024   # BR(H -> ~chi_2+ ~chi_1-)
     3.45394436E-03    2     1000022   1000022   # BR(H -> ~chi_10 ~chi_10)
     2.05450215E-02    2     1000023   1000023   # BR(H -> ~chi_20 ~chi_20)
     1.70481980E-02    2     1000022   1000023   # BR(H -> ~chi_10 ~chi_20)
     5.02912149E-02    2     1000022   1000025   # BR(H -> ~chi_10 ~chi_30)
     9.53431847E-03    2     1000022   1000035   # BR(H -> ~chi_10 ~chi_40)
     1.26080399E-01    2     1000023   1000025   # BR(H -> ~chi_20 ~chi_30)
     1.07917855E-02    2     1000023   1000035   # BR(H -> ~chi_20 ~chi_40)
#
#         PDG            Width
DECAY        36     4.18659669E+00   # A decays
#          BR         NDA      ID1       ID2
     7.39351413E-02    2           5        -5   # BR(A -> b       bb     )
     1.10082123E-02    2         -15        15   # BR(A -> tau+    tau-   )
     3.89188117E-05    2         -13        13   # BR(A -> mu+     mu-    )
     5.69693853E-05    2           3        -3   # BR(A -> s       sb     )
     2.97257115E-06    2           4        -4   # BR(A -> c       cb     )
     2.91137448E-01    2           6        -6   # BR(A -> t       tb     )
     7.22501054E-04    2          21        21   # BR(A -> g       g      )
     2.74013896E-07    2          22        22   # BR(A -> gam     gam    )
     7.14895327E-07    2          23        22   # BR(A -> Z       gam    )
     2.16216877E-02    2          23        25   # BR(A -> Z       h      )
     8.23703623E-02    2     1000024  -1000024   # BR(A -> ~chi_1+ ~chi_1-)
     1.37094594E-01    2     1000024  -1000037   # BR(A -> ~chi_1+ ~chi_2-)
     1.37094594E-01    2     1000037  -1000024   # BR(A -> ~chi_2+ ~chi_1-)
     4.93685140E-03    2     1000022   1000022   # BR(A -> ~chi_10 ~chi_10)
     3.59938221E-02    2     1000023   1000023   # BR(A -> ~chi_20 ~chi_20)
     2.64857131E-02    2     1000022   1000023   # BR(A -> ~chi_10 ~chi_20)
     1.31193609E-02    2     1000022   1000025   # BR(A -> ~chi_10 ~chi_30)
     4.30592271E-02    2     1000022   1000035   # BR(A -> ~chi_10 ~chi_40)
     1.60306412E-02    2     1000023   1000025   # BR(A -> ~chi_20 ~chi_30)
     1.05289993E-01    2     1000023   1000035   # BR(A -> ~chi_20 ~chi_40)
#
#         PDG            Width
DECAY        37     3.96408000E+00   # H+ decays
#          BR         NDA      ID1       ID2
     1.17594433E-04    2           4        -5   # BR(H+ -> c       bb     )
     1.16930900E-02    2         -15        16   # BR(H+ -> tau+    nu_tau )
     4.13401467E-05    2         -13        14   # BR(H+ -> mu+     nu_mu  )
     7.52561561E-07    2           2        -5   # BR(H+ -> u       bb     )
     2.87483170E-06    2           2        -3   # BR(H+ -> u       sb     )
     6.21401644E-05    2           4        -3   # BR(H+ -> c       sb     )
     3.68272128E-01    2           6        -5   # BR(H+ -> t       bb     )
     2.34982596E-02    2          24        25   # BR(H+ -> W+      h      )
     5.17846056E-10    2          24        36   # BR(H+ -> W+      A      )
     4.13647564E-02    2     1000024   1000022   # BR(H+ -> ~chi_1+ ~chi_10)
     1.37289300E-04    2     1000024   1000023   # BR(H+ -> ~chi_1+ ~chi_20)
     1.74841951E-01    2     1000024   1000025   # BR(H+ -> ~chi_1+ ~chi_30)
     1.60341217E-01    2     1000024   1000035   # BR(H+ -> ~chi_1+ ~chi_40)
     3.83787162E-02    2     1000037   1000022   # BR(H+ -> ~chi_2+ ~chi_10)
     1.81247890E-01    2     1000037   1000023   # BR(H+ -> ~chi_2+ ~chi_20)
#
#         PDG            Width
DECAY   1000039     0.00000000E+00   # gravitino decays
#          BR         NDA      ID1       ID2
