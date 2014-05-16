# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-3.62264486e+03
# dtanbeta/dmu=3.07693317e-02
# mu=8.27073698e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=4.43772897e-01
# dtanbeta/dm3sq=-9.37942045e-05
# m3sq= 8.79187592e+04
# SOFTSUSY3.4.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.0       # version number
Block MODSEL  # Select model
     1    1   # sugra
Block SMINPUTS             # Standard Model inputs
     1    1.27934000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.17200000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.25000000e+00   # mb(mb)
     6    1.73300000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    1.00000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    3.00000000e+02   # m0
     2    7.00000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.71152703e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.84947969e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03931127e+01   # MW
        25     1.17067670e+02   # h0
        35     1.00595670e+03   # H0
        36     1.00574925e+03   # A0
        37     1.00918465e+03   # H+
   1000021     1.56851468e+03   # ~g
   1000022     2.91962586e+02   # ~neutralino(1)
   1000023     5.52362296e+02   # ~neutralino(2)
   1000024     5.52454302e+02   # ~chargino(1)
   1000025    -8.53578078e+02   # ~neutralino(3)
   1000035     8.65066744e+02   # ~neutralino(4)
   1000037     8.65257047e+02   # ~chargino(2)
   1000001     1.45585623e+03   # ~d_L
   1000002     1.45384500e+03   # ~u_L
   1000003     1.45585264e+03   # ~s_L
   1000004     1.45384141e+03   # ~c_L
   1000005     1.33137028e+03   # ~b_1
   1000006     1.12000494e+03   # ~t_1
   1000011     5.57333133e+02   # ~e_L
   1000012     5.51506387e+02   # ~nue_L
   1000013     5.57327655e+02   # ~mu_L
   1000014     5.51500854e+02   # ~numu_L
   1000015     3.93492192e+02   # ~stau_1
   1000016     5.49632034e+02   # ~nu_tau_L
   2000001     1.39391535e+03   # ~d_R
   2000002     1.39934658e+03   # ~u_R
   2000003     1.39391163e+03   # ~s_R
   2000004     1.39934272e+03   # ~c_R
   2000005     1.38791967e+03   # ~b_2
   2000006     1.36509970e+03   # ~t_2
   2000011     4.00523820e+02   # ~e_R
   2000013     4.00508353e+02   # ~mu_R
   2000015     5.56881937e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05391045e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97667016e-01   # N_{1,1}
  1  2    -9.68350759e-03   # N_{1,2}
  1  3     6.19916779e-02   # N_{1,3}
  1  4    -2.69032971e-02   # N_{1,4}
  2  1     2.26501091e-02   # N_{2,1}
  2  2     9.79944386e-01   # N_{2,2}
  2  3    -1.62149718e-01   # N_{2,3}
  2  4     1.13593315e-01   # N_{2,4}
  3  1    -2.43873403e-02   # N_{3,1}
  3  2     3.51979248e-02   # N_{3,2}
  3  3     7.05139157e-01   # N_{3,3}
  3  4     7.07774775e-01   # N_{3,4}
  4  1    -5.96050000e-02   # N_{4,1}
  4  2     1.95898789e-01   # N_{4,2}
  4  3     6.87490560e-01   # N_{4,3}
  4  4    -6.96726373e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.72851048e-01   # U_{1,1}
  1  2    -2.31432147e-01   # U_{1,2}
  2  1     2.31432147e-01   # U_{2,1}
  2  2     9.72851048e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.86688380e-01   # V_{1,1}
  1  2    -1.62622390e-01   # V_{1,2}
  2  1     1.62622390e-01   # V_{2,1}
  2  2     9.86688380e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.39660658e-01   # F_{11}
  1  2     9.40548052e-01   # F_{12}
  2  1     9.40548052e-01   # F_{21}
  2  2    -3.39660658e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.88397609e-01   # F_{11}
  1  2     1.51888667e-01   # F_{12}
  2  1    -1.51888667e-01   # F_{21}
  2  2     9.88397609e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.00339207e-01   # F_{11}
  1  2     9.94953287e-01   # F_{12}
  2  1     9.94953287e-01   # F_{21}
  2  2    -1.00339207e-01   # F_{22}
Block gauge Q= 1.20066566e+03  # SM gauge couplings
     1     3.63001686e-01   # g'(Q)MSSM DRbar
     2     6.41211676e-01   # g(Q)MSSM DRbar
     3     1.04447278e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.20066566e+03  
  3  3     8.50555260e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.20066566e+03  
  3  3     1.33207708e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.20066566e+03  
  3  3     1.00279252e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.20066566e+03 # Higgs mixing parameters
     1     8.47927270e+02    # mu(Q)MSSM DRbar
     2     9.63615566e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43735195e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.04646571e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.20066566e+03  # MSSM DRbar SUSY breaking parameters
     1     2.97414568e+02      # M_1(Q)
     2     5.47621040e+02      # M_2(Q)
     3     1.52241159e+03      # M_3(Q)
    21     2.70177392e+05      # mH1^2(Q)
    22    -6.89511790e+05      # mH2^2(Q)
    31     5.49280258e+02      # meL(Q)
    32     5.49274690e+02      # mmuL(Q)
    33     5.47588989e+02      # mtauL(Q)
    34     3.94207648e+02      # meR(Q)
    35     3.94191922e+02      # mmuR(Q)
    36     3.89408237e+02      # mtauR(Q)
    41     1.40738376e+03      # mqL1(Q)
    42     1.40738010e+03      # mqL2(Q)
    43     1.29431157e+03      # mqL3(Q)
    44     1.35417944e+03      # muR(Q)
    45     1.35417551e+03      # mcR(Q)
    46     1.10797448e+03      # mtR(Q)
    47     1.34764615e+03      # mdR(Q)
    48     1.34764236e+03      # msR(Q)
    49     1.34074014e+03      # mbR(Q)
Block au Q= 1.20066566e+03  
  1  1    -1.54129775e+03      # Au(Q)MSSM DRbar
  2  2    -1.54129093e+03      # Ac(Q)MSSM DRbar
  3  3    -1.19522592e+03      # At(Q)MSSM DRbar
Block ad Q= 1.20066566e+03  
  1  1    -1.87853387e+03      # Ad(Q)MSSM DRbar
  2  2    -1.87852757e+03      # As(Q)MSSM DRbar
  3  3    -1.75708155e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.20066566e+03  
  1  1    -4.12752017e+02      # Ae(Q)MSSM DRbar
  2  2    -4.12744750e+02      # Amu(Q)MSSM DRbar
  3  3    -4.10545826e+02      # Atau(Q)MSSM DRbar
