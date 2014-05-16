# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.96083937e+03
# dtanbeta/dmu=-6.44328546e-01
# mu=8.68405915e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=9.69289142e-02
# dtanbeta/dm3sq=-2.21759538e-03
# m3sq= -3.06511887e+05
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
     3    4.00000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    3.60000000e+02   # m0
     2    6.00000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.76224225e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=6.82034205e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03891580e+01   # MW
        25     1.17909216e+02   # h0
        35     7.09001926e+02   # H0
        36     7.09009018e+02   # A0
        37     7.13866299e+02   # H+
   1000021     1.36311077e+03   # ~g
   1000022     2.50491293e+02   # ~neutralino(1)
   1000023     4.76950548e+02   # ~neutralino(2)
   1000024     4.77088166e+02   # ~chargino(1)
   1000025    -8.08654248e+02   # ~neutralino(3)
   1000035     8.16410178e+02   # ~neutralino(4)
   1000037     8.17389375e+02   # ~chargino(2)
   1000001     1.28924877e+03   # ~d_L
   1000002     1.28691129e+03   # ~u_L
   1000003     1.28921228e+03   # ~s_L
   1000004     1.28687473e+03   # ~c_L
   1000005     1.08432310e+03   # ~b_1
   1000006     9.34125660e+02   # ~t_1
   1000011     5.41168328e+02   # ~e_L
   1000012     5.35055678e+02   # ~nue_L
   1000013     5.41029421e+02   # ~mu_L
   1000014     5.34915296e+02   # ~numu_L
   1000015     2.51485757e+02   # ~stau_1
   1000016     4.86683485e+02   # ~nu_tau_L
   2000001     1.23763325e+03   # ~d_R
   2000002     1.24187350e+03   # ~u_R
   2000003     1.23756070e+03   # ~s_R
   2000004     1.24186891e+03   # ~c_R
   2000005     1.15070806e+03   # ~b_2
   2000006     1.15361505e+03   # ~t_2
   2000011     4.26418621e+02   # ~e_R
   2000013     4.26060743e+02   # ~mu_R
   2000015     5.13689748e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.60393024e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97789707e-01   # N_{1,1}
  1  2    -7.78812782e-03   # N_{1,2}
  1  3     6.25039270e-02   # N_{1,3}
  1  4    -2.11731961e-02   # N_{1,4}
  2  1     1.88492888e-02   # N_{2,1}
  2  2     9.84770812e-01   # N_{2,2}
  2  3    -1.47913541e-01   # N_{2,3}
  2  4     8.94021007e-02   # N_{2,4}
  3  1    -2.87906389e-02   # N_{3,1}
  3  2     4.21298172e-02   # N_{3,2}
  3  3     7.04606597e-01   # N_{3,3}
  3  4     7.07761063e-01   # N_{3,4}
  4  1    -5.68463146e-02   # N_{4,1}
  4  2     1.68495907e-01   # N_{4,2}
  4  3     6.91190558e-01   # N_{4,3}
  4  4    -7.00452167e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.77585652e-01   # U_{1,1}
  1  2    -2.10538102e-01   # U_{1,2}
  2  1     2.10538102e-01   # U_{2,1}
  2  2     9.77585652e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.91813916e-01   # V_{1,1}
  1  2    -1.27691644e-01   # V_{1,2}
  2  1     1.27691644e-01   # V_{2,1}
  2  2     9.91813916e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.35434852e-01   # F_{11}
  1  2     9.00220245e-01   # F_{12}
  2  1     9.00220245e-01   # F_{21}
  2  2    -4.35434852e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     8.10552567e-01   # F_{11}
  1  2     5.85665891e-01   # F_{12}
  2  1    -5.85665891e-01   # F_{21}
  2  2     8.10552567e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     3.19331280e-01   # F_{11}
  1  2     9.47643147e-01   # F_{12}
  2  1     9.47643147e-01   # F_{21}
  2  2    -3.19331280e-01   # F_{22}
Block gauge Q= 1.00889492e+03  # SM gauge couplings
     1     3.62517629e-01   # g'(Q)MSSM DRbar
     2     6.41504716e-01   # g(Q)MSSM DRbar
     3     1.05239343e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.00889492e+03  
  3  3     8.50573351e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.00889492e+03  
  3  3     4.90337555e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.00889492e+03  
  3  3     4.27722628e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.00889492e+03 # Higgs mixing parameters
     1     8.05566086e+02    # mu(Q)MSSM DRbar
     2     3.91934134e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43813920e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.36197984e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.00889492e+03  # MSSM DRbar SUSY breaking parameters
     1     2.54466508e+02      # M_1(Q)
     2     4.70473169e+02      # M_2(Q)
     3     1.32360810e+03      # M_3(Q)
    21    -1.47189824e+05      # mH1^2(Q)
    22    -6.37527768e+05      # mH2^2(Q)
    31     5.34816746e+02      # meL(Q)
    32     5.34675796e+02      # mmuL(Q)
    33     4.89462790e+02      # mtauL(Q)
    34     4.21534329e+02      # meR(Q)
    35     4.21172295e+02      # mmuR(Q)
    36     2.88577436e+02      # mtauR(Q)
    41     1.24858790e+03      # mqL1(Q)
    42     1.24855069e+03      # mqL2(Q)
    43     1.07783670e+03      # mqL3(Q)
    44     1.20437343e+03      # muR(Q)
    45     1.20436876e+03      # mcR(Q)
    46     9.40037320e+02      # mtR(Q)
    47     1.19903075e+03      # mdR(Q)
    48     1.19895700e+03      # msR(Q)
    49     1.09398872e+03      # mbR(Q)
Block au Q= 1.00889492e+03  
  1  1    -1.68706383e+03      # Au(Q)MSSM DRbar
  2  2    -1.68702371e+03      # Ac(Q)MSSM DRbar
  3  3    -1.17940886e+03      # At(Q)MSSM DRbar
Block ad Q= 1.00889492e+03  
  1  1    -1.95299970e+03      # Ad(Q)MSSM DRbar
  2  2    -1.95289754e+03      # As(Q)MSSM DRbar
  3  3    -1.64811980e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.00889492e+03  
  1  1    -6.66000592e+02      # Ae(Q)MSSM DRbar
  2  2    -6.65673910e+02      # Amu(Q)MSSM DRbar
  3  3    -5.56205594e+02      # Atau(Q)MSSM DRbar
