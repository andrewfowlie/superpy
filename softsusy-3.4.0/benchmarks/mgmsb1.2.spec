# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-1.22046100e+03
# dtanbeta/dmu=1.09184637e-01
# mu=3.14323732e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=2.85428995e-01
# dtanbeta/dm3sq=-1.22856144e-03
# m3sq= 1.34961341e+04
# SOFTSUSY3.4.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.0       # version number
Block MODSEL  # Select model
     1    2   # gmsb
Block SMINPUTS             # Standard Model inputs
     1    1.27934000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.17200000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.25000000e+00   # mb(mb)
     6    1.73300000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    1.50000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    4.00000000e+04   # lambda
     2    8.00000000e+04   # M_mess
     5    3.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.95939112e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04057403e+01   # MW
        25     1.11154896e+02   # h0
        35     3.86024389e+02   # H0
        36     3.85674883e+02   # A0
        37     3.94229006e+02   # H+
   1000021     9.46337205e+02   # ~g
   1000022     1.62902694e+02   # ~neutralino(1)
   1000023     2.68746924e+02   # ~neutralino(2)
   1000024     2.65481975e+02   # ~chargino(1)
   1000025    -3.15506724e+02   # ~neutralino(3)
   1000035     3.86477867e+02   # ~neutralino(4)
   1000037     3.86223841e+02   # ~chargino(2)
   1000039     7.52982857e-10   # ~gravitino
   1000001     9.02056618e+02   # ~d_L
   1000002     8.98670851e+02   # ~u_L
   1000003     9.02055482e+02   # ~s_L
   1000004     8.98669710e+02   # ~c_L
   1000005     8.58348068e+02   # ~b_1
   1000006     8.05237820e+02   # ~t_1
   1000011     2.62668404e+02   # ~e_L
   1000012     2.50063770e+02   # ~nue_L
   1000013     2.62667360e+02   # ~mu_L
   1000014     2.50062675e+02   # ~numu_L
   1000015     1.24150230e+02   # ~stau_1
   1000016     2.49563487e+02   # ~nu_tau_L
   2000001     8.66763157e+02   # ~d_R
   2000002     8.68146070e+02   # ~u_R
   2000003     8.66761569e+02   # ~s_R
   2000004     8.68145263e+02   # ~c_R
   2000005     8.73282735e+02   # ~b_2
   2000006     8.90069042e+02   # ~t_2
   2000011     1.30907873e+02   # ~e_R
   2000013     1.30903621e+02   # ~mu_R
   2000015     2.64600197e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.89224382e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.64851899e-01   # N_{1,1}
  1  2    -6.51685557e-02   # N_{1,2}
  1  3     2.20673704e-01   # N_{1,3}
  1  4    -1.26952699e-01   # N_{1,4}
  2  1    -2.33230555e-01   # N_{2,1}
  2  2    -6.35464531e-01   # N_{2,2}
  2  3     5.51919571e-01   # N_{2,3}
  2  4    -4.87004236e-01   # N_{2,4}
  3  1    -6.00513212e-02   # N_{3,1}
  3  2     8.12339580e-02   # N_{3,2}
  3  3     6.95699383e-01   # N_{3,3}
  3  4     7.11194243e-01   # N_{3,4}
  4  1    -1.05157782e-01   # N_{4,1}
  4  2     7.65074463e-01   # N_{4,2}
  4  3     4.03348821e-01   # N_{4,3}
  4  4    -4.90828519e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1    -5.68950014e-01   # U_{1,1}
  1  2     8.22372106e-01   # U_{1,2}
  2  1    -8.22372106e-01   # U_{2,1}
  2  2    -5.68950014e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1    -6.93734183e-01   # V_{1,1}
  1  2     7.20231132e-01   # V_{1,2}
  2  1    -7.20231132e-01   # V_{2,1}
  2  2    -6.93734183e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.62129331e-01   # F_{11}
  1  2     9.32127860e-01   # F_{12}
  2  1     9.32127860e-01   # F_{21}
  2  2    -3.62129331e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     5.69611424e-01   # F_{11}
  1  2     8.21914123e-01   # F_{12}
  2  1     8.21914123e-01   # F_{21}
  2  2    -5.69611424e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.52521912e-01   # F_{11}
  1  2     9.88300089e-01   # F_{12}
  2  1     9.88300089e-01   # F_{21}
  2  2    -1.52521912e-01   # F_{22}
Block gauge Q= 8.20582347e+02  # SM gauge couplings
     1     3.62662309e-01   # g'(Q)MSSM DRbar
     2     6.45029862e-01   # g(Q)MSSM DRbar
     3     1.06791209e+00   # g3(Q)MSSM DRbar
Block yu Q= 8.20582347e+02  
  3  3     8.68563997e-01   # Yt(Q)MSSM DRbar
Block yd Q= 8.20582347e+02  
  3  3     2.06184399e-01   # Yb(Q)MSSM DRbar
Block ye Q= 8.20582347e+02  
  3  3     1.52072546e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 8.20582347e+02 # Higgs mixing parameters
     1     3.06454260e+02    # mu(Q)MSSM DRbar
     2     1.45355741e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43925178e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.70065023e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 8.20582347e+02  # MSSM DRbar SUSY breaking parameters
     1     1.73047348e+02      # M_1(Q)
     2     3.27235279e+02      # M_2(Q)
     3     9.01438944e+02      # M_3(Q)
    21     5.57284369e+04      # mH1^2(Q)
    22    -8.42643334e+04      # mH2^2(Q)
    31     2.53786343e+02      # meL(Q)
    32     2.53785266e+02      # mmuL(Q)
    33     2.53454659e+02      # mtauL(Q)
    34     1.18735356e+02      # meR(Q)
    35     1.18730678e+02      # mmuR(Q)
    36     1.17286483e+02      # mtauR(Q)
    41     8.66831111e+02      # mqL1(Q)
    42     8.66829939e+02      # mqL2(Q)
    43     8.39044933e+02      # mqL3(Q)
    44     8.36112161e+02      # muR(Q)
    45     8.36111335e+02      # mcR(Q)
    46     7.79721678e+02      # mtR(Q)
    47     8.33361830e+02      # mdR(Q)
    48     8.33360200e+02      # msR(Q)
    49     8.30180749e+02      # mbR(Q)
Block au Q= 8.20582347e+02  
  1  1    -2.73808201e+02      # Au(Q)MSSM DRbar
  2  2    -2.73807813e+02      # Ac(Q)MSSM DRbar
  3  3    -2.57999753e+02      # At(Q)MSSM DRbar
Block ad Q= 8.20582347e+02  
  1  1    -2.91459634e+02      # Ad(Q)MSSM DRbar
  2  2    -2.91459091e+02      # As(Q)MSSM DRbar
  3  3    -2.85497897e+02      # Ab(Q)MSSM DRbar
Block ae Q= 8.20582347e+02  
  1  1    -2.73064985e+01      # Ae(Q)MSSM DRbar
  2  2    -2.73062992e+01      # Amu(Q)MSSM DRbar
  3  3    -2.72451564e+01      # Atau(Q)MSSM DRbar
