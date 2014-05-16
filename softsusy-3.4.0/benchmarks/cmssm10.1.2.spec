# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.92260744e+03
# dtanbeta/dmu=3.96548692e-02
# mu=6.69085543e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=4.44456467e-01
# dtanbeta/dm3sq=-1.54914238e-04
# m3sq= 4.69033460e+04
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
     1    1.37500000e+02   # m0
     2    5.50000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.81512818e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=8.79897917e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03955000e+01   # MW
        25     1.15524525e+02   # h0
        35     7.83714625e+02   # H0
        36     7.83462271e+02   # A0
        37     7.87806000e+02   # H+
   1000021     1.25180884e+03   # ~g
   1000022     2.25931445e+02   # ~neutralino(1)
   1000023     4.27307774e+02   # ~neutralino(2)
   1000024     4.27348721e+02   # ~chargino(1)
   1000025    -6.90897468e+02   # ~neutralino(3)
   1000035     7.03889518e+02   # ~neutralino(4)
   1000037     7.04170233e+02   # ~chargino(2)
   1000001     1.15134188e+03   # ~d_L
   1000002     1.14875403e+03   # ~u_L
   1000003     1.15133908e+03   # ~s_L
   1000004     1.14875122e+03   # ~c_L
   1000005     1.05451495e+03   # ~b_1
   1000006     8.81640184e+02   # ~t_1
   1000011     3.96770166e+02   # ~e_L
   1000012     3.88649496e+02   # ~nue_L
   1000013     3.96766595e+02   # ~mu_L
   1000014     3.88645853e+02   # ~numu_L
   1000015     2.44850679e+02   # ~stau_1
   1000016     3.87373090e+02   # ~nu_tau_L
   2000001     1.10217880e+03   # ~d_R
   2000002     1.10605651e+03   # ~u_R
   2000003     1.10217587e+03   # ~s_R
   2000004     1.10605352e+03   # ~c_R
   2000005     1.09829043e+03   # ~b_2
   2000006     1.09684400e+03   # ~t_2
   2000011     2.51946166e+02   # ~e_R
   2000013     2.51934753e+02   # ~mu_R
   2000015     3.97573385e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.06320008e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.96463721e-01   # N_{1,1}
  1  2    -1.47776634e-02   # N_{1,2}
  1  3     7.62531284e-02   # N_{1,3}
  1  4    -3.20489268e-02   # N_{1,4}
  2  1     3.28392245e-02   # N_{2,1}
  2  2     9.73559789e-01   # N_{2,2}
  2  3    -1.87206548e-01   # N_{2,3}
  2  4     1.26714759e-01   # N_{2,4}
  3  1    -3.04436012e-02   # N_{3,1}
  3  2     4.42376340e-02   # N_{3,2}
  3  3     7.04031972e-01   # N_{3,3}
  3  4     7.08135016e-01   # N_{3,4}
  4  1    -7.10972932e-02   # N_{4,1}
  4  2     2.23620189e-01   # N_{4,2}
  4  3     6.80792296e-01   # N_{4,3}
  4  4    -6.93873933e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.63578906e-01   # U_{1,1}
  1  2    -2.67424180e-01   # U_{1,2}
  2  1     2.67424180e-01   # U_{2,1}
  2  2     9.63578906e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.83388881e-01   # V_{1,1}
  1  2    -1.81511179e-01   # V_{1,2}
  2  1     1.81511179e-01   # V_{2,1}
  2  2     9.83388881e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.04714672e-01   # F_{11}
  1  2     9.14443019e-01   # F_{12}
  2  1     9.14443019e-01   # F_{21}
  2  2    -4.04714672e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.78893820e-01   # F_{11}
  1  2     2.04369492e-01   # F_{12}
  2  1    -2.04369492e-01   # F_{21}
  2  2     9.78893820e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.29005584e-01   # F_{11}
  1  2     9.91643867e-01   # F_{12}
  2  1     9.91643867e-01   # F_{21}
  2  2    -1.29005584e-01   # F_{22}
Block gauge Q= 9.53955617e+02  # SM gauge couplings
     1     3.62600631e-01   # g'(Q)MSSM DRbar
     2     6.42519992e-01   # g(Q)MSSM DRbar
     3     1.05630116e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.53955617e+02  
  3  3     8.58162300e-01   # Yt(Q)MSSM DRbar
Block yd Q= 9.53955617e+02  
  3  3     1.34713034e-01   # Yb(Q)MSSM DRbar
Block ye Q= 9.53955617e+02  
  3  3     1.00453869e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 9.53955617e+02 # Higgs mixing parameters
     1     6.85472679e+02    # mu(Q)MSSM DRbar
     2     9.66316781e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44005297e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.36453499e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 9.53955617e+02  # MSSM DRbar SUSY breaking parameters
     1     2.31219569e+02      # M_1(Q)
     2     4.28191262e+02      # M_2(Q)
     3     1.21653519e+03      # M_3(Q)
    21     1.31959187e+05      # mH1^2(Q)
    22    -4.53883374e+05      # mH2^2(Q)
    31     3.88977041e+02      # meL(Q)
    32     3.88973403e+02      # mmuL(Q)
    33     3.87875049e+02      # mtauL(Q)
    34     2.44085713e+02      # meR(Q)
    35     2.44073919e+02      # mmuR(Q)
    36     2.40491758e+02      # mtauR(Q)
    41     1.11150426e+03      # mqL1(Q)
    42     1.11150141e+03      # mqL2(Q)
    43     1.02471121e+03      # mqL3(Q)
    44     1.06947135e+03      # muR(Q)
    45     1.06946831e+03      # mcR(Q)
    46     8.80808075e+02      # mtR(Q)
    47     1.06436287e+03      # mdR(Q)
    48     1.06435990e+03      # msR(Q)
    49     1.05893476e+03      # mbR(Q)
Block au Q= 9.53955617e+02  
  1  1    -1.24227304e+03      # Au(Q)MSSM DRbar
  2  2    -1.24226749e+03      # Ac(Q)MSSM DRbar
  3  3    -9.59891365e+02      # At(Q)MSSM DRbar
Block ad Q= 9.53955617e+02  
  1  1    -1.51822639e+03      # Ad(Q)MSSM DRbar
  2  2    -1.51822125e+03      # As(Q)MSSM DRbar
  3  3    -1.41909502e+03      # Ab(Q)MSSM DRbar
Block ae Q= 9.53955617e+02  
  1  1    -3.27857351e+02      # Ae(Q)MSSM DRbar
  2  2    -3.27851496e+02      # Amu(Q)MSSM DRbar
  3  3    -3.26084015e+02      # Atau(Q)MSSM DRbar
