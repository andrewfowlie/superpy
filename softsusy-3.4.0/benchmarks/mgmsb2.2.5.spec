# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-3.32001799e+03
# dtanbeta/dmu=3.59844696e-02
# mu=8.34922844e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=2.79240679e-01
# dtanbeta/dm3sq=-1.79779409e-04
# m3sq= 5.87932478e+04
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
     1    1.60000000e+05   # lambda
     2    1.00000000e+14   # M_mess
     5    1.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.13224893e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03958348e+01   # MW
        25     1.15798105e+02   # h0
        35     1.05794785e+03   # H0
        36     1.05790418e+03   # A0
        37     1.06115333e+03   # H+
   1000021     1.19552388e+03   # ~g
   1000022     2.13897769e+02   # ~neutralino(1)
   1000023     4.12140706e+02   # ~neutralino(2)
   1000024     4.12272385e+02   # ~chargino(1)
   1000025    -8.31789455e+02   # ~neutralino(3)
   1000035     8.38638770e+02   # ~neutralino(4)
   1000037     8.39411223e+02   # ~chargino(2)
   1000039     3.76491429e+00   # ~gravitino
   1000001     1.46950089e+03   # ~d_L
   1000002     1.46750972e+03   # ~u_L
   1000003     1.46949507e+03   # ~s_L
   1000004     1.46750390e+03   # ~c_L
   1000005     1.29925295e+03   # ~b_1
   1000006     1.03325507e+03   # ~t_1
   1000011     7.09724719e+02   # ~e_L
   1000012     7.04972862e+02   # ~nue_L
   1000013     7.09709741e+02   # ~mu_L
   1000014     7.04957793e+02   # ~numu_L
   1000015     4.66290118e+02   # ~stau_1
   1000016     7.00042763e+02   # ~nu_tau_L
   2000001     1.32104329e+03   # ~d_R
   2000002     1.34815335e+03   # ~u_R
   2000003     1.32103465e+03   # ~s_R
   2000004     1.34814907e+03   # ~c_R
   2000005     1.33051128e+03   # ~b_2
   2000006     1.34216269e+03   # ~t_2
   2000011     4.82458906e+02   # ~e_R
   2000013     4.82414528e+02   # ~mu_R
   2000015     7.06079595e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.00914391e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.98006959e-01   # N_{1,1}
  1  2    -8.99168498e-03   # N_{1,2}
  1  3     5.94225488e-02   # N_{1,3}
  1  4    -1.92411215e-02   # N_{1,4}
  2  1     1.77840555e-02   # N_{2,1}
  2  2     9.89395604e-01   # N_{2,2}
  2  3    -1.26722027e-01   # N_{2,3}
  2  4     6.87138609e-02   # N_{2,4}
  3  1    -2.79397273e-02   # N_{3,1}
  3  2     4.16097905e-02   # N_{3,2}
  3  3     7.04708775e-01   # N_{3,3}
  3  4     7.07724197e-01   # N_{3,4}
  4  1    -5.37141459e-02   # N_{4,1}
  4  2     1.38867254e-01   # N_{4,2}
  4  3     6.95554477e-01   # N_{4,3}
  4  4    -7.02875982e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.83600470e-01   # U_{1,1}
  1  2    -1.80361067e-01   # U_{1,2}
  2  1     1.80361067e-01   # U_{2,1}
  2  2     9.83600470e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.95169270e-01   # V_{1,1}
  1  2    -9.81739463e-02   # V_{1,2}
  2  1     9.81739463e-02   # V_{2,1}
  2  2     9.95169270e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     1.87088433e-01   # F_{11}
  1  2     9.82343076e-01   # F_{12}
  2  1     9.82343076e-01   # F_{21}
  2  2    -1.87088433e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     4.10010702e-01   # F_{11}
  1  2     9.12080712e-01   # F_{12}
  2  1     9.12080712e-01   # F_{21}
  2  2    -4.10010702e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     7.87747956e-02   # F_{11}
  1  2     9.96892437e-01   # F_{12}
  2  1     9.96892437e-01   # F_{21}
  2  2    -7.87747956e-02   # F_{22}
Block gauge Q= 1.14416565e+03  # SM gauge couplings
     1     3.62806274e-01   # g'(Q)MSSM DRbar
     2     6.41730817e-01   # g(Q)MSSM DRbar
     3     1.05131142e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.14416565e+03  
  3  3     8.54040188e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.14416565e+03  
  3  3     1.96983343e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.14416565e+03  
  3  3     1.50685043e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.14416565e+03 # Higgs mixing parameters
     1     8.25430010e+02    # mu(Q)MSSM DRbar
     2     1.44798149e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43807352e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.16822997e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.14416565e+03  # MSSM DRbar SUSY breaking parameters
     1     2.17279027e+02      # M_1(Q)
     2     4.02026152e+02      # M_2(Q)
     3     1.09852146e+03      # M_3(Q)
    21     4.26650418e+05      # mH1^2(Q)
    22    -6.58236401e+05      # mH2^2(Q)
    31     7.06001177e+02      # meL(Q)
    32     7.05986125e+02      # mmuL(Q)
    33     7.01392160e+02      # mtauL(Q)
    34     4.77466346e+02      # meR(Q)
    35     4.77421533e+02      # mmuR(Q)
    36     4.63583911e+02      # mtauR(Q)
    41     1.43615823e+03      # mqL1(Q)
    42     1.43615225e+03      # mqL2(Q)
    43     1.29673327e+03      # mqL3(Q)
    44     1.31512797e+03      # muR(Q)
    45     1.31512356e+03      # mcR(Q)
    46     9.98572388e+02      # mtR(Q)
    47     1.28584027e+03      # mdR(Q)
    48     1.28583131e+03      # msR(Q)
    49     1.26949134e+03      # mbR(Q)
Block au Q= 1.14416565e+03  
  1  1    -1.01776165e+03      # Au(Q)MSSM DRbar
  2  2    -1.01775654e+03      # Ac(Q)MSSM DRbar
  3  3    -8.10600074e+02      # At(Q)MSSM DRbar
Block ad Q= 1.14416565e+03  
  1  1    -1.21491240e+03      # Ad(Q)MSSM DRbar
  2  2    -1.21490534e+03      # As(Q)MSSM DRbar
  3  3    -1.13693491e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.14416565e+03  
  1  1    -2.26843939e+02      # Ae(Q)MSSM DRbar
  2  2    -2.26836186e+02      # Amu(Q)MSSM DRbar
  3  3    -2.24476627e+02      # Atau(Q)MSSM DRbar
