# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.21776634e+03
# dtanbeta/dmu=-3.88388901e-01
# mu=6.67248959e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=9.48396464e-02
# dtanbeta/dm3sq=-1.31435471e-03
# m3sq= -2.09107597e+05
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
     1    1.20000000e+03   # m0
     2    4.50000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    2.08986760e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=3.61914103e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03908418e+01   # MW
        25     1.16883551e+02   # h0
        35     9.42583787e+02   # H0
        36     9.42686893e+02   # A0
        37     9.46262161e+02   # H+
   1000021     1.10074813e+03   # ~g
   1000022     1.88515051e+02   # ~neutralino(1)
   1000023     3.60117940e+02   # ~neutralino(2)
   1000024     3.60225266e+02   # ~chargino(1)
   1000025    -6.19932849e+02   # ~neutralino(3)
   1000035     6.30467870e+02   # ~neutralino(4)
   1000037     6.31565790e+02   # ~chargino(2)
   1000001     1.50618564e+03   # ~d_L
   1000002     1.50421681e+03   # ~u_L
   1000003     1.50613384e+03   # ~s_L
   1000004     1.50416494e+03   # ~c_L
   1000005     1.18351865e+03   # ~b_1
   1000006     9.83201296e+02   # ~t_1
   1000011     1.23341740e+03   # ~e_L
   1000012     1.23051511e+03   # ~nue_L
   1000013     1.23311570e+03   # ~mu_L
   1000014     1.23021278e+03   # ~numu_L
   1000015     9.95410495e+02   # ~stau_1
   1000016     1.13137954e+03   # ~nu_tau_L
   2000001     1.48466408e+03   # ~d_R
   2000002     1.48591561e+03   # ~u_R
   2000003     1.48456507e+03   # ~s_R
   2000004     1.48590886e+03   # ~c_R
   2000005     1.31116039e+03   # ~b_2
   2000006     1.21425600e+03   # ~t_2
   2000011     1.21108624e+03   # ~e_R
   2000013     1.21046943e+03   # ~mu_R
   2000015     1.13751873e+03   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.59229728e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.96172111e-01   # N_{1,1}
  1  2    -1.28121943e-02   # N_{1,2}
  1  3     8.20645273e-02   # N_{1,3}
  1  4    -2.72467784e-02   # N_{1,4}
  2  1     3.11340102e-02   # N_{2,1}
  2  2     9.75152975e-01   # N_{2,2}
  2  3    -1.88400859e-01   # N_{2,3}
  2  4     1.12305236e-01   # N_{2,4}
  3  1    -3.78232922e-02   # N_{3,1}
  3  2     5.54567550e-02   # N_{3,2}
  3  3     7.02897273e-01   # N_{3,3}
  3  4     7.08116777e-01   # N_{3,4}
  4  1    -7.23961178e-02   # N_{4,1}
  4  2     2.14095942e-01   # N_{4,2}
  4  3     6.80959582e-01   # N_{4,3}
  4  4    -6.96574316e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.63200557e-01   # U_{1,1}
  1  2    -2.68783718e-01   # U_{1,2}
  2  1     2.68783718e-01   # U_{2,1}
  2  2     9.63200557e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.86988307e-01   # V_{1,1}
  1  2    -1.60792045e-01   # V_{1,2}
  2  1     1.60792045e-01   # V_{2,1}
  2  2     9.86988307e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.95788423e-01   # F_{11}
  1  2     9.55253479e-01   # F_{12}
  2  1     9.55253479e-01   # F_{21}
  2  2    -2.95788423e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.84267117e-01   # F_{11}
  1  2     1.76686849e-01   # F_{12}
  2  1    -1.76686849e-01   # F_{21}
  2  2     9.84267117e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.51220707e-01   # F_{11}
  1  2     9.88500024e-01   # F_{12}
  2  1     9.88500024e-01   # F_{21}
  2  2    -1.51220707e-01   # F_{22}
Block gauge Q= 1.06740630e+03  # SM gauge couplings
     1     3.62258120e-01   # g'(Q)MSSM DRbar
     2     6.41713939e-01   # g(Q)MSSM DRbar
     3     1.05405454e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.06740630e+03  
  3  3     8.54178468e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.06740630e+03  
  3  3     5.16088204e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.06740630e+03  
  3  3     4.12825363e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.06740630e+03 # Higgs mixing parameters
     1     6.12070089e+02    # mu(Q)MSSM DRbar
     2     3.92071619e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43605024e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.07006209e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.06740630e+03  # MSSM DRbar SUSY breaking parameters
     1     1.90398248e+02      # M_1(Q)
     2     3.52809358e+02      # M_2(Q)
     3     9.97158048e+02      # M_3(Q)
    21     5.37923141e+05      # mH1^2(Q)
    22    -3.57904207e+05      # mH2^2(Q)
    31     1.23047107e+03      # meL(Q)
    32     1.23016953e+03      # mmuL(Q)
    33     1.13384691e+03      # mtauL(Q)
    34     1.20893262e+03      # meR(Q)
    35     1.20831563e+03      # mmuR(Q)
    36     1.00067197e+03      # mtauR(Q)
    41     1.47876436e+03      # mqL1(Q)
    42     1.47871165e+03      # mqL2(Q)
    43     1.16479953e+03      # mqL3(Q)
    44     1.46164057e+03      # muR(Q)
    45     1.46163369e+03      # mcR(Q)
    46     9.68011005e+02      # mtR(Q)
    47     1.45972968e+03      # mdR(Q)
    48     1.45962875e+03      # msR(Q)
    49     1.28507828e+03      # mbR(Q)
Block au Q= 1.06740630e+03  
  1  1    -1.35325121e+03      # Au(Q)MSSM DRbar
  2  2    -1.35321767e+03      # Ac(Q)MSSM DRbar
  3  3    -9.19232591e+02      # At(Q)MSSM DRbar
Block ad Q= 1.06740630e+03  
  1  1    -1.56517920e+03      # Ad(Q)MSSM DRbar
  2  2    -1.56509379e+03      # As(Q)MSSM DRbar
  3  3    -1.29428583e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.06740630e+03  
  1  1    -5.92607887e+02      # Ae(Q)MSSM DRbar
  2  2    -5.92303290e+02      # Amu(Q)MSSM DRbar
  3  3    -4.97393937e+02      # Atau(Q)MSSM DRbar
