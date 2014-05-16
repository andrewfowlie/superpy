# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.56562521e+03
# dtanbeta/dmu=-8.69395561e-01
# mu=7.54612172e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=9.72936511e-02
# dtanbeta/dm3sq=-3.00933225e-03
# m3sq= -2.59158982e+05
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
     1    3.30000000e+02   # m0
     2    5.00000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.87272935e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=7.06051903e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03906349e+01   # MW
        25     1.16991019e+02   # h0
        35     6.09531107e+02   # H0
        36     6.09545507e+02   # A0
        37     6.15189699e+02   # H+
   1000021     1.15306199e+03   # ~g
   1000022     2.07041789e+02   # ~neutralino(1)
   1000023     3.94704634e+02   # ~neutralino(2)
   1000024     3.94831637e+02   # ~chargino(1)
   1000025    -7.02935647e+02   # ~neutralino(3)
   1000035     7.10981947e+02   # ~neutralino(4)
   1000037     7.12167032e+02   # ~chargino(2)
   1000001     1.09880074e+03   # ~d_L
   1000002     1.09603236e+03   # ~u_L
   1000003     1.09876744e+03   # ~s_L
   1000004     1.09599896e+03   # ~c_L
   1000005     9.12841538e+02   # ~b_1
   1000006     7.76575972e+02   # ~t_1
   1000011     4.72695651e+02   # ~e_L
   1000012     4.65718746e+02   # ~nue_L
   1000013     4.72560378e+02   # ~mu_L
   1000014     4.65581581e+02   # ~numu_L
   1000015     2.08841587e+02   # ~stau_1
   1000016     4.19081038e+02   # ~nu_tau_L
   2000001     1.05637803e+03   # ~d_R
   2000002     1.05951758e+03   # ~u_R
   2000003     1.05631188e+03   # ~s_R
   2000004     1.05951342e+03   # ~c_R
   2000005     9.81127170e+02   # ~b_2
   2000006     9.89133496e+02   # ~t_2
   2000011     3.81824722e+02   # ~e_R
   2000013     3.81485016e+02   # ~mu_R
   2000015     4.49721044e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.61315363e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97115929e-01   # N_{1,1}
  1  2    -1.02792580e-02   # N_{1,2}
  1  3     7.15294536e-02   # N_{1,3}
  1  4    -2.31883046e-02   # N_{1,4}
  2  1     2.39132021e-02   # N_{2,1}
  2  2     9.82072716e-01   # N_{2,2}
  2  3    -1.61887375e-01   # N_{2,3}
  2  4     9.35618381e-02   # N_{2,4}
  3  1    -3.35081645e-02   # N_{3,1}
  3  2     4.93678544e-02   # N_{3,2}
  3  3     7.03707859e-01   # N_{3,3}
  3  4     7.07979708e-01   # N_{3,4}
  4  1    -6.37588035e-02   # N_{4,1}
  4  2     1.81632411e-01   # N_{4,2}
  4  3     6.88092483e-01   # N_{4,3}
  4  4    -6.99623626e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.73058746e-01   # U_{1,1}
  1  2    -2.30557318e-01   # U_{1,2}
  2  1     2.30557318e-01   # U_{2,1}
  2  2     9.73058746e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.91018962e-01   # V_{1,1}
  1  2    -1.33721419e-01   # V_{1,2}
  2  1     1.33721419e-01   # V_{2,1}
  2  2     9.91018962e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.73115530e-01   # F_{11}
  1  2     8.81000395e-01   # F_{12}
  2  1     8.81000395e-01   # F_{21}
  2  2    -4.73115530e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     8.12927559e-01   # F_{11}
  1  2     5.82364820e-01   # F_{12}
  2  1    -5.82364820e-01   # F_{21}
  2  2     8.12927559e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     3.55668182e-01   # F_{11}
  1  2     9.34612296e-01   # F_{12}
  2  1     9.34612296e-01   # F_{21}
  2  2    -3.55668182e-01   # F_{22}
Block gauge Q= 8.50961572e+02  # SM gauge couplings
     1     3.62132827e-01   # g'(Q)MSSM DRbar
     2     6.42336607e-01   # g(Q)MSSM DRbar
     3     1.06123658e+00   # g3(Q)MSSM DRbar
Block yu Q= 8.50961572e+02  
  3  3     8.56328253e-01   # Yt(Q)MSSM DRbar
Block yd Q= 8.50961572e+02  
  3  3     4.91160204e-01   # Yb(Q)MSSM DRbar
Block ye Q= 8.50961572e+02  
  3  3     4.27332013e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 8.50961572e+02 # Higgs mixing parameters
     1     6.99525487e+02    # mu(Q)MSSM DRbar
     2     3.92347808e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44010951e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     4.69191051e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 8.50961572e+02  # MSSM DRbar SUSY breaking parameters
     1     2.10482988e+02      # M_1(Q)
     2     3.90867401e+02      # M_2(Q)
     3     1.11771426e+03      # M_3(Q)
    21    -1.16846912e+05      # mH1^2(Q)
    22    -4.83190981e+05      # mH2^2(Q)
    31     4.67134701e+02      # meL(Q)
    32     4.66997452e+02      # mmuL(Q)
    33     4.23264884e+02      # mtauL(Q)
    34     3.77269158e+02      # meR(Q)
    35     3.76925410e+02      # mmuR(Q)
    36     2.51028413e+02      # mtauR(Q)
    41     1.06388358e+03      # mqL1(Q)
    42     1.06384958e+03      # mqL2(Q)
    43     9.10617504e+02      # mqL3(Q)
    44     1.02772040e+03      # muR(Q)
    45     1.02771617e+03      # mcR(Q)
    46     7.90834180e+02      # mtR(Q)
    47     1.02339875e+03      # mdR(Q)
    48     1.02333145e+03      # msR(Q)
    49     9.28907616e+02      # mbR(Q)
Block au Q= 8.50961572e+02  
  1  1    -1.48305137e+03      # Au(Q)MSSM DRbar
  2  2    -1.48301491e+03      # Ac(Q)MSSM DRbar
  3  3    -1.02166876e+03      # At(Q)MSSM DRbar
Block ad Q= 8.50961572e+02  
  1  1    -1.72568542e+03      # Ad(Q)MSSM DRbar
  2  2    -1.72559256e+03      # As(Q)MSSM DRbar
  3  3    -1.44972182e+03      # Ab(Q)MSSM DRbar
Block ae Q= 8.50961572e+02  
  1  1    -6.25682925e+02      # Ae(Q)MSSM DRbar
  2  2    -6.25368906e+02      # Amu(Q)MSSM DRbar
  3  3    -5.20590962e+02      # Atau(Q)MSSM DRbar
