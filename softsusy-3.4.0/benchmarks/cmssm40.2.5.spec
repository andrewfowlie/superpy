# SuperPy: Jacobian for naturalness priors.
# J = 1.47257624e+05
# b = 1.31647609e+05
# Mu = 9.11876000e+01
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
     1    7.50000000e+02   # m0
     2    6.50000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.79840900e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=8.97819232e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03886366e+01   # MW
        25     1.18443308e+02   # h0
        35     8.65273057e+02   # H0
        36     8.65265228e+02   # A0
        37     8.69291002e+02   # H+
   1000021     1.49024590e+03   # ~g
   1000022     2.73950154e+02   # ~neutralino(1)
   1000023     5.21376206e+02   # ~neutralino(2)
   1000024     5.21515841e+02   # ~chargino(1)
   1000025    -8.51062283e+02   # ~neutralino(3)
   1000035     8.59250067e+02   # ~neutralino(4)
   1000037     8.60102925e+02   # ~chargino(2)
   1000001     1.52213614e+03   # ~d_L
   1000002     1.52021001e+03   # ~u_L
   1000003     1.52209258e+03   # ~s_L
   1000004     1.52016638e+03   # ~c_L
   1000005     1.27098703e+03   # ~b_1
   1000006     1.08717259e+03   # ~t_1
   1000011     8.65185249e+02   # ~e_L
   1000012     8.61214048e+02   # ~nue_L
   1000013     8.64978780e+02   # ~mu_L
   1000014     8.61006715e+02   # ~numu_L
   1000015     6.09215716e+02   # ~stau_1
   1000016     7.89895057e+02   # ~nu_tau_L
   2000001     1.47304581e+03   # ~d_R
   2000002     1.47720366e+03   # ~u_R
   2000003     1.47296015e+03   # ~s_R
   2000004     1.47719805e+03   # ~c_R
   2000005     1.34730945e+03   # ~b_2
   2000006     1.31986691e+03   # ~t_2
   2000011     7.88645398e+02   # ~e_R
   2000013     7.88188348e+02   # ~mu_R
   2000015     8.03985206e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.59197685e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97966821e-01   # N_{1,1}
  1  2    -7.00710774e-03   # N_{1,2}
  1  3     5.97993158e-02   # N_{1,3}
  1  4    -2.09085031e-02   # N_{1,4}
  2  1     1.75984352e-02   # N_{2,1}
  2  2     9.84844451e-01   # N_{2,2}
  2  3    -1.46317435e-01   # N_{2,3}
  2  4     9.14489567e-02   # N_{2,4}
  3  1    -2.71348550e-02   # N_{3,1}
  3  2     3.94999021e-02   # N_{3,2}
  3  3     7.04907217e-01   # N_{3,3}
  3  4     7.07678792e-01   # N_{3,4}
  4  1    -5.49201052e-02   # N_{4,1}
  4  2     1.68736678e-01   # N_{4,2}
  4  3     6.91462989e-01   # N_{4,3}
  4  4    -7.00278980e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.78026115e-01   # U_{1,1}
  1  2    -2.08482417e-01   # U_{1,2}
  2  1     2.08482417e-01   # U_{2,1}
  2  2     9.78026115e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.91407083e-01   # V_{1,1}
  1  2    -1.30812827e-01   # V_{1,2}
  2  1     1.30812827e-01   # V_{2,1}
  2  2     9.91407083e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.67716988e-01   # F_{11}
  1  2     9.29937749e-01   # F_{12}
  2  1     9.29937749e-01   # F_{21}
  2  2    -3.67716988e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.13093173e-01   # F_{11}
  1  2     4.07750975e-01   # F_{12}
  2  1    -4.07750975e-01   # F_{21}
  2  2     9.13093173e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     2.37663360e-01   # F_{11}
  1  2     9.71347583e-01   # F_{12}
  2  1     9.71347583e-01   # F_{21}
  2  2    -2.37663360e-01   # F_{22}
Block gauge Q= 1.16443024e+03  # SM gauge couplings
     1     3.62625130e-01   # g'(Q)MSSM DRbar
     2     6.40749660e-01   # g(Q)MSSM DRbar
     3     1.04601378e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.16443024e+03  
  3  3     8.46875216e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.16443024e+03  
  3  3     4.95567015e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.16443024e+03  
  3  3     4.23709286e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.16443024e+03 # Higgs mixing parameters
     1     8.46596015e+02    # mu(Q)MSSM DRbar
     2     3.91663556e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43672839e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     9.57107132e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.16443024e+03  # MSSM DRbar SUSY breaking parameters
     1     2.77018429e+02      # M_1(Q)
     2     5.10792770e+02      # M_2(Q)
     3     1.42163473e+03      # M_3(Q)
    21     4.46851166e+04      # mH1^2(Q)
    22    -6.99015972e+05      # mH2^2(Q)
    31     8.60786195e+02      # meL(Q)
    32     8.60578691e+02      # mmuL(Q)
    33     7.91995865e+02      # mtauL(Q)
    34     7.85602254e+02      # meR(Q)
    35     7.85143917e+02      # mmuR(Q)
    36     6.22395668e+02      # mtauR(Q)
    41     1.47887380e+03      # mqL1(Q)
    42     1.47882876e+03      # mqL2(Q)
    43     1.25025188e+03      # mqL3(Q)
    44     1.43647434e+03      # muR(Q)
    45     1.43646857e+03      # mcR(Q)
    46     1.07896199e+03      # mtR(Q)
    47     1.43137625e+03      # mdR(Q)
    48     1.43128791e+03      # msR(Q)
    49     1.29697418e+03      # mbR(Q)
Block au Q= 1.16443024e+03  
  1  1    -1.78111537e+03      # Au(Q)MSSM DRbar
  2  2    -1.78107370e+03      # Ac(Q)MSSM DRbar
  3  3    -1.25136021e+03      # At(Q)MSSM DRbar
Block ad Q= 1.16443024e+03  
  1  1    -2.05354986e+03      # Ad(Q)MSSM DRbar
  2  2    -2.05344379e+03      # As(Q)MSSM DRbar
  3  3    -1.73151976e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.16443024e+03  
  1  1    -6.82580306e+02      # Ae(Q)MSSM DRbar
  2  2    -6.82247893e+02      # Amu(Q)MSSM DRbar
  3  3    -5.72739865e+02      # Atau(Q)MSSM DRbar
