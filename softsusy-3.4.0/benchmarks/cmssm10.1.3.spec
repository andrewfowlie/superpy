# SuperPy: Jacobian for naturalness priors.
# J = 2.95241733e+05
# b = 1.24438082e+05
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
     3    1.00000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    1.50000000e+02   # m0
     2    6.00000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.76385152e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.11810980e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03945125e+01   # MW
        25     1.16085532e+02   # h0
        35     8.48769229e+02   # H0
        36     8.48536746e+02   # A0
        37     8.52558307e+02   # H+
   1000021     1.35667045e+03   # ~g
   1000022     2.47689445e+02   # ~neutralino(1)
   1000023     4.68638801e+02   # ~neutralino(2)
   1000024     4.68701491e+02   # ~chargino(1)
   1000025    -7.45754754e+02   # ~neutralino(3)
   1000035     7.58164403e+02   # ~neutralino(4)
   1000037     7.58410228e+02   # ~chargino(2)
   1000001     1.24670821e+03   # ~d_L
   1000002     1.24433108e+03   # ~u_L
   1000003     1.24670519e+03   # ~s_L
   1000004     1.24432805e+03   # ~c_L
   1000005     1.14234954e+03   # ~b_1
   1000006     9.58226034e+02   # ~t_1
   1000011     4.31793771e+02   # ~e_L
   1000012     4.24323007e+02   # ~nue_L
   1000013     4.31789898e+02   # ~mu_L
   1000014     4.24319067e+02   # ~numu_L
   1000015     2.66973739e+02   # ~stau_1
   1000016     4.22941983e+02   # ~nu_tau_L
   2000001     1.19293706e+03   # ~d_R
   2000002     1.19736464e+03   # ~u_R
   2000003     1.19293391e+03   # ~s_R
   2000004     1.19736140e+03   # ~c_R
   2000005     1.18847632e+03   # ~b_2
   2000006     1.18197076e+03   # ~t_2
   2000011     2.74085051e+02   # ~e_R
   2000013     2.74072669e+02   # ~mu_R
   2000015     4.32299813e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05966217e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.96966160e-01   # N_{1,1}
  1  2    -1.26748387e-02   # N_{1,2}
  1  3     7.06798503e-02   # N_{1,3}
  1  4    -3.00363479e-02   # N_{1,4}
  2  1     2.86653150e-02   # N_{2,1}
  2  2     9.76134566e-01   # N_{2,2}
  2  3    -1.77579340e-01   # N_{2,3}
  2  4     1.21676566e-01   # N_{2,4}
  3  1    -2.80967645e-02   # N_{3,1}
  3  2     4.07347919e-02   # N_{3,2}
  3  3     7.04489163e-01   # N_{3,3}
  3  4     7.07987478e-01   # N_{3,4}
  4  1    -6.66884286e-02   # N_{4,1}
  4  2     2.12935047e-01   # N_{4,2}
  4  3     6.83494665e-01   # N_{4,3}
  4  4    -6.95015369e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.67359916e-01   # U_{1,1}
  1  2    -2.53406378e-01   # U_{1,2}
  2  1     2.53406378e-01   # U_{2,1}
  2  2     9.67359916e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.84730687e-01   # V_{1,1}
  1  2    -1.74084674e-01   # V_{1,2}
  2  1     1.74084674e-01   # V_{2,1}
  2  2     9.84730687e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.84606212e-01   # F_{11}
  1  2     9.23080745e-01   # F_{12}
  2  1     9.23080745e-01   # F_{21}
  2  2    -3.84606212e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.81375343e-01   # F_{11}
  1  2     1.92100069e-01   # F_{12}
  2  1    -1.92100069e-01   # F_{21}
  2  2     9.81375343e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.18042583e-01   # F_{11}
  1  2     9.93008534e-01   # F_{12}
  2  1     9.93008534e-01   # F_{21}
  2  2    -1.18042583e-01   # F_{22}
Block gauge Q= 1.03280996e+03  # SM gauge couplings
     1     3.62766868e-01   # g'(Q)MSSM DRbar
     2     6.42082769e-01   # g(Q)MSSM DRbar
     3     1.05212451e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.03280996e+03  
  3  3     8.55437331e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.03280996e+03  
  3  3     1.34162963e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.03280996e+03  
  3  3     1.00398894e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.03280996e+03 # Higgs mixing parameters
     1     7.40337327e+02    # mu(Q)MSSM DRbar
     2     9.65379787e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43912207e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     7.46252739e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.03280996e+03  # MSSM DRbar SUSY breaking parameters
     1     2.53207001e+02      # M_1(Q)
     2     4.67936400e+02      # M_2(Q)
     3     1.31909531e+03      # M_3(Q)
    21     1.56949245e+05      # mH1^2(Q)
    22    -5.28306184e+05      # mH2^2(Q)
    31     4.23672760e+02      # meL(Q)
    32     4.23668817e+02      # mmuL(Q)
    33     4.22477606e+02      # mtauL(Q)
    34     2.66158977e+02      # meR(Q)
    35     2.66146212e+02      # mmuR(Q)
    36     2.62265732e+02      # mtauR(Q)
    41     1.20403494e+03      # mqL1(Q)
    42     1.20403186e+03      # mqL2(Q)
    43     1.11028522e+03      # mqL3(Q)
    44     1.15796863e+03      # muR(Q)
    45     1.15796535e+03      # mcR(Q)
    46     9.54084783e+02      # mtR(Q)
    47     1.15234621e+03      # mdR(Q)
    48     1.15234300e+03      # msR(Q)
    49     1.14650293e+03      # mbR(Q)
Block au Q= 1.03280996e+03  
  1  1    -1.34291481e+03      # Au(Q)MSSM DRbar
  2  2    -1.34290883e+03      # Ac(Q)MSSM DRbar
  3  3    -1.03903288e+03      # At(Q)MSSM DRbar
Block ad Q= 1.03280996e+03  
  1  1    -1.63958179e+03      # Ad(Q)MSSM DRbar
  2  2    -1.63957626e+03      # As(Q)MSSM DRbar
  3  3    -1.53291558e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.03280996e+03  
  1  1    -3.56259932e+02      # Ae(Q)MSSM DRbar
  2  2    -3.56253602e+02      # Amu(Q)MSSM DRbar
  3  3    -3.54341104e+02      # Atau(Q)MSSM DRbar
