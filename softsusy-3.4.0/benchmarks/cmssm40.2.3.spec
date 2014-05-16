# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.73522587e+03
# dtanbeta/dmu=-6.03256556e-01
# mu=8.06383475e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=9.64454374e-02
# dtanbeta/dm3sq=-2.06542477e-03
# m3sq= -2.78268628e+05
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
     1    6.00000000e+02   # m0
     2    5.50000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.88704263e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=7.86821589e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03897995e+01   # MW
        25     1.17592526e+02   # h0
        35     7.31019979e+02   # H0
        36     7.31025955e+02   # A0
        37     7.35791157e+02   # H+
   1000021     1.27464309e+03   # ~g
   1000022     2.29843540e+02   # ~neutralino(1)
   1000023     4.38019454e+02   # ~neutralino(2)
   1000024     4.38150979e+02   # ~chargino(1)
   1000025    -7.50577825e+02   # ~neutralino(3)
   1000035     7.58905227e+02   # ~neutralino(4)
   1000037     7.59941039e+02   # ~chargino(2)
   1000001     1.28682549e+03   # ~d_L
   1000002     1.28451660e+03   # ~u_L
   1000003     1.28678715e+03   # ~s_L
   1000004     1.28447817e+03   # ~c_L
   1000005     1.06783522e+03   # ~b_1
   1000006     9.07726880e+02   # ~t_1
   1000011     7.03599221e+02   # ~e_L
   1000012     6.98786800e+02   # ~nue_L
   1000013     7.03420983e+02   # ~mu_L
   1000014     6.98607438e+02   # ~numu_L
   1000015     4.72304350e+02   # ~stau_1
   1000016     6.37453271e+02   # ~nu_tau_L
   2000001     1.24462197e+03   # ~d_R
   2000002     1.24793471e+03   # ~u_R
   2000003     1.24454620e+03   # ~s_R
   2000004     1.24792981e+03   # ~c_R
   2000005     1.14078537e+03   # ~b_2
   2000006     1.12642019e+03   # ~t_2
   2000011     6.35100516e+02   # ~e_R
   2000013     6.34701827e+02   # ~mu_R
   2000015     6.55539382e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.59913465e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97437285e-01   # N_{1,1}
  1  2    -8.97896633e-03   # N_{1,2}
  1  3     6.73111905e-02   # N_{1,3}
  1  4    -2.25265295e-02   # N_{1,4}
  2  1     2.15733182e-02   # N_{2,1}
  2  2     9.82863091e-01   # N_{2,2}
  2  3    -1.57135624e-01   # N_{2,3}
  2  4     9.39315292e-02   # N_{2,4}
  3  1    -3.11252527e-02   # N_{3,1}
  3  2     4.56185747e-02   # N_{3,2}
  3  3     7.04199700e-01   # N_{3,3}
  3  4     7.07850936e-01   # N_{3,4}
  4  1    -6.07015101e-02   # N_{4,1}
  4  2     1.78377320e-01   # N_{4,2}
  4  3     6.89115652e-01   # N_{4,3}
  4  4    -6.99726001e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.74611560e-01   # U_{1,1}
  1  2    -2.23902450e-01   # U_{1,2}
  2  1     2.23902450e-01   # U_{2,1}
  2  2     9.74611560e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.90932521e-01   # V_{1,1}
  1  2    -1.34360483e-01   # V_{1,2}
  2  1     1.34360483e-01   # V_{2,1}
  2  2     9.90932521e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.16083862e-01   # F_{11}
  1  2     9.09326245e-01   # F_{12}
  2  1     9.09326245e-01   # F_{21}
  2  2    -4.16083862e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     8.88645484e-01   # F_{11}
  1  2     4.58594815e-01   # F_{12}
  2  1    -4.58594815e-01   # F_{21}
  2  2     8.88645484e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     2.82873237e-01   # F_{11}
  1  2     9.59157303e-01   # F_{12}
  2  1     9.59157303e-01   # F_{21}
  2  2    -2.82873237e-01   # F_{22}
Block gauge Q= 9.81985216e+02  # SM gauge couplings
     1     3.62276455e-01   # g'(Q)MSSM DRbar
     2     6.41583979e-01   # g(Q)MSSM DRbar
     3     1.05451478e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.81985216e+02  
  3  3     8.52295843e-01   # Yt(Q)MSSM DRbar
Block yd Q= 9.81985216e+02  
  3  3     4.95397908e-01   # Yb(Q)MSSM DRbar
Block ye Q= 9.81985216e+02  
  3  3     4.24552002e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 9.81985216e+02 # Higgs mixing parameters
     1     7.46393564e+02    # mu(Q)MSSM DRbar
     2     3.92051622e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43871758e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.83734791e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 9.81985216e+02  # MSSM DRbar SUSY breaking parameters
     1     2.32721305e+02      # M_1(Q)
     2     4.30923109e+02      # M_2(Q)
     3     1.21817757e+03      # M_3(Q)
    21    -1.33708790e+04      # mH1^2(Q)
    22    -5.46547302e+05      # mH2^2(Q)
    31     6.99420159e+02      # meL(Q)
    32     6.99240820e+02      # mmuL(Q)
    33     6.40421647e+02      # mtauL(Q)
    34     6.32007153e+02      # meR(Q)
    35     6.31606905e+02      # mmuR(Q)
    36     4.89458644e+02      # mtauR(Q)
    41     1.24910242e+03      # mqL1(Q)
    42     1.24906278e+03      # mqL2(Q)
    43     1.05428767e+03      # mqL3(Q)
    44     1.21280744e+03      # muR(Q)
    45     1.21280241e+03      # mcR(Q)
    46     9.08903112e+02      # mtR(Q)
    47     1.20847782e+03      # mdR(Q)
    48     1.20839999e+03      # msR(Q)
    49     1.09287197e+03      # mbR(Q)
Block au Q= 9.81985216e+02  
  1  1    -1.58110853e+03      # Au(Q)MSSM DRbar
  2  2    -1.58107042e+03      # Ac(Q)MSSM DRbar
  3  3    -1.09681786e+03      # At(Q)MSSM DRbar
Block ad Q= 9.81985216e+02  
  1  1    -1.83179024e+03      # Ad(Q)MSSM DRbar
  2  2    -1.83169317e+03      # As(Q)MSSM DRbar
  3  3    -1.53899676e+03      # Ab(Q)MSSM DRbar
Block ae Q= 9.81985216e+02  
  1  1    -6.43164701e+02      # Ae(Q)MSSM DRbar
  2  2    -6.42844799e+02      # Amu(Q)MSSM DRbar
  3  3    -5.37278967e+02      # Atau(Q)MSSM DRbar
