# SuperPy: Jacobian for naturalness priors.
# J = 7.37059483e+04
# b = 8.15469777e+04
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
     1    6.00000000e+02   # m0
     2    5.00000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.95484135e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=6.15438300e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03905663e+01   # MW
        25     1.17107081e+02   # h0
        35     6.90304987e+02   # H0
        36     6.90309473e+02   # A0
        37     6.95366662e+02   # H+
   1000021     1.17095700e+03   # ~g
   1000022     2.08190433e+02   # ~neutralino(1)
   1000023     3.97021907e+02   # ~neutralino(2)
   1000024     3.97146958e+02   # ~chargino(1)
   1000025    -6.97379040e+02   # ~neutralino(3)
   1000035     7.05939302e+02   # ~neutralino(4)
   1000037     7.07079490e+02   # ~chargino(2)
   1000001     1.20234550e+03   # ~d_L
   1000002     1.19986771e+03   # ~u_L
   1000003     1.20230839e+03   # ~s_L
   1000004     1.19983050e+03   # ~c_L
   1000005     9.90190961e+02   # ~b_1
   1000006     8.35837858e+02   # ~t_1
   1000011     6.86803108e+02   # ~e_L
   1000012     6.81870419e+02   # ~nue_L
   1000013     6.86624716e+02   # ~mu_L
   1000014     6.81690839e+02   # ~numu_L
   1000015     4.70490494e+02   # ~stau_1
   1000016     6.20896580e+02   # ~nu_tau_L
   2000001     1.16514056e+03   # ~d_R
   2000002     1.16787318e+03   # ~u_R
   2000003     1.16506768e+03   # ~s_R
   2000004     1.16786847e+03   # ~c_R
   2000005     1.06486132e+03   # ~b_2
   2000006     1.05016419e+03   # ~t_2
   2000011     6.29304443e+02   # ~e_R
   2000013     6.28911665e+02   # ~mu_R
   2000015     6.38897826e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.60195080e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97049920e-01   # N_{1,1}
  1  2    -1.03714064e-02   # N_{1,2}
  1  3     7.22800980e-02   # N_{1,3}
  1  4    -2.36532960e-02   # N_{1,4}
  2  1     2.44482388e-02   # N_{2,1}
  2  2     9.81285459e-01   # N_{2,2}
  2  3    -1.64903437e-01   # N_{2,3}
  2  4     9.63742036e-02   # N_{2,4}
  3  1    -3.37028641e-02   # N_{3,1}
  3  2     4.95636751e-02   # N_{3,2}
  3  3     7.03696598e-01   # N_{3,3}
  3  4     7.07967977e-01   # N_{3,4}
  4  1    -6.44814518e-02   # N_{4,1}
  4  2     1.85781386e-01   # N_{4,2}
  4  3     6.87308913e-01   # N_{4,3}
  4  4    -6.99238070e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.71991886e-01   # U_{1,1}
  1  2    -2.35014411e-01   # U_{1,2}
  2  1     2.35014411e-01   # U_{2,1}
  2  2     9.71991886e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.90449508e-01   # V_{1,1}
  1  2    -1.37875932e-01   # V_{1,2}
  2  1     1.37875932e-01   # V_{2,1}
  2  2     9.90449508e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.28724525e-01   # F_{11}
  1  2     9.03435267e-01   # F_{12}
  2  1     9.03435267e-01   # F_{21}
  2  2    -4.28724525e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     8.95885782e-01   # F_{11}
  1  2     4.44284443e-01   # F_{12}
  2  1    -4.44284443e-01   # F_{21}
  2  2     8.95885782e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     2.91472636e-01   # F_{11}
  1  2     9.56579167e-01   # F_{12}
  2  1     9.56579167e-01   # F_{21}
  2  2    -2.91472636e-01   # F_{22}
Block gauge Q= 9.09329852e+02  # SM gauge couplings
     1     3.62083238e-01   # g'(Q)MSSM DRbar
     2     6.41966863e-01   # g(Q)MSSM DRbar
     3     1.05866877e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.09329852e+02  
  3  3     8.55106946e-01   # Yt(Q)MSSM DRbar
Block yd Q= 9.09329852e+02  
  3  3     4.96923586e-01   # Yb(Q)MSSM DRbar
Block ye Q= 9.09329852e+02  
  3  3     4.23632990e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 9.09329852e+02 # Higgs mixing parameters
     1     6.92874459e+02    # mu(Q)MSSM DRbar
     2     3.92246833e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43959952e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.09296265e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 9.09329852e+02  # MSSM DRbar SUSY breaking parameters
     1     2.10820269e+02      # M_1(Q)
     2     3.91199545e+02      # M_2(Q)
     3     1.11446203e+03      # M_3(Q)
    21     6.23654822e+03      # mH1^2(Q)
    22    -4.71933989e+05      # mH2^2(Q)
    31     6.83025339e+02      # meL(Q)
    32     6.82846000e+02      # mmuL(Q)
    33     6.24231033e+02      # mtauL(Q)
    34     6.26405332e+02      # meR(Q)
    35     6.26011157e+02      # mmuR(Q)
    36     4.86861935e+02      # mtauR(Q)
    41     1.16778473e+03      # mqL1(Q)
    42     1.16774630e+03      # mqL2(Q)
    43     9.77869097e+02      # mqL3(Q)
    44     1.13585015e+03      # muR(Q)
    45     1.13584528e+03      # mcR(Q)
    46     8.39353000e+02      # mtR(Q)
    47     1.13206956e+03      # mdR(Q)
    48     1.13199427e+03      # msR(Q)
    49     1.01996475e+03      # mbR(Q)
Block au Q= 9.09329852e+02  
  1  1    -1.47786169e+03      # Au(Q)MSSM DRbar
  2  2    -1.47782545e+03      # Ac(Q)MSSM DRbar
  3  3    -1.01691193e+03      # At(Q)MSSM DRbar
Block ad Q= 9.09329852e+02  
  1  1    -1.71621014e+03      # Ad(Q)MSSM DRbar
  2  2    -1.71611784e+03      # As(Q)MSSM DRbar
  3  3    -1.43757703e+03      # Ab(Q)MSSM DRbar
Block ae Q= 9.09329852e+02  
  1  1    -6.22555558e+02      # Ae(Q)MSSM DRbar
  2  2    -6.22242078e+02      # Amu(Q)MSSM DRbar
  3  3    -5.19361503e+02      # Atau(Q)MSSM DRbar
