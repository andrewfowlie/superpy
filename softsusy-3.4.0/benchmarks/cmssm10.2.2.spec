# SuperPy: Jacobian for naturalness priors.
# J = 4.23674343e+04
# b = 1.08397417e+05
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
     1    2.25000000e+02   # m0
     2    5.50000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.84678066e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=8.70325870e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03950919e+01   # MW
        25     1.15549973e+02   # h0
        35     8.03087862e+02   # H0
        36     8.02826535e+02   # A0
        37     8.07093716e+02   # H+
   1000021     1.25474379e+03   # ~g
   1000022     2.26322393e+02   # ~neutralino(1)
   1000023     4.28040730e+02   # ~neutralino(2)
   1000024     4.28082075e+02   # ~chargino(1)
   1000025    -6.90977640e+02   # ~neutralino(3)
   1000035     7.03991949e+02   # ~neutralino(4)
   1000037     7.04268105e+02   # ~chargino(2)
   1000001     1.16475381e+03   # ~d_L
   1000002     1.16219940e+03   # ~u_L
   1000003     1.16475092e+03   # ~s_L
   1000004     1.16219651e+03   # ~c_L
   1000005     1.06445356e+03   # ~b_1
   1000006     8.88690398e+02   # ~t_1
   1000011     4.34554129e+02   # ~e_L
   1000012     4.27134176e+02   # ~nue_L
   1000013     4.34549863e+02   # ~mu_L
   1000014     4.27129838e+02   # ~numu_L
   1000015     3.01537708e+02   # ~stau_1
   1000016     4.25662316e+02   # ~nu_tau_L
   2000001     1.11633696e+03   # ~d_R
   2000002     1.12016846e+03   # ~u_R
   2000003     1.11633395e+03   # ~s_R
   2000004     1.12016537e+03   # ~c_R
   2000005     1.11207648e+03   # ~b_2
   2000006     1.10554760e+03   # ~t_2
   2000011     3.08208830e+02   # ~e_R
   2000013     3.08196618e+02   # ~mu_R
   2000015     4.34977865e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.06183157e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.96456132e-01   # N_{1,1}
  1  2    -1.47835508e-02   # N_{1,2}
  1  3     7.63244103e-02   # N_{1,3}
  1  4    -3.21124300e-02   # N_{1,4}
  2  1     3.28907564e-02   # N_{2,1}
  2  2     9.73482576e-01   # N_{2,2}
  2  3    -1.87432032e-01   # N_{2,3}
  2  4     1.26961037e-01   # N_{2,4}
  3  1    -3.04492123e-02   # N_{3,1}
  3  2     4.42276178e-02   # N_{3,2}
  3  3     7.04034456e-01   # N_{3,3}
  3  4     7.08132931e-01   # N_{3,4}
  4  1    -7.11773893e-02   # N_{4,1}
  4  2     2.23957670e-01   # N_{4,2}
  4  3     6.80719694e-01   # N_{4,3}
  4  4    -6.93828105e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.63400272e-01   # U_{1,1}
  1  2    -2.68066999e-01   # U_{1,2}
  2  1     2.68066999e-01   # U_{2,1}
  2  2     9.63400272e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.83264861e-01   # V_{1,1}
  1  2    -1.82181812e-01   # V_{1,2}
  2  1     1.82181812e-01   # V_{2,1}
  2  2     9.83264861e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.97151363e-01   # F_{11}
  1  2     9.17753123e-01   # F_{12}
  2  1     9.17753123e-01   # F_{21}
  2  2    -3.97151363e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.82673477e-01   # F_{11}
  1  2     1.85345187e-01   # F_{12}
  2  1    -1.85345187e-01   # F_{21}
  2  2     9.82673477e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.28700828e-01   # F_{11}
  1  2     9.91683466e-01   # F_{12}
  2  1     9.91683466e-01   # F_{21}
  2  2    -1.28700828e-01   # F_{22}
Block gauge Q= 9.61499316e+02  # SM gauge couplings
     1     3.62546415e-01   # g'(Q)MSSM DRbar
     2     6.42427803e-01   # g(Q)MSSM DRbar
     3     1.05598973e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.61499316e+02  
  3  3     8.58048315e-01   # Yt(Q)MSSM DRbar
Block yd Q= 9.61499316e+02  
  3  3     1.34717940e-01   # Yb(Q)MSSM DRbar
Block ye Q= 9.61499316e+02  
  3  3     1.00441673e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 9.61499316e+02 # Higgs mixing parameters
     1     6.85470319e+02    # mu(Q)MSSM DRbar
     2     9.66211842e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43991018e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.67570075e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 9.61499316e+02  # MSSM DRbar SUSY breaking parameters
     1     2.31222590e+02      # M_1(Q)
     2     4.28205323e+02      # M_2(Q)
     3     1.21637152e+03      # M_3(Q)
    21     1.62290542e+05      # mH1^2(Q)
    22    -4.53186378e+05      # mH2^2(Q)
    31     4.27415559e+02      # meL(Q)
    32     4.27411212e+02      # mmuL(Q)
    33     4.26098711e+02      # mtauL(Q)
    34     3.01932424e+02      # meR(Q)
    35     3.01919948e+02      # mmuR(Q)
    36     2.98133766e+02      # mtauR(Q)
    41     1.12481307e+03      # mqL1(Q)
    42     1.12481012e+03      # mqL2(Q)
    43     1.03424045e+03      # mqL3(Q)
    44     1.08341273e+03      # muR(Q)
    45     1.08340959e+03      # mcR(Q)
    46     8.86461888e+02      # mtR(Q)
    47     1.07838690e+03      # mdR(Q)
    48     1.07838385e+03      # msR(Q)
    49     1.07279141e+03      # mbR(Q)
Block au Q= 9.61499316e+02  
  1  1    -1.24195615e+03      # Au(Q)MSSM DRbar
  2  2    -1.24195060e+03      # Ac(Q)MSSM DRbar
  3  3    -9.59582065e+02      # At(Q)MSSM DRbar
Block ad Q= 9.61499316e+02  
  1  1    -1.51788828e+03      # Ad(Q)MSSM DRbar
  2  2    -1.51788314e+03      # As(Q)MSSM DRbar
  3  3    -1.41875738e+03      # Ab(Q)MSSM DRbar
Block ae Q= 9.61499316e+02  
  1  1    -3.27865231e+02      # Ae(Q)MSSM DRbar
  2  2    -3.27859374e+02      # Amu(Q)MSSM DRbar
  3  3    -3.26091751e+02      # Atau(Q)MSSM DRbar
