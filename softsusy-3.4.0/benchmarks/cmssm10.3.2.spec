# SuperPy: Jacobian for naturalness priors.
# J = inf
# b = 1.06783422e+05
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
     1    3.50000000e+02   # m0
     2    5.25000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.92117624e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=7.57989534e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03950358e+01   # MW
        25     1.15295769e+02   # h0
        35     8.15774276e+02   # H0
        36     8.15463635e+02   # A0
        37     8.19751780e+02   # H+
   1000021     1.20846384e+03   # ~g
   1000022     2.16027970e+02   # ~neutralino(1)
   1000023     4.08482376e+02   # ~neutralino(2)
   1000024     4.08509197e+02   # ~chargino(1)
   1000025    -6.63267136e+02   # ~neutralino(3)
   1000035     6.76696162e+02   # ~neutralino(4)
   1000037     6.76978112e+02   # ~chargino(2)
   1000001     1.14868222e+03   # ~d_L
   1000002     1.14609594e+03   # ~u_L
   1000003     1.14867924e+03   # ~s_L
   1000004     1.14609296e+03   # ~c_L
   1000005     1.04385461e+03   # ~b_1
   1000006     8.66905732e+02   # ~t_1
   1000011     4.97827611e+02   # ~e_L
   1000012     4.91298990e+02   # ~nue_L
   1000013     4.97822069e+02   # ~mu_L
   1000014     4.91293380e+02   # ~numu_L
   1000015     3.96930775e+02   # ~stau_1
   1000016     4.89476428e+02   # ~nu_tau_L
   2000001     1.10428597e+03   # ~d_R
   2000002     1.10771825e+03   # ~u_R
   2000003     1.10428289e+03   # ~s_R
   2000004     1.10771505e+03   # ~c_R
   2000005     1.09953143e+03   # ~b_2
   2000006     1.08349668e+03   # ~t_2
   2000011     4.03403460e+02   # ~e_R
   2000013     4.03389608e+02   # ~mu_R
   2000015     4.97677731e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.06063531e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.96136224e-01   # N_{1,1}
  1  2    -1.60424047e-02   # N_{1,2}
  1  3     7.96362354e-02   # N_{1,3}
  1  4    -3.33666638e-02   # N_{1,4}
  2  1     3.54784157e-02   # N_{2,1}
  2  2     9.71771743e-01   # N_{2,2}
  2  3    -1.93399450e-01   # N_{2,3}
  2  4     1.30374896e-01   # N_{2,4}
  3  1    -3.17958951e-02   # N_{3,1}
  3  2     4.62063488e-02   # N_{3,2}
  3  3     7.03764009e-01   # N_{3,3}
  3  4     7.08216220e-01   # N_{3,4}
  4  1    -7.37761890e-02   # N_{4,1}
  4  2     2.30797082e-01   # N_{4,2}
  4  3     6.78948409e-01   # N_{4,3}
  4  4    -6.93050387e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.60861608e-01   # U_{1,1}
  1  2    -2.77028827e-01   # U_{1,2}
  2  1     2.77028827e-01   # U_{2,1}
  2  2     9.60861608e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.82270991e-01   # V_{1,1}
  1  2    -1.87466528e-01   # V_{1,2}
  2  1     1.87466528e-01   # V_{2,1}
  2  2     9.82270991e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.89100845e-01   # F_{11}
  1  2     9.21195165e-01   # F_{12}
  2  1     9.21195165e-01   # F_{21}
  2  2    -3.89100845e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.88022730e-01   # F_{11}
  1  2     1.54308406e-01   # F_{12}
  2  1    -1.54308406e-01   # F_{21}
  2  2     9.88022730e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.34543671e-01   # F_{11}
  1  2     9.90907665e-01   # F_{12}
  2  1     9.90907665e-01   # F_{21}
  2  2    -1.34543671e-01   # F_{22}
Block gauge Q= 9.39845261e+02  # SM gauge couplings
     1     3.62388175e-01   # g'(Q)MSSM DRbar
     2     6.42474280e-01   # g(Q)MSSM DRbar
     3     1.05746413e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.39845261e+02  
  3  3     8.59241726e-01   # Yt(Q)MSSM DRbar
Block yd Q= 9.39845261e+02  
  3  3     1.35029964e-01   # Yb(Q)MSSM DRbar
Block ye Q= 9.39845261e+02  
  3  3     1.00371917e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 9.39845261e+02 # Higgs mixing parameters
     1     6.57541069e+02    # mu(Q)MSSM DRbar
     2     9.66452462e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44003865e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.87276072e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 9.39845261e+02  # MSSM DRbar SUSY breaking parameters
     1     2.20343593e+02      # M_1(Q)
     2     4.08443517e+02      # M_2(Q)
     3     1.16409256e+03      # M_3(Q)
    21     2.20897259e+05      # mH1^2(Q)
    22    -4.15901218e+05      # mH2^2(Q)
    31     4.92130365e+02      # meL(Q)
    32     4.92124745e+02      # mmuL(Q)
    33     4.90430369e+02      # mtauL(Q)
    34     3.98915631e+02      # meR(Q)
    35     3.98901627e+02      # mmuR(Q)
    36     3.94663099e+02      # mtauR(Q)
    41     1.11001882e+03      # mqL1(Q)
    42     1.11001578e+03      # mqL2(Q)
    43     1.01395820e+03      # mqL3(Q)
    44     1.07208768e+03      # muR(Q)
    45     1.07208442e+03      # mcR(Q)
    46     8.63039472e+02      # mtR(Q)
    47     1.06750516e+03      # mdR(Q)
    48     1.06750202e+03      # msR(Q)
    49     1.06172291e+03      # mbR(Q)
Block au Q= 9.39845261e+02  
  1  1    -1.18992613e+03      # Au(Q)MSSM DRbar
  2  2    -1.18992081e+03      # Ac(Q)MSSM DRbar
  3  3    -9.18657241e+02      # At(Q)MSSM DRbar
Block ad Q= 9.39845261e+02  
  1  1    -1.45511106e+03      # Ad(Q)MSSM DRbar
  2  2    -1.45510612e+03      # As(Q)MSSM DRbar
  3  3    -1.35986740e+03      # Ab(Q)MSSM DRbar
Block ae Q= 9.39845261e+02  
  1  1    -3.13507086e+02      # Ae(Q)MSSM DRbar
  2  2    -3.13501472e+02      # Amu(Q)MSSM DRbar
  3  3    -3.11810469e+02      # Atau(Q)MSSM DRbar
