# SuperPy: Jacobian for naturalness priors.
# J = inf
# b = 3.41293034e+04
# Mu = 9.11876000e+01
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
     1    1.17000000e+05   # lambda
     2    1.30000000e+05   # M_mess
     5    1.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.00720872e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04003149e+01   # MW
        25     1.14214463e+02   # h0
        35     5.90783382e+02   # H0
        36     5.90516951e+02   # A0
        37     5.96220567e+02   # H+
   1000021     1.09537600e+03   # ~g
   1000022     1.90924072e+02   # ~neutralino(1)
   1000023     3.53353169e+02   # ~neutralino(2)
   1000024     3.52883759e+02   # ~chargino(1)
   1000025    -4.59626190e+02   # ~neutralino(3)
   1000035     4.91679266e+02   # ~neutralino(4)
   1000037     4.91554173e+02   # ~chargino(2)
   1000039     3.60477000e-09   # ~gravitino
   1000001     1.31061163e+03   # ~d_L
   1000002     1.30836935e+03   # ~u_L
   1000003     1.31060985e+03   # ~s_L
   1000004     1.30836757e+03   # ~c_L
   1000005     1.24250670e+03   # ~b_1
   1000006     1.15737163e+03   # ~t_1
   1000011     4.10082765e+02   # ~e_L
   1000012     4.02048609e+02   # ~nue_L
   1000013     4.10080995e+02   # ~mu_L
   1000014     4.02046805e+02   # ~numu_L
   1000015     1.98437503e+02   # ~stau_1
   1000016     4.01253215e+02   # ~nu_tau_L
   2000001     1.25097386e+03   # ~d_R
   2000002     1.25491737e+03   # ~u_R
   2000003     1.25097139e+03   # ~s_R
   2000004     1.25491610e+03   # ~c_R
   2000005     1.26202809e+03   # ~b_2
   2000006     1.27183996e+03   # ~t_2
   2000011     2.04125443e+02   # ~e_R
   2000013     2.04118252e+02   # ~mu_R
   2000015     4.10685897e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.33551108e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.89892316e-01   # N_{1,1}
  1  2    -2.78691699e-02   # N_{1,2}
  1  3     1.25570901e-01   # N_{1,3}
  1  4    -5.97365974e-02   # N_{1,4}
  2  1     8.98128282e-02   # N_{2,1}
  2  2     8.74544692e-01   # N_{2,2}
  2  3    -3.72523497e-01   # N_{2,3}
  2  4     2.97206127e-01   # N_{2,4}
  3  1    -4.43907702e-02   # N_{3,1}
  3  2     6.22240530e-02   # N_{3,2}
  3  3     7.01040334e-01   # N_{3,3}
  3  4     7.09013453e-01   # N_{3,4}
  4  1    -1.00380866e-01   # N_{4,1}
  4  2     4.80128169e-01   # N_{4,2}
  4  3     5.94979531e-01   # N_{4,3}
  4  4    -6.36710280e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     8.42114580e-01   # U_{1,1}
  1  2    -5.39298650e-01   # U_{1,2}
  2  1     5.39298650e-01   # U_{2,1}
  2  2     8.42114580e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.02318745e-01   # V_{1,1}
  1  2    -4.31069463e-01   # V_{1,2}
  2  1     4.31069463e-01   # V_{2,1}
  2  2     9.02318745e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.24833686e-01   # F_{11}
  1  2     9.74397154e-01   # F_{12}
  2  1     9.74397154e-01   # F_{21}
  2  2    -2.24833686e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     3.78493869e-01   # F_{11}
  1  2     9.25603798e-01   # F_{12}
  2  1     9.25603798e-01   # F_{21}
  2  2    -3.78493869e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     9.31178414e-02   # F_{11}
  1  2     9.95655095e-01   # F_{12}
  2  1     9.95655095e-01   # F_{21}
  2  2    -9.31178414e-02   # F_{22}
Block gauge Q= 1.18855422e+03  # SM gauge couplings
     1     3.63446808e-01   # g'(Q)MSSM DRbar
     2     6.43604183e-01   # g(Q)MSSM DRbar
     3     1.05287013e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.18855422e+03  
  3  3     8.58049592e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.18855422e+03  
  3  3     2.02200166e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.18855422e+03  
  3  3     1.51219730e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.18855422e+03 # Higgs mixing parameters
     1     4.51058240e+02    # mu(Q)MSSM DRbar
     2     1.44684313e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43363638e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     3.86170173e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.18855422e+03  # MSSM DRbar SUSY breaking parameters
     1     1.97798689e+02      # M_1(Q)
     2     3.70806087e+02      # M_2(Q)
     3     9.97189141e+02      # M_3(Q)
    21     1.43192806e+05      # mH1^2(Q)
    22    -1.75188883e+05      # mH2^2(Q)
    31     4.03313584e+02      # meL(Q)
    32     4.03311781e+02      # mmuL(Q)
    33     4.02759594e+02      # mtauL(Q)
    34     1.93654720e+02      # meR(Q)
    35     1.93647135e+02      # mmuR(Q)
    36     1.91311386e+02      # mtauR(Q)
    41     1.27538818e+03      # mqL1(Q)
    42     1.27538634e+03      # mqL2(Q)
    43     1.23139934e+03      # mqL3(Q)
    44     1.22177409e+03      # muR(Q)
    45     1.22177277e+03      # mcR(Q)
    46     1.13172751e+03      # mtR(Q)
    47     1.21657454e+03      # mdR(Q)
    48     1.21657197e+03      # msR(Q)
    49     1.21159957e+03      # mbR(Q)
Block au Q= 1.18855422e+03  
  1  1    -3.03358484e+02      # Au(Q)MSSM DRbar
  2  2    -3.03358055e+02      # Ac(Q)MSSM DRbar
  3  3    -2.85816733e+02      # At(Q)MSSM DRbar
Block ad Q= 1.18855422e+03  
  1  1    -3.22762951e+02      # Ad(Q)MSSM DRbar
  2  2    -3.22762354e+02      # As(Q)MSSM DRbar
  3  3    -3.16161514e+02      # Ab(Q)MSSM DRbar
Block ae Q= 1.18855422e+03  
  1  1    -3.17677014e+01      # Ae(Q)MSSM DRbar
  2  2    -3.17674668e+01      # Amu(Q)MSSM DRbar
  3  3    -3.16956619e+01      # Atau(Q)MSSM DRbar
