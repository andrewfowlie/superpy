# SuperPy: Jacobian for naturalness priors.
# J = 4.20537471e+04
# b = 9.97723793e+04
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
     1    7.50000000e+02   # m0
     2    3.50000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    2.31215977e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=8.47372536e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03987782e+01   # MW
        25     1.13054910e+02   # h0
        35     8.96427693e+02   # H0
        36     8.96374827e+02   # A0
        37     9.00038226e+02   # H+
   1000021     8.63752386e+02   # ~g
   1000022     1.41997466e+02   # ~neutralino(1)
   1000023     2.67296466e+02   # ~neutralino(2)
   1000024     2.67129098e+02   # ~chargino(1)
   1000025    -4.69573440e+02   # ~neutralino(3)
   1000035     4.86704574e+02   # ~neutralino(4)
   1000037     4.87085127e+02   # ~chargino(2)
   1000001     1.04977992e+03   # ~d_L
   1000002     1.04696076e+03   # ~u_L
   1000003     1.04977625e+03   # ~s_L
   1000004     1.04695708e+03   # ~c_L
   1000005     9.10481640e+02   # ~b_1
   1000006     7.24391615e+02   # ~t_1
   1000011     7.84050827e+02   # ~e_L
   1000012     7.79739006e+02   # ~nue_L
   1000013     7.84040050e+02   # ~mu_L
   1000014     7.79728174e+02   # ~numu_L
   1000015     7.53606525e+02   # ~stau_1
   1000016     7.76422364e+02   # ~nu_tau_L
   2000001     1.02993613e+03   # ~d_R
   2000002     1.03078653e+03   # ~u_R
   2000003     1.02993255e+03   # ~s_R
   2000004     1.03078259e+03   # ~c_R
   2000005     1.02305678e+03   # ~b_2
   2000006     9.40653400e+02   # ~t_2
   2000011     7.61568958e+02   # ~e_R
   2000013     7.61546676e+02   # ~mu_R
   2000015     7.81845779e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05405352e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.91946006e-01   # N_{1,1}
  1  2    -3.20280593e-02   # N_{1,2}
  1  3     1.13882117e-01   # N_{1,3}
  1  4    -4.52569229e-02   # N_{1,4}
  2  1     6.60357124e-02   # N_{2,1}
  2  2     9.54711705e-01   # N_{2,2}
  2  3    -2.44849980e-01   # N_{2,3}
  2  4     1.55606336e-01   # N_{2,4}
  3  1    -4.58679717e-02   # N_{3,1}
  3  2     6.71693159e-02   # N_{3,2}
  3  3     7.00264730e-01   # N_{3,3}
  3  4     7.09234601e-01   # N_{3,4}
  4  1    -9.78699956e-02   # N_{4,1}
  4  2     2.88076459e-01   # N_{4,2}
  4  3     6.60839359e-01   # N_{4,3}
  4  4    -6.86093842e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.35275225e-01   # U_{1,1}
  1  2    -3.53921254e-01   # U_{1,2}
  2  1     3.53921254e-01   # U_{2,1}
  2  2     9.35275225e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.74109952e-01   # V_{1,1}
  1  2    -2.26074770e-01   # V_{1,2}
  2  1     2.26074770e-01   # V_{2,1}
  2  2     9.74109952e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.08912664e-01   # F_{11}
  1  2     9.51090409e-01   # F_{12}
  2  1     9.51090409e-01   # F_{21}
  2  2    -3.08912664e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.98229894e-01   # F_{11}
  1  2     5.94733363e-02   # F_{12}
  2  1    -5.94733363e-02   # F_{21}
  2  2     9.98229894e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.99325302e-01   # F_{11}
  1  2     9.79933377e-01   # F_{12}
  2  1     9.79933377e-01   # F_{21}
  2  2    -1.99325302e-01   # F_{22}
Block gauge Q= 8.02091845e+02  # SM gauge couplings
     1     3.61745531e-01   # g'(Q)MSSM DRbar
     2     6.43484750e-01   # g(Q)MSSM DRbar
     3     1.06904006e+00   # g3(Q)MSSM DRbar
Block yu Q= 8.02091845e+02  
  3  3     8.69011348e-01   # Yt(Q)MSSM DRbar
Block yd Q= 8.02091845e+02  
  3  3     1.37635997e-01   # Yb(Q)MSSM DRbar
Block ye Q= 8.02091845e+02  
  3  3     9.97653216e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 8.02091845e+02 # Higgs mixing parameters
     1     4.61996583e+02    # mu(Q)MSSM DRbar
     2     9.68063933e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44043887e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     8.20640673e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 8.02091845e+02  # MSSM DRbar SUSY breaking parameters
     1     1.45252573e+02      # M_1(Q)
     2     2.70916262e+02      # M_2(Q)
     3     7.89951149e+02      # M_3(Q)
    21     5.83323823e+05      # mH1^2(Q)
    22    -1.98095990e+05      # mH2^2(Q)
    31     7.81417672e+02      # meL(Q)
    32     7.81406882e+02      # mmuL(Q)
    33     7.78205536e+02      # mtauL(Q)
    34     7.59369333e+02      # meR(Q)
    35     7.59347018e+02      # mmuR(Q)
    36     7.52709963e+02      # mtauR(Q)
    41     1.02501557e+03      # mqL1(Q)
    42     1.02501179e+03      # mqL2(Q)
    43     8.88829654e+02      # mqL3(Q)
    44     1.00893184e+03      # muR(Q)
    45     1.00892778e+03      # mcR(Q)
    46     7.09011823e+02      # mtR(Q)
    47     1.00708183e+03      # mdR(Q)
    48     1.00707814e+03      # msR(Q)
    49     9.99953421e+02      # mbR(Q)
Block au Q= 8.02091845e+02  
  1  1    -8.13408638e+02      # Au(Q)MSSM DRbar
  2  2    -8.13404969e+02      # Ac(Q)MSSM DRbar
  3  3    -6.24368812e+02      # At(Q)MSSM DRbar
Block ad Q= 8.02091845e+02  
  1  1    -9.98740059e+02      # Ad(Q)MSSM DRbar
  2  2    -9.98736653e+02      # As(Q)MSSM DRbar
  3  3    -9.32302526e+02      # Ab(Q)MSSM DRbar
Block ae Q= 8.02091845e+02  
  1  1    -2.11356060e+02      # Ae(Q)MSSM DRbar
  2  2    -2.11352218e+02      # Amu(Q)MSSM DRbar
  3  3    -2.10214432e+02      # Atau(Q)MSSM DRbar
