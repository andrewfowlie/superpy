# SuperPy: Jacobian for naturalness priors.
# J = 5.76933946e+05
# b = 4.26210973e+05
# Mu = 9.11876000e+01
# SOFTSUSY3.4.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.0       # version number
Block MODSEL  # Select model
     1    0   # nonUniversal
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
Block EXTPAR               # non-universal SUSY breaking parameters
     0    -1.00000000e+00  # Set MX=MSUSY
     1     7.50000000e+01  # M_1(MX)
     2     2.50000000e+03  # M_2(MX)
     3     2.50000000e+03  # M_3(MX)
     11    0.00000000e+00  # At(MX)
     12    0.00000000e+00  # Ab(MX)
     13    0.00000000e+00  # Atau(MX)
     23    2.50000000e+03  # mu(MX)
     26    2.50000000e+03  # mA(pole)
     31    1.00000000e+02  # meL(MX)
     32    1.00000000e+02  # mmuL(MX)
     33    2.50000000e+03  # mtauL(MX)
     34    1.00000000e+02  # meR(MX)
     35    1.00000000e+02  # mmuR(MX)
     36    2.50000000e+03  # mtauR(MX)
     41    2.50000000e+03  # mqL1(MX)
     42    2.50000000e+03  # mqL2(MX)
     43    2.50000000e+03  # mqL3(MX)
     44    2.50000000e+03  # muR(MX)
     45    2.50000000e+03  # mcR(MX)
     46    2.50000000e+03  # mtR(MX)
     47    2.50000000e+03  # mdR(MX)
     48    2.50000000e+03  # msR(MX)
     49    2.50000000e+03  # mbR(MX)
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=9.45792459e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04071188e+01   # MW
        25     1.16994146e+02   # h0
        35     2.50009591e+03   # H0
        36     2.49999179e+03   # A0
        37     2.50164491e+03   # H+
   1000021     2.64698268e+03   # ~g
   1000022     7.27758737e+01   # ~neutralino(1)
   1000023     2.45886238e+03   # ~neutralino(2)
   1000024     2.45904813e+03   # ~chargino(1)
   1000025    -2.53370976e+03   # ~neutralino(3)
   1000035     2.57991785e+03   # ~neutralino(4)
   1000037     2.58009371e+03   # ~chargino(2)
   1000001     2.60388992e+03   # ~d_L
   1000002     2.60285647e+03   # ~u_L
   1000003     2.60388992e+03   # ~s_L
   1000004     2.60285647e+03   # ~c_L
   1000005     2.58001945e+03   # ~b_1
   1000006     2.58165908e+03   # ~t_1
   1000011     2.39959487e+02   # ~e_L
   1000012     2.26865265e+02   # ~nue_L
   1000013     2.39959487e+02   # ~mu_L
   1000014     2.26865265e+02   # ~numu_L
   1000015     2.49898793e+03   # ~stau_1
   1000016     2.51704602e+03   # ~nu_tau_L
   2000001     2.58612155e+03   # ~d_R
   2000002     2.58540626e+03   # ~u_R
   2000003     2.58612155e+03   # ~s_R
   2000004     2.58540626e+03   # ~c_R
   2000005     2.60817059e+03   # ~b_2
   2000006     2.61253522e+03   # ~t_2
   2000011     1.30889660e+02   # ~e_R
   2000013     1.30889660e+02   # ~mu_R
   2000015     2.52258397e+03   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.04706433e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.99844849e-01   # N_{1,1}
  1  2    -4.94833270e-04   # N_{1,2}
  1  3     1.74529361e-02   # N_{1,3}
  1  4    -2.32967709e-03   # N_{1,4}
  2  1     7.44097429e-03   # N_{2,1}
  2  2     8.69037799e-01   # N_{2,2}
  2  3    -3.55756427e-01   # N_{2,3}
  2  4     3.43737254e-01   # N_{2,4}
  3  1    -1.06868463e-02   # N_{3,1}
  3  2     9.76352203e-03   # N_{3,2}
  3  3     7.06900303e-01   # N_{3,3}
  3  4     7.07165063e-01   # N_{3,4}
  4  1    -1.18617221e-02   # N_{4,1}
  4  2     4.94649100e-01   # N_{4,2}
  4  3     6.11084872e-01   # N_{4,3}
  4  4    -6.17864747e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     7.80493726e-01   # U_{1,1}
  1  2    -6.25163614e-01   # U_{1,2}
  2  1     6.25163614e-01   # U_{2,1}
  2  2     7.80493726e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     7.92622863e-01   # V_{1,1}
  1  2    -6.09712224e-01   # V_{1,2}
  2  1     6.09712224e-01   # V_{2,1}
  2  2     7.92622863e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.24829271e-01   # F_{11}
  1  2     9.05273489e-01   # F_{12}
  2  1     9.05273489e-01   # F_{21}
  2  2    -4.24829271e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     3.95959946e-01   # F_{11}
  1  2     9.18267783e-01   # F_{12}
  2  1     9.18267783e-01   # F_{21}
  2  2    -3.95959946e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     4.05971818e-01   # F_{11}
  1  2     9.13885596e-01   # F_{12}
  2  1     9.13885596e-01   # F_{21}
  2  2    -4.05971818e-01   # F_{22}
Block gauge Q= 2.50366726e+03  # SM gauge couplings
     1     3.64761070e-01   # g'(Q)MSSM DRbar
     2     6.36098330e-01   # g(Q)MSSM DRbar
     3     1.01384350e+00   # g3(Q)MSSM DRbar
Block yu Q= 2.50366726e+03  
  3  3     8.33687670e-01   # Yt(Q)MSSM DRbar
Block yd Q= 2.50366726e+03  
  3  3     1.25592748e-01   # Yb(Q)MSSM DRbar
Block ye Q= 2.50366726e+03  
  3  3     9.97753207e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 2.50366726e+03 # Higgs mixing parameters
     1     2.53486272e+03    # mu(Q)MSSM DRbar
     2     9.54823469e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43713089e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.15107085e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 2.50366726e+03  # MSSM DRbar SUSY breaking parameters
     1     7.50000000e+01      # M_1(Q)
     2     2.50000000e+03      # M_2(Q)
     3     2.50000000e+03      # M_3(Q)
    21    -3.41723281e+05      # mH1^2(Q)
    22    -6.24385011e+06      # mH2^2(Q)
    31     9.99999988e+01      # meL(Q)
    32     9.99999988e+01      # mmuL(Q)
    33     2.50000000e+03      # mtauL(Q)
    34     1.00000001e+02      # meR(Q)
    35     1.00000001e+02      # mmuR(Q)
    36     2.50000000e+03      # mtauR(Q)
    41     2.49999999e+03      # mqL1(Q)
    42     2.49999999e+03      # mqL2(Q)
    43     2.49999999e+03      # mqL3(Q)
    44     2.49999999e+03      # muR(Q)
    45     2.49999999e+03      # mcR(Q)
    46     2.49999999e+03      # mtR(Q)
    47     2.49999999e+03      # mdR(Q)
    48     2.49999999e+03      # msR(Q)
    49     2.49999999e+03      # mbR(Q)
Block au Q= 2.50366726e+03  
  1  1     6.10822209e-06      # Au(Q)MSSM DRbar
  2  2     6.10829474e-06      # Ac(Q)MSSM DRbar
  3  3     1.02279049e-05      # At(Q)MSSM DRbar
Block ad Q= 2.50366726e+03  
  1  1     2.09437893e-06      # Ad(Q)MSSM DRbar
  2  2     2.09445369e-06      # As(Q)MSSM DRbar
  3  3     3.43500190e-06      # Ab(Q)MSSM DRbar
Block ae Q= 2.50366726e+03  
  1  1     0.00000000e+00      # Ae(Q)MSSM DRbar
  2  2     1.21596833e-07      # Amu(Q)MSSM DRbar
  3  3     1.21904084e-07      # Atau(Q)MSSM DRbar
