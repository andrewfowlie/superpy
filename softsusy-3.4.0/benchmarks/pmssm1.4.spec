# SuperPy: Jacobian for naturalness priors.
# J = 4.08884361e+05
# b = 3.54215469e+05
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
     1     6.00000000e+02  # M_1(MX)
     2     2.50000000e+03  # M_2(MX)
     3     7.20000000e+02  # M_3(MX)
     11    0.00000000e+00  # At(MX)
     12    0.00000000e+00  # Ab(MX)
     13    0.00000000e+00  # Atau(MX)
     23    2.50000000e+03  # mu(MX)
     26    2.50000000e+03  # mA(pole)
     31    2.50000000e+03  # meL(MX)
     32    2.50000000e+03  # mmuL(MX)
     33    2.50000000e+03  # mtauL(MX)
     34    2.50000000e+03  # meR(MX)
     35    2.50000000e+03  # mmuR(MX)
     36    2.50000000e+03  # mtauR(MX)
     41    7.20000000e+02  # mqL1(MX)
     42    7.20000000e+02  # mqL2(MX)
     43    2.50000000e+03  # mqL3(MX)
     44    7.20000000e+02  # muR(MX)
     45    7.20000000e+02  # mcR(MX)
     46    2.50000000e+03  # mtR(MX)
     47    7.20000000e+02  # mdR(MX)
     48    7.20000000e+02  # msR(MX)
     49    2.50000000e+03  # mbR(MX)
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=7.86100708e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04032011e+01   # MW
        25     1.16632391e+02   # h0
        35     2.50009858e+03   # H0
        36     2.49999164e+03   # A0
        37     2.50164933e+03   # H+
   1000021     8.30991826e+02   # ~g
   1000022     5.91071784e+02   # ~neutralino(1)
   1000023     2.44521077e+03   # ~neutralino(2)
   1000024     2.44534814e+03   # ~chargino(1)
   1000025    -2.53398575e+03   # ~neutralino(3)
   1000035     2.57355045e+03   # ~neutralino(4)
   1000037     2.57363175e+03   # ~chargino(2)
   1000001     8.15371541e+02   # ~d_L
   1000002     8.11838163e+02   # ~u_L
   1000003     8.15371541e+02   # ~s_L
   1000004     8.11838163e+02   # ~c_L
   1000005     2.52132356e+03   # ~b_1
   1000006     2.52493548e+03   # ~t_1
   1000011     2.51893093e+03   # ~e_L
   1000012     2.51738994e+03   # ~nue_L
   1000013     2.51893093e+03   # ~mu_L
   1000014     2.51738994e+03   # ~numu_L
   1000015     2.49923910e+03   # ~stau_1
   1000016     2.51716819e+03   # ~nu_tau_L
   2000001     7.78980354e+02   # ~d_R
   2000002     7.75653449e+02   # ~u_R
   2000003     7.78980354e+02   # ~s_R
   2000004     7.75653449e+02   # ~c_R
   2000005     2.55087337e+03   # ~b_2
   2000006     2.55266414e+03   # ~t_2
   2000011     2.50366736e+03   # ~e_R
   2000013     2.50366736e+03   # ~mu_R
   2000015     2.52272970e+03   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.04712194e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.99801338e-01   # N_{1,1}
  1  2    -3.22074069e-04   # N_{1,2}
  1  3     1.89226325e-02   # N_{1,3}
  1  4    -6.25414722e-03   # N_{1,4}
  2  1     8.71619809e-03   # N_{2,1}
  2  2     8.83018030e-01   # N_{2,2}
  2  3    -3.37871994e-01   # N_{2,3}
  2  4     3.25646591e-01   # N_{2,4}
  3  1    -8.95268183e-03   # N_{3,1}
  3  2     9.78401644e-03   # N_{3,2}
  3  3     7.06920285e-01   # N_{3,3}
  3  4     7.07168886e-01   # N_{3,4}
  4  1    -1.55293751e-02   # N_{4,1}
  4  2     4.69236963e-01   # N_{4,2}
  4  3     6.21086275e-01   # N_{4,3}
  4  4    -6.27556651e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     8.29627681e-01   # U_{1,1}
  1  2    -5.58317034e-01   # U_{1,2}
  2  1     5.58317034e-01   # U_{2,1}
  2  2     8.29627681e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     8.40452070e-01   # V_{1,1}
  1  2    -5.41885891e-01   # V_{1,2}
  2  1     5.41885891e-01   # V_{2,1}
  2  2     8.40452070e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.59506199e-01   # F_{11}
  1  2     9.33142697e-01   # F_{12}
  2  1     9.33142697e-01   # F_{21}
  2  2    -3.59506199e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     3.97736205e-01   # F_{11}
  1  2     9.17499816e-01   # F_{12}
  2  1     9.17499816e-01   # F_{21}
  2  2    -3.97736205e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     4.08037834e-01   # F_{11}
  1  2     9.12965019e-01   # F_{12}
  2  1     9.12965019e-01   # F_{21}
  2  2    -4.08037834e-01   # F_{22}
Block gauge Q= 2.50363567e+03  # SM gauge couplings
     1     3.64599143e-01   # g'(Q)MSSM DRbar
     2     6.37149055e-01   # g(Q)MSSM DRbar
     3     1.04226634e+00   # g3(Q)MSSM DRbar
Block yu Q= 2.50363567e+03  
  3  3     8.30346658e-01   # Yt(Q)MSSM DRbar
Block yd Q= 2.50363567e+03  
  3  3     1.28666980e-01   # Yb(Q)MSSM DRbar
Block ye Q= 2.50363567e+03  
  3  3     9.96943961e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 2.50363567e+03 # Higgs mixing parameters
     1     2.53504991e+03    # mu(Q)MSSM DRbar
     2     9.54812252e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43771303e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.15715699e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 2.50363567e+03  # MSSM DRbar SUSY breaking parameters
     1     6.00000000e+02      # M_1(Q)
     2     2.50000000e+03      # M_2(Q)
     3     7.20000000e+02      # M_3(Q)
    21    -3.45287840e+05      # mH1^2(Q)
    22    -6.24770499e+06      # mH2^2(Q)
    31     2.50000000e+03      # meL(Q)
    32     2.50000000e+03      # mmuL(Q)
    33     2.50000000e+03      # mtauL(Q)
    34     2.50000000e+03      # meR(Q)
    35     2.50000000e+03      # mmuR(Q)
    36     2.50000000e+03      # mtauR(Q)
    41     7.19999998e+02      # mqL1(Q)
    42     7.19999998e+02      # mqL2(Q)
    43     2.50000000e+03      # mqL3(Q)
    44     7.19999998e+02      # muR(Q)
    45     7.19999998e+02      # mcR(Q)
    46     2.50000000e+03      # mtR(Q)
    47     7.19999998e+02      # mdR(Q)
    48     7.19999998e+02      # msR(Q)
    49     2.50000000e+03      # mbR(Q)
Block au Q= 2.50363567e+03  
  1  1     2.62246667e-06      # Au(Q)MSSM DRbar
  2  2     2.62248152e-06      # Ac(Q)MSSM DRbar
  3  3     3.86169834e-06      # At(Q)MSSM DRbar
Block ad Q= 2.50363567e+03  
  1  1     1.48617248e-06      # Ad(Q)MSSM DRbar
  2  2     1.48618480e-06      # As(Q)MSSM DRbar
  3  3     1.83664057e-06      # Ab(Q)MSSM DRbar
Block ae Q= 2.50363567e+03  
  1  1     0.00000000e+00      # Ae(Q)MSSM DRbar
  2  2     3.45922940e-08      # Amu(Q)MSSM DRbar
  3  3     3.45158408e-08      # Atau(Q)MSSM DRbar
