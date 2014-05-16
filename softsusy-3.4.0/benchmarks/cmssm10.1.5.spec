# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-3.62503300e+03
# dtanbeta/dmu=3.26163125e-02
# mu=8.27228868e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=4.44051861e-01
# dtanbeta/dm3sq=-9.93081060e-05
# m3sq= 8.23228935e+04
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
     1    1.75000000e+02   # m0
     2    7.00000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.67790022e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.86745030e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03934018e+01   # MW
        25     1.17040201e+02   # h0
        35     9.77200467e+02   # H0
        36     9.76997445e+02   # A0
        37     9.80509836e+02   # H+
   1000021     1.56472040e+03   # ~g
   1000022     2.91383739e+02   # ~neutralino(1)
   1000023     5.51268612e+02   # ~neutralino(2)
   1000024     5.51360734e+02   # ~chargino(1)
   1000025    -8.53787794e+02   # ~neutralino(3)
   1000035     8.65236861e+02   # ~neutralino(4)
   1000037     8.65432669e+02   # ~chargino(2)
   1000001     1.43572556e+03   # ~d_L
   1000002     1.43368235e+03   # ~u_L
   1000003     1.43572210e+03   # ~s_L
   1000004     1.43367889e+03   # ~c_L
   1000005     1.31647419e+03   # ~b_1
   1000006     1.10970143e+03   # ~t_1
   1000011     5.01771933e+02   # ~e_L
   1000012     4.95320226e+02   # ~nue_L
   1000013     5.01767461e+02   # ~mu_L
   1000014     4.95315699e+02   # ~numu_L
   1000015     3.11187549e+02   # ~stau_1
   1000016     4.93733038e+02   # ~nu_tau_L
   2000001     1.37264717e+03   # ~d_R
   2000002     1.37815906e+03   # ~u_R
   2000003     1.37264356e+03   # ~s_R
   2000004     1.37815535e+03   # ~c_R
   2000005     1.36710719e+03   # ~b_2
   2000006     1.35152498e+03   # ~t_2
   2000011     3.18456597e+02   # ~e_R
   2000013     3.18442305e+02   # ~mu_R
   2000015     5.01769953e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05491892e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97675633e-01   # N_{1,1}
  1  2    -9.67123647e-03   # N_{1,2}
  1  3     6.18910739e-02   # N_{1,3}
  1  4    -2.68196553e-02   # N_{1,4}
  2  1     2.25768096e-02   # N_{2,1}
  2  2     9.80061090e-01   # N_{2,2}
  2  3    -1.61744246e-01   # N_{2,3}
  2  4     1.13178382e-01   # N_{2,4}
  3  1    -2.43756963e-02   # N_{3,1}
  3  2     3.51992842e-02   # N_{3,2}
  3  3     7.05138296e-01   # N_{3,3}
  3  4     7.07775967e-01   # N_{3,4}
  4  1    -5.94932336e-02   # N_{4,1}
  4  2     1.95314457e-01   # N_{4,2}
  4  3     6.87596013e-01   # N_{4,3}
  4  4    -6.96795910e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.73075611e-01   # U_{1,1}
  1  2    -2.30486126e-01   # U_{1,2}
  2  1     2.30486126e-01   # U_{2,1}
  2  2     9.73075611e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.86845600e-01   # V_{1,1}
  1  2    -1.61665588e-01   # V_{1,2}
  2  1     1.61665588e-01   # V_{2,1}
  2  2     9.86845600e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.48680581e-01   # F_{11}
  1  2     9.37241619e-01   # F_{12}
  2  1     9.37241619e-01   # F_{21}
  2  2    -3.48680581e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.85072361e-01   # F_{11}
  1  2     1.72140767e-01   # F_{12}
  2  1    -1.72140767e-01   # F_{21}
  2  2     9.85072361e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.00666228e-01   # F_{11}
  1  2     9.94920253e-01   # F_{12}
  2  1     9.94920253e-01   # F_{21}
  2  2    -1.00666228e-01   # F_{22}
Block gauge Q= 1.18926641e+03  # SM gauge couplings
     1     3.63062041e-01   # g'(Q)MSSM DRbar
     2     6.41314549e-01   # g(Q)MSSM DRbar
     3     1.04483359e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.18926641e+03  
  3  3     8.50684685e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.18926641e+03  
  3  3     1.33202954e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.18926641e+03  
  3  3     1.00300029e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.18926641e+03 # Higgs mixing parameters
     1     8.48266627e+02    # mu(Q)MSSM DRbar
     2     9.63740313e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43752793e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     9.88726784e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.18926641e+03  # MSSM DRbar SUSY breaking parameters
     1     2.97409298e+02      # M_1(Q)
     2     5.47600235e+02      # M_2(Q)
     3     1.52266577e+03      # M_3(Q)
    21     2.13361241e+05      # mH1^2(Q)
    22    -6.91376607e+05      # mH2^2(Q)
    31     4.92904758e+02      # meL(Q)
    32     4.92900214e+02      # mmuL(Q)
    33     4.91524894e+02      # mtauL(Q)
    34     3.10275172e+02      # meR(Q)
    35     3.10260489e+02      # mmuR(Q)
    36     3.05789509e+02      # mtauR(Q)
    41     1.38742510e+03      # mqL1(Q)
    42     1.38742159e+03      # mqL2(Q)
    43     1.27994930e+03      # mqL3(Q)
    44     1.33325163e+03      # muR(Q)
    45     1.33324786e+03      # mcR(Q)
    46     1.09931685e+03      # mtR(Q)
    47     1.32659023e+03      # mdR(Q)
    48     1.32658657e+03      # msR(Q)
    49     1.31993129e+03      # mbR(Q)
Block au Q= 1.18926641e+03  
  1  1    -1.54178089e+03      # Au(Q)MSSM DRbar
  2  2    -1.54177407e+03      # Ac(Q)MSSM DRbar
  3  3    -1.19569030e+03      # At(Q)MSSM DRbar
Block ad Q= 1.18926641e+03  
  1  1    -1.87905475e+03      # Ad(Q)MSSM DRbar
  2  2    -1.87904845e+03      # As(Q)MSSM DRbar
  3  3    -1.75759887e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.18926641e+03  
  1  1    -4.12741912e+02      # Ae(Q)MSSM DRbar
  2  2    -4.12734646e+02      # Amu(Q)MSSM DRbar
  3  3    -4.10535591e+02      # Atau(Q)MSSM DRbar
