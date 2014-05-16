# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.55642809e+03
# dtanbeta/dmu=4.57789894e-02
# mu=6.45676397e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=2.79247849e-01
# dtanbeta/dm3sq=-3.09245296e-04
# m3sq= 3.01446421e+04
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
     1    1.20000000e+05   # lambda
     2    1.00000000e+14   # M_mess
     5    1.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=4.64761086e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03985751e+01   # MW
        25     1.13838254e+02   # h0
        35     8.07998333e+02   # H0
        36     8.07959256e+02   # A0
        37     8.12147235e+02   # H+
   1000021     9.21461851e+02   # ~g
   1000022     1.59245127e+02   # ~neutralino(1)
   1000023     3.06689193e+02   # ~neutralino(2)
   1000024     3.06784712e+02   # ~chargino(1)
   1000025    -6.43470999e+02   # ~neutralino(3)
   1000035     6.51881353e+02   # ~neutralino(4)
   1000037     6.52866133e+02   # ~chargino(2)
   1000039     2.82368571e+00   # ~gravitino
   1000001     1.12290984e+03   # ~d_L
   1000002     1.12025521e+03   # ~u_L
   1000003     1.12290536e+03   # ~s_L
   1000004     1.12025072e+03   # ~c_L
   1000005     9.91881263e+02   # ~b_1
   1000006     7.89270112e+02   # ~t_1
   1000011     5.37215736e+02   # ~e_L
   1000012     5.31012268e+02   # ~nue_L
   1000013     5.37204294e+02   # ~mu_L
   1000014     5.31000701e+02   # ~numu_L
   1000015     3.50840459e+02   # ~stau_1
   1000016     5.27233956e+02   # ~nu_tau_L
   2000001     1.01166212e+03   # ~d_R
   2000002     1.03130412e+03   # ~u_R
   2000003     1.01165546e+03   # ~s_R
   2000004     1.03130086e+03   # ~c_R
   2000005     1.01894619e+03   # ~b_2
   2000006     1.03463216e+03   # ~t_2
   2000011     3.64330030e+02   # ~e_R
   2000013     3.64296041e+02   # ~mu_R
   2000015     5.35189913e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.06305042e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.96602113e-01   # N_{1,1}
  1  2    -1.50955717e-02   # N_{1,2}
  1  3     7.72361128e-02   # N_{1,3}
  1  4    -2.43091571e-02   # N_{1,4}
  2  1     2.92496730e-02   # N_{2,1}
  2  2     9.83311101e-01   # N_{2,2}
  2  3    -1.58924804e-01   # N_{2,3}
  2  4     8.35861404e-02   # N_{2,4}
  3  1    -3.63849481e-02   # N_{3,1}
  3  2     5.45103398e-02   # N_{3,2}
  3  3     7.03020272e-01   # N_{3,3}
  3  4     7.08143527e-01   # N_{3,4}
  4  1    -6.78588293e-02   # N_{4,1}
  4  2     1.72916240e-01   # N_{4,2}
  4  3     6.88868628e-01   # N_{4,3}
  4  4    -7.00681930e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.73938101e-01   # U_{1,1}
  1  2    -2.26813966e-01   # U_{1,2}
  2  1     2.26813966e-01   # U_{2,1}
  2  2     9.73938101e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.92799573e-01   # V_{1,1}
  1  2    -1.19787344e-01   # V_{1,2}
  2  1     1.19787344e-01   # V_{2,1}
  2  2     9.92799573e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.42899516e-01   # F_{11}
  1  2     9.70051455e-01   # F_{12}
  2  1     9.70051455e-01   # F_{21}
  2  2    -2.42899516e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     5.09741386e-01   # F_{11}
  1  2     8.60327681e-01   # F_{12}
  2  1     8.60327681e-01   # F_{21}
  2  2    -5.09741386e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.05153018e-01   # F_{11}
  1  2     9.94456054e-01   # F_{12}
  2  1     9.94456054e-01   # F_{21}
  2  2    -1.05153018e-01   # F_{22}
Block gauge Q= 8.76516017e+02  # SM gauge couplings
     1     3.62244538e-01   # g'(Q)MSSM DRbar
     2     6.43150183e-01   # g(Q)MSSM DRbar
     3     1.06527454e+00   # g3(Q)MSSM DRbar
Block yu Q= 8.76516017e+02  
  3  3     8.63202471e-01   # Yt(Q)MSSM DRbar
Block yd Q= 8.76516017e+02  
  3  3     1.99575678e-01   # Yb(Q)MSSM DRbar
Block ye Q= 8.76516017e+02  
  3  3     1.50955756e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 8.76516017e+02 # Higgs mixing parameters
     1     6.37177140e+02    # mu(Q)MSSM DRbar
     2     1.45253959e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44133077e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.82096738e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 8.76516017e+02  # MSSM DRbar SUSY breaking parameters
     1     1.62306137e+02      # M_1(Q)
     2     3.02419766e+02      # M_2(Q)
     3     8.45276985e+02      # M_3(Q)
    21     2.42369729e+05      # mH1^2(Q)
    22    -3.94889938e+05      # mH2^2(Q)
    31     5.33709028e+02      # meL(Q)
    32     5.33697512e+02      # mmuL(Q)
    33     5.30193464e+02      # mtauL(Q)
    34     3.59426197e+02      # meR(Q)
    35     3.59391765e+02      # mmuR(Q)
    36     3.48789408e+02      # mtauR(Q)
    41     1.09610463e+03      # mqL1(Q)
    42     1.09610002e+03      # mqL2(Q)
    43     9.88866629e+02      # mqL3(Q)
    44     1.00535795e+03      # muR(Q)
    45     1.00535456e+03      # mcR(Q)
    46     7.62300720e+02      # mtR(Q)
    47     9.83660748e+02      # mdR(Q)
    48     9.83653829e+02      # msR(Q)
    49     9.71002108e+02      # mbR(Q)
Block au Q= 8.76516017e+02  
  1  1    -7.92290706e+02      # Au(Q)MSSM DRbar
  2  2    -7.92286675e+02      # Ac(Q)MSSM DRbar
  3  3    -6.28397950e+02      # At(Q)MSSM DRbar
Block ad Q= 8.76516017e+02  
  1  1    -9.48764818e+02      # Ad(Q)MSSM DRbar
  2  2    -9.48759229e+02      # As(Q)MSSM DRbar
  3  3    -8.87032607e+02      # Ab(Q)MSSM DRbar
Block ae Q= 8.76516017e+02  
  1  1    -1.73127897e+02      # Ae(Q)MSSM DRbar
  2  2    -1.73121879e+02      # Amu(Q)MSSM DRbar
  3  3    -1.71295545e+02      # Atau(Q)MSSM DRbar
