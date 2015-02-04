# SuperPy: Jacobian for naturalness priors.
# J = inf
# b = 1.20683278e+04
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
     1    3.50000000e+04   # lambda
     2    7.00000000e+04   # M_mess
     5    3.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=7.61727780e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04108408e+01   # MW
        25     1.10108118e+02   # h0
        35     3.39345359e+02   # H0
        36     3.38942713e+02   # A0
        37     3.48642950e+02   # H+
   1000021     8.39053578e+02   # ~g
   1000022     1.40798170e+02   # ~neutralino(1)
   1000023     2.32121627e+02   # ~neutralino(2)
   1000024     2.28419329e+02   # ~chargino(1)
   1000025    -2.81440057e+02   # ~neutralino(3)
   1000035     3.48501269e+02   # ~neutralino(4)
   1000037     3.48281199e+02   # ~chargino(2)
   1000039     5.80650000e-10   # ~gravitino
   1000001     7.99294026e+02   # ~d_L
   1000002     7.95446012e+02   # ~u_L
   1000003     7.99293018e+02   # ~s_L
   1000004     7.95444998e+02   # ~c_L
   1000005     7.59952922e+02   # ~b_1
   1000006     7.14432747e+02   # ~t_1
   1000011     2.31279539e+02   # ~e_L
   1000012     2.16884654e+02   # ~nue_L
   1000013     2.31278630e+02   # ~mu_L
   1000014     2.16883686e+02   # ~numu_L
   1000015     1.09009094e+02   # ~stau_1
   1000016     2.16446128e+02   # ~nu_tau_L
   2000001     7.68369036e+02   # ~d_R
   2000002     7.69254397e+02   # ~u_R
   2000003     7.68367626e+02   # ~s_R
   2000004     7.69253682e+02   # ~c_R
   2000005     7.74570074e+02   # ~b_2
   2000006     7.93598103e+02   # ~t_2
   2000011     1.16276744e+02   # ~e_R
   2000013     1.16273073e+02   # ~mu_R
   2000015     2.33641163e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -8.19974481e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.55294191e-01   # N_{1,1}
  1  2    -8.09318489e-02   # N_{1,2}
  1  3     2.47855174e-01   # N_{1,3}
  1  4    -1.39394611e-01   # N_{1,4}
  2  1    -2.62559216e-01   # N_{2,1}
  2  2    -6.53314769e-01   # N_{2,2}
  2  3     5.37964151e-01   # N_{2,3}
  2  4    -4.63505171e-01   # N_{2,4}
  3  1    -6.79272075e-02   # N_{3,1}
  3  2     9.21924623e-02   # N_{3,2}
  3  3     6.92496936e-01   # N_{3,3}
  3  4     7.12274132e-01   # N_{3,4}
  4  1    -1.17735134e-01   # N_{4,1}
  4  2     7.47081253e-01   # N_{4,2}
  4  3     4.11837807e-01   # N_{4,3}
  4  4    -5.08328300e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1    -5.81564299e-01   # U_{1,1}
  1  2     8.13500440e-01   # U_{1,2}
  2  1    -8.13500440e-01   # U_{2,1}
  2  2    -5.81564299e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1    -7.19857995e-01   # V_{1,1}
  1  2     6.94121364e-01   # V_{1,2}
  2  1    -6.94121364e-01   # V_{2,1}
  2  2    -7.19857995e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.95222434e-01   # F_{11}
  1  2     9.18585449e-01   # F_{12}
  2  1     9.18585449e-01   # F_{21}
  2  2    -3.95222434e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     5.99487715e-01   # F_{11}
  1  2     8.00383957e-01   # F_{12}
  2  1     8.00383957e-01   # F_{21}
  2  2    -5.99487715e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.73854954e-01   # F_{11}
  1  2     9.84771270e-01   # F_{12}
  2  1     9.84771270e-01   # F_{21}
  2  2    -1.73854954e-01   # F_{22}
Block gauge Q= 7.29473536e+02  # SM gauge couplings
     1     3.62409655e-01   # g'(Q)MSSM DRbar
     2     6.45744582e-01   # g(Q)MSSM DRbar
     3     1.07446075e+00   # g3(Q)MSSM DRbar
Block yu Q= 7.29473536e+02  
  3  3     8.72937134e-01   # Yt(Q)MSSM DRbar
Block yd Q= 7.29473536e+02  
  3  3     2.07499141e-01   # Yb(Q)MSSM DRbar
Block ye Q= 7.29473536e+02  
  3  3     1.52182160e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 7.29473536e+02 # Higgs mixing parameters
     1     2.71989952e+02    # mu(Q)MSSM DRbar
     2     1.45563671e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44072372e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.31665530e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 7.29473536e+02  # MSSM DRbar SUSY breaking parameters
     1     1.51185260e+02      # M_1(Q)
     2     2.86912197e+02      # M_2(Q)
     3     7.98376713e+02      # M_3(Q)
    21     4.26438153e+04      # mH1^2(Q)
    22    -6.75544397e+04      # mH2^2(Q)
    31     2.22469664e+02      # meL(Q)
    32     2.22468722e+02      # mmuL(Q)
    33     2.22179942e+02      # mtauL(Q)
    34     1.03663240e+02      # meR(Q)
    35     1.03659130e+02      # mmuR(Q)
    36     1.02392387e+02      # mtauR(Q)
    41     7.67260019e+02      # mqL1(Q)
    42     7.67258978e+02      # mqL2(Q)
    43     7.42570937e+02      # mqL3(Q)
    44     7.40596902e+02      # muR(Q)
    45     7.40596169e+02      # mcR(Q)
    46     6.90530914e+02      # mtR(Q)
    47     7.38232344e+02      # mdR(Q)
    48     7.38230896e+02      # msR(Q)
    49     7.35399742e+02      # mbR(Q)
Block au Q= 7.29473536e+02  
  1  1    -2.43801991e+02      # Au(Q)MSSM DRbar
  2  2    -2.43801644e+02      # Ac(Q)MSSM DRbar
  3  3    -2.29653916e+02      # At(Q)MSSM DRbar
Block ad Q= 7.29473536e+02  
  1  1    -2.59649048e+02      # Ad(Q)MSSM DRbar
  2  2    -2.59648562e+02      # As(Q)MSSM DRbar
  3  3    -2.54310929e+02      # Ab(Q)MSSM DRbar
Block ae Q= 7.29473536e+02  
  1  1    -2.38599126e+01      # Ae(Q)MSSM DRbar
  2  2    -2.38597384e+01      # Amu(Q)MSSM DRbar
  3  3    -2.38063737e+01      # Atau(Q)MSSM DRbar
