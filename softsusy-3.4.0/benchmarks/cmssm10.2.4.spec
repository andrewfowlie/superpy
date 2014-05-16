# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-3.39146965e+03
# dtanbeta/dmu=3.27695084e-02
# mu=7.75041734e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=4.43792460e-01
# dtanbeta/dm3sq=-1.07796703e-04
# m3sq= 7.40496900e+04
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
     1    2.75000000e+02   # m0
     2    6.50000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.75165623e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.44818526e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03935564e+01   # MW
        25     1.16614244e+02   # h0
        35     9.38838984e+02   # H0
        36     9.38617043e+02   # A0
        37     9.42287450e+02   # H+
   1000021     1.46434852e+03   # ~g
   1000022     2.70023877e+02   # ~neutralino(1)
   1000023     5.10929404e+02   # ~neutralino(2)
   1000024     5.11008653e+02   # ~chargino(1)
   1000025    -7.99937808e+02   # ~neutralino(3)
   1000035     8.11871574e+02   # ~neutralino(4)
   1000037     8.12084426e+02   # ~chargino(2)
   1000001     1.35934475e+03   # ~d_L
   1000002     1.35717957e+03   # ~u_L
   1000003     1.35934140e+03   # ~s_L
   1000004     1.35717621e+03   # ~c_L
   1000005     1.24287119e+03   # ~b_1
   1000006     1.04343957e+03   # ~t_1
   1000011     5.16397155e+02   # ~e_L
   1000012     5.10125576e+02   # ~nue_L
   1000013     5.16392079e+02   # ~mu_L
   1000014     5.10120440e+02   # ~numu_L
   1000015     3.62826904e+02   # ~stau_1
   1000016     5.08384511e+02   # ~nu_tau_L
   2000001     1.30193772e+03   # ~d_R
   2000002     1.30684072e+03   # ~u_R
   2000003     1.30193423e+03   # ~s_R
   2000004     1.30683712e+03   # ~c_R
   2000005     1.29649893e+03   # ~b_2
   2000006     1.27877514e+03   # ~t_2
   2000011     3.69700151e+02   # ~e_R
   2000013     3.69685762e+02   # ~mu_R
   2000015     5.16210930e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05586625e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97351625e-01   # N_{1,1}
  1  2    -1.10224764e-02   # N_{1,2}
  1  3     6.60466255e-02   # N_{1,3}
  1  4    -2.83916246e-02   # N_{1,4}
  2  1     2.53821714e-02   # N_{2,1}
  2  2     9.78163058e-01   # N_{2,2}
  2  3    -1.69511383e-01   # N_{2,3}
  2  4     1.17552830e-01   # N_{2,4}
  3  1    -2.61084388e-02   # N_{3,1}
  3  2     3.77553414e-02   # N_{3,2}
  3  3     7.04850583e-01   # N_{3,3}
  3  4     7.07869013e-01   # N_{3,4}
  4  1    -6.29589618e-02   # N_{4,1}
  4  2     2.04083489e-01   # N_{4,2}
  4  3     6.85630651e-01   # N_{4,3}
  4  4    -6.95914297e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.70262001e-01   # U_{1,1}
  1  2    -2.42057121e-01   # U_{1,2}
  2  1     2.42057121e-01   # U_{2,1}
  2  2     9.70262001e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.85721603e-01   # V_{1,1}
  1  2    -1.68383259e-01   # V_{1,2}
  2  1     1.68383259e-01   # V_{2,1}
  2  2     9.85721603e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.57314153e-01   # F_{11}
  1  2     9.33984259e-01   # F_{12}
  2  1     9.33984259e-01   # F_{21}
  2  2    -3.57314153e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.86875934e-01   # F_{11}
  1  2     1.61480308e-01   # F_{12}
  2  1    -1.61480308e-01   # F_{21}
  2  2     9.86875934e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.08381772e-01   # F_{11}
  1  2     9.94109346e-01   # F_{12}
  2  1     9.94109346e-01   # F_{21}
  2  2    -1.08381772e-01   # F_{22}
Block gauge Q= 1.12132635e+03  # SM gauge couplings
     1     3.62861324e-01   # g'(Q)MSSM DRbar
     2     6.41583006e-01   # g(Q)MSSM DRbar
     3     1.04797490e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.12132635e+03  
  3  3     8.52832115e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.12132635e+03  
  3  3     1.33666914e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.12132635e+03  
  3  3     1.00329685e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.12132635e+03 # Higgs mixing parameters
     1     7.94371474e+02    # mu(Q)MSSM DRbar
     2     9.64406373e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43812116e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     9.11704222e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.12132635e+03  # MSSM DRbar SUSY breaking parameters
     1     2.75277426e+02      # M_1(Q)
     2     5.07759777e+02      # M_2(Q)
     3     1.42089826e+03      # M_3(Q)
    21     2.31172278e+05      # mH1^2(Q)
    22    -6.06096706e+05      # mH2^2(Q)
    31     5.08677079e+02      # meL(Q)
    32     5.08671916e+02      # mmuL(Q)
    33     5.07110123e+02      # mtauL(Q)
    34     3.63435948e+02      # meR(Q)
    35     3.63421299e+02      # mmuR(Q)
    36     3.58968622e+02      # mtauR(Q)
    41     1.31370187e+03      # mqL1(Q)
    42     1.31369845e+03      # mqL2(Q)
    43     1.20807422e+03      # mqL3(Q)
    44     1.26444271e+03      # muR(Q)
    45     1.26443905e+03      # mcR(Q)
    46     1.03454091e+03      # mtR(Q)
    47     1.25841509e+03      # mdR(Q)
    48     1.25841154e+03      # msR(Q)
    49     1.25194139e+03      # mbR(Q)
Block au Q= 1.12132635e+03  
  1  1    -1.44229770e+03      # Au(Q)MSSM DRbar
  2  2    -1.44229130e+03      # Ac(Q)MSSM DRbar
  3  3    -1.11720349e+03      # At(Q)MSSM DRbar
Block ad Q= 1.12132635e+03  
  1  1    -1.75936325e+03      # Ad(Q)MSSM DRbar
  2  2    -1.75935733e+03      # As(Q)MSSM DRbar
  3  3    -1.64526125e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.12132635e+03  
  1  1    -3.84560744e+02      # Ae(Q)MSSM DRbar
  2  2    -3.84553943e+02      # Amu(Q)MSSM DRbar
  3  3    -3.82497655e+02      # Atau(Q)MSSM DRbar
