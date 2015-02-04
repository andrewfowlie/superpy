# SuperPy: Jacobian for naturalness priors.
# J = 1.21333401e+05
# b = 1.13710018e+05
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
     3    4.00000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    7.00000000e+02   # m0
     2    6.00000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.84442773e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=9.30313562e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03892150e+01   # MW
        25     1.18041012e+02   # h0
        35     8.07414348e+02   # H0
        36     8.07406064e+02   # A0
        37     8.11723029e+02   # H+
   1000021     1.38443709e+03   # ~g
   1000022     2.51959034e+02   # ~neutralino(1)
   1000023     4.79887074e+02   # ~neutralino(2)
   1000024     4.80022942e+02   # ~chargino(1)
   1000025    -8.00378721e+02   # ~neutralino(3)
   1000035     8.08686865e+02   # ~neutralino(4)
   1000037     8.09618659e+02   # ~chargino(2)
   1000001     1.41605868e+03   # ~d_L
   1000002     1.41397740e+03   # ~u_L
   1000003     1.41601728e+03   # ~s_L
   1000004     1.41393593e+03   # ~c_L
   1000005     1.17796254e+03   # ~b_1
   1000006     1.00406517e+03   # ~t_1
   1000011     8.05717973e+02   # ~e_L
   1000012     8.01474229e+02   # ~nue_L
   1000013     8.05521081e+02   # ~mu_L
   1000014     8.01276387e+02   # ~numu_L
   1000015     5.63303338e+02   # ~stau_1
   1000016     7.33671593e+02   # ~nu_tau_L
   2000001     1.37095587e+03   # ~d_R
   2000002     1.37464401e+03   # ~u_R
   2000003     1.37087449e+03   # ~s_R
   2000004     1.37463870e+03   # ~c_R
   2000005     1.25363443e+03   # ~b_2
   2000006     1.23007902e+03   # ~t_2
   2000011     7.35512161e+02   # ~e_R
   2000013     7.35076971e+02   # ~mu_R
   2000015     7.48885505e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.59483194e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97722076e-01   # N_{1,1}
  1  2    -7.90246293e-03   # N_{1,2}
  1  3     6.33733535e-02   # N_{1,3}
  1  4    -2.17262151e-02   # N_{1,4}
  2  1     1.94547134e-02   # N_{2,1}
  2  2     9.83847333e-01   # N_{2,2}
  2  3    -1.51729764e-01   # N_{2,3}
  2  4     9.29732102e-02   # N_{2,4}
  3  1    -2.90067456e-02   # N_{3,1}
  3  2     4.23533664e-02   # N_{3,2}
  3  3     7.04589784e-01   # N_{3,3}
  3  4     7.07755634e-01   # N_{3,4}
  4  1    -5.77129285e-02   # N_{4,1}
  4  2     1.73747427e-01   # N_{4,2}
  4  3     6.90300755e-01   # N_{4,3}
  4  4    -6.99975654e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.76351376e-01   # U_{1,1}
  1  2    -2.16189708e-01   # U_{1,2}
  2  1     2.16189708e-01   # U_{2,1}
  2  2     9.76351376e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.91117943e-01   # V_{1,1}
  1  2    -1.32985797e-01   # V_{1,2}
  2  1     1.32985797e-01   # V_{2,1}
  2  2     9.91117943e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.86654419e-01   # F_{11}
  1  2     9.22224680e-01   # F_{12}
  2  1     9.22224680e-01   # F_{21}
  2  2    -3.86654419e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.07908524e-01   # F_{11}
  1  2     4.19168357e-01   # F_{12}
  2  1    -4.19168357e-01   # F_{21}
  2  2     9.07908524e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     2.53748047e-01   # F_{11}
  1  2     9.67270349e-01   # F_{12}
  2  1     9.67270349e-01   # F_{21}
  2  2    -2.53748047e-01   # F_{22}
Block gauge Q= 1.07990028e+03  # SM gauge couplings
     1     3.62460504e-01   # g'(Q)MSSM DRbar
     2     6.41121252e-01   # g(Q)MSSM DRbar
     3     1.04982793e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.07990028e+03  
  3  3     8.49358380e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.07990028e+03  
  3  3     4.96030917e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.07990028e+03  
  3  3     4.23699153e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.07990028e+03 # Higgs mixing parameters
     1     7.95944676e+02    # mu(Q)MSSM DRbar
     2     3.91838294e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43758481e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     8.33216673e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.07990028e+03  # MSSM DRbar SUSY breaking parameters
     1     2.54876032e+02      # M_1(Q)
     2     4.70868871e+02      # M_2(Q)
     3     1.31978252e+03      # M_3(Q)
    21     3.00561352e+04      # mH1^2(Q)
    22    -6.19235193e+05      # mH2^2(Q)
    31     8.01543495e+02      # meL(Q)
    32     8.01345605e+02      # mmuL(Q)
    33     7.36149737e+02      # mtauL(Q)
    34     7.32531625e+02      # meR(Q)
    35     7.32095127e+02      # mmuR(Q)
    36     5.77375459e+02      # mtauR(Q)
    41     1.37567340e+03      # mqL1(Q)
    42     1.37563058e+03      # mqL2(Q)
    43     1.15995077e+03      # mqL3(Q)
    44     1.33677101e+03      # muR(Q)
    45     1.33676553e+03      # mcR(Q)
    46     9.99597971e+02      # mtR(Q)
    47     1.33211500e+03      # mdR(Q)
    48     1.33203103e+03      # msR(Q)
    49     1.20515188e+03      # mbR(Q)
Block au Q= 1.07990028e+03  
  1  1    -1.68092115e+03      # Au(Q)MSSM DRbar
  2  2    -1.68088128e+03      # Ac(Q)MSSM DRbar
  3  3    -1.17384139e+03      # At(Q)MSSM DRbar
Block ad Q= 1.07990028e+03  
  1  1    -1.94211331e+03      # Ad(Q)MSSM DRbar
  2  2    -1.94201178e+03      # As(Q)MSSM DRbar
  3  3    -1.63436373e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.07990028e+03  
  1  1    -6.62570212e+02      # Ae(Q)MSSM DRbar
  2  2    -6.62244128e+02      # Amu(Q)MSSM DRbar
  3  3    -5.54937920e+02      # Atau(Q)MSSM DRbar
