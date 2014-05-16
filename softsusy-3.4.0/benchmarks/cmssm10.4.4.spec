# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.65084833e+03
# dtanbeta/dmu=1.43734848e-02
# mu=6.11103996e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=4.40114730e-01
# dtanbeta/dm3sq=-6.17110736e-05
# m3sq= 1.45322868e+05
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
     1    1.05000000e+03   # m0
     2    5.00000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    2.05569689e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=3.56106574e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03955993e+01   # MW
        25     1.15542905e+02   # h0
        35     1.24696052e+03   # H0
        36     1.24690843e+03   # A0
        37     1.24960838e+03   # H+
   1000021     1.19892499e+03   # ~g
   1000022     2.07781375e+02   # ~neutralino(1)
   1000023     3.93581913e+02   # ~neutralino(2)
   1000024     3.93581505e+02   # ~chargino(1)
   1000025    -6.31593261e+02   # ~neutralino(3)
   1000035     6.46483078e+02   # ~neutralino(4)
   1000037     6.46624290e+02   # ~chargino(2)
   1000001     1.46189128e+03   # ~d_L
   1000002     1.45991330e+03   # ~u_L
   1000003     1.46188621e+03   # ~s_L
   1000004     1.45990822e+03   # ~c_L
   1000005     1.26992311e+03   # ~b_1
   1000006     1.01661732e+03   # ~t_1
   1000011     1.09845261e+03   # ~e_L
   1000012     1.09529475e+03   # ~nue_L
   1000013     1.09843777e+03   # ~mu_L
   1000014     1.09527987e+03   # ~numu_L
   1000015     1.05579452e+03   # ~stau_1
   1000016     1.09071843e+03   # ~nu_tau_L
   2000001     1.43350349e+03   # ~d_R
   2000002     1.43548547e+03   # ~u_R
   2000003     1.43349856e+03   # ~s_R
   2000004     1.43547999e+03   # ~c_R
   2000005     1.42384148e+03   # ~b_2
   2000006     1.29126513e+03   # ~t_2
   2000011     1.06599489e+03   # ~e_R
   2000013     1.06596418e+03   # ~mu_R
   2000015     1.09461150e+03   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.04674271e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.95660666e-01   # N_{1,1}
  1  2    -1.74526742e-02   # N_{1,2}
  1  3     8.42332514e-02   # N_{1,3}
  1  4    -3.54965129e-02   # N_{1,4}
  2  1     3.93254893e-02   # N_{2,1}
  2  2     9.67846033e-01   # N_{2,2}
  2  3    -2.05404206e-01   # N_{2,3}
  2  4     1.39773652e-01   # N_{2,4}
  3  1    -3.34166795e-02   # N_{3,1}
  3  2     4.83710544e-02   # N_{3,2}
  3  3     7.03487637e-01   # N_{3,3}
  3  4     7.08271637e-01   # N_{3,4}
  4  1    -7.74381710e-02   # N_{4,1}
  4  2     2.46230991e-01   # N_{4,2}
  4  3     6.75143701e-01   # N_{4,3}
  4  4    -6.91053263e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.55551792e-01   # U_{1,1}
  1  2    -2.94823292e-01   # U_{1,2}
  2  1     2.94823292e-01   # U_{2,1}
  2  2     9.55551792e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.79494068e-01   # V_{1,1}
  1  2    -2.01473004e-01   # V_{1,2}
  2  1     2.01473004e-01   # V_{2,1}
  2  2     9.79494068e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.33457897e-01   # F_{11}
  1  2     9.72366911e-01   # F_{12}
  2  1     9.72366911e-01   # F_{21}
  2  2    -2.33457897e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.99134915e-01   # F_{11}
  1  2     4.15863094e-02   # F_{12}
  2  1    -4.15863094e-02   # F_{21}
  2  2     9.99134915e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.38253151e-01   # F_{11}
  1  2     9.90396924e-01   # F_{12}
  2  1     9.90396924e-01   # F_{21}
  2  2    -1.38253151e-01   # F_{22}
Block gauge Q= 1.11563760e+03  # SM gauge couplings
     1     3.62448052e-01   # g'(Q)MSSM DRbar
     2     6.41730678e-01   # g(Q)MSSM DRbar
     3     1.05156080e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.11563760e+03  
  3  3     8.57554639e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.11563760e+03  
  3  3     1.35253472e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.11563760e+03  
  3  3     9.95273363e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.11563760e+03 # Higgs mixing parameters
     1     6.23972291e+02    # mu(Q)MSSM DRbar
     2     9.64101555e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43622675e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.58716966e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.11563760e+03  # MSSM DRbar SUSY breaking parameters
     1     2.10819942e+02      # M_1(Q)
     2     3.89854798e+02      # M_2(Q)
     3     1.10018517e+03      # M_3(Q)
    21     1.14794713e+06      # mH1^2(Q)
    22    -3.53858242e+05      # mH2^2(Q)
    31     1.09541001e+03      # meL(Q)
    32     1.09539516e+03      # mmuL(Q)
    33     1.09097365e+03      # mtauL(Q)
    34     1.06371481e+03      # meR(Q)
    35     1.06368408e+03      # mmuR(Q)
    36     1.05451095e+03      # mtauR(Q)
    41     1.42918342e+03      # mqL1(Q)
    42     1.42917822e+03      # mqL2(Q)
    43     1.24144275e+03      # mqL3(Q)
    44     1.40567542e+03      # muR(Q)
    45     1.40566979e+03      # mcR(Q)
    46     9.92264648e+02      # mtR(Q)
    47     1.40292120e+03      # mdR(Q)
    48     1.40291613e+03      # msR(Q)
    49     1.39319941e+03      # mbR(Q)
Block au Q= 1.11563760e+03  
  1  1    -1.11892469e+03      # Au(Q)MSSM DRbar
  2  2    -1.11891972e+03      # Ac(Q)MSSM DRbar
  3  3    -8.63683713e+02      # At(Q)MSSM DRbar
Block ad Q= 1.11563760e+03  
  1  1    -1.36814113e+03      # Ad(Q)MSSM DRbar
  2  2    -1.36813653e+03      # As(Q)MSSM DRbar
  3  3    -1.27848766e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.11563760e+03  
  1  1    -2.97143867e+02      # Ae(Q)MSSM DRbar
  2  2    -2.97138580e+02      # Amu(Q)MSSM DRbar
  3  3    -2.95567128e+02      # Atau(Q)MSSM DRbar
