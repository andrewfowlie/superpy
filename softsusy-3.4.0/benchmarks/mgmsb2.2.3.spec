# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.94129332e+03
# dtanbeta/dmu=4.02500392e-02
# mu=7.41144545e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=2.79156260e-01
# dtanbeta/dm3sq=-2.31179906e-04
# m3sq= 4.32705805e+04
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
     1    1.40000000e+05   # lambda
     2    1.00000000e+14   # M_mess
     5    1.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=7.48723509e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03968910e+01   # MW
        25     1.14910135e+02   # h0
        35     9.33705013e+02   # H0
        36     9.33658471e+02   # A0
        37     9.37316879e+02   # H+
   1000021     1.05931166e+03   # ~g
   1000022     1.86563482e+02   # ~neutralino(1)
   1000023     3.59542174e+02   # ~neutralino(2)
   1000024     3.59660285e+02   # ~chargino(1)
   1000025    -7.38363329e+02   # ~neutralino(3)
   1000035     7.45899288e+02   # ~neutralino(4)
   1000037     7.46761153e+02   # ~chargino(2)
   1000039     3.29430000e+00   # ~gravitino
   1000001     1.29687801e+03   # ~d_L
   1000002     1.29460127e+03   # ~u_L
   1000003     1.29687286e+03   # ~s_L
   1000004     1.29459611e+03   # ~c_L
   1000005     1.14632811e+03   # ~b_1
   1000006     9.11771434e+02   # ~t_1
   1000011     6.23595893e+02   # ~e_L
   1000012     6.18221863e+02   # ~nue_L
   1000013     6.23582673e+02   # ~mu_L
   1000014     6.18208537e+02   # ~numu_L
   1000015     4.08617193e+02   # ~stau_1
   1000016     6.13865083e+02   # ~nu_tau_L
   2000001     1.16705594e+03   # ~d_R
   2000002     1.19043640e+03   # ~u_R
   2000003     1.16704828e+03   # ~s_R
   2000004     1.19043263e+03   # ~c_R
   2000005     1.17522711e+03   # ~b_2
   2000006     1.18856938e+03   # ~t_2
   2000011     4.23406680e+02   # ~e_R
   2000013     4.23367460e+02   # ~mu_R
   2000015     6.20724509e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.02942275e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97453317e-01   # N_{1,1}
  1  2    -1.14287141e-02   # N_{1,2}
  1  3     6.70596140e-02   # N_{1,3}
  1  4    -2.14306751e-02   # N_{1,4}
  2  1     2.23959854e-02   # N_{2,1}
  2  2     9.86906494e-01   # N_{2,2}
  2  3    -1.40852390e-01   # N_{2,3}
  2  4     7.53299100e-02   # N_{2,4}
  3  1    -3.15831690e-02   # N_{3,1}
  3  2     4.71645655e-02   # N_{3,2}
  3  3     7.04035906e-01   # N_{3,3}
  3  4     7.07892259e-01   # N_{3,4}
  4  1    -5.98982841e-02   # N_{4,1}
  4  2     1.53819568e-01   # N_{4,2}
  4  3     6.92818198e-01   # N_{4,3}
  4  4    -7.01964871e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.79655921e-01   # U_{1,1}
  1  2    -2.00684521e-01   # U_{1,2}
  2  1     2.00684521e-01   # U_{2,1}
  2  2     9.79655921e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.94177796e-01   # V_{1,1}
  1  2    -1.07752075e-01   # V_{1,2}
  2  1     1.07752075e-01   # V_{2,1}
  2  2     9.94177796e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.11819052e-01   # F_{11}
  1  2     9.77308902e-01   # F_{12}
  2  1     9.77308902e-01   # F_{21}
  2  2    -2.11819052e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     4.58604623e-01   # F_{11}
  1  2     8.88640422e-01   # F_{12}
  2  1     8.88640422e-01   # F_{21}
  2  2    -4.58604623e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     9.01440634e-02   # F_{11}
  1  2     9.95928736e-01   # F_{12}
  2  1     9.95928736e-01   # F_{21}
  2  2    -9.01440634e-02   # F_{22}
Block gauge Q= 1.01069755e+03  # SM gauge couplings
     1     3.62544608e-01   # g'(Q)MSSM DRbar
     2     6.42385619e-01   # g(Q)MSSM DRbar
     3     1.05773395e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.01069755e+03  
  3  3     8.58251945e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.01069755e+03  
  3  3     1.98176381e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.01069755e+03  
  3  3     1.50812382e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.01069755e+03 # Higgs mixing parameters
     1     7.32111112e+02    # mu(Q)MSSM DRbar
     2     1.45008284e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43956251e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     9.10312729e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.01069755e+03  # MSSM DRbar SUSY breaking parameters
     1     1.89765299e+02      # M_1(Q)
     2     3.52259002e+02      # M_2(Q)
     3     9.72646547e+02      # M_3(Q)
    21     3.28158576e+05      # mH1^2(Q)
    22    -5.19248993e+05      # mH2^2(Q)
    31     6.20016381e+02      # meL(Q)
    32     6.20003089e+02      # mmuL(Q)
    33     6.15951620e+02      # mtauL(Q)
    34     4.18502625e+02      # meR(Q)
    35     4.18462971e+02      # mmuR(Q)
    36     4.06234082e+02      # mtauR(Q)
    41     1.26679570e+03      # mqL1(Q)
    42     1.26679040e+03      # mqL2(Q)
    43     1.14336534e+03      # mqL3(Q)
    44     1.16090475e+03      # muR(Q)
    45     1.16090085e+03      # mcR(Q)
    46     8.80885143e+02      # mtR(Q)
    47     1.13542195e+03      # mdR(Q)
    48     1.13541400e+03      # msR(Q)
    49     1.12090409e+03      # mbR(Q)
Block au Q= 1.01069755e+03  
  1  1    -9.06033522e+02      # Au(Q)MSSM DRbar
  2  2    -9.06028946e+02      # Ac(Q)MSSM DRbar
  3  3    -7.20223491e+02      # At(Q)MSSM DRbar
Block ad Q= 1.01069755e+03  
  1  1    -1.08312845e+03      # Ad(Q)MSSM DRbar
  2  2    -1.08312212e+03      # As(Q)MSSM DRbar
  3  3    -1.01316598e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.01069755e+03  
  1  1    -2.00105983e+02      # Ae(Q)MSSM DRbar
  2  2    -2.00099090e+02      # Amu(Q)MSSM DRbar
  3  3    -1.98003932e+02      # Atau(Q)MSSM DRbar
