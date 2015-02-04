# SuperPy: Jacobian for naturalness priors.
# J = inf
# b = 2.12834703e+05
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
     1    5.00000000e+02   # m0
     2    7.50000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.71306925e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=2.32293069e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03926120e+01   # MW
        25     1.17524018e+02   # h0
        35     1.13477705e+03   # H0
        36     1.13455838e+03   # A0
        37     1.13767487e+03   # H+
   1000021     1.68018660e+03   # ~g
   1000022     3.14733890e+02   # ~neutralino(1)
   1000023     5.95355213e+02   # ~neutralino(2)
   1000024     5.95456108e+02   # ~chargino(1)
   1000025    -9.05522884e+02   # ~neutralino(3)
   1000035     9.16743926e+02   # ~neutralino(4)
   1000037     9.16901657e+02   # ~chargino(2)
   1000001     1.59622561e+03   # ~d_L
   1000002     1.59440890e+03   # ~u_L
   1000003     1.59622151e+03   # ~s_L
   1000004     1.59440480e+03   # ~c_L
   1000005     1.45220829e+03   # ~b_1
   1000006     1.21849113e+03   # ~t_1
   1000011     7.07467656e+02   # ~e_L
   1000012     7.02805872e+02   # ~nue_L
   1000013     7.07459895e+02   # ~mu_L
   1000014     7.02798064e+02   # ~numu_L
   1000015     5.66700585e+02   # ~stau_1
   1000016     7.00255947e+02   # ~nu_tau_L
   2000001     1.53218980e+03   # ~d_R
   2000002     1.53794160e+03   # ~u_R
   2000003     1.53218559e+03   # ~s_R
   2000004     1.53793717e+03   # ~c_R
   2000005     1.52483403e+03   # ~b_2
   2000006     1.48162781e+03   # ~t_2
   2000011     5.74251796e+02   # ~e_R
   2000013     5.74232432e+02   # ~mu_R
   2000015     7.06027281e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05087066e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97909245e-01   # N_{1,1}
  1  2    -8.60378854e-03   # N_{1,2}
  1  3     5.86533203e-02   # N_{1,3}
  1  4    -2.57468586e-02   # N_{1,4}
  2  1     2.05330028e-02   # N_{2,1}
  2  2     9.81172503e-01   # N_{2,2}
  2  3    -1.56657233e-01   # N_{2,3}
  2  4     1.11073968e-01   # N_{2,4}
  3  1    -2.29156955e-02   # N_{3,1}
  3  2     3.29869743e-02   # N_{3,2}
  3  3     7.05374360e-01   # N_{3,3}
  3  4     7.07696081e-01   # N_{3,4}
  4  1    -5.68366482e-02   # N_{4,1}
  4  2     1.90100904e-01   # N_{4,2}
  4  3     6.88814425e-01   # N_{4,3}
  4  4    -6.97256000e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.74641699e-01   # U_{1,1}
  1  2    -2.23771218e-01   # U_{1,2}
  2  1     2.23771218e-01   # U_{2,1}
  2  2     9.74641699e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.87245548e-01   # V_{1,1}
  1  2    -1.59204986e-01   # V_{1,2}
  2  1     1.59204986e-01   # V_{2,1}
  2  2     9.87245548e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.05634973e-01   # F_{11}
  1  2     9.52148761e-01   # F_{12}
  2  1     9.52148761e-01   # F_{21}
  2  2    -3.05634973e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.93474572e-01   # F_{11}
  1  2     1.14053821e-01   # F_{12}
  2  1    -1.14053821e-01   # F_{21}
  2  2     9.93474572e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     9.29650919e-02   # F_{11}
  1  2     9.95669369e-01   # F_{12}
  2  1     9.95669369e-01   # F_{21}
  2  2    -9.29650919e-02   # F_{22}
Block gauge Q= 1.30491889e+03  # SM gauge couplings
     1     3.63073611e-01   # g'(Q)MSSM DRbar
     2     6.40701402e-01   # g(Q)MSSM DRbar
     3     1.04051838e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.30491889e+03  
  3  3     8.48204199e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.30491889e+03  
  3  3     1.32795009e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.30491889e+03  
  3  3     1.00129179e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.30491889e+03 # Higgs mixing parameters
     1     8.99472361e+02    # mu(Q)MSSM DRbar
     2     9.62629793e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43627041e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.32845513e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.30491889e+03  # MSSM DRbar SUSY breaking parameters
     1     3.19733699e+02      # M_1(Q)
     2     5.87648135e+02      # M_2(Q)
     3     1.62241669e+03      # M_3(Q)
    21     4.50344100e+05      # mH1^2(Q)
    22    -7.71702686e+05      # mH2^2(Q)
    31     7.00513033e+02      # meL(Q)
    32     7.00505177e+02      # mmuL(Q)
    33     6.98127723e+02      # mtauL(Q)
    34     5.69545975e+02      # meR(Q)
    35     5.69526456e+02      # mmuR(Q)
    36     5.63597957e+02      # mtauR(Q)
    41     1.54467524e+03      # mqL1(Q)
    42     1.54467105e+03      # mqL2(Q)
    43     1.41196913e+03      # mqL3(Q)
    44     1.48960733e+03      # muR(Q)
    45     1.48960282e+03      # mcR(Q)
    46     1.20031574e+03      # mtR(Q)
    47     1.48283993e+03      # mdR(Q)
    48     1.48283564e+03      # msR(Q)
    49     1.47496872e+03      # mbR(Q)
Block au Q= 1.30491889e+03  
  1  1    -1.63774566e+03      # Au(Q)MSSM DRbar
  2  2    -1.63773844e+03      # Ac(Q)MSSM DRbar
  3  3    -1.27125754e+03      # At(Q)MSSM DRbar
Block ad Q= 1.30491889e+03  
  1  1    -1.99454448e+03      # Ad(Q)MSSM DRbar
  2  2    -1.99453782e+03      # As(Q)MSSM DRbar
  3  3    -1.86593339e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.30491889e+03  
  1  1    -4.40695499e+02      # Ae(Q)MSSM DRbar
  2  2    -4.40687774e+02      # Amu(Q)MSSM DRbar
  3  3    -4.38352544e+02      # Atau(Q)MSSM DRbar
