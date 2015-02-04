# SuperPy: Jacobian for naturalness priors.
# J = 8.14580541e+04
# b = 1.60529126e+05
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
     1    9.50000000e+02   # m0
     2    4.50000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    2.12712602e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=2.87738757e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03964319e+01   # MW
        25     1.14836616e+02   # h0
        35     1.13080476e+03   # H0
        36     1.13075634e+03   # A0
        37     1.13370507e+03   # H+
   1000021     1.08804632e+03   # ~g
   1000022     1.85790043e+02   # ~neutralino(1)
   1000023     3.51522445e+02   # ~neutralino(2)
   1000024     3.51485226e+02   # ~chargino(1)
   1000025    -5.78814072e+02   # ~neutralino(3)
   1000035     5.94308914e+02   # ~neutralino(4)
   1000037     5.94507755e+02   # ~chargino(2)
   1000001     1.32508944e+03   # ~d_L
   1000002     1.32289103e+03   # ~u_L
   1000003     1.32508483e+03   # ~s_L
   1000004     1.32288642e+03   # ~c_L
   1000005     1.15058730e+03   # ~b_1
   1000006     9.19662154e+02   # ~t_1
   1000011     9.93641617e+02   # ~e_L
   1000012     9.90181452e+02   # ~nue_L
   1000013     9.93628121e+02   # ~mu_L
   1000014     9.90167914e+02   # ~numu_L
   1000015     9.55064135e+02   # ~stau_1
   1000016     9.86023237e+02   # ~nu_tau_L
   2000001     1.29958163e+03   # ~d_R
   2000002     1.30120236e+03   # ~u_R
   2000003     1.29957716e+03   # ~s_R
   2000004     1.30119739e+03   # ~c_R
   2000005     1.29082923e+03   # ~b_2
   2000006     1.17423386e+03   # ~t_2
   2000011     9.64495420e+02   # ~e_R
   2000013     9.64467503e+02   # ~mu_R
   2000015     9.90316954e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.04820187e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.94825098e-01   # N_{1,1}
  1  2    -2.08390912e-02   # N_{1,2}
  1  3     9.18749846e-02   # N_{1,3}
  1  4    -3.80492371e-02   # N_{1,4}
  2  1     4.56836494e-02   # N_{2,1}
  2  2     9.64568133e-01   # N_{2,2}
  2  3    -2.16153876e-01   # N_{2,3}
  2  4     1.44217968e-01   # N_{2,4}
  3  1    -3.66859265e-02   # N_{3,1}
  3  2     5.32860149e-02   # N_{3,2}
  3  3     7.02741105e-01   # N_{3,3}
  3  4     7.08498189e-01   # N_{3,4}
  4  1    -8.30070595e-02   # N_{4,1}
  4  2     2.57555138e-01   # N_{4,2}
  4  3     6.71558954e-01   # N_{4,3}
  4  4    -6.89770795e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.50506770e-01   # U_{1,1}
  1  2    -3.10703846e-01   # U_{1,2}
  2  1     3.10703846e-01   # U_{2,1}
  2  2     9.50506770e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.78081712e-01   # V_{1,1}
  1  2    -2.08221430e-01   # V_{1,2}
  2  1     2.08221430e-01   # V_{2,1}
  2  2     9.78081712e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.54743714e-01   # F_{11}
  1  2     9.67008604e-01   # F_{12}
  2  1     9.67008604e-01   # F_{21}
  2  2    -2.54743714e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.98928118e-01   # F_{11}
  1  2     4.62883874e-02   # F_{12}
  2  1    -4.62883874e-02   # F_{21}
  2  2     9.98928118e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.54482596e-01   # F_{11}
  1  2     9.87995510e-01   # F_{12}
  2  1     9.87995510e-01   # F_{21}
  2  2    -1.54482596e-01   # F_{22}
Block gauge Q= 1.01134642e+03  # SM gauge couplings
     1     3.62238815e-01   # g'(Q)MSSM DRbar
     2     6.42244256e-01   # g(Q)MSSM DRbar
     3     1.05665741e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.01134642e+03  
  3  3     8.60888043e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.01134642e+03  
  3  3     1.35946993e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.01134642e+03  
  3  3     9.95993287e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.01134642e+03 # Higgs mixing parameters
     1     5.71311033e+02    # mu(Q)MSSM DRbar
     2     9.65262214e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43744899e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.30542156e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.01134642e+03  # MSSM DRbar SUSY breaking parameters
     1     1.88858155e+02      # M_1(Q)
     2     3.50124766e+02      # M_2(Q)
     3     9.97548568e+02      # M_3(Q)
    21     9.38689987e+05      # mH1^2(Q)
    22    -2.98329407e+05      # mH2^2(Q)
    31     9.90764889e+02      # meL(Q)
    32     9.90751386e+02      # mmuL(Q)
    33     9.86734590e+02      # mtauL(Q)
    34     9.62267120e+02      # meR(Q)
    35     9.62239178e+02      # mmuR(Q)
    36     9.53907054e+02      # mtauR(Q)
    41     1.29502219e+03      # mqL1(Q)
    42     1.29501746e+03      # mqL2(Q)
    43     1.12437392e+03      # mqL3(Q)
    44     1.27399170e+03      # muR(Q)
    45     1.27398660e+03      # mcR(Q)
    46     8.98209226e+02      # mtR(Q)
    47     1.27154089e+03      # mdR(Q)
    48     1.27153628e+03      # msR(Q)
    49     1.26267545e+03      # mbR(Q)
Block au Q= 1.01134642e+03  
  1  1    -1.01828687e+03      # Au(Q)MSSM DRbar
  2  2    -1.01828232e+03      # Ac(Q)MSSM DRbar
  3  3    -7.84716396e+02      # At(Q)MSSM DRbar
Block ad Q= 1.01134642e+03  
  1  1    -1.24661982e+03      # Ad(Q)MSSM DRbar
  2  2    -1.24661561e+03      # As(Q)MSSM DRbar
  3  3    -1.16456459e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.01134642e+03  
  1  1    -2.68700618e+02      # Ae(Q)MSSM DRbar
  2  2    -2.68695807e+02      # Amu(Q)MSSM DRbar
  3  3    -2.67267234e+02      # Atau(Q)MSSM DRbar
