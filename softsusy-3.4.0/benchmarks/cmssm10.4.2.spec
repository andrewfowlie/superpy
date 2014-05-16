# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-2.19291541e+03
# dtanbeta/dmu=1.74938366e-02
# mu=5.07219789e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=4.39648877e-01
# dtanbeta/dm3sq=-9.36610273e-05
# m3sq= 9.26376279e+04
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
     1    8.50000000e+02   # m0
     2    4.00000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    2.21120337e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=2.24018433e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03976191e+01   # MW
        25     1.14019122e+02   # h0
        35     1.01399606e+03   # H0
        36     1.01394774e+03   # A0
        37     1.01720929e+03   # H+
   1000021     9.76358671e+02   # ~g
   1000022     1.63864696e+02   # ~neutralino(1)
   1000023     3.09432595e+02   # ~neutralino(2)
   1000024     3.09343067e+02   # ~chargino(1)
   1000025    -5.24850006e+02   # ~neutralino(3)
   1000035     5.41079153e+02   # ~neutralino(4)
   1000037     5.41354611e+02   # ~chargino(2)
   1000001     1.18774206e+03   # ~d_L
   1000002     1.18527047e+03   # ~u_L
   1000003     1.18773792e+03   # ~s_L
   1000004     1.18526632e+03   # ~c_L
   1000005     1.03079331e+03   # ~b_1
   1000006     8.22265547e+02   # ~t_1
   1000011     8.88837880e+02   # ~e_L
   1000012     8.85002903e+02   # ~nue_L
   1000013     8.88825739e+02   # ~mu_L
   1000014     8.84990714e+02   # ~numu_L
   1000015     8.54334649e+02   # ~stau_1
   1000016     8.81264561e+02   # ~nu_tau_L
   2000001     1.16508601e+03   # ~d_R
   2000002     1.16633159e+03   # ~u_R
   2000003     1.16508197e+03   # ~s_R
   2000004     1.16632713e+03   # ~c_R
   2000005     1.15725959e+03   # ~b_2
   2000006     1.05751618e+03   # ~t_2
   2000011     8.63017087e+02   # ~e_R
   2000013     8.62991978e+02   # ~mu_R
   2000015     8.86055611e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05047197e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.93659564e-01   # N_{1,1}
  1  2    -2.54589762e-02   # N_{1,2}
  1  3     1.01459066e-01   # N_{1,3}
  1  4    -4.12137132e-02   # N_{1,4}
  2  1     5.41960668e-02   # N_{2,1}
  2  2     9.60334005e-01   # N_{2,2}
  2  3    -2.29093687e-01   # N_{2,3}
  2  4     1.49457246e-01   # N_{2,4}
  3  1    -4.07296753e-02   # N_{3,1}
  3  2     5.93865840e-02   # N_{3,2}
  3  3     7.01719647e-01   # N_{3,3}
  3  4     7.08804532e-01   # N_{3,4}
  4  1    -8.96914239e-02   # N_{4,1}
  4  2     2.71263108e-01   # N_{4,2}
  4  3     6.66942035e-01   # N_{4,3}
  4  4    -6.88157029e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.43971987e-01   # U_{1,1}
  1  2    -3.30025587e-01   # U_{1,2}
  2  1     3.30025587e-01   # U_{2,1}
  2  2     9.43971987e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.76322770e-01   # V_{1,1}
  1  2    -2.16318859e-01   # V_{1,2}
  2  1     2.16318859e-01   # V_{2,1}
  2  2     9.76322770e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.79633211e-01   # F_{11}
  1  2     9.60106904e-01   # F_{12}
  2  1     9.60106904e-01   # F_{21}
  2  2    -2.79633211e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.98641652e-01   # F_{11}
  1  2     5.21042394e-02   # F_{12}
  2  1    -5.21042394e-02   # F_{21}
  2  2     9.98641652e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.74419978e-01   # F_{11}
  1  2     9.84671352e-01   # F_{12}
  2  1     9.84671352e-01   # F_{21}
  2  2    -1.74419978e-01   # F_{22}
Block gauge Q= 9.06827623e+02  # SM gauge couplings
     1     3.62006507e-01   # g'(Q)MSSM DRbar
     2     6.42823049e-01   # g(Q)MSSM DRbar
     3     1.06242240e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.06827623e+02  
  3  3     8.64664783e-01   # Yt(Q)MSSM DRbar
Block yd Q= 9.06827623e+02  
  3  3     1.36732936e-01   # Yb(Q)MSSM DRbar
Block ye Q= 9.06827623e+02  
  3  3     9.96780579e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 9.06827623e+02 # Higgs mixing parameters
     1     5.17374618e+02    # mu(Q)MSSM DRbar
     2     9.66569932e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43884108e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.04982123e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 9.06827623e+02  # MSSM DRbar SUSY breaking parameters
     1     1.66999414e+02      # M_1(Q)
     2     3.10475320e+02      # M_2(Q)
     3     8.94171876e+02      # M_3(Q)
    21     7.50482199e+05      # mH1^2(Q)
    22    -2.46326304e+05      # mH2^2(Q)
    31     8.86101186e+02      # meL(Q)
    32     8.86089035e+02      # mmuL(Q)
    33     8.82478909e+02      # mtauL(Q)
    34     8.60818560e+02      # meR(Q)
    35     8.60793422e+02      # mmuR(Q)
    36     8.53306541e+02      # mtauR(Q)
    41     1.16032505e+03      # mqL1(Q)
    42     1.16032079e+03      # mqL2(Q)
    43     1.00685836e+03      # mqL3(Q)
    44     1.14176968e+03      # muR(Q)
    45     1.14176510e+03      # mcR(Q)
    46     8.03810613e+02      # mtR(Q)
    47     1.13962038e+03      # mdR(Q)
    48     1.13961622e+03      # msR(Q)
    49     1.13161901e+03      # mbR(Q)
Block au Q= 9.06827623e+02  
  1  1    -9.16500929e+02      # Au(Q)MSSM DRbar
  2  2    -9.16496817e+02      # Ac(Q)MSSM DRbar
  3  3    -7.04980046e+02      # At(Q)MSSM DRbar
Block ad Q= 9.06827623e+02  
  1  1    -1.12355712e+03      # Ad(Q)MSSM DRbar
  2  2    -1.12355331e+03      # As(Q)MSSM DRbar
  3  3    -1.04923435e+03      # Ab(Q)MSSM DRbar
Block ae Q= 9.06827623e+02  
  1  1    -2.40109702e+02      # Ae(Q)MSSM DRbar
  2  2    -2.40105372e+02      # Amu(Q)MSSM DRbar
  3  3    -2.38821298e+02      # Atau(Q)MSSM DRbar
