# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-1.43086599e+03
# dtanbeta/dmu=9.06850960e-02
# mu=3.68686876e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=2.82197938e-01
# dtanbeta/dm3sq=-8.79659982e-04
# m3sq= 1.82960665e+04
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
     1    9.00000000e+04   # lambda
     2    1.00000000e+05   # M_mess
     5    1.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=4.58436368e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04035161e+01   # MW
        25     1.12330751e+02   # h0
        35     4.62225358e+02   # H0
        36     4.61891526e+02   # A0
        37     4.69151635e+02   # H+
   1000021     8.64571955e+02   # ~g
   1000022     1.45061841e+02   # ~neutralino(1)
   1000023     2.65582324e+02   # ~neutralino(2)
   1000024     2.64827272e+02   # ~chargino(1)
   1000025    -3.68263011e+02   # ~neutralino(3)
   1000035     4.01172207e+02   # ~neutralino(4)
   1000037     4.01215882e+02   # ~chargino(2)
   1000039     2.11776429e-09   # ~gravitino
   1000001     1.03056734e+03   # ~d_L
   1000002     1.02767032e+03   # ~u_L
   1000003     1.03056593e+03   # ~s_L
   1000004     1.02766891e+03   # ~c_L
   1000005     9.76792970e+02   # ~b_1
   1000006     9.13291009e+02   # ~t_1
   1000011     3.17799507e+02   # ~e_L
   1000012     3.07429997e+02   # ~nue_L
   1000013     3.17798149e+02   # ~mu_L
   1000014     3.07428595e+02   # ~numu_L
   1000015     1.52910634e+02   # ~stau_1
   1000016     3.06808492e+02   # ~nu_tau_L
   2000001     9.85001515e+02   # ~d_R
   2000002     9.87536569e+02   # ~u_R
   2000003     9.84999559e+02   # ~s_R
   2000004     9.87535570e+02   # ~c_R
   2000005     9.93416242e+02   # ~b_2
   2000006     1.00616337e+03   # ~t_2
   2000011     1.58936640e+02   # ~e_R
   2000013     1.58931151e+02   # ~mu_R
   2000015     3.19079969e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.59173549e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.83917914e-01   # N_{1,1}
  1  2    -4.32799989e-02   # N_{1,2}
  1  3     1.57612145e-01   # N_{1,3}
  1  4    -7.20471478e-02   # N_{1,4}
  2  1     1.21887636e-01   # N_{2,1}
  2  2     8.64324149e-01   # N_{2,2}
  2  3    -3.88777384e-01   # N_{2,3}
  2  4     2.94854737e-01   # N_{2,4}
  3  1    -5.62204366e-02   # N_{3,1}
  3  2     7.93596661e-02   # N_{3,2}
  3  3     6.97354144e-01   # N_{3,3}
  3  4     7.10097531e-01   # N_{3,4}
  4  1    -1.17848228e-01   # N_{4,1}
  4  2     4.94745037e-01   # N_{4,2}
  4  3     5.81126281e-01   # N_{4,3}
  4  4    -6.35319911e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     8.23299266e-01   # U_{1,1}
  1  2    -5.67607540e-01   # U_{1,2}
  2  1     5.67607540e-01   # U_{2,1}
  2  2     8.23299266e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.02288992e-01   # V_{1,1}
  1  2    -4.31131737e-01   # V_{1,2}
  2  1     4.31131737e-01   # V_{2,1}
  2  2     9.02288992e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.77255929e-01   # F_{11}
  1  2     9.60796102e-01   # F_{12}
  2  1     9.60796102e-01   # F_{21}
  2  2    -2.77255929e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     4.75880605e-01   # F_{11}
  1  2     8.79509892e-01   # F_{12}
  2  1     8.79509892e-01   # F_{21}
  2  2    -4.75880605e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.22955251e-01   # F_{11}
  1  2     9.92412216e-01   # F_{12}
  2  1     9.92412216e-01   # F_{21}
  2  2    -1.22955251e-01   # F_{22}
Block gauge Q= 9.37893365e+02  # SM gauge couplings
     1     3.62940586e-01   # g'(Q)MSSM DRbar
     2     6.44931689e-01   # g(Q)MSSM DRbar
     3     1.06542076e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.37893365e+02  
  3  3     8.66418658e-01   # Yt(Q)MSSM DRbar
Block yd Q= 9.37893365e+02  
  3  3     2.04715624e-01   # Yb(Q)MSSM DRbar
Block ye Q= 9.37893365e+02  
  3  3     1.51454054e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 9.37893365e+02 # Higgs mixing parameters
     1     3.59449508e+02    # mu(Q)MSSM DRbar
     2     1.45087756e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43642973e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     2.36651397e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 9.37893365e+02  # MSSM DRbar SUSY breaking parameters
     1     1.51689601e+02      # M_1(Q)
     2     2.86308742e+02      # M_2(Q)
     3     7.85300416e+02      # M_3(Q)
    21     8.47389560e+04      # mH1^2(Q)
    22    -1.13696833e+05      # mH2^2(Q)
    31     3.11339118e+02      # meL(Q)
    32     3.11337729e+02      # mmuL(Q)
    33     3.10913691e+02      # mtauL(Q)
    34     1.48353914e+02      # meR(Q)
    35     1.48348027e+02      # mmuR(Q)
    36     1.46540306e+02      # mtauR(Q)
    41     1.00167355e+03      # mqL1(Q)
    42     1.00167210e+03      # mqL2(Q)
    43     9.66802617e+02      # mqL3(Q)
    44     9.60964111e+02      # muR(Q)
    45     9.60963075e+02      # mcR(Q)
    46     8.89691747e+02      # mtR(Q)
    47     9.57082737e+02      # mdR(Q)
    48     9.57080707e+02      # msR(Q)
    49     9.53122825e+02      # mbR(Q)
Block au Q= 9.37893365e+02  
  1  1    -2.41594367e+02      # Au(Q)MSSM DRbar
  2  2    -2.41594022e+02      # Ac(Q)MSSM DRbar
  3  3    -2.27470337e+02      # At(Q)MSSM DRbar
Block ad Q= 9.37893365e+02  
  1  1    -2.57312744e+02      # Ad(Q)MSSM DRbar
  2  2    -2.57312263e+02      # As(Q)MSSM DRbar
  3  3    -2.51992432e+02      # Ab(Q)MSSM DRbar
Block ae Q= 9.37893365e+02  
  1  1    -2.43961793e+01      # Ae(Q)MSSM DRbar
  2  2    -2.43959987e+01      # Amu(Q)MSSM DRbar
  3  3    -2.43409012e+01      # Atau(Q)MSSM DRbar
