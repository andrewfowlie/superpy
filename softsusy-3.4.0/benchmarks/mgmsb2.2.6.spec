# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-3.50738081e+03
# dtanbeta/dmu=3.41942757e-02
# mu=8.81264120e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=2.79318258e-01
# dtanbeta/dm3sq=-1.60403278e-04
# m3sq= 6.74555117e+04
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
     1    1.70000000e+05   # lambda
     2    1.00000000e+14   # M_mess
     5    1.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.36625160e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03954587e+01   # MW
        25     1.16189296e+02   # h0
        35     1.11959886e+03   # H0
        36     1.11955634e+03   # A0
        37     1.12263586e+03   # H+
   1000021     1.26309122e+03   # ~g
   1000022     2.27572916e+02   # ~neutralino(1)
   1000023     4.38363474e+02   # ~neutralino(2)
   1000024     4.38499950e+02   # ~chargino(1)
   1000025    -8.78017656e+02   # ~neutralino(3)
   1000035     8.84575579e+02   # ~neutralino(4)
   1000037     8.85312020e+02   # ~chargino(2)
   1000039     4.00022143e+00   # ~gravitino
   1000001     1.55536684e+03   # ~d_L
   1000002     1.55349387e+03   # ~u_L
   1000003     1.55536069e+03   # ~s_L
   1000004     1.55348771e+03   # ~c_L
   1000005     1.37521033e+03   # ~b_1
   1000006     1.09372521e+03   # ~t_1
   1000011     7.52700332e+02   # ~e_L
   1000012     7.48204743e+02   # ~nue_L
   1000013     7.52684482e+02   # ~mu_L
   1000014     7.48188807e+02   # ~numu_L
   1000015     4.95092389e+02   # ~stau_1
   1000016     7.42989662e+02   # ~nu_tau_L
   2000001     1.39757457e+03   # ~d_R
   2000002     1.42654834e+03   # ~u_R
   2000003     1.39756544e+03   # ~s_R
   2000004     1.42654383e+03   # ~c_R
   2000005     1.40782669e+03   # ~b_2
   2000006     1.41879775e+03   # ~t_2
   2000011     5.11970798e+02   # ~e_R
   2000013     5.11923864e+02   # ~mu_R
   2000015     7.48687892e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.00210572e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.98214769e-01   # N_{1,1}
  1  2    -8.06665113e-03   # N_{1,2}
  1  3     5.62695205e-02   # N_{1,3}
  1  4    -1.83288258e-02   # N_{1,4}
  2  1     1.60195379e-02   # N_{2,1}
  2  2     9.90365101e-01   # N_{2,2}
  2  3    -1.20750749e-01   # N_{2,3}
  2  4     6.58756232e-02   # N_{2,4}
  3  1    -2.64274542e-02   # N_{3,1}
  3  2     3.93093951e-02   # N_{3,2}
  3  3     7.04963379e-01   # N_{3,3}
  3  4     7.07660226e-01   # N_{3,4}
  4  1    -5.11100745e-02   # N_{4,1}
  4  2     1.32539305e-01   # N_{4,2}
  4  3     6.96620149e-01   # N_{4,3}
  4  4    -7.03236419e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.85132517e-01   # U_{1,1}
  1  2    -1.71796173e-01   # U_{1,2}
  2  1     1.71796173e-01   # U_{2,1}
  2  2     9.85132517e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.95564674e-01   # V_{1,1}
  1  2    -9.40796495e-02   # V_{1,2}
  2  1     9.40796495e-02   # V_{2,1}
  2  2     9.95564674e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     1.76678670e-01   # F_{11}
  1  2     9.84268585e-01   # F_{12}
  2  1     9.84268585e-01   # F_{21}
  2  2    -1.76678670e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     3.87294407e-01   # F_{11}
  1  2     9.21956096e-01   # F_{12}
  2  1     9.21956096e-01   # F_{21}
  2  2    -3.87294407e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     7.40722496e-02   # F_{11}
  1  2     9.97252878e-01   # F_{12}
  2  1     9.97252878e-01   # F_{21}
  2  2    -7.40722496e-02   # F_{22}
Block gauge Q= 1.21064741e+03  # SM gauge couplings
     1     3.62925551e-01   # g'(Q)MSSM DRbar
     2     6.41435482e-01   # g(Q)MSSM DRbar
     3     1.04842906e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.21064741e+03  
  3  3     8.52151461e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.21064741e+03  
  3  3     1.96447814e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.21064741e+03  
  3  3     1.50626438e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.21064741e+03 # Higgs mixing parameters
     1     8.71564849e+02    # mu(Q)MSSM DRbar
     2     1.44703605e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43741099e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.30815583e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.21064741e+03  # MSSM DRbar SUSY breaking parameters
     1     2.31054096e+02      # M_1(Q)
     2     4.26885296e+02      # M_2(Q)
     3     1.16096881e+03      # M_3(Q)
    21     4.80639212e+05      # mH1^2(Q)
    22    -7.33080764e+05      # mH2^2(Q)
    31     7.48885986e+02      # meL(Q)
    32     7.48870061e+02      # mmuL(Q)
    33     7.44006483e+02      # mtauL(Q)
    34     5.06910708e+02      # meR(Q)
    35     5.06863336e+02      # mmuR(Q)
    36     4.92227250e+02      # mtauR(Q)
    41     1.52040306e+03      # mqL1(Q)
    42     1.52039675e+03      # mqL2(Q)
    43     1.37304601e+03      # mqL3(Q)
    44     1.39180545e+03      # muR(Q)
    45     1.39180077e+03      # mcR(Q)
    46     1.05712292e+03      # mtR(Q)
    47     1.36060895e+03      # mdR(Q)
    48     1.36059949e+03      # msR(Q)
    49     1.34335370e+03      # mbR(Q)
Block au Q= 1.21064741e+03  
  1  1    -1.07296522e+03      # Au(Q)MSSM DRbar
  2  2    -1.07295985e+03      # Ac(Q)MSSM DRbar
  3  3    -8.55312353e+02      # At(Q)MSSM DRbar
Block ad Q= 1.21064741e+03  
  1  1    -1.27995978e+03      # Ad(Q)MSSM DRbar
  2  2    -1.27995237e+03      # As(Q)MSSM DRbar
  3  3    -1.19804493e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.21064741e+03  
  1  1    -2.40132767e+02      # Ae(Q)MSSM DRbar
  2  2    -2.40124590e+02      # Amu(Q)MSSM DRbar
  3  3    -2.37634291e+02      # Atau(Q)MSSM DRbar
