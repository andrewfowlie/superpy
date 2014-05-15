* gm2_1L.F
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


        subroutine gm2_1L_1
	implicit none

#include "gm2.h"

	ComplexType tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
	ComplexType tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15
	ComplexType tmp16, tmp17, tmp18, tmp19, tmp20, tmp21, tmp22
	ComplexType tmp23, tmp24, tmp25, tmp26, tmp27, tmp28, tmp29
	ComplexType tmp30, tmp31

	tmp1 = CW*ZNeu(1,2) + SW*ZNeu(1,1)

	tmp2 = CW*ZNeu(2,2) + SW*ZNeu(2,1)

	tmp3 = CW*ZNeu(3,2) + SW*ZNeu(3,1)

	tmp4 = CW*ZNeu(4,2) + SW*ZNeu(4,1)

	tmp5 = 3 + 4*Log(MCha(1)) - 4*Log(MSf(1,1,2))

	tmp6 = 3 + 4*Log(MCha(2)) - 4*Log(MSf(1,1,2))

	tmp7 = -Log(MNeu(1)) + Log(MSf(1,2,2))

	tmp8 = -Log(MNeu(2)) + Log(MSf(1,2,2))

	tmp9 = -Log(MNeu(3)) + Log(MSf(1,2,2))

	tmp10 = -Log(MNeu(4)) + Log(MSf(1,2,2))

	tmp11 = 1 - 4*Log(MNeu(1)) + 4*Log(MSf(1,2,2))

	tmp12 = 1 - 4*Log(MNeu(2)) + 4*Log(MSf(1,2,2))

	tmp13 = 1 - 4*Log(MNeu(3)) + 4*Log(MSf(1,2,2))

	tmp14 = 1 - 4*Log(MNeu(4)) + 4*Log(MSf(1,2,2))

	tmp15 = -Log(MNeu(1)) + Log(MSf(2,2,2))

	tmp16 = -Log(MNeu(2)) + Log(MSf(2,2,2))

	tmp17 = -Log(MNeu(3)) + Log(MSf(2,2,2))

	tmp18 = -Log(MNeu(4)) + Log(MSf(2,2,2))

	tmp19 = 1 - 4*Log(MNeu(1)) + 4*Log(MSf(2,2,2))

	tmp20 = 1 - 4*Log(MNeu(2)) + 4*Log(MSf(2,2,2))

	tmp21 = 1 - 4*Log(MNeu(3)) + 4*Log(MSf(2,2,2))

	tmp22 = 1 - 4*Log(MNeu(4)) + 4*Log(MSf(2,2,2))

        tmp23 = MW*(CB*tmp1*
     &      (2*MNeu2(1)**3 + 3*tmp11*MNeu2(1)**2*MSf2(1,2,2) - 
     &        6*MNeu2(1)*MSf2(1,2,2)**2 + MSf2(1,2,2)**3)) - 
     &   3*CW*(MNeu(1)*(MNeu2(1) - MSf2(1,2,2))*
     &      (MNeu2(1)**2 + 4*tmp7*MNeu2(1)*MSf2(1,2,2) - 
     &        MSf2(1,2,2)**2)*ZNeuC(1,3))

        tmp24 = MW*(CB*tmp1*
     &      (2*MNeu2(1)**3 + 3*tmp19*MNeu2(1)**2*MSf2(2,2,2) - 
     &        6*MNeu2(1)*MSf2(2,2,2)**2 + MSf2(2,2,2)**3)) - 
     &   3*CW*(MNeu(1)*(MNeu2(1) - MSf2(2,2,2))*
     &      (MNeu2(1)**2 + 4*tmp15*MNeu2(1)*MSf2(2,2,2) - 
     &        MSf2(2,2,2)**2)*ZNeuC(1,3))

        tmp25 = MW*(CB*tmp2*
     &      (2*MNeu2(2)**3 + 3*tmp12*MNeu2(2)**2*MSf2(1,2,2) - 
     &        6*MNeu2(2)*MSf2(1,2,2)**2 + MSf2(1,2,2)**3)) - 
     &   3*CW*(MNeu(2)*(MNeu2(2) - MSf2(1,2,2))*
     &      (MNeu2(2)**2 + 4*tmp8*MNeu2(2)*MSf2(1,2,2) - 
     &        MSf2(1,2,2)**2)*ZNeuC(2,3))

        tmp26 = MW*(CB*tmp2*
     &      (2*MNeu2(2)**3 + 3*tmp20*MNeu2(2)**2*MSf2(2,2,2) - 
     &        6*MNeu2(2)*MSf2(2,2,2)**2 + MSf2(2,2,2)**3)) - 
     &   3*CW*(MNeu(2)*(MNeu2(2) - MSf2(2,2,2))*
     &      (MNeu2(2)**2 + 4*tmp16*MNeu2(2)*MSf2(2,2,2) - 
     &        MSf2(2,2,2)**2)*ZNeuC(2,3))

        tmp27 = MW*(CB*tmp3*
     &      (2*MNeu2(3)**3 + 3*tmp13*MNeu2(3)**2*MSf2(1,2,2) - 
     &        6*MNeu2(3)*MSf2(1,2,2)**2 + MSf2(1,2,2)**3)) - 
     &   3*CW*(MNeu(3)*(MNeu2(3) - MSf2(1,2,2))*
     &      (MNeu2(3)**2 + 4*tmp9*MNeu2(3)*MSf2(1,2,2) - 
     &        MSf2(1,2,2)**2)*ZNeuC(3,3))

        tmp28 = MW*(CB*tmp3*
     &      (2*MNeu2(3)**3 + 3*tmp21*MNeu2(3)**2*MSf2(2,2,2) - 
     &        6*MNeu2(3)*MSf2(2,2,2)**2 + MSf2(2,2,2)**3)) - 
     &   3*CW*(MNeu(3)*(MNeu2(3) - MSf2(2,2,2))*
     &      (MNeu2(3)**2 + 4*tmp17*MNeu2(3)*MSf2(2,2,2) - 
     &        MSf2(2,2,2)**2)*ZNeuC(3,3))

        tmp29 = MW*(CB*tmp4*
     &      (2*MNeu2(4)**3 + 3*tmp14*MNeu2(4)**2*MSf2(1,2,2) - 
     &        6*MNeu2(4)*MSf2(1,2,2)**2 + MSf2(1,2,2)**3)) - 
     &   3*CW*(MNeu(4)*(MNeu2(4) - MSf2(1,2,2))*
     &      (MNeu2(4)**2 + 4*tmp10*MNeu2(4)*MSf2(1,2,2) - 
     &        MSf2(1,2,2)**2)*ZNeuC(4,3))

        tmp30 = MW*(CB*tmp4*
     &      (2*MNeu2(4)**3 + 3*tmp22*MNeu2(4)**2*MSf2(2,2,2) - 
     &        6*MNeu2(4)*MSf2(2,2,2)**2 + MSf2(2,2,2)**3)) - 
     &   3*CW*(MNeu(4)*(MNeu2(4) - MSf2(2,2,2))*
     &      (MNeu2(4)**2 + 4*tmp18*MNeu2(4)*MSf2(2,2,2) - 
     &        MSf2(2,2,2)**2)*ZNeuC(4,3))

        tmp31 = MM2/(192.D0*MW*Pi**2*SW2)*
     &   (ELGF**3*((-(4.242640687119286D0*
     &              (MCha(1)*(MCha2(1) - MSf2(1,1,2))*
     &                (MCha2(1)**2 - 4*MCha2(1)*MSf2(1,1,2) + 
     &                  tmp5*MSf2(1,1,2)**2)*UChaC(1,2)*VChaC(1,1))
     &              ) + (-(4.242640687119286D0*
     &                 (MCha(1)*(MCha2(1) - MSf2(1,1,2))*
     &                   (MCha2(1)**2 - 4*MCha2(1)*MSf2(1,1,2) + 
     &                     tmp5*MSf2(1,1,2)**2)*UCha(1,2))) + 
     &              2*MW*(CB*
     &                 (MCha2(1)**3 - 6*MCha2(1)**2*MSf2(1,1,2) + 
     &                   3*
     &                    (1 + 4*Log(MCha(1)) - 4*Log(MSf(1,1,2)))*
     &                    MCha2(1)*MSf2(1,1,2)**2 + 
     &                   2*MSf2(1,1,2)**3)*VChaC(1,1)))*VCha(1,1))/
     &         (MCha2(1) - MSf2(1,1,2))**4 + 
     &        (-(4.242640687119286D0*
     &              (MCha(2)*(MCha2(2) - MSf2(1,1,2))*
     &                (MCha2(2)**2 - 4*MCha2(2)*MSf2(1,1,2) + 
     &                  tmp6*MSf2(1,1,2)**2)*UChaC(2,2)*VChaC(2,1))
     &              ) + (-(4.242640687119286D0*
     &                 (MCha(2)*(MCha2(2) - MSf2(1,1,2))*
     &                   (MCha2(2)**2 - 4*MCha2(2)*MSf2(1,1,2) + 
     &                     tmp6*MSf2(1,1,2)**2)*UCha(2,2))) + 
     &              2*MW*(CB*
     &                 (MCha2(2)**3 - 6*MCha2(2)**2*MSf2(1,1,2) + 
     &                   3*
     &                    (1 + 4*Log(MCha(2)) - 4*Log(MSf(1,1,2)))*
     &                    MCha2(2)*MSf2(1,1,2)**2 + 
     &                   2*MSf2(1,1,2)**3)*VChaC(2,1)))*VCha(2,1))/
     &         (MCha2(2) - MSf2(1,1,2))**4))/CB

        tmp31 = tmp31 - MM/(192.D0*CW2*MW*Pi**2*SW2)*
     &    (ELGF**3*((MM*(SW*tmp23*USf(1,1,2,2)*ZNeuC(1,1)) + 
     &            CW*MM*(tmp23*USf(1,1,2,2)*ZNeuC(1,2)) - 
     &            3*tmp1*(2*MW*(CB*SW*USf(1,2,2,2)*ZNeu(1,1)) + 
     &               CW*MM*(USf(1,1,2,2)*ZNeu(1,3)))*MNeu(1)*
     &             (MNeu2(1) - MSf2(1,2,2))*
     &             (MNeu2(1)**2 + 4*tmp7*MNeu2(1)*MSf2(1,2,2) - 
     &               MSf2(1,2,2)**2))*USfC(1,1,2,2) + 
     &         2*SW*USfC(1,2,2,2)*
     &          (3*CW*MM*(MNeu(1)*(MNeu2(1) - MSf2(1,2,2))*
     &               (MNeu2(1)**2 + 4*tmp7*MNeu2(1)*MSf2(1,2,2) - 
     &                 MSf2(1,2,2)**2)*USf(1,2,2,2)*ZNeu(1,1)*
     &               ZNeu(1,3)) - 
     &            3*MW*(CB*SW*MNeu(1)*(MNeu2(1) - MSf2(1,2,2))*
     &               (MNeu2(1)**2 + 4*tmp7*MNeu2(1)*MSf2(1,2,2) - 
     &                 MSf2(1,2,2)**2)*USf(1,1,2,2)*ZNeuC(1,1)**2)+
     &              (MM*((2*MW*
     &                     (CB*SW*
     &                       (2*MNeu2(1)**3 + 
     &                       3*tmp11*MNeu2(1)**2*MSf2(1,2,2) - 
     &                       6*MNeu2(1)*MSf2(1,2,2)**2 + 
     &                       MSf2(1,2,2)**3)*ZNeu(1,1)) + 
     &                    3*CW*
     &                     (MNeu(1)*(MNeu2(1) - MSf2(1,2,2))*
     &                       (MNeu2(1)**2 + 
     &                       4*tmp7*MNeu2(1)*MSf2(1,2,2) - 
     &                       MSf2(1,2,2)**2)*ZNeuC(1,3)))*
     &                  USf(1,2,2,2)) - 
     &               3*CW*MW*
     &                (CB*MNeu(1)*(MNeu2(1) - MSf2(1,2,2))*
     &                  (MNeu2(1)**2 + 
     &                    4*tmp7*MNeu2(1)*MSf2(1,2,2) - 
     &                    MSf2(1,2,2)**2)*USf(1,1,2,2)*ZNeuC(1,2)))
     &              *ZNeuC(1,1))))/(CB*(MNeu2(1) - MSf2(1,2,2))**4)

        tmp31 = tmp31 - MM/(192.D0*CW2*MW*Pi**2*SW2)*
     &    (ELGF**3*((MM*(SW*tmp24*USf(2,1,2,2)*ZNeuC(1,1)) + 
     &            CW*MM*(tmp24*USf(2,1,2,2)*ZNeuC(1,2)) - 
     &            3*tmp1*(2*MW*(CB*SW*USf(2,2,2,2)*ZNeu(1,1)) + 
     &               CW*MM*(USf(2,1,2,2)*ZNeu(1,3)))*MNeu(1)*
     &             (MNeu2(1) - MSf2(2,2,2))*
     &             (MNeu2(1)**2 + 4*tmp15*MNeu2(1)*MSf2(2,2,2) - 
     &               MSf2(2,2,2)**2))*USfC(2,1,2,2) + 
     &         2*SW*USfC(2,2,2,2)*
     &          (3*CW*MM*(MNeu(1)*(MNeu2(1) - MSf2(2,2,2))*
     &               (MNeu2(1)**2 + 4*tmp15*MNeu2(1)*MSf2(2,2,2) - 
     &                 MSf2(2,2,2)**2)*USf(2,2,2,2)*ZNeu(1,1)*
     &               ZNeu(1,3)) - 
     &            3*MW*(CB*SW*MNeu(1)*(MNeu2(1) - MSf2(2,2,2))*
     &               (MNeu2(1)**2 + 4*tmp15*MNeu2(1)*MSf2(2,2,2) - 
     &                 MSf2(2,2,2)**2)*USf(2,1,2,2)*ZNeuC(1,1)**2)+
     &              (MM*((2*MW*
     &                     (CB*SW*
     &                       (2*MNeu2(1)**3 + 
     &                       3*tmp19*MNeu2(1)**2*MSf2(2,2,2) - 
     &                       6*MNeu2(1)*MSf2(2,2,2)**2 + 
     &                       MSf2(2,2,2)**3)*ZNeu(1,1)) + 
     &                    3*CW*
     &                     (MNeu(1)*(MNeu2(1) - MSf2(2,2,2))*
     &                       (MNeu2(1)**2 + 
     &                       4*tmp15*MNeu2(1)*MSf2(2,2,2) - 
     &                       MSf2(2,2,2)**2)*ZNeuC(1,3)))*
     &                  USf(2,2,2,2)) - 
     &               3*CW*MW*
     &                (CB*MNeu(1)*(MNeu2(1) - MSf2(2,2,2))*
     &                  (MNeu2(1)**2 + 
     &                    4*tmp15*MNeu2(1)*MSf2(2,2,2) - 
     &                    MSf2(2,2,2)**2)*USf(2,1,2,2)*ZNeuC(1,2)))
     &              *ZNeuC(1,1))))/(CB*(MNeu2(1) - MSf2(2,2,2))**4)

        tmp31 = tmp31 - MM/(192.D0*CW2*MW*Pi**2*SW2)*
     &    (ELGF**3*((MM*(SW*tmp25*USf(1,1,2,2)*ZNeuC(2,1)) + 
     &            CW*MM*(tmp25*USf(1,1,2,2)*ZNeuC(2,2)) - 
     &            3*tmp2*(2*MW*(CB*SW*USf(1,2,2,2)*ZNeu(2,1)) + 
     &               CW*MM*(USf(1,1,2,2)*ZNeu(2,3)))*MNeu(2)*
     &             (MNeu2(2) - MSf2(1,2,2))*
     &             (MNeu2(2)**2 + 4*tmp8*MNeu2(2)*MSf2(1,2,2) - 
     &               MSf2(1,2,2)**2))*USfC(1,1,2,2) + 
     &         2*SW*USfC(1,2,2,2)*
     &          (3*CW*MM*(MNeu(2)*(MNeu2(2) - MSf2(1,2,2))*
     &               (MNeu2(2)**2 + 4*tmp8*MNeu2(2)*MSf2(1,2,2) - 
     &                 MSf2(1,2,2)**2)*USf(1,2,2,2)*ZNeu(2,1)*
     &               ZNeu(2,3)) - 
     &            3*MW*(CB*SW*MNeu(2)*(MNeu2(2) - MSf2(1,2,2))*
     &               (MNeu2(2)**2 + 4*tmp8*MNeu2(2)*MSf2(1,2,2) - 
     &                 MSf2(1,2,2)**2)*USf(1,1,2,2)*ZNeuC(2,1)**2)+
     &              (MM*((2*MW*
     &                     (CB*SW*
     &                       (2*MNeu2(2)**3 + 
     &                       3*tmp12*MNeu2(2)**2*MSf2(1,2,2) - 
     &                       6*MNeu2(2)*MSf2(1,2,2)**2 + 
     &                       MSf2(1,2,2)**3)*ZNeu(2,1)) + 
     &                    3*CW*
     &                     (MNeu(2)*(MNeu2(2) - MSf2(1,2,2))*
     &                       (MNeu2(2)**2 + 
     &                       4*tmp8*MNeu2(2)*MSf2(1,2,2) - 
     &                       MSf2(1,2,2)**2)*ZNeuC(2,3)))*
     &                  USf(1,2,2,2)) - 
     &               3*CW*MW*
     &                (CB*MNeu(2)*(MNeu2(2) - MSf2(1,2,2))*
     &                  (MNeu2(2)**2 + 
     &                    4*tmp8*MNeu2(2)*MSf2(1,2,2) - 
     &                    MSf2(1,2,2)**2)*USf(1,1,2,2)*ZNeuC(2,2)))
     &              *ZNeuC(2,1))))/(CB*(MNeu2(2) - MSf2(1,2,2))**4)

        tmp31 = tmp31 - MM/(192.D0*CW2*MW*Pi**2*SW2)*
     &    (ELGF**3*((MM*(SW*tmp26*USf(2,1,2,2)*ZNeuC(2,1)) + 
     &            CW*MM*(tmp26*USf(2,1,2,2)*ZNeuC(2,2)) - 
     &            3*tmp2*(2*MW*(CB*SW*USf(2,2,2,2)*ZNeu(2,1)) + 
     &               CW*MM*(USf(2,1,2,2)*ZNeu(2,3)))*MNeu(2)*
     &             (MNeu2(2) - MSf2(2,2,2))*
     &             (MNeu2(2)**2 + 4*tmp16*MNeu2(2)*MSf2(2,2,2) - 
     &               MSf2(2,2,2)**2))*USfC(2,1,2,2) + 
     &         2*SW*USfC(2,2,2,2)*
     &          (3*CW*MM*(MNeu(2)*(MNeu2(2) - MSf2(2,2,2))*
     &               (MNeu2(2)**2 + 4*tmp16*MNeu2(2)*MSf2(2,2,2) - 
     &                 MSf2(2,2,2)**2)*USf(2,2,2,2)*ZNeu(2,1)*
     &               ZNeu(2,3)) - 
     &            3*MW*(CB*SW*MNeu(2)*(MNeu2(2) - MSf2(2,2,2))*
     &               (MNeu2(2)**2 + 4*tmp16*MNeu2(2)*MSf2(2,2,2) - 
     &                 MSf2(2,2,2)**2)*USf(2,1,2,2)*ZNeuC(2,1)**2)+
     &              (MM*((2*MW*
     &                     (CB*SW*
     &                       (2*MNeu2(2)**3 + 
     &                       3*tmp20*MNeu2(2)**2*MSf2(2,2,2) - 
     &                       6*MNeu2(2)*MSf2(2,2,2)**2 + 
     &                       MSf2(2,2,2)**3)*ZNeu(2,1)) + 
     &                    3*CW*
     &                     (MNeu(2)*(MNeu2(2) - MSf2(2,2,2))*
     &                       (MNeu2(2)**2 + 
     &                       4*tmp16*MNeu2(2)*MSf2(2,2,2) - 
     &                       MSf2(2,2,2)**2)*ZNeuC(2,3)))*
     &                  USf(2,2,2,2)) - 
     &               3*CW*MW*
     &                (CB*MNeu(2)*(MNeu2(2) - MSf2(2,2,2))*
     &                  (MNeu2(2)**2 + 
     &                    4*tmp16*MNeu2(2)*MSf2(2,2,2) - 
     &                    MSf2(2,2,2)**2)*USf(2,1,2,2)*ZNeuC(2,2)))
     &              *ZNeuC(2,1))))/(CB*(MNeu2(2) - MSf2(2,2,2))**4)

        tmp31 = tmp31 - MM/(192.D0*CW2*MW*Pi**2*SW2)*
     &    (ELGF**3*((MM*(SW*tmp27*USf(1,1,2,2)*ZNeuC(3,1)) + 
     &            CW*MM*(tmp27*USf(1,1,2,2)*ZNeuC(3,2)) - 
     &            3*tmp3*(2*MW*(CB*SW*USf(1,2,2,2)*ZNeu(3,1)) + 
     &               CW*MM*(USf(1,1,2,2)*ZNeu(3,3)))*MNeu(3)*
     &             (MNeu2(3) - MSf2(1,2,2))*
     &             (MNeu2(3)**2 + 4*tmp9*MNeu2(3)*MSf2(1,2,2) - 
     &               MSf2(1,2,2)**2))*USfC(1,1,2,2) + 
     &         2*SW*USfC(1,2,2,2)*
     &          (3*CW*MM*(MNeu(3)*(MNeu2(3) - MSf2(1,2,2))*
     &               (MNeu2(3)**2 + 4*tmp9*MNeu2(3)*MSf2(1,2,2) - 
     &                 MSf2(1,2,2)**2)*USf(1,2,2,2)*ZNeu(3,1)*
     &               ZNeu(3,3)) - 
     &            3*MW*(CB*SW*MNeu(3)*(MNeu2(3) - MSf2(1,2,2))*
     &               (MNeu2(3)**2 + 4*tmp9*MNeu2(3)*MSf2(1,2,2) - 
     &                 MSf2(1,2,2)**2)*USf(1,1,2,2)*ZNeuC(3,1)**2)+
     &              (MM*((2*MW*
     &                     (CB*SW*
     &                       (2*MNeu2(3)**3 + 
     &                       3*tmp13*MNeu2(3)**2*MSf2(1,2,2) - 
     &                       6*MNeu2(3)*MSf2(1,2,2)**2 + 
     &                       MSf2(1,2,2)**3)*ZNeu(3,1)) + 
     &                    3*CW*
     &                     (MNeu(3)*(MNeu2(3) - MSf2(1,2,2))*
     &                       (MNeu2(3)**2 + 
     &                       4*tmp9*MNeu2(3)*MSf2(1,2,2) - 
     &                       MSf2(1,2,2)**2)*ZNeuC(3,3)))*
     &                  USf(1,2,2,2)) - 
     &               3*CW*MW*
     &                (CB*MNeu(3)*(MNeu2(3) - MSf2(1,2,2))*
     &                  (MNeu2(3)**2 + 
     &                    4*tmp9*MNeu2(3)*MSf2(1,2,2) - 
     &                    MSf2(1,2,2)**2)*USf(1,1,2,2)*ZNeuC(3,2)))
     &              *ZNeuC(3,1))))/(CB*(MNeu2(3) - MSf2(1,2,2))**4)

        tmp31 = tmp31 - MM/(192.D0*CW2*MW*Pi**2*SW2)*
     &    (ELGF**3*((MM*(SW*tmp28*USf(2,1,2,2)*ZNeuC(3,1)) + 
     &            CW*MM*(tmp28*USf(2,1,2,2)*ZNeuC(3,2)) - 
     &            3*tmp3*(2*MW*(CB*SW*USf(2,2,2,2)*ZNeu(3,1)) + 
     &               CW*MM*(USf(2,1,2,2)*ZNeu(3,3)))*MNeu(3)*
     &             (MNeu2(3) - MSf2(2,2,2))*
     &             (MNeu2(3)**2 + 4*tmp17*MNeu2(3)*MSf2(2,2,2) - 
     &               MSf2(2,2,2)**2))*USfC(2,1,2,2) + 
     &         2*SW*USfC(2,2,2,2)*
     &          (3*CW*MM*(MNeu(3)*(MNeu2(3) - MSf2(2,2,2))*
     &               (MNeu2(3)**2 + 4*tmp17*MNeu2(3)*MSf2(2,2,2) - 
     &                 MSf2(2,2,2)**2)*USf(2,2,2,2)*ZNeu(3,1)*
     &               ZNeu(3,3)) - 
     &            3*MW*(CB*SW*MNeu(3)*(MNeu2(3) - MSf2(2,2,2))*
     &               (MNeu2(3)**2 + 4*tmp17*MNeu2(3)*MSf2(2,2,2) - 
     &                 MSf2(2,2,2)**2)*USf(2,1,2,2)*ZNeuC(3,1)**2)+
     &              (MM*((2*MW*
     &                     (CB*SW*
     &                       (2*MNeu2(3)**3 + 
     &                       3*tmp21*MNeu2(3)**2*MSf2(2,2,2) - 
     &                       6*MNeu2(3)*MSf2(2,2,2)**2 + 
     &                       MSf2(2,2,2)**3)*ZNeu(3,1)) + 
     &                    3*CW*
     &                     (MNeu(3)*(MNeu2(3) - MSf2(2,2,2))*
     &                       (MNeu2(3)**2 + 
     &                       4*tmp17*MNeu2(3)*MSf2(2,2,2) - 
     &                       MSf2(2,2,2)**2)*ZNeuC(3,3)))*
     &                  USf(2,2,2,2)) - 
     &               3*CW*MW*
     &                (CB*MNeu(3)*(MNeu2(3) - MSf2(2,2,2))*
     &                  (MNeu2(3)**2 + 
     &                    4*tmp17*MNeu2(3)*MSf2(2,2,2) - 
     &                    MSf2(2,2,2)**2)*USf(2,1,2,2)*ZNeuC(3,2)))
     &              *ZNeuC(3,1))))/(CB*(MNeu2(3) - MSf2(2,2,2))**4)

        tmp31 = tmp31 - MM/(192.D0*CW2*MW*Pi**2*SW2)*
     &    (ELGF**3*((MM*(SW*tmp29*USf(1,1,2,2)*ZNeuC(4,1)) + 
     &            CW*MM*(tmp29*USf(1,1,2,2)*ZNeuC(4,2)) - 
     &            3*tmp4*(2*MW*(CB*SW*USf(1,2,2,2)*ZNeu(4,1)) + 
     &               CW*MM*(USf(1,1,2,2)*ZNeu(4,3)))*MNeu(4)*
     &             (MNeu2(4) - MSf2(1,2,2))*
     &             (MNeu2(4)**2 + 4*tmp10*MNeu2(4)*MSf2(1,2,2) - 
     &               MSf2(1,2,2)**2))*USfC(1,1,2,2) + 
     &         2*SW*USfC(1,2,2,2)*
     &          (3*CW*MM*(MNeu(4)*(MNeu2(4) - MSf2(1,2,2))*
     &               (MNeu2(4)**2 + 4*tmp10*MNeu2(4)*MSf2(1,2,2) - 
     &                 MSf2(1,2,2)**2)*USf(1,2,2,2)*ZNeu(4,1)*
     &               ZNeu(4,3)) - 
     &            3*MW*(CB*SW*MNeu(4)*(MNeu2(4) - MSf2(1,2,2))*
     &               (MNeu2(4)**2 + 4*tmp10*MNeu2(4)*MSf2(1,2,2) - 
     &                 MSf2(1,2,2)**2)*USf(1,1,2,2)*ZNeuC(4,1)**2)+
     &              (MM*((2*MW*
     &                     (CB*SW*
     &                       (2*MNeu2(4)**3 + 
     &                       3*tmp14*MNeu2(4)**2*MSf2(1,2,2) - 
     &                       6*MNeu2(4)*MSf2(1,2,2)**2 + 
     &                       MSf2(1,2,2)**3)*ZNeu(4,1)) + 
     &                    3*CW*
     &                     (MNeu(4)*(MNeu2(4) - MSf2(1,2,2))*
     &                       (MNeu2(4)**2 + 
     &                       4*tmp10*MNeu2(4)*MSf2(1,2,2) - 
     &                       MSf2(1,2,2)**2)*ZNeuC(4,3)))*
     &                  USf(1,2,2,2)) - 
     &               3*CW*MW*
     &                (CB*MNeu(4)*(MNeu2(4) - MSf2(1,2,2))*
     &                  (MNeu2(4)**2 + 
     &                    4*tmp10*MNeu2(4)*MSf2(1,2,2) - 
     &                    MSf2(1,2,2)**2)*USf(1,1,2,2)*ZNeuC(4,2)))
     &              *ZNeuC(4,1))))/(CB*(MNeu2(4) - MSf2(1,2,2))**4)

        tmp31 = tmp31 - MM/(192.D0*CW2*MW*Pi**2*SW2)*
     &    (ELGF**3*((MM*(SW*tmp30*USf(2,1,2,2)*ZNeuC(4,1)) + 
     &            CW*MM*(tmp30*USf(2,1,2,2)*ZNeuC(4,2)) - 
     &            3*tmp4*(2*MW*(CB*SW*USf(2,2,2,2)*ZNeu(4,1)) + 
     &               CW*MM*(USf(2,1,2,2)*ZNeu(4,3)))*MNeu(4)*
     &             (MNeu2(4) - MSf2(2,2,2))*
     &             (MNeu2(4)**2 + 4*tmp18*MNeu2(4)*MSf2(2,2,2) - 
     &               MSf2(2,2,2)**2))*USfC(2,1,2,2) + 
     &         2*SW*USfC(2,2,2,2)*
     &          (3*CW*MM*(MNeu(4)*(MNeu2(4) - MSf2(2,2,2))*
     &               (MNeu2(4)**2 + 4*tmp18*MNeu2(4)*MSf2(2,2,2) - 
     &                 MSf2(2,2,2)**2)*USf(2,2,2,2)*ZNeu(4,1)*
     &               ZNeu(4,3)) - 
     &            3*MW*(CB*SW*MNeu(4)*(MNeu2(4) - MSf2(2,2,2))*
     &               (MNeu2(4)**2 + 4*tmp18*MNeu2(4)*MSf2(2,2,2) - 
     &                 MSf2(2,2,2)**2)*USf(2,1,2,2)*ZNeuC(4,1)**2)+
     &              (MM*((2*MW*
     &                     (CB*SW*
     &                       (2*MNeu2(4)**3 + 
     &                       3*tmp22*MNeu2(4)**2*MSf2(2,2,2) - 
     &                       6*MNeu2(4)*MSf2(2,2,2)**2 + 
     &                       MSf2(2,2,2)**3)*ZNeu(4,1)) + 
     &                    3*CW*
     &                     (MNeu(4)*(MNeu2(4) - MSf2(2,2,2))*
     &                       (MNeu2(4)**2 + 
     &                       4*tmp18*MNeu2(4)*MSf2(2,2,2) - 
     &                       MSf2(2,2,2)**2)*ZNeuC(4,3)))*
     &                  USf(2,2,2,2)) - 
     &               3*CW*MW*
     &                (CB*MNeu(4)*(MNeu2(4) - MSf2(2,2,2))*
     &                  (MNeu2(4)**2 + 
     &                    4*tmp18*MNeu2(4)*MSf2(2,2,2) - 
     &                    MSf2(2,2,2)**2)*USf(2,1,2,2)*ZNeuC(4,2)))
     &              *ZNeuC(4,1))))/(CB*(MNeu2(4) - MSf2(2,2,2))**4)

	gm2_1L = tmp31/ELGF
	end
