
        dl22 = (MTy2*((-2 + 
     &          2*Log((MStop2(1)*MStop2(2))/MTy2**2))*MStop2(1)*
     &        MStop2(2) + MTy2*(-Xt2 + MStop2(4))))/
     &   (MStop2(1)*MStop2(2))

	if( imass .ne. 3 ) then

        dcH = 1/2.D0*(MTy2*
     &      (2*Log(MStop2(1))*
     &         ((2*MStop2(1)**2 - MSQ2*(MTy2 + 2*MStop2(1)))*
     &            MStop2(3)**2 + 
     &           MTy2*Xt2*(-2*MStop2(1)**2 + MSQ2*MStop2(4))) + 
     &        (MSQ2*Log(MSQ2)*MStop2(3)**3*
     &            (-2*MTy2*Xt2 + MStop2(3)**2 + 
     &              (2*MSQ2 + 2*MTy2 - MStop2(4))*
     &               (-2*MSQ2 + MStop2(4))) + 
     &           2*(-MSQ2 + MStop2(1))*
     &            (MTy2*(-MSQ2 + MStop2(2))*MStop2(3)*
     &               (MStop2(3)**2 + Xt2*(2*MSQ2 - MStop2(4))) + 
     &              Log(MStop2(2))*(-MSQ2 + MStop2(1))*
     &               (-2*MStop2(2)**2*
     &                  (-(MTy2*Xt2) + MStop2(3)**2) + 
     &                 MSQ2*
     &                  (2*MStop2(2)*MStop2(3)**2 + 
     &                    MTy2*(MStop2(3)**2 - Xt2*MStop2(4))))))/
     &         (-MSQ2 + MStop2(2))**2))/
     &    ((-MSQ2 + MStop2(1))**2*MStop2(3)**3)

	endif

	dcA = 0

	dlA = 0

	if( abs(MStop2(3)) .gt. 1D-14 ) then

        dcA = (MTy2*(2*Log(MStop2(2)/MStop2(1))*
     &        (MStop2(3)**2*MStop2(4) + 
     &          MTy2*(MStop2(3)**2 - 3*Xt2*MStop2(4))) + 
     &       (MStop2(3)**5 + 
     &          MStop2(3)*(3*MTy2*Xt2*MStop2(4)**2 + 
     &             MStop2(3)**2*
     &              (MTy2*(-2*Xt2 - MStop2(4)) - MStop2(4)**2)))/
     &        (MStop2(1)*MStop2(2))))/MStop2(3)**5

        dlA = (MTy2*(2*Log(MStop2(2)/MStop2(1))*
     &        (-(MTy2*Xt2) + MStop2(3)**2) + 
     &       (MTy2*MStop2(3)*(-MStop2(3)**2 + Xt2*MStop2(4)))/
     &        (MStop2(1)*MStop2(2))))/MStop2(3)**3

	endif

