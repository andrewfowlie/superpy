#ifndef _LOGINT_
#define _LOGINT_

#define I1(a,b) Jslog5( 4*(a),-4*(a),0.,0.,b )
#define I2(a,b) Jslog5( -(a),(a)+(b)-1.,a,-(a)+(b)-1.,1. )
#define I3(a,b) Jslog5( -(a),(a)-(b)+1., a,-(a)-(b)+1.,b )

//#define I_41(r1,r2,r3,r4) Jslog5(4*(r1)      ,-4*(r1)+(r2)-(r3)    , 4*(r4),-4*(r4)+(r2)-(r3),r3)
//#define I_42(r1,r2,r3,r4) Jslog5(-(r1)+2*(r4),(r1)-2*(r4)+(r2)-(r3), (r1)  ,-(r1)+(r2)-(r3)  ,r3)
#define I_43(r1,r2,r3,r4) Jslog5(-(r1)+2*(r4),(r1)-2*(r4)-(r2)+(r3), (r1)  ,-(r1)-(r2)+(r3)  ,r2)
#define I_5(r1,r2,r3,r4,r5) (Jflog5((r1)-2*(r5),r2,r3,r4,r5)-Jflog5(-(r1),r2,r4,r3,r5))

extern double  Jslog5(double An,double Bn,double Ad,double Bd,double C);
extern double  Jflog5(double r1,double r2,double r3,double r4,double r5);

#endif
