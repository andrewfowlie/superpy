#include"sf_epp.h"
#include "crt_util.h"
#include<math.h>
static double xin[2] ={0.938, 0.938}; /* Mass of proton [GeV] */
static int charge[2] ={1, 1};         /* Charge of proton [e] */
static double q2_max[2] ={2., 2.};    /* Max of photons Abs(q^2) [GeV^2] */
static double pt_min[2]={0.1,0.1};    /* Minimal Pt cut [GeV] of outgoing protons*/
static const double Maxpt_min = 1.;   /* Maximal of minimal  Pt cut [GeV], */
/* controls approximation of q2 as a function of pt  */

/* ********************************************************* */
/*  Equivalent photon approximation structure function.    * */
/*     Improved Weizsaecker-Williams formula               * */
/*   V.M.Budnev et al., Phys.Rep. 15C (1975) 181           * */
/* ********************************************************* */


int  p_epp__ (int*pNum){ if(pNum[0]==22 && pNum[1]==0) return 1; else return 0; }

void 
n_epp__ (int i, char *name)
{
  sprintf (name, "Proton.Photon(m=%.10G Ch=%d Q=%.8G Pt>%.9G)", xin[i - 1], charge[i - 1], q2_max[i - 1], pt_min[i-1]);
}


int 
r_epp__ (int i, char *name)
{ 
  double xin__, q2_max__, pt_min__;
  int charge__;
  if(4 != sscanf (name, "Proton.Photon(m=%lf%*[^=]%*c%d%*[^=]%*c%lf%*[^>]%*c%lf)", &xin__, &charge__, &q2_max__, &pt_min__))
   goto L10;

  if (xin__ <= 0.)
    goto L10;
  if (q2_max__ <= 0.)
    goto L10;
  if (pt_min__ < 0. || pt_min__ > Maxpt_min) goto L10;
  i--;
  xin[i] = xin__;
  q2_max[i] = q2_max__; /* [GeV^2]  */
  charge[i] = charge__;
  pt_min[i] = pt_min__; /* [GeV]    */
  return 1;
L10:return 0;
}

int i_epp__ (int i, double *be, double *mass)
{
  *be = 1.;
  *mass = xin[i - 1];
  return 1;
}


int mc_epp__(int i)
{ i--;

  switch(charge[i])
  { 
    case -1: return  -2212;
    case  1: return   2212;
  } 
}


int m_epp__ (int i,int*pString)
{
  void *pscr = NULL;
  i--;
  for (;;)
    {
      char strmen[] = "\050"
      " Incoming particle mass = XXX           "
      " Incoming particle charge = YYY         "
      " |Q^2|max = ZZZ                         "
      " Pt cut of outgoing proton = VVV        ";
      int mode;
      improveStr (strmen, "XXX", "%.10G GeV", xin[i]);
      improveStr (strmen, "YYY", "%d e", charge[i]);
      improveStr (strmen, "ZZZ", "%.8G GeV^2", q2_max[i]);
      improveStr (strmen, "VVV", "%.9G GeV", pt_min[i]);
      menu1 (38, 10, "", strmen, "n_sf_epp", &pscr, &mode);
      switch (mode)
	{
	case 0: return 1;
	case 1: correctDouble (40, 16, "Enter new value ", xin + i, 1);break;
	case 2:	charge[i]*=-1;break;
	case 3: correctDouble (40, 16, "Enter new value ", q2_max +i,1);break;
	case 4: 
	  do {correctDouble (40, 16, "Enter new value ", pt_min +i,1);}
          while((pt_min[0]>sqrt(q2_max[0])&& pt_min[0]>Maxpt_min) ||
                (pt_min[1]>sqrt(q2_max[1])&& pt_min[1]>Maxpt_min));
                           
	  /* kinematic and q2 approximation limitation */
	}
    }
    return 1;
}

static double phi(double x, double qq)
{
  double a = 7.16;
  double b = -3.96;
  double c = .028;
  double y,qq1,f;
  qq1=1+qq;
  y= x*x/(1-x);
  f=(1+a*y)*(-log(qq1/qq)+1/qq1+1/(2*qq1*qq1)+1/(3*qq1*qq1*qq1));
  f+=(1-b)*y/(4*qq*qq1*qq1*qq1);
  f+=c*(1+y/4)*(log((qq1-b)/qq1)+b/qq1+b*b/(2*qq1*qq1)+b*b*b/(3*qq1*qq1*qq1));
  return f;
}


double c_epp__ (int i, double x, double q)
{
#define alpha .0072992701
#define qz 0.71
  i--;
  { 
    double f, qmi, qma;
    
    qma=q2_max[i]/qz;
/*     x = omega/E = (E-E')/E  ; E,E' - incoming and outgoing protons energy
                              omega = E-E' - energy of emitted photon
*/
    qmi= xin[i]*xin[i]*x*x/(1-x)/qz; 
    qmi+=pt_min[i]*pt_min[i]/(1-x)/qz;

    f = alpha/M_PI*(phi(x,qma)-phi(x,qmi))*(1-x)/x;
    f *= charge[i]*charge[i];
    if (f < 0) f = 0;
    return f;
  }
}

