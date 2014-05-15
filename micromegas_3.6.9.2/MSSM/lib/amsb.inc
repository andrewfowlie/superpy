extern double xxxAMSBc(double m0,double m32, double tb, double sgn);

int xxxAMSB(double m0,double m32, double tb, double sgn)
{
  int err;
  err=xxxAMSBc(m0,m32,tb,sgn);
  assignValW("tb",tb);
  if(err==0) FillVal(1);
  return err;
}

int  xxxamsb_(double* m0,double* m32, double* tb, double* sgn)
{ return  xxxAMSB(* m0,* m32, *tb, *sgn  ); }
