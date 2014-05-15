#include "syst.h"
#include "polynom.h"
#include "tensor.h"
#include "spinor.h"

/*#include "reader_s.h"*/
#include "parser.h"
#include "sets.h"

#include "symb_tot.h"

/* Problem expected in the case S1==S2 */

void symb_clean(symb_data S)
{
  switch(S.type)
  {
     case errortp  : break;
     case numbertp :
     case polytp   : delpoly(&(S.expr.p));   break;
     case vectortp : 
     case indextp  : 
     case tenstp   : deltensor(&(S.expr.t)); break;
     case spintp   : delSpin(S.expr.s);      break;
     case etenstp  : delEtens(S.expr.et);    break;
  }               
}

symb_data symb_copy(symb_data S)
{
  symb_data ans;
  
  ans.type=S.type;

  switch(S.type)
  {
     case errortp  : ans.expr.p =NULL; break; 
     case numbertp :
     case polytp   : ans.expr.p =copypoly(S.expr.p); break;
     case vectortp : 
     case indextp  : 
     case tenstp   : ans.expr.t =copytens(S.expr.t);   break;
     case spintp   : ans.expr.s =copySpin(S.expr.s);   break;
     case etenstp  : ans.expr.et=copyEtens(S.expr.et); break;
  }
  return ans;               
}


symb_data symb_imult(symb_data S,int del,int factor)
{
  symb_data ans;
  if(S.type==errortp||S.type==indextp)
  { ans.type=errortp;ans.expr.p=NULL; return ans;}

  if(!del) ans=symb_copy(S); else { ans.type=S.type; ans.expr.p=S.expr.p; }

  switch(S.type)
  {
     case numbertp :
     case polytp   : multpolyint(&(ans.expr.p),factor); break;
     case vectortp : 
     case tenstp   : multtensint(&(ans.expr.t),factor); break;
     case spintp   : multSpinInt(&(ans.expr.s),factor); break;
     case etenstp  : multEtensInt(&ans.expr.et,factor); break;
  } 
  return ans;
}


symb_data symb_mult(symb_data S1, int del1, symb_data S2, int del2)
{
  symb_data ans={{NULL},errortp};
  
  if(S1.type==errortp || S2.type==errortp ||
     S1.type==indextp || S2.type==indextp) return ans;

  if(S1.expr.p==NULL)
  { if(del2) symb_clean(S2); S1.type=numbertp;  return  S1;}

  if(S2.expr.p==NULL)
  { if(del1) symb_clean(S1); S2.type=numbertp;  return  S2;}

 
  if(S1.type==S2.type)
  {  
     if(S1.expr.p==S2.expr.p)
     { S2=symb_copy(S2); 
       del1=(del1||del2);
       del2=1;
     } 

     switch(S1.type)
     { 
        case numbertp : 
        case polytp   : ans.expr.p=multtwopoly(S1.expr.p,S2.expr.p);   break;
        case tenstp   : ans.expr.t=multtwotens(S1.expr.t,S2.expr.t);   break;
        case spintp   : ans.expr.s=mult2Spin(S1.expr.s,S2.expr.s,0);   break;  
        case etenstp  : ans.expr.et=mult2Etens(S1.expr.et,S2.expr.et); break;
        default       : return ans;
     }
     ans.type=S1.type;
     if(del1) symb_clean(S1);
     if(del2) symb_clean(S2);
     
     return ans;
  } 

  
  if(S1.type< S2.type)
  {  symb_data S=S1;
     int del=del1;
     S1=S2; S2=S;
     del1=del2; del2=del;
  }

  if(S2.type==numbertp)
  { 
    long L=S2.expr.p->num;
    ans=symb_imult(S1,del1,L);
    if(del2) symb_clean(S2);
    return ans;
  }

  if(S2.type==vectortp||S2.type==spintp) return ans;
  
  if(del1) ans=S1; else ans=symb_copy(S1);
  
  switch(S1.type)
  {      
     case vectortp : 
     case tenstp   : multtenspoly(&ans.expr.t,S2.expr.p); break;
     case spintp   :
        switch(S2.type)
        { 
           case polytp   :multSpinPoly(&ans.expr.s,S2.expr.p); break;
           case tenstp   :multSpinTens(&ans.expr.s,S2.expr.t); break;
        } break;
     case etenstp  :
        switch(S2.type)
        { 
           case polytp   :multEtensPoly(&ans.expr.et,S2.expr.p); break;
           case tenstp   :multEtensTens(&ans.expr.et,S2.expr.t); break;
        } break;
  }
  
  if(del2) symb_clean(S2);
  return ans; 
}


symb_data  symb_typeUp(symb_data S, int del,int type)
{
   symb_data ans={{NULL},errortp};
   
   if(S.type==errortp || type<S.type) return ans;

   if(S.type==type){ if(del) ans=S; else ans=symb_copy(S); return ans; }

   switch(S.type)
   { case vectortp: case indextp: case spintp: case etenstp: return ans; }
   
   switch(type)
   { 
     case vectortp :
     case indextp  : return ans;
     case polytp   : if(del) ans=S; else ans=symb_copy(S); 
                     ans.type=polytp;
                     return ans;
     case tenstp   : ans.expr.t=newtensor1();
                     multtenspoly(&(ans.expr.t),S.expr.p); break;                     
     case spintp   : ans.expr.s=spin1(); 
                     if(S.type==tenstp) multSpinTens(&ans.expr.s,S.expr.t);
                     else               multSpinPoly(&ans.expr.s,S.expr.p);
                     break;
     case etenstp  : ans.expr.et=etens1();    
                     if(S.type==tenstp) multEtensTens(&ans.expr.et,S.expr.t);
                     else               multEtensPoly(&ans.expr.et,S.expr.p);
                     break;
   }

   ans.type=type; 
   if(del) symb_clean(S);   
   return ans;
}


static set symb_index(symb_data S)
{
  set index=set_constr(_E);
  tensor t; SpinTensor s; Etens e;
  int i,c;



  switch(S.type)
  {  case errortp: case numbertp: case polytp: case vectortp: case indextp: 
     return index;
  }

  if(S.type==errortp || S.expr.p==NULL)  return index;

  switch(S.type)
  {  case tenstp : for(i=0,t=S.expr.t; i<maxIndex; i++) 
                   if(t->tens[i]) set_add1(&index,i+1);
                   
                   break;

     case spintp : for(i=0,s=S.expr.s; i<s->l; i++) 
                   {  c=s->g[i];
                      if(c>0)  set_add1(&index,c);
                   }
                   for(i=0,t=s->tcoef; i<maxIndex; i++) 
                   if(t->tens[i]) set_add1(&index,i+1);
                   break;

     case etenstp: for(i=0,e=S.expr.et; i<4; i++) 
                   {  c=e->eps[i];
                      if(c>0)  set_add1(&index,c);
                   }
                   t=e->tcoef;
                   for(i=0; i<maxIndex; i++)  if(t->tens[i]) set_add1(&index,i+1);
  }
  return index;
} 


int symb_testindex(symb_data S)
{
  set index=set_constr(_E);
  tensor t; SpinTensor s; Etens e;
  int i,c;
  int first=1;


  switch(S.type)
  {  case errortp: case numbertp: case polytp: case vectortp: case indextp: 
     return 0;
  }

  if(S.type==errortp || S.expr.p==NULL)  return 0;
  
  switch(S.type)
  {  case tenstp : 
        
       for(t=S.expr.t;t;t=t->next)  
       {set   ind=set_constr(_E); 
         for(i=0; i<maxIndex; i++) if(t->tens[i]) set_add1(&ind,i+1);
         if(first) index=ind; else  if(!set_eq(index,ind)) return 1; 
         first=0;         
       }   return 0;

     case spintp :
     for(s=S.expr.s;s;s=s->next)
     { set  ind=set_constr(_E);
       for(i=0; i<s->l; i++) {  c=s->g[i];if(c>0)  set_add1(&ind,c);}
       t=s->tcoef;            
       for(i=0; i<maxIndex; i++)  if(t->tens[i]) set_add1(&ind,i+1);
       if(first) index=ind; else  if(!set_eq(index,ind)) return 1;    
       first=0;
     } return 0;                

     case etenstp:
     for(e=S.expr.et;e;e=e->next)
     { set  ind=set_constr(_E);
          
       for(i=0;i<4;i++) {  c=e->eps[i]; if(c>0)  set_add1(&ind,c);}
       t=e->tcoef;
       for(i=0; i<maxIndex; i++)  if(t->tens[i]) set_add1(&ind,i+1);
       if(first) index=ind; else  if(!set_eq(index,ind)) return 1; 
       first=0;
    } return 0;   
      
  }
  return 1;   
} 




symb_data symb_sum(symb_data S1, int del1, symb_data S2, int del2)
{
  symb_data ans={{NULL},errortp}, s1,s2;
  set index1,index2;
    
  if( S1.type==errortp || S2.type==errortp  
    ||S1.type==indextp || S2.type==indextp) return ans; 
 
  if(S1.expr.p==NULL)
  { if(del2) return S2; else return symb_copy(S2);}

  if(S2.expr.p==NULL)
  { if(del1) return S1; else return symb_copy(S1);}

 
  index1= symb_index(S1);
  index2= symb_index(S2); 
  
  if(!set_eq(index1,index2)) 
  {  
      printf("index problem\n");
      set_print(index1);
       set_print(index2);
        return ans;
         
  }

  if(S1.expr.p==S2.expr.p)
  { S2=symb_copy(S2); 
    del1=(del1||del2);
    del2=1;
  }

  if(S1.type<S2.type)
  {  symb_data S=S1;
     int del=del1;
     S1=S2; S2=S;
     del1=del2; del2=del;
  }
  
  if(S1.type>S2.type)                      
  {    
    s2=symb_typeUp(S2,del2,S1.type);
    if(s2.type==errortp) return ans; 
  } else 
  if(del2) s2=S2; else s2=symb_copy(S2);   
  if(del1) s1=S1; else s1=symb_copy(S1); 

  switch(s1.type)
  {      
     case numbertp : 
     case polytp   : sewpoly(&s1.expr.p,&s2.expr.p);   break;         
     case vectortp :
     case tenstp   : sewtens(&s1.expr.t,&s2.expr.t);   break;
     case spintp   : addSpin(&s1.expr.s,s2.expr.s);    break;
     case etenstp  : addEtens(&s1.expr.et,s2.expr.et); break;
   } 
 
   return s1;
}

symb_data symb_spur(symb_data S, int del)
{
   symb_data SS=symb_typeUp(S, del,spintp);
   Etens EE=calcspur(SS.expr.s);
   delSpin(SS.expr.s);
   SS.expr.et=EE;
   SS.type=etenstp;
   return SS;
}

int symb_iszero(symb_data S)
{
  if(S.type !=errortp && S.expr.p==NULL) return 1; else return 0; 
}

symb_data symb_real(symb_data S, int del)
{
  symb_data SS;
  if(del) SS=S; else SS=symb_copy(S);
  
  switch(SS.type)
  { case errortp: 
    case numbertp:
    case   polytp: return SS; 

    case vectortp:
    case tenstp:  tensRealPart(&(SS.expr.t)); return SS;
    case spintp:
          { SpinTensor s=SS.expr.s,pred=NULL,ss;
            while(s) 
            {
              if(s->g5) tensImPart(&(s->tcoef)); else tensRealPart(&(s->tcoef));
              ss=s;
              s=s->next;
              if(ss->tcoef==NULL) 
              { delunit(ss); 
                if(pred==NULL) SS.expr.s=s; else pred->next=s;
              } else pred=ss;
            }
          }  
          return SS;
    case etenstp:
          { Etens s=SS.expr.et,pred=NULL,ss;
            while(s) 
            {
              tensRealPart(&(s->tcoef));
              ss=s;
              s=s->next;
              if(ss->tcoef==NULL) 
              { delunit(ss); 
                if(pred==NULL) SS.expr.et=s; else pred->next=s;
              } else pred=ss;
            }
          }    
          return SS;
  }                      
}

symb_data symb_delEps(symb_data S, int del)
{
  symb_data SS;
  if(del) SS=S; else SS=symb_copy(S);
  
  switch(SS.type)
  { case errortp: 
    case numbertp:
    case polytp:
    case vectortp:
    case tenstp:  return SS;
    case spintp:
          { SpinTensor s=SS.expr.s,pred=NULL,ss;
            while(s) 
            {
              ss=s;
              s=s->next;
              if(ss->g5) 
              { 
                 deltensor(&(ss->tcoef)); delunit(ss); 
                 if(pred==NULL) SS.expr.s=s; else pred->next=s;
              } else pred=ss;
            }
          }  
          return SS;
    case etenstp:
          { Etens s=SS.expr.et,pred=NULL,ss;
            while(s) 
            { ss=s;
              s=s->next;
              if(ss->eps[0]!=X_MARK)
              { deltensor(&(ss->tcoef));
                delunit(ss); 
                if(pred==NULL) SS.expr.et=s; else pred->next=s;
              } else pred=ss;
            }
            ss=SS.expr.et;
            SS.type=tenstp;
            if(ss) {SS.expr.t=ss->tcoef;
            delunit(ss);}
          }    
          return SS;
  }                      
}
