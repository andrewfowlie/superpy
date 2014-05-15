
procedure initSum();
begin
   if( (length(inParticles) neq 2) or (length(inParticles) neq 2)) then
   begin
      write " Current summation is for 2->2 processes only";
      bye;
   end; 
 
   sum:=0;

   korder s;
   order s;

   let abs(s)=s;
   for all x let abs(x)^2=x^2;

   let p1.p2= ( s -p1.p1 -p2.p2)/2;
   let p1.p3= (-t +p1.p1 +p3.p3)/2;


   E1_cm:= (s+p1.p1-p2.p2)/(2*sqrt(s));
   E2_cm:= (s+p2.p2-p1.p1)/(2*sqrt(s));
   E3_cm:= (s+p3.p3-p4.p4)/(2*sqrt(s));
   E4_cm:= (s+p4.p4-p3.p3)/(2*sqrt(s));

   sum_mass:=p1.p1+p2.p2+p3.p3+p4.p4;
   sub_S:=0;
   if (p1.p1*p2.p2=0 and p3.p3*p4.p4 = 0) then 
      be_:= (s-(p1.p1) -(p2.p2))*(s-(p3.p3) -(p4.p4))/s^2
   else if  (p1.p1+p2.p2 =p3.p3+p4.p4 and p1.p1*p2.p2 =p3.p3*p4.p4) then
      be_:= ((s-p1.p1 -p2.p2)^2 - 4*p1.p1*p2.p2)/s^2 
   else  
   begin 
      beSquared:=16*(E1_cm^2-(p1.p1))*(E3_cm^2-(p3.p3))/s^2;
      Clear be_;
      if(p1.p1=p2.p2 and p3.p3=p4.p4 and (p1.p1=0 or p3.p3=0))
      then   sub_S:=2*sum_mass/(1-be_^2);
   end; 

   t_max:=(p1.p1-p3.p3-p2.p2+p4.p4)^2/(4*s) -(E1_cm^2 - (p1.p1))-
                                           (E3_cm^2 - (p3.p3)) +be_*s/2;

 
   t_min:=(p1.p1-p3.p3-p2.p2+p4.p4)^2/(4*s) -(E1_cm^2 - (p1.p1))-
                                           (E3_cm^2 - (p3.p3)) -be_*s/2;

   for all m,w let propDen(-p1-p2,m,w)= -1/(s-m^2);
   for all m,w let propDen( p1+p2,m,w)= -1/(s-m^2);

   for all m,w let propDen(-p1+p3,m,w)= -1/(t-m^2);
   for all m,w let propDen(p1-p3,m,w)=  -1/(t-m^2);

   for all m,w let propDen(-p2+p3,m,w)= -1/(sum_mass-s-t-m^2);
   for all m,w let propDen(p2-p3,m,w)=  -1/(sum_mass-s-t-m^2);

end;

procedure s_subst(s,x) ;
begin
   return (s where t=>x);
end; 
              

procedure addToSum();
begin

  int_fun:=int(numerator*denominator ,t);
  int_fun:=(int_fun where substitutions);
  
  sum:=sum+  totFactor* (s_subst(int_fun,t_max)- s_subst(int_fun,t_min));

end$

procedure finishSum();
begin

  on combinelogs;
  sum:=sum;
  off combinelogs;

  if not (sub_S=0) then  sum:=(sum where log(~x)=>log(sub({s=sub_S},x)));
  sum:= (sum where log(~x)=>-log(1/x) when  ordp(x,1/x));
  sum:=(sum/((s- p1.p1-p2.p2)^2 -4*p1.p1*p2.p2)/(16*pi) where substitutions);
end;

end;

