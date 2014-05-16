
procedure initSum;
begin

if( (length(inParticles) neq 2) or (length(inParticles) neq 2)) then
begin
   write " Current summation is for 2->2 processes only";
   bye;
end;
 

sum:=0;

let p1.p2= ( s -p1.p1 -p2.p2)/2;
let p1.p3= (-t +p1.p1 +p3.p3)/2;

operator sp,up,tp,up_,tp_;
 
for all m,w let propDen(-p1-p2,m,w)= -sp(m);
for all m,w let propDen( p1+p2,m,w)= -sp(m);


for all m,w let propDen(-p1+p3,m,w)= -tp(m);
for all m,w let propDen( p1-p3,m,w)= -tp(m);

for all m,w let propDen(-p2+p3,m,w)= -up(m);  
for all m,w let propDen( p2-p3,m,w)= -up(m);  

end;


Procedure addToSum;
begin
  scalar Ans, Ans2;
  Ans:=totFactor*numerator*denominator;
  if strcmp( PART(outParticles,1),PART(outParticles,2)) = t then
  begin
     Ans_2:= (Ans where {up(~x)=>up_(~x), tp(~y)=>tp_(~y),t=>t_});
     Ans_2:= (Ans_2 where { up_(~x)=>tp(~x),
              tp_(~y)=>up(~y),t_=>(p1-p4).(p1-p4)});
     Ans:=(Ans+Ans_2)/2;
     Clear Ans_2;
  end$

  Ans:=(Ans where
     s*sp(~m)=> 1+m^2*sp(m),
     t*tp(~m)=> 1+m^2*tp(m),
     t*up(~m)=>-1-(~m^2+s-p1.p1-p2.p2-p3.p3-p4.p4)*up(m),
     s*tp(~m1)*up(~m2)=>(p1.p1+p2.p2+p3.p3+p4.p4 -m1^2-m2^2)*tp(m1)*up(m2)
                        -tp(m1) -up(m2),
     tp(~m1)*tp(0)*m1^2=>(tp(m1) -tp(0)),
     up(~m1)*up(0)*m1^2=>(up(m1) -up(0)),
     sp(~m1)*sp(0)*m1^2=>(sp(m1) -sp(0)),
     substitutions);
     sum:=sum+Ans;
end;

symbolic procedure strcmp(s1,s2);
 explode s1=explode s2;
symbolic operator strcmp;


procedure finishSum(); begin  end;


END;
