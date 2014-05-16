
procedure initSum;
begin 
   kkk:=1;
   on rounded;
   precision 42;
   sum:=0;
   substitutions_:=append({propDen(~p,~m,~w)=>(m^2-p.p)},substitutions);

%   p1.p2:=8.000000001305593E+02;
%   p1.p3:=2.611210101477468E-07;
%   p1.p4:=3.969580164636039E+02;
%   p2.p3:=7.997770402860247E+02;
%   p2.p4:=1.123277225058524E-01;


p1.p2:=8.000000001306E+02;
p1.p3:=2.611210101477E-07;
p1.p4:=3.969580164636E+02;
%p1.p5:=4.030419836670E+02;
p2.p3:=7.997770402860E+02;
p2.p4:=1.123277225059E-01;
%p2.p5:=1.106321220289E-01;
%p3.p4:=3.968473843416E+02;
%p3.p5:=4.029296559444E+02;
%p4.p5:=2.117873545347E-01;


end;

procedure addToSum;
begin

%  ans:= ( (totFactor*numerator/denominator  where substitutions_) where
%parameters) ;

  ans:=( (  totFactor*numerator  where substitutions_) where parameters) ;


  sum:=sum+ans;

  write kkk,"  ", ans;
%  write ((denominator where substitutions_  ) where parameters);

  kkk:=kkk+1;
end;


symbolic procedure strcmp(s1,s2);
 explode s1=explode s2;
symbolic operator strcmp;


procedure POut_(j); begin  return mkid(p,j+length(inParticles)); end;
procedure sub_(A,B);begin  return sub(A,B); end;

procedure finishSum;
begin
  scalar i,j;
  vector q_;

  i:=1;
  while i<length(outParticles) do
  begin
     j:=i+1;
     while (j<=length(outParticles)) and  
            strcmp(part(outParticles,i),part(outParticles,j)) do
     begin
        if(j eq length(outParticles)) then 
          sum:= sum+ sub_({POut_(i)=POut_(j)},sum)
        else
          sum:= sum+ 
                    sub_({POut_(j)=POut_(i), q_=POut_(j)},
                       sub_({POut_(i)=q_},sum)
                       );
          sum:=(num(sum) where substitutions_ )/den(sum);  
                                            
        j:=j+1;
     end;  
     sum:=sum/(j-i);
     i:=j; 
  end;
  off nat;
end;

in"results/symb1.red"$

write "sum=",sum;

end;
