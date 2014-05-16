
procedure initSum;
begin 
   sum:=0;
   substitutions_:=append({propDen(~p,~m,~w)=>(m^2-p.p)},substitutions);
end;

procedure addToSum;
begin
  scalar nnn;
  sum:=sum+totFactor*numerator/denominator;
  nnn:=num(sum);
  sum:=(nnn  where substitutions_)/den(sum);
end;


symbolic procedure strcmp(s1,s2);
 explode s1=explode s2;
symbolic operator strcmp;


procedure POut_(j); begin  return mkid(p,j+length(inParticles)); end;
procedure sub_(A,B);begin  return sub(A,B); end;

procedure finishSum;
begin
  scalar i,j,nnn;
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
          nnn:=num(sum);
          sum:=(nnn where substitutions_ )/den(sum);  
                                            
        j:=j+1;
     end;  
     sum:=sum/(j-i);
     i:=j; 
  end;
  off nat;
end;

end;
