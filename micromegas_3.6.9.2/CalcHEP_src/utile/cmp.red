
symbolic procedure strcmp(s1,s2);
 explode s1=explode s2;
symbolic operator strcmp;


firstRun:=1;

procedure initSum;
begin
  if(firstRun=1) then 
  begin  
     diff_:=0; 

     if(length(inParticles)=2) and (length(outParticles)=2) then case_22:=1 
                                                            else case_22:=0;
     inParticles_cp:=inParticles;
     outParticles_cp:=outParticles; 
  end 
  else 
  begin
     if not ( (inParticles=inParticles_cp) and
             (outParticles=outParticles_cp) ) then 
     begin    
        out   "message"$
        write "Attempt to compare processes of different types"$
        shut  "message"$
        out T;
     end;
     diff_:=-diff_;
  end;

  if(case_22=1) then initSum22() else initSumDe();

end;

procedure addToSum;
begin 
  if case_22=1 then addToSum22() else addToSumDe();
end;

procedure finishSum();
begin
   if case_22=0 then diff_:=(diff_  where substitutions);
   if(firstRun=0) then
   begin
      if case_22=1 then 
        diff_:=(diff_ where tp(~m)=>1/(t-m^2),
                            sp(~m)=>1/(s-m^2),
                            up(~m)=>1/(-s-t+p1.p1+p2.p2+p3.p3+p4.p4 -m^2));
      out "message" $
      if( not (diff_ = 0)) then  write "Error"  else  write "Ok" $
      write "Check finished for ", inParticles, "->",outParticles $
      shut "message"$
      out T$
   end$
   firstRun:=0$
end$


%Default procedures

procedure AnsWhere(x,y);
begin
  Ans_:=(Ans_ where x=>y);    
end;

procedure POut_(j); begin  return mkid(p,j+length(inParticles)); end;

procedure initSumDe; 
begin 
   substitutions_:=append({propDen(~p,~m,~w)=>(m^2-p.p)},substitutions);
end;

procedure addToSumDe;
begin
  scalar Ans1,diff_1,i,j;
  vector q_;

  Ans_:=totFactor*numerator/denominator;
  i:=1;
  while i<length(outParticles) do
  begin
     j:=i+1;
     while (j<=length(outParticles)) and
            strcmp(part(outParticles,i),part(outParticles,j)) do
     begin
        Ans1:=Ans_;
        AnsWhere(POut_(i),q_);
        AnsWhere(q_,POut_(j));
        if not (j = length(outParticles)) then  AnsWhere(POut_(j),POut_(i));
        Ans_:= Ans1+Ans_;
        j:=j+1;
     end;
     Ans_:=Ans_/(j-i);
     i:=j;
  end;

  diff_:=diff_+ Ans_;
  clear Ans_;
  diff_1:=num(diff_);
  diff_:=(diff_1 where substitutions_)/den(diff_);
end;


% case 2->2 procedures

procedure initSum22; 
begin 
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


procedure addToSum22;
begin
  scalar Ans,Ans2;
  Ans:=totFactor*numerator*denominator;
  if strcmp( PART(outParticles,1),PART(outParticles,2)) = t THEN
  begin
     Ans_2:= (Ans where {up(~x)=>up_(~x), tp(~y)=>tp_(~y),t=>t_});
     Ans_2:= (Ans_2 where { up_(~x)=>tp(~x), tp_(~y)=>up(~y),t_=>(p1-p4).(p1-p4)});
     Ans:=(Ans+Ans_2)/2;
  end;
  Ans:=(Ans where 
     s*sp(~m)=> 1+m^2*sp(m),
     t*tp(~m)=> 1+m^2*tp(m),
     t*up(~m)=>-1-(~m^2+s-p1.p1-p2.p2-p3.p3-p4.p4)*up(m),
     s*tp(~m1)*up(~m2)=>(p1.p1+p2.p2+p3.p3+p4.p4 -~m1^2 -~m2^2)*tp(~m1)*up(~m2)
                        -tp(~m1) -up(~m2),
     tp(~m1)*tp(~m2)=>(tp(m1) -tp(m2))/(m1^2-m2^2),
     up(~m1)*up(~m2)=>(up(m1) -up(m2))/(m1^2-m2^2),
     sp(~m1)*sp(~m2)=>(sp(m1) -sp(m2))/(m1^2-m2^2), 
     substitutions);
                     
  diff_:=diff_+Ans;  
end;


% reading off results

in"results/symb1.red"$
in"results_/symb1.red"$

End$     
