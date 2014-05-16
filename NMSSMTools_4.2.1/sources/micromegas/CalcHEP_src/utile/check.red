
PathToResults:="results";

procedure in_(r)$ << in r$ >>$

symbolic procedure concatstr(x,y)$
compress  append( reverse cdr reverse explode x,cdr explode y)  $
symbolic operator concatstr$

procedure initSum(); begin  end;

procedure addToSum();
begin 
     FileName:=concatstr( concatstr(PathToResults,"/p"),
                       concatstr(DiagrNumber,".red"))$
     in_ FileName$

     diff_:=totFactor_*numerator_*denominator/denominator_
            -totFactor*numerator$    
%    Clear numerator_,numerator,totFactor_,totFactor, denominator_,denominator$       

     out  "message"$
     if (diff_=0) then  write DiagrNumber,  " OK" else
      Write  DiagrNumber, " Error " $
     out T$
     if (diff_=0) then  write DiagrNumber,  " OK" else
     Write  DiagrNumber, " Error " $
end$

procedure finishSum(); 
begin 
  out "message"$
  write "Check finished for ", inParticles, "->",outParticles $
  shut "message"$
  out T$
end;

in"results/symb1.red"$

End$     
