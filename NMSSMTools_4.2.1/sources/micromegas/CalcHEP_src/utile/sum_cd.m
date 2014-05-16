
initSum:=(sum=0;);

propDenSub[p_,m_,w_]:=((m^2-SC[p,p]));

addToSum:=(
   sum=Together[(sum+(totFactor*numerator/denominator) /. substitutions)];
   sum=Simplify[((Numerator[sum] /. propDen->propDenSub) /. substitutions)]
       /Denominator[sum];
);


POut$[j_]:=Module[
  {i0},
  i0=Length[inParticles];
  ToExpression["p"<>ToString[j+i0]]
];


finishSum:=Module[
  {i,j,q},
  For[i=1, i<Length[outParticles], i=j,
     For[j=i+1, j<=Length[outParticles] &&
                Part[outParticles,i]==Part[outParticles,j], j++,

        If[j == Length[outParticles], 
          sum=sum+(sum /.  POut$[i] ->POut$[j])
        ,
          sum=sum+(((sum /. POut$[i] -> q) /. POut$[j] -> POut$[i]) /. 
                                                      q-> POut$[j]);
        ];

        sum=Simplify[(Numerator[sum] /. propDen->propDenSub) /. substitutions]
            /Denominator[sum];        
     ];
     sum=sum/(j-i);
  ];
  sum=Simplify[sum];
];
