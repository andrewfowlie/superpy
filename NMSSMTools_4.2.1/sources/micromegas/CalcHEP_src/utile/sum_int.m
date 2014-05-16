
initSum:= Module[
   {MM1,MM2,MM3,MM4,tMid}, 

   If[Length[inParticles] !=2 || Length[outParticles] !=2,
   Print["Current summation is for 2->2 processes only"];
   Quit]; 
 
   sum=0;

   MM1=(SC[p1,p1] /. substitutions) ;
   MM2=(SC[p2,p2] /. substitutions) ;
   MM3=(SC[p3,p3] /. substitutions) ;
   MM4=(SC[p4,p4] /. substitutions) ;

   SC/: SC[p1,p2]=( s -MM1 -MM2)/2;
   SC/: SC[p1,p3]=(-t +MM1 +MM3)/2;

   E1cm= (s+MM1-MM2)/(2*Sqrt[s]); 
   E2cm= (s+MM2-MM1)/(2*Sqrt[s]);
   E3cm= (s+MM3-MM4)/(2*Sqrt[s]);
   E4cm= (s+MM4-MM3)/(2*Sqrt[s]);

   sumMass=MM1+MM2+MM3+MM4;

   beSquared =Simplify[16*(E1cm^2-MM1)*(E3cm^2-MM3)/s^2];
   subS=0;

   If[MM1*MM2==0 && MM3*MM4==0, 
      be$=Simplify[(s-MM1 -MM2)*(s-MM3 -MM4)/s^2],
      If[MM1+MM2==MM3+MM4 && MM1*MM2==MM3*MM4, 
        be$=Simplify[((s-MM1 -MM2)^2 - 4*MM1*MM2)/s^2], 
        Clear be$;        
        IF[MM1==MM2 && MM3==MM4 && (MM1==0 || MM3==0), 
           subS=Simplify[2*sumMass/(1-be$^2)]] 
      ] 
   ];

   tMid=(MM1-MM3-MM2+MM4)^2/(4*s)-(E1cm^2-MM1)-(E3cm^2-MM3);
   tMax=Simplify[tMid+be$*s/2];
   tMin=Simplify[tMid-be$*s/2];
   propDen[p_,m_,w_]:=((m^2- SC[p,p]) /. substitutions);
];


addToSum:=Module[
   {intFun,intFunL,intFun0}, 
   intFun=Integrate[Simplify[numerator/denominator /. substitutions] ,t];
   intFun0= intFun /. Log[x_] ->  0;
   intFunL= Simplify[intFun- intFun0];
  
   intFun0=Simplify[(intFun0 /. t->tMax)-(intFun0 /. t->tMin)];

   sum=sum + totFactor* Simplify[intFun0 + (intFunL /. Log -> LLog)];
];


Ordp[v1_,v2_]:=Module[
   {i,Res,s1,s2},
   s1=ToCharacterCode[ToString[InputForm[v1]]];
   s2=ToCharacterCode[ToString[InputForm[v2]]]; 

   For[i=1,i<=Length[s1]&&i<=Length[s2]&&Part[s1,i]==Part[s2,i],i++];

   If[i>Length[s1],Res=False,
      IF[i>Length[s2],Res=True,
        If[Part[s1,i]>Part[s2,i], Res=True,Res=False];
      ];
   ];
   Res
];


LLog[x_]:=Module[
   {lArg,lArg1,Ans}, 

   lArg= (x /. t->tMax)/(x /. t->tMin);
   If[subS !=0, lArg=(lArg /. s->subS)];
   lArg =Simplify[Together[lArg]]; 
   lArg1=Simplify[Together[1/lArg]];
   If[Ordp[lArg,lArg1],Ans=-mLog[lArg1], Ans= mLog[lArg]];
   Ans
];


finishSum:=(
   sum=Simplify[sum];
   sum=sum /. mLog -> Log;
   sum=(sum/((s- SC[p1,p1]-SC[p2,p2])^2 - 4*SC[p1,p1]*SC[p2,p2])/(16*pi));
);
