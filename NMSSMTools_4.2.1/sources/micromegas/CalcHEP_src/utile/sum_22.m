initSum:=(

   If[Length[inParticles] !=2 || Length[outParticles] !=2,
   Print["Current summation is for 2->2 processes only"];
   Quit ];

   MM$1=SC[p1,p1] /. substitutions;
   MM$2=SC[p2,p2] /. substitutions;
   MM$3=SC[p3,p3] /. substitutions;
   MM$4=SC[p4,p4] /. substitutions;

   massSum=Simplify[MM$1+MM$2+MM$3+MM$4];
   sum=0;
);

addToSum:=Module[
   {sqm},
   sqm=totFactor*numerator/denominator /.  propDen[p_,m_,w_]->(m^2-SC[p,p]);
   sqm=sqm /. substitutions;

   sqm=sqm /. {SC[p1,p2]->( s - MM$1 -MM$2)/2 ,
               SC[p1,p3]->(-t + MM$1 +MM$3)/2};            
   sum=sum+sqm;
];

finishSum:=Module[
   {sumU,u},
   If[Part[outParticles,1]==Part[outParticles,2], 
    sumU=sum /. t-> u;
    sum=(sum+ (sumU /. u-> (massSum -s-t)))/2; 
   ];
   sum=Apart[sum];
];
