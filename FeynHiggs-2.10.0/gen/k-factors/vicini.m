<< FormCalc`

foo = Log;
ifoo = Exp;
(*foo = ifoo = Identity;*)

Cfac[sqrts_, mh_, ktopNLO_, kbotReNLO_, kbotImNLO_, ktopNNLO__] :=
  N[{Log[sqrts], foo[mh], #}&/@
    {ktopNLO, kbotReNLO, kbotImNLO, ktopNLO ktopNNLO}]

Ctitles = {"cTopNLO", "cBotReNLO", "cBotImNLO", "cTopNNLO"};

Kfac[sqrts_, mh_, ktopNNLO_, ktopNLO_, kbotNLO_, ktopbotNLO_, __] :=
  N[{Log[sqrts], foo[mh], #}&/@
    {ktopNLO, ktopNNLO, kbotNLO, ktopbotNLO}]

Ktitles = {"kTopNLO", "kTopNNLO", "kBotNLO", "kTopBotNLO"};


get[f_, k_, t_] := MapThread[Rule, {t, get[f, k]}]

get[f_, k_] := Transpose[k@@@ get[f]]

(*get[f_] := ReadList[f, Number, RecordLists -> True]*)
get[f_] := DeleteCases[ToExpression[Import[f]], {}]

data[gghcTev] = get["vicini-data/cf-tev.dat", Cfac, Ctitles]

data[gghkTev] = get["vicini-data/kf-tev.dat", Kfac, Ktitles]

data[gghcLHC] = get["vicini-data/cf-lhc.dat", Cfac, Ctitles]

data[gghkLHC] = get["vicini-data/kf-lhc.dat", Kfac, Ktitles]



poly[vars__, n_Integer, s_String:"a"] :=
  MapIndexed[#1 ToSymbol[s, #2]&,
    Expand[Sum[Plus[vars]^i, {i, 0, n}]]]

fitfun[gghkLHC, "kTopBotNLO"] =
  poly[logm, logsqrts, 8]/(logm - foo[590])^3;

fitfun[gghcTev, "cBotReNLO"] =
  poly[logm, 8];

fitfun[gghcTev, "cBotImNLO"] =
  poly[logm, 6]/(logm - foo[320] + .0001 I)^2;

fitfun[gghcLHC | gghkLHC, _] = IndexIf[
  ifoo[logm] < 348, poly[logm, logsqrts, 8],
  poly[logm, logsqrts, 8, "b"]
] + mx Exp[logm] + sx Exp[logsqrts];

(*
fitfun[gghcLHC, _] = poly[logm, logsqrts, 8] +
  mx Exp[logm] + sx Exp[logsqrts];
*)

fitfun[gghcTev, _] = poly[logm, 8] + mx Exp[logm];

fitfun[gghkTev, _] = poly[logm, 8] + mx Exp[logm];


fitvars = {logsqrts, logm}

fit[li__][title:_ -> data_] :=
Block[ {fun},
  DATA[li, title] = data;
  FUN[li, title] = fun = dofit[Re[data], fitfun[li, title]];
  Print["fitfun = ", fun (*//Shallow*)];
  title -> fun
]


overlap = 10 (*10*)

dofit[d_, IndexIf[x_ < y_, a_, b__] + c_.] := IndexIf[ x < y,
  dofit[Select[d, ifoo[#[[2]]] <= y + overlap &], a + c],
  dofit[Select[d, ifoo[#[[2]]] >= y - overlap &], IndexIf[b] + c] ]

Off[FindFit::lstol, FindFit::sszero]

dofit[d_, fitfun_] :=
Block[ {thefit},
  Check[
    thefit = FindFit[d, fitfun,
      Complement[
        Cases[fitfun, x_Symbol /; Context[x] =!= "System`", {-1}],
        fitvars ],
      fitvars],
    Interrupt[] ];
  fitfun /. thefit
]


fplot[li__] := ffplot[li]/@ data[li]

ffplot[li__][title_ -> data_] := PLOT[li, title] = MapThread[
  Show[#2, DisplayFunction -> $DisplayFunction,
    PlotLabel -> title <> #1]&,
  { {"-Re", "-Im"},
    Transpose[logfffplot[fit[li][title -> data]]/@
      Split[data, #1[[1]] === #2[[1]] &]] } ]

fffplot[_ -> f_][li_] := theplot[
  {ifoo[#2], ##3}&@@@ li,
  f /. {logsqrts -> li[[1, 1]], logm -> foo[m]},
  {m, 75, 1050} ]

logfffplot[_ -> f_][li_] := theplot[
  {##2}&@@@ li,
  f /. logsqrts -> li[[1, 1]],
  {logm, foo[75], foo[1050]} ]

(* logfffplot[a_][b__] := (A=a; B={b}; Interrupt[]) *)

theplot[d_, f_, range_] := {
  theplot[d, f, range, Re],
  theplot[d, f, range, Im] }

theplot[d_, f_, range_, re_] := {
  ListPlot[{#1, re[#2]}&@@@ d,
    PlotStyle -> {Red},
    (*PlotRange -> All,*)
    DisplayFunction -> Identity],
  Plot[re[f], range,
    PlotStyle -> {Blue},
    PlotRange -> {Automatic, {Min[#], Max[#]}&[re[Last/@ d]]},
    DisplayFunction -> Identity] }


MkDir["f"];
MkDir["m"];

$DebugCmd = {"#ifdef DETAILED_DEBUG\n", "#endif\n", "DPROD ", " ENDL"};
SetOptions[PrepareExpr, DebugLines -> True];

fort[li_] :=
Block[ {file = ToString[li]},
  fplot[li];
  fits = fit[li]/@ data[li] /. {E^logm -> m, E^logsqrts -> sqrts};
  fits = fits /.
    (v_ -> IndexIf[c_, a_, b_]) -> IndexIf[c, v -> a, v -> b] //.
    {x___, IndexIf[c_, a1_, b1_], IndexIf[c_, a2_, b2_], y___} :>
      {x, IndexIf[c, {a1, a2}, {b1, b2}], y} /.
    (v_ -> a_) :> (v -> Re[a]) /; !FreeQ[a, Complex];
  Put[fits, ToFileName["m", file <> ".m"]];
  hh = OpenFortran[ToFileName["f", file <> ".f"]];
  WriteExpr[hh, fits, HornerStyle -> True];
  Close[hh]
]


fct := fort[gghcTev]

fkt := fort[gghkTev]

fcl := fort[gghcLHC]

fkl := fort[gghkLHC]

fall := {fct, fkt, fcl, fkl}

