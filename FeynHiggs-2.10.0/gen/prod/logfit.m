<< FormCalc`

debug = True

plotmax = 2000

get[tag_, file_String] := Select[getf[tag, file],
  !MemberQ[ {0., Indeterminate}, #[[-1]] ] &]

getf[tag_, file_] := getf[tag, file] = log[file, tag]/@
  ReadList[ToFileName["data", file], Number, RecordLists -> True]


foo = Log


log["ststh-lhc14.dat", _][{mh_, mst_, xs_?NumberQ}] :=
  {Sqrt[mh], mst, foo[10^-5 xs]} //N

log["tHm-lhc14.dat", _][{mh_, mt_, tb_, xsLO_, xsNLO_, ratio_}] :=
  {Sqrt[mh], tb, foo[1000 xsNLO]} //N

log["tHm2-lhc8.dat", _][{mh_, tb_, xsNLO_, xsNLOhi_, xsNLOlo_}] :=
  {Sqrt[mh], tb, foo[1000 xsNLO]} //N

log["tHm2-lhc8lo.dat", _][{mh_, tb_, xsNLO_, xsNLOhi_, xsNLOlo_}] :=
  {Sqrt[mh], tb, foo[1000 xsNLOlo]} //N

log["tHm2-lhc8hi.dat", _][{mh_, tb_, xsNLO_, xsNLOhi_, xsNLOlo_}] :=
  {Sqrt[mh], tb, foo[1000 xsNLOhi]} //N

log["bbh-tev.dat" | "bbh-lhc14.dat", _][{mh_, xs_?NumberQ}] :=
  {Sqrt[mh], foo[1000 xs]} //N

log["bbh-lhc.dat", _][{sqrts_, mh_, xs_?NumberQ}] :=
  {Sqrt[mh], sqrts, foo[1000 xs]} //N

log["ggh-tev.dat", _][{sqrts_, mh_, xsNNLO_, xsNNLL_, __}] :=
  {Sqrt[mh], foo[1000 xsNNLL]} //N

log["ggh-lhc.dat", _][{sqrts_, mh_, xsNNLO_, xsNNLL_, __}] :=
  {Sqrt[mh], Log[sqrts/1000], foo[1000 xsNNLL]} //N

log["ggh-smnlo-tev.dat", _][{sqrts_, mh_, xsLO_, xsNLO_, ___}] :=
  {Sqrt[mh], foo[1000 xsNLO]} //N

log["ggh-smnlo-lhc.dat", _][{sqrts_, mh_, xsLO_, xsNLO_, ___}] :=
  {Sqrt[mh], Log[sqrts/1000], foo[1000 xsNLO]} //N

log["hww-lo.dat", _][{mh_, xshww_, xshzz_, ___}] :=
  {mh, foo[xshww/1000]} //N

log["hzz-lo.dat", _][{mh_, xshww_, xshzz_, ___}] :=
  {mh, foo[xshzz/1000]} //N

log["hww-nlo.dat", _][{mh_, _, _, xshww_, xshzz_, ___}] :=
  {mh, foo[xshww/1000]} //N

log["hzz-nlo.dat", _][{mh_, _, _, xshww_, xshzz_, ___}] :=
  {mh, foo[xshzz/1000]} //N

log[file_, _][{mh_, xs_?NumberQ, ___}] := ({Sqrt[mh], foo[1000 xs]} //N) /;
  MemberQ[
    { "tth-lhc14.dat", "vbf-lhc14.dat", "wh-lhc14.dat", "zh-lhc14.dat",
      "tth-lhc8.dat",  "vbf-lhc8.dat",  "wh-lhc8.dat",  "zh-lhc8.dat",
      "tth-lhc7.dat",  "vbf-lhc7.dat",  "wh-lhc7.dat",  "zh-lhc7.dat"},
    file ]

log[__][{mh_, xs_?NumberQ, ___}] := {Sqrt[mh], foo[xs]} //N


warn[m|sqrtm] = MHiggs

warn[mst1] = MStop1

warn[other_] = other


poly[vars__, n_Integer, s_:"a"] := 
  MapIndexed[#1 ToSymbol[s, #2]&,
    Expand[Sum[Plus[vars]^i, {i, 0, n}]]]

(*fitfunc[tHm, _LHC] = poly[m, tb, 9]*)
fitfunc[tHm, _LHC] = If@@ {
  sqrtm < Sqrt[600.],    poly[sqrtm, tb, 9],
                         poly[sqrtm, tb, 8] }

fitfunc[tHm2|tHm2lo|tHm2hi, _LHC] = poly[sqrtm, tb, 2]

fitfunc[_StSth, _LHC] = poly[sqrtm, mst1, 2]/poly[sqrtm, mst1, 2, "b"]

fitfunc[_gghSMNLO, Tev] = Which@@ {
  sqrtm < Sqrt[2 173.3], poly[sqrtm, 5],
  sqrtm < Sqrt[800.],    poly[sqrtm, 7, "b"],
  True,                  poly[sqrtm, 4, "c"] }

(*fitfunc[_gghSMNLO, _LHC] = poly[sqrtm, logsqrts, 9]*)
(*fitfunc[_gghSMNLO, _LHC] = If@@ { sqrtm < Sqrt[2 173.3],
  poly[sqrtm, logsqrts, 7],
  poly[sqrtm, logsqrts, 7, "b"] }*)
fitfunc[_gghSMNLO, _LHC] = Which@@ {
  sqrtm < Sqrt[2 173.3], poly[sqrtm, logsqrts, 7],
  sqrtm < Sqrt[800.],    poly[sqrtm, logsqrts, 7, "b"],
  True,                  poly[sqrtm, logsqrts, 4, "c"] }

fitfunc[_gghSM, Tev] = If@@ { sqrtm < Sqrt[2 173.3],
  poly[sqrtm, 5],
  poly[sqrtm, 7, "b"] }

(*fitfunc[_gghSM, _LHC] = If@@ { sqrtm < Sqrt[2 173.3],
  poly[sqrtm, logsqrts, 3],
  poly[sqrtm, logsqrts, 5, "b"] }*)
fitfunc[_gghSM, _LHC] = Which@@ {
  sqrtm < Sqrt[2 173.3], poly[sqrtm, logsqrts, 3],
  sqrtm < Sqrt[800.],    poly[sqrtm, logsqrts, 5, "b"],
  True,                  poly[sqrtm, logsqrts, 2, "c"] }

fitfunc[_qqhSM, LHC[8]] = If@@ { sqrtm < Sqrt[290.],
  poly[sqrtm, 2],
  poly[sqrtm, 10, "b"] }

fitfunc[_qqhSM, LHC[14]] = poly[sqrtm, 4]

fitfunc[_bbhSM, _LHC] = poly[sqrtm, sqrts, 2]

fitfunc[_bbhSM, _] = poly[sqrtm, 3]

fitfunc[_tthSM, LHC[14]] = poly[sqrtm, 4]

fitfunc[_WhSM, LHC[8]] = poly[sqrtm, 4]

fitfunc[_ZhSM, LHC[8]] = poly[sqrtm, 4]

fitfunc[_WhSM, LHC[14]] = poly[sqrtm, 4]

fitfunc[GammaSM[_H0VV], _] = If@@ {
  m < 161., poly[m, 3] (*poly[m, 4]/poly[m, 4, "b"]*),
            poly[m, 4, "c"]/poly[m, 4, "d"] }

(*fitfunc[_hZZSM, _] = poly[m, 4]/poly[m, 4, "b"]*)

_fitfunc = poly[sqrtm, 2]

indeps = {m, sqrtm, sqrts, logsqrts, mst1, tb}


complete[d_] := d /; Length[ d[[1]] ] =!= 3

complete[d_] := Join[d, {Sqrt[2000], #, 10^-10}&/@ Union[#2&@@@ d]]


overlap = 10;
ran[-Infinity, b_] = #[[1]]^2 <= b^2 + overlap &;
ran[a_, Infinity] = #[[1]]^2 >= a^2 - overlap &
ran[a_, b_] = #[[1]]^2 >= a^2 - overlap && #[[1]]^2 <= b^2 + overlap &;

dofit[d_, If[x_ < y_, a_, b_] + r_.] := IndexIf@@ {
  x < y, dofit[Select[d, ran[-Infinity, y]], a + r],
         dofit[Select[d, ran[y, +Infinity]], b + r]}

dofit[d_, Which[x_ < y_, a_, x_ < z_, b_, True, c_] + r_.] := IndexIf@@ {
  x < y, dofit[Select[d, ran[-Infinity, y]], a + r],
  x < z, dofit[Select[d, ran[y, z]], b + r],
         dofit[Select[d, ran[z, Infinity]], c + r]}

dofit[d_, fitfun_] := fitfun /. FindFit[d, fitfun,
  Complement[
    Cases[fitfun, x_Symbol /; Context[x] =!= "System`", {-1}],
    fitvars ],  
  fitvars]


exp[i_IndexIf] := MapIf[exp, i]

fit[tag_][id_ -> file_] :=
Block[ {data, func, fitvars, thecat},
DATA=
  data = get[tag, file];
FUNC=
  func = fitfunc[id, tag];
FITVARS=
  fitvars = Select[indeps, !FreeQ[func, #]&];
  func = dofit[data, func];

  _thecat = 0;
  cat/@ data;
  _thecat =.;
  Cases[ DownValues[thecat], _[_[_[r___]], d_] :> 
    plot[debug, id, func, d, r] ];

  id -> {exp[func] /. 1. -> 1,
    Transpose[{fitvars, Min/@ #, Max/@ #}]& @
      Drop[Transpose[data], -1] }
]

fit[tag_][other_] = other


cat[{mh_, r___, xs_}] := thecat[r] = {thecat[r], {mh, xs}}


plot[True, id_, f_, d_, r___] :=
Block[ {dd = Level[d, {-2}], rul = Thread[Rest[fitvars] -> {r}]},
  If[ !FreeQ[f, sqrtm], dd = {#1^2, #2}&@@@ dd ];
  Show[
    ListPlot[Sort[dd],
      PlotStyle -> Red,
      DisplayFunction :> Identity],
    Plot[Evaluate[f /. rul /. sqrtm -> Sqrt[m]],
      {m, 30, plotmax},
      PlotStyle -> Blue,
      PlotRange -> All (* ({Min[#], Max[#]}&[Last/@ dd]) *),
      DisplayFunction :> Identity],
    DisplayFunction :> $DisplayFunction,
    PlotLabel -> ToString[id] <>
      ({" ", ToString[#1], "=", ToString[#2]}&@@@ rul),
    PlotRange -> {All, All} ]
] /; debug === True


normalize[expr_, vars_] :=
Block[ {num, den, d0},
  num = Numerator[expr];
  den = Denominator[expr];
  d0 = den /. Alternatives@@ vars -> 0;
  HornerForm[Distribute[num/d0], vars]/
    HornerForm[Distribute[den/d0], vars]
]


MkDir["f"];
MkDir["m"];

$DebugCmd = {"#ifdef DETAILED_DEBUG\n", "#endif\n", "DPROD ", " ENDL"};
SetOptions[PrepareExpr, DebugLines -> True];

write[file_, fits_] :=
Block[ {out},
  Put[fits, ToFileName["m", file <> ".m"]];
  out = OpenFortran[ToFileName["f", file <> ".f"]];
  writefit[out]/@ fits;
  Close[out];
  fits
]

writefit[out_][id_[h_, i_] -> r_] := (
  WriteString[out, "\n\tif( " <> ToFortran[h == i] <> " ) then"];
  writefit[out][id[h] -> r];
  WriteString[out, "\tendif\n"]
)

writefit[out_][id_ -> r_] := (
  WriteString[out, Apply[rangecheck[id], r[[2]], 1] <> "\n"];
  WriteExpr[out, id -> r[[1]]]
)

writefit[out_][other_String] := WriteString[out, other]


rangecheck[lhs_][var_Symbol, lo_, hi_] :=
  {"\n\tif( ", ToFortran[var], " .lt. ", ToFortran[lo],
(*
   " .or.\n     &      ", ToFortran[var], " .gt. ", ToFortran[hi],
*)
   " )\n     &    Warning('Extrapolating ",
     tof[lhs], " in ", tof[warn[var]], "')"}

rangecheck[_][__] = ""


tof[s_] := StringReplace[ToFortran[s /. h -> $h$],
  "$" ~~ x_ ~~ "$" :> "'//Digit(" <> x <> ")//'"]

(*
tof[s_Symbol] := ToFortran[s]

tof[s_[h_, ___]] := {ToFortran[s], "('//Digit(", ToFortran[h], ")//')"}
*)

nfitsTev := write["NHiggsProdFits-Tev", fit[Tev]/@ {
  bbhSM[h] -> "bbh-tev.dat",
  btagbhSM[h] -> "bbh-tev-tag.dat",
  gghSM[h] -> "ggh-tev.dat",
  gghSMNLO[h] -> "ggh-smnlo-tev.dat",
  tthSM[h] -> "tth-tev.dat",
  qqhSM[h] -> "vbf-tev.dat",
  WhSM[h] -> "wh-tev.dat",
  ZhSM[h] -> "zh-tev.dat" } /. {sqrts -> 2, logsqrts -> Log[2.]}]

nfitsLHC := write["NHiggsProdFits-LHC", fit[LHC[all]]/@ {
  bbhSM[h] -> "bbh-lhc.dat",
  gghSM[h] -> "ggh-lhc.dat",
  gghSMNLO[h] -> "ggh-smnlo-lhc.dat" }]

nfitsLHC7 := write["NHiggsProdFits-LHC7", fit[LHC[7]]/@ {
  tthSM[h] -> "tth-lhc7.dat",
  qqhSM[h] -> "vbf-lhc7.dat",
  WhSM[h] -> "wh-lhc7.dat",
  ZhSM[h] -> "zh-lhc7.dat" }]

nfitsLHC8 := write["NHiggsProdFits-LHC8", fit[LHC[8]]/@ {
  tthSM[h] -> "tth-lhc8.dat",
  qqhSM[h] -> "vbf-lhc8.dat",
  WhSM[h] -> "wh-lhc8.dat",
  ZhSM[h] -> "zh-lhc8.dat" }]

nfitsLHC14 := write["NHiggsProdFits-LHC14", fit[LHC[14]]/@ {
  btagbhSM[h] -> "bbh-lhc14-tag.dat",
  (*gghSM[h] -> "ggh-lhc.dat",*)
  tthSM[h] -> "tth-lhc14.dat",
  qqhSM[h] -> "vbf-lhc14.dat",
  WhSM[h] -> "wh-lhc14.dat",
  ZhSM[h] -> "zh-lhc14.dat",
  StSth[h] -> "ststh-lhc14.dat" }]

cfitsLHC14 := write["CHiggsProdFits-LHC14", fit[LHC[14]]/@ {
  tHm -> "tHm-lhc14.dat" } /. tb -> TBeff]

cfitsLHC8 := write["CHiggsProdFits-LHC8", fit[LHC[8]]/@ {
  tHm2 -> "tHm2-lhc8.dat",
  tHm2lo -> "tHm2-lhc8lo.dat",
  tHm2hi -> "tHm2-lhc8hi.dat" } /. tb -> TBeff]

decayfitsLO := write["NHiggsDecayFitsLO", fit[Decay]/@ {
  GammaSM[H0VV[h,3]] -> "hzz-lo.dat",
  GammaSM[H0VV[h,4]] -> "hww-lo.dat" }]

decayfitsNLO := write["NHiggsDecayFitsNLO", fit[Decay]/@ {
  GammaSM[H0VV[h,3]] -> "hzz-nlo.dat",
  GammaSM[H0VV[h,4]] -> "hww-nlo.dat" }]

t := write["test", fit[LHC[all]]/@ {
  gghSM[h] -> "ggh-lhc.dat" }]

tnlo := write["test", fit[LHC[all]]/@ {
  gghSMNLO[h] -> "ggh-smnlo-lhc.dat" }]

