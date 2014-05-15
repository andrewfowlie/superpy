(*

Syntax for the SLHA descriptor:

* blockname -> {members},
  defines a block,

* blockname :> {members}
  defines an overlapping block, i.e. does not take own storage,
  much like an equivalence statement in Fortran,

* members: either
  - a single symbol, "name",
  - an array of the form "name[{i1,[min,]max},...]",
  - an alias of the form "alias == name" 
    (takes up no separate storage),
  - an alias of the form "alias -> name"
    (has its own storage, used to construct overlapping blocks)

  aliases from a different block are written as "blockname**name".

*)



BeginPackage["GenSLHA`"]

Struct::usage = "Struct[ind][mem] defines a structure with members mem,
defined for index ranges ind."

SLHADefs::usage = "SLHADefs[desc] generates the definitions for
SLHADefs.h from the SLHA descriptions in desc."

SLHANames::usage = "SLHANames[desc] generates the data definitions for the
SLHA blocks in desc."

(* internal symbols *)

{i, Q, Severity, NLines, Code, Text, TextFlat,
  bracket, Slhadata, SlhaData}


Begin["`Private`"]

full = blockname <> "_" <> str[#1] &;

glob = "#define " <> #1 <> blockname <> " " <> str[#2] <> "\n" &

def = "#define " <> full[#2] <> " " <>
  str[data[#2, #3 + (offset += #1) - #1 + 1]] <> "\n" &

equ = "#define   " <> full[#1] <> " " <> #2 <> "\n" &

str = StringReplace[ToString[#, CForm], {" " -> "", "bracket" -> ""}] &

flat = ToExpression[ToString[#] <> "Flat"] &


data[n_] := SlhaData[n]

data[Q | _Q | Severity | NLines | _Code | _Text | _TextFlat, n_] := SlhaData[n]

data[_, n_] := SlhaData[n] /; blockname == "QExtPar"

data[name_, n_] := Slhadata[n]


block[name_ :> members_] :=
Block[ {offset = offset},
  block[name, members]
]

block[name_[_] -> members_] := block[name, members]

block[name_, members_] :=
Block[ {blockname = str[name], old, res},
  old = offset;
  res = member/@ members;
  glob["Offset", old] <>
  glob["Length", offset - old] <>
  "#define Block" <> blockname <> "(i) " <> str[data[old + i]] <> "\n" <>
  res <>
  "\n"
]


add[a_ == b_] := add[a] == b

add[a_ -> b_] := add[a] -> b

add[s_Symbol] := s[rr]

add[s_[j__]] := s[j, rr]


member[Struct[r__List][members__]] :=
Block[ {len, old, stride, rr, res},
  len = range[1, {}, 0][r][[1]];
  old = offset;
  member/@ {members};
  stride = offset - old;
  rr = Sequence@@ Insert[{r}, stride, {1, 1}];
  offset = old;
  res = member/@ add/@ {members};
  offset = old + stride len;
  res
]

member[s_Symbol] := def[1, s, 0]

member[s_Symbol[r_List]] := def@@ range[1, s[], 0][r]

member[s_Symbol[r__List]] :=
Block[ {d = def[0, #2, #3]},
  If[ Length[#4] > 1, d = {d, def[#1,
    Operate[flat, Prepend[Drop[#2, Length[#4]], i]], #3 - #5 + i - 1]} ];
  d
]&@@ range[1, s[], 0][r]

member[a_ == i_Integer] := equ[a, str[i]]

member[a_ == b_**s_] := equ[a, str[b] <> "_" <> str[s]]

member[a_ == s_] := equ[a, full[s]]

member[s_Symbol[r_List] -> b_] :=
  (offset += #1; member[#2 == b])&@@ range[1, s[], 0][r]

member[s_Symbol[r__List] -> b_] := (offset += #1;
  {member[#2 == b],
   member[flat[Head[#2]][i] ==
     (b /. x_[Sequence@@ #2, j__] :> flat[x][i, j])]})&@@ range[1, s[], 0][r]

member[s_Symbol -> b_] := (++offset; member[s == b])

member[other_] := (Message[SLHADefs::syntax, other, blockname]; Abort[])

SLHADefs::syntax = "Syntax error in `` in block ``."


range[len_, eq__][{stride_Integer, var__}, o___] :=
  {len, #2, #3, eq}&@@ range[stride, eq][{var}, o]

range[len_, lhs_, rhs_][{var_Symbol, min_:1, max_}, o___] :=
  range[len (max - min + 1), Append[lhs, var],
    rhs + bracket[len, var] - len min][o]

range[len_, eq__][] := {len, eq, eq}


bracket[1, s_] = s

bracket[n_, s_] = n bracket[s]


TimeStamp[] :=
  ToString[#3] <> " " <>
  {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
   "Sep", "Oct", "Nov", "Dec"}[[#2]] <> " " <>
  ToString[#1] <> " " <>
  ToString[#4] <> ":" <>
  StringTake["0" <> ToString[#5], -2]&@@ Date[]


SLHADefs::info = "Descriptor contains `` SLHA Blocks."


SLHADefs[desc_] :=
Block[ {offset = 0},
  Message[SLHADefs::info, Count[desc, _Rule]];
  "\
#if 0\n\
\tSLHADefs.h\n\
\t\tdeclarations for SLHALib data\n\
\t\tgenerated " <> TimeStamp[] <> "\n\
#endif\n\n\
#ifndef SLHADEFS_H\n\
#define SLHADEFS_H\n\n\
#define invalid (-999)\n\n\
" <> block/@ desc <> "\
#define nslhadata " <> str[offset] <> "\n\n\
#endif\n"
]


SLHANames[desc_, dname_, prefix_] :=
Block[ {blocks, labels, maxlen},
  {blocks, labels} = Transpose @
    Cases[desc, (b_[l_] -> _) :> {ToString[b], ToString[l]}];
  maxlen = Max[StringLength/@ blocks] + StringLength[prefix];
  MapIndexed[ {"#define ", #1, " ", ToString@@ #2, "\n"}&,
    labels ] <> "\n\
\tinteger nblocks\n\
\tparameter (nblocks = " <> str[Length[blocks]] <> ")\n\
\tcharacter*" <> str[maxlen] <> " " <> dname <> "(nblocks)\n" <>
  Apply[ {"\tdata ", dname, "(", #2, ") /\"",
    prefix, ToUpperCase[#1], "\"/\n"}&,
    Drop[Transpose[{blocks, labels}], -1], 1 ]
]


End[]

EndPackage[]

