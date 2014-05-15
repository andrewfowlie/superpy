
(* FeynCalc.m *)
(* this file should be called FeynCalc.m *)

(* ******************************************************************** *)
(* ********************** F E Y N C A L C 2.2beta6  ******************* *)
(*                           October 1994                             *)
(* ******************************************************************** *)

                 (* Copyleft     Rolf Mertig  1988 - 1994 *)


(*       NOTE: This version deviates somewhat from the manual!!!  
           A complete new version is being worked on, 
             together with a manual. 
                SUGGESTIONS and WISHES are most WELCOME !!!!!!!!
*)


(* ******************************************************************** *)
(* Mathematica 2.0 or higher required.                                  *) 
(* Version 2.2 works best *)
(* ******************************************************************** *)
(* In case you do not get to work what the manual claims is              *)
(* possible, feel free to write to the above adress or send e-mail to   *)
(* rolfm@rulgm4.leidenuniv.nl                                           *)
(* or to                                                                *)
(* rolfm@wri.com                                                        *)
(* ******************************************************************** *)
Print[" "];
WriteString["stdout", "FeynCalc 2.2beta6 (October1994) |"];
(* ******************************************************************** *)


(* ------------------------------------------------------------------*)
(* $BreitMaison is introduced in the context `Global *)
If[Global`$BreitMaison=!=True, Global`$BreitMaison = False];

$BreitMaison::usage=
"The setting of $BreitMaison determines whether the Breitenlohner-
Maison scheme is applied. If $BreitMaison=True, the so-called
naive gamma5 prescription is used, i.e. gamma5 anticommutes in
all dimensions.  The default is False.
$BreitMaison can only be set to True BEFORE loading FeynCalc.
The setting is then irreversible during the session. ";

(* ------------------------------------------------------------------*)

If[Attributes[Feyn`Symbols`$MU]==={Protected}, Unprotect[Feyn`Symbols`$MU]];
Clear["Feyn`Symbols`*"];
BeginPackage["Feyn`Symbols`"]

A0::usage=
"A0[m^2] is the Passarino - Veltman one point integral.";
 
A0ToB0::usage
"A0ToB0 is an option for A0. If set to True, A0[m^2] is expressed
by (1+ B0[0,m^2,m^2])*m^2;";
 
B0::usage=
"B0[pp,m1^2,m2^2] is the Passarino - Veltman two point integral.
All arguments are scalars and have dimension mass^2.";

B0Real::usage=
"B0Real is an option of B0 (default False). If set to True,
B0 is assumed to be real and
the relation B0[a,0,a] = 2 + B0[0,a,a]  is applied.";
 
B0Unique::usage=
"B0Unique is an option for B0. If set to True, B0[0,0,m2] is replcaed
by (B0[0,m2,m2]+1) and B0[m2,0,m2] simplifies to (B0[0,m2,m2]+2).";

B00::usage=
"B00[pp,m1^2,m2^2] is the Passarino - Veltman B00-function, i.e. the
coefficient function of g(mu nu). All arguments are scalars and have
dimension mass^2.";
 
B1::usage=
"B1[pp,m1^2,m2^2] is the Passarino - Veltman B1-function.
All arguments are scalars and have dimension mass^2.";
 
B11::usage=
"B11[pp,m1^2,m2^2] is the Passarino - Veltman B11-function, i.e.
the coefficient function of p(mu) p(nu).";

Born::usage=
"Born is an option of OneLoopSum. It may be set to the result of the
 (polarized) Born amplitude, which is multiplied with the result of
  OneLoopSum (though not at the end, but intermediately).";

BReduce::usage=
"BReduce is an option for B0, B00, B1, B11 determining whether
reductions to A0 and B0 will be done. ";
 
C0::usage=
"C0[p10, p12, p20, m1^2, m2^2, m3^2] is the scalar
Passarino - Veltman C0-function.  The convention for the arguments
is that if the denominator of the integrand has the form
([q^2-m1^2] [(q+p1)^2-m2^2] [(q+p2)^2-m3^2]),
the first three arguments of C0 are the scalar products
 p10 = p1^2, p12 = (p1-p2).(p1-p2), p20 = ,p2^2.";

CA::usage=
"CA is one of the Casimir operators of SU(N); CA = N"; 

CF::usage=
"CF is one of the Casimir operators of SU(N); CF = (N^2-1)/(2 N)"; 

CancelQ2::usage=
"CancelQ2 is an option for OneLoop. If set to True, cancellation
of all q^2's ( and (q^2)^(2n) ) with the first propagator via
q^2 ->( (q^2-m^2) + m^2 ) is performed.";

CancelQP::usage=
"CancelQP is an option for OneLoop. If set to True, cancelation of
all q.p's ( and (q.p)^(2n) ) is performed.";

ChangeDimension::usage=
"ChangeDimension[exp, dim] changes all LorentzIndex and Momenta in 
exp to dimension dim.";

ChargeConjugationMatrix::usage=
"ChargeConjugationMatrix denotes the charge conjugation matrix C.";

ChargeConjugationMatrixInverse::usage=
"ChargeConjugationMatrixInverse is the inverse of ChargeConjugationMatrix.";

ChiralityProjector::usage=
"ChiralityProjector[+1] is an input function for 1/2( 1 + gamma5 ).
ChiralityProjector[-1] is an input function for 1/2( 1 - gamma5 ).";

Chisholm::usage=
"Chisholm[x] substitutes products of three Dirac matrices or
slashes by the Chisholm identity.";

ChisholmSpinor::usage=
"ChisholmSpinor[x] uses for a DiraGamma between spinors
the Chisholm identity. As an optional second argument 1 or 2 may be given, indicating that ChisholmSpinor should only act upong the first resp. second part of a product of spinor chains.";

Collecting::usage=
"Collecting is an option of Contract2 and SquareAmplitude. 
If set to False, the result
is not collected w.r.t. to the K's.";

Collect2::usage=
"Collect2[expr, x] collects together terms which are not
free of any occurrence of x.
Collect2[expr, {x1, x2, ...}] collects together terms
 which are not free of any occurrence of x1, x2, ... .
The coefficients as well as the constant term are put over a 
common denominator and factored (depending on the option Factoring).";

Combine::usage=
"Combine[expr] puts terms in a sum over a common denominator, and cancels
factors in the result. Combine is similar to Together, but
accepts the option ProductExpand (default False) and works better than
Together on polynomials involving rationals with sums in the denominator.";

CombineGraphs::usage=
"CombineGraphs is an option for OneLoopSum.
Graphs indicated by its setting are summed before calculated with
OneLoop. Possible settings are:
CombineGraphs -> All;
CombineGraphs -> { {i, j}, {k, l, m, ...}, ... } which sums
the amplitudes number i and j into a single one as well as
number k, l, m, ... . The counting extends over all amplitudes
provided to OneLoopSum.
The default value is {}, if the setting is False, no summation
 will be done at the end of OneLoopSum.";

ComplexConjugate::usage=
"ComplexConjugate[expr] complex conjugates expr.
It operates on  Fermion-lines (products
of Spinor[..] .DiracMatrix[..] . Spinor[..]) and changes all
occuring LorentzIndex[mu] into LorentzIndex[ComplexIndex[mu]].
For taking the spin sum (i.e. constructing the traces) use
FermionSpinSum.";

ComplexIndex::usage=
"ComplexIndex is the head of a complex conjugate index.";

Contract::usage=
"Contract[expr] contracts pairs of Lorentz indices of metric tensors,
four-vectors and (depending on the optino EpsContract) of
Levi-Civita tensors in expr. For the contraction of Dirac matrices
with each other use DiracSimplify. \n \n
Contract[exp1, exp2] contracts (exp1*exp2), where exp1 and exp2 may be
larger products of sums of  metric tensors and 4-vectors.
Contract[exp1, exp2] should be used for polarization sums, where
exp2 should be the product (or expanded sum) of the polarization 
sums for the vector bosons.";

Contract2::usage=
"Contract2[expr] (still experimental).";

Contract3::usage=
"Contract3[expr] works like Contract, but more efficient if 
there are no Levi-Civita tensors around.";

D0::usage=
"D0[ p10, p12, p23, p30, p20, p13,  m1^2, m2^2, m3^2, m4^2 ] is the
 Passarino-Veltman D0-function. The convention for the arguments is
that if the denominator of the integrand has the form
( [q^2-m1^2] [(q+p1)^2-m2^2] [(q+p2)^2-m3^2] [(q+p3)^2-m4^2] ),
 the first six arguments of D0 are the scalar products
p10 = p1^2, p12 = (p1-p2)^2, p23 = (p2-p3)^2, p30 = p3^2,
p20 = p2^2, p13 = (p1-p3)^2.";

D0Convention::usage=
"D0Convention is an option for Write2. If set to 1, the convention for
the arguments of D0 is changed when writing a Fortran file with Write2:
The fifth and sixth argument of D0 are interchanged and the square root is
taken of the last four arguments.";

DB0::usage=
"DB0[p2,m1^2,m2^2] is the derivative of the two-point function
B0[p2,m1^2,m2^2] with respect to p2.";

DB1::usage=
"DB1[p2,m1^2,m2^2] is the derivative of B1[p2,m1^2,m2^2] with respect to p2.";

DenominatorOrder::usage=
"DenominatorOrder is an option for OneLoop, if set to True,
the FeynAmpDenominator's will be ordered canonically."

Dimension::usage=
"Dimension  is an option for DiracMatrix, DiracSlash, FourVector,
MetricTensor, OneLoop and ScalarProduct.";

DiracGamma::usage=
"DiracGamma[x, optdim]  is the head of all Dirac matrices and
Feynman slashes in the internal representation. Use DiracMatrix
and DiracSlash for input.";

DiracGammaT::usage=
"DiracGammaT[x] denotes the transpose of DiracGamma.";

DiracMatrix::usage=
"DiracMatrix[mu] is an input function for a Dirac matrix gamma(mu).
DiracMatrix[mu, nu, ...] or DiracMatrix[mu . nu . ...]
or DiracMatrix[mu] . DiracMatrix[nu] . ... is the input
for a product of (4-dimensional) Dirac matrices gamma(mu) gamma(nu) ... .
The dimension of the Dirac matrix may specified by the option
Dimension.
DiracMatrix[5], DiracMatrix[6] and DiracMatrix[7] are denote the
algebraic objects gamma5, gamma6 ( = 1/2 (1 + gamma_5) ) and
gamma7 ( = 1/2 (1 - gamma_5) ).";

DiracOrder::usage=
"DiracOrder[expr] orders the Dirac matrices in expr alphabetically.
DiracOrder[expr, orderlist] orders the Dirac matrices in expr according
to orderlist.";

DiracSimplify::usage=
"DiracSimplify[expr]   simplifies products of Dirac matrices
in expr. Double Lorentz indices and four vectors are contracted.
The Dirac equation is applied.
All DiracMatrix[5], DiracMatrix[6] and DiracMatrix[7] are moved to
the right. The order of the Dirac matrices is not changed.";

DiracSlash::usage=
"DiracSlash[p] is an input function for a Feynman slash.
A product of slashes may be entered by:
DiracSlash[p, q, ...] or
DiracSlash[p . q . s. ...] or DiracSlash[p] . DiracSlash[q] . ....";
 
DiracSpinor::usage=
"DiracSpinor[a,{1}] evaluates to Spinor[a, 1]";

DiracTrace::usage=
"DiracTrace[expr] is the head of Dirac Traces.
Whether the trace is  evaluated depends on the option
DiracTraceEvaluate.
The argument expr may be a product of Dirac matrices or slashes
separated by the Mathematica Dot \".\".";
 
DiracTraceEvaluate::usage=
"DiracTraceEvaluate   is an option for DiracTrace.
If set to False, DiracTrace remains unevaluated.";
 
Eps::usage=
"Eps[a, b, c, d] is the head of the totally antisymmetric epsilon
(Levi-Civita) tensor. The \"a,b, ...\" MUST have head LorentzIndex or
 FourVector.  Note that Eps has also an option Dimension (default 4).
For userfriendly input of Eps with Lorentz indices use LeviCivita.";

EpsAway::usage=
"EpsAway is an option for SquareAmplitude. If set to True all 
Levi-Civita tensors are replaced by 0 after squaring.";
 
EpsChisholm::usage=
"EpsChisholm[expr] substitutes for a gamma matrix contracted with
a Levi Civita tensor (Eps) the Chisholm identity.";
 
EpsContract::usage=
"EpsContract is an option  Contract specifying whether Levi-Civita
tensors Eps[...] will be contracted, i.e., products
of two  Eps are replaced via the determinant formula.";
 
EpsEvaluate::usage=
"EpsEvaluate[ expr ] applies total antisymmetry and
linearity (w.r.t. Momentum's) to all Levi-Civita tensors (Eps's)
in expr.";
 
EvaluateDiracTrace::usage=
"EvaluateDiracTrace[ expr ] evaluates all DiracTrace[...] in expr.";

Expand2::usage=
"Expand2[exp, x] expands all sums containing x.
Expand2[exp, {x1, x2, ...}]  expands all sums containing x1, x2, ....
Expand2[exp] is equivalent to Expand[x].";

Expanding::usage=
"Expanding is an option for Contract and DiracSimplify.
As option for Contract it  specifies whether expansion
will be done.  If set to False in DiracSimplify only a limitid
set of simplifications (multiplicative linearity etc.) is 
performed in DiracSimplify.";

 
ExpandScalarProduct::usage=
"ExpandScalarProduct[expr]  expands scalar products of sums of momenta
in expr.";

ExtraFactor::usage=
"ExtraFactor is an option for SquareAmplitude.";
 
ExtraVariables::usage=
"ExtraVariables is an option for OneLoopSum. It may be set to a
list of variables, which will be treated equivalently to the
scalar B0, C0, D0 w.r.t. to ordering.";

Factor1::usage=
"Factor1[poly] factorizes common terms  in the summands of poly.
It uses basically PolynomialGCD.";


Factor2::usage=
"Factor2[poly] factors a polynomial in a standard way.
Factor2 works better than Factor on polynomials involving rationals with
sums in the denominator. The maximum time in seconds Factor2 tries to factor,
is determined by the option FactorTime (default: 3600).";

Factor5::usage="experimental";

Factoring::usage=
"Factoring is an option for Collect2, Contract, Tr and  TraceEvaluate.
If set to True, the result will be factored, using Factor2.";

FactorFull::usage=
"FactorFull is an option of Factor2 (default False).
If set to False, products like
(a-b) (a+b) will be replaced by (a^2-b^2).";

FactorTime::usage=
"FactorTime is an option for Factor2. It denotes the maximum time
(in seconds) during which Factor2 tries to factor.";

FermionSpinSum::usage=
"FermionSpinSum[x] constructs the Traces out of squared ampliudes."; 

FeynAmp::usage=
"FeynAmp[name, q, amp]  is the head of the Feynman amplitude
given by the program FeynArts. \"name\" does the bookkeeping,
\"amp\"  is the analytical expression for the amplitude
and \"q \"is the integration variable.
In order to calculate the amplitute replace FeynAmp by OneLoop.";

FeynAmpList::usage=
"FeynAmpList[information][FeynAmp[...], FeynAmp[...], ...]
is a head from FeynArts denoting a collection of amplitudes.";

FeynAmpDenominator::usage=
"FeynAmpDenominator[ PropagatorDenominator[ ... ],
PropagatorDenominator[ ... ], ... ] is the head of
the denominators of the propagators, i.e. FeynAmpDenominator[ x ]
is the representation of 1/x.";

FeynCalcForm::usage=
"FeynCalcForm[expr] changes the printed output to a an easy to read form.
The default setting of $PrePrint is $PrePrint = FeynCalcForm, which
forces to display everything after applying FeynCalcForm. In order
to change to the normal (internal) Mathematica OutputForm do:
($PrePrint=.) .";

FinalFunction::usage=
"FinalFunction is an option for OneLoopSum. The function (head) indicated
acts upon the coeffients of the scalar integrals.";

FinalSubstitutions::usage=
"FinalSubstitutions  is an option for OneLoop and OneLoopSum and
Write2.
All substitutions indicated hereby are done at the
end of the calculation.";

FourVector::usage=
"FourVector[p, mu] is the input for a 4-dimensional Lorentz vector
p(mu).  It is transformed to the internal representation
Pair[Momentum[p], LorentzIndex[mu]]."

FreeQ2::usage=
"FreeQ2[expr, {form1, form2, ...}] yields True if expr does not contain any
occurence of form1, form2, ... and False otherwise.
FreeQ2[expr, form] is the same as FreeQ[expr, form].";

FUNCTION::usage=
"FUNCTION[exp, string]  is a head of an expression to be declared a
function (of type string), if used in Write2.";
 
FRH::usage="FRH[x] is equivalent to FixedPoint[ReleaseHold, x].";

GellMannMatrix::usage= "GellMannMatrix is replaced by SUNT.";

GellMannTrace::usage= "GellMannTrace is replaced by SUNTrace.";

GetOneLoopResult::usage=
"GetOneLoopResult[fname, {1,2, ...}] returns the sum of a
previously save results by OneLoop. If the files have names
as V1V1V1V1N1.m, V1V1V1V1N2.m, ..., the argument fname should
be V1V1V1V1, and the list denotes which N1, N2 extensions are
actually loaded.";

GraphName::usage=
"GraphName is the head of the first argument of OneLoop.
The arguments of GraphName should be , e.g., :
GraphName[ eezh, ... , N1 ],  where the last argument indicates
the running number of the amplitude. Upon setting the option
WriteOut of OneLoop to a different setting than False, e.g.
InputForm, a result file eezhN1.m is written out.";

InitialSubstitutions::usage=
"InitialSubstitutions is an option for OneLoop and OneLoopSum.
All substitutions indicated hereby will be performed at the
beginning of the calculation.";

IntermediateSubstitutions::usage=
"IntermediateSubstitutions is an options for OneLoop.
All substitutions indicated hereby will be performed
before tensorintegral decomposition.";

Isolate::usage=
"Isolate[expr] substitutes the abbreviation K[i] in HoldForm for
expr, if Length[expr]>0.
Isolate[expr, varlist] substitutes K[i] for all subsums in expr
which are free of any occurence of a member of the list varlist.
Instead of K any other head of the abbreviations may be specified
with the option IsolateHead.";

IsolateHead::usage=
"IsolateHead is an option for Isolate and Collect2.
Its default setting is K.";

IsolateSplit::usage=
"IsolateSplit is an option for Isolate. Its setting determines the
maximum number of characters of FortranForm[expr] which are
abbreviated by Isolate. If the  expression is larger than the
indicated maximum, it is split into smaller pieces and on each subsum
Isolate gets applied.
If no splitting is desired: IsolateSplit -> Infinity.";

K::usage=
"K[i] is the default setting of IsolateHead,
which is the head of abbreviations used by Isolate.
A K[i] returned by Isolate is given in HoldForm and can be
recovered by ReleasHold[K[i]].";

KeepOnly::usage=
"KeepOnly may be set to B0, C0, D0 keeping only the corresponding
coefficients. The default setting is False. If KeepOnly is set
to {} then the part of the amplitude which is not coefficient
of B0, C0, D0 is kept.";

LeptonSpinor::usage=
"LeptonSpinor[p, m] is the head of leptonic Dirac spinors.
Which of the spinors u, v,u_bar or v_bar
is understood, depends on the sign of the momentum argument and
the relative position of DiracSlash[p]: LeptonSpinor[sign p, mass]
is that spinor which yields  sign*mass*LeptonSpinor[p, mass] if
the Dirac equation is applied.";

LeviCivita::usage=
"LeviCivita[mu, nu, ro, si] is an input  function for the
totally antisymmetric Levi-Civita tensor. It transforms
to the internal representation Eps[ LorentzIndex[mu],  LorentzIndex[nu],
LorentzIndex[ro], LorentzIndex[si] ] (or with a second argument in 
LorentzIndex for the dimension, if the option Dimension of 
LeviCivita is changed) . For simplification of Levi-Civita
tensors use EpsEvaluate.";

LeviCivitaSign::usage=
"LeviCivitaSign is an option for DiracTrace and EpsChisholm. It determines
the sign in the result of a Dirac trace of four gamma matrices and gamma5.";

LorentzIndex::usage=
"LorentzIndex[mu] is the head of Lorentz indices.
The internal representation of a four-dimensional mu is LorentzIndex[mu].
For other than four dimensions enter LorentzIndex[mu, Dimension].
LorentzIndex[mu, 4] simplifies to LorentzIndex[mu].";

Mandelstam::usage=
"Mandelstam is an option for DiracTrace, OneLoop, OneLoopSum
and TrickMandelstam.  A typical setting is
Mandelstam -> {s, t, u, m1^2 + m2^2 + m3^2 + m4^2},
which stands for  s + t + u = m1^2 + m2^2 + m3^2 +  m4^2.
If other than scattering processes are calculated the setting
should be: Mandelstam -> {}.";
 
MetricTensor::usage=
"MetricTensor[mu, nu] is the input for a (4-dimensional) metric tensor.
For other than 4 dimensions change the options Dimension.";
 
Momentum::usage=
"Momentum is the head of a momentum in the internal representation.
A four-dimensional momentum p is represented by Momentum[p].
For other than four dimensions an extra argument must be given:
Momentum[q, dim].";

NF::usage=
"NF denotes the number of flavors.";

NumericalFactor::usage=
"NumericalFactor[expr] gives the numerical factor of expr.";

OneLoop::usage=
"OneLoop[name, q, amp ]
will calculate the one-loop Feynman amplitude amp.
The optional argument name indicates
the graph under consideration. q denotes the integration variable.";

OneLoopInfo::usage=
"OneLoopInfo[graphname] gives the last argument of FeynAmp[ ... ] generated
by FeynArts2.0.";

OneLoopResult::usage=
"OneLoopResult[ name ]  is the variable in the result file written out by
OneLoop to which the corresponding result is assigned; As name
the first argument of OneLoop is taken."

OneLoopSum::usage=
"OneLoopSum[ FeynAmp[ ... ], FeynAmp[ ... ] , ...]
will calculate a list of Feynman amplitudes by  replacing
FeynAmp  step by step by OneLoop.";

Pair::usage=
"Pair[a , b] is the head of a special pairing used in the internal
representation: a and b may have heads LorentzIndex or Momentum.
If both a and b have head
LorentzIndex, the metric tensor is understood.
If  a and b have head Momentum, a scalar product is meant
If one of  a and b has head LorentzIndex and the other
Momentum, a Lorentz vector (p_mu) is understood.";

PairCollect::usage=
"PairCollect is an option for DiracTrace specifying if
the result is collected with respect to Pair's.";

PartitHead::usage=
"PartitHead[expr, h] returns a list {ex1, h[ex2]} with ex1 freeof
expressions with head h, and h[ex2] having head h.";

PaVe::usage=
"PaVe[ i,j,... {p10,p12,...},{m1^2, mw^2, ...} ] denotes the invariant
(and scalar)
Passarino-Veltman integrals, i.e. the coefficient functions of
the tensor integral decomposition.  Joining plist and mlist gives the same
conventions as for A0, B0, C0, D0.  Automatic simlifications are
performed for the coefficient functions of two-point integrals and
for the scalar integrals.";

PaVeOrder::usage=
"PaVeOrder[expr] orders the arguments of all D0 in expr in a standard way.
PaVeOrder[expr, PaVeOrderList -> { {..., s, u, ...}, 
{... m1^2, m2^2, ...}, ...}] orders the arguments of all D0 in expr
according to the specified ordering lists.
The lists may contain only a subsequence of the D0-variables.";

PaVeOrderList::usage=
"PaVeOrderList is an option for PaVeOrder and PaVeReduce,
specifying in which order the arguments of D0 are to be permuted.";

PaVeReduce::usage=
"PaVeReduce[expr] reduces all Passarino-Veltman integrals
(i.e. all PaVe's) in expr down to scalar A0, B0, C0 and D0.";

Polarization::usage=
"Polarization[k] (= Polarization[k, I]) is the head of a 
4-dimensional polarization momentum with momentum k.
A slashed polarization vector (e1(k) slash) has to be entered
as  DiracSlash[Polarization[k]].
The internal representation for a polarization vector e1
corresponding to a boson with four momentum k is:
Momentum[ Polarization[ k, I ] ].
With this notation transversality of polarization vectors is
provided, i.e.  Pair[ Momentum[k],
Momentum[ Polarization[k, I] ] ] yields 0.
Polarization[k,-I] denotes the complex conjugate polarization.";

PolarizationSum::usage=
"PolarizationSum[ mu,nu, ... ] defines
(as abbreviations) different polarization sums.
PolarizationSum[mu, nu] = -g(mu nu);
PolarizationSum[mu, nu, k] = -g(mu nu) + k(mu) k(nu)/k^2;
PolarizationSum[mu, nu, k, n] = polarization sum for spin 1 fields;
(n = external momentum).
PolarizationSum[mu, nu, k, 0] is equivalent to -g(mu nu)";

PolarizationVector::usage=
"PolarizationVector[k, mu]  is an input function for a
polarization vector e(k)_mu. It is transformed to  the internal
representation Pair[Momentum[Polarization[p]], LorentzIndex[mu]].";

PostFortranFile::usage=
"PostFortranFile is an option for Write2 which may be set to a file
name (or a list of file names), which will be put at the end of the generated
Fortran file.";

Prefactor::usage=
"Prefactor is an option for OneLoop and OneLoopSum.
If set as option of OneLoop, the amplitude is multiplied by
Prefactor before calculation; if specified as option of OneLoopSum,
after calculation in the final result as a global factor.";

ProductExpand::usage=
"ProductExpand is an option for Combine. If set to False, no products are
expanded.";

PreFortranFile::usage=
"PreFortranFile is an option for Write2 which may be set to a file
name (or a list of file names),
which will be put at the beginning of the generated
Fortran file.";

PropagatorDenominator::usage=
"PropagatorDenominator[q, m] is a factor of the denominator of a
propagator.  If p is supposed to be D-dimensional enter:
PropagatorDenominator[Momentum[q, D], m].  What is meant is
1/(q^2-m^2).
PropagatorDenominator[p] evaluates to PropagatorDenominator[p, 0].";

PropagatorDenominatorExplicit::usage=
"PropagatorDenominatorExplicit[exp] inserts in exp all
PropagatorDenominator[p, m] to (ScalarProduct[p, p] - m^2)";

QuarkSpinor::usage=
"QuarkSpinor[p, m] is the head of hadronic Dirac spinors with
suppressed color index.
Which of the spinors u, v,u_bar or v_bar
is understood, depends on the sign of the momentum argument and
the relative position of DiracSlash[p]:
QuarkSpinor[sign p, mass]  is that spinor which yields
sign*mass*QuarkSpinor[sign p, mass] if the Dirac equation is applied.";

ReduceGamma::usage=
"ReduceGamma is an otpion for OneLoop. If set to True all
DiracMatrix6] and DiracMatrix[7] (i.e. all ChiralityProjector)
are reduced to Gamma5.";


ScalarProduct::usage=
"ScalarProduct[p, q] is the input for scalar product.
Expansion of sums of momenta in ScalarProduct is done with
ExpandScalarProduct. Scalar product may be set, e.g.,
ScalarProduct[a, b] = m^2; but a and b must not contain sums.
It is enouraged to always set ScalarProduct's BEFORE any
calculation. This improves the performance of FeynCalc.";

Schouten::usage="Schouten[expr] applies the Schouten identity on at most
42 terms in a sum. If Schouten should operate on larger expression you
can give a second argument, e.g.: Schouten[expr, 4711] which will work
on sums with less than 4711 terms.  \n
Schouten is also an option of Contract. It may be set to an integer
indicating the maximum number of terms onto which the function Schouten
will be applied.";

Scaling::usage=
"Scaling is an option of OneLoopSum.";

ReduceToScalars::usage=
"ReduceToScalars is an option for OneLoop  and OneLoopSum that
specifies whether the result will be reduced to scalar A0, B0, C0
and D0 scalar integrals.";

SelectGraphs::usage=
"SelectGraphs is an option for OneLoopSum indicating that only a
slected set of graphs of the list provided to OneLoopSum is to
be calculated.
Possible settings are: SelectGraphs -> { i, j,  ... }
or SelectGraphs -> { a, {b, c}, ...  }
which indicates the graphs to be taken from the list provided
to OneLoopSum. In the second setting the list {b, c} indicates that
all amplitudes from b to c should be taken.";

SetMandelstam::usage=
"SetMandelstam[s, t, u, p1, p2, p3, p4, m1, m2, m3, m4] defines the
Mandelstam variables  s=(p1+p2)^2, t=(p1+p3)^2, u=(p1+p4)^2 and sets
the pi on-shell: p1^2=m1^2, p2^2=m2^2, p3^2=m3^2, p4^2=m4^2.
Note that p1 + p1 + p2 + p3 + p4 = 0 is assumed. \n
SetMandelstam[x, {p1, p2, p3, p4, p5}, {m1, m2, m3, m4, m5}]
defines x[i, j] = (pi+pj)^2 and sets the pi on-shell.
The pi satisfy: p1 + p2 + p3 + p4 + p5 = 0.";

SetStandardMatrixElements::usage=
"SetStandardMatrixElements[ { expr1 }  -> { 1 },
{expr2} -> { 2 }, ... , enMomCon]
 is a function specifying which
abbreviations will be substituted for the standard
matrix elements made out of spinors, gamma matrices and
(eventually) scalar products containing the dependence on the
polarization vectors.  For enMomCon a rule indicating energy
momentum conservation may be given (e.g. {p4:>-p1-p2-p3}).
 Example: SetStandardMatrixElements[
{Spinor[p1].DiracSlash[k1].Spinor[p2] ScalarProduct[Polarization[k1],k2]}->
{1}, {p4:>-p1-p2-p3}].";

Negligible::usage=
"Negligible[me] is the head of small (negligible) masses.
This means that any mass with this head will be neglected if it
appears in a sum, but not as an argument of Passarino-Veltman  functions.";

NegligibleVariables::usage=
"NegligibleVariables is an option for OneLoop.
\"NegligibleVariables->{Melectron}\" i.e. will substitute \"Negligible[Melectron]\"
 for all Melectron's in the calculation.";

SpecificPolarization::usage=
"SpecificPolarization[exp, Polarization[r] -> {i, a, b},
Polarization[r2] -> {j, a2, b2}, ... ]
calculates exp with specifying the polarization of  the
momentum r to 0 (parallel), 1 (orthogonal), 2 (longitudinal),
\"+\" (righthanded), \"-\" (lefthanded).  The explicit polarization
vectors are constructed from r,a and b.";

Spinor::usage=
"Spinor[p, m, optarg] is the head of Dirac spinors.
Which of the spinors u, v,u_bar or v_bar
is understood, depends on the sign of the momentum (p) argument and
the relative position of DiracSlash[p]:
Spinor[sign p, mass]  is that spinor which yields
sign*mass*Spinor[p, mass] if the Dirac equation is applied.";

SpinorCollect::usage=
"SpinorCollect is an option for FermionSpinSum. If set to False the
 argument of FermionSpinSum has to be already collected w.r.t. Spinor.";

SpinSumExternalMomentum::usage=
"experimental";


SpinPolarizationSum::usage=
"SpinPolarizationSum is an option for SquareAmplitude. 
The set (pure) function  acts on the usual spin-sum.";

SquareAmplitude::usage=
"SquareAmplitude[amp] squares the amplitude amp. STILL EXPERIMENTAL!!!
Don' rely on it.";

StandardMatrixElement::usage=
"StandardMatrixElement[ ... ] is the head for matrixelemnt abbreviations.";

SU3Delta::usage= "SU3Delta is replaced by SUNDelta.";

SU3F::usage= "SU3F is replaced by SUNF.";

SUNNToCACF::usage= "SUNNToCACF is an option of SUNSimplify. If set to 
True, the Casismir operators CA (=N) and CF (=(N^2-1)/(2 N)) are introduced.";

SU3FToTraces::usage="SU3FToTraces is relaced by SUNFToTraces.";
 
SUNFToTraces::usage="SUNFToTraces is an option for SUNF and SUNSimplify, 
determining whether SUNF is expressed by traces.";

SUNFJacobi::usage="SUNFJacobi is an option for SUNSimplify, indicating 
whether the Jacobi identity should be used.";

SUNDelta::usage= "SUNDelta[i, j] is the Kronecker-delta for SU(N) with color
indices i and j."

SUNF::usage= "SUNF[i, j, k] are the structure constants of SU(N).";

SUNIndex::usage=
"SUNIndex is the head of SU(N)indices.";

SUNN::usage="SUNN denotes the number of colors. SUNDelta[a, a] yields
(SUNN^2 -1)."; 
 
SUNSimplify::usage="Simplifies products of SUNT (and complex conjugated)
matrices.";

SUNT::usage= "SUNT[a] is the generator of SU(N).";
 
SUNTrace::usage= "SUNTrace[expr] calculates the trace of expr.
All indices should occur twice and expr must be a product of SUNF's, SUNDelta's
and SUNT's.";

Trc::usage=
"Trc[exp] calculates the Dirac trace of exp.  Tr is identical to
DiracTrace, up to the default setting of DiracTraceEvaluate.";

TraceOfOne::usage=
"TraceOfOne is an option for Tr and DiracTrace.
Its setting determines the value of the unit trace.";

TrickMandelstam::usage=
"TrickMandelstam[expr, {s, t, u, m1^2 + m2^2 + m3^2 + m4^2}]
simplifies all sums in expr in such a way that one of the
Mandelstam variables s, t or u is eliminated by the
relation s + t + u = m1^2 + m2^2 + m3^2 + m4^2.
The trick is that the resulting sum has the most short number of terms.";

UVPartOnly::usage=
"
ROTTEN !!!! Don't use it.!!

UVPartOnly is an option for OneLoop. If set to True the result of 
OneLoop will be the coefficient of the UV-Divergence.";

Write2::usage=
" Write2[channel, val1 = expr1, val2 = expr2, ...] writes the settings
val1 = expr1, val2 = expr2 in sequence followed by a newline, to the
specified output channel.
Using Write of the Mathematica Version 2.0 (and 2.1) may on rare occasions
produce an output which cannot be read in identically. This bug is fixed by
Write2.";

WriteOutPaVe::usage=
"WriteOutPaVe is an option for PaVeReduce and OneLoopSum. If set to a string, the
results of all Passarino-Veltman PaVe's are written out.";

WriteOut::usage=
"WriteOut is an option for OneLoop. If set to True, the result of
OneLoop will be written to a file called \"name.res\", where name
is the first argument of OneLoop.";

$FeynC::usage = "If set to true the Thomas Grund - C -code for Collect2
is used (careful: Due to a bug in the Mathematica Kernel it is NOT
possible to send variables with local contexts to Collect2 and get them
back correctly).";                

$Kreimer::usage=
"experimental setup of the Kreimer-scheme for Gamma5.";

$Larin::usage=
"The deault is now (FeynCalc2.2beta1) that $Larin is False!!!!!!!!!!! \n
If set to True, the Larin-Gorishny-Akyampong-DelBurgo-scheme for
gamma5 in D-dimensions is used, i.e., before evaluating traces 
(but after moving gamma5 anticommuting in all dimensions to the 
right of the Dirac string) a product  gamma[mu].gamma5  is substituted
to  -I/6 Eps[mu,al,be,si] gamma[al,be,si], where all indices live in 
D-dimensions now. Especially the Levic-Civita tensor is taken to be
D-dimensional, i.e., contraction of two Eps's results in D's.
This has (FOR ONE AXIAL-VECTOR-CURRENT ONLY, it is not so clear if this
scheme also works for more than one fermion line involving gamma5)
the same effect as the Breitenlohner-Maison-'t Hooft-Veltman scheme.";


$LimitTo4::usage=
"$LimitTo4 is a global variable with default setting True.
If set to False no limit  Dimension -> 4 is performed after tensor
integral decomposition.";

$MemoryAvailable::usage=
"$MemoryAvailable is  a global variable which is set to an integer
n, where n is the available amount of main memory in MB.
The default is 16. It should be increased if possible.
The higher $MemoryAvailable can be,  the more intermediate steps do
not have to be repeated by FeynCalc.";

$MU::usage=
"$MU is the head for dummy indices which may be introduced by Chisholm.";

Protect[$MU];

$SpinorMinimal::usage=
"$SpinorMinimal is a global switch for an additional simplification
attempt in DiracSimplify for more than one Spinor-line.
The default is False, since otherwise it costs too much time.";

$ToughStuff::usage=
"$ToughStuff should be set to True if very lengthy 1-loop calculations
are done with OneLoopSum.";

$VeryVerbose::usage=
"$VeryVerbose  is a global variable with default setting
0. If set to 1, 2, ..., more and more
intermediate comments and informations
are displayed during calculations.";


Begin["Feyn`Calc`General`"];
(* *************************************************************************** *)
Factor1[x_] := Block[{factor1t1, factor1t2, factor1t3,mt,mi,m1,mp1,nx=x,iIii},
mt = (((# /. Plus -> mi /. mi -> Plus) /. m1 -> (-1)/.mp1 -> (-Plus[##]&)
      ) /. iIii -> I)&;
mi[y_, z__] := (m1 mp1[y,z] )/; (If[ Head[#] === Complex, False,
               If[# < 0, True, False]]& @ NumericalFactor[y]);
    nx = x /. Complex[0, b_] -> (b iIii);
    If[ Head[nx] =!= Plus, mt[nx /. Plus -> (Factor1[Plus[##]]&)],
        factor1t1 = Apply[ List, Map[# /. Plus -> factor1t3&, nx] ];
        factor1t2 = (PolynomialGCD @@ factor1t1) /. factor1t3 -> Plus;
        mt[(factor1t2 Apply[Plus,
                         Map[((# /. factor1t3 -> Plus) / factor1t2)&,
         factor1t1]])]
      ]               ];


(* dFRH *)
FRH[x_] := FixedPoint[ReleaseHold, x];
(* dCollect2 *)
Options[Collect2] = {Expanding -> True,
                     Factoring -> True, IsolateHead -> False } ;
Collect2[x_List, y__]      := Collect2[#, y]& /@ x;
Collect2[x_, y_, r___Rule] := Collect2[x, {y}, r] /; Head[y]=!=List;
(* Collect2[x_, y_List, ___]  := x /; FreeQ2[x, y]; *)

$FeynC = True;
Collect2[ expr_, vv_List,r___Rule ] := Block[
{exo,v,ru,nx,lk,fa,in,ih,pr,wr,tog,fr0,frx,lin,tv,mp,mp2,cd,i,co,ara,tim,
 new = 0,potrule,
 einss,re,compCON,ccflag = False, thc, ish},
{exo, fa, ih} = {Expanding, 
                  Factoring, IsolateHead} /. Join[{r}, Options[Collect2]];
v = Select[ vv, ((Head[#] =!= Plus) && (Head[#] =!= Times) && 
                 (!NumberQ[#]))& ];
If[Length[v] =!= Length[Union[v]], v = Union[v]];
(*check if the FastCalc-C-code by Th. Grund    is installed *)
If[(DownValues[Global`CCollect] =!= {}) && ($FeynC === True),
   potrule = { (a_Plus /; !FreeQ2[a, v])^n_ :>
                 Expand[a^n, y_ /; MemberQ[v, y]]
             };
   thc = Global`CCollect[expr /. potrule, v];
If[!FreeQ[thc, Global`Feyn], 
   input = expr;
   Print["CCt"];
$FeynC = False;
re = Collect2[expr, v, r];
$FeynC = True;
   ,
   If[fa =!= True, re = thc[[1]] + (thc[[2]] . thc[[3]]),
      print2["factoring constant part (C)"];
      re = Factor2[thc[[1]]];
      print2["factoring the rest (C)"];
      re = re + Map[Factor2, thc[[2]] . thc[[3]]];
     ];
      If[ih =!= False,
         ish[xx__] := Isolate[Plus[x], IsolateHead -> ih];
         re = Map[#/.Plus -> ish, re + nullll] /. nullll->0;
        ]

],
     

v = Select[ v, !FreeQ[expr, #]&];

If[Length[v] === 0,  re = expr,

tim = Timing[
nx = expr;
If[!FreeQ[nx, ComplexConjugate], 
   ccflag = True;
   nx = nx /. ComplexConjugate -> compCON;
   v = v /. ComplexConjugate -> compCON;
  ];

nx = nx/. HoldForm[k_[ii_]] -> lk[k][ii];
in = $VeryVerbose; If[!NumberQ[in], in = 0];
SetAttributes[{pr, wr}, HoldAll]; 
pr[i_, x__] := Print[x] /; in >= i;
wr[i_, x__] := WriteString["stdout", x] /; in >= i;
If[ fa === False,  tog[x_] := FixedPoint[ReleaseHold, x],
    fr0[x__] := Plus[x] /; !FreeQ2[{x}, v];
(*
    tog[x_]  := Factor2[FixedPoint[Combine[ReleaseHold[#]] &, x]];
*)
    tog[x_]  := Factor2[FRH[x]];
    frx[x__] := HoldForm[Plus[x]];

    nx = nx /. Plus -> fr0 /. fr0 -> frx 
  ];                                                   
If[exo === True, 
 wr[2,"expanding. "]; 
nx  = Expand2[nx, v]; (* lin denotes the part free of v *)
  ];
lin = Select[nx + ze1 + ze2, FreeQ2[#, v]&] /. ze1 -> 0 /. ze2 -> 0;
nx  = nx - lin;
If[fa === True,    wr[2, "linear part; LeafCount = ", LeafCount[lin]];
   lin = tog[lin]; wr[2, "; factored. "] 
  ];
tv = {}; (* tv is the list of all monomials to collect *)
mp[x_] := Block[{t1, t2},  (* "tv" is calculated as a side effect ! *)
                If[FreeQ2[x, v], x, t1 = Select[x t2, !FreeQ2[#, v]&];
                If[!MemberQ[tv, mp2[t1]], AppendTo[tv, mp2[t1]] ];
                (Select[x t2, FreeQ2[#, v]&]/t2) mp2[t1] ] ] (*endBlock*);
(*
                x/t1 mp2[t1] [] [] (*endBlock*);
*)
nx = (mp /@ (nx + ze) ) /. ze -> 0 /. mp -> mp2; pr[2,"length ",nx//Length,"."];
(* In case of denominators containing variables to be collected *)
cd[x_] := ((Numerator[#]/(Factor2[Denominator[#]] /.
   Plus-> (Collect2[Plus[##], v, r]&)))& @ x ) /; 
   (!FreeQ[Denominator[x], Plus]) && (!FreeQ2[Denominator[x], v]); 
If[Length[tv]>1, pr[2, "collecting ",Length[tv], " terms."]];
For[ i = 1, i <= Length[tv], i++, wr[2, "#",i];
     co = tog[ Coefficient[ nx, tv[[i]] ] ];
     If[Head[co] === Plus, co = tog[einss co] ];
     nx = nx /. tv[[i]] -> 0;
     If[ ih =!= False, 
         co = Isolate[co /. lk[ka_][j_] -> HoldForm[ka[j]],ara , 
                      IsolateHead -> ih];
         new = new + ((Isolate[FixedPoint[ReleaseHold, tv[[i]] /.
                              lk[ka_][j_] -> HoldForm[ka[j]]] /.  mp2 -> cd /. 
                              cd -> Identity, v, IsolateHead -> ih 
                              ] * co) /. einss -> 1) ,
         new = new + (((FixedPoint[ReleaseHold, tv[[i]]] /. mp2 -> cd /. 
                                cd -> Identity ) * co)/.einss->1)
   ]  ]    ][[1]];
wr[2,".\n"]; pr[2, "collected. time needed = ", tim //FeynCalcForm ];
If[ ih =!= False, 
    lin = Isolate[ FixedPoint[ReleaseHold, lin], v, IsolateHead->ih ],
    lin = FixedPoint[ReleaseHold, lin] ];
re = ((nx + new + lin) /. lk[ka_][j_] -> HoldForm[ka[j]] /. frx->Plus);
If[ccflag, re = re /. compCON -> ComplexConjugate];
  ](*endIf*);
    ] (* endifDownValues*);
re];
(* *************************************************************************** *)
(* *************************************************************************** *)
(* dCombine *)
Options[Combine] = {ProductExpand -> False };
Combine[x_, r___Rule] := Block[{h, new},
new = Together[ x /. Plus -> (If[FreeQ[{##}, _^_?Negative] && 
                                 FreeQ[{##}, Rational], h[##], Plus[##]
                                ]&) 
              ] /. h -> Plus;
If[(ProductExpand /. {r} /. Options[Combine]),
   new = ExpandNumerator[new//ExpandDenominator] 
  ];                       new];
(* *************************************************************************** *)
(* *************************************************************************** *)
SetAttributes[timeconstrained, HoldAll];
If[$OperatingSystem === "Unix", 
   timeconstrained[x__] := TimeConstrained[x],
    timeconstrained[x_,__] := x
  ];
(* dFactor5*)
Factor5[xy_] := If[Head[xy] === Plus,
                        Factor2[Map[# /. Plus -> plHF&, xy]] /. plHF -> Plus,
                        If[Head[xy] === Times, Map[Factor5, xy], xy]
                  ];

(* dFactor2 *)
Options[Factor2] = {FactorFull -> False, FactorTime -> 3600};
Factor2[x_, r___Rule] := Block[{fc,mt,mi,m1,mp1,cm,ff,pr,pp,h,fc5,tx,factor55,pl5,iIi},
fc5[y_] := If[Head[y] === Times, Map[fc5, y], 
                  If[Head[y] === Power, fc5[y[[1]]]^y[[2]],
(* change 1.9.93: put an Expand in ... *)
                     Factor[Expand[y]]]];
fc = mt[pp[ff[pp[Numerator[#]]] / ff[pp[Expand2[Denominator[#], Power]]]]]&;
mt = ((# /. Plus -> mi /. mi -> Plus) /. m1 -> (-1) /. mp1 -> (-Plus[##]&)) &;
mi[y_, z__] := (m1 mp1[y,z] )/; (If[ Head[#] === Complex, False,
               If[ # < 0, True, False] ]& @ NumericalFactor[y]);
If[$VersionNumber > 2.0, 
ff = timeconstrained[fc5[#], FactorTime /. {r} /. Options[Factor2], # ] &,
ff = timeconstrained[fc5[# /. Complex[0, b_] -> (b iIi) ]/. iIi -> I,
                         FactorTime /. {r} /. Options[Factor2], # ] &
  ];

pr = { fa_. pc[a_, b_]^n_. pc[a_, c_]^n_. :> 
            (fa pc[a^2, -b^2]^n) /; (((b + c) === 0) && IntegerQ[n]),
       fa_. pc[a_, b_]^n_. pc[c_, b_]^n_. :> 
            (fa pc[b^2, -a^2]^n) /; (((a + c) === 0) && IntegerQ[n]) };
pp = If[(FactorFull/. {r} /. Options[Factor2] ) =!= True,
        ( ((Numerator[#] /. Plus -> pc) //. pr) /. pc -> Plus /. pc -> Plus)/
        ( ((Denominator[#] /. Plus -> pc) //. pr) /. pc -> Plus /. pc -> Plus),
     #]&; 
(*
tx = vsu[x, 42000];
factor55[z_] := If[Head[z] === Plus,
                        fc[Map[# /. Plus -> pl5&, z]] /. pl5 -> Plus,
                        If[Head[z] === Times, Map[factor55, z], z]
                  ];
*)

(*
(fc[Combine[factor55[tx[[1]]], ProductExpand -> False]
   ] /. tx[[2]])
*)
(*
fc[Combine[x, ProductExpand -> False]]
*)
fc[x] ];
vsu[y_, les_:50000] := Block[{vv, xX, vs, iv, yr, vb, ly = LeafCount[y]},
               If[ly > les, print2["leafcount in vsu = ",ly];yr = y; vb = {},
                  vv = Variables[y];
                  vs = Table[vv[[iv]] -> ToExpression[StringJoin["xX",
                                       iv//ToString]], {iv, Length[vv]}];
                  vb = Map[Reverse, vs];
                  yr = y /. vs
                 ];
                    {yr, vb}];
                 
      
(* *************************************************************************** *)

(* *************************************************************************** *)
(* *************************************************************************** *)
(* dFreeQ2 *)
FreeQ2[_,{}]          := True;
FreeQ2[x_, y_]        := FreeQ[x, y] /; Head[y] =!= List;
FreeQ2[x_, {y_}]      := FreeQ[x, y];
FreeQ2[a_, {b_, c__}] := If[FreeQ[a, b], FreeQ2[a, {c}], False];
(* *************************************************************************** *)
(* *************************************************************************** *)
(* added August 1994; this fixes a terrible bug in Mathematica 
  (Versions 1.2 - 2.2)
*)
(* dExpand2 *)
Expand2[a_] := Expand[a];
Expand2[x_, a_ /; Head[a] =!= List] := Expand2[x, {a}];
Expand2[x_, l_List] := If[FreeQ[x, Plus], x,
  Block[{pl, t, plus},
(*
        SetAttributes[pl, HoldAll];
        SetAttributes[plus, HoldAll];
*)
        pl[y__] := If[FreeQ2[{Hold[y]}, l], plus[y], Plus[y]];
        t = Expand[x /. Plus -> pl]//.plus->Plus /. pl -> Plus;
        t
       ]                 ];

(* added by TH 23/12/94 *)
StandardOrder[e_]:=Block[{mi,ee},
  mi[y_, z__] := (m1 mp1[y,z] ) /; (If[ Head[#] === Complex, False,
                 If[ # < 0, True, False] ]& @ NumericalFactor[y]);
  ee=((e /. Complex -> mi /. mi -> Complex) /.
     m1 -> (-1) /. mp1 -> (-Complex[##]&));
  ee=((ee /. Plus -> mi /. mi -> Plus) /.
     m1 -> (-1) /. mp1 -> (-Plus[##]&))
];

(* dIsolate *)
 Options[Isolate] = {IsolateHead -> K, IsolateSplit -> 442};
 Isolate[y_HoldForm^n_., ___] := y^n;
(* this gives Problems if x has large HoldForm's ...
 Isolate[n_?NumberQ x_, y__]  := n Isolate[x, y];
*) (* better : *)
 Isolate[x_ /; NumericalFactor[x] =!=1, y__
        ] := NumericalFactor[x] Isolate[x/NumericalFactor[x], y]
 Isolate[x_?NumberQ, __]      := x;
 Isolate[x_, ___]             := x /; Length[x] < 1;     
 Isolate[ex_, r___Rule ]      := Isolate[ex, {}, r];
 Isolate[ex_, var_, r___Rule] := Isolate[ex, {var}, r]/;
                                (Head[var] =!= Rule) && Head[var] =!= List;

 Isolate[ exp_ /; Apply[Or[#===1, #===0]&, {NumericalFactor[exp]}], 
vars_List, ops___Rule ] := Block[{plush,vlist,res,split,kk},
      kk = IsolateHead/.{ops}/.Options[Isolate];
      split = IsolateSplit/.{ops}/.Options[Isolate];
      vlist = Flatten[{vars}];
(* This split-off is useful for various reasons (continuation lines, ...) *) 
      plush[x__] := If[ !FreeQ2[{x}, vlist], Plus[x],
         If[ (checkIsolate[x, split] === True ) && Length[{x}] > 4,
             Isolate[ Drop[Plus[x], Round[Length[Plus[x]]/2]] + 
               Isolate[Take[Plus[x], Round[Length[Plus[x]]/2]], vars, ops],
                                 vars, ops ], remIsolate[Plus[x], kk]
           ]          ];
          res = StandardOrder[exp] /. Plus -> plush /. plush -> Plus;
      If[Head[res] =!= HoldForm && vlist === {}, res = remIsolate[res, kk]];
                                            res](*endIsolate*);
(* three extra "global" functions *)
checkIsolate[x__, i_] := If[Head[i] === Integer,  (* LGF *)
 If[Length[Characters[ToString[FortranForm[Plus[x]]]]] > i, True, False],
 If[Head[i] === Complex, If[Length[{x}] > Im[i], True, False, False], False]];

tokIsolate[y_, ab_, uh_] := ab[ToExpression[ StringDrop[ToString[y],
                            StringLength[ToString[uh]]] ]];
remIsolate[x_, abb_] := remIsolate[x, abb] = Block[{temp, re, h, set},(*LGF*)
              If[ Head[abb]===Symbol, 
                  temp = tokIsolate[ uni[ToString[abb]]/.uni->Unique,abb, abb],
                  temp = tokIsolate[ Unique["dude"], abb, "dude" ] ];
              (* Need this in case old isolated K's are loaded *)
              If[ Head[temp]=!=abb,
                  re = Select[DownValues[abb], (#[[2]]===x) &];
                  If[Length[re] > 0, re = re[[1,1]] /. Literal -> HoldForm,
                     re = HoldForm @@ {temp}; Set@@{temp, x}
                    ],
                    re = HoldForm @@ {temp};
                    Set@@{temp, x} ];          re];
(* *************************************************************************** *)
(* *************************************************************************** *)
(* dNumericalFactor *)
NumericalFactor[x_] := 
If[NumberQ[x], x, If[Head[x] === Times, If[NumberQ[First[x]], First[x], 1], 1]];
(* *************************************************************************** *)
(* *************************************************************************** *)
(* dPartitHead *)
   PartitHead[x_, y_]      := {1, x} /; Head[x] === y;
   PartitHead[x_Times, y_] := {x, 1} /; FreeQ[x, y];
   PartitHead[x_, y_]      := {x, 0} /; FreeQ[x, y]; 
   PartitHead[x_Plus, y_]  := {#, x - #}& @ Select[x, FreeQ[#, y[___]]&];
   PartitHead[x_Times,y_]  := {x/#, #}& @ Select[x, If[Head[#]===y,True]&]; 
(* *************************************************************************** *)

(*dWrite2*)
 Options[Write2]={ FinalSubstitutions -> {},
                   FormatType -> InputForm, D0Convention -> 0, 
                   PageWidth    -> 50,
                   PreFortranFile -> "", PostFortranFile -> ""
                 };

 SetAttributes[ Write2, HoldAll ];
 Write2[f_String, x___, l_] := 
  Write2[f, Hold[x, l], dummyrule->False ]/; FreeQ[Hold[l], Rule];

 Write2[file_String, eeq__, opts___Rule] := Block[{j,vhf,vv,eq,k2str,ops,ide,
  aa0, be00, be11,be0, db0, ce0, de0, ansg,d0convention,oldopenops,pww,
  prefortran, postfortran, pagewidth,prerec,tostring,flag},

ops         = Join[{opts}, Options[Write2]];
{finsubst, pagewidth } = {FinalSubstitutions, PageWidth} /. ops;
{prefortran, postfortran}  = Flatten /@ {{PreFortranFile}, {PostFortranFile}} /. ops;
(* a modified Power function, for avoiding Fortran-Complications *)
pww[x_?NumberQ, 1/2]          := Sqrt[N[x]];
pww[x_?NumberQ, rat_Rational] := Power[N[x], N[rat]];
pww[x_,1/2]                   := Sqrt[x];
pww[x_, rat_Rational]         := Power[x,N[rat]];
pww[x_, he_]                  := (x^he) /; Head[he]=!=Rational;

If[finsubst === {},
   finsubst1 = FinalSubstitutions  /. Options[OneLoop];
   finsubst  =  FinalSubstitutions /. Options[OneLoopSum];
   If[finsubst==={}, finsubst = finsubst1];
   If[Head[finsubst]=!= List, finsubst = {}];
  ];

{aa0,be0,be1,be00,be11,db0,ce0,de0} = {A0,B0,B1,B00,B11,DB0,C0,D0}/.finsubst;
 (* allvar gives all Variables in HoldForm,( K[i] ) *)
allvar[y_] := Block[{arr={},ia,new,alt=Drop[#, -1]& /@ Position[y,HoldForm]},
                     For[ia = 1, ia <= Length[alt], ia++,
                         new = Part @@ Prepend[alt[[ia]], y];
                         If[!MemberQ[arr, new], AppendTo[arr,new] ] ];
               arr];
ide = {##}&;
eq = Flatten[{Hold[{eeq}]} /. Set -> Equal /. Hold -> ide];
mrel[x_] := MapAll[ReleaseHold, x];
(* vhf gives all "K" which are really present *)
vhf[n_. y_HoldForm]:= Block[{kk, qq}, kk = y[[1, 0]];
  (Table[ HoldForm @@ {qq[ii]}, {ii, y[[1,1]]} ] /. qq -> kk)/.finsubst
                           ] /; NumberQ[n];
vhf[y_] := Block[{te=y, var={}}, 
   While[!FreeQ[te, HoldForm], var = Union[ var, allvar[te] ];
         te = ReleaseHold[te//Variables]; var ]; 
   var/.finsubst];

If[(FormatType/.ops/.Options[Write2]) === FortranForm,
               oldopenops = Options[OpenWrite];
               SetOptions[OpenWrite, FormatType->FortranForm, 
                                     PageWidth-> pagewidth ];
                  WriteString[file,
"C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"];
                  Write[file];
                  WriteString[file,"C      ", ToString[Date[][[3]]],".",
                              ToString[Date[][[2]]], ".",ToString[Date[][[1]]]
                             ];
                  Write[file];
                  WriteString[file,
"C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"];
                  Write[file];
  ];
tostring = If[Head[#] === String, #, ToString[#]]&;
    For[ j=1, j<=Length[eq], j++,
         If[ !FreeQ[ eq[[j,2]], HoldForm ], vv = vhf[eq[[j,2]]], vv = {} ];
          eqj1 = eq[[j,1]]; 
          If[(FormatType/.ops/.Options[Write2]) === FortranForm,
             eqj2 = eq[[j,2]] /. Power->pww;
             If[Head[eqj1] === FUNCTION,
                WriteString[file, "      FUNCTION ", eqj1[[1]]//tostring,
                                  "()"];
                Write[file];
                eqj1 = eqj1[[1]]
               ];

               If[prefortran =!= {""},
                  prerec = Flatten[ReadList[#, Record]& /@ prefortran];
                  flag = False;
                  For[iir = 1, iir <= Length[prerec], iir++,
                      If[flag =!= True, 
                         If[!StringMatchQ[StringJoin @@ 
                                          Drop[prerec,iir-1], "*IMPLICIT*"],
                            If[Length[eqj1]===2, WriteString[file, eqj1[[2]]],
                               WriteString[file, "      REAL*8"]
                              ];
                            WriteString[file, "      ",eqj1, "\n"];
                            If[Length[vv] > 0,
                               vardec = StringJoin["      COMPLEX*16  ",
                                                   ToString[vv[[1,1,0]]],
                                                   "(",
                                                   ToString[Length[vv]+42],
                                                   ")\n"];
                               WriteString[file,vardec ];
                              ];
                            flag = True
                           ];
                        ];
                      WriteString[file, prerec[[iir]]];
                      Write[file];
                     ]
                 ],
              eqj2 = eq[[j,2]]
            ];

         Which[ 
               (FormatType/.ops/.Options[Write2]) === InputForm,
                     If[FreeQ[Streams[], file],
                        OpenWrite[file, FormatType -> InputForm]
                       ];
                     mal[x_]:=(True/;Length[Characters[ToString[ 
                                             InputForm[x]]]]<74
                              ) /; Length[x]<22;
                     If[ !FreeQ[ eqj2, HoldForm ],
                         For[iv=1, iv<=Length[vv], iv++,
                             If[mal[vv[[iv]]//ReleaseHold]=!=True,
                                WriteString[file, 
                                   ToString[vv[[iv]]], " = ( "],
                                WriteString[file, 
                                   ToString[vv[[iv]]], " = "] 
                               ];
                             Write[file, ReleaseHold[ vv[[iv]] ] ];
                             If[mal[vv[[iv]]//ReleaseHold]=!=True,
                                WriteString[file, "       );\n"]
                               ]
                            ] 
                       ];(* Write[file];*)
                     If[mal[eqj2]=!=True,
                        WriteString[file, eqj1//InputForm, " = ( " ],
                        WriteString[file, eqj1//InputForm, " = "]
                       ];
                     Write[ file, eqj2 ];
                     If[mal[eqj2]=!=True,
                        WriteString[file, "       );\n"]
                       ],

               (FormatType/.ops/.Options[Write2]) === FortranForm,
               oldopenops = Options[OpenWrite];
               SetOptions[OpenWrite, FormatType->FortranForm,
                                     PageWidth-> pagewidth
                         ];

               d0convention = D0Convention /. ops /. Options[Write2];
               If[ d0convention === 0, 
                   ansg[x_] := x/. 0 -> Null/.  finsubst
                 ];
If[ d0convention === 1, 
    ansg[v_. x_]:= (v x)/; FreeQ2[(v x)/.finsubst, {de0, ce0, be0, aa0, db0}];
    ansg[v_. x_] := Block[{args,t4,t5,t6,ll}, args = Apply[List, x];
         t4 = Take[args,4]; t5 = args[[5]]; t6 = args[[6]]; 
         ll = PowerExpand[ Sqrt[Take[args,-4]] ] /. finsubst;
         (v (Apply[de0, Join[t4, {t6,t5}, ll]]) /. 0 -> Null) ] /;
         ( (Head[x/.finsubst]===(de0)) && (Head[v/.finsubst] =!= (de0)) );
    ansg[v_. x_]:=Block[{args, mm}, args = List @@ x;
         mm = PowerExpand[ Sqrt[Take[args,-3]] ] /. finsubst;
         (v (ce0@@Join[Take[args,3], mm])/. 0 -> Null)] /;
          ( (Head[x/.finsubst]===(ce0)) && (Head[v/.finsubst] =!= (ce0)) );
    ansg[x_]:=Block[ {args, mm}, args = List@@x;
         mm = PowerExpand[ Sqrt[Take[args,-2]] ] /. finsubst;
         ((be0@@Join[{args[[1]]}, mm])/. 0 ->Null)]/;
                      Head[x/.finsubst]===(be0);
    ansg[v_. x_]:=Block[{args, mm}, args = List@@x;
         mm = PowerExpand[ Sqrt[Take[args,-2]] ] /. finsubst;
         (v (db0@@Join[{args[[1]]}, mm])/. 0 -> Null)] /; 
         ( (Head[x/.finsubst]===(db0)) && (Head[v/.finsubst] =!= (db0)) );
    ansg[x_]:=Block[{mm}, mm = PowerExpand[ Sqrt[x] ] /. finsubst;
                    aa0[mm] ] /; Head[x/.finsubst]===(aa0);
  ];
If[ !FreeQ[ eqj2, HoldForm ],
    For[iv=1, iv<=Length[vv], iv++,
        WriteString[file,"      ", (vv[[iv]])//FortranForm," = "];
(* assuming that a function Negligible in a Fortran file does not make sense ... *)
        If[!FreeQ2[{ be0, be1, be00, be11, db0, ce0, de0 }, 
                    Map[Head, Select[ Variables[ReleaseHold[ vv[[iv]] ] 
                                             ]/.finsubst,
                                   !FreeQ2[{be0, be1, be00, be11, db0, ce0, de0
                                           }, Head[#] ]& ]
                       ]
                  ],
        Apply[WriteString,{file, StringReplace[ToString[FortranForm[ansg[
                            ( ReleaseHold[vv[[iv]]]/.finsubst ) /.
                             Negligible->Identity/.
                               Power->pww ] /. ansg->Identity]],
                                           {"Null" -> "0D0"}]
                          }];
        Write[file], Write[file, ReleaseHold[vv[[iv]]]/.      
                           Negligible->Identity/.Power->pww ];
           ];
       ]
  ];
  WriteString[file, "      ",FortranForm[eqj1]," = "];
  Write[file, ansg[(eqj2/.Negligible-> Identity/. Power->pww )/.finsubst
                  ] /. ansg -> Identity];
  SetOptions @@ Prepend[oldopenops, OpenWrite]
          ](* endWhich *)
      ]; (* end j - loop *)

If[(FormatType/.ops/.Options[Write2]) === FortranForm,
   If[postfortran =!= {""},
      prerec = Flatten[ReadList[#, Record]& /@ postfortran];
      For[iir = 1, iir <= Length[prerec], iir++,
          WriteString[file, prerec[[iir]]]; Write[file]
         ]; Write[file]
     ]
 ];
Close @@ {file};
 (* for Fortran: check if no line is larger than 72 columns, and merge 
    lines, if possible  *)
If[(FormatType/.ops/.Options[Write2]) === FortranForm,
   rfile = ReadList[file, Record];
   OpenWrite @@ {file};
   If[ OddQ[ Length[rfile] ],
       rfile = Append[rfile, "                  "]
     ];
   For[ir = 1, ir < Length[rfile], ir = ir + 2,
          joinlabel = False;
          rf1 = StringReplace[rfile[[ir]],   {"\\"->""}];
          rf2 = StringReplace[rfile[[ir+1]], {"\\"->""}];
          If[(StringLength[rf1] > 72) || (StringLength[rf2]) > 72,
             Print["FORTRAN WARNING! Too long line encountered. Check "]
            ];
          If[ (StringLength[rf1] > 1) && (StringLength[rf2] > 7),
              If[StringTake[rf2, {6, 8}] === "-  ",
                 rf3 = StringDrop[rf2, 8];
                 If[StringTake[rf1, -1] === " ", rf1 = StringDrop[rf1,-1]];
                 If[StringLength[rf1] + StringLength[rf3] < 72,
                    joinlabel = True;
                    rf1 = StringJoin[rf1, rf3];
            ]   ]  ];
          If[joinlabel === True,
             WriteString[file, rf1, "\n" ],
             WriteString[file, rf1, "\n", rf2, "\n" ]
            ]
      ];
Close @@ {file};
 ];
file];
WriteString["stdout", "."];
End[];

(* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ *)
Begin["Feyn`Calc`Main`"];
(* changes to Mathematica - functions *)

Unprotect[ExpandAll];
ExpandAll[x_] := Expand[x] /; (!FreeQ[x,Hold]) || !FreeQ[x,HoldForm];
Protect[ExpandAll];
(* Fix a strange behaviour of NumberQ *)
Unprotect[NumberQ];
NumberQ[Power[a_?NumberQ,_]]:=True;
Protect[NumberQ]
Format[Continuation[_], StringForm] := "";
Format[Continuation[_]] := "";
StringBreak[_] := "";
SetOptions[ToString, PageWidth -> 62 ];
If[!NumberQ[Sqrt[3]], Unprotect[NumberQ];
   NumberQ[Power[a_?NumberQ,_]]:=True; Protect[NumberQ]
  ];

SetAttributes[timeconstrained, HoldAll];
If[$OperatingSystem === "Unix",
   timeconstrained[x__] := TimeConstrained[x],
    timeconstrained[x_,__] := x
  ];

 
(* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ *)

(* #################################################################### *)
(*                             Main10                                  *)
(* #################################################################### *)

$ToughStuff = False;
$VeryVerbose = 0;
(* ************************************************************************ *)
$MemoryAvailable = 16;

(* dotLindef : linearity of Dot-products *)
   dotLin[x_] := x /. Dot -> dotlindl /. dotlindl -> Dot;     
   dotlindl[] = 1;
   dotlindl[ a___, b_ c_, d___ ] := b dotlindl[a, c, d] /; 
                                        (noncommQ[b] === True);
   dotlindl[ a___, b_, d___ ]    := b dotlindl[a, d] /; 
                                        (noncommQ[b] === True);
   dotlindl[ a_Spinor, b___, c_Spinor, d_Spinor, e___, f_Spinor, g___]:=
   dotlindl[ a,b,c ] dotlindl[d,e,f,g];

(* relhdef *)
   relh[x_]:= FixedPoint[ReleaseHold,x];
(* fr567def, frlivcdef : two special FreeQ - checking functions *)
   fr567[x__]:=True/;FreeQ2[relh[{x}],
                            {DiracGamma[5],DiracGamma[6],DiracGamma[7]}];
   frlivc[x_?NumberQ]:=True;
   frlivc[x_,y__] := True/;FreeQ2[relh[{x,y}],{Momentum, LorentzIndex}];
(* Need it this way, since Eps may contain Momentum or LorentzIndex *)
   frlivc[x_]     := True/;(Head[x]=!=Momentum) && 
                           (Head[x]=!=LorentzIndex);

(* gamma67backdef: reinsertion of gamma6 and gamm7 *)
   gamma67back[x_] := x/.DiracGamma[6]->( 1/2 + DiracGamma[5]/2 )/.
                         DiracGamma[7]->( 1/2 - DiracGamma[5]/2 );

(* memsetdef : a dynamical memory dependent "Set" function *)
 SetAttributes[memset,HoldFirst];
   memset[x_,y_]:=Set[x,y]/;( (LeafCount[y]<50000) &&
                              (($MemoryAvailable- MemoryInUse[]/1000000.)>1.)
                            );
   memset[x_,y_]:=y/;((LeafCount[y]>=50000) || 
                      (!TrueQ[($MemoryAvailable- MemoryInUse[]/1000000.)>1.])
                     );

(* noncommQdef : checking non-commutativity *) 
   noncommQ[z_]:=memset[ noncommQ[z], TrueQ[noncQ[z]] ];
   noncQ[x_ ?NumberQ]:=True;
   noncQ[x_SUNTrace]:=True;

   noncQ[x_] := If[FreeQ2[relh[x],{DiracGamma,Spinor, SUNT,
                                      DiracGammaT, ChargeConjugationMatrix,
                                      ChargeConjugationMatrixInverse}],
                  True, False];

(* #################################################################### *)
(*                             Main11                                  *)
(* #################################################################### *)

(* print1def, print2def, print3def: print functions *)
 SetAttributes[{print1, print2, print3}, HoldAll];
 print1[x__]:=Print[x]/;$VeryVerbose>0;             
 print2[x__]:=Print[x]/;$VeryVerbose>1;
 print3[x__]:=Print[x]/;$VeryVerbose>2;

(* timefixdef : a more physics - like timing function *)
   timefix[n_]:= Which[ 0.<=n<0.5,   "time < 0.5 s",
                        0.5<=n<9.5,  N[n,2] "s",
                        9.5<=n<59.5,  N[n,2] "s",
                        59.5<=n<600,  N[n/60,2] "min",
                        600<=n<3570,  N[n/60,2] "min",
                        3569<n<36000, N[n/3600,2] "h",
                        36000<n,      N[n/3600,4] "h"
                      ];
(* FeynCalcFormdef + formatting stuff for FeynCalcForm: *)
   su3fuser[a_,b_,c_,___]:=fsu3U[a, b, c]/.fsu3U->"f";
   sumst[x_Plus]:=SequenceForm["(",x,")"];  sumst[y_]:=y;
   did[x_,___]:=x;
   fcdot[x__]:=Dot[x]/;FreeQ2[{x},{DiracGamma,Spinor,SUNT}];
   fcdot2[a_]:=a;
   Format[fcdot2[a_,b__]] := Infix[fcdot2[a,b], " ", 210, None];
   diF[x_Symbol -4]:=StringJoin[ToString[x],"-4"];   diF[x_]:=x;

   dea[yy__]:=Map[denfa,{yy}];
   denfa[_[x_,0]]:= SequenceForm["(",x^2,")"];
   denfa[_[x_,y_]]:= SequenceForm["(",x^2,"- ",y^2,")"];
   feynden[x__]:=1/fcdot2@@( dea @@ {x} );
   ditr[x_,___]:="tr"[x];
   fdprop[a_,b_]:=1/denfa[dudu[a,b]];
   compind[a_] := If[Head[a] === Symbol, StringJoin[ToString[a],"*"], a "*"];

(* change as a side effect the ordering Attribute of Plus and Times,
   but reinstall it again at the end. *)
 iDentity[a_,___] := a;
  epsd[a:_[_]..] := "eps"[a];
  epsd[a___, b_[c_,di_], d___] := Subscripted["eps"[di//diF]][a,b[c,di],d];
 FeynCalcForm[x_]:=Block[{xxxx},
                     Unprotect[Plus, Times];
                     ClearAttributes[Plus,Orderless];
                     ClearAttributes[Times,Orderless];
                    xxxx = (x/.
         (n_Real Second)->timefix[n]/.
         SUNTrace -> "tr" /.
         Dot->fcdot/.fcdot->fcdot2/.
         Eps->epsd/.
         Pair[ LorentzIndex[v_],LorentzIndex[w_] ]:>
         "g"[v, w]/.
         Pair[ LorentzIndex[v_,di_],LorentzIndex[w_,di_] ]:>
         Subscripted["g"[di//diF]][v, w]/.
         Pair[ Momentum[v_],Momentum[w_] ]:>
         SequenceForm@@Flatten[ {v//sumst ,{"."},w//sumst} ]/.
         Pair[ Momentum[v_,di_],Momentum[w_,di_] ]:>
         Subscripted[
         (SequenceForm@@Flatten[ {v//sumst ,{"."},w//sumst} ])[di//diF]]/.
         Pair[ Momentum[Polarization[v_,-I,sun___]],LorentzIndex[w_] ]:> 
           ("ep(*)"[v,w,sun]/.SUNIndex->iDentity)/.
         Pair[ Momentum[Polarization[v_,-I,___]],LorentzIndex[w_] ]:> 
           "ep(*)"[v,w]/.
         Pair[ Momentum[Polarization[v_,I,sun___]],LorentzIndex[w_] ]:> 
           ("ep"[v,w,sun]/.SUNIndex->iDentity)/.
         Pair[ Momentum[v_],LorentzIndex[w_] ]:>
         SequenceForm@@Flatten[ {sumst[v],"[",w,"]"} ]/.
         Polarization[ka_,-I,___]->"ep(*)"[ka]/.
         Polarization[ka_,I,___]->"ep"[ka]/.
         ComplexIndex -> compind /.
         Pair[ Momentum[v_,di_],LorentzIndex[w_,di_] ]:>
         Subscripted[ sumst[v][di//diF] ][w]/.
         DiracGamma[ LorentzIndex[v_] ]:>ToString["ga"[v]]/.
         DiracGamma[ LorentzIndex[v_,di_],di_]:>
         Subscripted[ ToString["ga"][di//diF] ][v]/.
         DiracGamma[ Momentum[v_] ]:>ToString["gs"[v]]/.
         DiracGamma[ Momentum[v_,di_],di_]:>
         Subscripted[ ToString["gs"][di//diF] ][v]/.
         DiracGamma[5]:>"ga[5]"/.DiracGamma[6]:>"ga[6]"/.
         DiracGamma[7]:>"ga[7]"/.  DiracGamma:>"ga"/.
         Literal[Spinor[p_,0,2,___] ]:>"u"[p/.Momentum->iDentity]/.
         Literal[Spinor[p_,0,1,___ ]]:>"u"[p/.Momentum->iDentity]/.
         Literal[Spinor[- p_, mass_,2, ___ ]]:>
         ToString[ "v"[p/.Momentum->iDentity,mass] ]/.
         Literal[Spinor[ p_,mass_,2, ___ ]]:> 
         ToString["u"[p/.Momentum->iDentity,mass] ]/.
         Literal[Spinor[ -p_, mass_,1,___ ]]:> 
         ToString[ "v"[p/.Momentum->iDentity,mass] ]/.
         Literal[Spinor[p_,mass_,1,___ ]]:> 
         ToString[ "u"[p/.Momentum->iDentity,mass] ]/.
         DiracTrace->ditr/.  
         SUNF -> su3fuser /. SUNDelta -> "d"/. SUNT -> "T" /.
         SUNIndex -> iDentity /.
         FeynAmpDenominator->feynden /.  PropagatorDenominator-> fdprop/.
          Momentum->did/.LorentzIndex->did);
          SetAttributes[Plus,Orderless];
          SetAttributes[Times,Orderless];
          Protect[Plus, Times];
         xxxx];
	(* inserted 22jun98 th: *)
Unprotect[$OutputForms];
AppendTo[$OutputForms, FeynCalcForm];
Protect[$OutputForms];
         

(* #################################################################### *)
(*                             Main12                                  *)
(* #################################################################### *)
(* Momentumdef,  LorentzIndexdef, Polarizationdef *)
(* for the moment *)
(*
 Polarization[a_, 1] := Polarization[a, I];
 Polarization[a_,-1] := Polarization[a,-I];
*)
 Polarization[-x_Symbol, I] := -Polarization[x,I];
 Polarization[-x_Symbol,-I] := -Polarization[x,-I];

 SetAttributes[ {DiracGamma, LorentzIndex, Momentum, Pair}, Constant ];
 Momentum[Momentum[x__]] := Momentum[x];
(* for the funny new FA2.0 - convention ... *)
 LorentzIndex[i_Integer] := LorentzIndex[i] = Unique["li"];
 LorentzIndex[LorentzIndex[in__]]:=LorentzIndex[in];
 Momentum[x_,4]     := Momentum[x];           
 LorentzIndex[x_,4] := LorentzIndex[x];       
 Momentum[0] = 0;
 Momentum[0,_] := 0;
 Momentum[_,0]     := 0;           
 LorentzIndex[_,0] := 0;       

(* #################################################################### *)
(*                             Main13                                  *)
(* #################################################################### *)
 Literal[ rev[yz__ /; FreeQ2[{yz}, {SUNT}] ] ]:=
Isolate[
Dot @@ (Reverse[FRH[{ yz }]]/.
                    DiracGamma[5]->(-DiracGamma[5])/.
                   {DiracGamma[6] :> DiracGamma[7], 
                    DiracGamma[7]:>DiracGamma[6]}/.
                     ChargeConjugationMatrix->
                    (-ChargeConjugationMatrix) /.
                     ChargeConjugationMatrixInverse ->
                    (-ChargeConjugationMatrixInverse)
                      
       ), IsolateHead->c$CC, IsolateSplit->Infinity
       ] /; Length[Position[{yz}, Spinor]] < 3;
 c$CCfrh /: HoldForm[c$CCfrh[ii_]] := c$CC[ii];

 ComplexIndex[ComplexIndex[x_]] := x;

 cLIndex[x_, dime___] := LorentzIndex[ComplexIndex[x], dime];
 cSIndex[x_] := SUNIndex[ComplexIndex[x]];
(* CHANGE JAN. 1993 *)
 Unprotect[Conjugate];
 Conjugate[x_Pair] := (x /. {Polarization[k_,a_,in___] :> 
                             Polarization[k,Conjugate[a],in] }
                      ) /;!FreeQ[x, Polarization];
 Protect[Conjugate];

conpa[x__] := Pair[x] /. {Polarization[k_, a_, in___] :>
                            Polarization[k, Conjugate[a], in]};
(* ComplexConjugatedef *)
 sunfcomp[a___] := SUNF @@ ({a}/.ComplexIndex -> Identity);
 sundcomp[a___] := SUNDelta @@ ({a}/.ComplexIndex -> Identity);
 ComplexConjugate[b_HoldForm] := b /; 
              FreeQ2[FRH[b], {Dot,LorentzIndex,SUNIndex,Complex}];
 ComplexConjugate[x_ /; (Head[x] =!= HoldForm)] := 
                  compcon[x] /. SUNF -> sunfcomp /.
                         SUNDelta -> sundcomp/.
                          compcon -> compcon2 /. compcon2 ->  
                              ComplexConjugate;

 compcon2[x_/;!FreeQ[x, HoldForm]] := compcon[FRH[x]];
 compcon[x_^n_?NumberQ] := compcon[x]^n;
 compcon[x_Plus] := compcon /@ x;
 compcon[x_Times]:= compcon /@ x;

 compcon[b_HoldForm] := b /; 
              FreeQ2[FRH[b], {Dot,LorentzIndex,SUNIndex,Complex}];
 compcon[x_ /; (Head[x] =!= Plus) && (Head[x] =!= Times)] := 
             Block[{nx=x},
                  If[!FreeQ[nx, SUNF], nx = Expand2[nx, SUNF]];
                    nx = (nx /. Dot -> rev /. rev -> Dot);
                    nx = nx //. c$CC -> c$CCfrh;
                    nx = nx /. 
                          LorentzIndex -> cLIndex /.
                          SUNIndex  -> cSIndex /.  
                          Complex[a_, b_] -> Complex[a, -b];
                   nx] /; FreeQ[x, HoldForm];

(* CAREFUL: Complex[a_, b_] -> Complex[a, -b] is only true if no complex
   variables are in denominators!!!!, (which is the case in HEP, unless you
   have width in the propagators ...)
*)

(* this might cause trouble ... (commented out 15.5.93) 
 Polarization/:Momentum[Polarization[k_,i___],di_Symbol
                               ]:=Momentum[Polarization[k,i]];
 Polarization/:Momentum[Polarization[k__], di_Symbol - 4 ]:= 0;
*)
(* #################################################################### *)
(*                             Main14                                  *)
(* #################################################################### *)
(* MetricTensordef *)
 Options[MetricTensor]={Dimension->4};
MetricTensor[x__]:=(*MetricTensor[x]=*)metricTensor[x];
 loin1[x_,___]:=x;
 metricTensor[a_ b_,opt___]:=metricTensor[a,b,opt];
 metricTensor[a_^2 ,opt___]:=metricTensor[a,a,opt];
 metricTensor[ x__ ]:=(metricTensor@@({x}/.LorentzIndex->loin1));
 metricTensor[x_,x_,op_:{}]:=(Dimension/.op/.Options[MetricTensor]);
 metricTensor[ x_, y_,op_:{} ] :=
    Pair[ LorentzIndex[x,Dimension/.op/.Options[MetricTensor] ],
          LorentzIndex[y,Dimension/.op/.Options[MetricTensor] ]
        ];                               
(* PolarizationVectordef *)
(* CHANGE JAN. 1993 *)
Polarization[k_]:=Polarization[k]=Polarization[k, I];
PolarizationVector[x__]:=PolarizationVector[x]=polarizationVector[x];

(* By default a second argument "I" is put into Polarization *)
(* This is changed to "-I" for conjugate polarization vectors *)
 polarizationVector[k_Polarization,mu_]:=
       FourVector[k, mu, Dimension -> 4 ];
 polarizationVector[k_Polarization,mu_,glu_]:=
       FourVector[k, SUNIndex[glu/.SUNIndex->Identity], 
                  mu, Dimension->4 ];
 polarizationVector[k_,mu_]:=
       FourVector[Polarization[k, I], mu, Dimension->4 ];

 polarizationVector[k_,mu_,glu_]:=
    If[FreeQ[glu, Blank],
       FourVector[Polarization[k, I,SUNIndex[glu/.SUNIndex->Identity]], 
                  mu, Dimension->4 ],
       FourVector[Polarization[k, I, glu], mu, Dimension -> 4]
      ];
(* #################################################################### *)
(*                             Main15                                  *)
(* #################################################################### *)
(*  DiracMatrixdef   *)
 Options[ DiracMatrix ] = {Dimension->4};
 Options[ DiracSlash  ] = {Dimension->4};

DiracMatrix[x__]:=DiracMatrix[x]=diracMatrix[x];

 If[Global`$BreitMaison === True,
    DiracMatrix[6] = 1/2 + 1/2 DiracGamma[5];
    ChiralityProjector[1] = 1/2 + 1/2 DiracGamma[5];
    DiracMatrix[7] = 1/2 - 1/2 DiracGamma[5];
    ChiralityProjector[-1] = 1/2 - 1/2 DiracGamma[5]
  ];

 diracMatrix[n_?NumberQ y_]:=n diracMatrix[y];
 DiracMatrix[n_?NumberQ y_,{}]:=n diracMatrix[y];
 diracMatrix[n_?NumberQ y_,opt_]:=n diracMatrix[y,opt];
 diracMatrix[x_,y_]:=diracMatrix[x].diracMatrix[y]/;(FreeQ[y,Rule]&&y=!={});
 diracMatrix[x_,y__,{}]:= diracMatrix[Dot[x,y]];
 diracMatrix[x_,y__,z_]:= diracMatrix[Dot[x,y],z]/;!FreeQ[z,Rule];
 diracMatrix[x_,y__,z_]:= diracMatrix[Dot[x,y,z]]/; FreeQ[z,Rule];
 diracMatrix[x_ y_Plus,opt_:{}]:= diracMatrix[Expand[x y],opt];
 diracMatrix[x_Plus,opt_:{}]:= diracMatrix[#,opt]& /@ x;
 diracMatrix[x_Dot,opt_:{}] :=  diracMatrix[#,opt]& /@ x;
 diracMatrix[n_Integer,___]:=DiracGamma[n];
 diracMatrix[5,opt_:{}]:=DiracGamma[5];
 diracMatrix[6,opt_:{}]:=DiracGamma[6];
 diracMatrix[7,opt_:{}]:=DiracGamma[7];
 diracMatrix["+"]=DiracGamma[6];  
 diracMatrix["-"]=DiracGamma[7];
If[Global`$BreitMaison=!=True,
   ChiralityProjector[1]=DiracGamma[6];
   ChiralityProjector[-1]=DiracGamma[7]
 ];
 diracMatrix[x_,op_:{}] := DiracGamma[LorentzIndex[ x,
            (Dimension/.op/.Options[DiracMatrix])  ] ,
               (Dimension/.op/.Options[DiracMatrix])
                                     ]/;(Head[x]=!=Dot && !IntegerQ[x]);
(*  DiracSlashdef   *)

DiracSlash[x__]:=DiracSlash[x]=diracSlash[x];
 
 ndot[]=1;
 ndot[a___,ndot[b__],c___]:=ndot[a,b,c];
 ndot[a___,b_Integer,c___]:=b ndot[a,c];
 ndot[a___,b_Integer x_,c___]:=b ndot[a,x,c];
 diracSlash[x_,y_]:=diracSlash[ndot[x,y]]/;(FreeQ[y,Rule]&&y=!={});
 diracSlash[x_,y__,{}]:=diracSlash[ndot[x,y]];
 diracSlash[x_,y__,z_]:=diracSlash[ndot[x,y],z]/;!FreeQ[z,Rule];
 diracSlash[x_,y__,z_]:=diracSlash[ndot[x,y,z]]/;FreeQ[z,Rule];
 diracSlash[x__]:= (diracSlash@@({x}/.Dot->ndot) )/;!FreeQ[{x},Dot];
 diracSlash[n_Integer x_ndot,opt_:{}]:=n diracSlash[x,opt];
 diracSlash[x_ndot,opt_:{}] := Expand[(diracSlash[#,opt]& /@ x)
                                     ]/.ndot->Dot;
(*   pull out a common numerical factor *)
 diracSlash[x_,op_:{}] := Block[{dtemp,dix,eins,numf,resd},
          dix = Factor2[ eins Expand[x]];
          numf = NumericalFactor[dix];
          resd = numf DiracGamma[ Momentum[Cancel[(dix/.eins->1)/numf],
            (Dimension/.op/.Options[DiracSlash])  ] ,
              (Dimension/.op/.Options[DiracSlash])
                                ]
                               ]/;((Head[x]=!=Dot)&&(Head[x]=!=ndot));

(* #################################################################### *)
(*                             Main16                                  *)
(* #################################################################### *)
(* FourVectordef, FourVectorDdef *)
Options[FourVector]={Dimension->4};

FourVector[x__]:=FourVector[x]=fourVector[x];

 fourVector[Literal[Dot[x_,y_]],op___]:=fourVector[x,y,op];
 fourVector[ x_Momentum,y___]:= fourVector[x[[1]],y];
 fourVector[ x_,y_LorentzIndex,op___]:= fourVector[x,y[[1]],op];
(*   pull out a common numerical factor *)
 fourVector[ x_,y_,opt_:{}]:=Block[{nx,numfa,one,result},
                                    nx = Factor2[one x];
                                    numfa = NumericalFactor[nx];
       result = numfa Pair[ LorentzIndex[y, Dimension/.opt/.
                                           Options[FourVector]],
                             Momentum[Cancel[nx/numfa]/.one->1,
                                      Dimension/.opt/.Options[FourVector]]
                          ]; result]/;!FreeQ[x,Plus];
 
 fourVector[ x_, y_,opt_:{} ] := Pair[
     LorentzIndex[y,Dimension/.opt/.Options[FourVector]],
     Momentum[x,Dimension/.opt/.Options[FourVector]] ]/;FreeQ[x,Plus];

(* ScalarProductdef *)
 Options[ScalarProduct]={Dimension->4};
 ScalarProduct/:Set[ScalarProduct[a__],z_]:= Block[{ste},
  ste = ScalarProduct[a];
          Set@@{ste/NumericalFactor[ste], z / NumericalFactor[ste]};
                                                  ];
 ScalarProduct[Literal[Dot[x_,y_]],op___]:=ScalarProduct[x,y,op];
 ScalarProduct[ x_, y_,opt_:{} ] := Pair[
     Momentum[x,Dimension/.opt/.Options[ScalarProduct]],
     Momentum[y,Dimension/.opt/.Options[ScalarProduct]] 
                                        ]/;FreeQ[{x,y}, Momentum];
 ScalarProduct[ x_, y_,opt_:{} ] := Pair[ x, y ]/;!FreeQ[{x,y}, Momentum];

(* ToFourDimensionsdef *)
 ToFourDimensions[x_]:=x/.Momentum[v_,_Symbol]->Momentum[v]/.
                          LorentzIndex[w_,_Symbol]->LorentzIndex[w]/.
                          Momentum[v_,4 - _Symbol]->0/.
                          LorentzIndex[w_,4 - _Symbol]->0;


(* #################################################################### *)
(*                             Main17                                  *)
(* #################################################################### *)

(*MomentumExpanddef*)
 MomentumExpand[x_]:=x/;FreeQ[x,Momentum]; 
   fourvecevlin[n_?NumberQ z_, dime___]  := n Momentum[z,dime];
   fourvecev[y_,di___] := memset[ fourvecev[y,di],
            ReleaseHold[Distribute[fourvecevlin[ 
                                  Expand2[relh[y],Momentum],Hold[di]]
                                  ]/.fourvecevlin->Momentum ]];
 MomentumExpand[x_] := x/.Momentum->fourvecev;

(* #################################################################### *)
(*                             Main18                                  *)
(* #################################################################### *)

(* ************************************************************** *)
 Pair[0,_]:=0;            (* don't reverse the order w.r.t.       *)
 Pair[Momentum[0],_]:=0;  (* SetAttributes!!!   (Bug in Math. ...)*)
(* ************************************************************** *)
 SetAttributes[{Pair,sCO,sceins,scev,sce,scevdoit,sczwei},Orderless];
(* ************************************************************** *)
 Pair[n_Integer x_,y_] := n Pair[x, y];                   (*Pairdef*)
 Pair[n_ x_Momentum, y_] := n Pair[x, y];
 Pair[y_, n_,x_Momentum] := n Pair[y, x];
 Pair[n_ x_LorentzIndex, y_] := n Pair[x, y];
 Pair[y_, n_,x_LorentzIndex] := n Pair[y, x];
 Pair[Momentum[ n_Integer x_,di___],y_] :=  n Pair[Momentum[x,di],y];
 Pair[Momentum[x_,___],Momentum[Polarization[x_, ___],___]]:=0;
 Pair[Momentum[x_,___],Momentum[Polarization[n_?NumberQ x_, ___],___]]:=0;

 Pair[Momentum[pi_,___],Momentum[Polarization[x_Plus, ki___],dii___]]:=
  scev[Momentum[x+pi,dii], Momentum[Polarization[x, ki],dii]]/;
        ( pi + Last[x] )===0;
 Pair[Momentum[pi_,___],Momentum[Polarization[x_Plus, ki___], dii___]]:=
  scev[Momentum[pi-x,dii], Momentum[Polarization[x, ki],dii]]/;
        ( pi - Last[x] )===0;

(* this is a convention ... *)
 Pair[Momentum[Polarization[x_,__],___], Momentum[Polarization[x_,__],___]
     ] = -1;

(* #################################################################### *)
(*                             Main19                                  *)
(* #################################################################### *)
                                     (*ExpandScalarProductdef*)
(*
 ExpandScalarProduct[x_] := Collect2[ MomentumExpand[x], Momentum ]/;
                            FreeQ[x, Pair] && (!FreeQ[x, Momentum]);
*)
(* this is essential for OneLoop !!!  (PropagatorDenominator)*)
 ExpandScalarProduct[x_] := If[LeafCount[x]<1000,
     Expand[ FixedPoint[pairexpand1,x, 3]//MomentumExpand],
             FixedPoint[pairexpand1,x, 3]//MomentumExpand
                              ];
(*
 pairexpand1[x_]:= Expand[  x/.Pair->scevdoit, Pair ];
*)
 pairexpand1[x_]:=  x/.Pair->scevdoit;
 pairexpand[x_] :=  x/.Pair->scev;  (*pairexpanddef*)
   scev[x_,y_]:= memset[ scev[x,y], scevdoit[x,y] ]; 
   scev[x_,y_]:= scevdoit[x,y]; 
   scevdoit[x_,y_] := Distribute[ sceins@@
                                 ( Expand[ MomentumExpand/@{x,y} ] )
                      ]/.sceins->sczwei/.sczwei->sCO/.sCO->Pair;
(* #################################################################### *)
(*                             Main20                                 *)
(* #################################################################### *)
(* ******************************************************************* *)
(* sCO has the actual contracting properties *)                (*sCOdef*)
(* ******************************************************************* *)
   sCO[ LorentzIndex[a_,di___], epsmu_ LorentzIndex[mu_, dimen___] ]:=
   ( epsmu /. LorentzIndex[mu,dimen]->LorentzIndex[a,di] ) /; 
   !FreeQ2[epsmu, {Eps, LorentzIndex[mu, dimen]}];
   
   sCO[ Momentum[x_,___],Momentum[Polarization[x_,___]]]:=0;
   sCO[ Momentum[x_,___],Momentum[Polarization[n_?NumberQ x_,___]]]:=0;
   sCO[Momentum[pi_,___],Momentum[Polarization[x_Plus, ki___]]]:=
    scev[Momentum[x+pi], Momentum[Polarization[x, ki]]]/;
        ( pi + Last[x] )===0;
   sCO[Momentum[pi_,___],Momentum[Polarization[x_Plus, ki___]]]:=
    scev[Momentum[pi-x], Momentum[Polarization[x, ki]]]/;
                ( pi - Last[x] )===0;

   sCO[ LorentzIndex[x_], LorentzIndex[x_] ]  := 4;
(*new ...*)
   sCO[ LorentzIndex[x_], LorentzIndex[x_,_Symbol] ]  := 4;

   sCO[ LorentzIndex[x_,di_], LorentzIndex[x_,di_] ] := di;
   sCO/: Literal[ sCO[LorentzIndex[z_,___],x_] ]^2 := sCO[x,x];
(*
   sCO/: Literal[ sCO[LorentzIndex[z_,___],x_] f_[a__][b___] ] :=
   (f[a][b]/.LorentzIndex[z,___]->x)/;
    !FreeQ[f[a][b]//Hold,LorentzIndex[z,___]];
*)
   sCO/: Literal[ sCO[LorentzIndex[z_,___],x_] f_[a__] ] :=
   (f[a]/.LorentzIndex[z,___]->x)/;
    !FreeQ[f[a],LorentzIndex[z,___]];

   sCO/: Literal[ Dot[A___, sCO[LorentzIndex[z_,___],x_],B___,
                  m_. f_[a__], c___ ] ] := 
   Dot[A,B,(m f[a]/.LorentzIndex[z,___]->x),c]/;
       !FreeQ[f[a], LorentzIndex[z,___]];

   sCO/: Literal[ Dot[A___, m_. f_[a__],B___, sCO[LorentzIndex[z_,___],x_],
                  c___ ] ] := 
   Dot[A.(m f[a]/.LorentzIndex[z,___]->x),B,c]/;
      !FreeQ[f[a]//Hold,LorentzIndex[z,___]];
(* #################################################################### *)
(*                             Main21                                 *)
(* #################################################################### *)
(* ******************************************************************** *)
(* definitions for dimension = D-4                                      *)
(* ******************************************************************** *)
   sCO[ _[_,_Symbol-4],_[_] ]:=0;
   sCO[ v_[x_,di_Symbol-4],w_[y_,di_Symbol] ] := sCO[v[x,di-4],w[y,di-4] ];
   sCO[ w_[y_,di_Symbol],v_[x_] ] := sCO[ v[x], w[y] ];
   sCO[ v_[x_], w_[y_,di_Symbol] ] := sCO[ v[x], w[y] ];
   sceins[0,_]:=0;                               (*sceinsdef*)
   sceins[a_LorentzIndex b_, c_] := b sceins[a, c];
   sceins[a_Momentum b_, c_] := b sceins[a, c];
   sczwei[ _[_],_[_,_Symbol-4] ]:=0;             (*sczweidef*)
   sczwei[ v_[x_,di_Symbol-4],w_[y_,di_Symbol] ]:=
			      sczwei[v[x, di-4], w[y, di-4]];
   sczwei[ w_[y_,di_Symbol],v_[x_] ]:=sczwei[ v[x],w[y] ];
   sce[x_,y_] := memset[sce[x, y],      (*scedef*)
                  Distribute[ sceins@@( Expand[ MomentumExpand/@{x,y} ])
                            ]/.sceins->sczwei/.sczwei->Pair
                       ];
   sCO[x_,y_] := memset[ sCO[x,y], 
                        Block[{sCOt=sce[x,y]},
                         If[ FreeQ[ sCOt, Pair ] ||
                              (Head[sCOt]=!=Plus)
                             , sCOt,Pair[x,y] 
                           ] ] ]/;FreeQ2[{x,y},{LorentzIndex}];
(* #################################################################### *)
(*                             Main22                                 *)
(* #################################################################### *)
(* LeviCivitadef*)
Options[LeviCivita] = { Dimension -> 4 };
 LeviCivita[ a__, ops___Rule ]:= EpsEvaluate[ Eps@@( LorentzIndex/@{a} ) 
				] /; frlivc[a] && FreeQ[{a},Rule] &&
         ( (Dimension /. {ops} /. Options[LeviCivita]) === 4);


 LeviCivita[ a__, ops___Rule ]:= EpsEvaluate[ Eps@@( 
   LorentzIndex[#, Dimension/.{ops}/.Options[LeviCivita]]& /@{a} ) 
				] /; frlivc[a] && FreeQ[{a},Rule] &&
         ( (Dimension /. {ops} /. Options[LeviCivita]) =!= 4);

(* ******************************************************************* *)
                                                         (*Epsdef*)
 Options[Eps] = {Dimension -> 4 };
 Eps[a___, lv1_[mu_,___], b___, lv1_[mu_,___],c___ ] := 0;
 Eps[x__] :=  0 /; !FreeQ[{x}, lv_[_,_Symbol -4]];
(*
 Eps[x__, ops___Rule] := 
       ( Eps@@( {x}//ToFourDimensions ) )/;
       (!FreeQ[{x},Momentum[_,_]] || !FreeQ[{x},LorentzIndex[_,_]]) &&
       ( (Dimension /. {ops} /. Options[Eps]) === 4);
*)

(* #################################################################### *)
(*                             Main23                                 *)
(* #################################################################### *)
 EpsEvaluate[x_] := x /; FreeQ[x,Eps];     (*EpsEvaluatedef*)
 EpsEvaluate[x_] := x//.Eps->epsev;        (*epscondef*)
   epscon/: epscon[a1_,a2_,a3_,a4_]^n_Integer?Positive :=  (   (
            ( -Det[{{sCO[a1,a1],sCO[a1,a2],sCO[a1,a3],sCO[a1,a4]},
                    {sCO[a2,a1],sCO[a2,a2],sCO[a2,a3],sCO[a2,a4]},
                    {sCO[a3,a1],sCO[a3,a2],sCO[a3,a3],sCO[a3,a4]},
                    {sCO[a4,a1],sCO[a4,a2],sCO[a4,a3],sCO[a4,a4]}}
                  ]//Expand 
            )/.sCO->Pair ) epscon[a1,a2,a3,a4]^(n-2) );
   epscon/: epscon[a1_,a2_,a3_,a4_] epscon[b1_,b2_,b3_,b4_] :=
            ( -Det[{{sCO[a1,b1],sCO[a1,b2],sCO[a1,b3],sCO[a1,b4]},
                    {sCO[a2,b1],sCO[a2,b2],sCO[a2,b3],sCO[a2,b4]},
                    {sCO[a3,b1],sCO[a3,b2],sCO[a3,b3],sCO[a3,b4]},
                    {sCO[a4,b1],sCO[a4,b2],sCO[a4,b3],sCO[a4,b4]}}
                  ]//Expand
            )/.sCO->Pair;                              (*epsevdef*)
   epsev[A__] := ( Expand /@ (Distribute[Dot[A]]//MomentumExpand) )/.
                 Dot->epsevlin/.epsevlin->epsevantilin;
   epsevlin[a___,b_ c_Momentum,d___] := b epsevlin[a,c,d];
   epsevlin[a___,b_ c_LorentzIndex,d___] := b epsevlin[a,c,d];
   epsevantilin[a__] := Signature[{a}] Eps@@Sort[{a}];
(* ************************************************************** *)
(* #################################################################### *)
(*                             Main24                                   *)
(* #################################################################### *)

(* #################################################################### *)
(*                             Main24A (07/1993)                        *)
(* #################################################################### *)
   Contract3[x_Plus] := Map[Contract3, x];
   Contract3[x_ /; Head[x] =!= Times] := Contract[x];
   Contract3[x_Times /; !FreeQ[x, LorentzIndex]] := 
    Block[{nx = x, nonli, lipa, nec = 0, ic,epli},
      nx = Contract[x, Expanding -> False];
        If[Head[nx] =!= Times, nec = Contract[nx],
           nonli = Select[nx, FreeQ[#, LorentzIndex]&];
           lipa  = Select[nx,!FreeQ[#, LorentzIndex]&];
           If[Head[lipa] =!= Times, epli = 1,
              epli  = Select[lipa, !FreeQ2[#,{Eps, DiracGamma}]&];
              lipa = lipa / epli;
             ];
           If[Head[lipa] =!= Times, 
              If[Head[lipa] === Plus,
                 nec = Contract3[lipa epli],
                 nec =  Contract[lipa epli]
                ],
              If[Length[lipa] < 2, nec = Contract[lipa epli],
                 nec = lipa[[1]] epli;
                 For[ic = 2, ic <= Length[lipa], ic++,
                     print2["ic = ", ic];
                     nec = Contract[nec, lipa[[ic]]];
                     nec = Collect2[nec, LorentzIndex, Factoring -> False];
                    ]; 
                ];
             ];
               nec = nec nonli;
          ];
            nec];
               
(* #################################################################### *)
(*                             Main24a (01/1993)                        *)
(* #################################################################### *)


Options[Contract2] = {Collecting -> False};

(* bb is assumed to be collected w.r.t. to LorentzIndex !!! *)
Contract2[a_, bb_, ops___Rule]:= Block[
        {sel, ct, rc, lco, lct, lastct, nop, b = bb, col, conT},
  col = Collecting /. {ops} /. Options[Contract2];
  If[Head[a] =!= Times, rc = Contract[a, b],
     lco[x_,y_] := If[Length[x]>Length[y], True, False];
     sel = Select[a, FreeQ[#, LorentzIndex]&];
     ct  = a/sel;
     nop = Select[ct, FreeQ[#, Plus]&];
     ct = ct/nop;
     If[Head[ct] =!= Times, rc = sel Contract[ct nop, b],
        ct = Sort[List @@ ct, lco];
        If[ nop =!= 1,
            lastct = contract21[b, nop],
            lastct = b nop
          ];
        lct = Length[ct];
        If[lct === 1, rc = sel contractLColl[ct[[1]], lastct] ];
        If[lct > 1,
           rc = sel Contract[Times @@ Take[ct, lct-1],
                             ct[[lct]], lastct ]
          ];
       ];
    ];
print2["lct = ",lct];
If[!FreeQ[rc, LorentzIndex],
   rc = Contract[rc, Expanding -> False];
  ];
If[!FreeQ[rc, LorentzIndex],
   print1["contracting agagin at the end of Contract2 "];
   rc = Contract[rc]
  ];
      rc];

Contract2[a_] := Block[{sel, ct, rc, lco, lct, lastct, nop},
  If[Head[a] =!= Times, rc = Contract[a],
     lco[x_,y_] := If[Length[x]>Length[y], True, False];
     sel = Select[a, FreeQ[#, LorentzIndex]&];
     ct  = a/sel;
     nop = Select[ct, FreeQ[#, Plus]&];
     ct = ct/nop;
     If[Head[ct] =!= Times, rc = sel Contract[ct nop],
        ct = Sort[List @@ ct, lco];
        If[ nop =!= 1, 
            lastct = contract21[Last[ct], nop],
            lastct = Last[ct] nop;
          ];
        lct = Length[ct];
        If[lct === 2, rc = sel contractLColl[ct[[1]], lastct] ];
        If[lct > 2, 
           rc = sel Contract[Times @@ Take[ct, lct-2], 
                             ct[[lct-1]], lastct ] 
          ];
       ];
    ];
print2["lct = ",lct];

If[!FreeQ[rc, LorentzIndex], 
   print1["contracting agagin at the end of Contract2 "];
   rc = Contract[rc] 
  ];
      rc];

Contract[a_, b_ /;Head[b] =!= Rule, c_ /; Head[c] =!= Rule, ops___Rule] := 
   Block[{lc, new = 0, i},
 print2["longcontract1"]; 
 lc = Contract[b, c, ops];
 print2["longcontract1done"]; 
   new = Contract[lc, a, ops];

(*
 If[Head[lc] =!= Plus,
    new = contractLColl[a, lc],
    For[i = 1, i <= Length[lc], i++,
        print2["X iC = ", i,"  out of ",Length[lc]];
        new = new + (Expand[ncc=Contract[lc[[i]], a] ] /. Pair->pair/.
                     pair2 -> Pair);
        print2["length of new = ",Length[new]];
       ] ];
*)

 new];

Contract[x_, y_ /; FreeQ2[y, {LorentzIndex,Eps}]] := Contract[x] y;
Contract[x_ /; FreeQ[x, DiracGamma],
         y_ /; !FreeQ[y, DiracGamma]] := Contract[y,x];

(*
Contract[a_, b_Times] :=  Contract[a b];
*)

Contract[a_, b_Times] := Block[{bb},
  If[MatchQ[b, Literal[Times__Pair]], contract21[ a, b ],
     bb = Collect2[b, LorentzIndex, Factoring -> False];
     If[Head[bb] === Plus,
        contractLColl[a, bb], 
        contract21[a, bb]
       ]
    ]                         ];

Contract[a_, b_ /; ((Head[b]=!=Times) && (Head[b] =!= Plus) && (Head[b] =!= Rule))
        ] := Contract[ a b ] /; (Head[a] =!= Plus);


Contract[a_, b_Plus, ops___Rule] := 
  If[(Collecting /. {ops} /. Options[Contract]) === True,
     contractLColl[a, Collect2[b, LorentzIndex]],
     contractLColl[a, b]
    ];

(* contractLColldef *)
contractLColl[a_, b_ /; Head[b] =!= Plus] := contract21[a, b];
contractLColl[lon_, shor_Plus] := Block[{neew = {}, long = lon,short = shor,tet},
    If[$VeryVerbose > 0, 
        WriteString["stdout","Long contraction ", Length[long], " * ",
                     Length[short], " \n "]
       ];
For[ij = 1, ij <= Length[short], ij++,
   If[$VeryVerbose > 2,
      WriteString["stdout"," | ", ij, "  XXX | "]
     ];
       tet = contract21[long, short[[ij]] ];

   If[!FreeQ[tet, LorentzIndex], 
      tet = tet /. Pair->pair /. pair2 -> Pair];
   If[!FreeQ[tet, LorentzIndex], 
      If[$VeryVerbose > 1, WriteString["stdout","expanding in contractLColl "]i];
      tet = Expand2[tet, LorentzIndex] /. Pair->pair /. pair2 -> Pair;
     ];
      If[Head[tet] === Plus,
         neew  = Join[neew, Apply[List, tet]],
         AppendTo[neew, tet]
        ];
   ];
                     Apply[Plus, neew]];

(* local easy contraction rules *) (* paird *)
fdi[]=4;
fdi[_Symbol] := 4;
fdi[xx_Symbol, xx_Symbol] := xx;
fdi[4, _Symbol] := 4;
fdi[_Symbol, 4] := 4;

SetAttributes[{pair, pair2}, Orderless];
pair[a_, b_] := pair[a, b] = pair2[a, b];
pair2[LorentzIndex[a_, di1___], LorentzIndex[a_, di2___]] := fdi[di1, di2];
pair2/: pair2[LorentzIndex[a_, dim1___], LorentzIndex[b_, dim2___]]^2 := 
        fdi[dim1, dim2];

pair2/: pair2[LorentzIndex[a_,de1___], Momentum[b_, de2___]]^2 := 
        Pair[Momentum[b, fdi[de1,de2]], Momentum[b,fdi[de1,de2]]];
pair2/: pair2[LorentzIndex[al_,di___], z_] pair2[LorentzIndex[al_,di2___], w_] := pair2[z, w];

(* contract21 can still have products in the first argument *)
(* by construction the second argument will always be either a 
   product or just Pair[  ,  ]
*)

contra3a[xx_, {pr_, prl__}] :=
      contra3a[contra3a[xx, {pr}], {prl}];

contra3b[xx_, {alien_ /; Head[alien] =!= pair2} ] := Expand2[xx alien, Pair];

contra3c[xx_, {Pair[LorentzIndex[mu_,di___], alpha_]} ] :=Block[{nxx},
    If[FreeQ[xx, LorentzIndex[mu,___]],
       nxx = Expand2[xx Pair[LorentzIndex[mu, di], alpha], Pair],
       If[Head[xx]===Plus, nxx = Apply[List, xx], nxx = {xx}];
print2["contra3c : Length of xx now ", Length[nxx]];
       nxx = nxx /. LorentzIndex[mu, ___] -> alpha;
       nxx = Apply[Plus, nxx];
      ];
         nxx];

contract21[xx_Plus, yy_] := contract22[xx, yy /. Pair -> pair] /. 
      contra4 -> contra3a /.  contra3a -> contra3b /.
      Pair -> pair /.  contra3b -> contra3c /. pair2 -> Pair;

list2[x_] := If[Head[x] === Times, List @@ x, {x}];
contract22[xx_, 0] := 0; 
contract22[xx_, yy_Pair] := contra3a[xx, {yy}] /.  contra3a -> contra3b /.
      Pair -> pair /.  contra3b -> contra3c /. pair2 -> Pair;
contract22[xx_, yy_Times]:= ( (yy/#) contra4[xx, list2[#]] )&[
                             Select[yy, !FreeQ[#, LorentzIndex]&]];


(*
contract21[xx_Plus, yy_] :=(iCcount=1; Apply[Plus,
                   Table[ (xx[[ii]] yy) /. Pair -> pair /. 
                                               pair2 -> Pair
                              , {ii, 1, Length[xx]}] ]);
*)
                                           
contract21[xx_ /;(Head[xx] =!= Plus) && (Head[xx] =!= Times), yy_] :=
  Contract[xx yy,Expanding -> False];

contract21[xxx_Times, yyy_] := ( (xxx/#) contit[#, yyy] )&[
                             Select[xxx, !FreeQ[#, LorentzIndex]&] ];
contit[xx_ , yy_] := 
  If[FreeQ[xx, LorentzIndex], 
     xx Contract[yy],
     If[Head[xx] =!= Times, Contract[xx yy],
        If[Length[xx] =!= 2, Contract[xx yy],
           If[(Head[xx[[1]]] === Plus) && (Head[xx[[2]]] === Plus),
iCcount = 1;
print2["contracting a product of a ",Length[xx[[1]]], " term sum  by a",
       Length[xx[[2]]], " term sum"];
(* that's the common situation !! *)
              Apply[ Plus, Flatten[Table[ (xx[[1, ii]] xx[[2, jj]] yy
                                               ) /. Pair -> pair /. 
                                                    pair2  -> Pair
                                              , {ii,1,Length[xx[[1]]]},
                                                 {jj,1,Length[xx[[2]]]}
                   ]              ]     ],
              Contract[xx yy]
             ] ] ] 
    ];


(* #################################################################### *)

(* coneinsdef    *)
   coneins[ x_ ]  := memset[coneins[x], x/.Pair->sCO/.sCO->Pair ]; 
(* contractlidef *)
   contractli[x_] := memset[contractli[x],x] /; FreeQ[x//Hold,LorentzIndex];
   contractli[x_] := Contract[ x, Expanding->True, Factoring->False,
                      EpsContract->False ];
   conall[ x_ ] := Contract[ x,                               (*conalldef*)
                   Expanding->True, EpsContract->True, Factoring->False ];
                                      (*Contractdef*)
(* CHANGE W.R.T. to MANUAL !!! *)
 Options[Contract] = { Collecting -> True,
                       Expanding->True, EpsContract->True, Factoring->False,
                        Schouten -> 0 };
 Contract[x_,opt___Rule] := x /; FreeQ2[ x,{LorentzIndex,Eps,Momentum} ];
 Contract[x_,opt___Rule] := Block[{ contractres,liplu,epscontractopt,
         lip,contractopt = Join[{opt},Options[Contract]]//Flatten,
           schout },
   contractexpandopt   = Expanding/.contractopt;
   contractepsopt      = EpsContract/.contractopt;
   contractfactoring   = Factoring/.contractopt;
   contractres = x/.Pair->sCO/.sCO->sceins/.sceins->sCO;
   If[ contractexpandopt === True,
       contractres = contractres /. 
                     {(yy_ /;(!FreeQ[yy, LorentzIndex] && 
                               FreeQ[yy,Eps] && !FreeQ[yy,Plus])
                      )^2 :>
                      ((Contract @@ {yy/.sCO->Pair, yy/.sCO->Pair} 
                      ) /. Pair -> sCO)
                     };
     ];
   schout = Schouten /. contractopt;
        If[ contractexpandopt === True,
            liplu[y__]:=lip[y]/;FreeQ[{y},LorentzIndex];
            liplu[y__]:=Plus[y]/;!FreeQ[{y},LorentzIndex];
            contractres = Expand[contractres/.Plus->liplu]/.lip->Plus
          ];
        If[ contractepsopt === True,
            If[ !FreeQ[contractres, Eps],
                contractres = Expand2[ contractres//EpsEvaluate, Eps];
                contractres = contractres/.Eps->epscon/.epscon->epsev
              ],
            contractres = contractres//EpsEvaluate
          ];
        If[ contractexpandopt=== True,
            contractres = Expand2[contractres, LorentzIndex] ];
        If[ (contractexpandopt===True) && (!FreeQ[contractres, Eps]) &&
            (contractepsopt===True),
            contractres = Expand2[ contractres, Eps ] ];
        If[ !FreeQ[ contractres,Eps ],
            contractres = contractres//EpsEvaluate//EpsEvaluate
          ];
        contractres = contractres/.Pair->sCO/.sCO->Pair;
       If[(contractepsopt===True) && 
          (!FreeQ2[contractres, {Eps,DiracGamma}]),
           contractres = doubleindex[Schouten[contractres, schout]]/.
                                     doubleindex->double2/.
                                     eepp->Eps/.double3->ident3
         ];
        If[ contractfactoring=== True,
            contractres = Factor2[ contractres ]
          ];
  contractres                 ](* EndContract *);

(* #################################################################### *)
(*                             Main25                                 *)
(* #################################################################### *)
   (*  doubleindexdef *)
  (* For canonizing dummy indices between Eps and gammas *)
  doubleindex[x_ ]:=x /;  FreeQ[x, Eps] || FreeQ[x, DiracGamma];
  doubleindex[x_Plus]:=doubleindex /@ x;
  double2[x_] := If[EvenQ[Length[Position[x, $MU]]],
                    double3a[x/.Eps->eepp/.$MU->Unique[DUmmY], 1] /. 
                     double3a-> double3,
                    double3a[x/.Eps->eepp, 1] /.  double3a-> double3
                   ];

  double3a[x_, i_] := double3a[x, i+1] /; !FreeQ[x, $MU[i]];

  double3[ m_. eepp[a1___, LorentzIndex[be_, di___], a2___], j_ ] := 
          (m/.be->$MU[j]) Eps[a1,LorentzIndex[$MU[j],di],a2]/; 
         (!FreeQ[m, LorentzIndex[be, di]]) && 
           FreeQ2[m, Select[{a1,a2}, Head[#]===LorentzIndex&]];

  double3[ m_. eepp[a1___, LorentzIndex[mu1_, di1___], a2___, 
                           LorentzIndex[mu2_, di2___],a3___], j_ ] :=
          (m/.mu1->$MU[j]/.mu2->$MU[j+1]) * 
           Eps[a1,LorentzIndex[$MU[j],di1],a2, 
                  LorentzIndex[$MU[j+1],di2],a3]/;
         (!FreeQ[m, LorentzIndex[mu1, di1]]) &&
         (!FreeQ[m, LorentzIndex[mu2, di2]]) &&
         (FreeQ2[m, Select[{a1,a2,a3}, Head[#]===LorentzIndex&]]);

  double3[ m_. eepp[a1___, LorentzIndex[mu1_, di1___], a2___,
                           LorentzIndex[mu2_, di2___], a3___,
                           LorentzIndex[mu3_, di3___], a4___ ], j_ ]:=
        (  (m/.mu1->$MU[j]/.mu2->$MU[j+1]/.mu3->$MU[j+2]) *
           Eps[a1,LorentzIndex[$MU[j],di1], a2, 
                  LorentzIndex[$MU[j+1],di2],a3,
                  LorentzIndex[$MU[j+2],di3], a4]
        )/; (!FreeQ2[m, LorentzIndex[mu1,di1]] && 
             !FreeQ2[m, LorentzIndex[mu2,di2]] &&
             !FreeQ2[m, LorentzIndex[mu3,di3]] 
            ) &&
           FreeQ2[m, Select[{a1,a2,a3,a4}, Head[#]===LorentzIndex&]];
  
  double3[ m_. eepp[LorentzIndexx[mu1_,di1___],LorentzIndex[mu2_,di2___],
                    LorentzIndex[mu3_,di3___],LorentzIndex[mu4_,di4___]], _]:=
       (  (m/.mu1->$MU[1]/.mu2->$MU[2]/.mu3->$MU[3]/.mu4->$MU[4]) *
          Eps[LorentzIndex[$MU[1],di1],  LorentzIndex[$MU[2],di2], 
              LorentzIndex[$MU[3],di3], LorentzIndex[$MU[4],di4] ]
       ) /; (!FreeQ2[m, {LorentzIndex[mu1,di1],LorentzIndex[mu2,di2],
                         LorentzIndex[mu3,di3],LorentzIndex[mu4,di4]}]
            );

(* #################################################################### *)
(*                             Main26                                   *)
(* #################################################################### *)
            
   (* Schoutendef *)
  Schouten[y_, 0] := y;
  Schouten[y_, oparg_:42] := FixedPoint[ schouten[#, oparg]&, y, 14];

  schouten[x_,optarg_:42]:=memset[schouten[x, optarg],
                      Block[{i,nx,temp0,temp,lind,liget,lisch,ltemp,ntemp,
                       schou,sor, all,result,epc, epsnterms,numberofli},
                      epc[a_, b_] := If[ Length[Position[ 
                                     PartitHead[a, Eps][[2]], LorentzIndex]]<
                                     Length[Position[  
                                     PartitHead[b, Eps][[2]], LorentzIndex]],
                                        True, False, False];
                      epsnterms[0] = 0;
                      epsnterms[a_] := Block[{tem},
                                              tem = Collect2[a, Eps,
                                                     Factoring ->False];
                                              If[Head[tem]===Plus,
                                                 Length[tem], 1]
                                            ];  
                                          
                        
                      nx = EpsEvaluate[ExpandScalarProduct[x]//Expand]//Expand;
                         liget[a_. Eps[x1_[y1_], x2_[y2_], x3_[y3_], 
                                       x4_[y4_]] * Pair[x5_[y5_], x6_[y6_]]
                              ] := {x1[y1],x2[y2],x3[y3],x4[y4],x5[y5],x6[y6]};
                         lisch[{a1_,a2_,a3_,a4_,a5_,a6_}]:=
                         Pair[a5,a6] Eps[a1,a2,a3,a4];
(* Split the sum into two parts *)
                      result = nx;
                         all  = PartitHead[nx, Eps];
                      If[ Head[all[[2]]]===Plus,
                          temp0 = PartitHead[all[[2]], Pair];
                          temp = temp0[[2]];
                      If[(Head[temp]===Plus) && Length[temp] > 1 && 
                         Length[temp] < optarg,
                         ltemp = Length[temp];
                         numberofli = Length[Position[temp[[1]],
                                                      LorentzIndex]];
                         i = 0;
temp = Catch[
                         While[i < Length[temp], i++;
                             If[!FreeQ2[temp[[i]],{Eps,Pair}],
                                lind = liget[ temp[[i]] ];
                                If[Length[lind]===6,
 (* create a list of 5 possible arrangements of terms of Schouten ident. *)
                                schou = lisch /@ Map[ 
                                    Append[#, Last[lind]]&,
                                    NestList[RotateLeft, Take[lind, 5], 4]
                                                    ];
                                sor = Sort[ schou, epc ]//Reverse;
                                Do[ If[FreeQ[temp, sor[[1]]], 
                                       sor = RotateLeft[sor]
                                      ], {6}];
print3["sor = ", sor];

                                If[!FreeQ[temp, sor[[1]]],
                                   ntemp = Expand[EpsEvaluate[temp/.sor[[1]]->
                                          (-Apply[Plus, Drop[sor,1]])
                                               ] ];
                                   If[(epsnterms[ntemp] < ltemp) || 
(* or all LorentzIndices are inside all Eps's *)
   (Union[Length[ Position[#, LorentzIndex] ]& /@
          Select[ Variables[ntemp], Head[#]===Eps& ]
         ]  === {numberofli}),         Throw[temp=ntemp]
                                     ]
                                  ]]]
                                ];
             temp];
     result = all[[1]]  + temp0[[1]] + temp
                            ] ];
         result] ];
                       
          
 
 (* #################################################################### *)
 (*                             Main27                                   *)
 (* #################################################################### *)

(* ************************************************************** *)
(* All you need for Dirac algebra                                 *)
(* ************************************************************** *)
   ds[x___] := memset[ ds[x], dr[x] ];    (*dsdef*)
(* ************************************************************** *)
(* The master function for Dirac-structures:  DiracSimplify       *)
(* ************************************************************** *)
 Options[diracSimplify] =
     {DiracInfo->False, DiracCanonicalFlag->False,
      InsideDiracTrace->False, DiracSubstitute67->False,
      Factoring -> False, DiracSimpCombine->False
     };                                         (*DiracSimplifydef*)
 dit[x_]:=DiracTrace[ diracSimplify[x,
                      Flatten[Prepend[{Options[DiracSimplify]},
                                      InsideDiracTrace->True]]
                                   ] ];
 Options[DiracSimplify] = {Expanding -> True};
 DiracSimplify[x_,y__, z___Rule]:=DiracSimplify[x.y, z];
 diracSimplify[z_]:=(Contract[z]/.DiracTrace->dit)/;!FreeQ[z,DiracTrace];

 DiracSimplify[a_, opts___Rule] := 
   If[ (Expanding /. {opts} /. Options[DiracSimplify]) === False,
       DiracEquation[
       dotLin[a] /. Eps -> ePSIN /. Dot -> drS /. drS -> ds /. 
                     dr -> drS /. drS -> ds /. dr -> Dot /. ePSIN-> Eps
                    ], 
       If[FreeQ[a, SUNT],
          oldDiracSimplify[a /.Dot->drS /. drS -> Dot] /. Pair -> sCO /.
          sCO -> Pair,
          oldDiracSimplify[a /.Dot->drS /. drS -> ds /. dr -> Dot
                          ] /. Pair -> sCO /. sCO -> Pair
         ]
     ];

 oldDiracSimplify[x_] := diracSimplify[x]/;FreeQ[x,Spinor];
 oldDiracSimplify[x_] := Block[{dre},
                     dre =  FixedPoint[ SpinorChainEvaluate, x, 142];
                     If[ !FreeQ[dre, Eps], 
                         dre = Contract[dre, EpsContract->True];
                         dre = FixedPoint[ SpinorChainEvaluate, dre, 142]
                         ,
                         If[!FreeQ[dre, LorentzIndex],
(*
                            dre = Contract[dre]
*)
                            dre = Contract[dre, Expanding -> False]
                           ];
                         dre = FixedPoint[ SpinorChainEvaluate, dre, 142];
                       ];
   If[!FreeQ[dre, LorentzIndex],
print2["contracting in oldDiracSimpify"];
      dre = Contract[dre];
print2["contracting in oldDiracSimpify done"];
     ];
 
   If[Length[DownValues[Feyn`Calc`OneLoop`spinorsandpairs]
            ] > 1,
      dre = (dre /. Dot -> Feyn`Calc`OneLoop`spinorsandpairs/.
             Feyn`Calc`OneLoop`spinorsandpairs->Dot
            )//dotLin
     ];
(*
   If[FreeQ[dre,Dot] || (!FreeQ[dre,StandardMatrixElement]), 
      dre = Expand2[dre, StandardMatrixElement] ];
*)
   If[!FreeQ[dre, DiracGamma], dre = Expand2[dre, DiracGamma]];
   If[LeafCount[dre] < 420, dre = Factor2[dre,FactorTime->10]];
       dre                 ] /; !FreeQ[x,Spinor];

 collone[x_,y_]:=(Print["collone"]; xx = x;Collect2[x,y, Factoring ->False]);

(* #################################################################### *)
(*                             Main28                                 *)
(* #################################################################### *)

(*diracSimplifydef *)
diracSimplify[x_,in___]:=x/;FreeQ2[x, {DiracGamma, ChargeConjugationMatrix,
                                       ChargeConjugationMatrixInverse}];
diracSimplify[x_,in___]:= memset[diracSimplify[x,in], Block[
       {diracopt,diracdt,diracndt=0,diraccanopt,diracpdt,diracgasu,
        diracldt,diracjj=0,info,diractrlabel,diracga67,diracsifac,
        diracpag,colle
       },
        (* There are several options *)
        diracopt      = Join[ Flatten[{in}],Options[diracSimplify] ];
        info         = DiracInfo/.diracopt;
        diraccanopt  = DiracCanonicalFlag/.diracopt;
        diractrlabel = InsideDiracTrace/.diracopt;
        diracga67    = DiracSubstitute67/.diracopt;
        diracgasu    = DiracSimpCombine/.diracopt;
        diracsifac   = Factoring/.diracopt;
        diracdt = dotLin[ x//DiracGammaExpand ];
print3["dir1"];
        If[ diracgasu === True,
            diracdt = contractli[diracgammacombine[diracdt/.Pair->sCO]
                                ]/.Dot->ds,
            diracdt = contractli[ diracdt ]/.Dot->ds
          ];
print3["dir2a"];
        diracdt = Expand2[ pairexpand[diracdt//fEx],dr ];
        If[diractrlabel===True,
(*
           diracdt = diracdt/.dr->trIC/.trI->ds//.
                     dr->drCOs/.drCO->trIC/.trI->ds;
*)
           diracdt = diracdt/.dr->trIC/.trI->ds;
              (* optimization *)
           colle[a_]:=If[ (Length[a]<20(*00*))||(Head[a]=!=Plus), a,   
                          Collect2[a, dr,Factoring -> False] ];
           dirfun[exp_]:=colle[ exp/.dr->drCOs/.drCO->dr/.
                                        dr->trIC/.trI->ds
                              ];
           diracdt = FixedPoint[ dirfun, diracdt ]/.drCOs->trIC/.trI->ds;
print3["dir2done"];
print["diracdt = ",diracdt];
           If[ FreeQ[ diracdt, dr ],
               diracdt = diracdt/.DiracGamma[_[__],___]->0;
               diracpag=PartitHead[diracdt,DiracGamma];
                   If[ diracpag[[2]] === DiracGamma[5], diracdt = 0 ];
                   If[ diracpag[[2]] === DiracGamma[6] ||
                       diracpag[[2]] === DiracGamma[7],
                       diracdt = 1/2  diracpag[[1]]
                     ]
             ]
          ];
print3["dir3"];
        If[FreeQ[diracdt,dr],
           diracndt=Expand[(diracdt/.sCO->scev)//DiracGammaExpand];
           If[ diracga67 === True, diracndt = Expand[ diracndt//gamma67back ] ],
         diracdt = Expand2[ diracdt,dr ];
         If[ Head[diracdt] === Plus, diracldt=Length[diracdt],
             If[ diracdt===0, diracldt = 0, diracldt = 1 ]
           ];
print3["in diracSimplify: working with ",diracldt," terms"];
      While[diracjj<diracldt,diracjj++;
            If[diracldt==1,
               diracpdt = diracdt, diracpdt = diracdt[[diracjj]]
              ];
            If[diractrlabel===True,
               diracpdt = diracpdt/.dr->trIC/.trI->ds//.dr->drCOs/.
                          drCO->trIC/.trI->ds;
               diracpdt = diracpdt//.dr->drCOs/.drCO->ds
              ];
(* maybe insert some TimeConstrained here later *)
print3["in diracSimplify: contraction done, expand now."];
       diracpdt = pairexpand[ diracpdt ]//Expand;
(*
            diracpdt = Expand2[ diracpdt,dr ];
*)
            If[diractrlabel===True,
               diracpdt = fEx[(diracpdt//DiracGammaExpand)/.dr->Dot]/.
                                    dr->trIC/.trI->ds//.dr->drCOs/.
                                    drCO->trIC/.trI->ds,
               diracpdt = fEx[DiracGammaExpand[diracpdt]/.dr->Dot]//.
                                    dr->drCOs/.drCO->ds
              ];
             If[ diracga67===True,
                 diracpdt = gamma67back[ diracpdt/.dr->dr67 ],
                 diracpdt = fEx[ diracpdt ]
               ];
             diracndt = diracndt + Expand[ diracpdt ];
             If[ info===True,
                 print["# ",diracjj," / ",diracldt," = ",
                        Length[diracndt] ]
               ]
           ];
   diracndt = diracndt/.dr->Dot/.sCO->scev;
   diracndt = Expand2[dotLin[diracndt],Dot];
   If[ (diraccanopt===True ), 
print3["diracordering in diracSimplify"];
        diracndt = DiracOrder[ diracndt ] ;
        diracndt = Expand2[dotLin[diracndt],Dot]
     ];
          ] (* If FreeQ[diracdt,dr] *);
print3["dir4"];
print["diracdt = ", diracdt ];
    diracndt = dotLin[diracndt];
   If[ diracsifac === True,
       diracndt = Factor2[ diracndt ] ];
print3["exiting diracSimplify"];
  diracndt]](* /; !FreeQ[x, DiracGamma]*);;  (* end of diracSimplify *)

(* #################################################################### *)
(*                             Main29                                   *)
(* #################################################################### *)

  (*DiracGammaExpanddef*)
 DiracGammaExpand[x_]:=DiracGammaExpand[relh[x]]/;FreeQ[x,Momentum]&&
                                                 !FreeQ[x,HoldForm];
 DiracGammaExpand[x_]:=x/.DiracGamma->gaev/.gaevlin -> DiracGamma;
   gaev[x__] := Apply[ gaevlin, ExpandAll[{x}//MomentumExpand] ];  (*gaevdef*)
   gaevlin[ n_Integer ]  := DiracGamma[n]; (* necessary !!!!!! *)
   gaevlin[x_ + y_,di___]      := gaevlin[x,di] + gaevlin[y,di];
   gaevlin[ z_. w_[y__],di___ ]    :=  z DiracGamma[ w[y],di ]/;
                                   (w===Momentum || w===LorentzIndex);

(* #################################################################### *)
(*                             Main30                                 *)
(* #################################################################### *)

 diracgammacombine[exp_]:=If[LeafCount[exp]<1000, 
                             DiracGammaCombine[exp], exp
                            ];

 DiracGammaCombine[exp_]:=exp//.gasumrules; (*DiracGammaCombinedef*)
(* in order to merge sums of DiracGamma's into one *)
   gasumrules =
    {n1_. DiracGamma[Momentum[x_,di___],di___] +
     n2_. DiracGamma[Momentum[y_,di___],di___]:>
             DiracGamma[ Momentum[n1 x + n2 y,di],di]/;
             (NumberQ[n1] && NumberQ[n2]),
     (n1_. DiracGamma[Momentum[x_,di___],di___] +
      n2_. DiracGamma[Momentum[x_,di___],di___] ):>
       (n1+n2) DiracGamma[Momentum[x,di],di],
     (n3_. Momentum[x_,di___] + n4_. Momentum[y_,di___]):>
       Momentum[ Expand[n3 x + n4 y],di]/;(NumberQ[n3]&&NumberQ[n4])
    };

(* ************************************************************** *)
(* all kinds of simplification functions                          *)
(* ************************************************************** *)
(* #################################################################### *)
(*                             Main31                                 *)
(* #################################################################### *)
            (* drdef *)
   dr[]=1;
   dr[a___,y_SUNT w_,b___]:= dr[a, y, w, b];
   dr[a___,y_ w_,b___]:=coneins[y ds[a,w,b]]/;(noncommQ[y]&&FreeQ[y,dr]);
   dr[a___,y_ ,b___]  :=coneins[y ds[a,b] ] /;(noncommQ[y]&&FreeQ[y,dr]);

If[ Global`$BreitMaison === True,
   dr[a__]:=( ds[a]/.DiracGamma[6]->(1/2 + DiracGamma[5]/2)/.
                     DiracGamma[7]->(1/2 - DiracGamma[5]/2)
            )/;(!FreeQ2[{a}, {DiracGamma[6], DiracGamma[7]}]) && 
                (Head[DiracGamma[6]]===DiracGamma);

   dr[b___,DiracGamma[7],DiracGamma[x_[c__],di___],d___ ] :=
     ( 1/2 ds[ b,DiracGamma[x[c],di],d ] -
       1/2 ds[ b,DiracGamma[5],DiracGamma[x[c],di],d ]
     )
  ];
   
   dr[b___,DiracGamma[5],DiracGamma[5],c___]:= ds[ b,c ];
   dr[b___,DiracGamma[5],DiracGamma[6],c___]:= ds[b,DiracGamma[6],c];
   dr[b___,DiracGamma[5],DiracGamma[7],c___]:=-ds[b,DiracGamma[7],c];
 
   dr[b___,DiracGamma[6],DiracGamma[x_[c__],di___],d___ ]:=
     ds[ b,DiracGamma[x[c],di], DiracGamma[7],d ];
 
   dr[b___,DiracGamma[6],DiracGamma[5],c___]:=ds[b,DiracGamma[6],c];
   dr[b___,DiracGamma[6],DiracGamma[7],c___ ] :=  0;
   dr[b___,DiracGamma[7],DiracGamma[6],c___ ] :=  0;
 
   dr[b___,DiracGamma[7],DiracGamma[x_[c__],di___],d___ ] :=
      ds[ b,DiracGamma[x[c],di],DiracGamma[6],d ];
 
   dr[b___,DiracGamma[7],DiracGamma[5],c___]:=-ds[b,DiracGamma[7],c];
 
   dr[b___,DiracGamma[6],DiracGamma[6],c___] :=
      ds[ b,DiracGamma[6],c ];

   dr[b___,DiracGamma[7],DiracGamma[7],c___] :=
      ds[ b,DiracGamma[7],c ];
 
   dr[b___,DiracGamma[5],c:DiracGamma[_[_]].. ,d___] :=
      (-1)^Length[{c}] ds[ b,c,DiracGamma[5],d];

(* o.k., some 4 years after the proposal of M.B., here it is: *)
   drS[b___,DiracGamma[7],DiracGamma[_[__],___] + (n_. mass_),
      xy:DiracGamma[_[__],___].. , DiracGamma[6], c___] := 
   (n mass drS[b, xy, DiracGamma[6], c]) /; NumberQ[n] &&
     OddQ[Length[{xy}]] && noncommQ[mass];
  
   drS[b___,DiracGamma[6],DiracGamma[_[__],___] + (n_. mass_ ),
      xy:DiracGamma[_[__],___].. , DiracGamma[7], c___] := 
   (n mass drS[b, xy, DiracGamma[7], c]) /; NumberQ[n] &&
     OddQ[Length[{xy}]] && noncommQ[mass];
 
   drS[b___,DiracGamma[6],DiracGamma[_[__],___] + (n_. mass_ ),
      xy:DiracGamma[_[__],___].. , DiracGamma[6], c___] := 
   (n mass drS[b, xy, DiracGamma[6], c]) /; NumberQ[n] &&
     EvenQ[Length[{xy}]] && noncommQ[mass];
 
   drS[b___,DiracGamma[7],DiracGamma[_[__],___] + (n_. mass_ ),
      xy:DiracGamma[_[__],___].. , DiracGamma[7], c___] := 
   (n mass drS[b, xy, DiracGamma[7], c]) /; NumberQ[n] &&
     EvenQ[Length[{xy}]] && noncommQ[mass];
 
   drS[b___,DiracGamma[6],DiracGamma[v_[w__],di___] + (n_. mass_ ),
       DiracGamma[6], c___] := 
   (n mass drS[b, DiracGamma[6], c] )/; NumberQ[n] && noncommQ[mass];
 
   drS[b___,DiracGamma[7],DiracGamma[v_[w__],di___] + (n_. mass_ ),
       DiracGamma[7], c___] := 
   (n mass drS[b, DiracGamma[7], c] )/; NumberQ[n] && noncommQ[mass];
 
   drS[b___,DiracGamma[6],DiracGamma[v_[w__],di___] + (n_. mass_ ),
       DiracGamma[7], c___] := 
   drS[b, DiracGamma[v[w],di], DiracGamma[7], c] /; NumberQ[n] &&
     noncommQ[mass];
  
   drS[b___,DiracGamma[7],DiracGamma[v_[w__],di___] + (n_. mass_),
       DiracGamma[6], c___] := 
   drS[b, DiracGamma[v[w],di], DiracGamma[6], c] /; NumberQ[n] &&
     noncommQ[mass];
  
   drS[b___,DiracGamma[6],DiracGamma[v_[w__],di___] + (n_. mass_ ),
       xy:DiracGamma[_[_]].. ,DiracGamma[7], c___] := 
   drS[b, DiracGamma[v[w],di], xy, DiracGamma[7], c] /; NumberQ[n] &&
          EvenQ[Length[{xy}]] && noncommQ[mass];
  
   drS[b___,DiracGamma[7],DiracGamma[v_[w__],di___] + (n_. mass_ ),
       xy:DiracGamma[_[__],___].. ,DiracGamma[6], c___] := 
   drS[b, DiracGamma[v[w],di], xy, DiracGamma[6], c] /; NumberQ[n] &&
          EvenQ[Length[{xy}]] && noncommQ[mass];
  
   drS[b___,DiracGamma[6],DiracGamma[v_[w__],di___] + (n_. mass_ ),
       xy:DiracGamma[_[__],___].. ,DiracGamma[6], c___] := 
   drS[b, DiracGamma[v[w],di], xy, DiracGamma[6], c] /; NumberQ[n] &&
          OddQ[Length[{xy}]] && noncommQ[mass];
  
   drS[b___,DiracGamma[7],DiracGamma[v_[w__],di___] + (n_. mass_),
       xy:DiracGamma[_[__],___].. ,DiracGamma[7], c___] := 
   drS[b, DiracGamma[v[w]], xy, DiracGamma[7], c] /; NumberQ[n] &&
          OddQ[Length[{xy}]] && noncommQ[mass];

(* NEW 10.9.93*)
   drS[b___, DiracGamma[LorentzIndex[c_,dI___],dI___],d_Plus, e___,
             DiracGamma[v_[w__],di___],
             DiracGamma[LorentzIndex[c_,dI___],dI___], f___] := (
   -drS[b, DiracGamma[LorentzIndex[c,dI], dI], d,e,
           DiracGamma[LorentzIndex[c,dI], dI], 
           DiracGamma[v[w], di], f
       ] + 2 sCO[LorentzIndex[c, dI], v[w]] *
             drS[b, DiracGamma[LorentzIndex[c, dI], dI], d, e, f])/;
     (And@@Map[ MatchQ[#, (n_. Pair[_,_] DiracGamma[__])/;NumberQ[n]
                      ]&, List@@d ]);

   drS[b___, DiracGamma[LorentzIndex[c_,dI___],dI___],d_Plus, 
             DiracGamma[LorentzIndex[c_,dI___],dI___], e___] := 
    (2-fdim[dI])  drS[b, d, e] /; 
     (And@@Map[ MatchQ[#, (n_. Pair[_,_] DiracGamma[__])/;NumberQ[n]
                      ]&, List@@d ]);
  
  If[ Global`$BreitMaison === False,
   dr[b___,DiracGamma[5],c:DiracGamma[_[__],_].. ,d___] :=
      (-1)^Length[{c}] ds[ b,c,DiracGamma[5],d ],
 
   dr[b___,DiracGamma[5],DiracGamma[x_[y__],d_Symbol -4] ,f___] :=
      ds[ b,DiracGamma[x[y],d-4],DiracGamma[5],f ];
   dr[b___,DiracGamma[5],DiracGamma[x_[y__],d_Symbol] ,f___] :=
      ( ds[b,DiracGamma[x[y],d-4],DiracGamma[5],f] -
        ds[b,DiracGamma[x[y],4],DiracGamma[5],f] )
    ];
 
 
 (* #################################################################### *)
 (*                             Main32                                 *)
 (* #################################################################### *)

(* gamma[mu] gamma[mu] ---> 4, etc. *)
   dr[b___,DiracGamma[LorentzIndex[c_]],
           DiracGamma[LorentzIndex[c_]],d___] := 4 ds[ b,d ];

   dr[b___,DiracGamma[LorentzIndex[c_,di_],di_],
           DiracGamma[LorentzIndex[c_,di_],di_],d___] := di ds[ b,d ];
 
   dr[b___,DiracGamma[LorentzIndex[c_,di_],di_],
           DiracGamma[LorentzIndex[c_,di_ -4],di_ -4],d___]:=(di-4) ds[ b,d ];
   
   dr[b___,DiracGamma[LorentzIndex[c_]],
           DiracGamma[LorentzIndex[c_,di_ -4],di_ -4],d___] := 0;
   
   dr[b___,DiracGamma[LorentzIndex[c_]],
           DiracGamma[LorentzIndex[c_,di_ ],di_ ],d___] := 4 ds[ b,d ];
 
fdim[]=4;    (* fdimdef *)
fdim[dimi_]:=dimi;
dcheck[dii_, diii__] := dimcheck[dii, diii] = 
If[Head[dii]===Symbol, True, If[Union[{dii, diii}]==={dii}, True, False]];
 
   dr[b___,DiracGamma[LorentzIndex[c_,dI___],dI___],
           DiracGamma[x_[y__],di1___],
           DiracGamma[LorentzIndex[c_,dI___],dI___],d___
     ] := ( (2-fdim[dI]) ds[b,DiracGamma[x[y],di1],d] ) /; dcheck[dI, di1];
 
   dr[b___,DiracGamma[LorentzIndex[c_,dI___],dI___],
           DiracGamma[x1_[y1__],d1___], DiracGamma[x2_[y2__],d2___],
           DiracGamma[LorentzIndex[c_,dI___],dI___],d___
     ] := (4 sCO[x1[y1],x2[y2]] ds[b,d] +
           (fdim[dI]-4) ds[b,DiracGamma[x1[y1],d1], DiracGamma[x2[y2],d2], d] 
          ) /; dcheck[dI, d1, d2];
 
   dr[b___,DiracGamma[LorentzIndex[c_,dI___],dI___],
           DiracGamma[x1_[y1__],d1___],
           DiracGamma[x2_[y2__],d2___],
           DiracGamma[x3_[y3__],d3___],
           DiracGamma[LorentzIndex[c_,dI___],dI___],d___
     ] := (-2 ds[b,DiracGamma[x3[y3],d3], DiracGamma[x2[y2],d2],
                   DiracGamma[x1[y1],d1],
               d] -
           (fdim[dI]-4) ds[b,DiracGamma[x1[y1],d1],
                             DiracGamma[x2[y2],d2],
                             DiracGamma[x3[y3],d3],
                         d]
          ) /; dcheck[dI, d1,d2,d3];
  dr[b___,DiracGamma[LorentzIndex[c_,dI___],dI___],
          DiracGamma[x1_[y1__],d1___],
          DiracGamma[x2_[y2__],d2___],
          DiracGamma[x3_[y3__],d3___],
          DiracGamma[x4_[y4__],d4___],
          DiracGamma[LorentzIndex[c_,dI___],dI___],d___
     ] := ( 2 ds[b,DiracGamma[x3[y3],d3],
                   DiracGamma[x2[y2],d2],
                   DiracGamma[x1[y1],d1],
                   DiracGamma[x4[y4],d4],
                 d] +
            2 ds[b,DiracGamma[x4[y4],d4],
                   DiracGamma[x1[y1],d1],
                   DiracGamma[x2[y2],d2],
                   DiracGamma[x3[y3],d3],
                 d] +
       (fdim[dI]-4) ds[b,DiracGamma[x1[y1],d1],
                         DiracGamma[x2[y2],d2],
                         DiracGamma[x3[y3],d3],
                         DiracGamma[x4[y4],d4],
                     d]
         ) /; dcheck[dI, d1,d2,d3,d4];
 dr[b___,DiracGamma[LorentzIndex[c_,dI___],dI___],
          DiracGamma[x1_[y1__],d1___],
          DiracGamma[x2_[y2__],d2___],
          DiracGamma[x3_[y3__],d3___],
          DiracGamma[x4_[y4__],d4___],
          DiracGamma[x5_[y5__],d5___],
          DiracGamma[LorentzIndex[c_,dI___],dI___],d___
    ] := ( 2 ds[b,DiracGamma[x2[y2],d2],
                  DiracGamma[x3[y3],d3],
                  DiracGamma[x4[y4],d4],
                  DiracGamma[x5[y5],d5],
                  DiracGamma[x1[y1],d1],
                d] -
           2 ds[b,DiracGamma[x1[y1],d1],
                  DiracGamma[x4[y4],d4],
                  DiracGamma[x3[y3],d3],
                  DiracGamma[x2[y2],d2],
                  DiracGamma[x5[y5],d5],
                d] -
           2 ds[b,DiracGamma[x1[y1],d1],
                  DiracGamma[x5[y5],d5],
                  DiracGamma[x2[y2],d2],
                  DiracGamma[x3[y3],d3],
                  DiracGamma[x4[y4],d4],
                d] -
      (fdim[dI]-4) ds[b,DiracGamma[x1[y1],d1],
                        DiracGamma[x2[y2],d2],
                        DiracGamma[x3[y3],d3],
                        DiracGamma[x4[y4],d4],
                        DiracGamma[x5[y5],d5],
                    d] ) /; dcheck[dI, d1,d2,d3,d4,d5];
 
   dr[b___,DiracGamma[Momentum[c_,dim1___],___],
           DiracGamma[Momentum[c_,dim2___],___],d___ ] :=
              scev[Momentum[c,dim1],Momentum[c,dim2]] ds[b,d];
 
   dr[ b___,DiracGamma[LorentzIndex[c_]],d:DiracGamma[_[_]].. ,
            DiracGamma[LorentzIndex[c_]],f___ ] :=
    -2 ds @@ Join[ {b},Reverse[{d}],{f} ] /; OddQ[Length[{d}]];
 
   dr[ b___,DiracGamma[Momentum[c__],dim___],
            DiracGamma[Momentum[x__],dii___],
            DiracGamma[Momentum[c__],di___],d___ ] := (
   2 scev[Momentum[c],Momentum[x]] ds[b,DiracGamma[Momentum[c],dim],d]
   - scev[Momentum[c],Momentum[c]] ds[b,DiracGamma[Momentum[x],dii],d]
                                                  );

(* #################################################################### *)
(*                             Main33                                 *)
(* #################################################################### *)

(* SUNstuff *)
   dr[ a___,b_,c:SUNT[i_].. ,d___] :=
     dr[ a, c, b, d ] /; FreeQ2[b, {SUNT}];
 
   Literal[dr[ a___,b_ dr[c:(SUNT[_])..], d___]]:=
     ( dr[c] dr[a, b, d] )/;FreeQ[{a, b, d}, SUNT];
 
   dr[ SUNT[i_], b___ ] :=
       (SUNT[i] ds[b]) /; FreeQ[{b}, SUNT];
 
   dr[ b__, SUNT[i_] ] :=
       (SUNT[i] ds[b]) /; FreeQ[{b}, SUNT];
 
   dr[ a__, b:SUNT[_].. ]:=(ds[b] ds[a])/; 
     FreeQ[{a}, SUNT];
 
   dr[ b:SUNT[_].., a__ ]:=(ds[b] ds[a])/; 
     FreeQ[{a}, SUNT];
(* #################################################################### *)
(*                             Main33a                                 *)
(* #################################################################### *)
   dr[ a___, ChargeConjugationMatrix, ChargeConjugationMatrix, b___ ] :=
     -dr[a, b];
   dr[ a___, ChargeConjugationMatrix, 
             ChargeConjugationMatrixInverse, b___ ] := dr[a, b];
   dr[ a___, ChargeConjugationMatrixInverse, 
             ChargeConjugationMatrix, b___ ] := dr[a, b];
   dr[ a___, ChargeConjugationMatrix, DiracGamma[5], b___ ] :=
     dr[a, DiracGammaT[5], ChargeConjugationMatrix, b];
   dr[ a___, ChargeConjugationMatrix, DiracGamma[6], b___ ] :=
     dr[a, DiracGammaT[6], ChargeConjugationMatrix, b];
   dr[ a___, ChargeConjugationMatrix, DiracGamma[7], b___ ] :=
     dr[a, DiracGammaT[7], ChargeConjugationMatrix, b];
   
   dr[ a___, ChargeConjugationMatrix, DiracGamma[x_], b___ ] :=
     -dr[a, DiracGammaT[x], ChargeConjugationMatrix, b] /; !NumberQ[x];

    
   



(* #################################################################### *)
(*                             Main34                                 *)
(* #################################################################### *)

   drCOs[x___] := memset[ drCOs[x], drCO[x] ];    (*drCOsdef*)
(* Dirac contraction rules *) (*drCOdef*)

   drCO[ b___,DiracGamma[lv_[c_,di_Symbol-4],di_Symbol-4], w___,
              DiracGamma[ww_[y__],dim___], 
              DiracGamma[lv_[c_,di_Symbol-4],di_Symbol-4], z___] :=
   (-drCO[ b, DiracGamma[lv[c,di-4],di-4],w,
             DiracGamma[lv[c,di-4],di-4],
             DiracGamma[ww[y],dim],z
        ] + 2 drCO[b, DiracGamma[ww[y],di-4], w,z] )/.drCO->ds;
       
       
   drCO[ b___,DiracGamma[LorentzIndex[c_]],d:DiracGamma[_[__]].. ,
         DiracGamma[x_[y__]],DiracGamma[LorentzIndex[c_]],f___ ] :=
       ( 2 ds @@ Join[ {b},Reverse[{d}],{DiracGamma[x[y]],f} ] +
         2 ds[ b,DiracGamma[x[y]],d,f ]
        ) /; OddQ[Length[{d}]];
 
   drCO[ b___,DiracGamma[c_, di___],d:DiracGamma[_[__],___].. ,
         DiracGamma[c_,dim___],f___
       ] :=
        Block[ {drCOij, drCOld = Length[{d}]},
     (-1)^drCOld scev[c,c] ds[b,d,f]
     + 2 Sum[(-1)^(drCOij+1) coneins[ Pair[c,{d}[[drCOij,1]] ]
            * ds@@Join[{b},Drop[{d},{drCOij,drCOij}],{DiracGamma[c,dim],f}]
                                    ],{drCOij,1,drCOld}
            ]
              ]/;((Length[{d}]>0)&&FreeQ[c,LorentzIndex]&&
                 (!NumberQ[c]) && !MatchQ[{di}, {_Symbol -4}]);

(* #################################################################### *)
(*                             Main35                                 *)
(* #################################################################### *)


   drCO[ b___,DiracGamma[LorentzIndex[c_,di_Symbol],di_Symbol],
         d:DiracGamma[_[_,dim___],dim___].. ,
         DiracGamma[LorentzIndex[c_,di_Symbol],di_Symbol],f___
       ]:=
   Block[{idrCO,jdrCO,lddrCO = Length[{d}]},
        (-1)^lddrCO ( di - 2 lddrCO ) ds[b,d,f] -
          4 (-1)^lddrCO  Sum[ (-1)^(jdrCO-idrCO) *
         coneins[ Pair[{d}[[idrCO,1]],{d}[[jdrCO,1]] ] *
                  ds@@Join[ {b},Drop[ Drop[{d},{idrCO,idrCO}],
                                     {jdrCO-1,jdrCO-1}
                                    ],{f}
                          ]
                ],
                       {idrCO,1,lddrCO-1},{jdrCO,idrCO+1,lddrCO}
                            ]/.Pair->scev
         ] /;(Length[{d}]>5);
 
   drCO[ b___,DiracGamma[lv_[c_,dim___],dim___],
              DiracGamma[vl_[x__],dii___],d___,
              DiracGamma[lv_[c_,di___],di___],f___
       ]:= -ds[b, DiracGamma[vl[x],dii], DiracGamma[lv[c,dim],dim], d,
                    DiracGamma[lv[c,di],di], f
                ] + 2 coneins[Pair[vl[x], lv[c,dim]] * 
                              ds[ b,d,DiracGamma[lv[c,di],di],f ]
                             ];
(* #################################################################### *)
(*                             Main36                                 *)
(* #################################################################### *)

                                                        (*dr67def*)
   dr67[ b___ ] := ds[ b ]/;FreeQ2[{b},{DiracGamma[6],DiracGamma[7]}];
   dr67[ b___,DiracGamma[6],z___ ] := 1/2 ds[b,z] +
                                      1/2 ds[ b,DiracGamma[5],z ];
   dr67[ b___,DiracGamma[7],z___ ] := 1/2 ds[b,z] -
                                      1/2 ds[ b,DiracGamma[5],z ];
 
   dIex[ a___,x_ + y_, b___] := ds[a,x,b] + ds[a,y,b];   (*dIexdef*)
                                                         (*dixdef*)

   dix[y_] :=  y/.dr->dIex/.dIex->ds;

(* #################################################################### *)
(*                             Main37                                   *)
(* #################################################################### *)

(* ************************************************************** *)
 SetAttributes[dr,Flat];   (* quite important!!! *)
(* ************************************************************** *)
 
(* This is the tricky function which does the expansion of the dr's *)
   fEx[z_]:=FixedPoint[ dix, z/.Dot->dr ];                (*fExdef*)
 
(* ************************************************************** *)
(* A function for bringing DiracGamma's into canonical order      *)
(* ************************************************************** *)
 

(* #################################################################### *)
(*                             Main38                                   *)
(* #################################################################### *)
 
DiracCanonical[ x_,y__ ]:=DiracCanonical[x.y];
   DiracCanonical[x_]:=memset[DiracCanonical[x],
         Block[{diraccanres},    (*DiracCanonicaldef*)
       diraccanres = x//.{
    de_[a___,DiracGamma[vl_[y__],di___],
             DiracGamma[lv_[z__],dim___],b___
       ] :>( (-ds[a,DiracGamma[lv[z],dim], DiracGamma[vl[y],di],b
                 ] +(  (2 sCO[vl[y],lv[z]] ds[a,b])/.sCO->scev  )
             )/.dr->de/.sCO->scev
           )/; !OrderedQ[{lv[y],vl[z]}]
                         };
(* change here in Expand : 24.5.93 *)
    diraccanres = Expand2[dotLin[ diraccanres ], DiracGamma
                        ] /. Pair -> sCO /. sCO->Pair;
(*
nonsense ....

    If[ Length[ diraccanres ]<42,
        If[ FreeQ[diraccanres,Dot],
            diraccanres = Factor2[ diraccanres ],
            diraccanres = Expand2[ diraccanres,Dot ]
          ]
      ];
*)
    diraccanres] ];

(* #################################################################### *)
(*                             Main39                                 *)
(* #################################################################### *)

(* ************************************************************** *)
(* and something for ordering it as desired                       *)
(* ************************************************************** *)
(*                                     DiracOrderdef *)
 DiracOrder[x_]:=DiracCanonical[x];
 DiracOrder[x_,y___,z_]:=DiracCanonical[x.y.z]/;Head[z]=!=List;
 DiracOrder[x_,y__,ord_List]:=DiracOrder[x.y,ord];
 DiracOrder[x_,ord_List]:=memset[DiracOrder[x,ord], 
                                 Block[
     {diracordrev=Reverse[ord], diracordz,
      diracordres=x,diracordi},
    Do[ diracordz = diracordrev[[diracordi]];
        diracordres = diracordres//.
            {de_[a___,DiracGamma[vl_[y__],di___],
             DiracGamma[lv_[diracordz0_,dime___],dim___],b___
       ] :>
   (  (-ds[a,DiracGamma[lv[diracordz0,dime],dim],
             DiracGamma[vl[y],di],b
          ]+
        ( 2 sCO[vl[y],lv[diracordz0,dime]] ds[a,b] )/.sCO->scev
      )/.dr->de
   ) /; !FreeQ[lv[diracordz0, dime], diracordz]
            }, {diracordi,1,Length[ord]}
      ];
      (Expand2[dotLin[diracordres], DiracGamma])/.Pair->sCO/.sCO->Pair]];


(* #################################################################### *)
(*                             Main40                                 *)
(* #################################################################### *)

(* ************************************************************** *)
(* these are the only definitions for the Head "DiracGamma";       *)
(* they are nessecary in order to                                 *)
(* perform all projections in the right way.                      *)
(* ************************************************************** *)
                                                 (* DiracGammadef *)

 DiracGamma[ x_ Momentum[pe_, di___],dii___] := 
  x DiracGamma[Momentum[pe, di], dii];
 DiracGamma[ x_ LorentzIndex[pe_, di___],dii___] := 
  x DiracGamma[LorentzIndex[pe, di], dii];
 DiracGamma[ x_[y_,___],4]:=DiracGamma[x[y] ];   (*  define 4 *)
 DiracGamma[ 5, __ ]:=DiracGamma[5];
 DiracGamma[ 6, __ ]:=DiracGamma[6];
 DiracGamma[ 7, __ ]:=DiracGamma[7];
 DiracGamma[ _, 0] := 0;
 DiracGamma[ 0 ] = 0;
 DiracGamma[ 0,_ ] := 0;
 DiracGamma[ LorentzIndex[x_],di_Symbol-4 ] := 0;(*    4, D-4 *)
 DiracGamma[ Momentum[x_],di_Symbol-4 ]     := 0;(*    4, D-4 *)
 DiracGamma[ Momentum[x_,di_Symbol -4]]     := 0;(*  D-4, 4   *)
 DiracGamma[ LorentzIndex[x_,di_Symbol -4]] := 0;(*  D-4, 4   *)
   
 DiracGamma[ LorentzIndex[x_,di_],di_Symbol-4 ]:= (*    D, D-4 *)
  DiracGamma[LorentzIndex[x,di-4],di-4];
 DiracGamma[ Momentum[x_,di_],di_Symbol-4 ]:=     (*    D, D-4 *)
  DiracGamma[ Momentum[x,di-4],di-4];
   
 DiracGamma[ LorentzIndex[x_,di_Symbol -4],di_Symbol ]:= (* D-4, D *)
  DiracGamma[LorentzIndex[x,di-4],di-4];
 DiracGamma[ Momentum[x_,di_Symbol -4],di_Symbol ]:=     (* D-4, D *)
  DiracGamma[Momentum[x,di-4],di-4];
   
 DiracGamma[ LorentzIndex[x_],di_Symbol]:=       (*    4, D *)
  DiracGamma[LorentzIndex[x]];
 DiracGamma[ n_. Momentum[x_],di_Symbol]:=       (*(n) 4, D *)
    ( n DiracGamma[Momentum[x]] ) /; NumberQ[n];
 
 DiracGamma[Momentum[x_,di_Symbol]] :=           (*    D, 4 *)
     DiracGamma[Momentum[x]];
 DiracGamma[LorentzIndex[x_,di_Symbol]] :=       (*    D, 4 *)
     DiracGamma[LorentzIndex[x]];
(* #################################################################### *)
(*                             Main40a                                  *)
(* #################################################################### *)
Unprotect[Transpose];
Transpose[x_] := Reverse[x] /.{ DiracGamma             :> DiracGammaT,
                       DiracGammaT            :> DiracGamma,
                      ChargeConjugationMatrix :> 
                   (- ChargeConjugationMatrix ),
                      ChargeConjugationMatrixInverse :>
                   (- ChargeConjugationMatrixInverse)
                     };
Protect[Transpose];
Unprotect[Inverse];
Inverse[x_] := x/.{ChargeConjugationMatrixInverse :> ChargeConjugationMatrix,
                   ChargeConjugationMatrix :> ChargeConjugationMatrixInverse};
(* #################################################################### *)
(*                             Main41                                  *)
(* #################################################################### *)
(*DiracSpinordef*)
DiracSpinor[a__,{1}]:=Spinor[a,1];

(* ************************************************************** *)
(* All kinds of stuff for pinors                                 *)
(* There is:  Spinor,             *)
(*             EpsChisholm                                        *)
(* ************************************************************** *)
(* Spinordef *)
 QuarkSpinor[p_, m_/;Head[m]=!=Pattern]     := QuarkSpinor[p,m]=Spinor[p,m,2];
 QuarkSpinor[n_. (p_ /;Head[p]=!=Pattern)]  := QuarkSpinor[n p]=Spinor[n p,0,2];
 LeptonSpinor[p_, m_/;Head[m]=!=Pattern]    := LeptonSpinor[p,m]=Spinor[p,m,1];
 LeptonSpinor[n_. (p_ /;Head[p]=!=Pattern)] := LeptonSpinor[n p]=Spinor[n p,0,1];

 Spinor[n_. x_, y___] :=( Spinor[n x,y]=Spinor[n Momentum[x], y] )/;
      (FreeQ2[{x,y},{Blank,BlankSequence,Momentum,HoldForm}]&& (n^2)===1 );

 Spinor[Momentum[x_,di_],m_,op___] := Spinor[Momentum[x],m,op];
 Spinor[kk_.+ n_. Momentum[ a_Plus ],m_, y___]:=
  Spinor[kk+ n Momentum[ a],m, y] = (
 Spinor[ MomentumExpand[kk + n Momentum[a]] ,m,y] );

 Spinor[n_. p_, _. Negligible[_], in___]:=Spinor[n p, 0, in];
 Spinor[n_. p_ /; FreeQ[p,Pattern]] := Spinor[n p,0,1];
 Spinor[p_, m_ /; FreeQ[m, Pattern]] := Spinor[p,m,1];
(* NOOOOO!!, no good for SUSY *)
(*
(* as a convention : *)
 Spinor[- Momentum[pe_],0,in__] := Spinor[ Momentum[pe],0,in ];
*)


(* #################################################################### *)
(*                             Main42                                   *)
(* #################################################################### *)
  eepsrule={ m_. Dot[ x___,DiracGamma[LorentzIndex[in_,diim___],___],y___] *
                    eeps[a___,LorentzIndex[in_,di___],b_,c___] :>
              -m Dot[x,DiracGamma[LorentzIndex[in,diim]],y] * 
               eeps[a,b,LorentzIndex[in,di],c]
           };
 epsspcrule0={  
            m_. Dot[ x___,DiracGamma[LorentzIndex[in_]],y___] *
               eeps[a_,b_,c_,LorentzIndex[in_]] :>
               (( - I m ( Dot[x,DiracGamma[a],DiracGamma[b],DiracGamma[c],
                           DiracGamma[5],y
                          ] -
                      scev[a,b] Dot[x,DiracGamma[c].DiracGamma[5],y] -
                      scev[b,c] Dot[x,DiracGamma[a].DiracGamma[5],y] +
                      scev[a,c] Dot[x,DiracGamma[b].DiracGamma[5],y]
                     )//SpinorChainEvaluate
                )//Expand
               ) /; FreeQ[{x,y}, Spinor[_, _Symbol,___]] 
             };
 epsspcrule={  
            m_. Dot[ x___,DiracGamma[LorentzIndex[in_]],y___] *
               eeps[a_,b_,c_,LorentzIndex[in_]] :>
               ( - I m ( Dot[x,DiracGamma[a],DiracGamma[b],DiracGamma[c],
                           DiracGamma[5],y
                          ] -
                      scev[a,b] Dot[x,DiracGamma[c].DiracGamma[5],y] -
                      scev[b,c] Dot[x,DiracGamma[a].DiracGamma[5],y] +
                      scev[a,c] Dot[x,DiracGamma[b].DiracGamma[5],y]
                     )//SpinorChainEvaluate
               )//Expand,
              DiracGamma[LorentzIndex[in_]] * 
              eeps[a_,b_,c_,LorentzIndex[in_]] :>
               ( - I ( Dot[x,DiracGamma[a],DiracGamma[b],DiracGamma[c],
                           DiracGamma[5],y
                          ] -
                      scev[a,b] Dot[x,DiracGamma[c].DiracGamma[5],y] -
                      scev[b,c] Dot[x,DiracGamma[a].DiracGamma[5],y] +
                      scev[a,c] Dot[x,DiracGamma[b].DiracGamma[5],y]
                     )//SpinorChainEvaluate
               )//Expand
            };

(* #################################################################### *)
(*                             Main421                                  *)
(* #################################################################### *)

(* Chisholmdef *)
 Chisholm[x_]:=FixedPoint[chish, x, 1];

 chish[x_]:=Block[{ sim, re, ind },
              ind = Unique[mu];
              re = x;
          sim[y_, index_]:=Contract[DiracSimplify[ 
                y/. {
                 Dot[a___, DiracGamma[lv1_[pe1_]],DiracGamma[lv2_[pe2_]],
                            DiracGamma[lv3_[pe3_]],b___ 
                      ] :>  (
               Pair[lv1[pe1], lv2[pe2]] Dot[a, DiracGamma[lv3[pe3]], b] -
               Pair[lv1[pe1], lv3[pe3]] Dot[a, DiracGamma[lv2[pe2]], b] +
               Pair[lv2[pe2], lv3[pe3]] Dot[a, DiracGamma[lv1[pe1]], b] +
               I Eps[ lv1[pe1],lv2[pe2],lv3[pe3],LorentzIndex[index] ]*
               Dot[a, DiracGamma[LorentzIndex[index]].
                      DiracGamma[5], b] 
                            )
                    }                            ], EpsContract->True
                                   ] // DiracOrder;                                          
             re = sim[re, ind];
          re];

(* #################################################################### *)
(*                             Main422                                  *)
(* #################################################################### *)


(*ChisholmSpinordef*)
 dsimp[x_]:=sirlin0[spcev0[x]];
 ChisholmSpinor[x_, choice_:0]:=memset[ChisholmSpinor[x,choice],
                             Block[{new=x, indi},
print3["entering ChisholmSpinor "];
  new = dotLin[new];
  If[choice===1, new = new/.{ Spinor[a__].b__ .Spinor[c__] * 
                              Spinor[d__].e__ .Spinor[f__]:>
                             nospinor[a].b.nospinor[c] * 
                              Spinor[d].e.Spinor[f]
                            }
    ];
  If[choice===2, new = new/.{ Spinor[a__].b__ .Spinor[c__] *
                              Spinor[d__].e__ .Spinor[f__]:>
                              Spinor[a].b.Spinor[c] *
                              nospinor[d].e.nospinor[f]
                            }
    ];

                    dsimp[Contract[dsimp[new/.{
      Literal[ Spinor[pe1_, m_, ql___] . DiracGamma[lv_[k_]] . h___ .
               Spinor[pe2_, m2_, ql___] ] :> Block[{indi},
                      indi = Unique[alpha];
     -1/Pair[pe1,pe2] ( Spinor[pe1, m, ql]. DiracGamma[pe1].
                        DiracGamma[lv[k]] . DiracGamma[pe2].h.
                        Spinor[pe2, m2, ql] -
                        Pair[pe1,lv[k]] Spinor[pe1, m, ql].
                            DiracGamma[pe2]. h . Spinor[pe2, m2, ql] -
                        Pair[lv[k],pe2] Spinor[pe1, m, ql].
                            DiracGamma[pe1] . h . Spinor[pe2, m2, ql]-
                          I Eps[pe1,lv[k],pe2,LorentzIndex[indi]] *
                        Spinor[pe1, m, ql].
                            DiracGamma[LorentzIndex[indi]].
                            DiracGamma[5].h.
                        Spinor[pe2, m2, ql]
                      )] }/.nospinor->Spinor], EpsContract->True] ] ]];
                              

   (*EpsChisholmdef*)
 EpsChisholm[x_] := x /; FreeQ[x, Eps] || FreeQ[x, DiracGamma]; 
 EpsChisholm[x_] := Block[{new=0, xx,i},
                          xx = Expand2[x,Eps];
                          If[ Head[xx]===Plus,
                              For[i=1, i<=Length[xx], i++,
                                  new = new + 
                                  ((xx[[i]]/.Eps->eeps)/.eepsrule/.
                                    eepsrule/.eepsrule/.eepsrule/.
                                    eepsrule/.
                                    epsspcrule0/.epsspcrule0/.
                                    epsspcrule/.epsspcrule)
                                 ],
                              new = (xx/.Eps->eeps)/.eepsrule/.
                                    eepsrule/.eepsrule/.eepsrule/.
                                    eepsrule;
                              new = new/.epsspcrule0/.epsspcrule0;
                              new = new/.epsspcrule/.epsspcrule
                            ];
                      new/.eeps->Eps] /; !FreeQ[x,Eps];
		      
(* #################################################################### *)
(*                             Main43                                   *)
(* #################################################################### *)

   spinlin[x_Plus]:=spinlin/@x;
   spinlin[a_] :=( (a/.Dot->dot)//.{
                 dot[x___,z_ b__,c___] :> z dot[x,b,c]/;noncommQ[z]===True,
                 dot[x___,z_ ,c___]    :> z dot[x,c]/;noncommQ[z]===True,
                 dot[x_Spinor,b___,c_Spinor,d_Spinor,e___,f_Spinor,g___]:>
                    dot[x,b,c] dot[d,e,f,g] }
                 )/.dot[]->1/.dot->Dot;
 SetAttributes[ SpinorChainEvaluate, Listable ];
                                           (*SpinorChainEvaluatedef*)
 SpinorChainEvaluate[y_]:=y /; FreeQ[y,Spinor];
 (* #################################################################### *)
 (*                             Main44                                   *)
 (* #################################################################### *)
 
 SpinorChainEvaluate[z_Plus]:= Block[{nz},
   nz = dotLin[z];
   If[Length[nz]>20, nz= Collect2[ nz, Spinor,Factoring -> False] ];
   If[Head[nz]=!=Plus, nz = SpinorChainEvaluate[nz],
      If[$sirlin =!= True, nz = Map[ spcev0, nz ],
         If[ FreeQ[nz, Literal[
                       Spinor[p1__] . 
                            (a__ /; FreeQ[{a}, DiracGamma[_,_]]
                            ) . Spinor[p2__] * 
                       Spinor[p3__] . (b__ /; FreeQ[{b}, DiracGamma[_,_]]
                            ) . Spinor[p4__] ]
                  ], nz = Map[ spcev0,nz ],
       nz = sirlin00[ Expand2[Map[ spcev0,z//sirlin0 ], Spinor] ]
           ] ] ];                  nz];
 SpinorChainEvaluate[x_]:= 
  If[$sirlin =!= True, Expand2[spcev0[x], Spinor],
  If[ FreeQ[x//dotLin, Literal[
                       Spinor[p1__] .
                            (a__ /; FreeQ[{a}, DiracGamma[_,_]]
                            ) . Spinor[p2__] *
                       Spinor[p3__] . (b__ /; FreeQ[{b}, DiracGamma[_,_]]
                            ) . Spinor[p4__] ]
           ],  
     Expand2[spcev0[x],Spinor],
     sirlin00[ Expand2[FixedPoint[spcev0, x//sirlin0, 3 ], Spinor] ]
    ]]/; !Head[x]===Plus;

(* Reference of Sirlin-relations: Nuclear Physics B192 (1981) 93-99; 
   Note that we take another sign in front of the Levi-Civita tensor
   in eq. (7), since we take (implicitly) \varepsilon^{0123} = 1
*)

(* This may be useful, but quit slow, therefore: 
*)
 (* #################################################################### *)
 (*                             Main441                                  *)
 (* #################################################################### *)

  $SpinorMinimal = False;

  sirlin00[x_]:= x/;($SpinorMinimal === False) || ($sirlin===False);
  sirlin00[x_]:=memset[sirlin00[x],
                     Block[{te, tg5, ntg5},
print3["sirlin001"];
                       te = sirlin0[x]//ExpandAll;
print3["sirlin002"];
                       If[FreeQ2[te,{DiracGamma[6],DiracGamma[7]}]&&
                          Head[te]===Plus && !FreeQ[te,DiracGamma[5]],
                          tg5 = Select[te, !FreeQ[#,DiracGamma[5]]& ];
                          ntg5 = te - tg5;
(*i.e. te = tg5 + ntg5 *)
                          test = Expand[tg5 + ChisholmSpinor[ntg5]];
                          If[nterms[test] < Length[te], te=test]
                         ];
print3["exiting sirlin00"];
                  te]] /; $SpinorMinimal ===  True;

(* ident3def *)

ident3[a_,_]:=a;
                  
 (* #################################################################### *)
 (*                             Main442                                  *)
 (* #################################################################### *)
 (* canonize different dummy indices *)  (*sirlin3def*)
 sirlin3a[x_]:=((sirlin3[Expand2[Contract[x],Spinor]/.
                         $MU->dum$y]/.dum$y->$MU)/.  sirlin3 -> Identity
	       )//Contract;
 sirlin3[a_Plus]:=sirlin3 /@ a;
 Literal[
 sirlin3[ m_. Spinor[p1__]. (ga1___) . 
	     DiracGamma[ LorentzIndex[la_] ]. (ga2___) .
	     Spinor[p2__] *
	     Spinor[p3__]. (ga3___) .
	     DiracGamma[ LorentzIndex[la_] ]. (ga4___) .
             Spinor[p4__]
      ]] := Block[{counter},
                   counter = 1;

             While[!FreeQ2[{m,ga1,ga2,ga3,a4}, 
                           {$MU[counter], dum$y[counter]} ],
                   counter = counter + 1
                  ];
       sirlin3[
         m Spinor[p1] . ga1 .
         DiracGamma[ LorentzIndex[$MU[counter]] ] . ga2 .  Spinor[p2] *
         Spinor[p3] . ga3 .  DiracGamma[ LorentzIndex[$MU[counter]] ] . ga4 .
         Spinor[p4] 
              ]  ] /; FreeQ[la, $MU];

 Literal[
 sirlin3[ m_. Spinor[p1__].(ga1___). 
             DiracGamma[ LorentzIndex[la_,di_],di_ ]. (ga2___) .
             Spinor[p2__] *
             Spinor[p3__].(ga3___). 
             DiracGamma[ LorentzIndex[la_,di_],di_ ]. (ga4___) .
             Spinor[p4__]
      ]] := ( m Spinor[p1] . ga1 . 
                 DiracGamma[ LorentzIndex[$MU[1], di],di ] . ga2 . Spinor[p2] *
                 Spinor[p3] . ga3 . DiracGamma[ LorentzIndex[$MU[1], di], di 
                                              ] . ga4 .
                 Spinor[p4]
              ) /; FreeQ2[{ga1,ga2,ga3,ga4}, DiracGamma[_,_]];

              
(* this is far from optimal, but for the moment sufficient *)
 $sirlin = True;


 (* #################################################################### *)
 (*                             Main443                                  *)
 (* #################################################################### *)

(* The Sirlin - identities are only valid in 4 dimensions and are only needed,
   if Dirac matrices are around 
*)
 sirlin0[x_]:=If[$sirlin=!=True, x,
                 If[ FreeQ2[x, {LorentzIndex, Momentum}],  x,
                     If[ FreeQ[x, Spinor], x,
                         If[ !FreeQ[x, DiracGamma[_,_]], 
                             sirlin3[x]/.sirlin3->Identity,
                             sirlin0doit[(x//sirlin2)/.sirlin2->Identity] 
                   ]   ]   ]
                ];

$sirlintime = 242;

 sirlin0doit[a_Plus]:=timeconstrained[
sirlin3a[Contract[
		   (Expand2[Map[sirlin1, a], Dot]/.
		    sirlin1->sirlin2) /. 
		   sirlin2 -> sirlin1/.sirlin1->sirlin2/.
                    sirlin2 -> Identity,EpsContract->True]
			 ] // spcev0,
                                     2 $sirlintime, a
                                    ];
 sirlin0doit[a_]:=timeconstrained[  
                    (sirlin3a[sirlin1[a]/.sirlin1->sirlin2/.
                        sirlin2 -> Identity
                       ] // spcev0),
                                  $sirlintime, a 
                                 ] /;Head[a]=!=Plus;

(*sirlin2def*)
 sirlin2[a_Plus]:=sirlin2/@a;

 Literal[
 sirlin2[m_. Spinor[pa__] . DiracGamma[Momentum[pj_]] .
                            DiracGamma[Momentum[pi_]] .
                            DiracGamma[LorentzIndex[mu_]].(vg5___).
             Spinor[pb__] *
             Spinor[Momentum[pi_],0,qf___] .
                    DiracGamma[LorentzIndex[mu_]] . (vg5___).
             Spinor[Momentum[pj_],0,qf___]
       ]] := (-sirlin2[ m Spinor[pa] . DiracSlash[pi,pj] .
                                       DiracMatrix[mu] . vg5 .
                          Spinor[pb] *
                          Spinor[Momentum[pi],0,qf] . 
                                       DiracMatrix[mu] . vg5 .
                          Spinor[Momentum[pj],0,qf]
                      ] + 
                2 m scev[Momentum[pi],Momentum[pj]] * 
                Spinor[pa] . DiracMatrix[mu] . vg5 .
                Spinor[pb] *
                          Spinor[Momentum[pi],0,qf] .
                                       DiracMatrix[mu] . vg5 .
                          Spinor[Momentum[pj],0,qf]
             )/; ({vg5}==={}) || ({vg5}==={DiracGamma[5]});
         

 Literal[
 sirlin2[m_. Spinor[pa__] . DiracGamma[Momentum[pi_]] .
                            DiracGamma[Momentum[pj_]] .
                            DiracGamma[LorentzIndex[mu_]].(vg5___).
             Spinor[pb__] *
             Spinor[Momentum[pi_],0,qf___] . 
                    DiracGamma[LorentzIndex[mu_]] . (vg5___).
             Spinor[Momentum[pj_],0,qf___]
       ]] :=(m scev[Momentum[pi], Momentum[pj]] * 
              Spinor[pa] . DiracMatrix[$MU[1]] .
              Spinor[pb] *
              Spinor[Momentum[pi],0,qf] . DiracMatrix[$MU[1]] . 
              Spinor[Momentum[pj],0,qf] +
             m scev[Momentum[pi], Momentum[pj]] *
              Spinor[pa] . DiracMatrix[$MU[1]]. DiracGamma[5] .
              Spinor[pb] *
              Spinor[Momentum[pi],0,qf] . DiracMatrix[$MU[1]] .
              DiracGamma[5] . Spinor[Momentum[pj],0,qf]
            ) /; ({vg5}==={}) || ({vg5}==={DiracGamma[5]});


 Literal[
 sirlin2[m_. Spinor[p1__]. (ga1___) .
	     DiracGamma[ LorentzIndex[la_] ].
	     DiracGamma[ LorentzIndex[nu_] ].
	     DiracGamma[6] . 
	     Spinor[p2__] *
	     Spinor[p3__]. (ga2___) .
	     DiracGamma[ LorentzIndex[la_] ].
	     DiracGamma[ LorentzIndex[nu_] ].
	     DiracGamma[7] .
	     Spinor[p4__] ]] :=  (
    m 4 Spinor[p1] . ga1 . DiracGamma[6] . Spinor[p2] *
        Spinor[p3] . ga2 . DiracGamma[7] . Spinor[p4] );

 Literal[
 sirlin2[m_. Spinor[p1__]. (ga1___) .
	     DiracGamma[ LorentzIndex[la_] ].
	     DiracGamma[ LorentzIndex[nu_] ].
	     DiracGamma[7] . 
	     Spinor[p2__] *
	     Spinor[p3__]. (ga2___) .
	     DiracGamma[ LorentzIndex[la_] ].
	     DiracGamma[ LorentzIndex[nu_] ].
	     DiracGamma[6] .
	     Spinor[p4__] ] ] :=  (
    m 4 Spinor[p1] . ga1 . DiracGamma[7] . Spinor[p2] *
        Spinor[p3] . ga2 . DiracGamma[6] . Spinor[p4] );
 (* #################################################################### *)
 (*                             Main444                                  *)
 (* #################################################################### *)


(* eq. (8) *)
 Literal[
 sirlin2[m_. Spinor[p1__]. (ga1___) . 
              DiracGamma[ LorentzIndex[mu_] ]. 
              DiracGamma[ lv_[rho_] ] .
              DiracGamma[ LorentzIndex[nu_] ]. (ga2___) .
            Spinor[p2__] * 
            Spinor[p3__]. (ga3___) . 
              DiracGamma[ LorentzIndex[mu_] ].
              DiracGamma[ lvt_[tau_] ] .
              DiracGamma[ LorentzIndex[nu_] ]. (ga4___) .
            Spinor[p4__]
       ]] := Block[{ii=1, ind, la, grho, gtau, gam5},
                    While[!FreeQ[{ga1,ro,ga2,ga3,tau,ga4}, $MU[ii]],
                          ii++];
             la = DiracGamma[LorentzIndex[$MU[ii]]];
             grho = DiracGamma[lv[rho]]; gtau = DiracGamma[lvt[tau]];
             gam5 = DiracGamma[5];
             Contract[
               2 m Pair[lv[rho], lvt[tau]] * 
                   Spinor[p1] . ga1 . la . ga2 .   Spinor[p2] * 
                   Spinor[p3] . ga3 . la . ga4 .   Spinor[p4] +
               2 m *
                   Spinor[p1] . ga1 . gtau . ga2 . Spinor[p2] *
                   Spinor[p3] . ga3 . grho . ga4 .   Spinor[p4] +
               2 m Pair[lv[rho], lvt[tau]] *
                   Spinor[p1] . ga1 . la . ga2 . gam5 . Spinor[p2] *
                   Spinor[p3] . ga3 . la . ga4 . gam5 . Spinor[p4] -
               2 m *
                   Spinor[p1] . ga1 . gtau . ga2 . gam5 . Spinor[p2] *
                   Spinor[p3] . ga3 . grho . ga4 . gam5 . Spinor[p4] 
                     ]
                   ];

(* eq. (12) of Sirlin *)

 Literal[
 sirlin2[m_. Spinor[p1__]. DiracGamma[ LorentzIndex[mu_] ].
                           DiracGamma[ lv_[rho_] ] .
                           DiracGamma[ LorentzIndex[sigma_] ].
                           DiracGamma[ lvt_[tau_] ].
                           DiracGamma[ LorentzIndex[nu_] ]. om_ .
             Spinor[p2__] *
             Spinor[p3__]. DiracGamma[ LorentzIndex[mu_] ].
                           DiracGamma[ lva_[alpha_] ] .
                           DiracGamma[ LorentzIndex[sigma_] ].
                           DiracGamma[ lvb_[beta_] ].
                           DiracGamma[ LorentzIndex[nu_] ]. om_ .
             Spinor[p4__]
       ]] := Contract[ m 16 Pair[lvt[tau],lvb[beta]] * 
                            Pair[lv[rho], lva[alpha]] * 
                           Spinor[p1] . DiracMatrix[mu] . om .
                           Spinor[p2] *
                           Spinor[p3] . DiracMatrix[mu] . om .
                           Spinor[p4]
                     ] /; (om===DiracGamma[6]) || 
                          (om===DiracGamma[7]);

(* eq. (13) of Sirlin *)
 Literal[
 sirlin2[m_. Spinor[p1__]. DiracGamma[ LorentzIndex[mu_] ].
                           DiracGamma[ lv_[rho_] ] .
                           DiracGamma[ LorentzIndex[sigma_] ].
                           DiracGamma[ lvt_[tau_] ].
                           DiracGamma[ LorentzIndex[nu_] ]. om1_ .
             Spinor[p2__] *
             Spinor[p3__]. DiracGamma[ LorentzIndex[mu_] ].
                           DiracGamma[ lva_[alpha_] ] .
                           DiracGamma[ LorentzIndex[sigma_] ].
                           DiracGamma[ lvb_[beta_] ].
                           DiracGamma[ LorentzIndex[nu_] ]. om2_ .
             Spinor[p4__]
       ]] :=(m 4 Spinor[p1] . DiracMatrix[mu].DiracGamma[lv[rho]].
                              DiracGamma[lv[beta]]. om1 .
                 Spinor[p2] *
                 Spinor[p3] . DiracMatrix[mu].DiracGamma[lva[alpha]].
                              DiracGamma[lvt[tau]]. om2 .
                                            Spinor[p4]
            ) /; ( (om1===DiracGamma[6])&& (om2===DiracGamma[7]) )||
                 ( (om1===DiracGamma[7])&& (om2===DiracGamma[6]) );
 (* #################################################################### *)
 (*                             Main445                                  *)
 (* #################################################################### *)


(* in case if no chiral projectors are present: *)
 Literal[
 sirlin2[m_. Spinor[p1__]. DiracGamma[ LorentzIndex[mu_] ].
                           DiracGamma[ lv_[rho_] ] .
                           DiracGamma[ LorentzIndex[sigma_] ].
                           DiracGamma[ lvt_[tau_] ].
                           DiracGamma[ LorentzIndex[nu_] ].
             Spinor[p2__] *
             Spinor[p3__]. DiracGamma[ LorentzIndex[mu_] ].
                           DiracGamma[ lva_[alpha_] ] .
                           DiracGamma[ LorentzIndex[sigma_] ].
                           DiracGamma[ lvb_[beta_] ].
                           DiracGamma[ LorentzIndex[nu_] ].
             Spinor[p4__]
       ]] := Block[{tmp,re},
                    tmp[ome1_,ome2_]:= sirlin2[ m Spinor[p1].
   DiracMatrix[mu].DiracGamma[lv[rho]].DiracMatrix[sigma].
   DiracGamma[lvt[tau]].DiracMatrix[nu].DiracGamma[ome1] .
   Spinor[p2] * 
   Spinor[p3].DiracMatrix[mu].DiracGamma[lva[alpha]].
   DiracMatrix[sigma].DiracGamma[lvb[beta]].DiracMatrix[nu].
   DiracGamma[ome2].  Spinor[p4]              ];
                   re = tmp[6,6] + tmp[6,7] + tmp[7,6] + tmp[7,7];
               re];

 (* #################################################################### *)
 (*                             Main446                                  *)
 (* #################################################################### *)

(* These are the ones calculated by FeynCalc  *)

 Literal[
sirlin2[
m_.  Spinor[pi__] . x1___ . DiracGamma[ LorentzIndex[mu_] ] . 
               DiracGamma[ LorentzIndex[nu_] ] . x2___ .
Spinor[pj__] *
Spinor[pk__] .  x3___ . DiracGamma[ vm_[a_] ] .  
                DiracGamma[ LorentzIndex[mu_] ] .
               DiracGamma[ LorentzIndex[nu_] ] . x4___ .
Spinor[pl__]
       ]] := Contract[ m (
2*Spinor[pi] . x1 . x2 . Spinor[pj]*
   Spinor[pk] . x3 . DiracGamma[vm[a]] . x4 . 
    Spinor[pl] + 
  2*Spinor[pk] . x3 . DiracGamma[LorentzIndex[al$mu]] . x4 .
    Spinor[pl]*
   Spinor[pi] . x1 . DiracGamma[vm[a]] . 
    DiracGamma[LorentzIndex[al$mu]] . x2 . Spinor[pj] - 
  2*Spinor[pi] . x1 . DiracGamma[5] . x2 .
    Spinor[pj]*
   Spinor[pk] . x3 . DiracGamma[vm[a]] . DiracGamma[5] . x4 .
    Spinor[pl] + 
  2*Spinor[pk] . x3 . DiracGamma[LorentzIndex[al$mu]] . 
    DiracGamma[5] . x4 .Spinor[pl]*
   Spinor[pi] . x1 .  DiracGamma[vm[a]] . 
    DiracGamma[LorentzIndex[al$mu]] . DiracGamma[5] . x2 . Spinor[pj]
             )];

sirlin2[ m_. *
Spinor[Momentum[pi_], 0, fq___] . DiracGamma[Momentum[pk_]] . 
     Spinor[Momentum[pj_], 0, fq___]*
    Spinor[Momentum[pl_], 0, fq___] . DiracGamma[Momentum[pj_]] . 
     Spinor[Momentum[pk_], 0, fq___] 
       ] := Contract[ m (
   -((Spinor[Momentum[pi], 0, fq] . DiracGamma[Momentum[pl]] . 
          Spinor[Momentum[pj], 0, fq]*
         Spinor[Momentum[pl], 0, fq] . DiracGamma[Momentum[pi]] . 
          Spinor[Momentum[pk], 0, fq]*Pair[Momentum[pj], Momentum[pk]])/
       Pair[Momentum[pi], Momentum[pl]]) + 
    (Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
        DiracGamma[5] . Spinor[Momentum[pj], 0, fq]*
       Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
        DiracGamma[5] . Spinor[Momentum[pk], 0, fq]*
       (-(Pair[Momentum[pi], Momentum[pl]]*
            Pair[Momentum[pj], Momentum[pk]]) + 
         Pair[Momentum[pi], Momentum[pk]]*
          Pair[Momentum[pj], Momentum[pl]] - 
         Pair[Momentum[pi], Momentum[pj]]*Pair[Momentum[pk], Momentum[pl]]))
      /(2*Pair[Momentum[pi], Momentum[pl]]) + 
    (Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
        Spinor[Momentum[pj], 0, fq]*
       Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
        Spinor[Momentum[pk], 0, fq]*
       (3*Pair[Momentum[pi], Momentum[pl]]*
          Pair[Momentum[pj], Momentum[pk]] + 
         Pair[Momentum[pi], Momentum[pk]]*Pair[Momentum[pj], Momentum[pl]] - 
         Pair[Momentum[pi], Momentum[pj]]*Pair[Momentum[pk], Momentum[pl]]))/
     (2*Pair[Momentum[pi], Momentum[pl]])
             ) ];
sirlin2[ m_. *
  Spinor[Momentum[pi_], 0, fq___] . DiracGamma[Momentum[pk_]] . 
     Spinor[Momentum[pj_], 0, fq___]*
    Spinor[Momentum[pl_], 0, fq___] . DiracGamma[Momentum[pi_]] . 
     Spinor[Momentum[pk_], 0, fq___] 
       ] := Contract[ m ( 
   -((Spinor[Momentum[pi], 0, fq] . DiracGamma[Momentum[pl]] . 
          Spinor[Momentum[pj], 0, fq]*
         Spinor[Momentum[pl], 0, fq] . DiracGamma[Momentum[pj]] . 
          Spinor[Momentum[pk], 0, fq]*Pair[Momentum[pi], Momentum[pk]])/
       Pair[Momentum[pj], Momentum[pl]]) + 
    (Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
        Spinor[Momentum[pj], 0, fq]*
       Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
        Spinor[Momentum[pk], 0, fq]*
       (Pair[Momentum[pi], Momentum[pl]]*Pair[Momentum[pj], Momentum[pk]] + 
         3*Pair[Momentum[pi], Momentum[pk]]*
          Pair[Momentum[pj], Momentum[pl]] - 
         Pair[Momentum[pi], Momentum[pj]]*Pair[Momentum[pk], Momentum[pl]]))
      /(2*Pair[Momentum[pj], Momentum[pl]]) + 
    (Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
        DiracGamma[5] . Spinor[Momentum[pj], 0, fq]*
       Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
        DiracGamma[5] . Spinor[Momentum[pk], 0, fq]*
       (-(Pair[Momentum[pi], Momentum[pl]]*
            Pair[Momentum[pj], Momentum[pk]]) + 
         Pair[Momentum[pi], Momentum[pk]]*Pair[Momentum[pj], Momentum[pl]] + 
         Pair[Momentum[pi], Momentum[pj]]*Pair[Momentum[pk], Momentum[pl]]))/
     (2*Pair[Momentum[pj], Momentum[pl]])
               ) ] /; First[
  Spinor[Momentum[pi], 0, fq] . DiracGamma[Momentum[pk]].
    Spinor[Momentum[pj], 0, fq]*
         Spinor[Momentum[pl], 0, fq] . DiracGamma[Momentum[pi]] .
          Spinor[Momentum[pk], 0, fq]]===
    Spinor[Momentum[pi], 0, fq] . DiracGamma[Momentum[pk]].
    Spinor[Momentum[pj], 0, fq];

sirlin2[ m_. *
  Spinor[Momentum[pi_], 0, fq___] . DiracGamma[Momentum[pk_]] . DiracGamma[5] . 
     Spinor[Momentum[pj_], 0, fq___]*
    Spinor[Momentum[pl_], 0, fq___] . DiracGamma[Momentum[pj_]] . 
         DiracGamma[5] . 
     Spinor[Momentum[pk_], 0, fq___] 
       ] := Contract[ m (
   Spinor[Momentum[pi], 0, fq] . DiracGamma[Momentum[pk]] . 
      Spinor[Momentum[pj], 0, fq]*
     Spinor[Momentum[pl], 0, fq] . DiracGamma[Momentum[pj]] . 
      Spinor[Momentum[pk], 0, fq] - 
    Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      Spinor[Momentum[pj], 0, fq]*
     Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      Spinor[Momentum[pk], 0, fq]*Pair[Momentum[pj], Momentum[pk]] + 
    Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      DiracGamma[5] . Spinor[Momentum[pj], 0, fq]*
     Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      DiracGamma[5] . Spinor[Momentum[pk], 0, fq]*
     Pair[Momentum[pj], Momentum[pk]]
             )      ];

sirlin2[ m_. *
  Spinor[Momentum[pi_], 0, fq___] . DiracGamma[Momentum[pk_]] . DiracGamma[5] . 
     Spinor[Momentum[pj_], 0, fq___]*
    Spinor[Momentum[pl_], 0,fq___]. DiracGamma[Momentum[pi_]] . DiracGamma[5] . 
     Spinor[Momentum[pk_], 0, fq___] 
       ] :=  Contract[ m (
   -(Spinor[Momentum[pi], 0, fq] . DiracGamma[Momentum[pk]] . 
        Spinor[Momentum[pj], 0, fq]*
       Spinor[Momentum[pl], 0, fq] . DiracGamma[Momentum[pi]] . 
        Spinor[Momentum[pk], 0, fq]) + 
    Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      Spinor[Momentum[pj], 0, fq]*
     Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      Spinor[Momentum[pk], 0, fq]*Pair[Momentum[pi], Momentum[pk]] + 
    Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      DiracGamma[5] . Spinor[Momentum[pj], 0, fq]*
     Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      DiracGamma[5] . Spinor[Momentum[pk], 0, fq]*
     Pair[Momentum[pi], Momentum[pk]]
              ) ];

sirlin2[ m_. *
  Spinor[Momentum[pi_], 0, fq___] . DiracGamma[Momentum[pl_]] . DiracGamma[5] . 
     Spinor[Momentum[pj_], 0, fq___]*
    Spinor[Momentum[pl_], 0, fq___] . DiracGamma[Momentum[pj_]] . DiracGamma[5] . 
     Spinor[Momentum[pk_], 0, fq___] 
       ] := Contract[ m (
   -(Spinor[Momentum[pi], 0, fq] . DiracGamma[Momentum[pl]] . 
        Spinor[Momentum[pj], 0, fq]*
       Spinor[Momentum[pl], 0, fq] . DiracGamma[Momentum[pj]] . 
        Spinor[Momentum[pk], 0, fq]) + 
    Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      Spinor[Momentum[pj], 0, fq]*
     Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      Spinor[Momentum[pk], 0, fq]*Pair[Momentum[pj], Momentum[pl]] + 
    Spinor[Momentum[pi], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      DiracGamma[5] . Spinor[Momentum[pj], 0, fq]*
     Spinor[Momentum[pl], 0, fq] . DiracGamma[LorentzIndex[la]] . 
      DiracGamma[5] . Spinor[Momentum[pk], 0, fq]*
     Pair[Momentum[pj], Momentum[pl]]
              ) ];

 (* #################################################################### *)
 (*                             Main447                                  *)
 (* #################################################################### *)

dig[LorentzIndex[a_,___]]:=a;
dig[Momentum[a_,___]]:=a;
dig[x_]:=x/;(Head[x]=!=LorentzIndex)&&(Head[x]=!=Momentum);
dig[n_?NumberQ]:={};
getV[x_List]:=Select[Flatten[{x}/.Dot->List]/.DiracGamma -> dige ,
		     Head[#]===dige&]/.dige->dig;

(* Get a list of equal gamma matrices *)
schnitt[x___][y___]:=Intersection[ 
                        Select[Flatten[{x}/.Dot->List],!FreeQ[#,LorentzIndex]&],
                        Select[Flatten[{y}/.Dot->List],!FreeQ[#,LorentzIndex]&]
                               ];

(* get a list of not equal slashes and matrices *)
comp[x___][y___]:=Select[ Complement[Flatten[Union[{x},{y}]/.Dot->List],
                             schnitt[x][y] ],
                          !FreeQ2[#, {LorentzIndex, Momentum}]&
                        ];
                 
(* sirlin1def *)
(* do some ordering with sirlin1 ... *)
   Literal[
   sirlin1[m_. Spinor[p1__]. (gam1__) . Spinor[p2__] *
               Spinor[p3__]. (gam2__) . Spinor[p4__]
          ]] :=  memset[sirlin1[m Spinor[p1].gam1.Spinor[p2] * 
                               Spinor[p3].gam2.Spinor[p4]
                              ],
Block[{schnittmenge, compmenge, result,order, orderl,orderr},
                      schnittmenge = schnitt[gam1][gam2];
                       compmenge   = comp[gam1][gam2];
                        leftind    = comp[gam1][schnittmenge];
                        rightind   = comp[gam2][schnittmenge];
print3["entering sirlin1"];
(* We need at least two dummy indices for the sirlin relations *)
                 If[ Length[schnittmenge] > 1,

(* Test for eq. (12) *)
    If[(Length[schnittmenge] === 3) && (Length[compmenge] > 3),
       orderl = Join[ Drop[leftind, {1,2}], {schnittmenge[[1]], 
                      leftind[[1]], schnittmenge[[2]], 
                      leftind[[2]], schnittmenge[[3]]}
                    ] // getV;
       orderr = Join[ Drop[rightind, {1,2}], {schnittmenge[[1]],
                      rightind[[1]], schnittmenge[[2]],
                      rightind[[2]], schnittmenge[[3]]}
                    ] // getV;
       result = 
       Expand[m Contract[
                 DiracOrder[ Spinor[p1].gam1.Spinor[p2], orderl ]*
                 DiracOrder[ Spinor[p3].gam2.Spinor[p4], orderr ] ]
             ]//sirlin2
       ];
                 
 
(* ... *)
 (* Test for eq. (8) *)                    
    If[(Length[schnittmenge] === 2) && (Length[compmenge] > 1),
       order = Join[{First[schnittmenge]}, compmenge,
                    {Last[schnittmenge]} ] // getV;
       result = sirlin2[ Expand[ m  DiracOrder[
                         Spinor[p1].gam1.Spinor[p2] *
                         Spinor[p3].gam2.Spinor[p4], order]
                                                ]//Contract
                       ]
       ];
                ];
           If[!ValueQ[result], 
              result = sirlin2[m *
                         Spinor[p1].gam1.Spinor[p2] *
                         Spinor[p3].gam2.Spinor[p4]
                                     ]
             ];
print3["exiting sirlin1"];
           result]] /; !FreeQ[{gam1}, LorentzIndex];
                        

(* #################################################################### *)
(*                             Main45                                   *)
(* #################################################################### *)

   spcev0[x_] := spcev000[x]/.spcev000->spcev0ev;
(*
   spcev000[ a_ b_ ] := a spcev000[b] /; noncommQ[a] === True;
*)
   spcev000[y_] := y /; noncommQ[y] === True;
   spcev000[y_Times] := Select[ y, FreeQ[#, Spinor]& ] spcev0ev[
                       Select[ y,!FreeQ[#, Spinor]& ]          ];
   spcev0ev[x_] := pairexpand[coneins[
                     Expand2[spinlin[x], Spinor]/.Dot->spcevs/.
                                     spcev->Dot
                                    ]
                             ](*//Expand*);

   spcevs[xx___] := memset[ spcevs[xx], FixedPoint[ spcev,{xx},4 ] ];
(*spcevsdef*)

  (*spcevdef*)
   spcev[y_List]:=spcev@@y;
   spcev[a___,b_,c___] := b spcev[a,c] /; noncommQ[b] === True; 
   spcev[] = 1;
   Literal[ spcev[x___,Spinor[a__],y___] ]:=
     Expand[ DiracCanonical[ DiracEquation[fEx[DiracGammaExpand[
                                               x.Spinor[a].y]]/.dr->Dot
                                          ] ] ]/; FreeQ[{x,y},Spinor];
   Literal[ spcev[x___,Spinor[a__],b___,Spinor[c__],y___] ]:=
      Block[ {spcevdi,spcevre,spcevj},
        spcevdi = diracSimplify[Dot[Spinor[a],b,Spinor[c]],
                                     InsideDiracTrace->False,
(*
                                     DiracCanonical->True,
*)
                                     DiracCanonical->False,
                                     DiracInfo->False,
                                     Factoring->False,
                                     DiracSimpCombine->True
                               ];
        spcevdi = Expand[ pairexpand[ spcevdi ] ];
        spcevdi = Expand2[ spcevdi,Dot ];
        If[ !(Head[spcevdi]===Plus),
            spcevre = spinlin[ spcevdi ];
            spcevre = DiracEquation[ spcevre ];
            (*spcevre = DiracCanonical[ spcevre ]*),
            spcevre = Sum[(* DiracCanonical[*)
                           DiracEquation[ spinlin[ spcevdi[[spcevj]] ] ]
                                        (* ]*),
                           {spcevj,1,Length[spcevdi]}
                         ]
          ];
        spcevre = dotLin[spcevs[x].spcevre.spcevs[y]];
        If[ !FreeQ[spcevre, SUNT],
            spcevre = (spcevre/.Dot->dr)/.dr->Dot
          ];
        spcevre//dotLin] /; FreeQ[{b}, Spinor];

(* #################################################################### *)
(*                             Main46                                   *)
(* #################################################################### *)

(* ************************************************************** *)
(*                  Using the Dirac equation                      *)
(* ************************************************************** *)
                                                (*DiracEquationdef*)
 DiracEquation[x_]:=(*DiracEquation[x]=*)diraceq[x];
 
    last[n_. Momentum[pe__]]:=Momentum[pe];
    last[x_Plus]:=PowerExpand[Sqrt[Last[x]^2]];

   diraceq[x_]:=x/;FreeQ[x,Spinor];
   diraceq[x_] := Expand2[ pairexpand[ x//.spCDieqRules ],Dot ];
   spCDieqRules = {
    Literal[doot_[ z___,Spinor[n_. Momentum[p_] + k_. ,m_, op___],
         DiracGamma[Momentum[p_,___],___],a___
       ]] :>(m/n doot[ z,Spinor[n Momentum[p] + k,m,op ],a ] -
             If[(k===0), 0 ,
                If[last[n Momentum[p] + k] =!= Momentum[p],0,
                   1/n doot[ z, Spinor[n Momentum[p] + k,m,op ],
                             DiracGamma[k], a ] 
                  ]
               ]
             )/; last[n Momentum[p]+k]===Momentum[p],
 
    Literal[doot_[ a___,DiracGamma[Momentum[p_,___],___],
         Spinor[n_. Momentum[p_] + k_. ,m_,op___],z___
       ]]  :>(m/n doot[ a,Spinor[ n Momentum[p] + k,m,op ],z ] -
              If[(k===0), 0 ,
                If[last[n Momentum[p] + k] =!= Momentum[p],0,
                   1/n doot[ a, DiracGamma[k],
                                Spinor[n Momentum[p] + k,m,op ],
                             z ]
                  ]
                ]
              ) /; last[n Momentum[p]+k]===Momentum[p],
 
    Literal[doot_[ a___,DiracGamma[Momentum[y__],___],
         DiracGamma[Momentum[y__],___],b___
       ]] :> scev[Momentum[y],Momentum[y]] doot[a,b],
 
    Literal[doot_[ z___,Spinor[n_. Momentum[p_] + k_. ,m_,op___],a___,
            DiracGamma[x_[y__],di___],
         DiracGamma[Momentum[p_,dim___],dim___],b___]] :>
      ( - doot[ z,Spinor[n Momentum[p]+k,m,op ],a,
               DiracGamma[Momentum[p,dim],dim],
               DiracGamma[x[y],di],b
             ]
          + 2 coneins[ Pair[x[y],Momentum[p,dim] ]  *
                       doot[ z,Spinor[n Momentum[p]+k,m,op],a,b ]
                     ]
      ) /; last[n Momentum[p]+k] === Momentum[p],
 
    Literal[doot_[ a___,DiracGamma[Momentum[p_,___],___],DiracGamma[5],
         Spinor[n_. Momentum[p_] + k_. ,m_,op___],z___ ]] :>
         (-m/n doot[a,DiracGamma[5],Spinor[n Momentum[p]+k,m,op],z ]-
          If[k===0, 0,
             If[last[n Momentum[p] + k] =!= Momentum[p],0,
                   1/n doot[ a, DiracGamma[k], DiracGamma[5],
                             Spinor[n Momentum[p] + k,m,op ], z]
               ]
            ]
         ) /; last[n Momentum[p]+k]===Momentum[p],
 
    Literal[doot_[ a___,DiracGamma[Momentum[p_,___],___],DiracGamma[6],
         Spinor[n_. Momentum[p_] + k_. ,m_,op___],z___ ]] :>
         (m/n doot[a,DiracGamma[7],Spinor[n Momentum[p]+k,m,op],z ]-
          If[k===0, 0,
             If[last[n Momentum[p] + k] =!= Momentum[p],0,
                   1/n doot[ a, DiracGamma[k], DiracGamma[6],
                             Spinor[n Momentum[p] + k,m,op ], z]
               ]
            ]
         ) /; last[n Momentum[p]+k]===Momentum[p],
 
    Literal[doot_[ a___,DiracGamma[Momentum[p_,___],___],DiracGamma[7],
         Spinor[n_. Momentum[p_] + k_. ,m_,op___],z___ ]] :>
         (m/n doot[a,DiracGamma[6],Spinor[n Momentum[p]+k,m,op],z ]-
          If[k===0, 0,
             If[last[n Momentum[p] + k] =!= Momentum[p],0,
                   1/n doot[ a, DiracGamma[k], DiracGamma[7],
                             Spinor[n Momentum[p] + k,m,op ], z]
               ]
            ]
         ) /; last[n Momentum[p]+k]===Momentum[p],
 
    Literal[doot_[ a___,DiracGamma[ Momentum[p_,dim___],dim___],
              DiracGamma[x_[y__],di___],b___,
              Spinor[n_. Momentum[p_] + k_. ,m_,op___],z___
       ]] :> (- doot[ a,DiracGamma[x[y],di],
                      DiracGamma[Momentum[p,dim],dim],b,
                    Spinor[n Momentum[p] + k,m,op ],z
                  ]
          + 2 coneins[ Pair[x[y],Momentum[p,dim]] *
                       doot[ a,b,Spinor[n Momentum[p] +k,m,op],z ]
                     ]
             ) /; last[n Momentum[p]+k]===Momentum[p]
            };


(* #################################################################### *)
(*                             Main47                                   *)
(* #################################################################### *)

(* ************************************************************** *)
(* Traces, and all that                                           *)
(* ************************************************************** *)
(* Trdef *)

 Options[ Trc ] = {Factoring->False,Mandelstam->{},
                   PairCollect->False,
                          DiracTraceEvaluate->True,
                          LeviCivitaSign -> (-1),
                         Schouten -> 442,
                       TraceOfOne -> 4
                  };
Trc[x__, rul___Rule] := DiracTrace @@ Join[{x}, Join[{rul},Options[Trc]]]; 

DiracTrace[0,___]:=0;
(*EvaluateDiracTracedef*)
 EvaluateDiracTrace[x_]:= x/.DiracTrace->Trc;
 Options[ DiracTrace ] = {Factoring->False,Mandelstam->{},
                          PairCollect->True,
                          DiracTraceEvaluate->False,
                          Schouten -> 442 (*,
                          LeviCivitaSign -> (-1) *)
                         };
 DiracTrace[a_ /; (FreeQ[a, DiracGamma] && !FreeQ[a, DiracGammaT]),b___] := 
    DiracTrace[(a//Transpose)//Reverse, b];
 DiracTrace[a___, x_,y_, z___]:=DiracTrace[a,x.y,z]/;FreeQ[y,Rule]&&
                                                     FreeQ[x,Rule];
 
                                                 (*DiracTracedef*)
 DiracTrace[x_,op___] := (diractraceev[dirspec[x] /.
                                       dirspec -> Identity,
                                  {op}] /. diractraceev->diractraceev2)/;
  ( DiracTraceEvaluate/.{op} /. (Join[{op},Options[DiracTrace]]//Flatten))===
     True;
(*  special case *)

dirli[LorentzIndex[xx_, ___],___] := xx;

dirspec[x_Dot] := Apply[dirspec2, x]/.dirspec2 -> Dot;
dirspec2[DiracGamma[LorentzIndex[a1_,dii_],dii_],
             DiracGamma[LorentzIndex[a2_,dii_],dii_],
             DiracGamma[LorentzIndex[a3_,dii_],dii_],
             a4:DiracGamma[LorentzIndex[_,dii_],dii_]..,
             DiracGamma[LorentzIndex[a1_,dii_],dii_],
             DiracGamma[LorentzIndex[a2_,dii_],dii_],
             DiracGamma[LorentzIndex[a3_,dii_],dii_],
             a4:DiracGamma[LorentzIndex[_,dii_],dii_]..
            ]:= dcs[dii]@@Join[{a1,a2,a3}, {a4}/.Diracgamma->dirli,
                              {a1,a2,a3}, {a4}/.Diracgamma->dirli
                             ];

dcs[dim_][x___] := dcs[dim][x] = (dics[dim][x] /. dics->dc);
dc[_][]=1; dics[_][]=1;
dics[dI_][a___, n_, n_, b___] := dI * dics[dI][a, b];
dics[dI_][a___, n_, z_, n_, b___ ] := (2-dI) * dics[dI][a, z, b];
dics[dI_][a___, n_, v_, w_, n_, b___
        ] := (dI-4) * dics[dI][a, v,w, b] + 4 (dics[dI]@@({a, b}/. v -> w));
dics[dI_][a___, n_, v_, w_, z_, n_, b___
        ] := (4-dI) * dics[dI][a, v,w,z, b] - 2 dics[dI][a, z,w,v,b];
dics[dI_][a___, n_, mu_, nu_, ro_,si_, n_, b___
        ] := (dI-4) * dics[dI][a, mu,nu,ro,si, b] +
             2 dics[dI][a, ro,nu,mu,si,b] + 2 dics[dI][a, si,mu,nu,ro,b];
dics[dI_][a___, n_, mu_, nu_, ro_, si_, de_, n_, b___
        ] := (4-dI) * dics[dI][a, mu,nu,ro,si,de, b] -
                 2 dics[dI][a, mu,de,nu,ro,si, b] -
                 2 dics[dI][a, mu,si,ro,nu,de, b] +
                 2 dics[dI][a, nu,ro,si,de,mu, b];
dicsav[dd_][x___] := dicsav[dd][x] = dics[dd][x];
dc[di_][a___, mu_, lim__, mu_, b___] :=
Expand[
Block[{m = Length[{lim}], i, j},
      (-1)^m ( (di-2 m) dicss[di][a,lim,b] -
      4 Sum[(-1)^(j-i) If[{lim}[[j]] === {lim}[[i]],
                           di (dicss[di] @@
                               Join[{a}, Delete[{lim}, {{i},{j}}], {b}]
                              ),
                              dicss[di] @@
                               (Join[{a}, Delete[{lim}, {{i},{j}}], {b}]
                                /. ({lim}[[j]]) -> ({lim}[[i]]))
                         ],
            {i,1,m-1}, {j,i+1,m}])
     ] /. dicss -> dicsav//. dics -> dcs];
(* ****************************************************** *)

                                                  (*diractraceevdef*)
 diractraceev[x_, opt___] := Block[{trfa = 1, enx = x},
   If[ Head[x] === Times,
       trfa = Select[x, FreeQ2[#, {DiracGamma, LorentzIndex, SUNIndex, Eps}]&];
       enx = x / trfa;
     ];
   If[!FreeQ[x, SUNT],
      enx = SUNSimplify[DiracTrace[enx,opt]] /.DiracTrace -> diractraceev2;
     ];
           diractraceev2[enx] trfa];

 diractraceev2[x_,opt_:{}]:= 
     ( TraceOfOne /. opt /.Options[Trc] /. Options[DiracTrace] ) * x /; 
         FreeQ[x,DiracGamma];
 
(* If no Dot's  but DiracGamma's are present *)
 diractraceev2[y_,opt_:{}]:=Block[
                              {diractrpa,diractrtemp,diractrresu,four},
  four = TraceOfOne/.opt /. Options[Trc] /. Options[DiracTrace];
  diractrtemp = Expand2[ conall[ y ]//gamma67back, DiracGamma ];
  If[ Head[diractrtemp]===Plus,
      diractrresu = Map[ Apply[Tr,Prepend[opt,#]]&,diractrtemp],
        diractrpa = PartitHead[ diractrtemp,DiracGamma ];
        diractrresu = diractrpa[[1]] four spursav[ diractrpa[[2]] ]
    ];
  diractrresu = Expand[diractrresu];
                  diractrresu] /;( FreeQ[y,Dot] && !FreeQ[y,DiracGamma]);
 
 
(* #################################################################### *)
(*                             Main48                                   *)
(* #################################################################### *)

 diractraceev2[nnx_,in_:{}]:= Block[{diractrjj,diractrlnx,diractrres,
                                    diractrny=0,mand,diractrfact,nx,
                                    diractrcoll,traceofone,schoutenopt,
                                    noepscontractlabel},
   opt = Join[ Flatten[{in}],Options[Trc], Options[DiracTrace] ];
   mand=Mandelstam/.opt;
   diractrfact=Factoring/.opt;
   diractrcoll=PairCollect/.opt;
   schoutenopt = Schouten /. opt;
   traceofone = TraceOfOne /.  opt;
   nx = Collect2[conall[nnx], Dot, Factoring -> False];
   If[ Head[nx]===Plus && Length[nx] > 100,
       diractrlnx = Length[nx]; diractrjj = 0;
       While[ diractrjj<diractrlnx,diractrjj++;
              diractrny = diractrny + Expand2[
                          ( Expand2[ diracSimplify[ nx[[diractrjj]],
                                     {InsideDiracTrace->True,
                                      Factoring->False,
                                      DiracCanonicalFlag->False}
                                                  ], Dot
                                   ]/.Dot->spursav/.Pair->scev/.
                                      Pair->sCO/.sCO->scev
                          )/.DiracGamma[5]->0/.
                             DiracGamma[6]->(1/2)/.DiracGamma[7]->(1/2),
                               DiracGamma
                              ];
            ],
noepscontractlabel= False;
If[Length[nx]>6,
   If[Length[Union[List@@nx]] === Length[nx],
      noepscontractlabel = True
     ]
  ];
(*
         nx = nx /. Dot -> spursav;
         If[!FreeQ[nx, spursav], nx = nx /. spursav -> Dot,
             nx = nx /. Pair -> sCO /. sCO  -> scev
           ];
*)
        If[FreeQ[nx,DiracGamma],
           diractrny = nx,
             diractrny = Expand[ (Expand2[ diracSimplify[ nx,
                                    {InsideDiracTrace->True,
                                     Factoring->False,
                                     DiracCanonicalFlag->False}
                                                        ],Dot
                                  ]/.Dot->spursav/.
                          Pair->scev/.Pair->sCO/.sCO->scev
                          )/.DiracGamma[5]->0/.DiracGamma[6]->(1/2)/.
                                               DiracGamma[7]->(1/2)
                        ];
         ];
     ];
 
   If[!FreeQ[diractrny, Eps], 
      If[Head[diractrny] === Plus,
         If[noepscontractlabel =!= True,
            diractrny = Contract[ diractrny, 
                        EpsContract -> True, Schouten->schoutenopt];
           ]
        ]
     ];
   If[ diractrfact===True, diractrres = Factor2[traceofone diractrny],
(* this  2@#$^^$#%^@*#$ mma!!!!; 
                           diractrres = Expand[ traceofone diractrny ]
*)
                           diractrres = traceofone diractrny
     ];
   If[ Length[ mand ] >0,
       diractrres = TrickMandelstam @@ Prepend[ {mand}, diractrres ]
     ];

   diractrpc[x__]:=Plus[x]/;FreeQ[{x},Pair];
   If[ diractrcoll===True,
   diractrpc[x__]:=collone[ Plus[x],Pair ];
       diractrres = diractrres/.Plus->diractrpc ];
                      diractrres]/;!FreeQ2[nnx,{Dot,DiracGamma}];
(* endof diractraceev1 *)
(* ************************************************************** *)
 
(* #################################################################### *)
(*                             Main49                                   *)
(* #################################################################### *)

 
spursavg[x___, LorentzIndex[a_, de___], LorentzIndex[a_, de___], y___] :=
  (de spursavg[x, y]) /. spursavg -> spug;
diracga[LorentzIndex[mu_, dii_]] := diracga[LorentzIndex[mu,dii],dii];
diracga[Momentum[p_, dii_]] := diracga[Momentum[p, dii],dii];
spug[x___] := spursav@@(Map[diracga, {x}] /. diracga -> DiracGamma);

(* calculation of traces (recursively) --  up to a factor of 4 *)
   spursav[x_DiracGamma,y_DiracGamma,r_DiracGamma,z_DiracGamma,
           DiracGamma[5]]:=Block[{dirsign},
        dirsign = LeviCivitaSign /. Options[Trc];
        dirsign I Apply[ Eps, {x,y,r,z}/.
                              DiracGamma[vl_[mp_,di___],di___]->vl[mp,di]
                 ]//EpsEvaluate];
 
   spursav[x___DiracGamma]:=memset[ spursav[x], spur[x] ];
 
   
   spur[]=1;
   spur[DiracGamma[5]]=0;
   spur[x_[y__],DiracGamma[5]]:=0;
   spur[a_[b__],x_[y__],DiracGamma[5]]:=0;
   spur[a_[b__],c_[d__],x_[y__],DiracGamma[5]]:=0;
   spur[a__] := (spur @@ Reverse[Transpose[{a}]]) /; 
                (!FreeQ[{a}, DiracGammaT]) && FreeQ[{a},DiracGamma];

(* This is a definition of   Trace( 1.2.3.4. gamma[5] ) *)
   spur[x_,y_,r_,z_,DiracGamma[5]]:=Block[{dirsign},
        dirsign = LeviCivitaSign /. Options[Trc];
        dirsign I Apply[Eps, {x,y,r,z}/.DiracGamma[vl_[mp_,dii___],___
                                                   ]->vl[mp,dii]
                       ]//EpsEvaluate];
 
   spur[m_,n_,r_,s_,l_,t_,DiracGamma[5]]:= Block[{dirsign, sres, ltr},
     If[($Kreimer === True) && (!OrderedQ[{m,n,r,s,l,t}]),
           Trc[1/(TraceOfOne/.Options[Trc]) DiracOrder[ m.n.r.s.l.t.DiracGamma[5] ] ],
        If[$Larin === True && !FreeQ[{m,n,r,s,l,t}, DiracGamma[LorentzIndex[_,_],_]],
           ltr[a1_, a2_, a3_, a4_, a5_][
                 DiracGamma[LorentzIndex[in_,di___], di___]
                                       ] :=  
    Block[{f1, f2, f3,drsi}, 
          drsi = LeviCivitaSign /. Options[Trc];
          drsi = drsi/(TraceOfOne/.Options[Trc]);
(*drsi is usually -1/4 *)
          {f1, f2, f3} = LorentzIndex[#, D]& /@ Unique[{"L","L","L"}];
          Trc[drsi I/6 Eps[LorentzIndex[in, di], f1, f2, f3] *
             a1.a2.a3.a4.a5.DiracGamma[f1, D] . DiracGamma[f2, D] .
                            DiracGamma[f3, D]
            ]
         ];
           Which[ MatchQ[t, DiracGamma[ LorentzIndex[__], ___]], 
                  ltr[m,n,r,s,l][t],
                  MatchQ[l, DiracGamma[ LorentzIndex[__], ___]], 
                  -ltr[m,n,r,s,t][l],
                  MatchQ[s, DiracGamma[ LorentzIndex[__], ___]], 
                  ltr[m,n,r,t,l][s],
                  MatchQ[r, DiracGamma[ LorentzIndex[__], ___]], 
                  -ltr[m,n,s,t,l][r],
                  MatchQ[n, DiracGamma[ LorentzIndex[__], ___]], 
                  ltr[m,r,s,t,l][n],
                  MatchQ[m, DiracGamma[ LorentzIndex[__], ___]], 
                  -ltr[n,r,s,t,l][m]
                ],
        dirsign = LeviCivitaSign /. Options[Trc];
       Expand[ + dirsign I (
        scev[ m//gc,n//gc ]  Apply[ Eps, {l,r,s,t}//gc ] -
        scev[ m//gc,r//gc ]  Apply[ Eps, {l,n,s,t}//gc ] +
        scev[ n//gc,r//gc ]  Apply[ Eps, {l,m,s,t}//gc ] +
        scev[ s//gc,l//gc ]  Apply[ Eps, {m,n,r,t}//gc ] +
        scev[ l//gc,t//gc ]  Apply[ Eps, {m,n,r,s}//gc ] +
        scev[ s//gc,t//gc ]  Apply[ Eps, {l,m,n,r}//gc ]
                                                       )//EpsEvaluate
                                               ] ] ] ];      (*spurdef*)

(* this trace is calculated via expressing  DiracMatrix[w1,w2,w3]
   by the Chisholm - identity;
   thus it is only valid in four dimensions and in the naive
   gamma5 prescription
*)
spur[w1_,w2_,w3_,w4_,w5_,w6_,w7_,w8_,DiracGamma[5]
    ]:= Block[{trsign,z1,z2,z3,z4,z5,z6,z7,z8},
{z1,z2,z3,z4,z5,z6,z7,z8} = 
{w1,w2,w3,w4,w5,w6,w7,w8} /.DiracGamma[vl_[mp_,dii___],___]->vl[mp,dii];
trsign = LeviCivitaSign /. Options[Trc];
(* trsign is usually  =  -1 *)
(* factor 4 is put later *)
trsign*I*(Eps[z5, z6, z7, z8]*Pair[z1, z4]*Pair[z2, z3] -
     Eps[z4, z6, z7, z8]*Pair[z1, z5]*Pair[z2, z3] -
     Eps[z5, z6, z7, z8]*Pair[z1, z3]*Pair[z2, z4] +
     Eps[z4, z6, z7, z8]*Pair[z1, z3]*Pair[z2, z5] +
     Eps[z5, z6, z7, z8]*Pair[z1, z2]*Pair[z3, z4] -
     Eps[z4, z6, z7, z8]*Pair[z1, z2]*Pair[z3, z5] +
     Eps[z3, z6, z7, z8]*Pair[z1, z2]*Pair[z4, z5] -
     Eps[z2, z6, z7, z8]*Pair[z1, z3]*Pair[z4, z5] +
     Eps[z1, z6, z7, z8]*Pair[z2, z3]*Pair[z4, z5] +
     Eps[z1, z2, z3, z8]*Pair[z4, z7]*Pair[z5, z6] -
     Eps[z1, z2, z3, z7]*Pair[z4, z8]*Pair[z5, z6] -
     Eps[z1, z2, z3, z8]*Pair[z4, z6]*Pair[z5, z7] +
     Eps[z1, z2, z3, z6]*Pair[z4, z8]*Pair[z5, z7] +
     Eps[z1, z2, z3, z7]*Pair[z4, z6]*Pair[z5, z8] -
     Eps[z1, z2, z3, z6]*Pair[z4, z7]*Pair[z5, z8] +
     Eps[z3, z4, z5, z8]*Pair[z1, z2]*Pair[z6, z7] -
     Eps[z2, z4, z5, z8]*Pair[z1, z3]*Pair[z6, z7] +
     Eps[z1, z4, z5, z8]*Pair[z2, z3]*Pair[z6, z7] +
     Eps[z1, z2, z3, z8]*Pair[z4, z5]*Pair[z6, z7] -
     Eps[z1, z2, z3, z5]*Pair[z4, z8]*Pair[z6, z7] +
     Eps[z1, z2, z3, z4]*Pair[z5, z8]*Pair[z6, z7] -
     Eps[z3, z4, z5, z7]*Pair[z1, z2]*Pair[z6, z8] +
     Eps[z2, z4, z5, z7]*Pair[z1, z3]*Pair[z6, z8] -
     Eps[z1, z4, z5, z7]*Pair[z2, z3]*Pair[z6, z8] -
     Eps[z1, z2, z3, z7]*Pair[z4, z5]*Pair[z6, z8] +
     Eps[z1, z2, z3, z5]*Pair[z4, z7]*Pair[z6, z8] -
     Eps[z1, z2, z3, z4]*Pair[z5, z7]*Pair[z6, z8] +
     Eps[z3, z4, z5, z6]*Pair[z1, z2]*Pair[z7, z8] -
     Eps[z2, z4, z5, z6]*Pair[z1, z3]*Pair[z7, z8] +
     Eps[z1, z4, z5, z6]*Pair[z2, z3]*Pair[z7, z8] +
     Eps[z1, z2, z3, z6]*Pair[z4, z5]*Pair[z7, z8] -
     Eps[z1, z2, z3, z5]*Pair[z4, z6]*Pair[z7, z8] +
     Eps[z1, z2, z3, z4]*Pair[z5, z6]*Pair[z7, z8])
              ] /; $Larin =!= True;

(* this trace has been calculated according to Larin,
   i.e. expression DiracMatrix[w8].DiracGamma[5] by
   (-I/6) LeviCivita[w8,mu,nu,la] DiracMatrix[mu,nu,la]
*)
spur[w1_,w2_,w3_,w4_,w5_,w6_,w7_,w8_,DiracGamma[5]
    ]:= Block[{trsign,z1,z2,z3,z4,z5,z6,z7,z8},
{z1,z2,z3,z4,z5,z6,z7,z8} = 
{w1,w2,w3,w4,w5,w6,w7,w8} /.DiracGamma[vl_[mp_,dii___],___]->vl[mp,dii];
trsign = LeviCivitaSign /. Options[Trc];
(* trsign is usually  =  -1 *)
(* factor 4 is put later *)
trsign*I*(Eps[z5, z6, z7, z8]*Pair[z1, z4]*Pair[z2, z3] - 
    Eps[z4, z6, z7, z8]*Pair[z1, z5]*Pair[z2, z3] + 
    Eps[z4, z5, z7, z8]*Pair[z1, z6]*Pair[z2, z3] - 
    Eps[z4, z5, z6, z8]*Pair[z1, z7]*Pair[z2, z3] - 
    Eps[z5, z6, z7, z8]*Pair[z1, z3]*Pair[z2, z4] + 
    Eps[z3, z6, z7, z8]*Pair[z1, z5]*Pair[z2, z4] - 
    Eps[z3, z5, z7, z8]*Pair[z1, z6]*Pair[z2, z4] + 
    Eps[z3, z5, z6, z8]*Pair[z1, z7]*Pair[z2, z4] + 
    Eps[z4, z6, z7, z8]*Pair[z1, z3]*Pair[z2, z5] - 
    Eps[z3, z6, z7, z8]*Pair[z1, z4]*Pair[z2, z5] + 
    Eps[z3, z4, z7, z8]*Pair[z1, z6]*Pair[z2, z5] - 
    Eps[z3, z4, z6, z8]*Pair[z1, z7]*Pair[z2, z5] - 
    Eps[z4, z5, z7, z8]*Pair[z1, z3]*Pair[z2, z6] + 
    Eps[z3, z5, z7, z8]*Pair[z1, z4]*Pair[z2, z6] - 
    Eps[z3, z4, z7, z8]*Pair[z1, z5]*Pair[z2, z6] + 
    Eps[z3, z4, z5, z8]*Pair[z1, z7]*Pair[z2, z6] + 
    Eps[z4, z5, z6, z8]*Pair[z1, z3]*Pair[z2, z7] - 
    Eps[z3, z5, z6, z8]*Pair[z1, z4]*Pair[z2, z7] + 
    Eps[z3, z4, z6, z8]*Pair[z1, z5]*Pair[z2, z7] - 
    Eps[z3, z4, z5, z8]*Pair[z1, z6]*Pair[z2, z7] + 
    Eps[z5, z6, z7, z8]*Pair[z1, z2]*Pair[z3, z4] - 
    Eps[z2, z6, z7, z8]*Pair[z1, z5]*Pair[z3, z4] + 
    Eps[z2, z5, z7, z8]*Pair[z1, z6]*Pair[z3, z4] - 
    Eps[z2, z5, z6, z8]*Pair[z1, z7]*Pair[z3, z4] + 
    Eps[z1, z6, z7, z8]*Pair[z2, z5]*Pair[z3, z4] - 
    Eps[z1, z5, z7, z8]*Pair[z2, z6]*Pair[z3, z4] + 
    Eps[z1, z5, z6, z8]*Pair[z2, z7]*Pair[z3, z4] - 
    Eps[z4, z6, z7, z8]*Pair[z1, z2]*Pair[z3, z5] + 
    Eps[z2, z6, z7, z8]*Pair[z1, z4]*Pair[z3, z5] - 
    Eps[z2, z4, z7, z8]*Pair[z1, z6]*Pair[z3, z5] + 
    Eps[z2, z4, z6, z8]*Pair[z1, z7]*Pair[z3, z5] - 
    Eps[z1, z6, z7, z8]*Pair[z2, z4]*Pair[z3, z5] + 
    Eps[z1, z4, z7, z8]*Pair[z2, z6]*Pair[z3, z5] - 
    Eps[z1, z4, z6, z8]*Pair[z2, z7]*Pair[z3, z5] + 
    Eps[z4, z5, z7, z8]*Pair[z1, z2]*Pair[z3, z6] - 
    Eps[z2, z5, z7, z8]*Pair[z1, z4]*Pair[z3, z6] + 
    Eps[z2, z4, z7, z8]*Pair[z1, z5]*Pair[z3, z6] - 
    Eps[z2, z4, z5, z8]*Pair[z1, z7]*Pair[z3, z6] + 
    Eps[z1, z5, z7, z8]*Pair[z2, z4]*Pair[z3, z6] - 
    Eps[z1, z4, z7, z8]*Pair[z2, z5]*Pair[z3, z6] + 
    Eps[z1, z4, z5, z8]*Pair[z2, z7]*Pair[z3, z6] - 
    Eps[z4, z5, z6, z8]*Pair[z1, z2]*Pair[z3, z7] + 
    Eps[z2, z5, z6, z8]*Pair[z1, z4]*Pair[z3, z7] - 
    Eps[z2, z4, z6, z8]*Pair[z1, z5]*Pair[z3, z7] + 
    Eps[z2, z4, z5, z8]*Pair[z1, z6]*Pair[z3, z7] - 
    Eps[z1, z5, z6, z8]*Pair[z2, z4]*Pair[z3, z7] + 
    Eps[z1, z4, z6, z8]*Pair[z2, z5]*Pair[z3, z7] - 
    Eps[z1, z4, z5, z8]*Pair[z2, z6]*Pair[z3, z7] + 
    Eps[z3, z6, z7, z8]*Pair[z1, z2]*Pair[z4, z5] - 
    Eps[z2, z6, z7, z8]*Pair[z1, z3]*Pair[z4, z5] + 
    Eps[z2, z3, z7, z8]*Pair[z1, z6]*Pair[z4, z5] - 
    Eps[z2, z3, z6, z8]*Pair[z1, z7]*Pair[z4, z5] + 
    Eps[z1, z6, z7, z8]*Pair[z2, z3]*Pair[z4, z5] - 
    Eps[z1, z3, z7, z8]*Pair[z2, z6]*Pair[z4, z5] + 
    Eps[z1, z3, z6, z8]*Pair[z2, z7]*Pair[z4, z5] + 
    Eps[z1, z2, z7, z8]*Pair[z3, z6]*Pair[z4, z5] - 
    Eps[z1, z2, z6, z8]*Pair[z3, z7]*Pair[z4, z5] - 
    Eps[z3, z5, z7, z8]*Pair[z1, z2]*Pair[z4, z6] + 
    Eps[z2, z5, z7, z8]*Pair[z1, z3]*Pair[z4, z6] - 
    Eps[z2, z3, z7, z8]*Pair[z1, z5]*Pair[z4, z6] + 
    Eps[z2, z3, z5, z8]*Pair[z1, z7]*Pair[z4, z6] - 
    Eps[z1, z5, z7, z8]*Pair[z2, z3]*Pair[z4, z6] + 
    Eps[z1, z3, z7, z8]*Pair[z2, z5]*Pair[z4, z6] - 
    Eps[z1, z3, z5, z8]*Pair[z2, z7]*Pair[z4, z6] - 
    Eps[z1, z2, z7, z8]*Pair[z3, z5]*Pair[z4, z6] + 
    Eps[z1, z2, z5, z8]*Pair[z3, z7]*Pair[z4, z6] + 
    Eps[z3, z5, z6, z8]*Pair[z1, z2]*Pair[z4, z7] - 
    Eps[z2, z5, z6, z8]*Pair[z1, z3]*Pair[z4, z7] + 
    Eps[z2, z3, z6, z8]*Pair[z1, z5]*Pair[z4, z7] - 
    Eps[z2, z3, z5, z8]*Pair[z1, z6]*Pair[z4, z7] + 
    Eps[z1, z5, z6, z8]*Pair[z2, z3]*Pair[z4, z7] - 
    Eps[z1, z3, z6, z8]*Pair[z2, z5]*Pair[z4, z7] + 
    Eps[z1, z3, z5, z8]*Pair[z2, z6]*Pair[z4, z7] + 
    Eps[z1, z2, z6, z8]*Pair[z3, z5]*Pair[z4, z7] - 
    Eps[z1, z2, z5, z8]*Pair[z3, z6]*Pair[z4, z7] + 
    Eps[z3, z4, z7, z8]*Pair[z1, z2]*Pair[z5, z6] - 
    Eps[z2, z4, z7, z8]*Pair[z1, z3]*Pair[z5, z6] + 
    Eps[z2, z3, z7, z8]*Pair[z1, z4]*Pair[z5, z6] - 
    Eps[z2, z3, z4, z8]*Pair[z1, z7]*Pair[z5, z6] + 
    Eps[z1, z4, z7, z8]*Pair[z2, z3]*Pair[z5, z6] - 
    Eps[z1, z3, z7, z8]*Pair[z2, z4]*Pair[z5, z6] + 
    Eps[z1, z3, z4, z8]*Pair[z2, z7]*Pair[z5, z6] + 
    Eps[z1, z2, z7, z8]*Pair[z3, z4]*Pair[z5, z6] - 
    Eps[z1, z2, z4, z8]*Pair[z3, z7]*Pair[z5, z6] + 
    Eps[z1, z2, z3, z8]*Pair[z4, z7]*Pair[z5, z6] - 
    Eps[z3, z4, z6, z8]*Pair[z1, z2]*Pair[z5, z7] + 
    Eps[z2, z4, z6, z8]*Pair[z1, z3]*Pair[z5, z7] - 
    Eps[z2, z3, z6, z8]*Pair[z1, z4]*Pair[z5, z7] + 
    Eps[z2, z3, z4, z8]*Pair[z1, z6]*Pair[z5, z7] - 
    Eps[z1, z4, z6, z8]*Pair[z2, z3]*Pair[z5, z7] + 
    Eps[z1, z3, z6, z8]*Pair[z2, z4]*Pair[z5, z7] - 
    Eps[z1, z3, z4, z8]*Pair[z2, z6]*Pair[z5, z7] - 
    Eps[z1, z2, z6, z8]*Pair[z3, z4]*Pair[z5, z7] + 
    Eps[z1, z2, z4, z8]*Pair[z3, z6]*Pair[z5, z7] - 
    Eps[z1, z2, z3, z8]*Pair[z4, z6]*Pair[z5, z7] + 
    Eps[z3, z4, z5, z8]*Pair[z1, z2]*Pair[z6, z7] - 
    Eps[z2, z4, z5, z8]*Pair[z1, z3]*Pair[z6, z7] + 
    Eps[z2, z3, z5, z8]*Pair[z1, z4]*Pair[z6, z7] - 
    Eps[z2, z3, z4, z8]*Pair[z1, z5]*Pair[z6, z7] + 
    Eps[z1, z4, z5, z8]*Pair[z2, z3]*Pair[z6, z7] - 
    Eps[z1, z3, z5, z8]*Pair[z2, z4]*Pair[z6, z7] + 
    Eps[z1, z3, z4, z8]*Pair[z2, z5]*Pair[z6, z7] + 
    Eps[z1, z2, z5, z8]*Pair[z3, z4]*Pair[z6, z7] - 
    Eps[z1, z2, z4, z8]*Pair[z3, z5]*Pair[z6, z7] + 
    Eps[z1, z2, z3, z8]*Pair[z4, z5]*Pair[z6, z7])
] /; $Larin === True;

   spur[x__,DiracGamma[6]]:=1/2 spur[x] + 1/2 spur[x,DiracGamma[5]];
   spur[x__,DiracGamma[7]]:=1/2 spur[x] - 1/2 spur[x,DiracGamma[5]];
 
 
   spur[x__]:=( DiracTrace@@ ( gamma67backj[ {x} ] )
              ) /; !FreeQ2[{x},{DiracGamma[6],DiracGamma[7]}];
 
(*
   gc[x_]:=ToFourDimensions[ x/.DiracGamma->gach ];         (*gcdef*)
*)
   gc[x_]:=x/.DiracGamma->gach;
   gach[ex_,___]:=ex /; Length[ex]>0;                     (*gachdef*)
   gach[n_Integer]=DiracGamma[n];
 
   spur[y__] :=Block[ {spx,le=Length[{y}],tempres,i,spurjj,tempr,
                       temp2 = 0,fi,spt, resp,scx,dirsign},
                spx = ( {y}//DiracGammaExpand )/.DiracGamma->gach;
                scx[a_,b_]:=scev[spx[[a]],spx[[b]]];
 
                resp =
     Catch[
    If[ OddQ[le] && fr567[spx], Throw[0] ];
    If[ le===2,Throw[scev[spx[[1]],spx[[2]]]/.Pair->sCO/.sCO->Pair] ];
    If[ le===4,
        Throw[(scx[1,2] scx[3,4]-scx[1,3] scx[2,4]+scx[1,4] scx[2,3]
              )//Expand
             ]
      ];
     If[ le===6,
         Throw[(
          scx[1,6] scx[2,5] scx[3,4] - scx[1,5] scx[2,6] scx[3,4] -
          scx[1,6] scx[2,4] scx[3,5] + scx[1,4] scx[2,6] scx[3,5] +
          scx[1,5] scx[2,4] scx[3,6] - scx[1,4] scx[2,5] scx[3,6] +
          scx[1,6] scx[2,3] scx[4,5] - scx[1,3] scx[2,6] scx[4,5] +
          scx[1,2] scx[3,6] scx[4,5] - scx[1,5] scx[2,3] scx[4,6] +
          scx[1,3] scx[2,5] scx[4,6] - scx[1,2] scx[3,5] scx[4,6] +
          scx[1,4] scx[2,3] scx[5,6] - scx[1,3] scx[2,4] scx[5,6] +
          scx[1,2] scx[3,4] scx[5,6]
                )//Expand
              ]
       ];
 
   If[FreeQ[spx,DiracGamma[5]],
        For[i=2, i<le+1, i++,
            temp2 += ((-1)^i) coneins[
                     scev[spx[[1]],spx[[i]]] spt@@Rest[Drop[spx,{i,i}]]
                                      ]
           ]; Throw[Expand[ temp2/.spt->spursavga/.spursavga->spug] ],
 
    If[Global`$BreitMaison===True,
        dirsign = LeviCivitaSign /. Options[Trc];
    fi = Table[LorentzIndex[ Unique[] ],{spurjj,1,4}];
Print["check5"];
    Throw[ DiracTrace @@
           ( {spx}/.DiracGamma[5]->
             (dirsign I/24 (DiracGamma[fi[[1]]].DiracGamma[fi[[2]]].
                    DiracGamma[fi[[3]]].DiracGamma[fi[[4]]]
                   ) (Eps@@fi)
             )
           )
         ]
      ,
       If[$Kreimer === True, NochNichtFertig,
          If[$Larin === True,
             {fi1, fi2, fi3} = LorentzIndex[#,D]& /@ Unique[{"a","b","c"}];
              drsi = LeviCivitaSign /. Options[Trc];
              drsi = drsi/(TraceOfOne/.Options[Trc]);
             (*drsi is usually -1/4 *)
             temp2 = spx /. {a___, lomo_[mUU_,di___], DiracGamma[5]} :>
                     Trc[ drsi I/6 Eps[lomo[mUU,di], fi1, fi2, fi3] *
                         Dot @@ Map[DiracGamma[#,D]&, {a,fi1,fi2,fi3}]];
         ]  ];

    If[$Larin === False && $Kreimer === False,
    fi = LorentzIndex[Unique[]];
    temp2 =
    scev[spx[[le-3]],spx[[le-2]]] spt@@Append[Drop[
                                              Drop[spx,{le-3,le-2}],-1
                                                  ], DiracGamma[5]]-
    scev[spx[[le-3]],spx[[le-1]]] spt@@Append[Drop[
                                              Drop[spx,{le-3,le-3}],-2
                                                  ], DiracGamma[5]]+
    scev[spx[[le-2]],spx[[le-1]]] spt@@Append[Drop[spx,-3],
                                              DiracGamma[5]
                                             ] +
    ( I Eps[spx[[le-3]],spx[[le-2]],spx[[le-1]],fi] *
       spt @@ Append[Drop[spx,-4],fi]
    );
      ];

    temp2 = temp2/.spt->spursavga /. spursavga -> spug;
 
    Throw[conall[temp2] ]
     (*if Global`$BreitMaison===True*)
      ]
    ];
       ]//pairexpand//Expand ;
            resp];

(* #################################################################### *)
(*                             Main50                                   *)
(* #################################################################### *)

(* ************************************************************** *)
(* Properties and special cases of traces (up to a factor 4) *)
   tris[x___] := tris[x] = trI[x];                  (*trisdef*)
   trI[a_+b_] := tris[a] + tris[b];                  (*trIdef*)
   trI[] = 1;
   trI[ DiracGamma[5] ] = 0;
   trI[ DiracGamma[6] ] = 1/2;
   trI[ DiracGamma[7] ] = 1/2;
 
 
   trI[ a:DiracGamma[_[__]].. ,DiracGamma[n_] ] := 0 /;
      (OddQ[Length[{a}]]&&(n==5 || n==6 || n==7));
 
   If[ Global`$BreitMaison === False,
       trI[ a:DiracGamma[_[__],___].. ,DiracGamma[n_] ] := 0 /;
          (OddQ[Length[{a}]]&&(n==5 || n==6 || n==7))
    ];
 
   trI[ d:DiracGamma[__].. ] := 0/;(OddQ[Length[{d}]] && fr567[ d ]);
 
   trI[ d:DiracGamma[_[__],___].. ,DiracGamma[5] ] := 0/;Length[{d}]<4;
 
   trI[x_] :=  x /; FreeQ[ {x},DiracGamma ];
 
   trI[ DiracGamma[a_[b__],___],DiracGamma[c_[d__],___],
        DiracGamma[6] ] := 1/2 scev[ a[b],c[d] ];
 
   trI[ DiracGamma[a_[b__],___],DiracGamma[c_[d__],___],
        DiracGamma[7] ] := 1/2 scev[ a[b],c[d] ];
 
   trI[ x__ ] := spursav[ x ]/;( Length[{x}]<11 && fr567[x]);
 
(* #################################################################### *)
(*                             Main51                                   *)
(* #################################################################### *)

                       (*trICdef*)
   $Larin = False;
   $Kreimer = False;
(* cyclic property *)
   trIC[y___]:=If[$Kreimer =!= True,
                  tris @@ cyclic[y],
                  tris[y]
                 ];
   cyclic[x__]:=RotateLeft[{x},Position[{x},First[Sort[{x}]]][[1,1]]];
   cyclic[]:={};


(*
 Options[ TrickMandelstam ] = {CollectMandelstam->False,
                               CutTrickMandelstam->2342};
*)

(* #################################################################### *)
(*                             Main52                                   *)
(* #################################################################### *)

(*ntermsdef*)
  nterms[x_Plus]:=Length[x];
   nterms[x_]:=Block[{ntermslex = Expand[x]},
                     If[ Head[ntermslex]===Plus,
                         ntermslex = Length[ntermslex],
                         If[x===0, ntermslex = 0, ntermslex = 1]
                       ];
           ntermslex];
(*nsortQdef*)
   nsortQ[x_,y_]:=True/;nterms[x]<=nterms[y];
   nsortQ[x_,y_]:=False/;nterms[x]>nterms[y];

(*combinedef*)
 combine[x_]:=Combine[x];



  (*TrickMandelstamdef*)
  TrickMandelstam[x_, man_List]:=x /; Length[man]=!=4;

(* prefer m^2-M^2 instead of (m-M)*(m+M) *)
   factor3[x_]:=factor3[x]=Factor2[x, FactorFull -> False];

(* #################################################################### *)
(*                             Main53                                   *)
(* #################################################################### *)
   TrickMandelstam[ y_, __ ] := factor3[y] /; FreeQ[y,Plus];
   TrickMandelstam[x_,es_,te_,uu_, mas_]:=TrickMandelstam[x, {es,te,uu,mas}];
   TrickMandelstam[x_List,y__]:=Map[TrickMandelstam[#,y]&, x];
   TrickMandelstam[a_ , {es_, te_, uu_, mm_}]:=Block[{tres},
       tres = trickmandelstam[a//factor3, {es,te,uu,mm}];
       If[LeafCount[tres]<2000, tres = Cancel[tres]];
          tres];

   trickmandelstam[ yy_Times, args_List ] := Map[ TrickMandelstam[ #,args ]&,yy ];
   trickmandelstam[ yy_Power, args_List ] := TrickMandelstam[yy[[1]],args ]^yy[[2]];
   trickmandelstam[ y_, args_List ]:= Block[{nulLl},
    trickmandelstam[nulLl+y,args]/.nulLl->0] /; (Head[y]=!=Times) &&
     (Head[y]=!=Power) && (Head[y]=!=Plus);

   drickstu[exp_,{},___] := exp;   
   drickstu[exp_,{s_,t_,u_,ma_},___] := exp /; !FreeQ[{s,t,u},Plus];

   short1[x_Plus,es_,te_,uu_,ma_]:=(Sort[{x, Expand[ x/.te->(ma-es-uu) ],
                                         Expand[ x/.uu->(ma-te-es) ]
                                     },nsortQ]//First );
   short1[a_ b_,c__]:=short1[a,c] short1[b,c];
   short1[a_^n_,c__]:=short1[a,c]^n;

   short1[x_,__]:=x/;(Head[x]=!=Plus) && (Head[x]=!=Times) && 
                     (Head[x]=!=Power);
 trickmandelstam[x_Plus,man_List]:=Block[{tricktemp,merk,nx=x,plusch, plusch0},
   merk[y_]:=memset[ merk[y],drickstu[y,man] ];
   plusch0[z__]:= Plus[z] /; !FreeQ[{z},plusch0];
(* This is for arguments of D0, etc. ... *)
   plusch[z__]:=drickstu[Plus[z],man]/;(Length[{z}]===(Length[Plus@@man]-1))&&
                                       FreeQ[{z},Plus];
   plusch[z__]:=(factor3 /@ Collect2[ Plus[z], Take[man, 3] ] ) /; 
                 Length[{z}]=!=(Length[Plus@@man]-1);
   tricktemp = merk[ nx ];
   (tricktemp/.Plus->plusch0/.plusch0->plusch /.
   (*drickstu->drickback/.*)plusch->Plus)]/;(Length[man]===4 || man==={}) &&
                                                  Head[x]=!=Times;

   drickback[x_,__]:=x;

   drickstu[ x_Plus,{s_,t_,u_,m_}  ]:=
     Block[{result,tristemp,eM,otherv,nuLL,trickman},
(* Check if an overall factorization is possible *)
      tristemp = factor3[ x/.s->(m-t-u) ];
  If[Head[tristemp]=!=Plus, 
     result = TrickMandelstam[tristemp,{s,t,u,m}],
         otherv = Complement[ Variables[tristemp], Variables[s+t+u+m] ];
(* The simplifications cannot occur outside certain coefficients *)
         If[ otherv =!= {},
             result = factor3/@ (Collect2[ eM tristemp, Append[otherv,eM] ]
                                );
             result = Map[short1[#,s,t,u,m]&, result+nuLL]/.nuLL->0/.eM->1;
             result = Map[factor3, result]
            ,
             result = short1[tristemp, s,t,u,m]
           ]
    ];
    result];
(* ------------------------------------------------------------ *)

(* #################################################################### *)
(*                             Main54                                   *)
(* #################################################################### *)
(*PolarizationSumdef*)

Options[PolarizationSum] = {Dimension -> 4};

PolarizationSum[mu_,nu_, ops___Rule]:= 
    -MetricTensor[mu, nu, Dimension -> 
          (Dimension /. {ops} /. Options[PolarizationSum])];  
PolarizationSum[mu_,nu_,k_, ops___Rule]:= 
   -MetricTensor[mu,nu, Dimension ->
          (Dimension /. {ops} /. Options[PolarizationSum])] +
     FourVector[k,mu] FourVector[k,nu]/
       Factor2[ExpandScalarProduct[ScalarProduct[k,k]]];
PolarizationSum[mu_,nu_,k_,0, ops___Rule] := -MetricTensor[mu, nu,
 Dimension -> (Dimension /. {ops} /. Options[PolarizationSum])];
PolarizationSum[mu_,nu_,k_,n_, ops___Rule]:=Collect2[
  -MetricTensor[mu,nu, Dimension ->
          (Dimension /. {ops} /. Options[PolarizationSum])
               ] - FourVector[k,mu] FourVector[k,nu]/
               Factor2[scev[Momentum[k],Momentum[n]]^2] *
               Factor2[scev[Momentum[n],Momentum[n]]] +
            ( FourVector[n,mu] FourVector[k,nu] +
              FourVector[n,nu] FourVector[k,mu] )/
               Factor2[scev[Momentum[k],Momentum[n]]], Pair];
(* #################################################################### *)
(*                             Main55                                   *)
(* #################################################################### *)
     (*sma*)
small2/: small2[x_]^n_ := small2[x^2] /; n > 0;
small2/: small2[_] a_ :=0;
small3/: small3[_] + a_ :=a;
small4[x_^m_]:=Negligible[x]^m;
   sma[x_]:=x/;FreeQ[x,Negligible];
   sma[x_]:=x/.Negligible->small2/.small2->small3/.
                         small3->small4/.small4->Negligible;


setit[a_,b_,___]:=set[a,sma[(b//Expand)]]/.set->Set;
(* SetMandelstamdef *)
SetMandelstam[s_,t_,u_, { {p1_, m12_}, {p2_, m22_} } ->
                        { {p3_, m32_}, {p4_, m42_} }] :=
SetMandelstam[s,t,u, p1, p2, -p3, -p4, Sqrt[m12], Sqrt[m22],
                                       Sqrt[m32], Sqrt[m42]];
SetMandelstam[s_,t_,u_,p1_,p2_,p3_,p4_,m1_,m2_,m3_,m4_]:=Block[
      {settemp,oldmem,setvars,sol,pp1, pp2, pp3, pp4},
      (* By definition the momenta p1, p2, p3, p4 are 4-dimensional, thus: *)
      (* note that p1, p2, p3, p4 may have have a minus - sign *)
      {pp1, pp2, pp3, pp4} = #/NumericalFactor[#] & /@ {p1, p2, p3, p4};
      { setdel[ Momentum[pp1, _Symbol], Momentum[pp1] ],
        setdel[ Momentum[pp2, _Symbol], Momentum[pp2] ],
        setdel[ Momentum[pp3, _Symbol], Momentum[pp3] ],
        setdel[ Momentum[pp4, _Symbol], Momentum[pp4] ]
      } /. setdel -> SetDelayed;
      oldmem = $MemoryAvailable;
      $MemoryAvailable = 0;
      settemp = {ScalarProduct[p1,p1] == m1^2,
                 ScalarProduct[p2,p2] == m2^2,
                 ScalarProduct[p3,p3] == m3^2,
                 ScalarProduct[p4,p4] == m4^2,
                 ScalarProduct[p1,p2] == sma[1/2 s - 1/2 m1^2 - 1/2 m2^2],
                 ScalarProduct[p1,p3] == sma[1/2 t - 1/2 m1^2 - 1/2 m3^2],
                 ScalarProduct[p1,p4] == sma[1/2 u - 1/2 m1^2 - 1/2 m4^2],
                 ScalarProduct[p2,p3] == sma[1/2 u - 1/2 m2^2 - 1/2 m3^2],
                 ScalarProduct[p2,p4] == sma[1/2 t - 1/2 m2^2 - 1/2 m4^2],
                 ScalarProduct[p3,p4] == sma[1/2 s - 1/2 m3^2 - 1/2 m4^2]
                }//pairexpand//Expand;
     setvars = Select[Variables[settemp/.Equal->List], SameQ[Head[#],Pair]&];     sol=Solve[ settemp,setvars ]/.Rule->setit;
     $MemoryAvailable = oldmem;
     sol
   ];

(* #################################################################### *)
(*                             Main56                                   *)
(* #################################################################### *)

SetMandelstam[x_, pl_List, ml_List]:=Block[
          {settemp,oldmem,setvars,sol,n=Length[ml], psu,pkl,sq2,eqq,ppl},
      oldmem = $MemoryAvailable;
      $MemoryAvailable = 0;
     ppl = #/NumericalFactor[#] & /@ pl;
     Table[ setdel[ Momentum[ppl[[ij]], _Symbol], Momentum[ppl[[ij]]] ],
            {ij, Length[ppl]} ] /. setdel -> SetDelayed;

      settemp = Join[ Table[ScalarProduct[pl[[i]], pl[[i]]] == ml[[i]]^2,
                            {i,1,n}],
                      Table[ScalarProduct[pl[[j]], pl[[j+1]]] == 
                    sma[1/2 x[j,j+1] - 1/2 ml[[j]]^2 - 1/2 ml[[j+1]]^2],
                            {j,1,n-1}],
                     {ScalarProduct[ pl[[1]],pl[[n]] ] ==
                    sma[1/2 x[1,n] - 1/2 ml[[1]]^2 - 1/2 ml[[n]]^2]}
                    ]//ExpandScalarProduct//Expand;


     setvars = Select[Variables[settemp/.Equal->List], SameQ[Head[#],Pair]&];   
  sol=Solve[ settemp,setvars ]/.Rule->setit;

  sq2[y_]:=ScalarProduct[y, y]//ExpandScalarProduct//Expand;
  pkl = {};
  For[ k=1, k<=n, k++,
       For[ l=k+1, l<= n, l++,
            npk = ScalarProduct[ pl[[k]], pl[[l]] ]//pairexpand;
            If[ (Head[npk] === Pair) || (Head[-npk]=== Pair),
                AppendTo[pkl,{pl[[k]], pl[[l]]}] 
              ]
          ]
     ];          
            
  psu = Plus@@pl;
  enm[a_]:=Expand[ - Apply[ Plus, Drop[pl,{a,a}] ]  ];
 (* p46 *)

Do[
  eqq = {sq2[psu] == 0};
  eqq = Join[ eqq, Table[ sq2[pl[[ii]] +  pl[[n]]] -
                          sq2[enm[ii] + pl[[n]]] ==0 , {ii, 2,n-3}]
            ];
  For[ j1 = 1, j1<n-2, j1++,
       For[ j2 = j1 + 2, j2<n, j2++,
            If[ EvenQ[j2-j1],
                AppendTo[ eqq, sq2[pl[[j1]] + pl[[j2]]] -
                               sq2[pl[[j1]] + enm[j2] ] ==0
                        ],
                AppendTo[ eqq, sq2[pl[[j1]] + pl[[j2]]] -
                               sq2[enm[j1]  + pl[[j2]]] ==0
                        ]
     ]    ]   ];
  var =  pairexpand[ScalarProduct@@#&/@pkl];
  var = Select[ Variables[var], Head[#]===Pair&];
  If[Length[var] > 0,
  nsol = Solve[ eqq, var ];
  nsol = nsol /. Rule -> setit
    ]
, {2}];

  $MemoryAvailable = oldmem;
  MapAll[ Expand, Append[sol, nsol]//Flatten ]
   ];

(* Eq. ... R.M. Dissertation *)
(* Ex.:  SpecificPolarization[ ..., Polarization[r] -> {0,a,b} ... ] *)
(* 0 :  Parallel *)
(* 1 :  Orthogonal *)
(* 2 :  Longitudinal *)
(* "+":  Righthanded*)
(* "-" :  Lefthanded*)

powerexpand[x_]:=x;
(* eq. 5.23 of A.Denner *)
standmat[]=1;
standmat[Spinor[c_. Momentum[pe1_],0,i_] ,
                       DiracGamma[ch_] ,
                       Spinor[Momentum[pe2_],0,i_]
                      ] := Sqrt[
    Expand[ ExpandScalarProduct[ 2 ScalarProduct[pe1, pe2]
          ]                    ]] /; (ch === 6) || (ch === 7) && (c^2)===1;

standmat[Spinor[c_. Momentum[pe1_],0,i_] ,
                       DiracGamma[Momentum[ka_]] ,
                       DiracGamma[ch_] ,
                       Spinor[Momentum[pe2_],0,i_]
                      ] := Sqrt[
    Expand[ ExpandScalarProduct[ 4 ScalarProduct[pe1, ka] *
            ScalarProduct[pe2, ka] - 2 ScalarProduct[pe1, pe2] *
            ScalarProduct[ka, ka]
         ]                    ]] /; ((ch === 6) || (ch === 7)) && ((c^2)===1) &&
                                    FreeQ[ka, Polarization];

(* SpecificPolarizationdef *)
SpecificPolarization[exp_]:=exp/.Dot->standmat/.standmat->Dot;
SpecificPolarization[exp_, polrule__Rule] := 
Block[{new, polrule2, i,mrule,standcalc,lk,factti},
new = exp /. Dot->standmat /. standmat-> Dot;
new = Contract[ dotLin[new], Expanding -> False ];
If[Head[new]===Polarization, new = Momentum[new]];
polrule2 = {polrule};
mrule[a_,b_]:=Momentum[a]->b;
If[FreeQ[polrule2, Momentum], 
   polrule2 =  polrule2/.Rule->mrule
  ];
pmc[xx_]=xx /. {n_."+" :> n "-",  "-" m_. :> "+" m};

        lk/: HoldForm[lk[x_]]:=lk[x];
(* Loop over the specific polarizations *)
For[i = 1, i <= Length[polrule2], i++,
If[ MemberQ[{0,1,2,"+","-"}, polrule2[[i,2,1]]],
oldnew = new;
    new = ( new /. polrule2[[i,1]] -> pols[ polrule2[[i,1]], 
                                            polrule2[[i,2]]
                                          ] 
          );
If[oldnew===new, 
    new = ( new /. Conjugate[polrule2[[i,1]]] -> pols[ polrule2[[i,1]], 
                                            polrule2[[i,2]]//pmc
                                          ] 
          )
  ];
 ](* *);
print1["substituted "];
    new = Contract[ dotLin[ new ], Expanding -> False ];
    new = ExpandScalarProduct[new];
    If[!FreeQ2[new, {LorentzIndex, Eps}],
        new = Contract[ new, EpsContract->True ] // sma;
      ];
    If[!FreeQ[new, DiracGamma],
       new = Collect2[ new, {Spinor,DiracGamma,Momentum,
                             Polarization}, IsolateHead->FF,
                             IsolateSplit -> Infinity,Factoring -> True];
        new = DiracSimplify[ new ] //Contract// ExpandScalarProduct;
        If[!FreeQ2[new, {LorentzIndex, Eps}],
            new  = Contract[ new, EpsContract->True ] // sma;
          ];
        new = new /.  Dot->standmat /. standmat-> Dot;
print1["simplified"];
      ];
        new = FixedPoint[ #/.FF->lk/.lk->FF &, new];
(*exper *) ];
print1["last ", LeafCount[new]];
(*
new  = Contract[ new, EpsContract->True ]//ExpandScalarProduct;
*)
If[!FreeQ[new, Eps],
print1["epsstuff"];
   If[!FreeQ[new, DiracGamma],
       new = DiracSimplify[Collect2[ new, DiracGamma, Factoring -> True]];
     ];
   new  = Contract[ new, EpsContract->True ] //ExpandScalarProduct; 
   If[!FreeQ2[new, {Eps, LorentzIndex}],
       print1["chisholm"];
       new = Collect[new, Eps]//EpsChisholm;
       print1["chisholm done"];
     ];
   new = new /.  Dot->standmat /. standmat-> Dot;
  ];

If[ LeafCount[new] < 1000, 
    print2["factoring"];
    new = Factor2[new],
    If[!FreeQ[new, Dot],
       print1["collecting now"];
(*
       factti = FactorTime /. Options[Factor2];
       SetOptions[Factor2, FactorTime -> 42];
*)
       new = Collect2[ new, Dot, Factoring -> False];
(*
       SetOptions[Factor2, FactorTime -> factti];
*)
     ];
  ];
print1["exiting"];
new];

pols[ Momentum[Polarization[r_, ii___]], {spe_,a_,b_}]:= 
pols[ Momentum[Polarization[r, ii]], {spe,a,b}]= 
Block[{sc,resu, dummymu, star, npc, spec},
spec = spe;
sc[x_, y_]:= sc[x,y]=ExpandScalarProduct[ ScalarProduct[x,y] ]//Factor2;
(* Remember that internally a complex conjugate polarization vector is 
   denoted by Polarization[x, -I]
*)
star[]=1;
star[ I]=1;
star[-I]=-1;

(* This is for the possible change of definition of "+" and "-" *)
If[!NumberQ[spec], npc = NumericalFactor[spec]; spec=spec/npc, npc =1];

resu = Momentum[Polarization[r,ii]];
resu = Which[
             (spec === "+") || (spec === "-"), 
             npc 1/Sqrt[2] pols[Momentum[Polarization[r, ii]], {0, a,b}
                               ]  + 
(* is this right ???? *)
                 npc star[ii] ToExpression[StringJoin[spec,"1"]]* 
             I/Sqrt[2] pols[Momentum[Polarization[r, ii]], {1,a,b}
                               ],
             spec === 0, 
       powerexpand[ 1/Sqrt[  Expand[   
         sc[a,b] (2 sc[a,r] sc[b,r] - sc[r,r] sc[a,b]) + 
         sc[a,a] sc[b,b] sc[r,r] - sc[b,b] sc[a,r]^2 - 
         sc[a,a] sc[b,r]^2      ] 
                          ] 
                  ] * 
       powerexpand[ 1/Sqrt[  Expand[ 
          ( (sc[b,r] + sc[a,r])^2 - sc[r,r] sc[a+b, a+b] )
              ] ]       ] * 
      ( Momentum[b] ( sc[a,a+b] sc[r,r] - 
                       sc[a,r] (sc[b,r] + sc[a,r])
                     ) - 
        Momentum[a] ( sc[b,a+b] sc[r,r] -
                       sc[b,r] (sc[b,r] + sc[a,r])
                     ) +
        Momentum[r] ( sc[a,b] ( sc[a,r] - sc[b,r] ) +
                       sc[b,b] sc[a,r] - sc[a,a] sc[b,r])
      ),

     spec === 1,
       dummymu = Unique[$MU];
       powerexpand[ 1/Sqrt[ Expand[ 
         sc[a,b] (2 sc[a,r] sc[b,r] - sc[r,r] sc[a,b]) +
         sc[a,a] sc[b,b] sc[r,r] - sc[b,b] sc[a,r]^2 -
         sc[a,a] sc[b,r]^2  
              ] ] ] * LorentzIndex[dummymu] *
       EpsEvaluate[
       Eps[ LorentzIndex[dummymu], Momentum[b], Momentum[a],
            Momentum[r]
          ]       ],

     spec === 2,
       powerexpand[ 1/Sqrt[ Expand[
         sc[r,r] ( (sc[b,r] + sc[a,r])^2 - 
                      sc[r,r] sc[b+a,b+a] ) 
              ] ] ] *
       ( Momentum[r] ( sc[b,r] + sc[a,r]) -
         Momentum[a+b] sc[r,r]
       )
      ];
resu = Collect2[ resu, Momentum, Factoring -> True];
resu];
 
End[];
WriteString["stdout", "."];
$PrePrint = FeynCalcForm;
(* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ *)
Begin["Feyn`Calc`Main`"];
(* ******** *)

ChangeDimension[x_, diim_] := Block[{xx, dirGG, DIRGAMM,pAiR, ld , md,spi4},
If[ diim === 4, xx = x /. LorentzIndex[a_,___] -> LorentzIndex[a] /.
                     Momentum[b_,___] -> Momentum[b],
ld[a_,___] := LorentzIndex[a,diim];
md[a_,___] := Momentum[a,diim];
spi4[a_, b__] := Spinor[ ToFourDimensions[a], b];
xx = x /. Pair -> pAiR /. DiracGamma->DIRGAMM;
xx = xx /. Momentum -> md /. LorentzIndex -> ld;
dirGG[aa_,___] := DIRGAMM[aa, diim];
xx = xx /. DIRGAMM -> dirGG/.  DIRGAMM -> DiracGamma /. pAiR -> Pair;
xx = xx /. Spinor -> spi4;
  ];
                                 xx];

 duM=Unique[System`C]; 
(*
 duM=Unique[dUU]; 
*)
gamcount[x_] := Length[Position[x, DiracGamma]];
specsir[a_] := Block[{tem, nosp},
                      nosp = Select[a, FreeQ[#, Spinor]&]/.$MU-> duM;
                      tem = (a/nosp)/.$MU-> duM;
                      If[!FreeQ[tem, duM[2]],
                         tem = Contract[DiracOrder[tem,
                                        {duM[1], Momentum[_], duM[2]}
                                                  ]//DiracSimplify]
                        ];
                      If[Length[Position[tem, Spinor]] > 2, 
                         print2["entering specsir2 with ",tem//FeynCalcForm];
                         tem = specsir2[tem];
                         print2["exiting specsir2 with", tem//FeynCalcForm];
                        ];
                     Expand[tem nosp]
                    ];
ChisholmSave[x_]    := memset[ChisholmSave[x], Chisholm[x]];
EpsChisholmSave[x_] := memset[EpsChisholmSave[x], EpsChisholm[x]];
specsir2[x_] := specsir2[x] = 
  Block[{ste, ste1, ste2,ste1r,ste2r,duMy},
        ste = x /. $MU -> duMy;
   (* there are two possibilities ... *)
   ste1r ={Literal[Spinor[pe1__]] . g1__ . Literal[Spinor[pe2__]] *
           Literal[Spinor[pe3__]] . g2__ . Literal[Spinor[pe4__]] *
           Literal[Spinor[pe5__]] . g3__ . Literal[Spinor[pe6__]] :>
  (Expand[DiracSimplify[
      DiracOrder[DiracSimplify[Contract[Spinor[pe1] . g1 . Spinor[pe2] *
                    EpsChisholmSave[ Spinor[pe3] . g2 . Spinor[pe4] *
                                 ChisholmSave[Spinor[pe5] . g3 . Spinor[pe6]]
                        ] ] ], {$MU[1], Momentum[_], $MU[2]} ] 
           ]] /. $MU -> duMy
    )/; 
      (Length[{g1}] < 3) && (Length[{g2}] < 3) && (Length[{g3}] > 2)
          };

ste2r = {Literal[Spinor[pe1__]] . g1__ . Literal[Spinor[pe2__]] *
         Literal[Spinor[pe3__]] . g2__ . Literal[Spinor[pe4__]] *
         Literal[Spinor[pe5__]] . g3__ . Literal[Spinor[pe6__]] :>
  (Expand[DiracSimplify[
      DiracOrder[DiracSimplify[
           Contract[Spinor[pe3] . g2 . Spinor[pe4] *
                    EpsChisholm[ Spinor[pe1] . g1 . Spinor[pe2] *
                                 Chisholm[Spinor[pe5] . g3 . Spinor[pe6]]
                        ] ] ],{$MU[1], Momentum[_], $MU[2]} ] 
            ]] /. $MU -> duMy
   )/; 
      (Length[{g1}] < 3) && (Length[{g2}] < 3) && (Length[{g3}] > 2)
        };
ste1 = ste/.ste1r;
If[(!FreeQ[ste1, DiracMatrix[duMy[_]].DiracGamma[Momentum[_]].
                DiracMatrix[duMy[_]]
         ]) && (ste1=!= ste), 
   ste1 = ste1 /. ste2r /. ste1r
  ];
If[(ste1 =!= 0) && (!FreeQ[ste1, DiracMatrix[duMy[_]].DiracGamma[Momentum[_]].
                        DiracMatrix[duMy[_]]
                       ])  , 
  ste2 = ste /. ste2r;

If[(!FreeQ[ste2, DiracMatrix[duMy[_]].DiracGamma[Momentum[_]].
                DiracMatrix[duMy[_]]
         ]) && (ste2=!= ste), 
    ste2 = ste2 /. ste1r /. ste2r
  ],
 ste2 = ste1
  ];
If[ste1 =!= 0,
  ste1 = DiracSimplify[DiracOrder[DiracSimplify[ste1]]]//Expand;
  ste2 = DiracSimplify[DiracOrder[DiracSimplify[ste2]]]//Expand;
  If[!FreeQ[ste1, DiracMatrix[duMy[_]].DiracGamma[Momentum[_]].
                DiracMatrix[duMy[_]]],
     ste1 = First[Sort[Select[{ste1, ste2}, FreeQ[#,RuleDelayed]&], 
             (LeafCount[{##}[[1]]] < LeafCount[{##}[[2]]] )& ] ]
    ];
  ];
ste1];
(* end of specsir2 *)

(* **************************************************************** *)
(* SquareAmplitudedef *)
(* **************************************************************** *)
(* careful: this function is still under development!!! *)

(* ----------------------------------------------------------------- *)
susa[xxx_] := susa[xxx] = SUNSimplify[xxx, Factoring->True];
(* /. SUNTrace -> ST; *)
susm[nuLL]:=0;
susm[xx_Plus]                   := Map[susm, xx];
susm[xx_ /; FreeQ[xx,SUNIndex]] := xx;
susm[xx_ /;( (Head[xx] =!= Plus) && (!FreeQ[xx, SUNIndex]))
    ]                           := Expand[
          Select[xx, FreeQ[#1, SUNIndex] & ]*
    susa[Select[xx, !FreeQ[#1, SUNIndex] & ]] (*, SUNIndex*)];

Options[SquareAmplitude] = {
                  Dimension -> 4, 
                  EnergyMomentumConservation -> {},
                  EpsAway -> False,
                  ExtraFactor -> 1,
                  Collecting -> False,
                  Factoring -> False,
                  FinalSubstitutions -> {},
                  InitialSubstitutions -> {},
                  IntermediateSubstitutions -> {},
                  IsolateHead -> K,
                  IsolateSplit -> 4711 I,
                  Mandelstam -> {Global`S, Global`T},
                  SelectGraphs -> All,
                  SpinPolarizationSum -> 1,
                  SpinSumExternalMomentum -> Automatic,
                  WriteOut -> False
                            };

SquareAmplitude[ FeynAmpList[he__, procname_Rule, process_Rule][amps__],
    opts___Rule] := Block[
                           {sli = {}, amp = {amps}, proctype, dim, nam,
                            colorpart = 1,pluiso,gluON,den1,numfaN,enmomsubst,
                            enmomback,inisubst,pli,
                            prn,mand,factoring,exmom,exm, es,te, sel, scalP,
         p1,p2,k1,k2,k3,k4,m12,m22,m32,m42,m52,m62,
                            finsubst,proc,extrafact},
    prnam = procname[[2]];
prn       = StringJoin[ToString[prnam], ".sq"];
dim       = Dimension   /. {opts} /. Options[SquareAmplitude];
mand      = Mandelstam  /. {opts} /. Options[SquareAmplitude];
inisubst  = InitialSubstitutions /.  {opts} /. Options[SquareAmplitude];
collecting= Collecting /. {opts} /. Options[SquareAmplitude];
factoring = Factoring   /. {opts} /. Options[SquareAmplitude];
finsubst  = FinalSubstitutions /. {opts} /. Options[SquareAmplitude];
isolhead  = IsolateHead /. {opts} /. Options[SquareAmplitude];
isolsplit = IsolateSplit/. {opts} /. Options[SquareAmplitude];
exmom     = SpinSumExternalMomentum /. {opts} /. Options[SquareAmplitude];
extrafact = ExtraFactor /. {opts} /. Options[SquareAmplitude];

{es, te} = {mand[[1]], mand[[2]]};
If[Length[mand === 3], uu = mand[[3]]];

sel = SelectGraphs /. {opts} /. Options[SquareAmplitude];
If[sel === All, sli = {amps}, For[jj = 1, jj <= Length[amp], jj++, 
                                  If[MemberQ[sel, jj], AppendTo[sli, amp[[jj]]]] ] 
  ];

plsI[xx__] := Isolate[Plus[xx], {DiracGamma, Spinor, LorentzIndex, SUNIndex},
                     IsolateHead -> isolhead, IsolateSplit->555I];

(* *********************************** *)
{es, te} = {mand[[1]], mand[[2]]};
If[Length[mand===3], uu = mand[[3]]];
proc  = Last[process];
proctype = Map[Length, proc];
Which[ 
       proctype === (2 -> 1),
       proc = proc[[1]] -> Join[proc[[2]], {{0,0,0}, {0,0,0}, {0,0,0}}],
       proctype === (1 -> 2),
       proc = Join[proc[[1]], {{0,0,0}}] -> Join[proc[[2]], {{0,0,0}, {0,0,0}, {0,0,0}}],
       proctype === (2 -> 2),
       proc = proc[[1]] -> Join[proc[[2]], {{0,0,0}, {0,0,0}}],
       proctype === (2 ->3),
       proc = proc[[1]] -> Join[proc[[2]], {{0,0,0}}]
     ];

    p1m1 = {proc[[1, 1, 2]], proc[[1, 1, 3]]^2};
    p2m2 = {proc[[1, 2, 2]], proc[[1, 2, 3]]^2};
    k1m3 = {proc[[2, 1, 2]], proc[[2, 1, 3]]^2};
    k2m4 = {proc[[2, 2, 2]], proc[[2, 2, 3]]^2};
    k3m5 = {proc[[2, 3, 2]], proc[[2, 3, 3]]^2};
    k4m6 = {proc[[2, 4, 2]], proc[[2, 4, 3]]^2};
 fields = { proc[[1,1,1]], proc[[1,2,1]], 
            proc[[2,1,1]], proc[[2,2,1]], proc[[2,3,1]], proc[[2,4,1]]
          };

(* this has to be thought over again ... *)
extrafact = extrafact ( NF^0 (*( (Length[Position[fields, Global`F[I]]] + 
                             Length[Position[fields, Global`F[-I]]])/2
                             )*)
                      );
extrafact = ExpandScalarProduct[extrafact]//Expand;

 {p1, p2, k1, k2, k3, k4} = #[[1]]& /@ {p1m1, p2m2, k1m3, k2m4, k3m5, k4m6};
 {m12, m22, m32, m42, m52, m62} = #[[2]]& /@ {p1m1, p2m2, k1m3, k2m4, k3m5, k4m6};
If[MemberQ[{1 -> 2, 2 -> 1, 2 -> 2}, proctype],
    enmomsubst = k1 -> p1+p2-k2-k3-k4;
    enmomback = {p1+p2-k2-k3-k4->k1,
                 Momentum[p1]+Momentum[p2]-Momentum[k2]-
                 Momentum[k3]-Momentum[k4] -> Momentum[k1],
                -Momentum[p1]-Momentum[p2]+Momentum[k2]+
                 Momentum[k3]+Momentum[k4] -> (-Momentum[k1])
                },
     enmomsubst = {};
     enmomback  = {};
  ];

scalP[a_, b_, c_] := If[dim =!= 4,
 Apply[ Set, {ScalarProduct[a,b, Dimension -> dim], c//Expand}];
 Apply[ Set, {ScalarProduct[a,b, Dimension -> 4], c//Expand} ];
 Apply[ Set, {pair2[Momentum[a, dim], Momentum[b, dim]], c//Expand} ];
 Apply[ Set, {pair2[Momentum[a ], Momentum[b ]], c//Expand} ]
 , 
  Apply[ Set, {ScalarProduct[a,b, Dimension -> 4], c//Expand} ];
  Apply[ Set, {pair2[Momentum[a], Momentum[b]], c//Expand}];
                       ];
sCP[a_] := ScalarProduct[a, a, Dimension -> dim];
sCP[a_,b_] := ScalarProduct[a, b, Dimension -> dim];

If[Head[Pair[Momentum[p1], Momentum[p1]]] === Pair,
If[proctype  === (1 -> 2),
     scalP[p1, p1, m12];  scalP[k1, k1, m32];  scalP[k2, k2, m42];
     scalP[k1, k2, 1/2 (sCP[p1] - sCP[k1] - sCP[k2])];
     scalP[p1, k1, 1/2 (sCP[p1] - sCP[k2] + sCP[k1])];
     scalP[p1, k2, 1/2 (sCP[p1] - sCP[k1] + sCP[k2])];
  ];

If[ proctype  === (2 -> 1),
   scalP[p1,p1, m12];
   scalP[p2,p2, m22];
   scalP[k1,k1, m32];
   scalP[p1,p2, 1/2 (sCP[k1] - sCP[p1] - sCP[p2])];
   scalP[p1,k1,-1/2 (sCP[p2] - sCP[p1] - sCP[k1])];
   scalP[p2,k1,-1/2 (sCP[p1] - sCP[p2] - sCP[k1])];
  ];

uUu = Expand[m12 + m22 + m32 + m42 - es - te];
If[ proctype === (2 -> 2),
    scalP[p1, p1, m12];
    scalP[p1, p2, - m12/2 - m22/2 + es/2];
    scalP[p1, k1,   m12/2 + m32/2 - te/2];
    scalP[p1, k2,   m12/2 + m42/2 - uUu/2];
    scalP[p2, p2, m22];
    scalP[p2, k1,  m22/2 + m32/2 - uUu/2];
    scalP[p2, k2,  m22/2 + m42/2 - te/2];
    scalP[k1, k1, m32];
    scalP[k1, k2, -m32/2 - m42/2 + es/2];
    scalP[k2, k2, m42];
  ];

If[ proctype === (2 ->3),
    scalP[p1, p1, m12];
    scalP[p2, p2, m22];
    scalP[k1, k1, m32];
    scalP[k2, k2, m42];
    scalP[k3, k3, m52];
    scalP[p1, p2, Global`P[p1, p2]/2 - m12/2 - m22/2];
    scalP[p1, k1,-Global`P[p1,-k1]/2 + m12/2 + m32/2];
    scalP[p1, k2,-Global`P[p1,-k2]/2 + m12/2 + m42/2];
    scalP[p1, k3,-Global`P[p1,-k3]/2 + m12/2 + m52/2];
    scalP[p2, k1,-Global`P[p2,-k1]/2 + m22/2 + m32/2];
    scalP[p2, k2,-Global`P[p2,-k2]/2 + m22/2 + m42/2];
    scalP[p2, k3,-Global`P[p2,-k3]/2 + m22/2 + m52/2];
    scalP[k1, k2, Global`P[k1, k2]/2 - m32/2 - m42/2];
    scalP[k1, k3, Global`P[k1, k3]/2 - m32/2 - m52/2];
    scalP[k2, k3, Global`P[k2, k3]/2 - m42/2 - m52/2];
  ];
];

(*
If[(k3 =!= 0) && (k4 =!= 0),
    SetMandelstam[Global`P,  {p1, p2, -k1, -k2, -k3, -k4},
                             {p1m1[[2]], p2m2[[2]], k1m3[[2]],
                              k2m4[[2]], k3m5[[2]], k4m6[[2]]}
                 ]
  ];
*)

(* *********************************** *)

(* sum the amplitudes *)
amp = Sum[ PropagatorDenominatorExplicit[sli[[i,2]], {p1,p2,-k1,-k2,-k3,-k4}], {i,1, Length[sli]}
         ] /. inisubst;
amp0 = amp;
If[dim =!= 4,
   amp = ChangeDimension[amp, dim];
   extrafact = ChangeDimension[extrafact, dim];
  ];

print1["extrafact = ",extrafact];

If[!FreeQ[amp, DiracGamma],
   amp = DiracSimplify[amp /. enmomsubst /. 
                       Pair -> pPpPpP, Expanding -> False
                      ] /. enmomback /. pPpPpP->Pair,
     amp = ExpandScalarProduct[Contract[amp/.enmomsubst]] /. enmomback
  ];

print1["checK"];
sund[xxxx_] := 1;
sund[glu_, xxxx_] := SUNDelta[glu, SUNIndex[xxxx]];
gluON[]=Sequence[];
mU = Unique[System`C]; sU = Unique[System`C]; gL = Unique[System`C];
(*
mU = Unique[mUUn]; sU = Unique[sUUu]; gL = Unique[glLLu];
*)

gluON[_SUNIndex]=gL;

If[$VeryVerbose > 0, Print["Length of amp = " , Length[amp]]];
(* a list of rules for substituting indices for the polarization momenta *)
momtolor = {Momentum[Polarization[p1, _, glui___SUNIndex], ___] :>
                       (LorentzIndex[mU[1, gluON[glui]],dim] sund[glui, sU[1]]),
            Momentum[Polarization[p2, _, glui___SUNIndex], ___] :>
                       (LorentzIndex[mU[2, gluON[glui]],dim] sund[glui, sU[2]]),
            Momentum[Polarization[k1, _, glui___SUNIndex], ___] :>
                       (LorentzIndex[mU[3, gluON[glui]],dim]  sund[glui, sU[3]]),
            Momentum[Polarization[k2, _, glui___SUNIndex], ___] :>
                       (LorentzIndex[mU[4, gluON[glui]],dim]  sund[glui, sU[4]]),
            Momentum[Polarization[k3, _, glui___SUNIndex], ___] :>
                       (LorentzIndex[mU[5, gluON[glui]],dim]  sund[glui, sU[5]]),
            Momentum[Polarization[k4, _, glui___SUNIndex], ___] :>
                       (LorentzIndex[mU[6, gluON[glui]],dim]  sund[glui, sU[6]])
           };

pairdel[ c___,  n_ b_LorentzIndex, d___] :=
pairdel[c,n b, d] = n pairdel[c,b, d];

dirdel[n_ b_LorentzIndex, d___] := dirdel[n b, d] = n DiracGamma[b, d];

If[!FreeQ[amp, SUNIndex],
   (* collecting w.r.t. to SUNIndex *)
   print1["collecting w.r.t. SUNIndex"];
   If[ $FeynContract === True, 
       pair2PAIR[a_, b_] := If[FreeQ[{a,b}, SUNIndex], Global`PAIR[a,b], Pair[a, b]];
       print1["substing PAIR"];
       amp = amp /. Pair -> pair2PAIR;
amp00 = amp;
       print1["substing PAIR done"];
       amp = Global`CCollect[amp, {Pair, SUNF, SUNT, SUNDelta}];
       print1["back from CCollect"];
       amp = amp /. Global`PAIR -> Pair,
       amp = Collect2[amp, SUNIndex,Factoring -> False];
     ];
   amp11=amp;
   print1["SUNSimplifying"];
   amp = SUNSimplify[amp, SUNFToTraces -> False];
  ];

oldamp = amp;
amp = amp /. momtolor;

If[!FreeQ[amp, SUNIndex],
print2["pairdel"];
amp = amp/. Pair -> pairdel /. pairdel->Pair/.
            DiracGamma->dirdel/.dirdel->DiracGamma;
amp = ChangeDimension[amp, dim];
  print2["contracting"];

amp = Contract[ amp/.enmomsubst, Expanding -> False ];
amp = ExpandScalarProduct[amp] /. enmomback;
(*
amp = ChangeDimension[amp, dim];
*)

   print1["collecting w.r.t. SUNIndex"];
If[ $FeynContract === True,
       print1["substing PAIR( 2 )"];
       amp = amp /. Pair -> pair2PAIR;
       print1["substing PAIR( 2 ) done"];
amp00 = amp;
       amp = Global`CCollect[amp, {Pair, SUNF, SUNT, SUNDelta}];
       print1["back from CCollect"];
       amp = amp /. Global`PAIR -> Pair,
       amp = Collect2[amp, SUNIndex,Factoring -> False];
  ];

amp22 = amp;

amppp=amp;

If[Head[amp] === Times,
   colorpart = Select[amp, !FreeQ[#, SUNIndex]&];
   amp = Select[amp, FreeQ[#, SUNIndex]&],
    colorpart = 1;
   iii=0;
(*
   amp = Map[#/.Plus-> plUUU&, amp + nuLL1 + nuLL2];
   amp = Sum[(Print[iii++];susa[amp[[ij]]]), {ij, Length[amp]}];
   amp = Expand2[amp /. nuLL1 -> 0 /. nuLL2->0,SUNIndex] /. plUUU -> Plus;
*)
  ];
If[!FreeQ[colorpart, SUNIndex], print2["colorpart = ",colorpart]];
  ];

amp44=amp;
print2["Length of amp44 = ",Length[amp44]];
print2["leafcount of amp44 = ",LeafCount[amp44]];
If[Head[amp44] === Plus,
   namp = 0; lamp44 = Length[amp44];
   For[ia = 1, ia <= lamp44, ia++,
       print2["contracting amp44[[",ia,"]]  out of ", lamp44];
       namp = namp + ExpandScalarProduct[Contract[amp44[[ia]]]];
      ];
   amp = namp,
   amp = Contract[amp]//ExpandScalarProduct;
  ];
       
amp45 = amp;
print2["Length of amp45 = ",Length[amp45]];

nuLL = Unique[System`C];
(*
nuLL = Unique[nuLlllllL];
*)
Collect3[x_ /; (Head[x] =!= Plus), _]:=x;
Collect3[su_Plus, a_Symbol] := Block[{suntr, ll, j=0, vp, new = 0,cc,
                                      temp = su},
print1["this cridder is  ", Length[temp], "  terms long"];
While[ (ll = Length[temp]) > 0,
       j++;
       print2["j = ",j, "    length of temp = ",ll];
       If[Head[temp] === Plus,
          part1 = First[temp], part1 = temp
         ];
       vp = Select[part1 DUM, !FreeQ[#, a]&];
       print2["vp = ", vp//FeynCalcForm];
       cc = Coefficient[temp, vp];
       new = new + vp cc;
       temp = temp - ((#vp& /@ (cc + nuLL))/.nuLL->0)
     ];
new];

If[ 
      $FeynContract === True,
       print1["substing PAIR( 3 )"];
       amp = amp /. Pair -> pair2PAIR;
       print1["substing PAIR( 3 ) done"];
amp00 = amp;
       amp = Global`CCollect[amp, {Pair, SUNF, SUNT, SUNDelta}];
       print1["back from CCollect"];
       amp = amp /. Global`PAIR -> Pair,
       amp = Collect2[amp, SUNIndex,Factoring -> False];
  ];

amp55 = amp;

(* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC *)
(* construct the tensor from the polarization sums *)
(* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC *)

pli = {p1, p2, k1, k2, k3, k4};
pli2 = ExpandScalarProduct[ScalarProduct[#, #]]& /@ pli;

prod = extrafact;
exm = {0,0,0,0,0,0};
If[exmom =!= Automatic, 
   If[ Head[exmom] === List, exm = exmom];
   If[ Head[exmom] === Symbol, exm = {exmom, exmom,exmom,exmom,exmom,exmom}];
  ];
mom4set[xp_Symbol] := Apply[SetDelayed, {Momentum[xp,___Symbol], Momentum[xp]}];
Map[mom4set, exm];
print2["exm = ",exm];
print2["fields = ", fields];
For[ij = 1, ij < 7, ij++,
    If[pli =!= 0,
   (* check for photon *)
    If[(fields[[ij]] === Global`V[1]) &&  !FreeQ[amp, mU[ij]],
       prod = prod (- MetricTensor[mU[1], ComplexIndex[mU[1]],
                           Dimension -> dim])
      ];
   (* check for gluon*)
    If[(fields[[ij]] === Global`G[1]) && !FreeQ[amp, mU[ij,gL]],
       If[exmom === Automatic,
          Which[
                proctype === (2->3),
                exm = {pli[[2]], pli[[1]], pli[[5]], pli[[3]], pli[[4]]},
                proctype === (2->2),
                exm = {pli[[2]], pli[[1]], pli[[4]], pli[[3]]},
                proctype === (2->1),
                exm = {pli[[3]], pli[[1]], pli[[2]]},
                proctype === (1->2),
                exm = {pli[[3]], pli[[1]], pli[[2]]}
               ];
          ];
       prod = prod PolarizationSum[mU[ij, gL], ComplexIndex[mU[ij, gL]],
                                   pli[[ij]], exm[[ij]], Dimension -> dim];
       ];
   (* check for massive Vectorbosons *)
    If[MatchQ[fields[[ij]], Global`V[iii_ /; iii=!=1]]  &&  !FreeQ[amp, mU[ij]],
       prod = prod PolarizationSum[mU[ij], ComplexIndex[mU[ij]],
                         pli[[ij]], Dimension -> dim]
      ];
   (* check for ghosts *)
    If[MatchQ[fields[[ij]], Global`U[_]],
Print["GHOSTCHECK"];
       prod = I prod 
      ];
  ];
   ];
print1["polarization sums =  ", prod//FeynCalcForm];
If[dim =!= 4,
   prod = prod/.{Pair[LorentzIndex[a_], LorentzIndex[b_]]:>
                 Pair[LorentzIndex[a, dim], LorentzIndex[b, dim]],
                 Pair[LorentzIndex[a_], Momentum[b_]]:>
                 Pair[LorentzIndex[a, dim], Momentum[b, dim]],
                 Pair[Momentum[a_], Momentum[b_]]:>
                 Pair[Momentum[a, dim], Momentum[b,dim]],
                 DiracGamma[LorentzIndex[a_]] :>
                 DiracGamma[LorentzIndex[a,dim],dim],
                 DiracGamma[Momentum[a_]] :>
                 DiracGamma[Momentum[a, dim], dim]
                };
  ];
amp = Contract[ amp, Expanding -> False ];
(*
pluli[x__] := Collect2[ Plus[x], LorentzIndex , Factoring -> False ];
*)
pluli[x__] := Collect2[ Plus[x], LorentzIndex , Factoring -> True ];
amp66 = amp;
print1[Length[amp66]];

frhDot[xx__] := FRH[Dot[xx]];
amp = ((# /. Plus -> plsI&) /@ amp)/.Dot -> frhDot;

print1["leafcount = ",LeafCount[amp]];
print1["combining"];
amp1 = Combine[amp];
print1["combining done "];
den1 = Denominator[amp1];

print2["denominator = ",den1//FRH];
print3["LeafCount = ",LeafCount[amp1]];

colorfactor = susa[colorpart ComplexConjugate[colorpart]];
If[colorfactor =!= 1, 
  print1["the global colorfactor of the squared amplitudes is ", colorfactor];
  ];

nam = Numerator[amp1];
numfaN = NumericalFactor[nam];
nam = nam/numfaN;
den1 = den1/numfaN;

(*
If[FreeQ2[nam, {SUNIndex, DiracGamma}], 
   nam = pluli[nam](*,
   nam = Map[#/.Plus->pluli&, nam]; *)
  ];
*)

(*
Global`$ZWISCHEN = True;
If[Global`$ZWISCHEN === True, nam >> "nam.s"; den1 >> "den1.s"; prod >>"prod.s" ]; 
*)
prod = Collect2[prod, LorentzIndex];

amp = SquareAmplitude[nam, ExtraFactor -> (prod colorfactor), 
                            IsolateHead -> isolhead,
                            Mandelstam -> {}, Factoring -> False,
                            IntermediateSubstitutions -> enmomsubst,
                            FinalSubstitutions -> enmomback,
                            EnergyMomentumConservation -> 
                            {-p1,-p2,k1,k2,k3,k4 }
                     ];
(* ++++++++++++++++++++++++++++++++++ *)
print1["collecing "];
amp = amp  /. enmomsubst;
amp77 = amp;
$CheckCollect = False;
If[$CheckCollect === True,
   amp = Collect2[ amp, isolhead ]/.Plus -> plsI;
  ];

amp = amp / den1  ComplexConjugate[1/den1];
If[collecting === True, 
   amp = Collect2[ amp, isolhead ]
  ];
If[factoring =!= True,
   pluisol[xx__] := Isolate[Plus[xx],
                     IsolateHead -> isolhead, IsolateSplit->888];
   amp = Isolate[amp /. Plus -> pluisol, IsolateHead -> isolhead, 
                 IsolateSplit -> 888],

   amp = amp /. finsubst;
   amp = Collect2[ amp, isolhead ]//FRH//Factor2;
   amp = amp /. finsubst;
   amp = Factor2[amp];
   If[(proctype === (2->2)) && (Length[mand] === 4),
      amp = TrickMandelstam[amp, mand]
     ];
  ];
amp];


(* XFX *)
SquareAmplitude[ exp_ /; FreeQ[exp, FeynAmp], ops___Rule ] := Block[
{ts,sps,amp2,cts,plu3,neamp,intermed,amp, 
 amps1, amps2,pamp,lis, saveminamp, file,finsubst2, isolhead2,
 nwres, mand2, extrafact2,factoring2,epsAway},

ts = exp;
enmo = EnergyMomentumConservation /. {ops} /. Options[SquareAmplitude];
intermed = IntermediateSubstitutions /. {ops} /.
              Options[SquareAmplitude];
isolhead2 = IsolateHead /. {ops} /. Options[SquareAmplitude];
saveminamp = WriteOut    /. {ops} /. Options[SquareAmplitude];
extrafact2 = ExtraFactor /. {ops} /. Options[SquareAmplitude];
factoring2 = Factoring   /. {ops} /. Options[SquareAmplitude];
finsubst2  = FinalSubstitutions  /. {ops} /. Options[SquareAmplitude];
epsAway    = EpsAway /. {ops} /. Options[SquareAmplitude];

If[StringQ[saveminamp], file = FileNames @@ {saveminamp}, file = {}];

plu3[y__] := Isolate[Plus[y], 
              {DiracGamma, LorentzIndex,Spinor,Eps,SUNIndex},
                  IsolateHead -> isolhead2,
                     IsolateSplit -> 444I];

etl[x__] := plu3[x] /; FreeQ[{x}, Eps];
nuLL1 = Unique[System`C]; nuLL2 = Unique[System`C];
(*
nuLL1 = Unique[cCU]; nuLL2 = Unique[cCU];
*)
etl[y__] := Block[{ee},
                   ee = PartitHead[Plus[y], Eps];
                   (Collect2[ee[[2]]+nuLL1 + nuLL2, Eps,Factoring->False
                            ] /. nuLL1 -> 0 /. nuLL2 -> 0(*/.Plus->plu2*)) + 
(ee[[1]]/.Plus->plu3)
(*
(Factor2[ee[[1]], FactorTime->242]/.Plus->plu3)
*)
                 ] /; !FreeQ[{y}, Eps];
(* ++++++++++++++++++++++++++++++++++++ *)

print2["enmo = ", enmo];

If[file =!= {},
   print2["loading previous result ", file];
   (Get @@ {saveminamp});
   ts = Global`$TSAMP;
lis = Length[ts];
nuLlLL = Unique[System`C];
(* 
 nuLlLL = Unique[cCUn]; 
 *)
(*
print2["colling2insquareamplitudeagain"];
    ts = ((#/.Plus->etl)& /@ (ts + nuLlLL)) /. nuLlLL->0
*)
 ,

 ts = PropagatorDenominatorExplicit[ts, enmo];

print2["diracsimplifyinginsquareamplitude"];
ts = DiracSimplify[ts];
ots = ts;
print2[LeafCount[ots]];

If[Head[ts] === Plus,
print2["colling2insquareamplitude"];
    ts = Collect2[ts, {Spinor,SUNIndex}, Factoring -> False];
ots2 = ts;
(*
    ts = Collect2[ts, {Spinor,SUNIndex}, Factoring -> True];
*)
  ];

(* sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss *)
spmin[xxxx_, enmomcon_List] := Block[{xx=xxxx,
enms = (Map[Momentum, enmomcon]//ExpandScalarProduct), smmin,sol, sba, check, res=0, simP},
If[Length[enmomcon] === 0, res = xxxx,
               sol[pe_,___] := Solve[(Plus@@enms)==0, 
                        PowerExpand[Sqrt[pe^2]]][[1,1]];
               sba[pe_,___] := {Reverse[sol[pe]], Map[-#&,Reverse[sol[pe]]]};
              check[pe_, __, a_List]:= !FreeQ[a, Last[sol[pe][[2]]]/
                                        NumericalFactor[Last[sol[pe][[2]]]]];
Literal[smmin[Spinor[pe1__], a___, Spinor[pe2__]]] := 
 smmin[Spinor[pe1] . a . Spinor[pe2]] = 
 If[FreeQ[{a}, Momentum],
    Spinor[pe1] . a . Spinor[pe2],
  If[ check[pe1, {a}], 
      DiracSimplify[(Spinor[pe1]/.sol[pe1]) . a . Spinor[pe2]] /. sba[pe1],
      If[check[pe2, {a}],
         DiracSimplify[Spinor[pe1] . a . (Spinor[pe2]/.sol[pe2])] /. sba[pe2], 
         Spinor[pe1] . a . Spinor[pe2]
        ],
      Spinor[pe1] . a . Spinor[pe2]
    ]
   ];
specsu[x_]:=specsu2[Expand2[x, Spinor]];
specsu2[x_Plus]:=Map[specsu2, x];
simP[xy_] := (*simP[xy]=*)
specsu[
FixedPoint[Expand[DiracSimplify[
             Contract[#(* /.Dot->smmin /. smmin -> Dot*)]]]&, xy, 7]
      ] /. specsu2 -> specsir;

xx = DiracSimplify[xx];
xx = Collect2[xx, Spinor, Factoring->False];
xx = Collect2[specsu[xx]/.specsu2 -> specsir, Spinor,Factoring->False];
If[Head[xx] =!= Plus,
   res = simP[xx],
   For[iji = 1, iji <= Length[xx], iji++,
      print2["iji = ",iji, " out of ", Length[xx]];
      res = res + (nP = simP[xx[[iji]]]);
       print2["nP = ",nP//FeynCalcForm];
     ];
  ];
];
res];
(* sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss *)

If[(!FreeQ[ts, Spinor]) && (enmo=!= {}),
   print2["minimizing SME's"];
   ts = spmin[ts, enmo];
print2["colling2insquareamplitudeagain"];
    ts = (Collect2[ts, {Spinor, SUNIndex} ] + nuLlLL) /. nuLlLL->0;
  ];
print1["ts//Length = ",ts//Length];


ts = ((#/.Plus->etl)& /@ (ts + nuLlLL)) /. nuLlLL->0;

If[StringQ[saveminamp],
   Write2 @@ {saveminamp, Global`$TSAMP == FRH[ts]}
  ];
];
(* ++++++++++++++++++++++++ *)

Factor3[x_] := Factor2[x, FactorTime -> 422];

If[ValueQ[Global`ENM], enmLI = Global`ENM, enmLI = {}];
(*
epsi[x_] := x /. Eps -> epsimP;
epsimP[x__] := EpsEvaluate[Eps @@ ({x} /. enmLI)];
*)

extrafact2 = Collect2[ extrafact2 , LorentzIndex];

(* for the non-fermionic case *)
contractP[a_,b_] := Block[{tem, tem1, tem2},
print2["entering contractP"];
                           tem = Contract2[a, b, Collecting -> False];
print2["exiting contractP"];
                     tem];

ijj=1;
(*
mulL[a_,a_]:=1;
mulL[a_,b_]:=2;
*)
mulL[__]:=1;
neamp = 0;
print2["ts[[1]] = ", ts[[1]]//FeynCalcForm];

(* Select w.r.t. to the structure of the amplitude *)
If[!FreeQ[ts, Spinor],

cts  = ComplexConjugate[ts//FRH] /. Plus -> plu3;
cts = ((#/.Plus->etl)& /@ (cts + nuLlLL)) /. nuLlLL->0;
print2["length of cts = ", cts//Length];

If[Head[ts] =!= Plus,
   neamp = FermionSpinSum[susm[ts cts],
                          ExtraFactor -> extrafact2] /. enmLI,
(* somehow very very weird ... *)
For[i = 1, i<=Length[ts], i++,
(*
    For[j = i, j <= Length[ts], j++,
 (* ] *)
*)
    For[j = 1, j <= Length[ts], j++,
print2["i = ",i, "  j = ",j, " out of ",Length[ts]];
        calc = mulL[i, j] (FermionSpinSum[susm[ts[[i]] cts[[j]]],
                                         ExtraFactor -> extrafact2] /. enmLI);
If[(epsAway===True) && ((Length[Position[calc, DiracGamma[5]]] +
                         Length[Position[calc, DiracGamma[6]]] +
                         Length[Position[calc, DiracGamma[7]]]
                       )) === 1 ,
   calc = calc /. Literal[x:Dot[DiracGamma[__]..]] :> 
                   (Dot[x] /. DiracGamma[7] -> (1/2) /.
                              DiracGamma[6] -> (1/2) /.
                              DiracGamma[5] -> 0
                   );
   calc = DiracSimplify[calc, Expanding -> False];
  ];
        If[calc===0&&lis>2, print2["this one vanishes ..."],
print2["LeafCount newtrace = ", LeafCount[calc]];
           neamp = neamp + calc;
           print2["Length of neamp = ",Length[neamp]];
          ];
       ]
   ]
  ],

pluIS[xx__] := If[FreeQ[{xx}, SUNIndex], pluIS2[xx], Plus[xx]];
(*
ts = FRH[ts] /. Plus -> pluIS;
*)
cts = ComplexConjugate[ts//FRH] /. Plus -> plu3;
If[Head[ts] === Plus,
   samp = susm[Expand[susa[ts cts /. Plus -> pluIS] ]],
   samp = susm[susa[ts cts/. Plus -> pluIS ]]
  ];
samp = samp /. pluIS2 -> Plus;
(* WWW *)
   
If[Head[samp] === Times,
   sampfa = Select[samp, FreeQ[#, LorentzIndex]&];
   print2["sampfa = ",sampfa];
   newsamp = samp/sampfa,
   sampfa = 1; newsamp = samp
  ];

print2["sampfa = ",sampfa];
print2["extrafact2 = ",extrafact2];
ti = Timing[
If[Head[newsamp]===Plus,
   samp = {};
   For[iis = 1, iis<=Length[newsamp], iis++,
       print2["contract samp ", iis," out of ",Length[newsamp]];
       nsamp = contractP[newsamp[[iis]], extrafact2];
       AppendTo[samp, nsamp];
      ];
   samp = Apply[Plus, samp]
   ,
   samp = contractP[newsamp, extrafact2];
  ];
   print2["Length of samp = ",samp//Length];
           ];
 neamp = samp sampfa;

If[$VeryVerbose >0, Print["time needed for polarization sums = ",
   ti[[1]]]//FeynCalcForm];

print1["    collecting now "];

];
   
print2["after expanding the length of the amplitude is ", 
       Length[neamp + neuladsg]-1];

If[FreeQ[neamp, DiracTrace], amps2 = neamp,
(*
If[Head[neamp] === Times,  
   pamp  = Select[neamp, !FreeQ2[#, {DiracGamma,Complex}]&];
   amps2 = neamp / pamp,   
    pamp = neamp; amps2 = 1
  ];
*)

SetOptions[DiracTrace, PairCollect -> True];
SetOptions[Tr, PairCollect -> True];
mom4[xyx_,___] := Momentum[xyx];

lastsimp[yyy_] := Block[{rel, mul = 1, yy = yyy},
                   If[Head[yy] === Times,
                      mul = Select[yy, FreeQ2[#,{Eps, LorentzIndex}]&];
                      yy = yy/mul
                     ];            
                       rel =  Contract[yy/.intermed, 
                                   EpsContract -> True] ;
                       If[!FreeQ[rel,Pair],
                          rel = rel /. Momentum  -> mom4; 
                          rel = ExpandScalarProduct[rel];
                         ];
print2["beforefactor3"];
                       If[LeafCount[rel]<9200, rel = Factor3[rel]];
print2["afterfactor3"];
                       If[epsAway===True, rel = rel /. Eps[__]->0];
                       If[rel ===0, print2["TRS00000000"],
                       rel = rel mul;
                       If[FreeQ[rel, Eps], 
                          rel = rel/.Plus->plu3,
                          rel = rel /. Plus-> etl
                         ];
                         ];
If[rel =!= 0, print2["TRESULT = ", rel//FeynCalcForm]];
                   rel];


TRS[xyz__] := memset[TRS[xyz] , Trc[xyz]];
epsimp[x__] := If[!FreeQ[{x},Eps], Schouten[Plus[x], 4444], Plus[x]];

If[!FreeQ[ neamp, DiracGamma],
   If[Head[neamp]=!=Plus,
      nwres = { SUNSimplify[neamp/.DiracTrace->TRS]},
      nwres = {};
      For[ij = 1, ij <= Length[neamp], ij++, print2["ij = ",ij, " out of ", Length[neamp](*, 
                  "|| length of nwres = ",Length[nwres]*)];
          AppendTo[nwres, (ww=lastsimp[neamp[[ij]] /. DiracTrace -> TRS])];
         ], 
      nwres = {SUNSimplify[neamp/.DiracTrace->TRS]}
     ];
  ];
print2["changing list to sum"];
amps2 = (Plus@@nwres);
];

If[factoring2 === True,  
   amps2 = Factor2[ExpandScalarProduct[
              (amps2//FRH) /. intermed] /. Momentum -> mom4],
   amps2 = Expand[ExpandScalarProduct[
              (amps2//FRH) /. intermed] /. Momentum -> mom4]
  ];
amps2 = ExpandScalarProduct[amps2/.intermed] /. finsubst2;
Print["amps = ",amps];
If[Length[mand2] ===4, amps2 = TrickMandelstam[amps2, mand2] ];

amps2
];

(* ************************************************************************ *)


propdenex[a_, m_] := propdenex[a, m] = 1/Factor2[
          FixedPoint[ExpandScalarProduct, ScalarProduct[a, a] - m^2] ];

Literal[propdp[su_ /; Length[su]<=2][a_, m_]] := propdenex[a,m];

propdp[su_][a_, m_] := 
 Block[{na = ExpandScalarProduct[Momentum[a]], i, vn, tte = {na}},
        vn = Variables[na];
tte = Join[tte,
       Table[ Expand[na /. (Solve[ su==0, vn[[i]] ][[1]])], {i, Length[vn]}]
          ];
tte = Sort[Union[tte], nsortQ][[1]];
propdenex[tte, m]
] /; Length[su] > 2;

fdpsave[x__] := fdpsave2 @@ ({x}/.PropagatorDenominator->pdsave);

PropagatorDenominatorExplicit[x_ /;FreeQ[x, PropagatorDenominator], ___] := x;
PropagatorDenominatorExplicit[x_] := x /. FeynAmpDenominator -> fdpsave /. 
       PropagatorDenominator -> propdenex /. pdsave -> PropagatorDenominator/.
          fdpsave2 -> FeynAmpDenominator;
PropagatorDenominatorExplicit[x_, {}] := PropagatorDenominatorExplicit[x];
PropagatorDenominatorExplicit[x_, {en__}] := x /. 
FeynAmpDenominator -> fdpsave /.
PropagatorDenominator -> propdp[Plus @@ ExpandScalarProduct[ Momentum /@ {en}]
                               ] /. pdsave -> PropagatorDenominator/.
                                     fdpsave2 -> FeynAmpDenominator;

(* FermionSpinSumdef *)
 If[$Context === "Global`", dotLin[x__] := Feyn`Calc`Main`dotLin[x]];
 Literal[trsimp[a_. DiracGamma[_,___]]] := 0 /; FreeQ[a, DiracGamma];
 Literal[trsimp[Dot[d__]]] := Trc[Dot[d] ] /; Length[{d}] < 4;
Options[FermionSpinSum] = {SpinPolarizationSum -> Identity,
                           SpinorCollect -> False,
                           ExtraFactor -> 1};
FermionSpinSum[x_,ops___] := Block[{spsf,spir,spir2,dirtri, nx,nnx, is=1, sufu,exf,
                                     plsp,lis, cOL, spinorCollect},
            nx = x;
            If[!FreeQ[x, Spinor], 
               spsf = SpinPolarizationSum /. {ops} /. Options[FermionSpinSum];
               exf = ExtraFactor/. {ops} /. Options[FermionSpinSum];
               spinorCollect= SpinorCollect/. {ops} /. Options[FermionSpinSum];

(* FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF *)
(* fermion polarization sums *)
   spir = { (* ubar u , vbar v *)
            Spinor[s_. Momentum[pe1_], arg__ ]^2 :>
           (dotLin[spsf[ (DiracSlash[pe1] + s First[{arg}]) ] ]),
           (Spinor[s_. Momentum[pe1_], arg__] . dots___ ) *
           (dots2___ . Spinor[s_. Momentum[pe1_], arg__ ] )  :>
            dots2 . dotLin[spsf[(DiracSlash[pe1] + s First[{arg}])]] . dots
          };

   spir2 = Spinor[s_. Momentum[pe_], arg__] . dots___ .
           Spinor[s_. Momentum[pe_], arg__] :> DiracTrace[(
           dotLin[spsf[(DiracSlash[pe] + s First[{arg}])]] . dots
                                               )] /; FreeQ[{dots}, Spinor] ;

   dirtri = DiracTrace[n_. a_Dot] DiracTrace[m_. b_Dot] :>
             DiracTrace[ DiracTrace[n a] m b] /; Length[a] <= Length[b] &&
               Head[n] =!= Dot && Head[m] =!= Dot;
(* FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF *)
uNi=Unique[System`C];
kK = Unique[System`C];
(*
uNi=Unique[cCC];
kK = Unique[ckC];
*)
pliui[xxx__] := pliui[xxx] = If[Length[{xxx}] < 5,
                                Isolate[Factor2[Plus[xxx]]],
                                Isolate[Plus[xxx], IsolateHead -> kK,
                                IsolateSplit -> 444I]
                               ];
$isoFlag = True;

cOL[xy_] := Block[{temP = xy, nodot = 0, ntemP},
print3["entering cOL"];
temp0 = temP;
                   If[Head[temP] === Plus,
                      nodot = Select[temP, FreeQ[#, Dot]&];
                      temP = temP - nodot;
                      nodot = nodot/. {a_ DiracGamma[5] :> 0 /; 
                               FreeQ[a, DiracGamma]};
                     ];
(*
If[$isoFlag,
                   If[Head[temP] === Plus,
                      temP = Map[#/.Plus -> pliui&, temP],
                      temP = temP /. Plus -> pliui;
                     ];
  ];
*)
temp1 = temP;
                   temP = Collect2[temP, DiracGamma, Factoring -> False 
                                  ];
print3["collected in cOL"];
(*
                   If[Head[temP] === Plus,
(* this step by step factorization is ESSENTIAL!!!!! *)
(* because this bloody Mma is not able to it directly ...*)
                      ntemP = 0;  lntemP = Length[temP];
                      For[ijn = 1, ijn <= Length[temP], ijn++,
                          print3["ijn = ",ijn,"(", lntemP,")"];
If[$isoFlag,
                          ntemP = ntemP + Factor2[FRH[Factor2[temP[[ijn]]]]]
           ,              ntemP = ntemP + Factor2[temP[[ijn]]]
  ];
            

                         ];
                      temP = ntemP ,
If[$isoFlag,
                      temP = Factor2[Factor2[temP]//FRH],
                      temP = Factor2[temP]
  ]
                     ];
*)
print3["exiting cOL"];
                   temP + nodot
                 ];

dirtracesep[xy_] := If[Head[xy] =!= Times, DiracTrace[xy//cOL],
    SUNSimplify[Select[xy, FreeQ2[#, {DiracGamma, LorentzIndex, Eps}]&]] * 
       DiracTrace[Select[xy,!FreeQ2[#, {DiracGamma, LorentzIndex, Eps}]&]//cOL]
                      ];

sufu[xyx_] := Block[{tsuf,spif,mulEx,mulEx2,epSimp,doT, xx=xyx, memm},
 print2[is++, "out of", lis, " Mem = [",memm = N[MemoryInUse[]/10^6,3],"]"];
If[(memm > $MemoryAvailable) &&IntegerQ[is/10], 
   print2["sharing "];Share[]; print2[" done"]];
                       spif = Select[xx, !FreeQ[#, Spinor]&];
                       tsuf = xx / spif;
(*
epSimp[xxx_] := If[FreeQ[xxx, Eps], xxx, DiracSimplify[xxx]];
*)
epSimp[xxx_] := DiracSimplify[DiracOrder[xxx] /. Dot -> doT /.
    {doT[a__, DiracGamma[5]] :> 0 /; Length[{a}] < 4,
     doT[DiracGamma[5]] :> 0,
     doT[a__DiracGamma] :> 0 /; FreeQ2[{a}, {DiracGamma[5], DiracGamma[6], 
              DiracGamma[7]}] && OddQ[Length[{a}]]
    } /. doT -> Dot, Expanding -> False    ];
                        
Literal[mulEx[mul_. DiracTrace[xy_]]] := 
  If[!FreeQ[(exf tsuf), LorentzIndex], 
                            dirtracesep[DiracSimplify[
     Contract[mul xy tsuf, exf], Expanding -> False, EpsContract ->
          False]//epSimp],
     dirtracesep[DiracSimplify[ mul exf tsuf xy,Expanding -> False]//epSimp]
    ];

mulEx2[ xy_ ] := mul tsuf exf xy;

(mulEx[(((spif)//.spir//.spir2//.dirtri) /. DiracTrace->trsimp/.
              trsimp->DiracTrace /. $MU->uNi )
      ] /.mulEx->mulEx2 
)                  ]; (* endofsufu *)

(*
 nx = Expand2[nx, Spinor];
*)
plsphold = Unique[System`C];
plsp[xyx__] := If[FreeQ[{xyx}, Spinor], plsphold[xyx], Plus[xyx]];

If[spinorCollect === True,
print2["collectinsufu"];
 nx = Collect2[nx /. Plus -> plsp, Spinor, Factoring -> False
              ] /. plsphold -> Plus;
print2["collectinsufudone"];
  ];

onx = nx;
               If[!FreeQ[ nx, $MU], nx = nx /. $MU->Unique[System`C]];
  If[Head[nx] === Plus,  
     lis = Length[nx];
     nnx = 0;
     For[iin = 1, iin <= lis , iin++, 
         nnx = nnx + sufu[nx[[iin]]];
        ];
     nx = nnx,
     lis = 1; nx = sufu[nx]
    ];

(* in case somthing went wrong .. *)
If[nx =!= 0 && FreeQ[nx, DiracTrace], Print[MIST];Dialog[]; nx = x exf];
              ] (* endIfFreeQ[x, Spinor]*);
(*
If[!FreeQ2[exf, {LorentzIndex, Eps}],  
    nx = Contract[ nx exf ], nx = nx exf ];
*)
nx];
           
                  
WriteString["stdout", "#"];
End[];

(* ******************************************************************* *)
(*                             oneloop.m                               *)
(* ******************************************************************* *)
Begin["Feyn`Calc`OneLoop`"];

If[$VersionNumber<2.0, 
   Print["Sorry, ... FeynCalc requires Mathematica 2.0 or 
          higher, bye bye ..."] ; Quit[] 
  ];


$LimitTo4 = True;
(* ********************************************************************* *)
(*                          oneloop10                                    *)
(* ********************************************************************* *)
   expandscalarproduct[x_]:=Feyn`Calc`Main`pairexpand[x];
   momentumexpand[x_]:=Feyn`Calc`Main`MomentumExpand[x];
   coneins[ x_ ] := Feyn`Calc`Main`memset[coneins[x], 
       x/.Pair->Feyn`Calc`Main`sCO/. Feyn`Calc`Main`sCO->Pair];
   diracgammacombine[x_]:=Feyn`Calc`Main`DiracGammaCombine[x];
   dotlin[x_] :=Feyn`Calc`Main`dotLin[x];
   print1[x__]:=Feyn`Calc`Main`print1[x];
   print2[x__]:=Feyn`Calc`Main`print2[x];
   print3[x__]:=Feyn`Calc`Main`print3[x];
(* contractlidef *)
   contractli[x_] := Contract[ x, Expanding->True, Factoring->False,
                                  EpsContract->False 
                             ] // Expand;
   conall[ x_ ] := Contract[ x,                               (*conalldef*)
                   Expanding->True, EpsContract->True, Factoring->False 
                           ] // Expand;

 collin[ expr_, varh_, expa_]:=Block[{nx,se,i,co,new, null,res,frex,totvarh},
      new = 0;
      nx = expr;
      If[ expa =!= True,
          frex[yy__]:=Plus[yy]/;!FreeQ2[{yy}, varh];
          nx = nx/.Plus -> frex
        ];
      nx = Expand[nx];
      totvarh = Select[ Variables[ nx ], (!FreeQ2[#,varh])&];
      If[(Length[nx] < 42) && (LeafCount[nx]<5000), 
         res = Collect[ nx, totvarh ],
          For[ i = 1, i <= Length[totvarh], i++,
               co = Coefficient[ nx, totvarh[[i]] ];
               nx = nx - Expand[ co totvarh[[i]] ];
               new = new + totvarh[[i]] co ];
        res = nx + new
        ];
      If[ expa =!= True, res = res/.frex -> Plus ];
      res];



   spinorchainevaluate[x_]:=DiracOrder[
      Expand2[Contract[ Feyn`Calc`Main`SpinorChainEvaluate[x] ], DiracGamma],{q}
                                      ];

(* ********************************************************************* *)
(*                          oneloop11                                    *)
(* ********************************************************************* *)


Options[OneLoop]={
                  DenominatorOrder           -> False,
                  Dimension                  -> D,
                  FinalSubstitutions         -> {}, 
                  Factoring                  -> False,
                  FormatType                 -> InputForm,
                  InitialSubstitutions       -> {},
                  IntermediateSubstitutions  -> {},
                  IsolateHead                -> False,
                  Mandelstam                 -> {},  
                  CancelQ2                   -> True,  
                  CancelQP                   -> False,
                  Prefactor                  -> 1,
                  ReduceGamma                -> False,
                  ReduceToScalars            -> False,
                  NegligibleVariables             -> {}, 
                  UVPartOnly                 -> False,
                  WriteOut                   -> True ,
                  WriteOutPaVe               -> False
                 };
(* setting WriteOut to "" retrieves also previously calculated results *)

(* OneLoopdef *)
OneLoop[qq_,amp_]:=OneLoop[False, qq,amp];

OneLoop[qq_,amp_,opt1_,opt___]:=OneLoop[False, qq,amp,opt1,opt]/;
				!FreeQ[opt1,Rule];

OneLoop[grname_,q_,integ_,opts___] := Block[
       {oneamp=integ,iv,onemandel,denorder,isolatehead,tric,
        var2,smallv,finsubst,fact,q2cancel,qpcancel,pvabbrev,
        writeout,prop,dnq,dfsor,dfsorb, denomOrder,dfli,
        vcid,intcan,de, oneoptga67,vin3,q2canc, sqB,usual,tostandmat,
        vint4,ret,vret,fmas,vva,smalist,isol,i, $higherpoint,
        pvlist, pva,pvar,arglist,npref, nfact, nre,formattype,pLu,
        collpav,simpit,prefactor,in3p, in3q, resin3q,newprefactor,
        newnewprefactor,
        newret, facto,ret2,defs,dim,breakdown,options, name=grname,DDim,
        newoneamp,ip,lenneu, lenneu2, neuamp,paone,paone2,oneselect,fsub,
        intermedsubst,writeoutrecover = False, oldfile, ftemp,
        writeoutpav, uvpart,to4dim,tone,ltone,holdoneamp,
         lssu
        },

options = {opts};
(* for FA2.0 *)
If[Length[options] > 0, 
   If[ MatchQ[options[[1]] ,(a_List -> b_List)],
       FA2Info = options[[1]];
       options = Rest[options]
     ]
  ];

( oneopt    = Join[ options,Options[OneLoop] ];
 onemandel  = Mandelstam/.oneopt;
 denorder   = DenominatorOrder/.oneopt;
 dim        = Dimension/.oneopt;
If[ dim===True, dim = D ];
 formattype = FormatType/.oneopt;
 isolatehead= IsolateHead/.oneopt;
 prefactor  =  Prefactor/.oneopt ;
 smallv     =  Flatten[{ NegligibleVariables/.oneopt }];
 inisubs    =  InitialSubstitutions/.oneopt;
 finsubst   = FinalSubstitutions/.oneopt;
 intermedsubst = IntermediateSubstitutions/.oneopt;
 fact       = Factoring/.oneopt;
 q2cancel   = CancelQ2/.oneopt;
 qpcancel   = CancelQP/.oneopt;
 writeoutpav = WriteOutPaVe /. oneopt;
 uvpart = UVPartOnly/.oneopt;


 reducegamma67 = ReduceGamma/.oneopt;
 writeout   = WriteOut/.oneopt;
If[writeout === True || writeout === " ", writeout = ""; writeoutrecover = False];
If[writeout === True || writeout === "" , writeoutrecover = True];
 breakdown  = ReduceToScalars/.oneopt;
 If[(breakdown===True) && ($LimitTo4 === False),
    SetOptions[PaVeReduce, Dimension -> dim]
   ];
);

If[uvpart === True, 
   breakdown = True; 
   $LimitTo4=True; 
   SetOptions[PaVeReduce, Dimension -> True];
   dim = 4
 ];

 If[(breakdown===True) && ( (writeoutpav===False) || (writeoutpav===True) ),
    writeoutpav = ""
   ];

(* print1def, print2def, print3def: print functions *)
 SetAttributes[{print1, print2, print3}, HoldAll];
 print1[x__]:=Print[x]/;$VeryVerbose>0;
 print2[x__]:=Print[x]/;$VeryVerbose>1;
 print3[x__]:=Print[x]/;$VeryVerbose>2;

tim = Timing[
If[ StringQ[ name ] && StringQ[writeout], name= StringJoin[writeout, name]];
   If[ Head[name]===GraphName, 
        name = StringJoin @@ (ToString/@{First[name], Last[name]}),
        If[(name=!=False) && (Head[name]=!=String), name=ToString[name]]
     ];
   If[ Head[name]===String,
       Which[ formattype === InputForm,   name = StringJoin[name, ".m"],
              formattype === FortranForm, name = StringJoin[name, ".for"],
              formattype === MacsymaForm, name = StringJoin[name, ".mac"],
              formattype === MapleForm,   name = StringJoin[name, ".map"]
     ]      ];

If[StringQ[name],
oldfile = FileNames @@ {name};
If[oldfile =!= {},
   print1["oldfile  =", oldfile];
   ftemp =( Get @@ {name} );
   If[ValueQ[OneLoopResult[grname]] && FreeQ[ftemp, FeynAmpDenominator],
      oneamp = OneLoopResult[grname]
     ];
  ]
];

oneampstart = oneamp;

(* in case oneamp has no FeynAmpDenominator: write oneamp out *)
If[ FreeQ[ oneamp , FeynAmpDenominator], 
    oneampresult = oneamp,

(* A0[0] = 0 *)

FeynAmpDenominator[PropagatorDenominator[Momentum[q], 0]]:=0;

(* ********************************************************************* *)
(*                          oneloop12                                    *)
(* ********************************************************************* *)


If[ q2cancel === True, q2cancel = 1 ];
(* * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(* Starting  the game:  *)
(* * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * *)

(*epschisholmdef*)
   epschisholm[x_]:=x/;FreeQ[x, Eps] || FreeQ[x, DiracGamma];
   epschisholm[x_Plus]:=Map[epschisholm,x];
   epschisholm[x_]:=If[reducegamma67===True, 
                       EpsChisholm[x],
                       EpsChisholm[x]/.DiracGamma[5] ->
                                      (DiracGamma[6]-DiracGamma[7])
                      ] /; Head[x]=!=Plus;


   If[ FreeQ[oneamp, Spinor[p1__] . a__ . Spinor[p2__] *
                     Spinor[p3__] . b__ . Spinor[p4__]
            ], Feyn`Calc`Main`$sirlin = False
     ];
print3["$sirlin = ", Feyn`Calc`Main`$sirlin];

pri[iii_]:=print["check ",iii,"  ",oneamp//FeynCalcForm];
(* smallv *)
   smav=Table[smallv[[iv]]->Negligible[ smallv[[iv]] ], 
              {iv,1,Length[smallv]} ];
(* do the initial substitutions *)
   oneamp = ( oneamp/.smav )/.inisubs;
(* neglect any small variable in the numerators of the fermion propagators  *)
   smalldirac /: smalldirac[_] + DiracGamma[a__] := DiracGamma[a];
(* and in the spinors *)
   smalldirac /: Spinor[ pe_, smalldirac[_], op___] := Spinor[pe,0,op];
   oneamp = oneamp /. Negligible -> smalldirac /. smalldirac -> Negligible;
(* put heads on the momenta in the denominators *)
   denf[x_,y_]:=denprop[Momentum[x],y]/;FreeQ[x,Momentum];
   denf[x_,y_]:=denprop[x,y]/;!FreeQ[x,Momentum];
(* if a propagator has no integration momentum *)
   denprop[a_, b_] := (1/TrickMandelstam[(Pair[a,a]//ExpandScalarProduct)
                        - b^2, onemandel] /. Negligible[_]->0)/; FreeQ[a, q];
(* ********************************************************************* *)
(*                          oneloop13                                    *)
(* ********************************************************************* *)

(* ONEAMPCHANGE: denominators *)

 If[ Cases[x, Literal[Spinor[a_,_,_] . (___) . Spinor[b_,_,_] * 
                       Spinor[c_,_,_] . (___) . Spinor[d_,_,_] ] 
          ] =!= {},
     $fourfermion = True
   ];
oneamp00 = oneamp;
   oneamp = oneamp/.PropagatorDenominator -> denf/.denprop -> 
                                             PropagatorDenominator;

(*
(* do a partial fraction decomposition of propagators with the same 
   momenta (sometimes necessary for general gauge) *)
Literal[
feynampdenpartfrac[ a___, PropagatorDenominator[ qpe1_, 0], b___,
                          PropagatorDenominator[ qpe1_, 0], c___
       ]          ] := ( 1/(Global`$ZeroMass1^2 - Global`$ZeroMass2^2) feynampdenpartfrac[a, 
                          PropagatorDenominator[ qpe1, Global`$ZeroMass1], b, c] -
                         1/(Global`$ZeroMass1^2 -  Global`$ZeroMass2^2) feynampdenpartfrac[a,  
                          PropagatorDenominator[ qpe1,  Global`$ZeroMass2], b, c]                         
                       );
*)
Literal[
feynampdenpartfrac[ a___, PropagatorDenominator[ qpe1_, ma1_ ], b___,
                          PropagatorDenominator[ qpe1_, ma2_ ], c___
       ]          ] := ( 1/(ma1^2 - ma2^2) feynampdenpartfrac[a, 
                          PropagatorDenominator[ qpe1, ma1 ], b, c] -
                         1/(ma1^2 - ma2^2) feynampdenpartfrac[a,  
                          PropagatorDenominator[ qpe1, ma2 ], b, c]                         
                       ) /; (ma1^2 - ma2^2)=!=0;

print2["length of oneamp = ", Length[oneamp]];
   oneamp = Expand[ oneamp /. FeynAmpDenominator -> feynampdenpartfrac /.
                              feynampdenpartfrac -> FeynAmpDenominator,
                    FeynAmpDenominator ];


(* reduce N-point to 4-point,  still temporary  *)
SetAttributes[ Y, Orderless ];
smn[a_]:=0;
smdaemon[x_]:=x/.Negligible->smn;
Y[{pi_, mi_}, {pj_, mj_}]:=Expand[smdaemon[ mi^2 + mj^2 - 
                                   ScalarProduct[pi - pj, pi - pj]//
                                   ExpandScalarProduct
                                          ]
                                 ];

(* pli is the list of momenta (without q), being one less than mli *)

(* ********************************************************************* *)
(*                          oneloop14                                    *)
(* ********************************************************************* *)
(* Reducing N-=point to 4 point *)
(* ********************************************************************* *)

gr454[mli_List, pli_List]:= gr454[mli, pli] = 
   Block[{pl=Prepend[pli, 0], n = Length[mli]},
   (* create a list of the columns of eq. (4.54) *) 
 columnli = Prepend[ Table[ Y[{pl[[i]], mli[[i]]}, {pl[[j]], mli[[j]]}], 
                            {j,1,n}, {i,1,n}],
                     Array[#^0&, n]
                   ];
    columnli];

(*define this such that sgr[{...},{...}, -1  ] is the Y_{ij} determinant *)
sgr[mli_, pli_, i_] := sgr[mli, pli, i] = 
                        Factor2[Drop[gr454[mli, pli], {i+2,i+2}]//Det];

redamp[x_Plus,qu_]:=redamp[#,qu]&/@x;

redamp[x_,qu_] := Block[{pl, ml, fm, pp,res=x},
                         fd = PartitHead[x, FeynAmpDenominator];
                         If[ Length[fd[[2]]] > 4, 
                             fm[pe_,___]:=pe;
                             pp[a_, b_]:= Expand[(a/.Momentum->fm) - qu];
                             mm[a_, b_]:= b;
                         pl = List @@ ( fd[[2]]/.PropagatorDenominator->pp );
                         pl = Rest[pl];
                         ml =  List @@ ( fd[[2]]/.PropagatorDenominator->mm );
(*Extract the determinant *)
                        res = {1/sgr[ml, pl, -1] ,
                                Sum[ (-1)^j subdethold[sgr[ml,pl, j]] fd[[1]] *
                                     Drop[fd[[2]], {j+1,j+1}], 
                                     {j, 0, Length[fd[[2]]]-1}
                                   ]
                               }
                           ];
                         res] /; Head[x] =!= Plus;
    subdethold[x_Plus]:=x/;Length[x] < 4;
    subdethold[x_]:=ReleaseHold[Isolate[x, IsolateSplit->Infinity,
                                   IsolateHead->SUB
                                ] /. SUB -> SUBDET
                               ];

   $higherpoint = False;
   fdhigh[xx__] := If[Length[Union[{xx}]] < 5, 0, FeynAmpDenominator[xx]];
   If[ !FreeQ[oneamp /. FeynAmpDenominator -> fdhigh, FeynAmpDenominator],
      $higherpoint = True;
      namp = redamp[ oneamp,q ];
      prefactor = prefactor namp[[1]];
      oneamp = namp[[2]]
     ];

(* Put here the i pi^2 from the integrals *)
   newprefactor = prefactor I Pi^2;
(* ONEAMPCHANGE: extract coupling constants *)
If[Head[oneamp] === Times, 
   oneselect = Select[oneamp,
                FreeQ2[#, {Pair, PropagatorDenominator,Eps,
                           dim,
                           Momentum,LorentzIndex, SU3f,SUNF,
                           SU3delta, GellMannMatrix,
                           DiracGamma, Spinor}
                      ]&
                     ];
   oneamp = oneamp/oneselect;
   newprefactor =  Factor2[ newprefactor oneselect ];
   ];
 If[(!FreeQ[ newprefactor, dim ]) || (!FreeQ[newprefactor, LorentzIndex])
     || (!FreeQ[newprefactor, DiracGamma]),
    newnewprefactor = Select[newprefactor, 
                             !FreeQ2[#, {dim, LorentzIndex, DiracGamma}]&];
    If[!FreeQ[newnewprefactor, DiracGamma], 
       oneamp = newnewprefactor . oneamp,
       oneamp = oneamp newnewprefactor
      ];
    newprefactor = smalld[newprefactor / newnewprefactor + nUUUl
                         ]/. nUUUl -> 0;
    If[newprefactor === 0, oneamp = 0];
   ];

(* ********************************************************************* *)
(*                          oneloop15                                    *)
(* ********************************************************************* *)


(* change the dimension ,   to4dimdef*)
(* THIS IS TOTALLY WRONG !!!!!!!!!!!!!!!
   to4dim[x_]:=x/.{Momentum[v_,_]:>Momentum[v],
                   LorentzIndex[w_,_]:>LorentzIndex[w]} /. dim -> 4;
*)
(* therefore back again *)

   to4dim[x_]:=x/.{Momentum[v_,_]:>Momentum[v],
                   LorentzIndex[w_,_]:>LorentzIndex[w]};
oneamp0 = oneamp;
print2["check"];
(* ONEAMPCHANGE: make dimensions right *)
   If[ dim===4, oneamp = SUNSimplify[to4dim[ oneamp ], SUNFToTraces -> False],
       oneamp = oneamp + null;
       newone = 0;
       For[ i=1, i<=Length[oneamp], i++,
            newone = newone + SUNSimplify[
                     (oneamp[[i]]/.
            {Pair[LorentzIndex[a_],LorentzIndex[b_]]:>
             Pair[LorentzIndex[a,dim],LorentzIndex[b,dim]],
             Pair[LorentzIndex[a_],Momentum[p_]]:>
             Pair[LorentzIndex[a,dim],Momentum[p,dim] ]/;!FreeQ[p,q],
             Pair[Momentum[qu1_],Momentum[qu2_]]:>
             Pair[Momentum[qu1,dim],Momentum[qu2,dim]],
             DiracGamma[LorentzIndex[a_]]:>DiracGamma[LorentzIndex[a,dim],dim],
             dirga_[a_]:>DiracGamma[a/.Momentum[pe_]->Momentum[pe,dim],dim]/;
                         (dirga===DiracGamma) && !FreeQ[a,Momentum],
             PropagatorDenominator[mom_, mas_]:>
             PropagatorDenominator[mom/.Momentum[ppe_]->Momentum[ppe,dim],mas]
            }),     SUNFToTraces -> False]
          ];
       oneamp = newone/.null -> 0
     ];
   print3[" oneamp = ", oneamp];
If[(!FreeQ[oneamp, SUNF]) || (!FreeQ[oneamp, SUNDelta]) || 
    (!FreeQ[oneamp, SUNT]),
   oneamp = oneamp /. SUNF -> sUNF /. SUNDelta -> sUNDelta;
   AppendTo[finsubst, {sUNF -> SUNF, sUNDelta -> SUNDelta}];
  ];

(* ********************************************************************* *)
(*                          oneloop16                                    *)
(* ********************************************************************* *)

(* in order to get scalar products in the right dimension *)
   Feyn`Calc`Main`scevdoit[ Momentum[x_,_], Momentum[y_,_] ] := 
   Pair[Momentum[x], Momentum[y]] /; FreeQ[{x, y}, q];
(* for propagators without the integration variable q *)
    dnq[a___,PropagatorDenominator[pe_,ma_],b___] := 
   expandscalarproduct[ 
   (1/Factor2[ smalld[TrickMandelstam[ expandscalarproduct[Pair[pe,pe]]-ma^2,
                    onemandel ]//Expand] ] ) dnq[a,b] ] /; FreeQ[pe,q];
   prode[a_,b_] := PropagatorDenominator[a//momentumexpand,b] /; !FreeQ[a,q];
   prode[ppe_,mm_]:=
   expandscalarproduct[
   (1/Factor2[ smalld[TrickMandelstam[ expandscalarproduct[Pair[ppe,ppe]]-mm^2,
                    onemandel ]//Expand] ] ) ]  /; FreeQ[ppe,q];

(* ONEAMPCHANGE: extract denominators without q's *)
   oneamp = oneamp/.FeynAmpDenominator->dnq/.dnq->FeynAmpDenominator/.
                    PropagatorDenominator->prode;

   dfsor[ve_,ma_]:=defs[ma][ve]//momentumexpand;
   dfsorb[a_][b_]:=PropagatorDenominator[b,a];
(* denomOrder orders the PropagatorDenominator         , denomOrderdef *)
(* the eventually necessary translation of q is done with intcan.       *)
   denomOrder[ a__ ]:=Block[{ dfli={a},ru },
(* this is checked always *)
     ru={PropagatorDenominator[c_. Momentum[ce_. q ,di___],ma0_]:>
         PropagatorDenominator[Sqrt[c^2] Momentum[Sqrt[ce^2] q,di],ma0],
         PropagatorDenominator[x_+be_ Momentum[q,di___],mA_] :> 
         PropagatorDenominator[Expand[-x+Momentum[q,di] ],mA]/;(be===-1)
        };
     dfli = dfli//.ru;
     If[!MatchQ[ First[dfli], PropagatorDenominator[Momentum[q,___],_] ],
        denorder = True ];
     If[ denorder === True,
         dfli = Sort[ dfli/.PropagatorDenominator->dfsor ]/.defs->dfsorb;
(* For boxes: bring a evtl. 0 at position 3 *)
     If[ Length[ dfli ] === 4, dfli = RotateRight[dfli,2] ] ];
     If[ Length[dfli]===3 && denorder === True,
         dfli = dfli/.{den1_, PropagatorDenominator[pe2_,ma1_],
                              PropagatorDenominator[pe3_,ma2_]}:>
                      {den1,  PropagatorDenominator[pe3,ma2],
                              PropagatorDenominator[pe2,ma1]}/;
                      Length[ pe2 ] > Length[ pe3 ]
       ];
     If[ denorder === True,
         dfli = dfli/.{den1_,PropagatorDenominator[pe2_,ma2_],den3_,
                             PropagatorDenominator[pe3_,ma2_]}:>
                      {den1, PropagatorDenominator[pe3,ma2],den3,
                             PropagatorDenominator[pe2,ma2]}/;
                             Length[pe2]>Length[pe3];
      print3["after denomOrdering : ",dfli];
       ];
      FeynAmpDenominator@@dfli] (* endBlock *);
(* ********************************************************************* *)
(*                          oneloop17                                    *)
(* ********************************************************************* *)

   vcid[pe_,___]:=pe;
(* for translating the integration variable q in the first propagator *)
(* intcandef *)
   intcan[x_]:=x/;FreeQ[x,q];
   intcan[x_Plus]:=intcan/@x;
   intcan[ a_ b_ ]:=a intcan[b]/; FreeQ[a,q];
   intcan[any_. FeynAmpDenominator[ 
               PropagatorDenominator[p_ + Momentum[q,dim], m1_],dfrest___ ] 
         ]:=denomExpand[ ( any FeynAmpDenominator[ 
                 PropagatorDenominator[ p + Momentum[q,dim],m1], dfrest] 
                         )/.q->( q - (p/.Momentum->vcid) )
                       ];
(*  bringing the denominator in a canonical form                         *)

(* ONEAMPCHANGE : fixing all denominators , evtl. ordering *)
   oneamp = denomExpand[ oneamp ];
       oneamp = denomExpand[ oneamp ]/.FeynAmpDenominator->denomOrder;
       oneamp = (intcan[oneamp//momentumexpand]//momentumexpand)/.
                intcan->Identity ;
   
   oneamp = Map[ (Numerator[#]/Factor2[Denominator[#]])&, oneamp + null ];
   oneamp = oneamp/.null->0;
   print3["oneamp after ordering = ",oneamp];

(* ONEAMPCHANGE : bringing all denominator in standard order *)
   consum[x_]:=conall[x]/;Head[x]=!=Plus;    (* consumdef *)
   consum[x_Plus]:=Block[{nx=0,i},
                         For[i=1, i<=Length[x], i++,
                             print2["contracting # ",i," out of ",Length[x]];
                             nx = nx + conall[x[[i]]] 
                            ];
                         nx];
  (* consum2def *)
   consum2[x_]:=smalld[expandscalarproduct[conall[x]]]/;Head[x]=!=Plus; 
   consum2[x_Plus]:=Block[{nx=0,i},
                         For[i=1, i<=Length[x], i++,
                             print2["contracting # ",i," out of ",Length[x]];
                             nx = nx + (Expand[
                                 expandscalarproduct[conall[x[[i]]]]
                                             ]//smalld)
                            ];
                      nx];
(* ********************************************************************* *)
(*                          oneloop18                                    *)
(* ********************************************************************* *)

   oneamp = Contract[ oneamp, Expanding -> False];
oneampc = oneamp;
   oneamp = DiracSimplify[oneamp, Expanding -> False];
oneampd = oneamp;
   oneamp=collin[ consum[oneamp]//smalld, FeynAmpDenominator, False ];
(* ONEAMPCHANGE : contracting Lorentz indices and expanding *)
(* now a canonized form is achieved *)
(* and all Lorentzindices which are not part of a DiracGamma contracted *)

   If[!FreeQ[oneamp,DiracGamma],print1["Simplification of Dirac structures"]];
   If[!FreeQ[oneamp,DiracTrace], print1["and calculation of Traces"] ];

(* ONEAMPCHANGE : contracting Lorentz indices and trace calculation *)
   
   If[ !FreeQ[oneamp, DiracTrace],
       neuamp = 0;
       oneamp = collin[ oneamp, DiracTrace, False];
       If[ Head[oneamp] =!= Plus, oneamp = EvaluateDiracTrace[oneamp],
           For[i=1, i<=Length[oneamp], i++,
               print1["calculating trace # ",i,"  out of ",Length[oneamp]];
timi = Timing[
               neuamp = neuamp + Expand[ oneamp[[i]]//EvaluateDiracTrace ]
             ][[1]];
 print1["time needed = ",timi];
              ];
           oneamp = Collect2[neuamp, FeynAmpDenominator, Factoring -> False]
         ]
     ];

(* ONEAMPCHANGE : bringing gamma5, gamma6 and gamma7 into standard way according to 
                  the setting of the option ReduceGamma
*)
   If[reducegamma67=!=True,
      oneamp = oneamp /.(* DiracGamma[5]->(DiracGamma[6]-DiracGamma[7])/.*)
                 {Literal[ Spinor[p1__].(a___).Spinor[p2__]] :>
                   (Spinor[p1].a.DiracGamma[6].Spinor[p2] +
                    Spinor[p1].a.DiracGamma[7].Spinor[p2]
                   ) 
                 };
    ];

   If[reducegamma67===True,
      oneamp = oneamp /. DiracGamma[6] -> (1/2 + DiracGamma[5]/2)/.
                         DiracGamma[7] -> (1/2 - DiracGamma[5]/2)
     ];
(* ********************************************************************* *)
(*                          oneloop19                                    *)
(* ********************************************************************* *)

   print1["contraction, etc. "];
   diracorder[x_]:= x /; FreeQ[ x, DiracGamma ];
   diracorder[x_]:=Collect2[DiracOrder[x(*, {q}*)], Dot,Factoring->False
                           ] /; Head[x]=!=Plus;
   diracorder[x_Plus]:=Collect2[ Map[DiracOrder[#(*, {q}*)]&, x], Dot,
                                  Factoring->False];

   $Test = True;
   paex[x_,y_]:=If[!FreeQ[{x,y}, Eps], 
                   ExpandScalarProduct[Pair[x,y]],
                   Pair[x,y]
                  ];

   epsit[x_]:=If[$Test === True, 
                 epschisholm[Collect2[x/.Pair->paex, {Eps,Spinor}, 
                                                Factoring -> False ]], 
                 x
                ];
   simpit[x_ /; FreeQ[x, DiracGamma]] := (print3["checkkkk"];consum2[x]);
   simpit[x_] := diracorder[ Collect2[ consum2[x]//epsit, Dot, 
                                        Factoring->False
                                     ] // DiracSimplify
                           ] /; !FreeQ[x, DiracGamma];
    
(* we want to keep the distinction in different graphs *)
(* ONEAMPCHANGE *)
oneamp1 = oneamp; print3["oneamp1 = ",oneamp1];
   neuamp = 0;
   If[ Head[oneamp] === Plus, lenneu = Length[oneamp], lenneu = 1];
       For[ ip=1, ip <= lenneu, ip++,
            print1["working with part # ",ip," out of ",lenneu];
            If[ lenneu > 1, 
                parf = PartitHead[oneamp[[ip]], FeynAmpDenominator]//smalld,
                parf = PartitHead[oneamp, FeynAmpDenominator]//smalld
              ];
            neuamp = neuamp + ((Expand[ parf[[2]] simpit[parf[[1]]] 
                                      ]) //smalld)
          ];
print2["Length of neuamp = ", Length[neuamp]];
If[MemoryInUse[]> $MemoryAvailable 10^6,
timsh = Timing[
Share[]
             ][[1]];
print1["time needed for Share = ",timsh];

print1[N[ MaxMemoryUsed[]/10^6,2 ], "MB used"];
  ];

If[ FreeQ[neuamp, Dot],  oneamp = neuamp,
   oneamp = Collect2[neuamp, Dot, Factoring->False];
  ];
oneamp2 = oneamp; print3["oneamp2 = ",oneamp2];

(* ONEAMPCHANGE *)
   If[ intermedsubst =!= {},
       oneamp = oneamp /. intermedsubst /. intermedsubst;
          neuamp = 0;
   If[ Head[oneamp] === Plus, lenneu = Length[oneamp], lenneu = 1];
       For[ iip=1, iip <= lenneu, iip++,
            print1["working with (substituted) part # ",iip," out of ",lenneu];
            If[ lenneu > 1,
                parf = PartitHead[oneamp[[iip]], FeynAmpDenominator],
                parf = PartitHead[oneamp, FeynAmpDenominator]
              ];
            neuamp = neuamp + ((Expand[ parf[[2]] simpit[parf[[1]]]
                                      ]) //smalld)
          ];
      ];
If[MemoryInUse[]>10^7,
timsh = Timing[
Share[]
             ][[1]];
print1["time needed for Share = ",timsh];

print1[N[ MaxMemoryUsed[]/10^6,2 ], "MB used"];
  ];


   oneamp = neuamp;
       
   print1["check1"];
   print3["neuamp = ",neuamp];

(* Zwischenspiel... *)
  smadot[]=1;
  standard /: standard[a_] standard[b_] := standard[a b];
  smadot[x___]:=standard[ dotlin[Dot[x]]/.dim->4 ] /.  
                standard -> StandardMatrixElement;
  If[StandardMatrixElement =!= Identity,
     StandardMatrixElement[x_Plus]:=Map[StandardMatrixElement,x]
    ];
  
(* ********************************************************************* *)
(*                          oneloop20                                    *)
(* ********************************************************************* *)

(* tostandmatdef *)
  tostandmat[xy_]:=Block[{te, standa},
                          standa[a__]:= dotlin[Dot[a]] /; 
                                           !FreeQ2[{a}, {Dot,Pair,Spinor}] ;
                          standa[a_,b__]:=StandardMatrixElement[a,b]/;
                                            FreeQ2[a, {Dot,Pair,Spinor}];
			  te = tempstandmat[Expand[xy]]//Expand;
			  If[ !FreeQ[te, DiracGamma[6]],
			      te = te/.StandardMatrixElement -> standa/.
                                   standa -> spinorsandpairs/.
                                   spinorsandpairs->smadot/.
                                   DiracGamma[6] -> 
				  (1/2 + DiracGamma[5]/2);
			      te = spinorchainevaluate[te]//Expand;
			      te = tempstandmat[te]//Expand
                            ];
			  If[ !FreeQ[te, DiracGamma[7]],
			      te = te/.StandardMatrixElement -> standa/.
                                   standa -> spinorsandpairs/.
                                   spinorsandpairs->smadot/.
                                   DiracGamma[7] -> 
				  (1/2 - DiracGamma[5]/2);
			      te = DiracSimplify[te]//DiracOrder//Expand;
			      te = tempstandmat[te//Contract]//Expand
                            ];
                          te = te /.Dot->spinorsandpairs/.
                                   spinorsandpairs->smadot;
                      te] /; reducegamma67 === True;
            
    tostandmat[xy_]:=Block[{te, standa, pluh, pluh2},
                          standa[xx_]:=xx/;!FreeQ2[xx,{Dot,Pair,Polarization}];
                          standa[a_,b__]:=StandardMatrixElement[a,b];
                          pluh[x__]:=Plus[x] /; !FreeQ2[{x}, 
                          {Spinor,GellMannMatrix, SUNIndex, Polarization,Pair}];
                          te = xy /. Plus -> pluh /. pluh->pluh2;
                          te = tempstandmat[Expand[te]]//Expand;
                          te = te /. pluh2 -> Plus;
                              te = te/.StandardMatrixElement -> standa/.
                                   DiracGamma[5] ->
                                  (DiracGamma[6] - DiracGamma[7]);
			      te = te//dotlin;
    (* 1 = gamma6 + gamma7 *)
                              If[!FreeQ[te, Spinor],
                                 te = te //. {
                                 Literal[
                                 Spinor[p1__].(a___).Spinor[p2__] ]:>
                                 (Spinor[p1].a.DiracGamma[6].Spinor[p2] +
                                  Spinor[p1].a.DiracGamma[7].Spinor[p2]
                                 ) /; FreeQ2[{a}, {DiracGamma[6],
                                                   DiracGamma[7]} ]
                                                  };
                                 te = Expand[DiracSimplify[te]//DiracOrder,
                                             Spinor]//Contract;
                                 te = Expand2[ te, Pair ]
                                ];
                              te = tempstandmat[te];
                         te = te /.Dot->spinorsandpairs/.       
                                     spinorsandpairs->smadot;
                      te] /; reducegamma67 =!= True;

(* ********************************************************************* *)
(*                          oneloop21                                    *)
(* ********************************************************************* *)

(* ONEAMPCHANGE : spinor stuff and matrixelements *)
   print1[" Dirac-Algebra again"];
   print2["before spinorch: oneamp = ",oneamp//Length];
(*
   oneamp =  conall[spinorchainevaluate[oneamp//smalld ]]//
            expandscalarproduct;
*)
   If[ Head[oneamp]===Plus,
       print1["substituting"];
       If[ Length[oneamp]<10, 
            oneamp = oneamp/.{DiracGamma[Momentum[pe_,di_Symbol],di_Symbol]:>
                              DiracGamma[Momentum[pe]]/;FreeQ[pe,q] };
               oneamp =  conall[oneamp//smalld ]// expandscalarproduct,
            newamp = 0;
            For[ij = 1, ij <= Length[oneamp], ij++,
                If[ IntegerQ[ij/100], print2["ij = ",ij] ];
                temp =  (oneamp[[ij]]/.
                  {DiracGamma[Momentum[pe_,di_Symbol],di_Symbol]:>
                   DiracGamma[Momentum[pe]]/;FreeQ[pe,q] });
                temp = conall[temp//smalld] // expandscalarproduct;
                newamp = newamp + Expand2[temp, q]
               ];
            oneamp = newamp
         ];
     ];
              
(*
   oneamp = Collect2[oneamp/.subdethold->Identity, q, Factoring-> False];
*)
print1["collect w.r.t. ", q];
   oneamp = Collect2[oneamp, q, Factoring -> False];

        q2canc[aq_Plus]:=Block[{nq2, iq, la},
             nq2 = 0; la = Length[aq];
             For[iq=1, iq<=la, iq++, 
                 print2["cancel q2 : iq= ",iq,"  / ",la];
                 nq2 +=Expand[ q2canc[aq[[iq]]] ]//smalld
                ];         nq2]; 
   q2canc[x_]:=x/;FreeQ[x,q];
   q2canc[]=1;
   q2canc[bbb_]:=( conall[ expandscalarproduct[ Q2Cancel[bbb,q2cancel,q] 
                       ]           ]
               )/;!(Head[bbb]===Plus);

(* ONEAMPCHANGE : cancelling q^2 *)
   If[ IntegerQ[q2cancel] && 
       !FreeQ[  oneamp, Pair[ Momentum[q,___],Momentum[q,___] ]  ],
       print1["cancelling q^2"];
       If[ Head[oneamp]===Plus,
           in3p = Select[ oneamp, FreeQ[#, Pair[ Momentum[q,___],
                                              Momentum[q,___]] ]& 
                        ];
           in3q =  Collect2[ oneamp - in3p, {Pair,q}, Factoring -> False]; 
	   resin3q = q2canc[ in3q ]; oneamp = resin3q + in3p,
           oneamp = q2canc[ oneamp ]
         ] 
     ];
(* ********************************************************************* *)
(*                          oneloop22                                    *)
(* ********************************************************************* *)

(* ONEAMPCHANGE : cancelling q p 's *)
  If[(qpcancel === True) && (oneamp =!= 0),
print1["cancelling qp's"];
     qpcanc0[aa_ b_]:=(aa qpcanc0[b])/;FreeQ[aa,q];
     qpcanc[]=1;
     qpcanc[aq_Plus]:=Block[{nq2=0, iq, la = Length[aq]},
             For[iq=1, iq<=la, iq++, 
                print2["cancel qp : iq= ",iq,"  / ",la];

tocancel =  qpcanc0[aq[[iq]]];	
print3["tocancel = ",tocancel];
                 nq2 +=Expand[ qpcanc0[aq[[iq]]]/.qpcanc0->qpcanc ]//smalld;
                ];         nq2];
     qpcanc[b_] := QPCancel[b,Momentum[q, dim] ] /;!(Head[b]===Plus);
     oneamp = qpcanc[Collect2[oneamp, q, Factoring->False] ];
(* order the denominators again *)
     If[ denorder === True,
         oneamp = denomExpand[ oneamp ]/.FeynAmpDenominator->denomOrder
       ];
      oneamp = (intcan[oneamp//momentumexpand]//ExpandScalarProduct)/.
                intcan->Identity ;
    ];
oneamp1 = oneamp;
print1["simplifying again"];
If[!FreeQ[oneamp, Spinor],
   oneamp = Map[spinorchainevaluate, oneamp +  nUUl]/.nUUl->0;
   oneamp = Expand2[ ExpandScalarProduct[ oneamp ], q]
  ];

(* ONEAMPCHANGE : cancelling q^2 *)
If[ IntegerQ[q2cancel],
For[iiq = 1, iiq<5, iiq++,
   print2["iiq = ",iiq];
   If[ IntegerQ[q2cancel] && 
       !FreeQ[  oneamp, Pair[ Momentum[q,___],Momentum[q,___] ]  ],
       oneamp = Collect2[oneamp, q, Factoring-> False];
       print1["cancelling q^2"];
       If[ Head[oneamp]===Plus,
           in3p = Select[ oneamp, FreeQ[#, Pair[ Momentum[q,___],
                                              Momentum[q,___]] ]& 
                        ];
           in3q =  Collect2[ oneamp - in3p, {Pair,q}, Factoring->False]; 
	   resin3q = q2canc[ in3q ]; 
           oneamp = resin3q + in3p,
           oneamp = q2canc[ oneamp ]
         ] ] 
  ]];

   print1["collecting w.r.t.", q, "  ", Length[oneamp], " terms"];

oneamp = FixedPoint[ ReleaseHold, smalld[oneamp] ];

timecoll = $FactorTime;
$FactorTime=360;
tim=Timing[

If[$ToughStuff =!=True, 
   oneamp = Collect2[ oneamp, q ],
   oneamp = Collect2[ oneamp, q, Factoring -> False];
  ];

oneamp2 = oneamp;
   If[ Head[oneamp] === Plus, 
       newoneamp = 0; 
       lenneu2 = Length[oneamp];
       If[ dim=!=4, qdi = {q,dim}, qdi = {q} ];
       For[ii=1, ii<=lenneu2, ii++,
           print2["isolating; ii  = ", ii, " out of ", lenneu2];
           newoneamp = newoneamp + Isolate[oneamp[[ii]], qdi,
                                           IsolateHead ->FC, 
                                           IsolateSplit->Infinity]
          ];
       oneamp = newoneamp
     ];
         ];

$FactorTime=  timecoll;
   print1["time for collect2  and isolate : ", tim[[1]]//FeynCalcForm];
   print1["isolated and collected, before tensorintegraldecomposition:/n 
           length of isolated graph = ",Length[oneamp]];
   print1["memory in use = ",MemoryInUse[]];
   print3["oneamp = ", oneamp];
   
(* ONEAMPCHANGE : tensorintegraldecomposition *)
   If[!FreeQ2[oneamp, {Pair[x_,y_Plus] , DiracGamma[x_Plus]}],
 print1["cheCK"];
      oneamp = DiracSimplify[oneamp]//expandscalarpoduct
     ];
sirlinsave = Feyn`Calc`Main`$sirlin;
Feyn`Calc`Main`$sirlin = False;
print1["q = ",q];
   oneamp = tensint[ oneamp,dim,q,{Mandelstam->onemandel} ];
   oneamp = oneamp /.{B1 :> bB1, B00 :> bB00, B11 :> bB11};
   oneamp = oneamp//smalld;
   oneamp = oneamp /.{bB1 :> B1, bB00 :> B00, bB11 :> B11};
Feyn`Calc`Main`$sirlin = sirlinsave;
   oneamp = FixedPoint[ReleaseHold, oneamp];
   print1["after tensint "];
   print3["after tensint : oneamp = ", oneamp];
   If[!FreeQ[oneamp, Spinor],
      oneamp = Collect2[oneamp, Spinor, Factoring -> False],
      If[!FreeQ[oneamp, Pair],
(*
         oneamp = Collect2[oneamp, {Pair, Polarization},   
                           Factoring -> False]
Change10/94
*)
         If[LeafCount[oneamp] > 100000,
            print2["specialexpanding"];
            holdoneamp = Hold @@ {oneamp};
            tone = 0; ltone = Length[oneamp];
            For[ijt = 1, ijt <= ltone, ijt++,
                print2["ijt = ",ijt, " out of ",ltone];
                tone = tone + Expand[holdoneamp[[1,ijt]]];
               ];
            Clear[holdoneamp];
            oneamp = tone;
           ];
            
,
         oneamp = Expand[oneamp]
        ]
     ];

(* StandardMatrixElement[ Spinor[] ... ] -> Spinor[] ... *)

standback[x_]:= x /; !FreeQ[ x, Spinor ];
oneamp = oneamp /. StandardMatrixElement -> standback /.
                   standback -> StandardMatrixElement;

oneampsave = oneamp;
print3["oneampsave = ", oneampsave];

(* ********************************************************************* *)
(*                          oneloop23                                    *)
(* ********************************************************************* *)

   print1["after collecting "];

(* ONEAMPCHaNGE : spinor stuff and matrixelements *)
   If[ (!FreeQ[ oneamp,Spinor ]),
           If[Length[oneamp]>42, 
              oneamp = Isolate[ oneamp, {Spinor, DiracGamma, Pair},
                                IsolateHead-> FCT, 
                                IsolateSplit -> Infinity];
             ];
           print1["length of oneamp now: ",Length[oneamp]];
           If[(reducegamma67 === True) ||
               (!FreeQ[oneamp, Literal[
                       Spinor[p1__] . a___ . Spinor[p2__] *
                       Spinor[p3__] . b___ . Spinor[p4__] ]
                      ]),
               oneamp = oneamp/. DiracGamma[7]->(1/2 - DiracGamma[5]/2)/.
                                 DiracGamma[6]->(1/2 + DiracGamma[5]/2);
             ];

       oneamp = Expand2[DiracSimplify[oneamp],Dot]//DiracOrder//
                 expandscalarproduct;
       oneamp = DiracSimplify[oneamp];

           If[(reducegamma67 =!= True)  && 
              (!FreeQ[oneamp, Literal[
                      Spinor[p1__] . a___ . Spinor[p2__] *
                      Spinor[p3__] . b___ . Spinor[p4__] ]
                     ]),
               oneamp = oneamp /.
                 {Literal[ Spinor[p1__].(a___).Spinor[p2__]] :>
                   (Spinor[p1].a.DiracGamma[6].Spinor[p2] +
                    Spinor[p1].a.DiracGamma[7].Spinor[p2]
                   ) 
                 };
               oneamp = DiracSimplify[oneamp]
             ];

       If[$higherpoint===False, 
          oneamp = FixedPoint[ReleaseHold, oneamp]
         ];
       print1["after spinorchainevaluate"]
     ];

   oneamp = oneamp /.{B1 :> bB1, B00 :> bB00, B11 :> bB11};
   oneamp = oneamp//smalld;
   oneamp = oneamp /.{bB1 :> B1, bB00 :> B00, bB11 :> B11};

(* If something changed, we have to order again *)
   If[oneamp =!= oneampsave,
      If[!FreeQ[oneamp, Spinor],
         oneamp = Collect2[oneamp, Spinor, Factoring -> False],
         If[!FreeQ[oneamp, Pair],
            oneamp = Collect2[oneamp, {Pair, Polarization},   
                              Factoring -> False],
            oneamp = Expand[oneamp]
           ]
        ]
     ];
   print1["after Collecting "];
   
standmatd[xxx__]:=StandardMatrixElement[dotdotlin[xxx]];
standmatd[]=1;
(* this may take a lot of time ... *)
  If[  (StandardMatrixElement =!= Identity),
  If[ !FreeQ2[ oneamp, {Dot, DiracGamma, Polarization,
                        SU3f, SU3delta, GellMannMatrix
                        } ] && FreeQ[ oneamp, Spinor ],
      If[Head[oneamp]===Plus,
         oneamp = Sum[ (oneamp[[iii]] spip[])/.
                       spip -> spinorsandpairs /.
                       spinorsandpairs -> standmatd, 
                       {iii,1,Length[oneamp]}],
         oneamp = Expand2[oneamp spip[], Polarization]/.
                  spip->spinorsandpairs /. spinorsandpairs -> 
                  standmatd
        ]
    ]];
  
  If[ (StandardMatrixElement =!= Identity) &&
     (!FreeQ2[oneamp, {Spinor, Polarization, SUNIndex}]),
     print1["before tostandmat"];
     If[ Head[oneamp] === Plus, 
         paone = Select[oneamp, !FreeQ2[#, {Spinor, Polarization, SUNIndex}]&];
         paone = {oneamp - paone, paone},
         paone = {0, oneamp}
       ];
     paone2 = Isolate[Collect2[paone[[2]], {Spinor,Polarization,SUNIndex},
                               Factoring -> False], {Spinor,Polarization,SUNIndex},
                               IsolateHead -> tempFC,IsolateSplit->Infinity];
     oneamp = FixedPoint[ReleaseHold, tostandmat[paone2] ] + paone[[1]];
     print1["after tostandmat"];
  print3["after tostandmat : oneamp = ",oneamp];
    ];
(* ********************************************************************* *)
(*                          oneloop25                                    *)
(* ********************************************************************* *)


(* ONEAMPCHANGE : spinor stuff and matrixelements *)
   If[ !FreeQ[oneamp, StandardMatrixElement], 
       print1["collecting w.r.t. standard matrixelements "];
       oneamp = collin[oneamp/.inisubs, StandardMatrixElement, True 
                      ];
       print1["collecting done"];
     ];

(* ONEAMPCHANGE : inserting the subdeterminants again *)
   If[!FreeQ[ oneamp, SUBDET ], 
      oneamp = oneamp /. SUBDET -> SUB;
      print1["subdeterminants reinserted"]
     ];

   tric[y_Plus]:=tric /@ y;
   tric[x_]:=TrickMandelstam[x, onemandel]/;Length[onemandel]>0;
   tric[x_]:=x/;onemandel==={};
   fma[xx_]:=True/;Head[xx]===StandardMatrixElement;
   vva[xx_]:=True/;!FreeQ[{xx},PaVe];

If[LeafCount[oneamp]<1000, 
   voneampsma =  Variables[oneamp]/.Negligible->Identity;
   smalist = Select[ voneampsma, fma ];
   pvlist = Select[ voneampsma, vva ]//Union;
   arglist = Union[pvlist/.PaVe->pvar],
   arglist = {}
  ];

   collpav[x_Symbol]:=x;
   collpav[a_StandardMatrixElement b_] := a collpav[b];
   collpav[x_]:=tric[ Collect2[x,{A0,B0,B1,B00,B11,C0,D0,PaVe},
                           Factoring -> True] ];
(* ********************************************************************* *)
(*                          oneloop26                                    *)
(* ********************************************************************* *)

(* substituting the final substitutions *)
fsub[x_, su_]:=Block[{nx=x,ij}, 
                         For[ij=1, ij<=Length[su], ij++,
                             nx=nx/.su[[ij]]
                            ];nx
               ];

   oneamp = fsub[ oneamp, finsubst ];
   If[!FreeQ[oneamp, SUNIndex],
      oneamp = SUNSimplify[oneamp, SUNFToTraces -> True];
     ];

(* getting a common factor *)

       oneamp = oneamp + nUl1 + nUl2;
   If[ fact===True,
       factor3[x_]:=Factor2[x, FactorFull -> False];
       npref[0]=0;
       npref[w_ v_]:=(factor3[w]/.Plus->pluS) npref[v]/;
         FreeQ2[w,{A0, B0,B1, C0, D0, B00, B11, PaVe}];
       oneamp = factor3[(npref /@ Map[factor3,oneamp ])/.nUl1->0/.nUl2->0
                          ]/.pluS->Plus/.npref->collpav ,
       oneamp = Map[collpav, oneamp] /. nUl1 ->0 /.nUl2 -> 0;
     ];

(* Isolating *)
   If[ isolatehead=!=False, isol[x_]:=Isolate[x,IsolateHead->isolatehead];
       sh[he_][x__]:=isol[he[x]];
       scaliso = {A0->sh[A0], B0->sh[B0], B1->sh[b1], B00->sh[B00], 
                  B11->sh[B11], C0->sh[C0], D0->sh[D0]};
       isoplu[x__]:=isol[Plus[x]];
       oneamp = isol[  oneamp/.D0->sh[D0]/.C0->sh[C0]/.
                          B11->sh[B11]/.B00->sh[B00]/.B1->sh[b1]/.
                          B0->sh[B0]/.A0->sh[A0]/.Plus->isoplu/.
          isoplu->Plus ];
       isol[x_]:=x
     ];
       
(* putting everything together again, including the prefactor *)
   oneampresult = oneamp;
   oneampresult = fsub[newprefactor, finsubst] oneampresult;

   print3["oneampresult = ",oneampresult];

   If[ isolatehead=!=False, oneampresult = isol[oneampresult] ];

](* end of If FreeQ oneamp, FeynAmpDenomiantor *) ;
 
(* ********************************************************************* *)
(*                          oneloop27                                    *)
(* ********************************************************************* *)

(* writing the result in a file specified by grname *)
   If[ (writeout=!=False )&&  StringQ[ name ] && 
       (!ValueQ[OneLoopResult[grname]]),
       If[ formattype===FortranForm,
            wri[ name, Hold[grname  = oneampresult], FormatType -> FortranForm
               ]/.wri -> Write2,
           wri[ name, Hold[OneLoopResult[ grname ] = oneampresult], 
           FormatType->formattype]/.wri -> Write2;
           If[Length[FA2Info] > 0,
              wri[ name, Hold[OneLoopInfo[ grname ] = FA2Info], 
              FormatType->InputForm]/.wri -> Write2;
             ];
         ]
     ];
   If[ (!ValueQ[ OneLoopResult[ grname ]]) && (grname=!=False), 
       set[ OneLoopResult[grname], oneampresult ]/.set->Set 
     ];
   If[ (!ValueQ[ OneLoopInfo[ grname ]]) && (grname=!=False), 
       set[ OneLoopInfo[grname], FA2Info]/.set->Set 
     ];

(* ********************************************************************* *)
(*                          oneloop28                                    *)
(* ********************************************************************* *)


(* now it's done and depending on $VeryVerbose some printing may be done *)
(* --------------------------------------------------------------------- *)
(* some cosmetics for print1 *)
(* --------------------------------------------------------------------- *)
(* only if no abbreviations are introduced: print eventually *)
(* for printing purposes abbreviations are useful,but  *)
(* this may actually under certain circumstances be incorrect! *)
(* Though the result returned by OneLoop is of course correct *)
If[ isolatehead===False,
    PaVeAbbreviate[x_]:=x/.PaVe->paVeabbrevi/.paVeabbrevi->PaVe;
    paVeabbrevi[x__,{y__},{m1_,m2_,m3_}]:=
    ToExpression[ StringJoin@@Prepend[Map[ToString,{x}],"C"] ];
    paVeabbrevi[x__,{y__},{m1_,m2_,m3_,m4_}]:=
    ToExpression[ StringJoin@@Prepend[Map[ToString,{x}],"D"] ];
    pva[xxx_]:= xxx//PaVeAbbreviate;
   pvar[xx__,li1_{},li2_List]:="As"[li2[[1]]]/;Length[li2]===1;
   pvar[xx__,li1_List,li2_List]:="Bs"@@Flatten[{li1,li2}]/;Length[li2]===2;
   pvar[xx__,li1_List,li2_List]:="Cs"@@Flatten[{li1,li2}]/;Length[li2]===3;
   pvar[xx__,li1_List,li2_List]:="Ds"@@Flatten[{li1,li2}]/;Length[li2]===4;

print2[" "];
print2[" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"];
print2[" "];
If[ grname=!=False, 
   If[ $VeryVerbose >0, Print[" The result for ",grname, " is ",
        Shallow[ oneampresult,{6, 20} ]
                                ],
   print3[oneampresult//FeynCalcForm] 
     ];
  ];

If[ Length[ arglist ]>0 && !breakdown===True, print2[ "Arguments:  ",
arglist//Union ]
  ];
print2[" "];
print2[" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"];
print2[" "];
If[MemoryInUse[]>10^7,
timsh = Timing[
Share[]
             ][[1]];
print1["time needed for Share = ",timsh];
         
print1[N[ MaxMemoryUsed[]/10^6,2 ], "MB used"];
  ];
  ] (* endIf IsolateHead *);

(* ----------------------------------------------------------------- *)
]; (* end Timing *)
If[grname=!=False,
   print1["CPU - time for ",grname, ": ", tim[[1]]//FeynCalcForm],
   print1["CPU - time : ", tim[[1]]//FeynCalcForm]
  ];


oneampresult]/;FreeQ[q,Rule];

(* ******************************************************************* *)



(* ********************************************************************* *)
(*                          oneloop30                                    *)
(* ********************************************************************* *)


(* OneLoopSumdef *)

Options[OneLoopSum]={ 
                      Born               -> 1,           
                      CombineGraphs      -> {},
                      Dimension          -> True,
                      ExtraVariables      -> {},
                      FinalFunction       ->Identity,
                      FinalSubstitutions -> {},
                      FormatType         -> InputForm,
                      InitialSubstitutions  -> {},
                      IntermediateSubstitutions -> {},
                      IsolateHead        -> K,
                      KeepOnly           -> False,
                      Mandelstam         -> {},
                      Prefactor          -> 1,
                      ReduceToScalars    -> True,
                      Scaling            -> {},
                      SelectGraphs       -> All,
                      WriteOutPaVe       -> False
                    };

SetAttributes[OneLoopSum, HoldAll];
OneLoopSum[ex_, ops___]:= Block[{mand,reduce,na,i,exx,nsres,sres,
            len0,len,jj, vars, varpave,isolatehead,alreadysummed,
            filename, formattype,selectgraphs,combinegraphs,pres,
            aa0,bb0,cc0,dd0,ddb0,finfunc,inisuB,extravars,
            acdc, lnw, nres, j, nvd, set, npi3, newpa, simp,nplin,
             mandelspec,sumcol,colll,prefactor,feynli,dims,fsub,keeponly,
            d0multiply, c0multiply, d0scalIsolate,c0scalIsolate,vsm,
           intermedsub,isol2,combinelist,np,nnp,npavopt,iiv,lres,mmsu,iim,
           feynAmpden1, feynAmpden2,masss1, masss2, fim,checklabel },

exx = ex;
sres = 0;

SetOptions[DiracTrace,  DiracTraceEvaluate -> False ];
exx = FixedPoint[ ReleaseHold, exx ];
opsli         = Join[ {ops}//ReleaseHold, Options[OneLoopSum] ];

born          = Born /. opsli;
combinegraphs = CombineGraphs/.opsli;
extravars     = ExtraVariables/.opsli;
isolatehead   = IsolateHead/.opsli;
If[isolatehead === True, isolatehead = IsolateHead/.Options[Isolate]];
dims           = Dimension /.opsli;
finalsubst    = FinalSubstitutions/.opsli;
finfunc       = FinalFunction/.opsli;
formattype    = FormatType/.opsli;
inisuB        = InitialSubstitutions /. opsli;
intermedsub   = IntermediateSubstitutions /. opsli;
keeponly      = KeepOnly /. opsli;
mandelspec    = Mandelstam/.opsli;
prefactor     = Prefactor/.opsli;
reduce        = ReduceToScalars/.opsli;
scaling       = Scaling /. opsli;
writeoutpave  = WriteOutPaVe/.opsli;

If[intermedsub =!= {},
   SetOptions[ OneLoop, IntermediateSubstitutions -> intermedsub ]
  ];

(* ********************************************************************* *)
(*                          oneloop31                                    *)
(* ********************************************************************* *)

If[inisuB =!= {}, exx = exx /. inisuB];
If[ !FreeQ[exx, FeynAmp],
    (* For FA2.0 *)
    If[ (Head[Head[exx]] =!= FeynAmpList) &&
        (Length[exx] === 1              ) && 
        (Head[Head[exx[[1]]]] === FeynAmpList),
        exx = exx[[1]]
      ];

(* bringing selectgraphs and combinegraphs in a standard form *)
selectgraphs  = SelectGraphs/.Join[ {ops}, Options[OneLoopSum] ];
print1["selectghraphs = ",selectgraphs];

If[ selectgraphs =!= All, 
    selectgraphs  = selectgraphs //. { a___, {i_, j_}, b___} :>
                                     { a, i, j, b} /; ( (j-i)^2 === 1);
    selectgraphs  = selectgraphs //. { a___, {i_Integer, j_Integer}, b___ } :>
                                     {a, Range[i,j],b};
    selectgraphs  = Flatten[ {selectgraphs} ] // Sort;
  ]; 

len0 = Length[exx];
If[ combinegraphs === All, combinegraphs = Range[len0] ];
If[ selectgraphs === All, selectgraphs = Range[len0] ];
  
combinegraphs = combinegraphs //. { a___, {b___, {i_,j_}, c___}, d___} :>
                                  { a, Join[{b}, Range[i,j], {c}], d};
If[ MatchQ[combinegraphs, {i___Integer}], combinegraphs={combinegraphs}];
combinegraphs = Map[ Sort, combinegraphs ];
If[combinegraphs =!= False, 
   ncombine = {};
   For[ ii = 1, ii <= Length[combinegraphs], ii++,
        If[ Length[ Intersection[selectgraphs, combinegraphs[[ii]]] ] > 0,
            AppendTo[ncombine, Intersection[selectgraphs,combinegraphs[[ii]]]]
      ]   ];
    combinegraphs = ncombine;
    combinelist   = Flatten[combinegraphs];
  ];

  print1["Selecting graphs # ",selectgraphs];
  print1["Combining graphs # ", combinegraphs];

If[ Head[exx]=!=FeynAmp,
    feynli[___][xx__]:={xx};
    feynli[xx__]:={xx};
    exx = {List@@(exx/.FeynAmpList->feynli)}//Flatten,
    exx = {exx}
  ];

If[ (Head[selectgraphs]===List && Length[selectgraphs]>0) ||
    (Head[combinegraphs]===List && Length[combinegraphs]>0),
    sumit[exli_, numli_List]:=
    Block[{nam,lq = First[exli][[2]]},
           nam = GraphName[exli[[1,1,1]],
                           ToExpression[ StringJoin@@ Map[ ToString,
                                         Table[First[ exli[[numli[[jj]]]]]//
                                               Last, {jj,Length[numli]} ]
                          ]            ]      ];
           amps = Sum[ exli[[numli[[ij]]]]//Last, {ij,Length[numli]} ];
           FeynAmp[nam,lq,amps]
         ];
    nex = {};
    alreadysummed = {};
    For[ ii = 1, ii <= len0, ii++,
      If[ MemberQ[selectgraphs, ii] && FreeQ[alreadysummed,ii], 
         If[ !FreeQ[combinelist, ii],
             sumli = sumit[exx, Select[combinegraphs, !FreeQ[#,ii]&][[1]]];
             AppendTo[nex, sumli];
             AppendTo[alreadysummed, Select[combinegraphs, !FreeQ[#,ii]&][[1]]],
            AppendTo[nex, exx[[ii]]];
            AppendTo[alreadysummed, ii]
           ]
        ]
       ];
    exx = nex
  ];
len = Length[exx];
fim[a_] := If[a==={}, False, a[[1]]];

For[ i=1, i<=len, i++,
     na = exx[[i,1]];
     If[ na=!=False, print1["Calculating ",na, " ; # ",i," out of ",len ] ];
     checklabel = False;
     If[(i > 1) && !FreeQ[exx[[i]], DiracTrace] && 
        (Last[exx[[i-1]]] =!= 0) && (Last[exx[[i]]] =!= 0),
        feynAmpden1 = Select[Last[exx[[i-1]]], !FreeQ[#, FeynAmpDenominator]&];
        feynAmpden2 = Select[Last[exx[[i]]], !FreeQ[#, FeynAmpDenominator]&];
        If[ (Head[feynAmpden2] === FeynAmpDenominator) && 
            (Head[feynAmpden1] === FeynAmpDenominator) && 
            Length[feynAmpden1] > 2, 
            masss1 = Union[ #[[2]]& /@ feynAmpden1];
            masss2 = Union[ #[[2]]& /@ feynAmpden2];
            If[ (Length[masss1] === 1) && (Length[masss2] === 1) && 
               (((!FreeQ[NegligibleVariables/.Options[OneLoop], masss1//fim]) &&
                 (!FreeQ[NegligibleVariables/.Options[OneLoop], masss2//fim])
                ) ||
                ((FreeQ[NegligibleVariables/.Options[OneLoop], masss1//fim]) &&
                 (FreeQ[NegligibleVariables/.Options[OneLoop], masss2//fim])
               )),
                masss1 = masss1[[1]];
                masss2 = masss2[[1]];
                raTio = (Last[exx[[i-1]]] /. masss1 -> masss2)/Last[exx[[i]]];
                If[NumberQ[raTio],
                   pres = OneLoop[exx[[i, 1]], exx[[i, 2]], 
                                  (OneLoopResult[exx[[i-1, 1]]]/ raTio ) /. masss1 -> masss2];
                   checklabel = True
                  ]
        ] ] 
       ];
     If[ (i > 2 ) && !FreeQ[exx[[i]], DiracTrace],
        feynAmpden1 = Select[Last[exx[[i-2]]], !FreeQ[#, FeynAmpDenominator]&];
        feynAmpden2 = Select[Last[exx[[i]]], !FreeQ[#, FeynAmpDenominator]&];
        If[ (Head[feynAmpden2] === FeynAmpDenominator) &&
            (Head[feynAmpden1] === FeynAmpDenominator) &&
            (Length[feynAmpden1] === Length[feynAmpden2]) &&
            (Length[feynAmpden1] > 2),
            masss1 = Select[Union[ #[[2]]& /@ feynAmpden1], FreeQ[#, 0]&];
            masss2 = Select[Union[ #[[2]]& /@ feynAmpden2], FreeQ[#, 0]&];
            If[ (Length[masss1] === Length[masss2]) && (Length[masss1]>0),
                (((!FreeQ[NegligibleVariables/.Options[OneLoop], masss1//fim]) &&
                 (!FreeQ[NegligibleVariables/.Options[OneLoop], masss2//fim])
                ) ||
                ((FreeQ[NegligibleVariables/.Options[OneLoop], masss1//fim]) &&
                 (FreeQ[NegligibleVariables/.Options[OneLoop], masss2//fim])
               )),
                masssub =  Table[ masss1[[iji]] -> masss2[[iji]], 
                                  {iji, Length[masss1]}];
                raTio = (Last[exx[[i-2]]] /. masssub)/Last[exx[[i]]];
                If[NumberQ[raTio],
                   pres = OneLoop[exx[[i, 1]], exx[[i, 2]],
                                  (OneLoopResult[exx[[i-2, 1]]]/ raTio ) /.
                          masssub];
                   checklabel = True
                  ]
        ] ] 
       ];  
  
     If[ checklabel === False,  
         pres = exx[[i]]/.FeynAmp -> OneLoop
       ];
     If[!FreeQ[pres, StandardMatrixElement],
        pres = Expand2[pres, StandardMatrixElement]
       ];
     sres =  sres + pres 
   ]
,
 (*If FreeQ FeynAmp*)

(* ********************************************************************* *)
(*                          oneloop32                                    *)
(* ********************************************************************* *)

sres = exx;

If[ keeponly =!= False,
    SetOptions[B0, B0Unique -> True];
    SetOptions[A0, A0ToB0   -> True];
    Which[ keeponly === D0,
           B0[__] := 0; C0[__] := 0; PaVe[__,{_,_,_},{_,_,_}] := 0,
           keeponly === C0,
           D0[__] := 0; B0[__] := 0,
           keeponly === B0,
           C0[__] := 0; D0[__] := 0,
           keeponly === {},
           C0[__] := 0;  D0[__] :=0
        ]
  ];


dsi[x__]:=dsi[x]=DiracOrder[ DiracSimplify[ Dot[x] ] ];
sres = sres /. Dot -> dsi;
(*
If[!FreeQ[sres, StandardMatrixElement], 
   sres = Expand2[ sres, StandardMatrixElement]
  ]
*)
] (*If FreeQ FeynAmp*);

If[combinegraphs === False, nres = sres,

If[ reduce === True,
    mand = Mandelstam/.Options[OneLoop],
    mand = {}
  ];
If[ (mand === {}) && (mandelspec=!={}), mand = mandelspec ];
print1["mand = ",mand];

   If[ Length[mand]===4,
       mansu={mand[[3]]->( mand[[4]] - mand[[1]] - mand[[2]] )},
       mansu = {}
     ];

   collp[x_]:= Block[{temp,ntemp,iit},
                     temp = x/.mansu;
  If[reduce =!= False,
  print2["collecting w.r.t. PaVe " ];
    If[!FreeQ[temp, PaVe],
       If[ Head[temp] =!= Plus,
           temp = Collect2[temp, {A0,B0,C0,D0,PaVe}, Factoring -> True],
           ntemp = 0;
           For[iit = 1, iit <= Length[temp], iit++,
               print2["collecting #  ",iit,  " out of ",Length[temp]];
               ntemp = ntemp + Collect2[temp[[iit]],  {A0,B0,C0,D0,PaVe}, 
                                        Factoring-> True];
              ];
           temp = Collect2[ntemp, {A0,B0,C0,D0,PaVe}, Factoring -> True]
         ];
    ];
print2["PaVe-collection done"];
    ];
                     temp
                  ];

If[FreeQ[sres, StandardMatrixElement],
   vsm = {},
   print1["searching StandardMatrixElement"];
   vsm = Variables[ sres /. {a_StandardMatrixElement _. :> a} ];
   vsm = Select[vsm, (Head[#] === StandardMatrixElement) &];
  ];

If[(!ValueQ[$SMECollect]) || ($SMECollect === True), 
If[vsm=!={}, print1["collect with respect to StandardMatrixElement"];
    nsres = 0;
   For[ij=1, ij<=Length[vsm], ij++,
 print1["ij = ",ij, "  out of ",Length[vsm]];
       dif =  D[ sres, vsm[[ij]] ];
       nsres = nsres + collp[ dif ] vsm[[ij]];
       sres = sres /. vsm[[ij]] -> 0
      ];
   sres = nsres + sres;
   print1["collecting done"],
sres = collp[sres]
  ];
  ];
npavopt = PaVeOrderList/.Options[PaVeOrder];
paveorder[xxx_]:=PaVeOrder[xxx, PaVeOrderList -> npavopt];
(* insert here eventually previously calculated PaVe's *)
sres = sres /. PaVe->pavesave /. pavesave -> PaVe;

(* ********************************************************************* *)
(*                          oneloop33                                    *)
(* ********************************************************************* *)

 pavvar[y_]:=Block[{alt,arr,ia,new},
                   alt = Drop[#,-1]& /@ Position[y,PaVe];
                   arr = {};
                   For[ia=1, ia<=Length[alt], ia++,
                       new = Part @@ Prepend[alt[[ia]], y];
                       If[!MemberQ[arr, new], AppendTo[arr,new] ]
                      ];
              arr];

If[ !FreeQ[sres, PaVe],
    varpave = {};
    If[Head[sres]===Plus,
       lres = Length[sres];
       For[iiv = 1, iiv <= lres, iiv++,
           print2["searching for PaVe;  iiv = ",iiv, " out of ", lres];
           varpave = Union[varpave, pavvar[sres[[iiv]]]];
          ],
       varpave = pavvar[sres]
      ];
    varpave = FixedPoint[ ReleaseHold, varpave ];
    lenpa = Length[varpave];

pavit[xXX_PaVe, dir_, prev_:False] := Block[{nx, file, temp, set,xxx},
   paV[xy__, p_List, m_List] := PaVe[xy,C,p,C,m];
   xxx = paV@@xXX;
   nx = StringReplace[ ToString[InputForm[xxx], PageWidth -> 222],
                       {", "->"","^"->"","{"->"", "/" -> "",
                       "}"->"", "["->"", "]"->"", "*" -> "", " " -> "" }
                     ];
                      nx = StringJoin[dir, nx, ".s"];
                      print1["nx = ",nx];
                      file = FileNames @@ {nx};
                      print1["file  =", file];
                      If[file =!= {},
                         temp =( Get @@ {nx} ) // paveorder;
 (* If something went wrong in writing the file *)                        
                         If[ Head[temp]=!=Plus, file = {} ]
                        ];
                      If[(file ==={}) && (keeponly === False),
tim = Timing[
                        If[prev === False,
                           temp = PaVeReduce[xXX, Dimension -> dims]//paveorder,
                           temp = paveorder[prev]
                          ];
            ][[1]];
print1["Time needed = ",tim//FeynCalcForm];
                         OpenWrite @@ {nx};
                         WriteString @@ {nx, "( "};
                         Write @@ {nx, temp};
                         WriteString @@ {nx, "  ) "};
                         Close @@ {nx}
                        ];
                           temp] (* pavitend *);


If[ reduce === True,
    For[ j=1,j<=lenpa,j++,
          print1["working with # ",j," out of ",lenpa ];
          print1["calculating ",InputForm[ varpave[[j]] ] ];

If[writeoutpave===True, writeoutpave = ""];
If[ !StringQ[writeoutpave],
tii=Timing[
          nvd = PaVeReduce[ varpave[[j]], IsolateHead ->False, 
                            Dimension -> dims
                          ] // paveorder
          ];
print1[tii[[1]]," needed"]
  ];

(* ********************************************************************* *)
(*                          oneloop34                                    *)
(* ********************************************************************* *)
SQR[xxx_]:=PowerExpand[Sqrt[xxx]];
              If[ StringQ[writeoutpave], 
                  (* Check if the difference w.r.t. the previous PaVe
                     is only in the mass arguments. *)
                  nvd = False;
                  If[ j > 1, 
                      If[ (Take[varpave[[j-1]],{-2,-2}]=== 
                           Take[varpave[[j]], {-2,-2}]) && 
                          ((FreeQ[{Last[varpave[[j-1]]], 
                                  Last[varpave[[j]]]}, Negligible]
                          ) || 
                          ( (Union[(!FreeQ[#,Negligible])& /@ 
                                   Last[varpave[[j-1]]]] === {True}) &&
                            (Union[(!FreeQ[#,Negligible])& /@ 
                                   Last[varpave[[j]]]] === {True})  
                          )) ,
                          mmsu = Table[SQR[Last[varpave[[j-1]]][[iim]]] ->
                                       SQR[Last[varpave[[j]]][[iim]]],
                                       {iim, Length[Last[varpave[[j]]]]}
                                      ];
                          
                          If[(varpave[[j-1]] /. mmsu) === varpave[[j]],
                             nvd = pavit[varpave[[j]], writeoutpave,
                                         (varpave[[j-1]]/.PaVe->pavesave
                                         ) /. mmsu
                                        ];
                  ] ] ];
                  If[nvd === False,
                     nvd = pavit[varpave[[j]], writeoutpave]
                    ]
                ];

              set[ varpave[[j]]/.PaVe->pavesave ,nvd ]/.set->Set
       ];
sres = sres /. PaVe->pavesave /. pavesave -> PaVe;
  ]

]; (* If !FreeQ[ sres, PaVe ] *)

   
   isol2[a_?NumberQ] := a;
   isol2[a_?NumberQ b_]:=a isol2[b];
   isol2[isol2[a_]]:=isol2[a];
   acdc = Join[extravars,{A0,B0,B1,B00,B11,DB0,C0,D0,PaVe}];
   acdc = Union[acdc, acdc /. finalsubst];
print1["acdc = ",acdc];
(* partdef *)
   part[a_Times]:=Block[{pAA}, 
                        pAA = Select[a, !FreeQ2[#, acdc]&];
                        pAA part[a/pAA]
                       ] /; !FreeQ2[a, acdc];
                  
   tog[y_] := Combine[ReleaseHold[y]/.mansu/.scaling/.scaling];
   If[ Length[mand]===4 ,
        sumcol[xx_]:=xx/.Plus->colll/.colll->Plus;
       colll[yy__]:=isol2[Collect2[Plus[yy], Variables[Take[mand/.
                                    scaling/.scaling,3]], 
                                   Factoring->True]
                         ] /; FreeQ2[{yy},acdc];
(*simpdef*)
       simp[y_]:=sumcol[Factor2[ y /. scaling /. scaling, FactorFull->False 
                               ]//smalld ],
       colll[yy__]:=isol2[Plus[yy]];
       sumcol[xx_]:=xx/.Plus->colll/.colll->Plus;
       simp[y_]:=Factor2[ y/. scaling /. scaling, FactorFull -> False]//smalld  
     ];
   If[ Length[mandelspec] === 4,
       simp[y_]:=sumcol[TrickMandelstam[
                             Factor2[ (y /. scaling /. scaling)//smalld,
                                      FactorFull -> False],
                                         mandelspec/.scaling]
                      ]
     ];
 born = simp[born /. scaling /. scaling];
   
 lnw = Length[sres];
 If[ Head[sres]=!=Plus, lnw = 1 ];
 nres = 0;
 print1["substituting "];
 sres =  sres/.mansu ;
 print1["done"];
(* here we have the loop over the StandardMatrixElement *)
  If[ FreeQ[sres, StandardMatrixElement],  lnw = 1];
  For[jj=1, jj<=lnw, jj++, 
      If[lnw === 1, 
         If[  FreeQ[sres, StandardMatrixElement],
              newpa = {sres, 1},
              newpa =  PartitHead[ Expand2[sres, StandardMatrixElement], 
                                   StandardMatrixElement ]
           ],
              newpa = PartitHead[ sres[[jj]],StandardMatrixElement ]
        ];
       print1[" # ",jj," out of ",lnw,"  ", newpa[[2]] ];
                (* Collect wrt. all the scalar integrals *)
If[MemoryInUse[] > 10 10^6, 
   print1["share "]; tis = Timing[Share[]][[1]];
   print1["time needed for share = ",tis];
  ];

print1["Shallow  ", Shallow[newpa[[1]]]];
np = newpa[[1]];
print1["leafcount of np = ",LeafCount[np]];
If[$ToughStuff === True || (LeafCount[np]>100000), 
  oldnp = np;
  np = Collect2[np, acdc, Factoring -> False]
  ,
tinp = Timing[
                np = Collect2[ np, acdc, Factoring -> True];
                nnp = paveorder[np];
                If[ np =!= nnp, 
                    np = Collect2[nnp, acdc, Factoring -> True];
                  ];
             ];
print1["timing for collecting = ",tinp[[1]]];
  ];
zero[__]:=0;
                (* combine the terms without PaVe's *)
                nplin = np/.A0->zero;
		nplin = nplin/.B0->zero;
		nplin = nplin/.C0->zero;
		nplin = nplin/.D0->zero;
		nplin = nplin/.DB0->zero;
		nplin = nplin/.B1->zero;
		nplin = nplin/.B00->zero;
		nplin = nplin/.B11->zero;
                nplin = nplin/.PaVe->zero;
                If[extravars =!= {},
                   For[iext = 1, iext <= Length[extravars], iext++,
                       nplin = nplin /. extravars[[iext]] -> 0
                      ]
                  ];
If[ nplin === 0 , print1["nplin = 0"],
    print1["leafcount of nplin = ",LeafCount[nplin]];
  ];
                If[nplin =!= 0, 
                   If[ keeponly === False,
                       np = np - nplin ,
                       If[ keeponly === {}, np = 0,
                           If[keeponly === B0 || keeponly === C0 || 
                              keeponly === D0, 
                              np = np - nplin;
                              nplin = 0 ]
                     ]   ]
                  ];
                pres = 0;
                If[ Head[np]===Plus,
                    lnp = Length[np];
       print1["combining coefficients of A0, B0, ..."];
                 (* This ist the loop of A0B0C0D0 *)
                    (* putting it now over a common denominator *)
                    For[i3=1, i3<=lnp, i3++,
                        print1["i3 = ",i3,"   out of ",lnp, "  
LeafCount = ",LeafCount[np[[i3]]] ];
                        npi3 = finfunc @@ {part[ finfunc[np[[i3]]] 
                                               ]/.part->simp};
                        npi3 = part[npi3 born]/.part->simp;
print2["npi3 = ",npi3];
                        pres = pres + npi3 
                       ], 
                    pres = part[np ]/.part->simp;
                    If[born=!=1,
                       pres = part[ pres born ] /. part -> simp 
                      ]
                  ];
print2["factoring nplin, LeafCount =  ", LeafCount[nplin]];
nplin = Cancel[simp[ Factor2[finfunc[nplin]]  ] * 
               simp[ Factor2[finfunc[born]]]//finfunc ];
print2["nplin = ",nplin];
                 nres = nres +  newpa[[2]] (pres + nplin)
      ] (* endFor*);

(* ********************************************************************* *)
(*                          oneloop35                                    *)
(* ********************************************************************* *)


fsub[x_]:=Block[{nx=x,su,ij}, su = finalsubst;
                         For[ij=1, ij<=Length[su], ij++,
                             nx=nx/.su[[ij]]
                            ];nx
               ];
mand=fsub[mand];
nres = fsub[ FixedPoint[tog, prefactor, 5] nres];
check = nres;

    {aa0, bb0, bb1, bb00, bb11, ddb0, cc0, dd0} = 
    {A0, B0,   B1,  B00,  B11,  DB0,  C0,  D0} // fsub;

If[ isolatehead=!=False,

    print1["isolating now "];

    plupp0[x__]:=Plus[x] /; !FreeQ[{x},plupp0];
    plupp1[x__]:=Factor2[ains Plus[x]];
    isolfact[x_]:= isol2[x/.Plus->plupp0/.plupp0->plupp1]/.ains->1;
    If[ Length[mand]===4,
        isolmand[x_]:= Isolate[x, {mand[[1]],mand[[2]],mand[[3]]},
                                       IsolateHead->isolatehead],
        isolmand[x_]:=Isolate[x,IsolateHead->isolatehead ]
      ];

    isolate0[x_]:=Isolate[x, IsolateHead->isolatehead ];

    isc[x_][y__]:=isol1[x][(TrickMandelstam[fsub[x[y]]/.dd0->D0,mand
                                           ]//paveorder)/.
                            D0 -> dd0, IsolateHead->isolatehead];
If[Length[mand]===4,
    isol1[_][x_, y_]:= Isolate[x, y] /; FreeQ2[x, Take[mand, 3]]
  ];
(* for scaling *)
   d0multiply = (D0 /. scaling) /. D0 -> 1;
   c0multiply = (C0 /. scaling) /. C0 -> 1;
   db0multiply = (DB0 /. scaling) /. DB0 -> 1;
   d0scalIsolate[x_,he_] := Isolate[x d0multiply, he];
   c0scalIsolate[x_,he_] := Isolate[x c0multiply, he];
   db0scalIsolate[x_,he_] := Isolate[x db0multiply, he];

    nres = isolate0[( nres )/.
                           dd0    -> isc[dd0]/.  
                           cc0    -> isc[cc0]/.
                           bb11   -> isc[bb11]/. bb00  -> isc[bb00]/.
                           bb1    -> isc[bb1]/.  bb0   -> isc[bb0]/. 
                           ddb0   -> isc[ddb0]/.
                           aa0    -> isc[aa0]/.  PaVe -> isc[PaVe]/.
                           isol1[dd0] -> d0scalIsolate/.
                           isol1[cc0] -> c0scalIsolate/.
                           isol1[ddb0] -> db0scalIsolate/.
                           isol1[bb1] -> Isolate/.
                           isol1[bb00] -> Isolate/.
                           isol1[bb11] -> Isolate/.
                           isol1[bb0] -> Isolate/.
                           isol1[aa0] -> Isolate/.
                           isol1[PaVe] -> Isolate/.
                           isol2 -> isolfact/.
                           isol2 -> isolmand
                   ],
(* ********************************************************************* *)
(*                          oneloop36                                    *)
(* ********************************************************************* *)

(* If isolatehead .. *)
(*Only if the option Factoring of OneLoop is True, factor also here *)
    If[(Factoring/.Options[OneLoop]) === True,
       specrule = {(a_Symbol - b_Symbol) (a_Symbol+b_Symbol)->(a^2-b^2)};
       factor3[x_]:=Factor2[x]/.specrule;
       isol22[a_ b_]:=isol22[a] isol22[b];
       isol22[a_]:=a/;Head[a]=!=Plus;
       isc2[x_][y__]:=(TrickMandelstam[x[y]/.dd0->D0,mand]//paveorder)/.
                     D0 -> dd0;
       nres = nres/.dd0 -> isc2[dd0]/.cc0    -> isc2[cc0]/.
                    bb11-> isc2[bb11]/. bb00 -> isc2[bb00]/.
                    bb1 -> isc2[bb1]/. bb0   -> isc2[bb0]/.
                    dbb0 -> isc2[ddb0] /.
                    aa0 -> isc2[aa0]/.  PaVe -> isc2[PaVe];
       nres = nres/.isol2 -> isol22;
       nres = Map[factor3, nres + nuLL]/.specrule;
       colp[x__]:=Map[TrickMandelstam[#,mand]&, 
                      Collect2[Plus[x], {aa0,bb0,bb00,bb11,bb1,ddb0,
                                         cc0,dd0,PaVe}, 
                               Factoring -> True]
                     ];
       nres = factor3[ Map[(#/.Plus->hoLdP)&, nres]/.nuLL->0/.
                       hoLdP[0]->0 ] /. hoLdP -> colp,
       nres = nres/.isol2->Identity
      ];
       nres = nres/.isol22->Identity ;
   disc[y__]:=TrickMandelstam[ D0[y],mand ]//paveorder;
   cisc[y__]:=TrickMandelstam[ C0[y],mand ]//paveorder;
   dbisc[y__]:=TrickMandelstam[ DB0[y],mand ]//paveorder;
   nres = nres /. D0->disc /. C0->cisc /. DB0 -> dbisc;
  ];
 ](*If combinegraphs ... *);
print2["The result of OneLoopSum is ", nres];
nres];
(*endOneLoopSum *)


(* ******************************************************************* *)

     (*smallddef *)
small2/: small2[x_]^n_ := small2[x^2] /; n > 0;
small2/: small2[_] a_ :=0;
small3/: small3[_] + a_ :=a;
small4[x_^m_]:=Negligible[x]^m;
   smalld[x_]:=x/;FreeQ[x,Negligible];
   smalld[x_]:=x/.Negligible->small2/.small2->small3/.
                         small3->small4/.small4->Negligible;

(* ******************************************************************* *)
(* ********************************************************************* *)
(*                          oneloop37                                    *)
(* ********************************************************************* *)

                                               (*spinorsandpairsdef*)
   dotdotlin[x___]:=dotlin[Dot[x]];
(*tempstandmatdef*)
   standma[x_]:=dotdotlin[x]/;!FreeQ2[x, {Polarization, Dot}];
   tempstandmat[x_]:=Block[{ttt},
                     ttt = x;
If[StandardMatrixElement =!= Identity,
                     If[(*FreeQ[ttt, Spinor] && *)
                        !FreeQ2[ttt,{SUNF,SUNDelta, SUNT, Pair}],
If[LeafCount[ttt]>500, print2["expanding in tempstandmat"]];
                        ttt = Expand[ttt spinorsandpairs[]];
If[LeafCount[ttt]>500, print2["expanding in tempstandmat done"]];
                        ttt = ttt /. spinorsandpairs -> dotsp
                       ];
                       
 
                     If[(Length[DownValues[spinorsandpairs]]>1) ||
                        ValueQ[StandardMatrixElement],
                     ttt = x/.Dot->spinorsandpairs/.
                           spinorsandpairs->StandardMatrixElement/.
                           StandardMatrixElement->standma/.standma->
                           StandardMatrixElement;
                       ];
(*
                    If[$fourfermion =!= True,
                       If[FreeQ[ttt, Polarization] && !FreeQ[ttt,Spinor],
                          ttt = SpecificPolarization[ttt]
                         ]
                      ];
*)
  ];
                       ttt];

 

   spinorsandpairs[a_,b__] := dotdotlin[a,b]//spinorsandpairs;
   dotsp[]=1;
   dotsp[x_]:=x;

   spinorsandpairs/: spinorsandpairs[x___] Pair[ Momentum[a__],
                                               Momentum[b__]
                                             ]^n_. :=
    spinorsandpairs[dotsp[x] Pair[Momentum[a],Momentum[b]]^n]/;
     !FreeQ[{a,b},Polarization];
   spinorsandpairs/: spinorsandpairs[x___] Eps[w__] :=
   spinorsandpairs[dotsp[x] Eps[w]]/;!FreeQ[{w}, Polarization];

   spinorsandpairs/:
 spinorsandpairs[x___] a_SUNT:= spinorsandpairs[dotsp[x], a];

   spinorsandpairs/:
 spinorsandpairs[x___] a_SUNF:= spinorsandpairs[dotsp[x] a];

   spinorsandpairs/:
 spinorsandpairs[x___] a_SUNDelta:= spinorsandpairs[dotsp[x] a];

 spinorsandpairs/: spinorsandpairs[x___] spinorsandpairs[y___] :=
                      spinorsandpairs[dotsp[x] dotsp[y]];

(* ********************************************************************* *)
(*                          oneloop38                                    *)
(* ********************************************************************* *)

(* *************************************************************** *)
(* Tensorintegraldecomposition *)
(* *************************************************************** *)
  (* PropagatorDenominatordef *)

(* This is needed, since the famous FeynArts creates tadpoles with 
   photon propagators ...
*)
Literal[PropagatorDenominator[0,0]]:=Block[{},
Print["!@#&**T(&()*&G_*&+#@U(BJ)IM:LM:LM<WPYOU_()FJ_jdfmbw34059hj"];
Print["nononononononononononononononononononononononononononononon"];
Print[" "];Quit[]
                                 ];
(* Tadpoles ( here was a bug found by G.Weiglein and M.Boehm (yeah!!!!) ...) *)
 PropagatorDenominator[0, m_] := -1/m^2;
 PropagatorDenominator[x_]:=PropagatorDenominator[x,0]/;FreeQ[x,Pattern];
 PropagatorDenominator[x_, y_] := (PropagatorDenominator[x, y] = 
         PropagatorDenominator[
  momentumexpand[ Momentum[x] ],y])/;FreeQ2[x,{Momentum, Pattern}];
  momentumexpandsave[a_] := momentumexpandsave[a] = momentumexpand[a];
  PropagatorDenominator[x_, y_] := (PropagatorDenominator[x, y] = 
  PropagatorDenominator[x//momentumexpand, y] )/; 
   (FreeQ[{x,y}, Pattern] ) && ( momentumexpandsave[x] =!= x);

 denomExpand[y__] := y/.FeynAmpDenominator->denexp; (*denomExpanddef*)
 denexp[z__] := Expand //@ momentumexpand[ FeynAmpDenominator[z] ];
(* *************************************************************** *)
(* suind substitutes for qu dummy indices for the tensor integral  *)
(* decomposition, the first argument of suind is a sumand *)
(* *************************************************************** *)
  (* Hier das RICHTIGE  suind ::: *)
   suind[ y_,qu_,dim_,md_] := Block[ {i,res=y, posli,
                                      currentposli},          (*suinddef*)
           posli = Position[y,Momentum[qu,___] ];
           If[posli=!={},
              For[i=1, i <= Length[posli],i++,
                  currentposli = Position[res, Momentum[qu,___] ];
                  res = ReplacePart[res, LorentzIndex[md[i],dim], 
                                    currentposli[[1]] ]
                 ]
             ];
                           res];  


(* *************************************************************** *)
(* for the divergent parts  "epsilon - substitution"               *)
(* *************************************************************** *)
   epst[x__] := If[uvpart === True, 0, epst2[x] ];
   epst2[gr_,x_,4,resid_] := to4dim[ gr ];       (*epstdef*)
   epst2[gr_,x_,d_Symbol,_] := to4dim[ gr ]/;FreeQ[x,d];
   epst2[gr_,x_,d_Symbol,resid_] := Block[{epstresul,epstin,epseps,epx},
       epx = to4dim[x];        
    epstin = Expand[(epx/.d->(4-epseps))-(epx/.d->4)
                  ]//spinorchainevaluate;
       epstresul = to4dim[gr + Normal[ Series[epstin,{epseps,0,1}]
                                     ]/.epseps->resid ]//Expand;
                               epstresul];

(* ********************************************************************* *)
(*                          oneloop39                                    *)
(* ********************************************************************* *)

(* *************************************************************** *)
(* A useful evaluation function ( for  tensint )                   *)
(* *************************************************************** *)
   to4d2[x_]:= x /;  $LimitTo4 =!= True;
   to4d2[x_]:=(x/.{Momentum[fope_,_]:>Momentum[fope]/;FreeQ[fope,q],
                   LorentzIndex[muu_,_] :> LorentzIndex[muu]}
              ) /; $LimitTo4 === True;
   dirsim[a_ b_]:=a dirsim[b] /; FreeQ2[a, {Spinor,DiracGamma}];
 SetAttributes[eval,Listable];                            (*evaldef*)
   eval[evy_] := Feyn`Calc`Main`memset[ eval[evy],   
                 Block[{evalte,nul1,nul2,ie,neval,nt},
       evalte = to4d2[ evy/.NonCommutativeMultiply->Times ];
       If[  !FreeQ[ evalte, LorentzIndex ], 
            evalte = contractli[ evalte ];   
         ];
       evalte = FixedPoint[ReleaseHold,evalte]//to4d2;
       evalte = Expand2[ dotlin[ evalte ]//expandscalarproduct, Dot ];
       If[(Length[evalte]>5) && !FreeQ[evalte, Dot] , 
          evalte = Collect2[ evalte, Dot, Factoring -> False]
         ];
       If[  !FreeQ[ evalte, DiracGamma], 
            evalte = Map[ dirsim, evalte + nul1 ]/.
                    dirsim->DiracSimplify/.nul1 -> 0;
            If[LeafCount[evalte]>100, 
               evalte = Collect2[evalte,Dot, Factoring -> False];
              ];
            evalte = evalte + nul1 + nul2;
            neval = 0;
            For[ie=1, ie<=Length[evalte], ie++,
                If[Length[evalte[[ie]]]>0, 
                   print3["ie = ",ie, " out of ",Length[evalte]]
                  ];
                nt = Contract[DiracOrder[evalte[[ie]]]
                             ]//expandscalarproduct;
                nt = DiracSimplify[nt]//DiracOrder;
                nt = Expand2[nt, Dot];
                print3["length of nt = ",nt//Length];
                neval = neval + nt
               ];
            evalte = neval/.nul1->0/.nul2->0;
          ];
       If[  !FreeQ[evalte, Eps], 
            evalte = evalte//EpsEvaluate//epschisholm;
               (* The default is that Eps's will be contracted away!*)
            evalte = Expand[ conall[ evalte ]//expandscalarproduct];
            evalte = epschisholm[ evalte ]//DiracSimplify//DiracOrder;
            evalte = Contract[ evalte ];
            evalte = Expand[ evalte//expandscalarproduct ]
         ];
       evalte = tempstandmat[ evalte ];
       evalte = Expand[evalte];
                 evalte]];

(* ********************************************************************* *)
(* ********************************************************************* *)
(*                          oneloop40                                    *)
(* ********************************************************************* *)
   (*tensintdef*)
tensint[x_,dim_,q_,options___] := (*tensint[x,dim,q,options]=*)
  Block[{tensj,tensi,tensic,tensg=0,mandel,tensx=x,tensdnp,tensdnp1,
         tenslnt,tensldn,tensqmax,tenslep,tensdnqq,tensdnqqb,
         tensqc,tensjq,tensfq,ltx
        },
print2["entering tensint "];
print3["entering tensint",q,"  dimension  ",dim,"   ",x//FeynCalcForm];
(* diracSimplify must have been used previously)              *)
 
   mandel =  Mandelstam/.Join[ options,Options[ tensint ] ];
(* tensor integral decomposition *)
   ltx = nterms[tensx];

   If[  Head[tensx]===Plus, tenslnt = Length[tensx], tenslnt = 1  ];
 
(* The tensj - loop runs over all different q-monomials *)
    For[ tensj=1, tensj <= tenslnt, tensj++,
         print1["tensorintegral # ",tensj," / ",tenslnt ];
         tensqc[tensj][any_]:=0;
         If[ tenslnt===1,
             tensdnp = PartitHead[ tensx,FeynAmpDenominator ],
             tensdnp = PartitHead[ tensx[[tensj]],FeynAmpDenominator ]
           ];
(* Collect according to the number of q's *)
   tensdnp1 = Collect2[ tensdnp[[1]],q, Factoring -> False];
   pairpow/: pairpow[a___,Momentum[q,di___],b___]^n_Integer?Positive :=
            (pairpow[a,Momentum[q,di],b]^(n-1))**pairpow[a,Momentum[q,di],b];
oten1=tensdnp1;
   tensdnp1= tensdnp1/.Pair->pairpow/.pairpow->Pair;
 
   If[ Head[tensdnp1]===Plus, tensldn = Length[tensdnp1], tensldn = 1];
   tensqmax[tensj]=0;
   For[ tensic=1, tensic <= tensldn, tensic++,
        If[   tensldn===1, 
            tenslep = Length[Position[tensdnp1,q ] ];
            tensqc[tensj][tenslep] += suind[ tensdnp1,q,dim,mud ],
             tenslep = Length[ Position[tensdnp1[[tensic]],q ] ];
             tensqc[tensj][tenslep] += suind[tensdnp1[[tensic]],q,dim,mud]
          ];
        If[  tenslep > tensqmax[tensj], tensqmax[tensj] = tenslep  ]
      ](*tensic - loop*);
                         
   For[ tensjq=0, tensjq <= tensqmax[tensj], tensjq++,
        tdenlen = Length[ tensdnp[[2]] ];
        print2["Tensorintegral (N = ",tdenlen,
               ") : # of q's = ",tensjq,
               " decomposing ", Length[tensqc[tensj][tensjq]]," term(s)"
              ];
        tensg += tdec[ tensqc[tensj][tensjq], tensdnp[[2]],q,
                       tensjq,dim,mud,mandel
                     ]/.NonCommutativeMultiply->Times
      ](*tensjq - loop*)
       ](*tensj - loop*);
  Expand[tensg]] (* end tensint *);
 
(* ********************************************************************* *)
(*                          oneloop41                                    *)
(* ********************************************************************* *)

(* *************************************************************** *)
(* tensor integrals; "qn" denotes the number of "q's"  *)
(* *************************************************************** *)
pavremember[x__] := Feyn`Calc`Main`memset[pavremember[x], PaVeReduce[x]];

   tdec[ expr_,props_,Q_,qn_ ,di_,mudu_,mand_]:=          (*tdecdef*)
   Feyn`Calc`Main`memset[ tdec[ expr,props,Q,qn,di,mudu,mand],
   Block[{spl0, mande,tensps={},tdecnew,tdec0j,tdectij, 
           tensdf2,tensdf1, pav0,
           tdecex = expr/.Pair-> Feyn`Calc`Main`sCO,   
           tdi,tdecti,tdectj,tdectk,tdectl,tdectm,tdecr=0,tdecpl,tdecml,tdeclpl,
           rul,spl,add 
          },
   print3["entering tdec with ", expr//FeynCalcForm];
   tensdf2[_,b_]:=b;
   tensdf1[a_,_]:= Expand[ to4dim[
                           momentumexpand[a-Momentum[Q]]]
                         ];
   If[ tdecex===0, tdecr = 0,
   add[gra_,pva_,exp_]:= Block[{addre, pv = pva},
      If[breakdown === True, 
opv = pv;
         pv = pavremember[pv, WriteOutPaVe -> writeoutpav ];
         If[uvpart === True,
            pv = pv /. C0[__] -> 0 /. D0[__]->0 /. B0[__]-> UVDELTA;
            If[FreeQ2[pv, {B1,B00,B11}], 
print2["pv = ",pv//FeynCalcForm];
              pv = Factor2[D[ pv, UVDELTA]],
print2["pv = ",pv//FeynCalcForm];
              Print["problems with uvcheck in OneLoop!!!", Dialog[]]
              ];
print2["uvcheck ",pv//FeynCalcForm];
           ];
        ];
(*XXX*)
      If[ $LimitTo4 === True,
          addre = gra+( Expand[
                         expandscalarproduct[ pv to4dim[ exp/.di->4 ] ]
                              ]),
          addre = gra+( Expand[
                         expandscalarproduct[ pv exp]
                     ])
        ];              addre];

(* calculate the List of scalar products needed as arguments *)

(* get the list of p's from the propagators *)
   tdecpl = Drop[ Expand[ props//momentumexpand]/.
                  PropagatorDenominator->tensdf1/.
                  FeynAmpDenominator->List,1 ]//diracgammacombine;
   tdecpl = Expand[tdecpl];
   print2["tdecpl = ",tdecpl];
   tdecml = props/.PropagatorDenominator->tensdf2/.
                   FeynAmpDenominator->List;           
   print3["tdecml = ",tdecml];
   tdecml = #^2& /@ tdecml;
   print3["tdecml = ",tdecml];
   tdeclpl = Length[ tdecpl ];
(* D_0, D_mu, D_munu and D_munuro are 4-dimensional, if $LimitTo4 is True *)
  If[ ($LimitTo4===True) && (tdeclpl===3) && (qn<4), 
      tdecex = tdecex/.di->4; tdi=4,
      tdi=di     
   ];
(* calculation of (N (N-1)/2) scalar pipj - arguments *)
   spl0[a_,b_,man_]:=(TrickMandelstam@@
                      Prepend[{man},Expand[Pair[a-b,a-b] ]//expandscalarproduct]
                     )//smalld;
   spl[aa__]:=spl0[aa,mand]//ExpandAll;
   Which[ tdeclpl == 1, tensps = { spl[ tdecpl[[1]],0 ] },
          tdeclpl == 2, tensps = { spl[ tdecpl[[1]],0 ],
                                   spl[ tdecpl[[2]],tdecpl[[1]] ],
                                   spl[ tdecpl[[2]],0] },
          tdeclpl === 3, tensps = { spl[ tdecpl[[1]],0 ],
                                    spl[ tdecpl[[2]],tdecpl[[1]] ],
                                    spl[ tdecpl[[3]],tdecpl[[2]] ],
                                    spl[ tdecpl[[3]],0 ],
                                    spl[ tdecpl[[2]],0 ],
                                    spl[ tdecpl[[3]],tdecpl[[1]] ]}
(*,
          tdeclpl === 4, tensps = { pl[ tdecpl[[1]],0 ],
                                    spl[ tdecpl[[2]],tdecpl[[1]] ],
                                    spl[ tdecpl[[3]],tdecpl[[2]] ],
                                    spl[ tdecpl[[3]],0 ],
 ...                                  
                                   
}
*)
        ];
(* scalar integrals *)
   If[ qn==0,
       tdecnew = eval[ tdecex ];
      If[ $LimitTo4 === True,
       Which[ tdeclpl == 0,                        (* e A0 = 2m^2 *)
              tdecr = epst[ tdecr,tdecnew,tdi, 2 tdecml[[1]] ],
              tdeclpl == 1,
              tdecr = epst[ tdecr,tdecnew,tdi, 2 ]; (* e B0 = 2 *)
            ] 
        ];
(* if the option DenominatorOrder is True, then order here again *)
       If[ denomOrder === True, 
           pav0 = PaVeOrder[PaVe[0,tensps,tdecml]],
           pav0 = PaVe[0,tensps,tdecml]
         ];
       tdecr = add[ tdecr, pav0, tdecnew ]
     ];

   If[ qn==1,
       tdecnew = Table[  eval[ tdecex/.LorentzIndex[mudu[1],___]->
                               tdecpl[[tdecti]] ], {tdecti,1,tdeclpl}
                      ];
       If[  ($LimitTo4 === True) && (tdeclpl === 1),   (* e B1 = -1 *)
            tdecr = epst[ tdecr,tdecnew[[1]], tdi,-1 ]
         ];
            For[ tdectj=1,tdectj<=tdeclpl,tdectj++,
                 tdecr = add[ tdecr, PaVe[tdectj,tensps,tdecml],
                              tdecnew[[tdectj]]
               ]            ]
     ];

   If[ qn==2,
       tdecnew = eval[ tdecex/.LorentzIndex[mudu[1],dime___]->
                               LorentzIndex[mudu[2],dime]
                     ];
      If[$LimitTo4 === True,
       Which[ 
             tdeclpl == 0,        (* e A00 = m^4/2 *)
                            tdecr = epst[  tdecr,tdecnew,tdi,
                                           tdecml[[1]]^2/2
                                        ],
             tdeclpl == 1,        (* e B00  *)
                            tdecr = epst[  tdecr,tdecnew,tdi,
                                           (-1/3 spl[tdecpl[[1]],0] +
                                             tdecml[[1]] +
                                             tdecml[[2]] )/2
                                        ] ,
             tdeclpl == 2,       (* e C00 = 1/2 *)
                            tdecr = epst[ tdecr, tdecnew, tdi, 1/2 ]
            ] 
        ];
       tdecr = add[ tdecr, PaVe[0,0,tensps,tdecml], tdecnew ];

       tdecnew = Table[{Sort[{tdecti,tdectj}],
                 tdecex/.LorentzIndex[mudu[1],___]->tdecpl[[tdecti]]/.
                         LorentzIndex[mudu[2],___]->tdecpl[[tdectj]]
                       },{tdectj,1,tdeclpl},{tdecti,1,tdeclpl}
                      ];
       tdecnew = eval[ Flatten[ tdecnew,1 ] ];
       If[ ($LimitTo4 === True) && (tdeclpl == 1),  (* e B11 = 2/3 *)
           tdecr = epst[ tdecr,tdecnew[[1,2]],tdi, 2/3 ]
         ];
       For[ tdectj=1,tdectj<=Length[tdecnew],tdectj++,
            tdecr = add[ tdecr,
                    PaVe@@Join[tdecnew[[tdectj,1]],{tensps},{tdecml}],
                         tdecnew[[tdectj,2]]
                       ];            
            If[$LimitTo4 === True, tdecr = tdecr /.tdi->4];
          ]
     ];

   If[ qn == 3,               (* The  00i - terms *)
       tdecnew ={};
        For[ tdectij = 1, tdectij <= tdeclpl, tdectij++,
            tdecnew = Append[ tdecnew, eval[
        conall[ tdecex/.NonCommutativeMultiply->Times/.
                {LorentzIndex[mudu[1],dime___]->
                 LorentzIndex[mudu[2],dime],
                 LorentzIndex[mudu[3],___]->tdecpl[[tdectij]]}]
      + conall[ tdecex/.NonCommutativeMultiply->Times/.
                {LorentzIndex[mudu[2],dime___]->
                 LorentzIndex[mudu[3],dime],
                 LorentzIndex[mudu[1],___]->tdecpl[[tdectij]]}]
      + conall[ tdecex/.NonCommutativeMultiply->Times/.
                {LorentzIndex[mudu[1],dime___]->
                 LorentzIndex[mudu[3],dime],
                 LorentzIndex[mudu[2],___]->tdecpl[[tdectij]]}]
                                          ]
                            ]
           ];
       If[ ($LimitTo4 === True) && (tdeclpl==2),  (* C001 = -1/6 *)
           tdecr = epst[ tdecr,tdecnew[[1]], tdi,-1/6 ];
                tdecr = epst[ tdecr,tdecnew[[2]], tdi,-1/6 ]
         ];
       For[ tdec0j = 1, tdec0j <= tdeclpl, tdec0j++,
            tdecr = add[  tdecr, PaVe[0,0,tdec0j,tensps,tdecml],
                          tdecnew[[tdec0j]]  ]
          ];

 
       For[ tdecti = 1, tdecti <= tdeclpl, tdecti++,
            For[ tdectj = 1, tdectj <= tdeclpl, tdectj++,
                 For[ tdectk=1, tdectk <= tdeclpl, tdectk++,
                      tdecr = add[ tdecr,
                              PaVe@@Join[Sort[{tdecti,tdectj,tdectk}],
                                              {tensps},{tdecml}
                                             ],
                           (eval[ tdecex/.LorentzIndex[mudu[1],___]->
                                          tdecpl[[tdecti]]/.
                                          LorentzIndex[mudu[2],___]->
                                          tdecpl[[tdectj]]/.
                                          LorentzIndex[mudu[3],___]->
                                          tdecpl[[tdectk]]
                                ]
                           )
                                ];                
                      If[ $LimitTo4 === True, tdecr = tdecr/.tdi->4 ]
                    ]        ]                   ] 
   ];

(* D_munurosi *)
 If[ qn == 4,    (* the 0000 - terms *)
     tdecnew = eval[(tdecex/.LorentzIndex[mudu[1],dime___] -> 
                             LorentzIndex[mudu[2],dime]/.
                             LorentzIndex[mudu[3],dime___] ->
                             LorentzIndex[mudu[4],dime]
                    ) +
                     (tdecex/.LorentzIndex[mudu[1],dime___] -> 
                              LorentzIndex[mudu[3],dime]/.
                              LorentzIndex[mudu[2],dime___] -> 
                              LorentzIndex[mudu[4],dime]
                     ) +
                     (tdecex/.LorentzIndex[mudu[1],dime___] -> 
                              LorentzIndex[mudu[4],dime]/.
                              LorentzIndex[mudu[2],dime___] -> 
                              LorentzIndex[mudu[3],dime] )
                   ];
     If[ ($LimitTo4 === True) && (tdeclpl == 3), 
         tdecr = epst[ tdecr, tdecnew, tdi, 1/12 ] 
       ];   (* e D0000 = 1/12 *)
     tdecr = add[ tdecr, PaVe[0,0,0,0, tensps, tdecml], tdecnew ];
    
                 (* the 00ij - terms *)
     tdecnew = Table[{Sort[{0,0,tdecti,tdectj}],   
                     (tdecex/.LorentzIndex[mudu[1],___] -> LorentzIndex[mudu[2],___]/. 
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectj]]
                     ) +
                     (tdecex/.LorentzIndex[mudu[2],___] -> LorentzIndex[mudu[3],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectj]]
                     ) +
                     (tdecex/.LorentzIndex[mudu[1],___] -> LorentzIndex[mudu[3],___]/. 
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectj]]
                     ) +
                     (tdecex/.LorentzIndex[mudu[1],___] -> LorentzIndex[mudu[4],___]/. 
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdectj]]
                     ) +
                     (tdecex/.LorentzIndex[mudu[2],___] -> LorentzIndex[mudu[4],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdectj]]
                     ) +
                     (tdecex/.LorentzIndex[mudu[3],___] -> LorentzIndex[mudu[4],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdectj]]
                     )},{tdectj,1,tdeclpl},{tdecti,1,tdeclpl}
                      ];

     tdecnew = eval[ Flatten[ tdecnew,1 ] ];
       For[ tdectj=1,tdectj<=Length[tdecnew],tdectj++,
            tdecr = add[ tdecr,
                    PaVe@@Join[tdecnew[[tdectj,1]],{tensps},{tdecml}],
                         tdecnew[[tdectj,2]]
                      ];
            If[$LimitTo4 ===True, tdecr = tdecr/.tdi->4]
          ];
     
       For[ tdecti = 1, tdecti <= tdeclpl, tdecti++,
            For[ tdectj = 1, tdectj <= tdeclpl, tdectj++,
                 For[ tdectk=1, tdectk <= tdeclpl, tdectk++,
                      For[ tdectl=1, tdectl <= tdeclpl, tdectl++,
                      tdecr = add[ tdecr,
                              PaVe@@Join[Sort[{tdecti,tdectj,tdectk,tdectl}],
                                              {tensps},{tdecml}
                                             ],
                           (eval[ tdecex/.LorentzIndex[mudu[1],___]->
                                          tdecpl[[tdecti]]/.
                                          LorentzIndex[mudu[2],___]->
                                          tdecpl[[tdectj]]/.
                                          LorentzIndex[mudu[3],___]->
                                          tdecpl[[tdectk]]/.
                                          LorentzIndex[mudu[4],___]->
                                          tdecpl[[tdectl]]
                                        ] 
                                   )
                                ];
            If[$LimitTo4 ===True, tdecr = tdecr/.tdi->4]
                ] ] ] ]

](*end if qn ==4 *);
If[qn === 5,
(*
Global`te=tdecex;
*)
(* here was a BUG indicated by Schorsch Weiglein  ... *)
   (* the 0000i - terms *)
lssu[exxp_, tdectii_][lmu1_, lmu2_, lmu3_, lmu4_, lmu5_] :=
      (exxp /. 
       LorentzIndex[mudu[lmu1], dime___] ->
       LorentzIndex[mudu[lmu2], dime]/.
       LorentzIndex[mudu[lmu3], dime___] ->
       LorentzIndex[mudu[lmu4], dime]/.
       LorentzIndex[mudu[lmu5], dime___] -> tdecpl[[tdectii]]
      ) +
      (exxp/.
       LorentzIndex[mudu[lmu1], dime___] ->
       LorentzIndex[mudu[lmu3], dime]/.
       LorentzIndex[mudu[lmu2], dime___] ->
       LorentzIndex[mudu[lmu4], dime]/.
       LorentzIndex[mudu[lmu5], dime___] -> tdecpl[[tdectii]]
      ) +
      (exxp/.
       LorentzIndex[mudu[lmu1], dime___] ->
       LorentzIndex[mudu[lmu4], dime]/.
       LorentzIndex[mudu[lmu2], dime___] ->
       LorentzIndex[mudu[lmu3], dime]/.
       LorentzIndex[mudu[lmu5], dime___] -> tdecpl[[tdectii]]
      );

   tdecnew =  Table[ lssu[tdecex, tdecti][1,2,3,4,5] +
                     lssu[tdecex, tdecti][2,3,4,5,1] +
                     lssu[tdecex, tdecti][3,4,5,1,2] +
                     lssu[tdecex, tdecti][4,5,1,2,3] +
                     lssu[tdecex, tdecti][5,1,2,3,4] , 
                    {tdecti,1,tdeclpl}
                   ];

   tdecnew = eval[ Flatten[ tdecnew,1 ] ];

       For[ tdec0j = 1, tdec0j <= tdeclpl, tdec0j++,
            tdecr = add[  tdecr, PaVe[0,0,0,0,tdec0j,tensps,tdecml],
                          tdecnew[[tdec0j]]  ]
          ];

   (* the 00ijk - terms *)
   tdecnew = Table[{Sort[{0,0,tdecti,tdectj,tdectk}],
                     (tdecex/.LorentzIndex[mudu[1],___] -> LorentzIndex[mudu[2],___]/.
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[5],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[1],___] -> LorentzIndex[mudu[3],___]/.
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[5],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[1],___] -> LorentzIndex[mudu[4],___]/.
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[5],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[1],___] -> LorentzIndex[mudu[5],___]/.
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[2],___] -> LorentzIndex[mudu[3],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[5],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[2],___] -> LorentzIndex[mudu[4],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[5],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[2],___] -> LorentzIndex[mudu[5],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[3],___] -> LorentzIndex[mudu[4],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[5],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[3],___] -> LorentzIndex[mudu[5],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[4],___] -> tdecpl[[tdectk]]
                     ) + 
                     (tdecex/.LorentzIndex[mudu[4],___] -> LorentzIndex[mudu[5],___]/.
                              LorentzIndex[mudu[1],___] -> tdecpl[[tdecti]]/.
                              LorentzIndex[mudu[2],___] -> tdecpl[[tdectj]]/.
                              LorentzIndex[mudu[3],___] -> tdecpl[[tdectk]]
                     )},{tdectk,1,tdeclpl},{tdectj,1,tdeclpl},{tdecti,1,tdeclpl}
                  ];

     tdecnew = eval[ Flatten[ tdecnew,1 ] ];
 
     tdecnew = eval[ Flatten[ tdecnew,1 ] ];
       For[ tdectj=1,tdectj<=Length[tdecnew],tdectj++,
            tdecr = add[ tdecr,
                    PaVe@@Join[tdecnew[[tdectj,1]],{tensps},{tdecml}],
                         tdecnew[[tdectj,2]]
                      ];
            If[$LimitTo4 ===True, tdecr = tdecr/.tdi->4]
          ];

       For[ tdecti = 1, tdecti <= tdeclpl, tdecti++,
            For[ tdectj = 1, tdectj <= tdeclpl, tdectj++,
                 For[ tdectk=1, tdectk <= tdeclpl, tdectk++,
                      For[ tdectl=1, tdectl <= tdeclpl, tdectl++,
                           For[ tdectm=1, tdectm <= tdeclpl, tdectm++,
                                tdecr = add[ tdecr,
                                 PaVe@@Join[Sort[{tdecti,tdectj,tdectk,tdectl,tdectm}],
                                                 {tensps},{tdecml}
                                           ],
                                       (eval[ tdecex/.LorentzIndex[mudu[1],___]->
                                          tdecpl[[tdecti]]/.
                                          LorentzIndex[mudu[2],___]->
                                          tdecpl[[tdectj]]/.
                                          LorentzIndex[mudu[3],___]->
                                          tdecpl[[tdectk]]/.
                                          LorentzIndex[mudu[4],___]->
                                          tdecpl[[tdectl]]/.
                                          LorentzIndex[mudu[5],___]->
                                          tdecpl[[tdectm]]
                                           ]
                                      )
                                         ];
            If[$LimitTo4 ===True, tdecr = tdecr/.tdi->4]
                ] ] ] ] ];

  ];
    
  tdecr = Expand[tdecr];

(* end if tdecex == 0 *) ];
print3["exiting tdec with ",tdecr];
tdecr]]; (*tdec end *)

(* ************************************************************** *)

(* ********************************************************************* *)
(*                          oneloop42                                    *)
(* ********************************************************************* *)

(* ---------------------------------------------------------------- *)
(* ---------------------------------------------------------------- *)
   nterms[x_Plus]:=Length[x];    (*ntermsdef *)
   nterms[x_]:=Block[{ntermslex = Expand[x]},
                     If[ Head[ntermslex]===Plus,
                         ntermslex = Length[ntermslex],
                         If[x===0, ntermslex = 0, ntermslex = 1]
                       ];
           ntermslex];
(* ------------------------------------------------------------ *)
(* ************************************************************************* *)
(* ************************************************************************* *)
(*             for cancelling q^2's and q.p's                         *)
(* ************************************************************************* *)
(* ************************************************************************* *)

(* ********************************************************************* *)
(*                          oneloop43                                    *)
(* ********************************************************************* *)


   dndummy[x__]:= dummy FeynAmpDenominator[x];
   vcid[x_,___]:=x;

(* For cancelling  q^2 's *)       (* Q2Canceldef *)

 Q2Cancel[yy_,q2c_,Q_]:=(*Feyn`Calc`Main`memset[ Q2Cancel[yy,q2c,Q],
                        *) Block[ {q2rule,q2res,q2sim,ey},

print3["entering Q2cancel ",yy," ",q2c];
      ey  = Expand2[ expandscalarproduct[ yy ],Q ];

    (* A simplification function *)
    If[ FreeQ[ey, Pair[Momentum[Q,di___],Momentum[Q,di___]] ],
        q2res = ey,
(* The rules of the game *)
   q2rule = {(
             (Pair[Momentum[Q,di___],Momentum[Q,di___]]^n_:1) FeynAmpDenominator[
              PropagatorDenominator[Momentum[Q,dime___],m1_],a___,
                    PropagatorDenominator[p_,m2_] ,y___ ] fa_.
             ) :> q2sim[ ( (fa Pair[ Momentum[Q,di],
                                     Momentum[Q,di]
                                   ]^(n-1) FeynAmpDenominator[
                             PropagatorDenominator[p,m2],a,y ]
                           )/.Q->(2 Q-(p/.Momentum->vcid/.vcid->Momentum))
                          ) + fa Pair[Momentum[Q,di],
                                      Momentum[Q,di]
                                     ]^(n-1) m1^2 FeynAmpDenominator[
                                   PropagatorDenominator[Momentum[Q,dime],m1],a,
                                      PropagatorDenominator[p,m2],y
                                                       ]
                       ] /; Length[{a}]===(q2c-1)
            };
   q2sim[x_]:= Expand[ expandscalarproduct[ spinorchainevaluate[x] 
                                       ]//denomExpand ] // smalld;
    q2res = ((ey//denomExpand)/.FeynAmpDenominator->dndummy/.
             dndummy->FeynAmpDenominator
            )//.q2rule;
    q2res = q2res/.dummy->1;
     q2res = q2sim[q2res]
      ](* end of if freeq  ey ... *);
    print3["exiting Q2Cancel ",q2res];
                                   q2res](*]*);

(* ************************************************************************* *)
(* ************************************************************************* *)

   QPCancel[eyin_, Q_Momentum]:=                     (* QPCanceldef *)
   Feyn`Calc`Main`memset[ QPCancel[ey,Q],
       Block[{qprule, qpsim, qpupec, qpre, q, translate, extractm, tran,
              feyndunique, fdunique, ey, eyfa = 1, qprf,nuLlLLL},
   nuLlLLL = Unique[System`C];
   q = Q[[1]];
   If[Head[eyin] === Times, 
      ey   =  Select[eyin, !FreeQ[#, q]&];
      eyfa =  eyin / ey,
      ey   =  eyin
     ];
   extractm[a_, ___] := a;
   tran[a_, x_, y_] := a /. (Rule @@ ( {x, y} /. Momentum -> extractm ));

  translate[x_] := Expand2[x, q] /. 
         { FeynAmpDenominator[a___, PropagatorDenominator[-Q+pe_., m0_], b___]:>
           FeynAmpDenominator[a, PropagatorDenominator[ Q - pe, m0], b],
          ( fa_. FeynAmpDenominator[PropagatorDenominator[Q + pe_, m0_], 
                                 PropagatorDenominator[Q + ka_, m1_], z___]
          )        :>
        tran[fa FeynAmpDenominator[PropagatorDenominator[Q + pe, m0], 
                                   PropagatorDenominator[Q + ka, m1], z],
                          Q, Q - pe ],
          ( fa_. FeynAmpDenominator[PropagatorDenominator[Q + pe_, m0_]] ):>
                tran[fa FeynAmpDenominator[PropagatorDenominator[Q + pe, m0]],
                     Q, Q - pe]
         };
(* create a unique order of the propagators using translational invariance *)
$FeynUnique = True;
   feyndunique[xx_] := If[ $FeynUnique =!= True, xx,
                           Map[ fdunique,
                    Collect2[xx, FeynAmpDenominator, Factoring -> False] +
                               nuLlLLL] /. nuLlLLL -> 0
                         ];
fdunique[xx_] := xx /; FreeQ[xx, FeynAmpDenominator];
Literal[fdunique[efc_. FeynAmpDenominator[a__] ]] :=
      Block[{old, new, icou = 0, i},
       new = fdun1[efc FeynAmpDenominator[a]];
print2["checkunique"];
       new[[1]]];
(* generate n integrals (for an n-point function) *)
   fdun1[0]=0;
   Literal[fdun1[ fc_. FeynAmpDenominator[a__] ]] := 
       Block[{ rl = {}, al,prd4, nal = {}, feynampdenO, fdsort},
   fdsort[x_] := FeynAmpDenominator[x];
   fdsort[x_, y__] := FeynAmpDenominator @@ Join[{x}, Sort[{y}]];
   prd4[xy_ , ma_] := PropagatorDenominator[
              ChangeDimension[xy , 4] /. Momentum[Q[[1]]] ->  Q, ma];
   al  = feynampdenO[a] /. PropagatorDenominator -> prd4;
(* seperately: FeynAmpDenominator-List *)
(* and coefficient list *)
     For[i = 1, i <= Length[al], i++,
         If[!MatchQ[al[[i, 1]], Q + pe_],
            AppendTo[nal, {fc, al /. FeynAmpDenominator -> fdsort}],
            AppendTo[nal, tran[ {fc, Prepend[Drop[al, {i,i}], al[[i]]]},
                                 Q,  2 Q - al[[i,1]] 
                                 ] /. PropagatorDenominator -> prd4 ];
           ];
        ];
      nal  = nal/. feynampdenO -> FeynAmpDenominator;
      nal = Join[nal, tran[nal, Q, -Q]] //.  
             FeynAmpDenominator[w___, PropagatorDenominator[-Q + pe_. , m0_], b___]:>
             FeynAmpDenominator[w, PropagatorDenominator[ Q - pe, m0], b];
      nal = nal /. FeynAmpDenominator -> fdsort;
      nal = Sort[Union[nal], OrderedQ[{#1[[2]], #2[[2]]}]&];
      nal = Map[(#[[1]] #[[2]])&, nal];
      nal];

(* A simplification function, which also translates (evtl.)  *)
   qpsim[x_]:=  Expand[ DiracSimplify[ 
                        x // smalld // translate 
                      ]//ExpandScalarProduct]  // smalld;

   qprule = {

(* -------------------------------------------------------------------- *)
(* hm, this is a good idea, but it does not work somehow ... 
   commented out 7.7.93 
             (
    Pair[p_Momentum, Momentum[q, di___]]  fa_. *
     FeynAmpDenominator[ a___, PropagatorDenominator[Q + p1_. , m1_], b___,
                               PropagatorDenominator[Q + p2_Plus, m2_], c___ ]
             ) :>  (
    (Pair[p2, Momentum[q, di]] fa *
      FeynAmpDenominator[a, PropagatorDenominator[Q + p1, m1], b, 
                            PropagatorDenominator[Q + p2, m2], c] - 
     qpsim[(Pair[p2, Momentum[q, di]] - 
            Pair[p, Momentum[q, di]])  fa *
      FeynAmpDenominator[a, PropagatorDenominator[Q + p1, m1], b,
                            PropagatorDenominator[Q + p2, m2], c]
          ]
    )              ) /; (First[p2] === p) || (First[p2] === (-p)),
*)
(* -------------------------------------------------------------------- *)

             (
   Pair[p_, Momentum[q, ___]]^n_.  fa_. * 
     FeynAmpDenominator[ a___, PropagatorDenominator[Q + p1_. , m1_], b___, 
                               PropagatorDenominator[Q + p2_, m2_], c___ ] 
             ):> 
       qpsim[ fa/2  Pair[p, Momentum[q]]^(n-1) Cancel[(p1 - p2)/p] ( 
                    FeynAmpDenominator[ a, b, 
                                       PropagatorDenominator[Q + p2, m2], c ] - 
                    FeynAmpDenominator[ a, PropagatorDenominator[Q + p1, m1], 
                                        b, c ] + 
                    FeynAmpDenominator[ a, PropagatorDenominator[Q + p1, m1], 
                                        b, PropagatorDenominator[Q + p2, m2], 
                                        c ] *
                    (-Pair[p1, p1] + Pair[p2, p2] - m2^2 + m1^2)
                                               )
            ] /; (p === p1 - p2) || (p === p2 - p1),

(* ************************************************************************* *)

   (* like Q2Cancel *)
   (Pair[Q, Q]^n_. FeynAmpDenominator[ a___, PropagatorDenominator[Q, m1_], 
                        y___ ] fa_. ) :> 
       (qpsim[  fa Pair[Q, Q]^(n-1) FeynAmpDenominator[ a, y ]
               + fa Pair[Q, Q]^(n-1) m1^2  * 
                   FeynAmpDenominator[ a, PropagatorDenominator[Q, m1], y ]
              ]
        ) /; Length[{a, y}] > 0,

(* ************************************************************************* *)
     FeynAmpDenominator[ PropagatorDenominator[Q, 0].. ] fa_ :> 0 /; FreeQ[{fa}, q],

     Pair[a_, Q]  FeynAmpDenominator[ PropagatorDenominator[Q, _] ] fa_ :> 
        0 /; FreeQ[{a,fa}, Q],
     Pair[Q, a_]  FeynAmpDenominator[ PropagatorDenominator[Q, _] ] fa_ :> 
        0 /; FreeQ[{a,fa}, Q]

(* ************************************************************************* *)

          } (* endqprule *);

   qprf[xy_] := If[Head[xy]===Times, 
                   Expand2[ (Select[xy, !FreeQ[#, q]&] /. qprule) *
                            Select[xy,  FreeQ[#, q]&], q
                         ],
                   xy /. qprule
                  ];
                   
   qpupec[ xxx_ ]:= If[Head[xxx] === Plus, Map[qprf, xxx], qprf[xxx]];
   qpupec[ xxx_ ]:= xxx/.qprule;

   qpre = FixedPoint[ qpupec,
                      Expand2[ey//translate//ExpandScalarProduct , q]
                    ];
   print3["uniqueness"];
   qpre = ExpandScalarProduct[feyndunique[qpre] // DiracSimplify];
   print2["againg in QPC"];
   qpre = FixedPoint[ qpupec, Expand2[qpre, q] ] // feyndunique //
                    DiracSimplify // ExpandScalarProduct;
   If[eyfa =!= 1, qpre = Expand2[qpre eyfa, q]];
   print2["againg in QPC done"];

print3["exiting qpc with ",qpre // FeynCalcForm];
      qpre]];

(* ************************************************************************* *)

(* ************************************************************************* *)
(* ************************************************************************* *)


(* ********************************************************************* *)
(*                          oneloop44                                    *)
(* ********************************************************************* *)

(* SetStandardMatrixElementdef *)
Options[SetStandardMatrixElements] = {WriteOut -> False};
 SetStandardMatrixElements[rx_List,en_:{}, op___Rule]:=
 Block[{links={},nmat,mat,ix,i,ii,j,sup,newli= {}, ops, enm,
          savmem,neweq,mati,set,isos,isolspc,nullll,x,x2,sumand,
          mati1,set2, filename},
  ops = {op}; enm = en;
  If[{ops}==={} || (!FreeQ[enm, WriteOut]), ops=en; enm = {}];
  filename = WriteOut /. ops /. Options[SetStandardMatrixElements];
  If[StringQ[filename], file = FileNames @@ {filename}];
     If[ValueQ[file] && (file =!= {}), 
        temp = Get @@ file;
        temp = Select[temp, !FreeQ[#,Spinor]&] 
       ];
     If[Length[temp]>0, temp/.Literal->Identity/.RuleDelayed->Set;
 print2["loading old matrixelementdefinitions"]
       ,

  savmem=$MemoryAvailable;
  $MemoryAvailable=0;
  x = {};
  For[ix = 1, ix <= Length[rx], ix ++,
      If[ FreeQ[rx[[ix,1]], Plus], 
          x = Prepend[x,{rx[[ix,1]],
                          StandardMatrixElement@@Flatten[{rx[[ix,2]]}]
                        } ],
          x = Append[x, {rx[[ix,1]],
                          StandardMatrixElement@@Flatten[{rx[[ix,2]]}]
                        } ]
        ]
     ];
  x = Flatten[x, 1];
 If[ Cases[x, Literal[Spinor[a_,_,_] . (___) . Spinor[b_,_,_] * 
                       Spinor[c_,_,_] . (___) . Spinor[d_,_,_] ] 
          ] =!= {},
    x = x /. DiracGamma[6] -> (1/2 + 1/2 DiracGamma[5]);
    x = x /. DiracGamma[7] -> (1/2 - 1/2 DiracGamma[5])
   ];


print2[Length[x]];
print2["enm = ",enm];
print2["ops= ",ops];
  If[ enm==={},
      mat = DiracSimplify[ x ]//Expand//DiracOrder//Contract,
      mat = DiracSimplify[ x/.enm ]//Expand//DiracOrder//Contract
    ];
       
   nmat = Expand[ expandscalarproduct[mat] ]//smalld;
   nmat = nmat /. Dot -> spinorsandpairs;
   isos[b_spinorsandpairs]:=b;
   isos[a_?NumberQ]:=a;
   isos[a_ b_spinorsandpairs]:= b isos[a];
   isolspc[xx_]:=Map[ isos, 
             collin[ xx, spinorsandpairs, True] + nullll
                    ]/.nullll->0;
   mat = Table[ {isolspc[ nmat[[2 ii - 1]] ],
                nmat[[2 ii]]},{ii,1,Length[nmat]/2}
              ];
   pat[x_,_]:=x;
(* Need this for pattern in the SME's, e.g., 4-fermion processes *)
   set2[a_, b_]:=Set @@ {a, b/.Pattern->pat};
 For[ i=1, i<=Length[mat],i++,
print2["i = ",i," out of ", Length[mat]];
      mat = Expand[ (mat//expandscalarproduct)/.spinorsandpairs->Dot ];
      mat = mat /. Dot -> spinorsandpairs;
      mati1 = Expand[isolspc[mat[[i,1]]]];
      For[j=1,j<=nterms[mati1],j++,
          If[ nterms[mati1] === 1,
              sumand = mati1,
              sumand = mati1[[j]]
            ];
        print2["sumand = ",sumand];
          If[!(FreeQ[sumand,spinorsandpairs]),
             sup = PartitHead[ sumand,spinorsandpairs]//Expand;
            If[ (!MemberQ[ links,sup[[2]]/.spinorsandpairs->Dot/.
                                           Pair -> bier ]) &&
                Head[ sup[[2]] ] === spinorsandpairs,
print2["o.k1"];
                links =  Append[ links, sup[[2]]/.spinorsandpairs->Dot/.
                                                  Pair -> bier
                               ]//Expand;
                neweq = set[ sup[[2]], Together/@(collin[
                             Expand[
                              ((mat[[i,2]]-mati1+sumand )/sup[[1]]
                                    )/.isos->Identity
                                   ],spinorsandpairs ,True
                                            ]
                                     )
                           ];
    (* Avoid things like  a=a *)
    If[ (neweq[[1]] - neweq[[2]]) =!= 0, 
print2["setting"];
        newli = Append[ newli,neweq/.set->set2 ]//Expand;
    j = nterms[mati1]+1
      ];
  ]    
 ]
     ]
      ] (* i - loop *);
 
print2["
Solving the system of linear equations for standard
matrix elements"];

(*
(* This takes care of the fact that  1 = Gamma6 + Gamma7 *)
If[ (!FreeQ[ rx, DiracGamma[6] ]) || (!FreeQ[ rx, DiracGamma[7] ]),
    SetStandardMatrixElements[
      rx/.DiracGamma[6]->(1/2 + DiracGamma[5]/2)/.
          DiracGamma[7]->(1/2 - DiracGamma[5]/2), enm]
  ];
*)
 $MemoryAvailable=savmem;
If[StringQ[filename], Put @@ {DownValues[spinorsandpairs], filename}]
];

newli];
(* ------------------------------------------------------------ *)

GetOneLoopResult[x_, li_List] := Block[{name, list,new, lenli},
name = ToString[x];
list = Table[StringJoin[name, "N",li[[i]]//ToString, ".m"], {i, Length[li]}];
none[__,y_] := ToExpression[ StringJoin[ Rest[ Characters[ToString[y]] ] ]];
new = 0;
lenli = Length[li];
For[j = 1, j <= lenli, j++, 
    print1["loading #  ",j, " out of ",lenli];
    Get[list[[j]]];
   ];
For[j = 1, j <= lenli, j++,
    print1["summing # ",j, " out of ",lenli];
    new = new + Expand2[ DownValues[OneLoopResult][[j,2]],
                        StandardMatrixElement ];
   ];
new];
WriteString["stdout", "."];
End[];

Begin["Feyn`Calc`PaVe`"];
(* ***************************************************************** *)
(*                          pave10                                   *)
(* ***************************************************************** *)
SetAttributes[memsetP, HoldFirst];
collect2P[x__] := Collect2[x, Factoring -> False];
collect3[x__]  := Collect2[x, Factoring -> True];
combine[x__]   :=Combine[x, ProductExpand -> False ];
memsetP[x_,y_] :=Feyn`Calc`Main`memset[x, y];
nterms2[x_]    :=Feyn`Calc`Main`nterms[x];
(* ***************************************************************** *)
(*                          pave11                                   *)
(* ***************************************************************** *)
 Negligible[0] = 0;                                       (*Negligibledef*)
 Negligible[x_^pow_] := Negligible[x]^pow;
(* ***************************************************************** *)
(*                          pave12                                   *)
(* ***************************************************************** *)
(* Symmetry in the indices *)
 PaVe[i_,j__,  pl_List, ml_List ]  := PaVe @@Join[Sort[{i,j}],{pl,ml}]/;
                                      !OrderedQ[{i,j}];
(* Special cases of PaVe: *)
 PaVe[0, {}, {x_}]      := A0[x];
 PaVe[0, {p2}, {x_,y_}] := B0[p2,x,y];
(* A0def *)
 Options[A0] = {A0ToB0 -> True};    
 A0[Negligible[_]^_. , ___] := 0;  (* In dimensional regularization: A0(0)=0 *)
 A0[small[_]^_. ,___] := 0;
 A0[0,___] = 0;
 A0[mm_, op___]:=(mm + mm B0[0,mm,mm])/; ( A0ToB0/.{op}/.Options[A0] ) &&
                                       !( BReduce/.Options[B0] );
(* there is no tensorial 1-point function *)
 PaVe[_,{},{_}] := 0;
(* but a non-zero coefficient of g_munu *)
 PaVe[0,0,{},{m2_}] := (m2/4 A0[m2] + m2^2/8) /; $LimitTo4 === True;
(* ***************************************************************** *)
(*                          pave13                                   *)
(* ***************************************************************** *)
(* B0def*)
 Options[B0] = { BReduce -> False, B0Unique -> False, B0Real -> False };
 PaVe[0, {p_}, {m1_, m2_}] := B0[p, m1, m2];
 B0[pe_,me2_,me1_,opt___]:=
    B0 @@ Prepend[ {me1,me2,opt}, Expand[pe]] /; !OrderedQ[{me2,me1}];
 B0[Negligible[pp_]^j_., Negligible[a_]^n_., Negligible[b_]^m_.] := B0[pp^j, a^n, b^m];
 B0[0, Negligible[a_]^n_., Negligible[b_]^m_.] := B0[0, a^n, b^m];
 bop[x___] := BReduce/.Flatten[ Join[{x},Options[B0]] ];
 nos[x_] := True/;(x=!=0) && FreeQ[x,Negligible]&&FreeQ[x,small];

 nos[x_] := If[(x =!= 0) && FreeQ[x, Negligible] && FreeQ[x, small],
               True, False];

(* Ansgar does not like this, therefore the default of B0Unique is False ...*)
 B0[0,0,mm_,opt___]:=( B0[0,mm,mm] + 1 ) /; nos[mm] && 
          ( (B0Unique/.{opt}/.Options[B0]) === True );
 B0[mm_,0,mm_,opt___]:=( B0[0,mm,mm] + 2 ) /; 
          ( (B0Unique/.{opt}/.Options[B0]) === True ) && 
          ( (B0Real/.{opt}/.Options[B0]) === True );

 B0[Negligible[pp_]^n_., m1_, m2_, opt___] := B0[0, m1, m2, opt]/;nos[m1]||nos[m2];
 B0[pp_, Negligible[m1_]^n_., m2_, opt___] := B0[pp, 0, m2, opt]/;nos[pp]||nos[m2];
 B0[pp_, m1_, Negligible[m2_]^n_.,opt___ ] := B0[pp, m1, 0, opt]/;nos[pp]||nos[m1];
(* ***************************************************************** *)
(*                          pave14                                   *)
(* ***************************************************************** *)
 (* DB0def *)
 If[$VersionNumber > 2.2, Unprotect[Derivative]];
 Derivative[1, 0, 0][B0][pp_,m02_,m12_]=DB0[pp,m02,m12];
 (* also DB0 is symmetric in its mass arguments *)
  DB0[pe_,me2_,me1_,opt___]:=
     DB0 @@ Prepend[ {me1,me2,opt}, Expand[pe]] /; !OrderedQ[{me2,me1}];
 sull[_] := 0;
 DB0[pe_,me1_,me2_] := (DB0 @@ ({pe,me1,me2}/.Negligible->sull)) /; nos[pe] &&
                             !FreeQ[{me1,me2}, Negligible];
 (* Derivative of B1 *)
 bop1[x___] := BReduce/.Flatten[ Join[{x},Options[DB1]] ];
 Options[DB1]={ BReduce -> False };
 Derivative[1, 0, 0][B1][pp_,m02_,m12_]=DB1[pp,m02,m12];
 DB1[m_, m_, 0, opt___] := (- DB0[m,m,0] + 1/2/m) /; bop1[opt];
DB1[pp_, m02_, m12_, opt___]:=(
- (m12 - m02)/(2 pp^2) ( B0[pp,m02,m12] - B0[0,m02,m12] ) +          
(m12 - m02 - pp)/(2 pp) DB0[pp,m02,m12] ) /; nos[pp] && bop1[opt];

(* special cases *)
 B0[kl_, kmm_, mm_, opt___ ] :=  (A0[mm]/mm) /; nos[mm] && bop[opt] &&
                                         (((kl E+kmm)/.Negligible[_]->0)===0);
 B0[kl_, m_, kmm_, opt___ ] :=  (A0[mm]/mm) /; nos[mm] && bop[opt] &&
                                         (((kl E+kmm)/.Negligible[_]->0)===0);

 B0[0,mm_,mm_,opt___]        := (A0[mm]/mm - 1)/;nos[mm] && bop[opt];
 B0[Negligible[_]^n_.,mm_,mm_,opt___] := (A0[mm]/mm - 1)/;nos[mm] && bop[opt];
(* smaddef *)
  smanull[_]:=0;
  smad[x_]:=Block[{nx=Factor2[x]/.Negligible->smanull},
                   Factor2[Numerator[nx]]/ Factor2[Denominator[nx]]
                 ];
(* ***************************************************************** *)
(*                          pave15                                   *)
(* ***************************************************************** *)
(* B1def *)
 Options[B1] = {BReduce->True};
 PaVe[1,{pp_},{mm1_,mm2_}]  := B1[pp, mm1, mm2];
 B1[a_,b_,c_,ops___Rule] :=  bb1[a, b, c] /; 
                 ((BReduce/.{ops}/.Options[B1])===True) && 
                  Head[bb1[a,b,c]] =!= bb1;
(* Special cases, if photon and fermionic small masses are present *)
 bb1[Negligible[me_]^n_., Negligible[me_]^n_., Negligible[mla_]^m_.]:=
   ( -1/2 B0[Negligible[me]^n, Negligible[me]^n, 0] - 1/2 )/; TrueQ[mla < me];

 bb1[Negligible[me_]^n_., Negligible[mla_]^n_., Negligible[me_]^m_.]:=
   (1/2 - 1/2 B0[Negligible[me]^n,0 ,Negligible[me]^n]) /; TrueQ[mla < me];

(* other special cases of B1 *)

(* B1( p,m,m ) = -1/2 B0( p,m,m )  *)
 bb1[pp_,mm_,mm_] := -1/2 B0[pp,mm,mm];
 bb1[mm_, mm_, 0]:= -1/2 B0[mm, mm, 0] - 1/2;
 bb1[mm_, 0, mm_]:= 1/2 - B0[mm,0,mm]/2;
 bb1[0,0,mm_]:=-1/2 B0[0,0,mm]+1/4;
 bb1[Negligible[_]^n_.,0,mm_]:=( -1/2 B0[0,0,mm] + 1/4 )/;nos[mm];
 bb1[0,Negligible[_]^n_.,mm_]:=( -1/2 B0[0,0,mm] + 1/4 )/;nos[mm];
 bb1[0,mm_,0]           :=( -1/2 B0[0,0,mm] - 1/4 )/;nos[mm];

 bb1[Negligible[_]^n_.,Negligible[_]^n_.,mm_]:=( -1/2 B0[0,0,mm] + 1/4 )/;nos[mm];
 bb1[Negligible[_]^n_.,mm_,Negligible[_]^n_.]:=( -1/2 B0[0,0,mm] - 1/4 )/;nos[mm];
 bb1[Negligible[_]^n_.,mm_,0]:=( -1/2 B0[0,0,mm] - 1/4 )/;nos[mm];
(* B1 in general *)
 bb1[pp_,ma0_,ma1_ ]:=(smad[ma1-ma0]/(2 pp) (B0[pp,ma0,ma1] -
                                            B0[0,ma0,ma1]) - 
                        1/2 B0[pp,ma0,ma1]
                      ) /; nos[pp];
(* ***************************************************************** *)
(*                          pave16                                   *)
(* ***************************************************************** *)
(* B00def *)
 PaVe[0,0,{p_},{m1_,m2_}]  := B00[p,m1,m2] /; $LimitTo4 === True;;
 Options[B00]={BReduce->True};
 B00[x__,  BReduce->True]:= b00[x] /; $LimitTo4 === True;
 B00[x__,  BReduce->True]:= PaVeReduce[PaVe[0,0,{First[{x}]},Rest[{x}]]] /; 
         $LimitTo4 === False;
 B00[x_,y_,z_]:= b00[x,y,z]/;( BReduce/.Options[B00] )===True && 
                             $LimitTo4 === True;
 B00[x_,y_,z_]:= B00[x,y,z, BReduce->True]/;( BReduce/.Options[B00] )===True && 
                             $LimitTo4 === False;
 b00[0,mm_,mm_] := mm / 2 ( B0[0,mm,mm] + 1 )/;nos[mm];
 b00[Negligible[_]^n_.,mm_,mm_] := mm / 2 ( B0[0,mm,mm] + 1 )/;nos[mm];
 b00[pp_,mm_,mm_]     :=  1/6 ( A0[mm]+B0[pp,mm,mm] smad[2 mm - pp/2] +
                                  smad[2 mm - pp/3]) /;nos[pp];
 b00[pp_,mm1_,mm2_]         :=  ( 1/6 ( A0[mm2]+
                               (B1[pp,mm1,mm2] ) smad[pp-mm2+mm1] )+
                                     smad[mm1/3] B0[pp,mm1,mm2] +
                                smad[ 1/6 ( mm1 + mm2 - pp/3 ) ] );

(* ***************************************************************** *)
(*                          pave17                                   *)
(* ***************************************************************** *)
(*                              B11def                                *)
(* ****************************************************************** *)
 Options[B11]={BReduce->True};
 PaVe[1,1,{pp_},{mm1_,mm2_}]  := B11[pp,mm1,mm2] /; $LimitTo4 === True;
 B11[pe_, mm1_, mm2_,  BReduce->True]     := 
     b11[pe, mm1, mm2] /; ($LimitTo4 === True) && 
                          (nos[pe] || ( (!nos[pe]) && (mm1 === mm2))) ;

 B11[x__,  BReduce->True]     := 
   PaVeReduce[PaVe[1,1,{{x}[[1]]}, Rest[{x}] ]] /; 
      ($LimitTo4 === False) && (x[[1]] =!=0);
 B11[x_,y_,z_]:= b11[x,y,z]/; (( BReduce/.Options[B11] )===True )&& 
                             ($LimitTo4 === True ) &&
                             (nos[x] || ( (!nos[x]) && (y === z)));
 B11[x_,y_,z_]:= B11[x,y,z, BReduce->True]/;( BReduce/.Options[B11] )===True && 
                             ($LimitTo4 === False) && nos[x];
 b11[ 0,mm1_,mm1_ ] := 1/3 * B0[ 0,mm1,mm1 ];
 b11[ Negligible[_]^n_.,mm1_,mm1_ ] := 1/3 * B0[ 0,mm1,mm1 ];
 b11[ pp_,mm_,mm_]:= ( 1/(3pp) ( A0[mm]+B0[pp,mm,mm] smad[pp-mm]-   
                                         smad[mm - pp/6] )) /;nos[pp];
 b11[ pp_,m1_,m2_ ] := ( 1/(3 pp) ( A0[m2] - smad[2 (pp-m2 + m1)]*
                       (PaVe[1,{pp},{m1,m2}]) - smad[m1] B0[pp,m1,m2] - 
                        smad[ 1/2 (m1 + m2 - pp/3 )]) )/;nos[pp];
(* ***************************************************************** *)
(*                          pave18                                   *)
(* ***************************************************************** *)
(* Symmetries of the PaVe - integrals *)
(* ****************************************************************** *)
(* Notation :   p10 = p1^2;  p12 = (p1-p2)^2;  etc.                   *)
(* ****************************************************************** *)
(* C2 --> C1, C22 --> C11,  C002 --> C001, C222 --> C111,   *)
(* if p10=p20  and  m2=m3    *)
PaVe[2,{p10_,p12_,p10_},{m1_,m2_,m2_}]  :=PaVe[1,{p10,p12,p10},{m1,m2,m2}];
PaVe[2,2,{p10_,p12_,p10_},{m1_,m2_,m2_}]:=PaVe[1,1,{p10,p12,p10},{m1,m2,m2}];
PaVe[0,0,2,{p10_,p12_,p10_},{m1_,m2_,m2_}]:=PaVe[0,0,1,{p10,p12,p10},{m1,m2,m2}];
PaVe[1,2,2,{p10_,p12_,p10_},{m1_,m2_,m2_}]:=PaVe[1,1,2,{p10,p12,p10},{m1,m2,m2}];
PaVe[2,2,2,{p10_,p12_,p10_},{m1_,m2_,m2_}]:=PaVe[1,1,1,{p10,p12,p10},{m1,m2,m2}];
(* a special case *)
 PaVe[ 2,{p10_, pp_,pp_},{m_,m_,m2_} ]:=
  - 2 PaVe[1,{p10,pp,pp},{m,m,m2}] - PaVe[0,{p10,pp,pp},{m,m,m2}];
(* *********************************************************************** *)
(*  D's: The argument list is (in general) : p10, p12, p23, p30, p20, p13  *)
(* *********************************************************************** *)
 pav[{a__},pl_List,ml_List]:=PaVe[a,pl,ml];
(*  1 <---> 2;   p20=p10,  p23=p13 , m3 = m2  *)
 PaVe[x__,{p10_,p12_,p13_,p30_,p10_,p13_},{m1_,m2_,m2_,m4_}]:=
  pav[{x} /. {1:>2, 2:>1}, {p10,p12,p13,p30,p10,p13},{m1,m2,m2,m4} ]/;
   Count[{x}, 2] > Count[{x}, 1];

(*  1 <---> 3;   p10=p30,  p12=p23 , m2 = m4  *)     
 PaVe[x__,{p10_,p12_,p12_,p10_,p20_,p13_},{m1_,m2_,m3_,m2_}]:=
    pav[{x} /. {1:>3, 3:>1}, {p10,p12,p12,p10,p20,p13},{m1,m2,m3,m2} ]/;
     Count[{x}, 3] > Count[{x}, 1];

(*  2 <---> 3;   p30=p20,  p13=p12 , m3 = m4  *)
 PaVe[x__,{p10_,p12_,p23_,p20_,p20_,p12_},{m1_,m2_,m3_,m3_}]:=
  pav[{x} /. {2:>3, 3:>2}, {p10,p12,p23,p20,p20,p12},{m1,m2,m3,m3}]/;
     Count[{x}, 3] > Count[{x}, 2];

(* in order to canonize the C0's  (args:   p1^2, (p2-p1)^2, p2^2)  *)
 PaVe[0, {p10_, p12_, p20_}, {m1_, m2_, m3_}] := cord[p10, p12, p20,m1,m2,m3];
 PaVe[0, {p10_, p12_, p23_, p30_, p13_, p20_}, {m1_, m2_, m3_, m4_}]:=
   D0[p10, p12, p23, p30, p13, p20, m1, m2, m3, m4](*//PaVeOrder*);


 cord[a_,b_,c_, m1_,m2_,m3_]:=
     C0@@( Sort[{ {a,b,c, m1,m2,m3}, {c,b,a, m1,m3,m2},
                  {a,c,b, m2,m1,m3}, {b,c,a, m2,m3,m1},
                  {c,a,b, m3,m1,m2}, {b,a,c, m3,m2,m1} } ][[1]] );

   cord[C0[six__],{}]:=cord[six];
   cord[C0[te__], argu_List ]:= Block[{int, puref, arg, smalist, six,
                                       varg, sma, pw},
       six =  {te}/. smallLL->sma;
       If[FreeQ[six, sma],
          arg = argu,
          smalist = Select[Variables[six/.Power->pw], (!FreeQ[#, sma])&]/.pw->Power;
          If[!FreeQ[smalist, Power],
             arg = (argu/.smallLL->Identity) /.
                   Map[(#[[1,1]] -> (#[[1]]) )&, smalist ],
             arg = argu/.smallLL->sma
            ];
         ];
       varg = Variables[arg];
       For[iv=1,iv<=Length[varg],iv++,
           If[(!FreeQ[six, varg[[iv]]^2]) && FreeQ[arg,varg[[iv]]^2],
              arg = arg/.varg[[iv]]->(varg[[iv]]^2)];
          ];
       puref = func[Apply[or,(stringmatchq[slot[1], #]& /@ tomatch[arg])
                         ]]/.slot->Slot/.func->Function/.or->Or/.
                          stringmatchq->StringMatchQ;
       int = Select[ tostring /@ (oldper@@six),
                     func[ stringmatchq[slot[1],tomatch[arg]]
                         ]/.slot->Slot/.func->Function/.
                           stringmatchq->StringMatchQ
                          ];
       If[Length[int] === 0, int = six,int=ToExpression[int[[1]]]];
       int/.sma->smallLL] /; Length[{te}]===6 && Length[argu]>0;

(* ***************************************************************** *)
(*                          pave19                                   *)
(* ***************************************************************** *)
(* D0def *)
 Options[PaVeOrder] = {PaVeOrderList -> {}};

(* smallLL is intermediately introduced for Negligible *)
(* PaVeOrderdef *)
 PaVeOrder[expr_,opt___Rule]:=Block[{new, dordering, opli, cordering,
             be0, be1, aa0, be11, be00, dordering0,j,nulL },
   opli = PaVeOrderList/.{opt}/. Options[PaVeOrder];
   If[opli === False, new = expr,
   new = expr/.B0->be0/.B1->be1/.A0->aa0/.B00->be00/.B11->be11/.
         Negligible->smallLL/. 0->nulL;
   opli = opli /. 0 -> nulL /. Negligible->smallLL;
   dordering0[ten__]:=(D0@@(oldper[ten][[1]]));
   If[ Length[opli]>0, 
       If[ Head[opli[[1]]]=!=List, opli = {opli}];
       If[expr=!=(D0@@opli[[1]]),
          new = new /. D0 -> dordering0;
          For[j=1, j<=Length[opli], j++,  
              dordering[j][ten10__]:= D0 @@ dord[ D0[ten10], opli[[j]]];
              cordering[j][six06__]:= C0 @@ cord[ C0[six06], opli[[j]]];
              new = new/.D0->dordering[j]/.C0 -> cordering[j]
             ]
        ],
       new = new /. D0 -> dordering0 /. C0 -> cord;
      ];
   new = new(*/.C0->cord*) /. nulL -> 0 /. smallLL -> Negligible;
   new = new/.be0->B0/.be1->B1/.aa0->A0/.be00->B00/.be11->B11;
    ];
    new];

(* Make use of the nice new StringReplace *)
   
   tostring = ToString[InputForm[#], PageWidth -> 4711]&;
   tomatch[{li:{__}..}]:= tomatch /@ {li};
   tomatch[{li__}]:=StringReplace[tostring[{li}],{"{"->"*","}"->"*"}]/;
                                  Head[{li}[[1]]]=!=List;
   dord[D0[ten__],{}]:=dord[D0[ten]];
   dord[D0[te__], argu_List ]:= Block[{int, puref, arg, smalist, ten, 
                                       varg, sma, pw},
       ten =  {te}/. smallLL->sma;
       If[FreeQ[ten, sma], 
          arg = argu,
          smalist = Select[Variables[ten/.Power->pw], (!FreeQ[#, sma])&]/.pw->Power;
          If[!FreeQ[smalist, Power], 
             arg = (argu/.smallLL->Identity) /. 
                   Map[(#[[1,1]] -> (#[[1]]) )&, smalist ],
             arg = argu/.smallLL->sma
            ];
         ];
       varg = Variables[arg];
       For[iv=1,iv<=Length[varg],iv++,
           If[(!FreeQ[ten, varg[[iv]]^2]) && FreeQ[arg,varg[[iv]]^2], 
              arg = arg/.varg[[iv]]->(varg[[iv]]^2)];
          ];
       puref = func[Apply[or,(stringmatchq[slot[1], #]& /@ tomatch[arg])
                         ]]/.slot->Slot/.func->Function/.or->Or/.
                          stringmatchq->StringMatchQ;
       int = Select[ tostring /@ (oldper@@ten), 
                     func[ stringmatchq[slot[1],tomatch[arg]]
                         ]/.slot->Slot/.func->Function/.
                           stringmatchq->StringMatchQ 
                          ];
       If[Length[int] === 0, int = ten,int=ToExpression[int[[1]]]];
       int/.sma->smallLL] /; Length[{te}]===10 && Length[argu]>0;

(* If no ordering list is given, a standard representative is returned *)
 dord[D0[ten__]]:=(oldper[ten][[1]])/;Length[{ten}]===10;

oldper[a_,b_,c_, m1_,m2_,m3_] := 
Sort[{ {a,b,c, m1,m2,m3}, {c,b,a, m1,m3,m2},
                  {a,c,b, m2,m1,m3}, {b,c,a, m2,m3,m1},
                  {c,a,b, m3,m1,m2}, {b,a,c, m3,m2,m1} } 
    ];

(* This list has been calculated with FeynCalc! *)
oldper[p10_,p12_,p23_,p30_,p20_,p13_,m0_,m1_,m2_,m3_]:=Sort[{
   {p10, p12, p23, p30, p20, p13, m0, m1, m2, m3},
   {p10, p13, p23, p20, p30, p12, m0, m1, m3, m2},
   {p20, p12, p13, p30, p10, p23, m0, m2, m1, m3},
   {p20, p23, p13, p10, p30, p12, m0, m2, m3, m1},
   {p30, p13, p12, p20, p10, p23, m0, m3, m1, m2},
   {p30, p23, p12, p10, p20, p13, m0, m3, m2, m1},
   {p10, p20, p23, p13, p12, p30, m1, m0, m2, m3},
   {p10, p30, p23, p12, p13, p20, m1, m0, m3, m2},
   {p12, p20, p30, p13, p10, p23, m1, m2, m0, m3},
   {p12, p23, p30, p10, p13, p20, m1, m2, m3, m0},
   {p13, p30, p20, p12, p10, p23, m1, m3, m0, m2},
   {p13, p23, p20, p10, p12, p30, m1, m3, m2, m0},
   {p20, p10, p13, p23, p12, p30, m2, m0, m1, m3},
   {p20, p30, p13, p12, p23, p10, m2, m0, m3, m1},
   {p12, p10, p30, p23, p20, p13, m2, m1, m0, m3},
   {p12, p13, p30, p20, p23, p10, m2, m1, m3, m0},
   {p23, p30, p10, p12, p20, p13, m2, m3, m0, m1},
   {p23, p13, p10, p20, p12, p30, m2, m3, m1, m0},
   {p30, p10, p12, p23, p13, p20, m3, m0, m1, m2},
   {p30, p20, p12, p13, p23, p10, m3, m0, m2, m1},
   {p13, p10, p20, p23, p30, p12, m3, m1, m0, m2},
   {p13, p12, p20, p30, p23, p10, m3, m1, m2, m0},
   {p23, p20, p10, p13, p30, p12, m3, m2, m0, m1},
   {p23, p12, p10, p30, p13, p20, m3, m2, m1, m0}       }];

(* ***************************************************************** *)
(*                          pave20                                   *)
(* ***************************************************************** *)
(* These are the ultimate formulas for the reduction of coefficient 
   functions. The reference is: Techniques for the calculation of 
   electroweak radiative correctinos at the one-loop level and ...,
   by A. Denner (slightly rewritten by R.M), 
   to appear in Fortschritte der Physik, 1993 
*) 
(* ***************************************************************** *)
(* Notation :    pij = (pi - pj)^2;  where p0 = (0,0,0,0), i.e., 
                 p10 = p1^2, etc.  *)
(* ***************************************************************** *)
PaVeBr[i__, p_List, m_List] := tT[Length[m]][i][Join[p, m]];
(*breakdowndef*)
breakdown[x_]:= If[FreeQ[x,PaVe], x,
                   FixedPoint[(#/.T->tT)&,x/.PaVe->PaVeBr] ];
drop[] = {}; (* i.e. no index is an empty list *)
drop[x__] := Drop[{x},-1];
(* A Kronecker delta *)
delt = If[ #1 === #2, 1, 0]&;
(* ***************************************************************** *)
(*                          pave21                                   *)
(* ***************************************************************** *)

(* This is only valid for UV - Divergences !! *)
$epsilon /: $epsilon^n_Integer?Positive := 0;
$epsilon /: $epsilon A0[mm_] := 2 mm;
$epsilon /: $epsilon B0[_, _, _] := 2;
$epsilon /: Literal[$epsilon B1[_, _, _]] := -1;
$epsilon /: $epsilon C0[__] :=0;
$epsilon /: $epsilon D0[__] :=0;
$epsilon /: $epsilon T[3][1][_] :=0;
$epsilon /: $epsilon T[3][2][_] :=0;
$epsilon /: $epsilon T[4][ij__][_] :=0 /; Length[{ij}] < 4;

(* ***************************************************************** *)
(*                          pave22                                   *)
(* ***************************************************************** *)
(* Things with head Negligible are discarded in this function *)
SetAttributes[demon, Listable];
Literal[demon[demon[x_]]] := demon[x];
null[_]:=0;
demon[x_]:=memsetP[
               demon[x],
                 Block[{nx=x, den}, (* This is quite tricky ... *)
                       den =  Factor2[Denominator[nx]];
                       If[ Head[den]=== Plus,
                           nx = nx /. Negligible->null,
                           (* Now:if there is something small in the numerator*)
                           nx = Factor2[nx] /. Negligible->null ];
                        Factor2[Numerator[nx]]/ Factor2[Denominator[nx]]
                      ]
                  ];
(* Remember:  pij = (pi - pj)^2, i.e, p12 = 1/2 ( p10+p20-p12 ), 
	      with  p0 = (0,0,0,0)
*)
kinmainv[2, {p10_, p12_, p20_, m02_,m12_,m22_}] := 
 1/demon[ p10 p20 -  (1/2 (p10+p20-p12) )^2 ] *
  demon[{ {p20, -1/2 (p10+p20-p12)}, {-1/2 (p10+p20-p12), p10} }];

(* Calculate determinants only once. *)
det[x_List] := det[x] = demon[ demon[Det[demon[x]]] // Factor2 ];

(* ***************************************************************** *)
(*                          pave23                                   *)
(* ***************************************************************** *)

kinmainv[3, {p10_, p12_, p23_, p30_, p20_, p13_,m02_, m12_, m22_, m32_}] :=
  Block[{p1p2, p1p3, p2p3},
   {p1p2, p1p3, p2p3} = demon[{p10-p12+p20, p10-p13+p30, p20-p23+p30}/2];
   1/det[{ {p10, p1p2, p1p3}, {p1p2, p20, p2p3}, {p1p3, p2p3, p30} }] *
   demon[
        {{-p2p3^2 + p20*p30, p1p3*p2p3 - p1p2*p30, p1p2*p2p3 - p1p3*p20},
         {p1p3*p2p3 - p1p2*p30, -p1p3^2 + p10*p30, p1p2*p1p3 - p10*p2p3},
         {p1p2*p2p3 - p1p3*p20, p1p2*p1p3 - p10*p2p3, -p1p2^2 + p10*p20}}
 ]      ];

(* Xinvdef  !!!!! No precaution is taken for 0 determinants !!!! *)
(* Only for B1 and B11 the Xinv[1][1,1] will not be used *)
Xinv[1][1,1][{pp_,_,_}] := 1/pp;
Xinv[2][i_, j_][a_]:= Xinv[2][i,j][a] = kinmainv[2, a][[i,j]];
Xinv[3][i_, j_][a_]:= Xinv[3][i,j][a] = kinmainv[3, a][[i,j]];

(* we put the Negligible - demon here *) 
f[1][{pp_, m02_,m12_}] := demon[ pp - m12 + m02 ]; 
f[1][{p10_,p12_,p20_, m02_,m12_,m22_}] := demon[ p10 - m12 + m02 ]; 
f[2][{p10_,p12_,p20_, m02_,m12_,m22_}] := demon[ p20 - m22 + m02 ];
f[1][{p10_, p12_, p23_, p03_, p20_, p13_, m02_, m12_, m22_, m32_}] :=
 demon[ p10 - m12 + m02 ];
f[2][{p10_, p12_, p23_, p03_, p20_, p13_, m02_, m12_, m22_, m32_}] :=
 demon[ p20 - m22 + m02 ];
f[3][{p10_, p12_, p23_, p03_, p20_, p13_, m02_, m12_, m22_, m32_}] :=
 demon[ p03 - m32 + m02 ];

T[0][][___] := 0; (* in dimensional regularization *)
T[1][][{mm_}] := A0[mm];
T[2][][{pp_, m12_, m22_}] := B0[pp, m12, m22];
T[3][][{p10_,p12_,p20_, m02_,m12_,m22_}] := C0[p10,p12,p20,m02,m12,m22];
T[4][][{p10_, p12_, p23_, p03_, p20_, p13_, m02_, m12_, m22_, m32_}] :=
                   D0[p10,p12,p23,p03,p20,p13, m02,m12,m22,m32];

(* ***************************************************************** *)
(*                          pave24                                   *)
(* ***************************************************************** *)

(* The translated argument lists obtained by canceling  *)
c[0][{pp_, m12_, m22_}] := {m22};
c[1][{pp_, m12_, m22_}] := {m12};
c[0][{p10_,p12_,p20_, m02_,m12_,m22_}] := {p12, m22, m12};
c[1][{p10_,p12_,p20_, m02_,m12_,m22_}] := {p20, m02, m22};
c[2][{p10_,p12_,p20_, m02_,m12_,m22_}] := {p10, m02, m12};
c[0][{p10_, p12_, p23_, p03_, p20_, p13_, m02_, m12_, m22_, m32_}] :=
     {p13, p12, p23, m32, m12, m22};
c[1][{p10_, p12_, p23_, p03_, p20_, p13_, m02_, m12_, m22_, m32_}] :=
     {p20, p23, p03, m02, m22, m32};
c[2][{p10_, p12_, p23_, p03_, p20_, p13_, m02_, m12_, m22_, m32_}] :=
     {p10, p13, p03, m02, m12, m32};
c[3][{p10_, p12_, p23_, p03_, p20_, p13_, m02_, m12_, m22_, m32_}] :=
     {p10, p12, p20, m02, m12, m22};

(* getmdef *)
getm[{_}] := 0; getm[{_,_,_}] := 1; getm[{_,_,_,_,_,_}] := 2;
getm[{_,_,_,_,_,_,_,_,_,_}] := 3
(* ***************************************************************** *)
(*                          pave25                                   *)
(* ***************************************************************** *)
pluep2[x__]:=Plus[x]/;!FreeQ2[{x}, {tT,B0,B1,B00,B11,C0,D0,T}];

(* equation  (4.18) of A. Denners review *)
tT[N_Integer][0,0, i___Integer][a_List] := Block[{P, M, k, epsi },
  P = 2 + Length[{i}]; M = getm[a];
       1/(2 + P - M) (R[N, 0, 0][i][a] - Sum[ R[N, k][k, i][a], {k, M}])+
       Expand[( 1/(2 + P - M) epsi/(2 + P - M)* 
                     (R[N, 0, 0][i][a] - Sum[ R[N, k][k, i][a], {k, M}])
              )/.Plus->pluep2 ]/.epsi->$epsilon/.pluep2->Plus
                                                ] /; $LimitTo4 === True;

(* $dIM gets defined from the option of PaVeReduce *)
tT[N_Integer][0,0, i___Integer][a_List] := Block[{P, M, k, epsi },
  P = 2 + Length[{i}]; M = getm[a];
       1/($dIM + P -2 - M) (R[N, 0, 0][i][a] - Sum[ R[N, k][k, i][a], {k, M}])
                                                ] /; $LimitTo4 =!= True;
(* two special cases *)
R[n0__][j__][a_] := ( R[n0][Sequence @@ Sort[{j}]][a] ) /; !OrderedQ[{j}];
R[n0__][{}][a_] :=  R[n0][][a];

T[n0__][j__][a_] := ( T[n0][Sequence @@ Sort[{j}]][a] ) /; !OrderedQ[{j}];


(* special B-stuff*)
tT[2][1][{p10_,m02_,m12_}]   := B1[p10,  m02, m12];
(* XXX *)
tT[2][1,1][{p10_,m02_,m12_}] := B11[p10, m02, m12] /; 
                  ($LimitTo4===True) || (p10 === 0);
tT[2][0,0,1][{p10_,0,0}]   := 
  If[$LimitTo4 === True,
      p10/36 + (p10*B0[p10, 0, 0])/24, -(p10*B0[p10, 0, 0])/(8*(1 - $dIM))
    ];
tT[2][1,1,1][{p10_,0,0}]   := 
  If[$LimitTo4 === True,
     -1/12 - B0[p10, 0, 0]/4, ((2 + $dIM)*B0[p10, 0, 0])/(8*(1 - $dIM))
    ];

(* (4.18), *)
tT[N_Integer][k_Integer,i___Integer][a_List] := 
(*tT[N][k,i][a] =*) Block[ {P, M ,r, kp },     
  P = 1 + Length[{i}]; M = getm[a];
   Sum[ Xinv[M][k, kp][a]  ( R[N,kp][i][a] - 
       Sum[ delt[ kp,{i}[[r]] ] * (T[N]@@Join[{0,0}, Delete[{i},r]])[a],
             {r, P-1} ]     ), {kp, M} ]              ];

(* no M's in i *)
R[N_Integer, 0, 0][i___Integer][a_List] := Block[{q,P,M},
      q = Length[{i}]; P = 2 + q; M = getm[a];
      demon[a[[-N]]] T[N][i][a]  + T[N-1][i][ c[0][a] ]
                                                 ] /; FreeQ[{i}, getm[a]];

R[N_Integer,0,0][i___Integer, mm:(_Integer)..][a_List] := Block[
     {q,M,P,j,k},
      q = Length[{i}]; M = getm[a]; P = Length[{mm}] + 2 + q;
      demon[ a[[-N]] ] T[N][i, mm][a] + 
(* here was the tough bug found by Ralph Schuster ... *)
      (-1)^(P - q) ( T[N - 1][i][c[0][a]] + Sum[
      Binomial[P - 2 - q, j] * Sum @@ Prepend[ Array[List[k[#], M - 1]&, j],
       (T[N-1]@@Join[{i}, Array[k,j]])[c[0][a]]
                                             ], {j,P - 2 - q}
                                                  ] )          ] /; 
                                          ({mm}[[1]] === getm[a]);

(* ***************************************************************** *)
(*                          pave26                                   *)
(* ***************************************************************** *)

(* 4.19 , no M's*)
R[N_Integer, k_Integer][i___Integer][a_List] := Block[{q,P,M},
     q = Length[{i}]; P = 1 + q; M = getm[a];
     1/2( (T[N - 1] @@ til[i][k])[ c[k][a] ] theta[k, i] -
     f[k][a] T[N][i][a]  -  T[N - 1][i][c[0][a]] 
        )                                             ] /;
       FreeQ[{i}, getm[a]];

R[N_Integer,k_Integer][i___Integer, mm:(_Integer)..][a_List]:=Block[
      {q, P, M, kk, j},
      q = Length[{i}]; P = Length[{mm}] + 1 + q; M = getm[a];
       1/2( (T[N-1] @@ til[i,mm][k])[c[k][a]] theta[k, i, mm] -
       f[k][a] T[N][i,mm][a] -(-1)^(P - 1 - q) (
          T[N - 1][i][c[0][a]] + 
          Sum[ 
          Binomial[P - 1 - q, j]  Sum @@ Prepend[Array[List[kk[#], M -1 ]&,j],
                                 (T[N - 1]@@Join[{i}, Array[kk, j]])[c[0][a]]
                                                ], {j, P - 1 - q}
          ]))                                                        ] /;
                                          ({mm}[[1]] === getm[a]);

(* 4.20 *)
(* thetadef *)
theta[k_Integer, i___Integer] := 1 /; FreeQ[{i}, k];
theta[k_Integer, i___Integer] := 0 /;!FreeQ[{i}, k];

tm[a_Integer, b_Integer]:= a /; a<=b;
tm[a_Integer, b_Integer]:= (a-1) /; a > b;
til[][_]={};
til[x__][k_]:= Map[tm[#,k]&, {x}]; 

(* ***************************************************************** *)
(*                          pave27                                   *)
(* ***************************************************************** *)
(* PaVeReducedef *)
(* ********************************************************************** *)
(* Decomposition down to scalar integrals *)
(* ********************************************************************** *)


cancel[x_]:=Cancel[x/.Plus->pll]/.pll->Plus;
PaVeReduce[x_, y___Rule]:= Block[{op, wriout, nnx = x},
        op = Join[{y}, Options[PaVeReduce]];
        wriout = WriteOutPaVe /. op;
        If[!FreeQ[nnx, StandardMatrixElement], 
           nnx = Expand2[nnx, StandardMatrixElement];
          ];
        If[StringQ[wriout] && (Head[x] === PaVe),
           nnx  = pavitp @@ Join[{nnx, wriout}, op],
           nnx = pavereduce[nnx, y]
          ];
              nnx];

pavitp[xXX_PaVe, dir_,opts___] := Block[{nx, file, temp, set,xxx},
   paV[xy__, p_List, m_List] := PaVe[xy,C,p,C,m];
   xxx = paV@@xXX;
   nx = StringReplace[ ToString[InputForm[xxx], PageWidth -> 222],
                       {", "->"","^"->"","{"->"", "/" -> "",
                       "}"->"", "["->"", "]"->"", "*" -> "", " " -> "" }
                     ];
                      nx = StringJoin[dir, nx, ".s"];
    If[Streams[nx] === {},
                      file = FileNames @@ {nx};
                      If[file =!= {},
                         print2P["file  =", file];
                         temp =( Get @@ {nx} ) // PaVeOrder;
 (* If something went wrong in writing the file *)
                         If[ Head[temp]=!=Plus, file = {} ]
                        ];
                      If[file ==={} ,
                          temp = PaVeReduce[xXX, WriteOutPaVe->False,opts
                                           ]//PaVeOrder;
                         print2P["writing result to ",nx];
                         OpenWrite @@ {nx, FormatType -> InputForm };
                         WriteString @@ {nx, "( "};
                         Write @@ {nx, temp};
                         WriteString @@ {nx, "  ) "};
                         Close @@ {nx}
                        ],
        temp = PaVeReduce[xXX, WriteOutPaVe->False,opts]//PaVeOrder;
      ];
                           temp]; 

pavereduce[0,___]:=0;
pavereduce[w_,___]:=w/;nterms2[w]===1 && FreeQ[w,PaVe];
pavereduce[ a_ b_,ops___ ]:=cancel[ a pavereduce[ b,ops] ]/;
                           FreeQ[a,PaVe]&&!FreeQ[a,StandardMatrixElement];
pavereduce[ w_(*Plus*),ops___ ]:=
     Block[{mpa,nw,nn,pre,re=0,nulll,op = Flatten[{ops}]},
            mpa = w/.StandardMatrixElement[__]->0;
            nw = w - mpa;
            re = pavereduce[mpa, ops];
            If[Head[nw] === Plus, nn = Length[nw], nn = 1];
            For[ij=1,ij<=nn,ij++,
                print2P["breaking down # ",ij," / ",nn];
                If[ nn===1, pre = PartitHead[ nw,StandardMatrixElement ],
                            pre = PartitHead[ nw[[ij]],StandardMatrixElement ]
                  ];
                re = re + pre[[2]] pavereduce[ pre[[1]],ops ]
               ];
        re]/;!FreeQ[w, StandardMatrixElement];

pavereduce[pvli_List,op___]:=Block[{i,set,le=Length[pvli],npvli},
                                   npvli = {};
                            Do[ print2P[" Working with # ",i," out of ",le];
                                npvli=Append[ npvli,
                                       pavereduce[pvli[[i]],op] ],
                                {i,le}
                              ];
                            npvli];
 
(* ********************************************************************** *)
 
Options[ PaVeReduce ] = { Dimension -> True,
                          IsolateHead->False,
                          Mandelstam->{},
                          PaVeOrderList -> {},
                          WriteOutPaVe -> ""
                       };

(* This default setting of Dimension results --- together with 
   $LimitTo4 = True --- into dimensional regularization 
   (for the ultraviolett divergencies ) with 
   the limit Dimension -> 4 being taken.
   If the option Dimension is set to some explicit variable (d for instance),
   no limit is taken and d occurs in the result.
*)
(* ***************************************************************** *)
(* ***************************************************************** *)
(*                          pave28                                   *)
(* ***************************************************************** *)

pavereduce[brex_,optis___]:=Block[{sq,t,tt,ma,rest,lin,
                          result,var,ij,lra,nra,des,hed,mand3,
                          isolating,kkk, ktri,trick,dimen,
                          tog,su,msu,isok,isocc,is,pvs,pl2,paveorderli,
                          il,pail, breakx,  rlin,plup,plusu,iit,
                          colstu,cofun2,ir,nresult,cofun,newt, tvarS
                        },
  { paveorderli , mand,isok, dimen}= { PaVeOrderList, Mandelstam,
                                       IsolateHead, Dimension
(*, FinalSubstitutions*)}/.
     Join[ {optis},Options[PaVeReduce] ];

(* a little bit fishy, since this yields a side effect, but it is only 
   used once in tT *)
(*
If[(dimen =!= True) && ($LimitTo4 =!= True), $dIM = dimen];
*)
If[(dimen =!= True) && ($LimitTo4 =!= True), $dIM = dimen, $dIM = D];

If[!FreeQ[ brex, PaVe], breakx = collect2P[brex, PaVe], breakx=brex ];
isolateP[x__]:=Isolate[x, IsolateHead -> isok];

tri[xx_  yy_]:= tri[xx] tri[yy];
tri[any_ xx_]:= ( any tri[xx] )/;FreeQ[any,Plus] || Head[any]===PaVe;
tri[any_ ]:=any /;FreeQ[any,Plus] || Head[any]===PaVe;
 
mand = mand/.Negligible->Identity;
If[ mand==={}, 
    If[($LimitTo4 === False ) && (Head[brex] === PaVe),
       tvarS = Variables[ Join @@ Take[brex, -2] ];
       trick[z_] := Collect2[z, tvarS],
       trick[z_]:=z
      ],
    trick[z_]:=trick[z]=TrickMandelstam[z,mand]//Factor2
  ];
 
msu = {};
   
pl2[x__]:=kkk[ Plus[x] ]/; FreeQ2[ {x},{A0,B0,B1,B00,B11,C0,D0,PaVe} ];
backpc[a_,b_,c_,d_,e_,f_]:=PaVeOrder[C0[a,b,c,d,e,f],
                                     PaVeOrderList -> paveorderli];
backpd[a_] := D0[a];
backpd[a_,b_,c_,d_,e_,f_,m1_,m2_,m3_,m4_]:=
 PaVeOrder[ D0[a,b,c,d,e,f,m1,m2,m3,m4], PaVeOrderList -> paveorderli];
print2P["starting pavebr ing"];
pluep[yy__]:=Plus[yy]/;!FreeQ2[{yy}, {$epsilon,A0,B0,B1,B00,B11}];
tim = Timing[
t =  breakdown[ (breakx/.msu) ];
            ];
t = t/.C0->backpc/.D0->backpd;
  
If[ !FreeQ[t, HoldForm], t = FixedPoint[ReleaseHold,t] ];
t = collect3[ t, {A0,B0,B1,B00,B11,C0,D0,PaVe,Dot}];
If[ !FreeQ[t, $epsilon], 
     t = Expand[t/.Plus->pluep]/.$epsilon->0/.pluep->Plus
  ];

result = t;
 
(* ***************************************************************** *)
(*                          pave29                                   *)
(* ***************************************************************** *)


(* get the "linear" part *)
print2P["check4 ", MemoryInUse[]," MB used"];
 cofun[0]=0;
 cofun[a_ b_]:=a cofun[b]/;!FreeQ2[a, {A0,B0,B1,B00,B11,C0,D0,PaVe} ];
 cofun2[0]=0;
 cofun2[y_]:=y/;FreeQ[y,Plus];
 cofun2[yy_]:= trick[yy];
  
print2P["check7"];
 
If[ isok=!=False && Head[result]=!=PaVe,
    result = cofun/@( result + nuLL );
    isolatefirst[a_ b_]:=(a isolatefirst[b])/;
                           FreeQ2[a,{A0,B0,B1,B00,B11,C0,D0,PaVe}];
    isolatefirst[a_]:=a /; FreeQ2[a,{A0,B0,B1,B00,B11,C0,D0,PaVe}];
    result = isolatefirst /@ result;
    isolatetri[a_]:=isolateP[ trick[a] ];
    result = (result/.nuLL->0)/.isolatefirst->isolatetri/.cofun->cofun2
  ];
   
If[ isok=!=False, result = isolateP[ result ], result = trick /@ result ];
result]/;FreeQ[brex,StandardMatrixElement];
(* **************************************************************** *)
print2P[x__]:=Feyn`Calc`Main`print2[x];
WriteString["stdout", "."];
End[];
(* ***********************************************************************  *)
(*                             SUN.m                                        *)
(* ***********************************************************************  *)
(* similar to GellMann.m, but use T_a matrices instead of GellMann matrices *)
(* Remember: T_a = 1/2 lambda_a *)
(* ***********************************************************************  *)
(* SU(N) algebra *)
(* ***********************************************************************  *)
 Begin["Feyn`Calc`SUN`"];
(* ***********************************************************************  *)
(* ******************************** SUNF is totally antisymmetric *********)
(* SUNF       structure constant f_{abc}
   SUNDelta   SU(N) delta function
   SUNT       Generators of SU(N)
   SUNN         "N" = Number of colors
*)

SetAttributes[SUNIndex, Flat];
(* change the default ... *)
Options[SUNF]={SUNFToTraces -> False};
SUNF[a___,x_,b___] := SUNF[a, SUNIndex[x], b] /; FreeQ2[x, 
                               {SUNIndex, Rule, Pattern}];

SUNF /: SUNF[a___, SUNIndex[ComplexIndex[b_]], c___] *
        SUNF[d___, SUNIndex[ComplexIndex[b_]], e___] := 
        Block[{dum = Unique[System`C]},
              SUNF[a, SUNIndex[dum], c] SUNF[d, SUNIndex[dum], e]
             ];

(*
SUNF /: m_. SUNF[a___, SUNIndex[ComplexIndex[b_]], c___
                ]:= m SUNF[a, SUNIndex[b], c] /; FreeQ[m, ComplexIndex[b]];
*)

Literal[SUNF[a___,x_,b___,x_,c___]]:=0;
Literal[SUNF[a___,x_,y_,b___]]:=-SUNF[a,y,x,b]/;!OrderedQ[{x,y}]/;Length[{a}]<2;
(* ******************************** SUNDelta is symmetric *************** *)
 SetAttributes[ SUNDelta, Orderless ];
(********************************** For SU(3) **************************** *)
 SUNDelta[a_ /; FreeQ[a, SUNIndex] , b_ /; FreeQ[b, SUNIndex]] := 
      SUNDelta[SUNIndex[a], SUNIndex[b]];
 Literal[SUNDelta[i_SUNIndex, i_SUNIndex]] := SUNN^2 - 1;
SUNDelta /: Literal[SUNDelta[i_SUNIndex, j_SUNIndex]^2] := (SUNN^2 - 1) /; 
    i =!= j;
(************************** Define the contraction property for the SUNDelta *)
 SUNDelta/: SUNDelta[i_SUNIndex, j_SUNIndex] y_[z__] := 
              ( y[z] /. i -> j )/; !FreeQ[y[z]//Hold, i] &&
                 FreeQ[y[z], SUNDelta[__]^n_Integer?Negative];

(* change SUNT[a. b] to SUNT[a].SUNT[b] (etc) *)
 SUNT[Literal[Dot[x__]]]:=Dot@@( SUNT/@{x} );
 SUNT[x_,y__]:=SUNT[x.y];
 SUNT[b_ /; (Head[b] =!= SUNIndex)] := SUNT[SUNIndex[b]];

 SUNT/: Literal[SUNT[x_SUNIndex ]. SUNT[x_SUNIndex ]] := 1/2 (SUNN^2-1)/SUNN;

(* change SUNT' which are multiplied with each other to lambdaT's *)
 lambdaT[1]=1;
 gm2lambdaT[x__]:= (gmlin@@( {x}/.SUNT->lambdaT ) )/.gmlin->Dot;
(************************* linearity ************************************* *)
 gmlin/: Literal[gmlin[gmlin[x__]]] := gmlin[x];
 gmlin[ a___, b_ c_, d___ ] := b gmlin[a,c,d]/;FreeQ[b, lambdaT] && 
                                 Feyn`Calc`Main`noncommQ[b];
 gmlin[ a___, b_ , d___ ]   := b gmlin[a,d]/;FreeQ[b, lambdaT] &&
                                 Feyn`Calc`Main`noncommQ[b];
 gmlin[]=1;

gellm1[x_, y__] := gellm1[x. y];
gellm2[x_, y__] := gellm2[x. y];
(************************* cyclicity ************************************* *)
 gmcyc[x__]:=gellm1 @@ First[ NestList[RotateLeft,{x},Length[{x}]-1]//Sort ];
(******************** define the properties of trace of Gell-Mann matrices *)
 gellm2[ ] = gmcyc[ ] = SUNN;         (* unit trace  *)
(************************ each single Gell-Mann matrix has vanishing trace *)
 gellm2[ lambdaT[_] ] := 0;
(************************ Cvitanovic - rules ***************************** *)
 gellm2[Literal[Dot[a___, lambdaT[i_SUNIndex], b___, 
                          lambdaT[i_SUNIndex], c___]]]:=
        1/2 gmcyc[b] gmcyc[a, c] - 1/2/SUNN gmcyc[a, b, c]; 

 gellm2/: gellm2[Literal[Dot[a___, lambdaT[i_SUNIndex], b___]]]^2 :=
          1/2 gmcyc[a, b, a, b] - 1/2/SUNN gmcyc[a, b]^2;

 gellm2/: gellm2[Literal[Dot[a___, lambdaT[i_SUNIndex], b___]]] * 
	  gellm2[Literal[Dot[c___, lambdaT[i_SUNIndex], d___]]]:=
          1/2 gmcyc[a, d, c, b] - 1/2/SUNN gmcyc[a, b] gmcyc[c, d];

 Literal[SUNF[i_,j_,k_,SUNFToTraces -> False]]:= SUNF[i,j,k];
 Literal[SUNF[i_,j_,k_,op_:{}]]:= 2 I (SUNTrace[ SUNT[i,k,j] ] - 
			      SUNTrace[ SUNT[i,j,k] ]
			     )/;
                         (SUNFToTraces/.Flatten[Join[{op},Options[SUNF]]]);
   
 f2tr[i_,j_,k_,___]:= 2 I (gmcyc @@ lambdaT/@{i,k,j} - gmcyc @@ lambdaT/@{i,j,k});
(*********** do the application of the Cvitanovic - relations step by step *)
 cvit[x_Plus]:=cvit/@x;
 cvit[x_]:= cvit[x]=ExpandAll[ x/.gellm1->gellm2 ];
(*********** this is the function which puts everything together ********* *)
 SUNTrace[n_]:=  SUNN n /; FreeQ[n,SUNT] && FreeQ[n, Pattern];
 SUNTrace[n_ y_]:=n SUNTrace[y]/;FreeQ[n,SUNT];
 Literal[SUNTrace[Dot[a:(SUNT[SUNIndex[ComplexIndex[_]]])..]]] := 
   SUNTrace[ Dot @@ Reverse[{a}/.ComplexIndex->Identity] ]; 
 Literal[SUNTrace[1 ..]]= SUNN;
 Literal[SUNTrace[SUNT[x_SUNIndex ]]] := 0;
 Literal[SUNTrace[SUNT[x_SUNIndex] . SUNT[y_SUNIndex]]] := 1/2 SUNDelta[x, y];
 Literal[SUNTrace[ x_, y__ ] ]:= SUNTrace[x.y];
 Literal[SUNTrace[n_?NumberQ y_] ] := n SUNTrace[y];
 SUNTrace[ SUNTrace[x_] y_. ]  := SUNTrace[x] SUNTrace[y];
 SUNTrace/:  Literal[SUNTrace[(A___).  SUNT[x_SUNIndex].B___]] *
		  Literal[SUNTrace[(a___). SUNT[x_SUNIndex].b___]] := 
              FixedPoint[cvit,  (gmcyc[A.SUNT[x].B] * 
		                 gmcyc[a.SUNT[x].b])/.
				 SUNTrace->gellm1/.
				 Dot->gm2lambdaT/.gellm1->gellm2/.SUNF->f2tr 
	 	        ]/.lambdaT->SUNT/.gellm2->SUNTrace;
 SUNTrace /: Literal[SUNTrace[x_]^2] := SUNTrace[x] * SUNTrace[x];

 SUNTrace[ expr_ ]:= fixgell[expr]/;NumberQ[fixgell[expr]] && FreeQ[expr, Pattern];
 SUNTrace[ expr_ ]:= (fixgell[expr (*/.SUNT[1]->1*) ]/.
                           gellm2->SUNTrace)/;
			   (expr=!=(fixgell[expr]/.gellm2->Identity)) && 
                           FreeQ[expr, Pattern];
 gellm1[x_Plus]:=gellm1 /@ x;
 gellm1/: gellm1[x_ y_] := x gellm1[y]/;FreeQ[x,lambdaT];
 gellm1/: gellm1[x_Dot gellm1[y___]]:=gellm1[y] gellm1[x];
 gellex[z_]:=gellm1[ExpandAll[z]];
 fixgell[x_]:=fixgell[x]=
 FixedPoint[cvit,       ( gellm1[ ExpandAll[x/.SUNTrace->gellex/.
				         Dot->gm2lambdaT/.SUNF->f2tr] ]
	                )/.gellm1->gellm2, 9
           ]/.lambdaT->SUNT;
(* *********************************************************************** *)

(*Renamed*)
SetAttributes[Rename, Listable];

Rename[exp_] := Block[{new=exp,old,subst,uh,ne,uhc=0, uuh, suI,
                       suh= {}, sub={}, sam, dummy = Unique[System`C]}, 
 sam[iii_, uuh_ ] := (AppendTo[sub, SUNIndex[uuh] -> SUNIndex[iii]]
                     ) /; FreeQ[sub, SUNIndex[uuh]];
uuh[] := SUNIndex[dummy[uhc++]];
suii[0] = 0;
suii[a_Plus] := suii /@ a;
suii[y_] := Block[{ste, ste2, dumy},
                  ste = Select[y dumy, !FreeQ2[#, {SUNIndex}]&];
                  ste2 = Select[y dumy, FreeQ2[#, {SUNIndex}]&];
                  (ste2 suI[ste])/.dumy->1](* /. suI -> Identity*);
subst = {suI[SUNT[SUNIndex[ii_]] ff_] :>
           (sam[ii, uh[]]; (SUNT[SUNIndex[ii]] (ff /. suI -> Identity)) /. 
           SUNIndex[ii]-> uh[])/;
           (!FreeQ[ff,SUNIndex[ii]]) && FreeQ[ii,dummy[_]],

        Dot[ A___,SUNT[SUNIndex[ii_]], B___, 
                  SUNT[SUNIndex[ii_]], Z___]      :> 
           (sam[ii, uh[]];
            Dot[ A,SUNT[SUNIndex[ii]],B,SUNT[SUNIndex[ii]],Z ] /.
               SUNIndex[ii] -> uh[] 
           ) /; FreeQ[ii, dummy[_]] , (* 
        suI[SUNTrace[Dot[ A___,SUNT[SUNIndex[jj_]], B___]] ff_]   :>
           (sam[jj, uh[]]; (sTr[Dot[ A,SUNT[SUNIndex[jj]],B ]] ff
                           ) /. SUNIndex[jj] -> uh[] /. sTr -> SUNTrace
           )/; !FreeQ2[ff,{SUNIndex[jj]}] && FreeQ[jj,dummy[_]], *)
        suI[Dot[ A___,SUNT[SUNIndex[jj_]], B___] ff_]   :>
           (sam[jj, uh[]]; (Dot[ A,SUNT[SUNIndex[jj]],B ] ff
                           ) /. SUNIndex[jj] -> uh[] 
           )/; !FreeQ[ff,SUNIndex[jj]] && FreeQ[jj,dummy[_]],
        suI[Literal[SUNF[A___,SUNIndex[ij_],B___]] *
            Literal[SUNF[V___,SUNIndex[ij_],W___]] ff_] :> 
           (sam[ij, uh[]]; (SUNF[A, SUNIndex[ij], B] *
                            SUNF[V, SUNIndex[ij], W] ff
                           )/. SUNIndex[ij] -> uh[]
           )/; FreeQ[ij, dummy[_]],
        suI[Literal[SUNF[A___,SUNIndex[ij_],B___]] ff_] :> 
           (sam[ij, uh[]];(SUNF[A,SUNIndex[ij],B] ff
                          )/. SUNIndex[ij] -> uh[]
           )/; !FreeQ[ff, SUNIndex[ij]] && FreeQ[ij, dummy[_]]
        };
(* CHANGE 28.6. 93 *)
suff[x_] := FixedPoint[(uh[] = uuh[]; (suii[#]//.subst)/.suI->Identity)&, x, 42];
(*
suff[x_] := FixedPoint[(uh[] = uuh[]; (suii[#]/.subst)/.suI->Identity)&, x, 42];
*)
If[$VeryVerbose>1,Print["expanding w.r.t. SUNF done "]];
new = Expand2[new, SUNIndex];
If[$VeryVerbose>1,Print["expanding w.r.t. SUNF done "]];
new = suff[new];
new = backsubfun[new, sub];
new];

backsubfun[xxx_, {}]:=xxx;
backsubfun[xxx_, {a_ -> b_, c___Rule}] := 
 If[FreeQ[xxx, b], backsubfun[xxx /. a -> b,{c}], backsubfun[xxx, {c}]];
(* *********************************************************************** *)

Literal[sunTRACEcyc[Dot[z:SUNT[_]..]]] :=
 sunTRACE[Dot@@RotateLeft[{z}, Position[{z},Last[Sort[{z}]]][[1,1]]]];

(*SUNSimplifydef*)
SetAttributes[SUNSimplify, Listable];
Options[SUNSimplify] = {
                        Expanding    -> False,
                        Factoring    -> False,
                        SUNFJacobi   -> False,
                        SUNNToCACF   -> True, 
                        SUNFToTraces -> True};
SUNSimplify[x_, opts___Rule] := Block[{af, temp = x, sft, sunf, expan,sunsi,
                                       jac, expanding, factoring,ntemp,dotT},
If[(!FreeQ[x, SUNIndex]) || (!FreeQ[x, SUNN]),

sunsi = {Literal[Dot[xx___, SUNT[a_SUNIndex], SUNT[b_SUNIndex],
                  SUNT[a_SUNIndex], yy___]] :>
         ((-1)/(2 SUNN) dotT[xx, SUNT[b], yy]),
         Literal[SUNF[xx_SUNIndex, yy_SUNIndex, zz_SUNIndex] * 
                 SUNT[zz_SUNIndex] . SUNT[yy_SUNIndex]] :>
          I/2 SUNN SUNT[xx],
         Literal[SUNF[xx_SUNIndex, yy_SUNIndex, zz_SUNIndex] * 
                 SUNT[yy_SUNIndex] . SUNT[zz_SUNIndex]] :>
         -I/2 SUNN SUNT[xx],
          Literal[SUNF[xx_SUNIndex, yy_SUNIndex, zz_SUNIndex] *
                 SUNT[xx_SUNIndex] ] :>
         I SUNT[zz] . SUNT[yy] - I SUNT[yy] . SUNT[zz],
          Literal[SUNF[xx_SUNIndex, yy_SUNIndex, zz_SUNIndex] *
                 SUNT[zz_SUNIndex] ] :>
         I SUNT[yy] . SUNT[xx] - I SUNT[xx] . SUNT[yy]
         }  /. dotT[]->1 /. dotT->Dot;
expanding = Expanding /. {opts} /. Options[SUNSimplify];
factoring = Factoring/. {opts} /. Options[SUNSimplify];
af = SUNNToCACF /. {opts} /. Options[SUNSimplify];
sft = SUNFToTraces/. {opts} /. Options[SUNSimplify];
jac = SUNFJacobi /. {opts} /. Options[SUNSimplify];
If[sft === True, sunf[a__] := SUNF[a, SUNFToTraces -> True], 
                 sunf[a__] := SUNF[a, SUNFToTraces -> False]
  ];

If[!FreeQ[temp, Dot], temp = DiracSimplify[temp, Expanding -> False]];
temp = Rename[temp /. sunsi];
If[!FreeQ[temp, DiracTrace],
   temp = DiracSimplify[temp, Expanding -> False] /.
          {Literal[DiracTrace[Dot[xx__SUNT] dd_, dops___]] :>
                 SUNTrace[Dot[xx]] DiracTrace[dd, dops],
           Literal[DiracTrace[SUNT[xx_]  dd_, ___]] :> 0 /; FreeQ[dd,SUNIndex],
           Literal[DiracTrace[SUNT[xx_] . dd_, ___]] :> 0 /; FreeQ[dd,SUNIndex]
          } 
  ];

If[FreeQ2[temp, {SUNTrace}] && sft === False,
   expan = Identity, expan = Expand(*All*)[#, SUNIndex]&
  ];

temp = FixedPoint[expan, temp /. SUNTrace -> sunTRACEcyc /.
                        Dot -> gm2lambdaT /.lambdaT -> SUNT /.
                         sunTRACE -> SUNTrace /. SUNF -> sunf
             ];
If[jac === True && !FreeQ[temp, SUNF],
   temp = temp /. { SUNF[SUNIndex[a_], SUNIndex[c_], SUNIndex[e_]]*
                    SUNF[SUNIndex[b_], SUNIndex[d_], SUNIndex[e_]] :>
                    (SUNF[SUNIndex[a], SUNIndex[b], SUNIndex[e]]* 
                     SUNF[SUNIndex[c], SUNIndex[d], SUNIndex[e]] +
                      SUNF[SUNIndex[b], SUNIndex[c], SUNIndex[e]]*
                      SUNF[SUNIndex[a], SUNIndex[d], SUNIndex[e]]
                     ) /; Sort[{ {a,c,e, b,d,e}, {a,b,e, c,d,e}, {b,c,e, a,d,e} 
                               }][[1]] === {a,c,e, b,d,e}
                  };
  ];

If[!FreeQ[temp, ComplexIndex],
   temp = temp /. Dot -> dot /. 
  {
   (SUNT[SUNIndex[a_]] /; FreeQ[{a}, ComplexIndex]) * 
    SUNT[SUNIndex[ComplexIndex[b_]]] :>
   SUNTrace[SUNT[SUNIndex[a]] . SUNT[SUNIndex[b]] ], 

   (SUNT[SUNIndex[ComplexIndex[a_]]] ) * 
   dot[(b:SUNT[SUNIndex[_]]..) /; FreeQ[{b}, ComplexIndex]]  :> 
   SUNTrace[SUNT[SUNIndex[a]] . b ], 

   (SUNT[SUNIndex[a_ /; FreeQ[{a}, ComplexIndex]] ]) * 
   dot[b:SUNT[SUNIndex[ComplexIndex[_]]]..]  :> 
   SUNTrace[SUNT[SUNIndex[a]] . b ], 

   dot[(a:SUNT[SUNIndex[_]]..) /; FreeQ[{a}, ComplexIndex]] * 
   dot[b:SUNT[SUNIndex[ComplexIndex[_]]]..] :>
   SUNTrace[Dot @@ Join[{a}, Reverse[{b}]/.ComplexIndex -> Identity]] 
  } /. dot -> Dot;
 ];

If[!FreeQ[temp, SUNIndex], temp = temp /. sunsi];

If[expanding === True, temp = Expand[temp],
   If[LeafCount[temp] < 242 && Head[temp] === Plus, 
      ntemp = Expand[temp];
      If[LeafCount[ntemp] < LeafCount[temp], temp = ntemp]
     ];
  ];
      
If[factoring === True, temp = Factor2[temp]];


If[ af === True, temp = temp /. (1-SUNN^2) -> (-CF 2 CA) /. 
                              SUNN -> CA /. (-1 + CA^2)->(2 CA CF)];
 ](*thatsthemainIf*);
temp];
(* *********************************************************************** *)
WriteString["stdout", ".|\n"];
End[];
EndPackage[];
(* *********************************************************************** *)
(* end of feyncalc XXX END of FeynCalc XXX End of ... *)
(* *********************************************************************** *)
