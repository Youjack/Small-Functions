(* ::Package:: *)

$LoadAddOns = {"FeynArts"};
BeginPackage["YoujackOneLoop`",{"FeynCalc`"}];
Needs["FeynArts`"];
$FAVerbose = 0;

SetScatOptions::usage =
  "Set options for calculating scattering amplitudes:
  amplitudes are M without an \[ImaginaryI];
  external lines are not truncated.";

Set1PIOptions::usage =
  "Set options for calculating 1PIs:
  amplitudes are M without an \[ImaginaryI];
  external lines are truncated;
  reducible diagrams are excluded.";

AddMassRegulator::usage = 
  "AddMassRegulator[amp, m] add mass m to every massless propagator in amp.";

Z2d::usage =
  "Z2d[expr, {{Z1,d1},{Z2,d2},...}] replaces Zi in expr with di,
  and expand di to first order.";

DRCoupling::usage =
  "DRCoupling[expr, g, dimg] multiplies g by ((4\[Pi]\[ExponentialE]^-\[Gamma])^(-1/2)\[Mu])^dimg.";
DRCoupling\[Prime]::usage =
  "DRCoupling\[Prime][expr, g, dimg] multiplies g by \[Mu]^dimg.";

DRExpand::usage =
  "DRExpand[expr, dim] replaces D\[RightArrow]dim-2\[Epsilon] in expr,
  expands expr with respect to \[Epsilon],
  and returns {\[Epsilon] dependent part, \[Epsilon] independent part}.";

DRExpandScaleless::usage =
  "DRExpandScaleless[expr, dim, zero\[CapitalDelta]] replaces D \[RightArrow] dim-2\[Epsilon] in expr,
  then replaces zero\[CapitalDelta]^-Epsilon \[RightArrow] Epsilon(1/EpsilonUV - 1/EpsilonIR).
  The return format is the same as DRExpand.";

FCFPOneLoop::usage =
  "FCFPOneLoop[int, l, x] applies `FCFeynmanParametrize` to compute one-loop tensor integral.
  (the final result have an implicit Feynman-parameters integral with Gamma[b])";

FPOneLoopDenom::usage =
  "FPOneLoopDenom[denom, l, x, \[CapitalDelta]] Feynman-parametrizes denom,
  and return {l\[RightArrow]shifted l, \[CapitalDelta]\[RightArrow]..., power of denom}.";

ReduceOneLoopNumr::usage =
  "ReduceOneLoopNumr[numr, lRule] shifts loop momentum in numr according to lRule,
  reduces the loop-momentum tensors to scalars,
  and returns {{a1 in (l^2)^a1, prefactor1},...}.";

StdOneLoop::usage =
  "StdOneLoop[b, reducedNumr, \[CapitalDelta]] gives the standard one-loop integral.
  (the final result have an implicit Feynman-parameters integral with Gamma[b])";

FPOneLoop::usage =
  "FPOneLoop[ampLoop, l, x, \[CapitalDelta]]
  (the final result have an implicit Feynman-parameters integral with Gamma[b])";

FPIntegrate::usage =
  "FPIntegrate[b, expr, x] integrates over b Feynman parameters x[i] in expr and multiplies it by Gamma[n].";

Begin["`Private`"];

$lorentzIndices = {Global`\[Mu], Global`\[Nu], Global`\[Rho], Global`\[Sigma]};

(* set options *)
(
  SetOptions[Paint, Numbering -> Simple, SheetHeader -> None, ColumnsXRows -> {4, 1}];
  SetOptions[Render, ImageSize -> {4 128, 128}];

  SetOptions[CreateFeynAmp, PreFactor -> -I, GaugeRules -> {_FAGaugeXi -> 1}];

  SetOptions[FCFAConvert, List -> False, SMP -> False,
    LorentzIndexNames -> $lorentzIndices,
    ChangeDimension -> D,
    UndoChiralSplittings -> True, Contract -> True,
    DropSumOver -> True];
);
SetScatOptions[] := (
  SetOptions[CreateTopologies,
    ExcludeTopologies -> {Tadpoles, WFCorrections, WFCorrectionCTs}];
  SetOptions[CreateFeynAmp, Truncated -> False];
);
Set1PIOptions[] := (
  SetOptions[CreateTopologies,
    ExcludeTopologies -> {Tadpoles, WFCorrections, WFCorrectionCTs, Reducible}];
  SetOptions[CreateFeynAmp, Truncated -> True];

  $KeepLogDivergentScalelessIntegrals = True;
  SetOptions[FCFeynmanParametrize, FeynmanIntegralPrefactor -> "Textbook"];
);

(* functions for lopp calculation *)

AddMassRegulator[amp_, m_] :=
  ReplaceAll[amp, PropagatorDenominator[p_, 0] :> PropagatorDenominator[p, m]];

Z2d[expr_, ZList_List] := Module[{\[Alpha]},
  ReplaceAll[expr, (#[[1]] -> 1 + \[Alpha] #[[2]]) & /@ ZList] //
    Normal@Series[#, {\[Alpha], 0, 1}] & //
    ReplaceAll[\[Alpha] -> 1] //
    Collect[#, (#[[2]]) & /@ ZList, Simplify] &
];

DRCoupling[expr_, g_, dimg_] :=
  ReplaceAll[expr, g -> g ((4 Pi E^-EulerGamma)^(-1/2) ScaleMu)^dimg];
DRCoupling\[Prime][expr_, g_, dimg_] :=
  ReplaceAll[expr, g -> g ScaleMu^dimg];

DRExpand[expr_, dim_] :=
  FCReplaceD[expr, D -> dim - 2 Epsilon] //
  Normal@Series[#, {Epsilon,0,0}] & //
  {SelectNotFree2[#, Epsilon], SelectFree2[#, Epsilon]} & // Simplify;

DRExpandScaleless[expr_, dim_, zero\[CapitalDelta]_] :=
  expr // Collect[#, zero\[CapitalDelta], Simplify] & //
  If[Head@# === Plus, List@@#, List@#] & // Map[
    Switch[
      SelectNotFree[#, zero\[CapitalDelta]] /. {zero\[CapitalDelta]^n_ :> n, 1 -> 0} /.
        D -> 4 - 2 Epsilon // Simplify,
      0, {0,#},
      -Epsilon,
        SelectFree[#, zero\[CapitalDelta]] //
        FCReplaceD[#, D -> dim - 2 Epsilon] & // Simplify //
        Normal@Series[Epsilon #, {Epsilon,0,0}] & //
        { # / Epsilon, - # / EpsilonIR } &,
      _, {0,0}
    ]&
  ] // Apply[Plus];

FCFPOneLoop[int_, l_, x_] := Module[
  {tensorIntList, FPList},
  (* decompose `int` into {{int1,factor1},{int2,factor2},...} *)
  tensorIntList = int //
    DiracSimplify // Simplify //
    Expand2[#, l] & //
    If[Head@# === Plus, List@@#, List@#] & //
    Uncontract[#, l, Pair -> {l}] & //
    ({SelectNotFree[#, l], SelectFree[#, l]} &) /@ # &;
  (* apply `FCFeynmanParametrize` to `inti` and divide it by Gamma[n] *)
  FPList = MapAt[
    (#[[1]] #[[2]]) & @
      FCFeynmanParametrize[#, {l}, Names -> x] &,
    #, 1] & /@ tensorIntList;
  (#[[1]] #[[2]] // Contract) & /@ FPList // Apply[Plus] // #/Gamma[n] & // Simplify
];

FPOneLoopDenom[denom_, l_, x_, \[CapitalDelta]_] := Module[
  { momPart, massPart,
    b, xSum,
    lPart, \[CapitalDelta]mom, \[CapitalDelta]mass },
  { momPart, massPart } = FCI[denom] // Apply[List, #, {0,1}] & // Transpose;
  b = Length@momPart;
  DataType[x, FCVariable] = True;
  xSum = Sum[
    DataType[x[i], FCVariable] = True;
    x[i], {i, 1, b}];
  (* Feynman parametrizes `momPart` and separates it into `lPart` and `\[CapitalDelta]mom` *)
  { lPart, \[CapitalDelta]mom } = MapIndexed[SPD[#1] Apply[x][#2] &, momPart] //
    Apply[Plus] // CompleteSquare[#, l] & // ReplaceAll[xSum -> 1] //
    {SelectNotFree2[#, l], SelectFree2[#, l]} &;
  (* Feynman parametrizes `massPart` to `\[CapitalDelta]mass` *)
  \[CapitalDelta]mass = MapIndexed[#1^2 Apply[x][#2] &, massPart] // Apply[Plus];
  (* return the results : ( l^2 - \[CapitalDelta] )^b *)
  { l -> l - (lPart[[1]] // SelectFree2[#, l] & //
      ReplaceAll[Momentum[p_,___] -> p] // Simplify),
    \[CapitalDelta] -> - \[CapitalDelta]mom + \[CapitalDelta]mass //
      ExpandScalarProduct // FullSimplify[#, Assumptions -> xSum == 1] &,
    b }
];

ReduceOneLoopNumr::notsupp = "Not supported tensor `1`.";
ReduceOneLoopNumr[numr_, lRule_Rule] := Module[
  { l = lRule[[1]] },
  FCI[numr] //
  SUNSimplify //
  (* translate `l` *)
  FCReplaceMomenta[#, {lRule}] & //
  (* decompose `numr` into {tensor of `l`, prefactor} *)
  Append[
    (* non-const part *)
    # // FCTraceExpand // DiracGammaExpand // DotExpand // ExpandScalarProduct // Expand //
      SelectNotFree2[#, l] & //
      DiracSimplify // Simplify //
      Uncontract[#, l, Pair -> {l}] & //
      Collect2[#, l] & //
      If[Head@# === Plus, List@@#, List@#] & //
      Map[{SelectNotFree[#, l], SelectFree[#, l] // Simplify} &],
    (* const part, not simplified *)
    { 1, FCReplaceMomenta[#, {l -> 0}] }
  ] & //
  (* reduce tensors of `l` and return the power `a` of '(l^2)^a' *)
  Map@ReplaceAll@{
    (* Note the args in `Pair` are in alphabetical order! *)
    {0, _} -> Nothing, (* this case happens in scalar integrals *)
    {1, pref_} :> {0, pref},
    {Pair[LorentzIndex[__], M_], _} -> Nothing,
    {Pair[LI1_, M1_] Pair[LI2_, M2_], pref_} :> {1, Pair[LI1,LI2]/D pref // Contract},
    {Pair[p1__] Pair[p2__] Pair[p3__], _} -> Nothing,
    {n_, pref_} :> (Message[ReduceOneLoopNumr::notsupp, n]; {n, pref})
  } //
  (* collect terms of the same power of `l` *)
  GatherBy[#, First] & //
  Map[{#[[1, 1]], Plus@@(Transpose[#][[2]]) // Simplify} &]
];

StdOneLoop[b_, reducedNumr_, \[CapitalDelta]_] :=
  FCI[reducedNumr] //
  If[Head[#[[1]]] =!= List, {#}, #] & //
  (* calculate the standard one-loop integral *)
  Map@Replace[
    {a_, pref_} :>
      I (-1)^(a-b) / (4Pi)^(D/2) \[CapitalDelta]^(D/2+a-b) *
      Gamma[-D/2-a+b] Gamma[D/2+a] /( Gamma[b] Gamma[D/2] ) pref
  ] //
  (* gather the results *)
  Map[Apply[Times]] // Apply[Plus] // Simplify;

FPOneLoop[ampLoop_, l_, x_, \[CapitalDelta]_] := Module[
  { denom, numr,
    FPDenom, reducedNumr },
  { denom, numr } = {
    SelectNotFree[#, FeynAmpDenominator],
    SelectFree[#, FeynAmpDenominator]
  } & @ Collect[FCI[ampLoop], FeynAmpDenominator[__]];
  FPDenom = FPOneLoopDenom[denom, l, x, \[CapitalDelta]];
  reducedNumr = ReduceOneLoopNumr[numr, FPDenom[[1]]] //
    Map[MapAt[Calc /* Simplify, #, 2] &];
  { StdOneLoop[FPDenom[[3]], reducedNumr, \[CapitalDelta]],
    FPDenom[[2]] }
];

FPIntegrate[b_, expr_, x_, opts___?OptionQ] := Module[
  { xSum, int },
  xSum[n_] := Plus @@ Table[x[i], {i, 1, n}];
  int = expr /. x[b] -> 1 - xSum[b - 1];
  For[i = b - 1, i > 0, i--,
    int = Integrate[int, {x[i], 0, 1 - xSum[i - 1]}, opts]
  ];
  Gamma[b] int
];

End[];

EndPackage[];
