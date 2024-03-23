(* ::Package:: *)

$LoadAddOns = {"FeynArts"};
BeginPackage["YoujackOneLoop`",{"FeynCalc`"}];
Needs["FeynArts`"];
$FAVerbose = 0;

Set1PIOptions::usage =
  "amplitudes are M without an \[ImaginaryI];
  external lines are truncated.";

Z2d::usage =
  "Z2d[expr, {{Z1,d1},{Z2,d2},...}] replaces Zi in expr with di,
  and expand di to first order."

DRCoupling::usage =
  "DRCoupling[expr, g] replace the coupling g with 4\[Pi].";

DRExpand::usage =
  "DRExpand[expr, dim] replaces D\[RightArrow]dim-2\[Epsilon] in expr,
  expands expr with respect to \[Epsilon],
  and puts terms of different order into a List.";

FCFPOneLoop::usage =
  "FCFPOneLoop[int, l, x] applies `FCFeynmanParametrize` to compute one-loop tensor integral.";

FPOneLoopDenom::usage =
  "FPOneLoopDenom[denom, l, x, \[CapitalDelta]] Feynman-parametrizes denom,
  and return {l\[RightArrow]shifted l, \[CapitalDelta]\[RightArrow]..., power of denom}.
  (denom is usually ampLoop[[-1]])";

ReduceOneLoopNumr::usage =
  "ReduceOneLoopNumr[numr, lRule] shifts loop momentum in numr using lRule,
  reduces the loop-momentum tensors to scalars,
  and return {{loop-momentum scalar1,prefactor1},...}.
  (numr is usually Delete[ampLoop,-1])";

StdOneLoop::usage =
  "StdOneLoop[b, reducedNumr, \[CapitalDelta]] gives the standard one-loop integral.";

FPOneLoop::usage = 
  "FPOneLoop[ampLoop, l, x, \[CapitalDelta]]"

Begin["`Private`"];

Set1PIOptions[] := (
  $KeepLogDivergentScalelessIntegrals = True;

  SetOptions[CreateTopologies,
    ExcludeTopologies -> {Tadpoles, WFCorrections, WFCorrectionCTs, Reducible}];
  SetOptions[Paint, Numbering -> Simple, SheetHeader -> None, ColumnsXRows -> {2, 1}];
  SetOptions[Render, ImageSize -> {2 128, 128}];

  SetOptions[CreateFeynAmp, PreFactor -> -I, Truncated -> True];
  SetOptions[FCFAConvert, List -> False, SMP -> True,
    LorentzIndexNames -> {Global`\[Mu], Global`\[Nu], Global`\[Rho], Global`\[Sigma]},
    ChangeDimension -> D,
    UndoChiralSplittings -> True, Contract -> True];
  SetOptions[FCFeynmanParametrize, FeynmanIntegralPrefactor -> "Textbook"];
);

Z2d[expr_, ZList_List] := Module[{\[Alpha]},
  ReplaceAll[expr, (#[[1]] -> 1 + \[Alpha] #[[2]]) & /@ ZList] //
    Normal@Series[#, {\[Alpha], 0, 1}] & //
    ReplaceAll[\[Alpha] -> 1] //
    Collect[#, (#[[2]]) & /@ ZList, Simplify] &
];

DRCoupling[expr_, g_] :=
  ReplaceAll[expr, g -> g ((4 Pi E^-EulerGamma)^(-1/2) ScaleMu)^Epsilon];

DRExpand[expr_, dim_] :=
  FCReplaceD[expr, D -> dim - 2 Epsilon] //
  Series[#, {Epsilon, 0, 0}] & //
  Table[#[[3,n]] Epsilon^(#[[4]] + n - 1), {n,1,Length@#[[3]]}] & //
  Simplify;

FCFPOneLoop[int_, l_, x_] := Module[
  {tensorIntList, FPList},
  (* decompose `int` into {{int1,factor1},{int2,factor2},...} *)
  tensorIntList = int //
    DiracSimplify // Simplify //
    Expand2[#, l] & //
    If[Head@# === Plus, List@@#, List@#] & //
    Uncontract[#, l, Pair -> {l}] & //
    ({SelectNotFree[#, l], SelectFree[#, l]} &) /@ # &;
  (* apply `FCFeynmanParametrize` to `inti` *)
  FPList = MapAt[
    (#[[1]] #[[2]]) & @
      FCFeynmanParametrize[#, {l}, Names -> x] &,
    #, 1] & /@ tensorIntList;
  (#[[1]] #[[2]] // Contract) & /@ FPList // Plus@@# & // Simplify
];

FPOneLoopDenom[denom_, l_, x_, \[CapitalDelta]_] := Module[
  { momPart, massPart,
    b, xSum,
    lPart, \[CapitalDelta]mom, \[CapitalDelta]mass },
  { momPart, massPart } = FCI[denom] // Apply[List, #, {0,1}] & // Transpose;
  b = Length@momPart;
  xSum = Plus @@ Table[
    DataType[x[i], FCVariable] = True;
    x[i], {i, 1, b}];
  (* Feynman parametrizes `momPart` and separates it into `lPart` and `\[CapitalDelta]mom` *)
  { lPart, \[CapitalDelta]mom } = MapIndexed[SPD[#1] Apply[x][#2] &, momPart] //
    Apply[Plus] // CompleteSquare[#, l] & // ReplaceAll[xSum -> 1] //
    {SelectNotFree[#, l], SelectFree[#, l]} &;
  (* Feynman parametrizes `massPart` to `\[CapitalDelta]mass` *)
  \[CapitalDelta]mass = MapIndexed[#1^2 Apply[x][#2] &, massPart] // Apply[Plus];
  (* return the results : ( l^2 - \[CapitalDelta] )^b *)
  { l -> l - (lPart[[1]] // SelectFree[#, l] & // Simplify),
    \[CapitalDelta] -> - \[CapitalDelta]mom + \[CapitalDelta]mass //
      ExpandScalarProduct // FullSimplify[#, Assumptions -> xSum == 1] &,
    b }
];

ReduceOneLoopNumr[numr_, lRule_] := Module[
  { l = lRule[[1]] },
  numr //
    (* translate `l` *)
    ReplaceAll[lRule] //
    (* decompose `numr` into {power of `l`,prefactor} *)
    DiracSimplify // Simplify //
    Expand2[#, l] & // If[Head@# === Plus, List@@#, List@#] & //
    Uncontract[#, l, Pair -> {l}] & //
    {SelectNotFree[#, l], SelectFree[#, l]} & /@ # & //
    (* reduce tensors of `l` *)
    Map@ReplaceAll@{
      {Pair[idx1_, __] Pair[idx2_, __], pref_} :> {Pair[idx1, idx2]/D FCI@SPD[l], pref},
      {Pair[LorentzIndex[__], __], __} -> Nothing
    } //
    (* collect terms of the same power of `l` *)
    Map[Apply[Times]/*Contract] // Apply[Plus] // Simplify //
    Collect2[#, l] & // If[Head@# === Plus, List@@#, List@#] & //
    {SelectNotFree[#, l], SelectFree[#, l] // Simplify} & /@ # &
];

StdOneLoop[b_, reducedNumr_, \[CapitalDelta]_] :=
  FCI[reducedNumr] //
    (* calculate the standard one-loop integral *)
    Map@Replace@{
      {Pair[__], pref_} :> {
        I (-1)^(b+1) \[CapitalDelta]^(D/2-b+1)/(4\[Pi])^(D/2) D/2 Gamma[-D/2+b-1]/Gamma[b],
        pref},
      {1, pref_} :> {
        I (-1)^(b  ) \[CapitalDelta]^(D/2-b  )/(4\[Pi])^(D/2) Gamma[-D/2+b]/Gamma[b],
        pref}
    } //
    (* multiply by Gamma[b] in Feynman parametrization *)
    Map[Apply[Times]] // Apply[Plus] // Gamma[b] # & // Simplify;

FPOneLoop[ampLoop_, l_, x_, \[CapitalDelta]_] := Module[
  {FPdenom, reducedNumr},
  FPdenom = FPOneLoopDenom[ampLoop[[-1]], l, x, \[CapitalDelta]];
  reducedNumr = ReduceOneLoopNumr[Delete[ampLoop,-1], FPdenom[[1]]];
  { StdOneLoop[FPdenom[[3]], reducedNumr, \[CapitalDelta]],
    FPdenom[[2]] }
];

End[];

EndPackage[];
