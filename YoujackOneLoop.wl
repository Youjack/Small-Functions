(* ::Package:: *)

$LoadAddOns = {"FeynArts"};
BeginPackage["YoujackOneLoop`",{"FeynCalc`"}];
Needs["FeynArts`"];
$FAVerbose = 0;

Set1PIOptions::usage =
  "amplitudes are M without an \[ImaginaryI];
  external lines are truncated.";

AddMassRegulator::usage = 
  "AddMassRegulator[amp, m] add mass m to every massless propagator in amp.";

Z2d::usage =
  "Z2d[expr, {{Z1,d1},{Z2,d2},...}] replaces Zi in expr with di,
  and expand di to first order.";

DRCoupling::usage =
  "DRCoupling[expr, g] replace the coupling g with 4\[Pi].";

DRExpand::usage =
  "DRExpand[expr, dim] replaces D\[RightArrow]dim-2\[Epsilon] in expr,
  expands expr with respect to \[Epsilon],
  and puts terms of different order into a List.";

FCFPOneLoop::usage =
  "FCFPOneLoop[int, l, x] applies `FCFeynmanParametrize` to compute one-loop tensor integral.
  (the final result have an implicit Feynman-parameters integral with Gamma[n])";

FPOneLoopDenom::usage =
  "FPOneLoopDenom[denom, l, x, \[CapitalDelta]] Feynman-parametrizes denom,
  and return {l\[RightArrow]shifted l, \[CapitalDelta]\[RightArrow]..., power of denom}.";

ReduceOneLoopNumr::usage =
  "ReduceOneLoopNumr[numr, lRule] shifts loop momentum in numr using lRule,
  reduces the loop-momentum tensors to scalars,
  and return {{loop-momentum scalar1,prefactor1},...}.";

StdOneLoop::usage =
  "StdOneLoop[b, reducedNumr, \[CapitalDelta]] gives the standard one-loop integral.
  (the final result have an implicit Feynman-parameters integral with Gamma[n])";

FPOneLoop::usage = 
  "FPOneLoop[ampLoop, l, x, \[CapitalDelta]]
  (the final result have an implicit Feynman-parameters integral with Gamma[n])";

FPIntegrate::usage =
  "FPIntegrate[b, expr, x] integrates over b Feynman parameters x[i] in expr and multiplies Gamma[n].";

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

AddMassRegulator[amp_, m_] :=
  ReplaceAll[amp, PropagatorDenominator[p_, 0] :> PropagatorDenominator[p, m]];

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

ReduceOneLoopNumr::notsupp = "Not supported tensor `1`.";
ReduceOneLoopNumr[numr_, lRule_] := Module[
  { l = lRule[[1]] },
  FCI[numr] //
    (* translate `l` *)
    ReplaceAll[lRule] //
    (* decompose `numr` into {power of `l`,prefactor} *)
    DiracSimplify //
    Uncontract[#, l, Pair -> {l}] & //
    Collect2[#, l] & // If[Head@# === Plus, List@@#, List@#] & //
    Map[{SelectNotFree[#, l], SelectFree[#, l] // Simplify} &] //
    (* reduce tensors of `l` *)
    Map@ReplaceAll@{
      (* Note the args in `Pair` are in alphabetical order! *)
      {1, pref_} :> {1, pref},
      {Pair[LorentzIndex[__], M_], _} -> Nothing,
      {Pair[LI1_, M1_] Pair[LI2_, M2_], pref_} :> {FCI@SPD[l], Pair[LI1,LI2]/D pref // Contract},
      {n_, _} :> (Message[ReduceOneLoopNumr::notsupp, n]; Abort[])
    } //
    (* collect terms of the same power of `l` *)
    GatherBy[#, First] & //
    Map[{#[[1, 1]], Plus@@(Transpose[#][[2]]) // Simplify} &]
];

StdOneLoop::notsupp = "Not supported numerator `1`.";
StdOneLoop[b_, reducedNumr_, \[CapitalDelta]_] :=
  FCI[reducedNumr] //
    (* calculate the standard one-loop integral *)
    Map@Replace@{
      {1, pref_} :>
        I (-1)^(b  ) \[CapitalDelta]^(D/2-b  )/(4\[Pi])^(D/2) Gamma[-D/2+b]/Gamma[b] pref,
      {Pair[__], pref_} :>
        I (-1)^(b+1) \[CapitalDelta]^(D/2-b+1)/(4\[Pi])^(D/2) D/2 Gamma[-D/2+b-1]/Gamma[b] pref,
      {n_, _} :> (Message[StdOneLoop::notsupp, n]; Abort[])
    } //
    (* gather the results *)
    Map[Apply[Times]] // Apply[Plus] // Simplify;

FPOneLoop[ampLoop_, l_, x_, \[CapitalDelta]_] := Module[
  { denom, numr,
    FPdenom, reducedNumr },
  { denom, numr } = {
    SelectNotFree[#, FeynAmpDenominator],
    SelectFree[#, FeynAmpDenominator]
  } & @ FCI[ampLoop];
  FPdenom = FPOneLoopDenom[denom, l, x, \[CapitalDelta]];
  reducedNumr = ReduceOneLoopNumr[numr, FPdenom[[1]]];
  { StdOneLoop[FPdenom[[3]], reducedNumr, \[CapitalDelta]],
    FPdenom[[2]] }
];

FPIntegrate[b_, expr_, x_] := Module[
  { xSum, int },
  xSum[n_] := Plus @@ Table[x[i], {i, 1, n}];
  int = expr /. x[b] -> 1 - xSum[b - 1];
  For[i = b - 1, i > 0, i--,
    int = Integrate[int, {x[i], 0, 1 - xSum[i - 1]}]
  ];
  Gamma[b] int
];

End[];

EndPackage[];
