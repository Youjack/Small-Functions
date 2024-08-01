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

(* auxiliary functions *)

LogCombine::usage = 
  "LogCombine[expr, factor:1] combines all logs in expr
  (and takes the factor in the coefficient into the log)."

AddMassRegulator::usage = 
  "AddMassRegulator[amp, m] add mass m to every massless propagator in amp.";

Z2d::usage =
  "Z2d[expr, {{Z1,d1},{Z2,d2},...}] replaces Zi in expr with di,
  and expand di to first order.";

BetaFromZ::usage =
  "BetaFromZ[ZgMS, g, dimg] calculates the beta function of dimensionless g from ZgMS.
  (dimg is the dimension of g in D-dimensional spacetime)";

(* dimensional regularization *)

DRCoupling::usage =
  "DRCoupling[expr, g, dimg] multiplies g by (\[ExponentialE]^\[Gamma]\[Mu]/(4\[Pi]))^(dimg/2).
  (can only be applied to dimensionless couplings)";
DRCoupling\[Prime]::usage =
  "DRCoupling\[Prime][expr, g, dimg] multiplies g by \[Mu]^dimg.
  (can only be applied to dimensionless couplings)";

DRAddScaleMu::usage =
  "DRAddScaleMu[expr, dim, \[CapitalDelta], eps:\[Epsilon]] replaces D\[Rule]dim-2\[Epsilon] in expr,
  then replaces \[CapitalDelta]^(-\[Epsilon])\[Rule](\[ExponentialE]^\[Gamma] \[Mu]^2/(4\[Pi]\[CapitalDelta]))^\[Epsilon]";

DRExpand::usage =
  "DRExpand[expr, dim] replaces D\[Rule]dim-2\[Epsilon] in expr,
  expands expr with respect to \[Epsilon],
  and returns {\[Epsilon] dependent part, \[Epsilon] independent part}.";

DRExpandScaleless::usage =
  "DRExpandScaleless[expr, dim, zero\[CapitalDelta]] replaces D \[Rule] dim-2\[Epsilon] in expr,
  then replaces \[Mu]^(2\[Epsilon])zero\[CapitalDelta]^(-\[Epsilon]) \[Rule] Epsilon(1/\[Epsilon]UV - 1/\[Epsilon]IR).
  The return format is the same as DRExpand.";

(* loop calculations *)

FCFPOneLoop::usage =
  "FCFPOneLoop[int, l, x] applies `FCFeynmanParametrize` to compute one-loop tensor integral.
  (the final result have an implicit Feynman-parameters integral with Gamma[b])";

FPOneLoopDenom::usage =
  "FPOneLoopDenom[denom, l, x, \[CapitalDelta]] Feynman-parametrizes denom (an FAD),
  and return {l\[Rule]shifted l, \[CapitalDelta]\[Rule]..., power of denom}.";

ReduceOneLoopNumr::usage =
  "ReduceOneLoopNumr[numr, lRule] shifts loop momentum in numr according to lRule,
  reduces the loop-momentum tensors to scalars,
  and returns {{a1 in (l^2)^a1, prefactor1},...}.";

StdOneLoop::usage =
  "StdOneLoop[b, FPNumr, \[CapitalDelta]] gives the standard one-loop integral.
  (the final result have an implicit Feynman-parameters integral with Gamma[b])";

FPOneLoop::usage =
  "FPOneLoop[ampLoop, l, x, \[CapitalDelta]]
  (the final result have an implicit Feynman-parameters integral with Gamma[b])";

FPIntegrate::usage =
  "FPIntegrate[b, expr, x] integrates over b Feynman parameters x[i] in expr
  and multiplies it by Gamma[n].";

(* PlusDistribution *)

PlusDistributionExplicit::usage =
  "PlusDistributionExplicit[expr, x] explicitly calculates PlusDistributions of variable x in expr";

PlusDistributionExpand::usage =
  "PlusDistributionExpand[expr, x] expands PlusDistributions of variable x in expr
  with in terms of PlusDistribution[Log[1-x]^n/(1-x)].";

PlusDistributionIntegrate::usage =
  "PlusDistributionIntegrate[expr, x] integrates expr (may involve PlusDistributions) w.r.t x from 0 to 1."

Begin["`Private`"]; (*----------------------------------------------------------------------------*)

$lorentzIndices = {Global`\[Mu], Global`\[Nu], Global`\[Rho], Global`\[Sigma],
  Global`\[Kappa], Global`\[Lambda]};

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

(* auxiliary fucntions *)

PlusToList[sum_] := 
  If[sum === 0, {}, If[Head@sum === Plus, List@@sum, List@sum]];
TimesToList[prod_] := 
  If[prod === 0, {}, If[Head@prod === Times, List@@prod, List@prod]];

LogCombine[expr_, factor_ : 1] := Module[
  {noLog, haveLog, logList},
  {noLog, haveLog} = expr // PowerExpand // Collect[#, Log[_]] & // PlusToList //
    {SelectFree[#, Log], SelectNotFree[#, Log]} &;
  (* decomposd c Log[X] into {c,X} *)
  logList = haveLog // Map[{SelectFree[#, Log], Identity @@ SelectNotFree[#, Log]} &];
  (Plus @@ noLog) + If[logList === {}, 0,
    logList[[1, 1]] / factor *
    Log[Times @@ Map[#[[2]]^(factor #[[1]] / logList[[1, 1]] // Simplify) &, logList]]]
];

AddMassRegulator[amp_, m_] :=
  FCI@amp /. PropagatorDenominator[p_, 0] :> PropagatorDenominator[p, m];

Z2d[expr_, ZList_List] := Module[{c},
  expr //
  ReplaceAll[(#[[1]] -> 1 + c #[[2]]) & /@ ZList] //
  Normal@Series[#, {c, 0, 1}] & //
  ReplaceAll[c -> 1] //
  Collect[#, (#[[2]]) & /@ ZList, Simplify] &
];

BetaFromZ[ZgMS_, g_, dimg_] := Module[
  {dimgg, C1},
  dimgg = ((dimg /. D -> D - 2Epsilon) - dimg // Simplify) / Epsilon;
  C1 = g SeriesCoefficient[ZgMS - 1, {Epsilon, Infinity, 1}];
  - Epsilon dimgg g - dimgg C1 + dimgg g D[C1, g]
]

(* dimensional regularization *)

DRCoupling[expr_, g_, dimg_] :=
  expr /. g -> g ((E^EulerGamma ScaleMu^2)/(4 Pi))^(dimg/2);
DRCoupling\[Prime][expr_, g_, dimg_] :=
  expr /. g -> g (ScaleMu^2)^(dimg/2);

DRAddScaleMu[expr_, dim_, Delta_, eps_:Epsilon] :=
  FCReplaceD[expr, D -> dim - 2 eps] // Simplify //
  (* to deal with nonsymbol Deltas, Exponent is not suitable here *)
  ReplaceAll[{
    Delta^(     b_ * eps) /; b<0 :>           ((E^EulerGamma ScaleMu^2)/(4 Pi Delta))^(-b eps),
    Delta^(a_ + b_ * eps) /; b<0 :> Delta^a * ((E^EulerGamma ScaleMu^2)/(4 Pi Delta))^(-b eps)
  }];

DRExpand[expr_, dim_] :=
  FCReplaceD[expr, D -> dim - 2 Epsilon] //
  Normal@Series[#, {Epsilon,0,0}] & //
  {SelectNotFree2[#, Epsilon], SelectFree2[#, Epsilon]} & // Simplify;

DRExpandScaleless[expr_, dim_, Delta_] :=
  expr // Collect[#, Delta, Simplify] & //
  PlusToList // Map[
    Switch[
      Exponent[#, Delta] /. D -> 4 - 2 Epsilon // Simplify,
      0, {0,#}, (* do not change terms without Delta *)
      -Epsilon,
        SelectFree[#, Delta] //
        FCReplaceD[#, D -> dim - 2 Epsilon] & // Simplify //
        Normal@Series[ScaleMu^(-2Epsilon) Epsilon #, {Epsilon,0,0}] & //
        { # / Epsilon, - # / EpsilonIR } &,
      _, {0,0}
    ] &
  ] // Apply[Plus];

(* loop calculations *)

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

FPOneLoopDenom[denom_, l_, x_, Delta_] := Module[
  { momPart, massPart,
    b, xSum,
    lPart, DeltaMom, DeltaMass },
  { momPart, massPart } = FCI[denom] // Apply[List, #, {0,1}] & // Transpose;
  b = Length@momPart;
  DataType[x, FCVariable] = True;
  xSum = Sum[
    DataType[x[i], FCVariable] = True;
    x[i], {i, 1, b}];
  (* Feynman parametrizes `momPart` and separates it into `lPart` and `DeltaMom` *)
  { lPart, DeltaMom } = MapIndexed[SPD[#1] Apply[x][#2] &, momPart] //
    Apply[Plus] // CompleteSquare[#, l] & // ReplaceAll[xSum -> 1] //
    {SelectNotFree2[#, l], SelectFree2[#, l]} &;
  (* Feynman parametrizes `massPart` to `DeltaMass` *)
  DeltaMass = MapIndexed[#1^2 Apply[x][#2] &, massPart] // Apply[Plus];
  (* return the results : ( l^2 - Delta )^b *)
  { l -> l - (lPart[[1]] // SelectFree2[#, l] & //
      ReplaceAll[Momentum[p_,___] :> p] // Simplify),
    Delta -> - DeltaMom + DeltaMass //
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
    (* non-const part; could be empty *)
    # // FCTraceExpand // DiracGammaExpand // DotExpand // ExpandScalarProduct // Expand //
      SelectNotFree2[#, l] & //
      DiracSimplify // Simplify //
      Uncontract[#, l, Pair -> {l}] & //
      Collect2[#, l] & // PlusToList //
      Map[{SelectNotFree[#, l], SelectFree[#, l] // Simplify} &],
    (* const part, not simplified *)
    { 1, FCReplaceMomenta[#, {l -> 0}] }
  ] & //
  (* reduce tensors of `l` and return the power `a` of '(l^2)^a' *)
  Map@ReplaceAll@{
    (* Note the args in `Pair` are in alphabetical order! *)
    {1, pref_} :> {0, pref},
    {Pair[LorentzIndex[__], M_], _} -> Nothing,
    {Pair[LI1_, M1_] Pair[LI2_, M2_], pref_} :> {1, Pair[LI1,LI2]/D pref // Contract},
    {Pair[p1__] Pair[p2__] Pair[p3__], _} -> Nothing,
    {Pair[LI1_, M1_] Pair[LI2_, M2_] Pair[LI3_, M3_] Pair[LI4_, M4_], pref_} :> {2,
      (Pair[LI1,LI2] Pair[LI3,LI4] + Pair[LI1,LI3] Pair[LI2,LI4] + Pair[LI1,LI4] Pair[LI3,LI2]) /
      (D(D+2)) pref // Contract},
    {Pair[p1__] Pair[p2__] Pair[p3__] Pair[p4__] Pair[p5__], _} -> Nothing,
    {n_, pref_} :> (Message[ReduceOneLoopNumr::notsupp, n]; {n, pref})
  } //
  (* collect terms of the same power of `l` *)
  GatherBy[#, First] & //
  Map[{#[[1, 1]], Plus@@(Transpose[#][[2]]) // Simplify} &]
];

StdOneLoop[b_, FPNumr_, Delta_] :=
  FCI[FPNumr] //
  If[Head[#[[1]]] =!= List, {#}, #] & //
  (* calculate the standard one-loop integral *)
  Map@Replace[
    {a_, pref_} :>
      I (-1)^(a-b) / (4Pi)^(D/2) Delta^(D/2+a-b) *
      Gamma[-D/2-a+b] Gamma[D/2+a] /( Gamma[b] Gamma[D/2] ) pref
  ] //
  (* gather the results *)
  Map[Apply[Times]] // Apply[Plus] // Simplify;

FPOneLoop[ampLoop_, l_, x_, Delta_] := Module[
  { denom, numr,
    FPDenom, FPNumr },
  { denom, numr } = {
    SelectNotFree[#, FeynAmpDenominator],
    SelectFree[#, FeynAmpDenominator]
  } & @ Collect[FCI[ampLoop], FeynAmpDenominator[__]];
  FPDenom = FPOneLoopDenom[denom, l, x, Delta];
  FPNumr = ReduceOneLoopNumr[numr, FPDenom[[1]]] //
    Map[MapAt[Calc /* Simplify, #, 2] &];
  { StdOneLoop[FPDenom[[3]], FPNumr, Delta],
    FPDenom[[2]] }
];

FPIntegrate[b_, expr_, x_, opts___?OptionQ] := Module[
  { xSum, int },
  xSum[n_] := Plus @@ Table[x[i], {i, 1, n}];
  int = expr /. x[b] -> 1 - xSum[b - 1];
  Do[
    int = Integrate[int, {x[i], 0, 1 - xSum[i - 1]}, opts],
    {i, b-1, 1, -1}];
  Gamma[b] int
];

(* PlusDistribution *)

(
  PlusDistribution[Log[x_ * (1 - x_)]/(1 - x_)] =.;
  PlusDistribution[PlusDistribution[F_]] := PlusDistribution[F];
  PlusDistribution[f_ * PlusDistribution[F_]] := PlusDistribution[f * F];
  PlusDistribution[DiracDelta[_]] := 0;
  PlusDistribution[f_ * DiracDelta[_]] := 0;
);

PlusDistributionExplicit[expr_, x_] := Module[{xx},
  expr /. PlusDistribution[F_] :> F - DiracDelta[1-x] Integrate[F /. x -> xx, {xx, 0, 1}]
];

PlusDistributionExpand[expr_, x_] :=
  expr /. PlusDistribution[F_] :> Module[
    { FList,
      InvF, LogF, PD,
      xx, havePD, noPD},
    FList = F // 
      ReplaceAll[log_Log /; !FreeQ[log,x] :> PowerExpand@log] // Expand // PlusToList //
      ReplaceAll@{(1-x)^-1 -> InvF, (x-1)^-1 -> -InvF, Log[1-x] -> LogF} //
      ReplaceAll@{LogF^n_ InvF :> PD[n], LogF InvF -> PD[1], InvF -> PD[0]};
    havePD = SelectNotFree[FList, PD] //
      Map[{SelectNotFree[#, PD] /. PD[n_] :> n, SelectFree[#, PD]} &] //
      Map[If[(#[[2]] /. x -> 1) === 0,
        PlusDistributionExplicit[PlusDistribution[#[[2]] Log[1-x]^#[[1]] / (1-x)], x],
        #[[2]] PlusDistribution[Log[1-x]^#[[1]] / (1-x)] - DiracDelta[1-x] Integrate[
          Log[1-xx]^#[[1]] / (1-xx) *
          ((#[[2]] /. x -> xx) - (#[[2]] /. x -> 1)),
          {xx, 0, 1}]
      ] &];
    noPD = SelectFree[FList, PD] //
      Map[PlusDistributionExplicit[PlusDistribution[#], x] &];
    Plus @@ havePD + Plus @@ noPD
  ] //
  Collect[#, {PlusDistribution[_], EpsilonIR, DiracDelta[_], Log[_]}, Simplify] &;

PlusDistributionIntegrate[expr_, x_] :=
  Collect[expr, {PlusDistribution[_], DiracDelta[_]}] // PlusToList //
  Map[ReplaceAll@{
    f_ * PlusDistribution[F_] :> Integrate[F (f - (f /. x -> 1)), {x,0,1}],
    PlusDistribution[_] -> 0,
    f_ * DiracDelta[1-x] :> (f /. x -> 1),
    DiracDelta[1-x] -> 1,
    f_ :> Integrate[f, {x,0,1}]
  }] // Apply[Plus] // Simplify;

End[];

BeginPackage["YoujackOneLoop`FormFactorTools`"]; (*===============================================*)
Needs["FeynCalc`"];

DecomposeOSVertex::usage =
  "DecomposeOSVertex[expr, spinorL, spinorR, qin] contracts expr as spinorL.expr.spinorR,
  simplies it with Gordon identity, then removes the spinors,
  and returns {terms with \[Gamma], terms with \[Sigma], other terms}.
  (qin is the momentum of the incoming photon)";

Begin["`Private`"];

DecomposeOSVertex[expr_, spinorL_, spinorR_, qin_] := Module[
  { momL, momR, m, pp },
  momL = FCI[spinorL][[1]]; momR = FCI[spinorR][[1]];
  m = FCI[spinorL][[2]];
  expr //
    spinorL.#.spinorR & // DiracSimplify // Simplify //
    ReplaceAll@{Spinor[__].Spinor[__] -> 1, Spinor[__] -> Sequence[]} //
    ReplaceAll@{
      momR -> Momentum[(pp - qin)/2, D], -momR -> -Momentum[(pp - qin)/2, D],
      momL -> Momentum[(pp + qin)/2, D], -momL -> -Momentum[(pp + qin)/2, D]
    } // ExpandScalarProduct // Simplify //
    ReplaceAll[FCI@FVD[pp, LI_] :>
      2 m GAD[LI] - I DiracSigma[GAD[LI], GSD[qin]] (* Gordon identity *)
    ] // ExpandScalarProduct //
    { SelectNotFree2[SelectFree2[#, DiracSigma], DiracGamma] // Simplify,
      SelectNotFree2[#, DiracSigma] // Simplify,
      SelectFree2[#, DiracGamma, DiracSigma] // Simplify } &
];

End[];

EndPackage[];

BeginPackage["YoujackOneLoop`DiracGammaTools`"]; (*===============================================*)
Needs["FeynCalc`"];

DiracDecompose::usage =
  "DiracDecompose[expr] decomposes the product of Dirac gamma matrices in terms of the conventional basis.";

Begin["`Private`"];

DiracDecompose[expr_] := Module[
  {LI1, LI2},
                                  1/4 DiracTrace[                              expr] +
      GA[LI1]                     1/4 DiracTrace[GA[LI1]                     . expr] +
      GA5                         1/4 DiracTrace[GA5                         . expr] +
      GA[LI1].GA5                 1/4 DiracTrace[GA5.GA[LI1]                 . expr] +
  1/2 DiracSigma[GA[LI1],GA[LI2]] 1/4 DiracTrace[DiracSigma[GA[LI1],GA[LI2]] . expr] //
    DiracSubstitute67 // DiracSimplify[#, ToDiracGamma67->False] & // FCE //
    ReplaceAll[GA[LI1_?(#=!=5&)].GA[LI2_?(#=!=5&)] :>
      ToDiracSigma[GA[LI1].GA[LI2], GA[LI1], GA[LI2]]] //
    FullSimplify // Expand
];

End[];

EndPackage[];

EndPackage[];
