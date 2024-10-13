BeginPackage["YoujackPackage`"];

(* Styles *)
ItalicStyle;
PlainStyle;
NumStyle;

(* LinearFitPlot *)
LinearFitPlot::usage =
  "LinearFitPlot[{labels},
  {data}, IncludeConstantBasis->False,
  range, axesQty, axesDim]
  = {TableForm@info, figure}";

(* Wolfram Player *)
WolframPlayer::usage = "WolframPlayer[expr]\nWolframPlayer[expr, name]";

(* InteractiveMapAt *)
InteractiveMapAt::usage = "";

Begin["`Private`"];

(* Styles *)
$YoujackPlotColor = ColorData[97,"ColorList"];
(* SetAttributes[YoujackMathForm, HoldAll]; *)
(* YoujackMathForm[expr_] := Style[TraditionalForm@HoldForm@expr, FontFamily->"Times"] *)
ItalicStyle[expr_] := Style[expr, Italic, FontFamily->"Times"];
PlainStyle[expr_] := Style[expr, Plain, FontFamily->"Times"];
NumStyle[num_?NumberQ] := PlainStyle@num;

(* LinearFitPlot *)
LinearFitPlot[labels_?ListQ,
  data_?ListQ, includeConstantBasis_?OptionQ,
  range_, axesQty_?ListQ, axesDim_?ListQ] := Module[
    {
      len,
      model, func, info,
      color, listPlot, plot, axesLabel, figure
    },

    len = Length@labels;

    (* Fit *)
    model = Table[LinearModelFit[data[[n]], x,x, includeConstantBasis], {n,1,len}];
    func = Table[model[[n]]@"BestFit", {n,1,len}];
    info = Table[{
      labels[[n]],
      func[[n]] /. {x->"x"},
      Row@{Superscript["R",2],"=",model[[n]]@"RSquared"}
    }, {n,1,len}];

    (* Plot *)
    color = Table[$YoujackPlotColor[[n]], {n,1,len}];
    listPlot = ListPlot[data, PlotLegends->labels, PlotStyle->color];
    plot = Plot[func, {x,range[[1]],range[[2]]}, PlotStyle->color];
    axesLabel = Table[DisplayForm@Row@{
      ItalicStyle@axesQty[[i]], PlainStyle@Row@{" (",axesDim[[i]],")"}
    }, {i,1,2}];
    figure = Show[listPlot, plot,
      PlotRange->All, PlotRangePadding->None, AxesOrigin->{0,0},
      AxesLabel->axesLabel, GridLines->Automatic];

    (* Return *)
    {TableForm@info, figure}
  ];

(* WolframPlayer *)
WolframPlayer[expr_, name_?StringQ] :=
  With[
    { dir = Evaluate@FileNameJoin@{$UserDocumentsDirectory, "Wolfram Player", name<>".cdf"} },
    Export[dir, #, "CDF"]& @ Notebook[{Cell[BoxData@ToBoxes@expr,"Output"]}, WindowSize->All];
    StartProcess@{"WolframPlayer", "\""<>dir<>"\""};
    ToString@Head@expr<>" in Wolfram Player"
  ];
WolframPlayer[expr_] := WolframPlayer[expr, CreateUUID["CDFOutput-"]];

(* InteractiveMapAt *)
InteractiveMapAt::wrongsize = "The argument list has a wrong size.";
InteractiveMapAt::backat0 = "Already back to depth\[Hyphen]0.";
InteractiveMapAt[fSeq___, OptionsPattern[{Print -> True}]][expr_] := Module[
  {
    fList, fListLen,
    exprList, posList, exprP, fexprP, posString,
    LocalVarColor = RGBColor[0.235, 0.49, 0.568],
    HeadColor = RGBColor[1., 0.72, 0.]
  },
  If[OddQ[fListLen = Length[fList = {fSeq}]],
    Message[InteractiveMapAt::wrongsize]; Abort[]];
  exprList = {expr}; posList = {};
  Do[
    Which[
      fList[[2i-1]] === 0, (
        posString = {"Depth[", Length@posList, "] = "};
        exprP = exprList[[-1]];
        fexprP = fList[[2i]][exprP];
        exprList[[-1]] = fexprP;
      ),
      fList[[2i-1]] === Back, If[posList =!= {},
        posString = {"Back to Depth[", Length@posList - 1, "] = "};
        exprP = Head[exprList[[-2]]][
          Delete[exprList[[-2]], List /@ posList[[-1]]],
          exprList[[-1]]
        ];
        fexprP = fList[[2i]][exprP];
        exprList = Delete[exprList, -1]; exprList[[-1]] = fexprP;
        posList = Delete[posList, -1],
        (* posList === {} *)
        Message[InteractiveMapAt::backat0]; Continue[]
      ],
      True, (
        posString = {
          "Depth[", Length@posList + 1, "] \[Congruent] ",
          "Depth[", Length@posList, "]",
          "\[LeftDoubleBracket]", fList[[2 i - 1]], "\[RightDoubleBracket] = "};
        exprP = exprList[[-1, fList[[2i-1]]]];
        fexprP = fList[[2i]][exprP];
        AppendTo[exprList, fexprP];
        AppendTo[posList, fList[[2i-1]]];
      )
    ];
    If[OptionValue[Print], Print @@
      (posString // Map[Style[#,LocalVarColor]&]) ~Join~
      {Style[Head@exprP,HeadColor], Style[List@@exprP,Black]} ~Join~
      If[fList[[2 i]] === Identity, {}, {
        Style[" \[Rule] ",LocalVarColor],
        Style[Head@fexprP,HeadColor], Style[List@@fexprP,Black]
      }]
    ],
    {i, fListLen/2}];
  Do[
    exprList[[i]] = Head[exprList[[i]]][
      Delete[exprList[[i]], List /@ posList[[i]]],
      exprList[[i + 1]]
    ],
    {i, Length@posList, 1, -1}];
  exprList[[1]]
];

End[];

EndPackage[];
