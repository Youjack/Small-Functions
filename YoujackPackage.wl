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
    range_, axesQty_?ListQ, axesDim_?ListQ] :=
    Module[
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

End[];

EndPackage[];
