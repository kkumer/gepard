(* $Id:*)

(*     GeParD - Mathematica interface    *)
(*     ==============================    *)


Print["GeParD - Mathematica interface (2008-05-15)"];


If[$VersionNumber<5.999,  (* Mathematica 5.*)
BeginPackage["gepard`", "Format`", "NumericalMath`NLimit`", "Graphics`Graphics`",
                        "Utilities`FilterOptions`"], 
                     (* else Mathematica 6.* *)
BeginPackage["gepard`", "Format`", "NumericalCalculus`", "Utilities`FilterOptions`"]]



GepardInit::usage = "GepardInit[] performs initialization of GeParD parameters. 
It reads GEPARD.INI for default parameters. Parameters SPEED, P, SCHEME, ANSATZ,
DATFILE and OUTFILE can be overriden by specifying corresponding options 
e.g. GepardInit[P->0, DATFILE->\"dvcs\"]."

GepardInitInternal::usage = "MathLink function ..."

GepardFit::usage = "GepardFit[{par1, par2, ...}, options] does the complete
fitting of data to variable parameter {par1, par2, ...}, using ansatz depending
on the complete set of (variable and fixed) parameters specified in array
Parameters. Options are the same as for GepardInit"

PrettyStatus::usage = "PrettyStatus[] prints out the present status of fit and parameters in
the \"pretty\" form."

JustCali3::usage = "JustCali3[{par1, par2, ...}, options] is a simplified version of GepardFit which
doesn't do fitting but just printing out of initial value of chi-square. Useful for debugging.
Options are the same as for GepardInit"

CompileMoments::usage = "CompileMoments creates compiled Mathematica objects corresponding to 
moments of GPDs and PDFs as defined in GPDMom[j,t,xi] and PDFMom[j,t,xi]. These are then
available as ANSATZ 'MMA'."

SpliceToFortran::usage = "SpliceToFortran[path-to-src_, substitutions_] creates Fortran subroutine
SPLICE in file splice.f. This subroutine calculates moments of GPDs and PDFs as
defined in Mathematica functions GPDMom[j,t,xi] and PDFMom[j,t,xl]. path-to-src is string
with path to gepard/src directory, and substitutions_ are list of replacement rules which should
take care of some differences in names of Mathematica and Fortran variables and functions.
After recompiling gepard.exe this ansatz is available as ANSATZ 'SPLICE'."

Parameters::usage = "Parameters is an array of the form {{id1, par1, val1, step1, low1, hi1},
{id2, par2, ...}, ...}. Here idn is the parameter number, par1 is parameter name, val1 is
its starting value, etc. as expected by Minuit. See Minuit manual."

SaveParameters::usage = "SaveParameters[] writes present values of fitting
parameters (which are only recorded in the internal Minuit-Gepard memory)
into the list Parameters, preparing it thus for next fitting procedure."

ParameterID::usage = "ParameterID[symbol] returns parameter number (id) of a parameter with
name symbol, as specified by array Parameters."

ParameterName::usage = "ParameterName[id] returns name string corresponding to parameter number id, as specified by array Parameters."

GPD::usage = "GPD[{val1, val2, ...}, t, xi] is a function that is contacted by MathLink Minuit
interface and should give GPD's for .... GPD[flavor, x, t, xi] gives value of x-space
GPDs (flavor: 1=Q, 2=G). It relies on moments GPDcurrent[j,t,xi] which have to be set up e.g.
by GPDcurrent[j_, t_, xi_] = GPDMom[j, t, xi] /. PAR[n_] :> First[MinuitGetParameter[n]]"

PDF::usage = "PDF[{val1, val2, ...}, t, xi] is a function that is contacted by MathLink Minuit
interface and should give PDF's for .... "

slope::usage = "slope[flavor, x] gives GPD slope at x for flavor 1=Q or 2=G"

plotslopes::usage = "plots slope[flavor, x] for flavor 1=Q or 2=G"

plotPDFs::usage = "plots GPD[flavor, x, 0, 0] for flavor 1=Q or 2=G"

ChiSquareProbability::usage = "ChiSquareProbability[d, chisq] gives probability that
chi-square for a system with d degrees of freedom will be larger than or equal to chisq"

ReducedChiSquareProbability::usage = "ReducedChiSquareProbability[d, chisq] gives probability that
reduced chi-square (chi-square divided by number of degrees of freedom) for a system with
d d.o.f will be larger than or equal to chisq"

GPDMom::usage = "GPDMom ... should be defined by user"
PDFMom::usage = "PDFMom ... should be defined by user"

GPDMomCompiledQ::usage = "... produced by CompileMoments."
GPDMomCompiledG::usage = "... produced by CompileMoments."
GPDcurrent::usage = "GPDcurrent[j_, t_, xi_] = GPDMom[j, t, xi] /. PAR[n_] :> First[MinuitGetParameter[n]] is how it could be set."
jValues::usage = "jValues is list of complex coordinates on Mellin-Barnes contour 
set by call to MinuitInit[1] or GepardFit."
PAR::usage = "PAR[n] is the symbol for n-th fitting parameter."

MinuitInit::usage = "MinuitInit"
MinuitSetParameter::usage = "MinuitSetParameter"
MinuitGetParameter::usage = "MinuitGetParameter"
MinuitCommand::usage = "MinuitCommand[command_] executes MINUIT command, which should be specified as a string."
MinuitStatus::usage = "MinuitStatus"
GetChiSquares::usage = "GetChiSquares"
PrintMinuitCommand::usage = "MinuitCommand[command_] executes MINUIT
command and reads MINUIT output from file fit.mnt and prints it. 
MinuitCommand[command_, fontsize_] does the printing with font
size fontsize (values 4-10 are nice, 5 is default)."

SPEED::usage = "GeParD parameter"
P::usage = "GeParD parameter"
SCHEME::usage = "GeParD parameter"
ANSATZ::usage = "GeParD parameter"

DATFILE::usage = "GeParD parameter. Filename without extension."
OUTFILE::usage = "GeParD parameter. Filename without extension."

GPDQsplice::usage = "Fortran format function to be spliced into splice.f"
GPDGsplice::usage = "Fortran format function to be spliced into splice.f"
PDFQsplice::usage = "Fortran format function to be spliced into splice.f"
PDFGsplice::usage = "Fortran format function to be spliced into splice.f"

FCM::usage = "Fortran format array holding the expressions for GPDs"

j::usage = "j - complex conformal moment"
t::usage = "t - Mandelstam variable. In Fortran represented as DEL2"
xi::usage = "xi - DVCS kinematical variable."
xbj::usage = "xbj - DIS kinematical variable x_Bjorken."

lobj::usage = "lobj - MathLink link"

cffH::usage = "cffH[xi, t, q2, q02, options] returns singlet CFF H(xi, t, q2, q02). 
Options are same as for GepardInit."

cffE::usage = "cffE[xi, t, q2, q02, options] returns singlet CFF E(xi, t, q2, q02). 
Options are same as for GepardInit."

F2::usage = "F2[xbj, q2, q02, options] returns singlet F2(xbj, q2, q02). 
Options are same as for GepardInit."

cffHInternal::usage = "MathLink function ..."
cffEInternal::usage = "MathLink function ..."
F2Internal::usage = "MathLink function ..."
cInt::usage = "MathLink function ..."
MinuitContour::usage = "MathLink function ..."

AllParameterValues::usage = "AllParameterValues[] returns complete list that corresponds to 
COMMON block PAR. For internal use."

PlotMinuitContour::usage = "PlotMinuitContour[par1, par2, npts, options] plots two-parameter correlation contour produced by Minuit's MNCONT. par1 and par2 are symbols from Parameters list, npts is number of points on the contour and options are same as for ListPlot. See also PlotMinuitContourFixed and PlotMinuitContourFixedAll."

PlotMinuitContourFixed::usage = "PlotMinuitContourFixed[par1, par2, npts, options] plots two-parameter correlation contour produced by Minuit's MNCONT. par1 and par2 are symbols from Parameters list, npts is number of points on the contour and options are same as for ListPlot. Other parameters are temporarily fixed and then released for faster plotting. See also PlotMinuitContour and PlotMinuitContourFixedAll"

PlotMinuitContourFixedAll::usage = "PlotMinuitContourFixed[par1, par2, npts, options] plots two-parameter correlation contour produced by Minuit's MNCONT. par1 and par2 are symbols from Parameters list, npts is number of points on the contour and options are same as for ListPlot. Other parameters are temporarily fixed and then released for faster plotting. Additionally values of par1 and par2 are returned to old values before plotting. See also PlotMinuitContour and PlotMinuitContourFixed"
Begin["`Private`"]

AllParameterValues[] := Block[{iv = Table[0, {70}]}, 
(iv[[#[[1]]]] = #[[3]]) & /@ Parameters; iv]

SpliceToFortran[path_String, tfor_] := 
  Block[{}, Off[AssignFunction::undef];
    GPDQsplice = 
      FortranAssign[FCM[1], GPDMom[j, t, xi][[1]] /. tfor, 
        AssignOptimize -> False, AssignMaxSize -> Infinity, 
        AssignPrecision -> Infinity]; 
    GPDGsplice = 
      FortranAssign[FCM[2], GPDMom[j, t, xi][[2]] /. tfor, 
        AssignOptimize -> False, AssignMaxSize -> Infinity, 
        AssignPrecision -> Infinity];
    PDFQsplice = 
      FortranAssign[FCM[1], PDFMom[j, t, xi][[1]] /. tfor, 
        AssignOptimize -> False, AssignMaxSize -> Infinity, 
        AssignPrecision -> Infinity]; 
    PDFGsplice = 
      FortranAssign[FCM[2], PDFMom[j, t, xi][[2]] /. tfor, 
        AssignOptimize -> False, AssignMaxSize -> Infinity, 
        AssignPrecision -> Infinity];
    Splice[StringJoin[path, "/splice_template.f"], 
      StringJoin[path, "/splice.f"], FormatType -> OutputForm]]

lobj = Install["gepard.exe"]  (* Installing C and Fortran routines. *)

defaultopts = {SPEED -> -1, P -> -1, SCHEME -> "DFLT", ANSATZ -> "DFLT", 
  DATFILE -> "DFLT", OUTFILE -> "DFLT"};
Options[GepardInit] = defaultopts;
Options[cffH] = defaultopts;
Options[cffE] = defaultopts;
Options[F2] = defaultopts;

GepardInit[(opts___)?OptionQ] := GepardInitInternal @@ ( 
  {SPEED, P, SCHEME, ANSATZ, DATFILE, OUTFILE} 
        /. {opts} /. Options[GepardInit] )

cffH[(xi_)?NumericQ, (t_)?NumericQ, (q2_)?NumericQ, (q02_)?NumericQ, (opts___)?OptionQ] := cffHInternal @@ ( {xi, t, q2, q02,
		SPEED, P, SCHEME, ANSATZ} /. {opts} /. Options[cffH] )

cffE[(xi_)?NumericQ, (t_)?NumericQ, (q2_)?NumericQ, (q02_)?NumericQ, (opts___)?OptionQ] := cffEInternal @@ ( {xi, t, q2, q02,
		SPEED, P, SCHEME, ANSATZ} /. {opts} /. Options[cffE] )

F2[(xbj_)?NumericQ, (q2_)?NumericQ, (q02_)?NumericQ, (opts___)?OptionQ] := F2Internal @@ ( {xbj, q2, q02, SPEED, P, SCHEME, ANSATZ} /. {opts} /. Options[F2] )

GepardFit[pars_, (opts___)?OptionQ] := Block[{varpars = ParameterID /@ pars, 
      allpars = First[Transpose[Parameters]], status, ierr, chis}, 
      fixedpars = Complement[allpars, varpars]; 
      GepardInit[opts];
      jValues = MinuitInit[1]; (MinuitSetParameter @@ #1 & ) /@ Parameters; 
      MinuitCommand[StringJoin["fix ", StringJoin[
         Flatten[Table[{" ", ToString[fixedpars[[n]]]}, 
           {n, Length[fixedpars]}]]]]]; 
      Print[First[AbsoluteTiming[ierr = MinuitCommand["migrad"]]]]; 
      If[ierr == 4, 
       Print["Abnormal termination (e.g. MIGRAD not converged)"]; Abort[], 
       Null]; MinuitCommand["cali 3"]; status = MinuitStatus[]; 
      GPDcurrent[j_, t_, xi_] = GPDMom[j, t, xi] /. 
        PAR[n_] :> First[MinuitGetParameter[n]]; 
      Print[StringJoin["Quality of covariance matrix = ", 
        ToString[Last[status]]]]; 
      Print["  ----    Total and partial chi-squares  ----- "]; 
       chis = GetChiSquares[];
       Print[ChiTable[chis]];
      Print["  ----    Parameter status :       ----- "]; 
      ParameterStatus = Select[Table[Join[Parameters[[n]], 
          (MinuitGetParameter /@ Transpose[Parameters][[1]])[[n]]], 
         {n, Length[Parameters]}], Last[#1] != 0 & ]]

PrettyStatus[] := Block[{status, chis},
       MinuitCommand["cali 3"]; status = MinuitStatus[]; 
      GPDcurrent[j_, t_, xi_] = GPDMom[j, t, xi] /. 
        PAR[n_] :> First[MinuitGetParameter[n]]; 
      chis = GetChiSquares[];
      Print[StringJoin["Quality of covariance matrix = ", 
        ToString[Last[status]]]]; 
       Print[ChiTable[chis]];
      Print["  ----    Parameter status :       ----- "]; 
      ParameterStatus = Select[Table[Join[Parameters[[n]], 
          (MinuitGetParameter /@ Transpose[Parameters][[1]])[[n]]], 
         {n, Length[Parameters]}], Last[#1] != 0 & ]]
 
JustCali3[pars_, (opts___)?OptionQ] := Block[{varpars = ParameterID /@ pars, 
      allpars = First[Transpose[Parameters]], status, ierr}, 
     fixedpars = Complement[allpars, varpars]; 
      GepardInit[opts]; jValues = MinuitInit[1]; 
      (MinuitSetParameter @@ #1 & ) /@ Parameters; 
      MinuitCommand[StringJoin["fix ", StringJoin[
         Flatten[Table[{" ", ToString[fixedpars[[n]]]}, 
           {n, Length[fixedpars]}]]]]]; ierr = MinuitCommand["cali 3"]; 
      status = MinuitStatus[]; Print[StringJoin["\!\(\[Chi]\^2\) = ", 
        ToString[First[status]]]]; GPDcurrent[j_, t_, xi_] = 
       GPDMom[j, t, xi] /. PAR[n_] :> First[MinuitGetParameter[n]]; 
      ]

plotPDFs[] := LogLinearPlot[{GPD[1, x, 0, 0], GPD[2, x, 0, 0]}, 
     {x, 1/10000, 0.1}, PlotRange -> All, AxesLabel -> 
      {"x", "xf(x)"}, PlotStyle -> 
      {{Thickness[0.01], RGBColor[0, 0, 1]}, {Thickness[0.01], 
        RGBColor[1, 0, 0]}}]

plotslopes[] := LogLinearPlot[{slope[1, x], slope[2, x]}, {x, 1/10000, 0.01}, 
     PlotRange -> All, AxesLabel -> {"x", "B(x)"}, 
     PlotStyle -> {{Thickness[0.01], RGBColor[0, 0, 1]}, 
       {Thickness[0.01], RGBColor[1, 0, 0]}}]
 
slope[f_Integer, x_] := - ND[Log[GPD[f, x, -t, 0]], t, 0]
 
GPD[f_Integer, (x_)?NumericQ, (t_)?NumericQ, (xi_)?NumericQ] := 
    Chop[(1/(2*Pi))*NIntegrate[GPDcurrent[0.5 + I*y, t, xi][[f]]/
        x^(0.5 + I*y), {y, -10, -2, -0.5, 0, 0.5, 2, 10}]]
 
GPD[pars_, t_, xi_] := 
    Module[{j, upars = Join[{t, xi}, pars[[Transpose[Parameters][[1]]]]]}, 
     Flatten[Table[j = jValues[[jind]]; args = Join[{Re[j], Im[j]}, upars]; 
        {GPDMomCompiledQ @@ args, GPDMomCompiledG @@ args}, 
       {jind, Length[jValues]}]]]

PDF[pars_, t_, xi_] := 
    Module[{j, upars = Join[{t, xi}, pars[[Transpose[Parameters][[1]]]]]}, 
     Flatten[Table[j = jValues[[jind]]; args = Join[{Re[j], Im[j]}, upars]; 
        {PDFMomCompiledQ @@ args, PDFMomCompiledG @@ args}, 
       {jind, Length[jValues]}]]]
 
(*FIXME: PAR(18) and PAR(28) hard-wired here!! *)

CompileMoments[] := Block[{},
        GPDMomCompiledQ = 
            Compile[Evaluate[{rej, imj, t, xi}~Join~
                  Union[Cases[GPDMom[j, t, xi], PAR[n_], Infinity]~Join~{PAR[18], PAR[28]}]], 
              Evaluate[GPDMom[rej + I imj, t, xi][[1]]]];
        GPDMomCompiledG = 
            Compile[Evaluate[{rej, imj, t, xi}~Join~
                  Union[Cases[GPDMom[j, t, xi], PAR[n_], Infinity]~Join~{PAR[18], PAR[28]}]], 
              Evaluate[GPDMom[rej + I imj, t, xi][[2]]]];
        PDFMomCompiledQ = 
            Compile[Evaluate[{rej, imj, t, xi}~Join~
                  Union[Cases[GPDMom[j, t, xi], PAR[n_], Infinity]~Join~{PAR[18], PAR[28]}]], 
              Evaluate[PDFMom[rej + I imj, t, xi][[1]]]];
        PDFMomCompiledG = 
            Compile[Evaluate[{rej, imj, t, xi}~Join~
                  Union[Cases[GPDMom[j, t, xi], PAR[n_], Infinity]~Join~{PAR[18], PAR[28]}]], 
              Evaluate[PDFMom[rej + I imj, t, xi][[2]]]];
]


ParameterID[symb_Symbol] := Select[Parameters, #1[[2]] == ToString[symb] & ][[
     1,1]]

ParameterName[id_Integer] := ToString[Select[Parameters, #1[[1]] == id & ][[
     1,2]]]

ChiSquareProbability[d_, chisq_] := Gamma[d/2, chisq/2]/Gamma[d/2]
 
ReducedChiSquareProbability[d_, chisq_] := ChiSquareProbability[d, d chisq]

rdfile[fname_] := Block[{str, 
    outp = {}, ln = Null}, str = OpenRead[fname]; While[
        ln =!= EndOfFile, AppendTo[outp, ln = Read[str, String]]]; Close[str];
       outp]

PrintMinuitCommand[comm_?StringQ, fsize_:5] := 
  Module[{before, after, fname="fit.mnt"}, 
    before = rdfile[fname]; MinuitCommand[comm]; after = rdfile[fname]; 
    StylePrint[TableForm[Take[after, Length[before] - Length[after]]], 
      FontFamily->"Bitstream vera sans mono", 
      FontSize -> fsize, CellLabel -> fname]]

(* Nice printing any table *)
PrettyTable[m_List, opts___] := DisplayForm[StyleBox[
      GridBox[m, Evaluate[FilterOptions[GridBox, opts]], GridFrame -> 2, 
        RowLines -> {1, 0}, ColumnLines -> {1, 1, 0}], 
      Evaluate[FilterOptions[StyleBox, opts]], Background -> GrayLevel[0.9]]]

(* Nice printing chi-squares table *)
ChiTable[chis_List] := PrettyTable[
      (* Calculating chi-square probabilities and coloring bad ones *)
       chisa = {#[[1]], #[[2]], (prb = ChiSquareProbability[#[[2]], #[[1]]]; 
        prbpr = NumberForm[prb, {4, 3}]; 
        If[prb < 0.2, DisplayForm[StyleBox[prbpr, FontColor -> Hue[1]]], 
          prbpr])} & /@ chis;
        chis1=Transpose[chisa]; 
        chis2={Map[NumberForm[#, {5, 1}] &, chis1[[1]], {1}], chis1[[2]], chis1[[3]]};
    Transpose[Prepend[Transpose[Prepend[chis2, 
           (* labels of columns : *)
 {"    total     "}~Join~Table["part-"<>ToString[k], {k, 1, Length[First[chis2]] - 1}]
               ]], 
           (* labels of rows : *)
    {"", "\!\(\[Chi]\^2\)", "d.o.f", "probabilities"}
    ]]]

(* Plotting two-parameter correlation ellipse *)

PlotMinuitContour[par1_Symbol, par2_Symbol, npts_Integer, (opts___)?OptionQ] :=
  Module[{list},
    list = MinuitContour[ParameterID[par1], ParameterID[par2], npts]; 
    AppendTo[list, First[list]]; (* closing contour *)
    ListPlot[list, Join[
      {AxesLabel -> {par1,par2}, PlotJoined->True, 
      PlotStyle -> {Hue[1], Thickness[0.01]}}, {opts}]]
        ]

(* Same, but with all the other parameters fixed *)

PlotMinuitContourFixed[par1_Symbol, par2_Symbol, npts_Integer, (opts___)?OptionQ] :=
  Module[{allvarpars, contpars, strpars},
    contpars = {ParameterID[par1], ParameterID[par2]};
    allvarpars = Transpose[Select[Table[Join[Parameters[[n]], 
     (MinuitGetParameter /@ Transpose[Parameters][[1]])[[n]]], 
            {n, Length[Parameters]}], 
                Last[#1] != 0 &]
                       ][[1]];
    otherpars = Complement[allvarpars, contpars];
    strpars = StringJoin[(ToString[#] <> " ") & /@ otherpars];
    MinuitCommand["fix " <> strpars];
    PlotMinuitContour[par1, par2, npts, opts];
    MinuitCommand["release " <> strpars];
        ]

(* Same, but with all the other parameters fixed, and contour parameters   *)
(*  returned to previous value after plotting.                             *)

PlotMinuitContourFixedAll[par1_Symbol, par2_Symbol, npts_Integer, (opts___)?OptionQ] :=
  Module[{allvarpars, contpars, strpars, oldvals},
    contpars = {ParameterID[par1], ParameterID[par2]};
    oldvals = First[Transpose[MinuitGetParameter /@ contpars]];
    allvarpars = Transpose[Select[Table[Join[Parameters[[n]], 
     (MinuitGetParameter /@ Transpose[Parameters][[1]])[[n]]], 
            {n, Length[Parameters]}], 
                Last[#1] != 0 &]
                       ][[1]];
    otherpars = Complement[allvarpars, contpars];
    strpars = StringJoin[(ToString[#] <> " ") & /@ otherpars];
    MinuitCommand["fix " <> strpars];
    PlotMinuitContour[par1, par2, npts, opts];
    MinuitCommand["release " <> strpars];
    MinuitCommand["set param " <> ToString[contpars[[1]]] <> " " <> ToString[oldvals[[1]]]];
    MinuitCommand["set param " <> ToString[contpars[[2]]] <> " " <> ToString[oldvals[[2]]]];
        ]

SaveParameters[] := Block[{},
Parameters = Map[{#[[1]], #[[2]], 
  MinuitGetParameter[#[[1]]][[1]], #[[4]], #[[5]], #[[6]]}&, Parameters];]

End[ ]

EndPackage[ ]
