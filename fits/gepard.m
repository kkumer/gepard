(* $Id:*)

(*     GeParD - Mathematica interface    *)
(*     ==============================    *)


Print["GeParD - Mathematica interface (2008-01-15)"];


BeginPackage["gepard`", "Format`", "NumericalMath`NLimit`", "Graphics`Graphics`"]


GepardInit::usage = "GepardInit[] performs initialization of GeParD parameters. 
It reads GEPARD.INI for default parameters. Parameters SPEED, P, SCHEME and ANSATZ 
can then be overriden by specifying corresponding options e.g. GepardInit[P->0, ANSATZ->\"MMA\"]."

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

ParameterID::usage = "ParameterID[symbol] returns parameter number (id) of a parameter with
name symbol, as specified by array Parameters."

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
MinuitCommand::usage = "MinuitCommand[command_] executes MINUIT command"
MinuitStatus::usage = "MinuitStatus"
GetChiSquares::usage = "GetChiSquares"
PrintMinuitCommand::usage = "MinuitCommand[command_, file_] executes MINUIT
command and reads MINUIT output from file and prints it. 
MinuitCommand[command_, file_, fontsize_] does the printing with font
size fontsize (values 4-10 are nice)."

SPEED::usage = "GeParD parameter"
P::usage = "GeParD parameter"
SCHEME::usage = "GeParD parameter"
ANSATZ::usage = "GeParD parameter"

GPDQsplice::usage = "Fortran format function to be spliced into splice.f"
GPDGsplice::usage = "Fortran format function to be spliced into splice.f"
PDFQsplice::usage = "Fortran format function to be spliced into splice.f"
PDFGsplice::usage = "Fortran format function to be spliced into splice.f"

FCM::usage = "Fortran format array holding the expressions for GPDs"

j::usage = "j - complex conformal moment"
t::usage = "t - Mandelstam variable. In Fortran represented as DEL2"
xi::usage = "xi - DVCS kinematical variable."

lobj::usage = "lobj - MathLink link"

cffH::usage = "cffH[xi, t, q2, q02, options] returns singlet CFF H(xi, t, q2, q02). 
Other parameters are read from  GEPARD.INI and additionally parameters SPEED, P, SCHEME and ANSATZ 
can be set by specifying corresponding options e.g. GepardInit[P->0, ANSATZ->\"MMA\"]."

cffE::usage = "cffE[xi, t, q2, q02, options] returns singlet CFF E(xi, t, q2, q02). 
Other parameters are read from  GEPARD.INI and additionally parameters SPEED, P, SCHEME and ANSATZ 
can be set by specifying corresponding options e.g. GepardInit[P->0, ANSATZ->\"MMA\"]."

cffHInternal::usage = "MathLink function ..."
cffEInternal::usage = "MathLink function ..."
cInt::usage = "MathLink function ..."

AllParameterValues::usage = "AllParameterValues[] returns complete list that corresponds to 
COMMON block PAR. For internal use."

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

Options[GepardInit] = {SPEED -> -1, P -> -1, SCHEME -> "DFLT", ANSATZ -> "DFLT"};
Options[cffH] = {SPEED -> -1, P -> -1, SCHEME -> "DFLT", ANSATZ -> "DFLT"};
Options[cffE] = {SPEED -> -1, P -> -1, SCHEME -> "DFLT", ANSATZ -> "DFLT"};

GepardInit[(opts___)?OptionQ] := GepardInitInternal @@ ( {SPEED, P, SCHEME, ANSATZ} 
        /. {opts} /. Options[GepardInit] )

cffH[(xi_)?NumericQ, (t_)?NumericQ, (q2_)?NumericQ, (q02_)?NumericQ, (opts___)?OptionQ] := cffHInternal @@ ( {xi, t, q2, q02,
		SPEED, P, SCHEME, ANSATZ} /. {opts} /. Options[cffH] )

cffE[(xi_)?NumericQ, (t_)?NumericQ, (q2_)?NumericQ, (q02_)?NumericQ, (opts___)?OptionQ] := cffEInternal @@ ( {xi, t, q2, q02,
		SPEED, P, SCHEME, ANSATZ} /. {opts} /. Options[cffE] )

GepardFit[pars_, (opts___)?OptionQ] := Block[{varpars = ParameterID /@ pars, 
      allpars = First[Transpose[Parameters]], status, ierr, chis, dof, 
      probchi}, fixedpars = Complement[allpars, varpars]; 
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
      chis = GetChiSquares[];
      dof = chis[[1,2]]; probchi = ChiSquareProbability[dof, 
        chis[[1,1]]];  
      Print[StringJoin["quality of covariance matrix = ", 
        ToString[Last[status]]]]; 
      Print[Transpose[chis]];
      Print[StringJoin["Probability of this and larger total chi-square = ", 
        ToString[probchi]]]; 
      Print["  ----    Parameter status :       ----- "]; 
      ParameterStatus = Select[Table[Join[Parameters[[n]], 
          (MinuitGetParameter /@ Transpose[Parameters][[1]])[[n]]], 
         {n, Length[Parameters]}], Last[#1] != 0 & ]]

PrettyStatus[] := Block[{tchis},
       MinuitCommand["cali 3"]; status = MinuitStatus[]; 
      GPDcurrent[j_, t_, xi_] = GPDMom[j, t, xi] /. 
        PAR[n_] :> First[MinuitGetParameter[n]]; 
      chis = GetChiSquares[];
      dof = chis[[1,2]]; probchi = ChiSquareProbability[dof, 
        First[status]]; 
      Print[StringJoin["quality of covariance matrix = ", 
        ToString[Last[status]]]]; 
      Print[Transpose[chis]];
      Print[StringJoin["Probability of this and larger total chi-square = ", 
        ToString[probchi]]]; 
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
      Print["  ----    Parameter status :       ----- "]; 
      ParameterStatus = Select[Table[Join[Parameters[[n]], 
          (MinuitGetParameter /@ Transpose[Parameters][[1]])[[n]]], 
         {n, Length[Parameters]}], Last[#1] != 0 & ]]

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

ChiSquareProbability[d_, chisq_] := Gamma[d/2, chisq/2]/Gamma[d/2]
 
ReducedChiSquareProbability[d_, chisq_] := ChiSquareProbability[d, d chisq]

rdfile[fname_] := Block[{str, 
    outp = {}, ln = Null}, str = OpenRead[fname]; While[
        ln =!= EndOfFile, AppendTo[outp, ln = Read[str, String]]]; Close[str];
       outp]

PrintMinuitCommand[comm_?StringQ, fname_?StringQ, 
    fsize_:8] := Module[{before, after}, before = rdfile[fname]; 
  MinuitCommand[
    comm]; after = rdfile[fname]; StylePrint[TableForm[
      Take[after, Length[before] - 
    Length[after]]], "Input", FontSize -> fsize, CellLabel -> fname]]

End[ ]

EndPackage[ ]
