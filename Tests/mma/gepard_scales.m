(* $Id: gepard.m,v 1.5 2006-04-06 16:10:22+02 kkumer Exp kkumer $ *)

(*               gepard                     *)
(*               ======                     *)


(* Various comments and checks of these expressions are in notebook gepard.nb *)

(*
Print["gepard: calculating stuff related to GEneralized PARton Distributions Ver 0.8 (2006-03-24)"];
*)

BeginPackage["gepard`"]

(* 0. Exported functions *)

CFF::usage = "CFF[xi, Delta2, Q2, {Q20}}, {{FF, flavor}, ansatz}, {NF, approx, {scheme}}, {{c, accurracy, LAMBDA}, eval}]
returns Compton form factor for kinematical parameters xi, Delta^2, Q^2, (Q^2)0, 
where Q is \"caligraphic Q\".  FF should be one of {HJi, HtJI, EJi, EtJi}. 
This and ansatz are just labels denoting corresponding input GPD that has to be
provided as global function Fcm[].
Argument approx is one of {LO, NLO, NNLO}.
Argument eval is one of {GEP, GEPI, GEPfix}, standing for exact routines, 2dim j-r interpolation,
and j interpolation with R fixed. Before using the last two, one needs to create 
interpolation functions by calling interpol[] or interpolfix[].
Scheme must be bar2CS or MS (MS stands for MSbar and is available only up to NLO and
only for eval=GEP). 
Also flavor must be SP since only singlet case is implemented. 
{c, accurracy, LAMBDA<=10) are parameters of inverse Mellin integration. See: Fcm"

Fcm::usage = "Fcm[{j, xi, Delta2, Q2}, {{FF, flavor}, ansatz}, {appr, scheme, rf, rr}]
is conformal moment of GPD function FF and it has to be provided by user. FF is one of 
{HJi, HtJI, EJi, EtJi}.  Ansatz stands for set of parameters describing the shape.
Argument approx is one of {LO, NLO, NNLO}. rf = Q^2/mu_{f}^2, rr = Q^2/mu_{r}^2."

cDVCS::usage = "cDVCS[{j, rf, rr, R, {sector, flavor}, {nf, order, scheme}, {c, eval}]
is term of order <order> in expansion in astrong/(2 Pi) of
moment j of Wilson coefficient vector {C_Q, C_G} multiplied by evolution operator \\Epsilon(R),
with R=astrong(Q2)/atrong(Q02). rf is Q^2/mu_f^2, and rr is Q^2/mu_r^2. Argument
eval is one of {GEP, GEPI, GEPfix}, standing for exact routines, 2dim j-r interpolation,
and j interpolation with R fixed. Before using the last two, one needs to create 
interpolation functions by calling interpol[] or interpolfix[].
Scheme must be bar2CS or MS (MS stands for MSbar and is available only up to NLO and
only for eval=GEP). 
Also flavor must be SP since only singlet case is implemented.
\n
cDVCS[{j, rf, rr}, {V, SP}, {nf, order, scheme}, {GEP}] gives the same, but before
evolution, and cDVCS[{j, rf, rr}, {V, SP}, {nf, scheme}, {GEP}] gives a list with
all three orders=0,1,2. Note that here eval=GEP.
\n
This is a low level function and a call to it MUST BE PRECEDED BY A CALL TO
common[j+1, nf] to initialize common blocks!!
"

RAope::usage = " RAope[{j, R}, {V, SP}, {nf, scheme}, {GEP}]
is moment j of evolution operator \\Epsilon(R) with R=astrong(Q2)/atrong(Q02).
\n
This is a low level function and a call to it MUST BE PRECEDED BY A CALL TO
common[j+1, nf] to initialize common blocks!!
"

Fpw::usage = "Fpw[{j, xi, Delta2, Q2, {Q02}}, {{FF, flavor}, ansatz}, {nf, {appr, order}, 
{scheme, rf_:1, rr_:1}}, {c, eval}] is j-th partial wave of CFF. Arguments have same
meaning as for CFF."

interpol::usage = "interpol[{npoints, rf, rr, {Rstart, Rend, Rstep}}, {V, SP}, {NF, scheme}, {c_:0.5, LAM_:10}]
creates 2D interpolation fuction for evolved Wilson coefficients interpolating on 6 x npoints
points along inverse Mellin integration path (which are more dense around zero), and
for scale ratio R points defined by iterator {Rstart, Rend, Rstep}. rf is Q^2/mu_f^2, and
rr is Q^2/mu_r^2. {V, SP} has to be just that since only unpolarized (vector) singlet
case is implemented. scheme has to be bar2CS for the same reason. Inverse Mellin
integration goes from c-i*LAM to c+i*LAM. npoints=30, Rstep=0.05 is 10^-5 accurate"

interpolfix::usage = "interpolfix[{npoints, rf, rr, Q2, Q02}, {V, SP}, {NF, scheme}, {c_:0.5, LAM_:10}] 
creates interpolation fuction for evolved Wilson coefficients
interpolating on 6 x npoints
points along inverse Mellin integration path (which are more dense around zero), for
fixed scale ratio R=astrong(Q2)/atrong(Q02). rf is Q^2/mu_f^2, and
rr is Q^2/mu_r^2. {V, SP} has to be just that since only unpolarized (vector) singlet
case is implemented. scheme has to be bar2CS for the same reason. c is parameter
of inverse Mellin integration."


common::usage = "common[n, NF] initializes values in memory for DIS anomalous dimensions
and Wilson coefficients of complex moment n=j+1. This has to be called before calling
low-level functions, but CFF[] takes care of it itself."

AS2P::usage = "AS2P[mu2, NF, approx] returns alpha_strong(mu2)/(2 Pi) evolved from
fixed boundary condition alpha_strong(2.5)/Pi = 0.1 by direct integration of RGE.
approx is 0-3, corresponding to LO-NNNLO. Functions NP2O and O2NP may be useful."

NP2O::usage = "NP2O[LO] = 0, ..., NP2O[NNLO] = 2, ..."

O2NP::usage = "O2NP[0] = LO, ..., O2NP[2] = NNLO, ..."

V::usage = "label"
SP::usage = "label"
HJi::usage = "label"
toySP::usage = "label"
bar2CS::usage = "label"
MS::usage = "label"
LO::usage = "label"
NLO::usage = "label"
NNLO::usage = "label"
GEP::usage = "GEP is evaluation label, denoting evaluation using exact Fortran routines. See: CFF"
GEPI::usage = "GEPI is evaluation label, denoting evaluation using interpolation. See: CFF"
GEPIfix::usage = "GEPIfix is evaluation label, denoting evaluation using interpolation with
fixed scales. See: CFF"
CFFold::usage = "Old CFF: integration from -LAM to LAM."
CFFmon::usage = "CFF with Reap[ ... Sow[y]] for investigating integrand calls."
CFFphi::usage = "CFF with slanted integration contour for faster convergence"
CFFDM::usage = "CFF a la Dieter"

Begin["`Private`"]

Install["gepard.exe"] (* Installing Fortran routines. (Why are they not private?) *)


(* Transforming LO->0, NLO->1, ... and vice versa *)
NP2O[appr_] := Count[Characters[ToString[appr]], "N"]
O2NP[ord_] := ToExpression[StringJoin[Table["N", {ord}]] <> "LO"]

(* Transforming HJi, EJi -> V, and HtJi, EtJI -> A *)
GPD2Sector[gpd_] := Switch[gpd, HJi, V, EJi, V, HtJi, A, EtJi, A]

(* Choosing Tan for V and -Cot for A case *)
MinCot[z_] := -Cot[z]
GPD2Trigonometric[gpd_] := Switch[gpd, HJi, Tan, EJi, Tan, HtJi, MinCot, 
     EtJi, MinCot]
 

(* Complete j-th partial wave, up to some order appr=LO, NLO or NNLO *)
(* ----------------------------------------------------------------- *)

(* Calls input GPD Fcm[] which has to be already defined somewhere. *)

Fpw[{(j_)?NumericQ, (xi_)?NumericQ, Delta2_, (Q2_)?NumericQ, 
      {(Q02_)?NumericQ}}, {{F_, flavor_ /; flavor == SP}, argAns__}, 
      {nf_, LO, {scheme_, rf_:1, rr_:1}}, 
      {c___, eval_ /; eval == GEP || eval==GEPI || eval==GEPIfix}] := 
    Block[{R = AS2P[Q2/rf, nf, 0]/AS2P[Q02, nf, 0]}, 
     If[eval == GEP, common[j + 1, nf],Null]; 
    cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] . RAope[
     {j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] . 
     Fcm[{j, xi, Delta2, Q02}, {{F, flavor}, argAns}, {LO, scheme, rf, rr}]]

Fpw[{(j_)?NumericQ, (xi_)?NumericQ, Delta2_, (Q2_)?NumericQ, 
      {(Q02_)?NumericQ}}, {{F_, flavor_ /; flavor == SP}, argAns__}, 
      {nf_, NLO, {scheme_, rf_:1, rr_:1}}, 
      {c___, eval_ /; eval == GEP || eval==GEPI || eval==GEPIfix}] := 
    Block[{R = AS2P[Q2/rf, nf, 1]/AS2P[Q02, nf, 1]}, 
     If[eval == GEP, common[j + 1, nf],Null]; 
    (cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] . RAope[
     {j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] + 
 (* NLO *)
    AS2P[Q2/rr, nf, 1] cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, 
    {eval}][[2]] . RAope[{j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] +
     AS2P[Q2/rf, nf, 1] cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, 
    {eval}][[1]] . RAope[{j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[2]]) . 
  Fcm[{j, xi, Delta2, Q02}, {{F, flavor}, argAns}, {NLO, scheme, rf, rr}]
]

Fpw[{(j_)?NumericQ, (xi_)?NumericQ, Delta2_, (Q2_)?NumericQ, 
      {(Q02_)?NumericQ}}, {{F_, flavor_ /; flavor == SP}, argAns__}, 
      {nf_, NNLO, {scheme_, rf_:1, rr_:1}}, 
      {c___, eval_ /; eval == GEP || eval==GEPI || eval==GEPIfix}] := 
    Block[{R = AS2P[Q2/rf, nf, 2]/AS2P[Q02, nf, 2]}, 
     If[eval == GEP, common[j + 1, nf],Null]; 
  (  cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] . RAope[
     {j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] + 
 (* NLO *)
    AS2P[Q2/rr, nf, 2] cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, 
    {eval}][[2]] . RAope[{j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] +
     AS2P[Q2/rf, nf, 2] cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, 
    {eval}][[1]] . RAope[{j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[2]] +
 (* NNLO *)
    AS2P[Q2/rr, nf, 3]^2 cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, 
    {eval}][[3]] . RAope[{j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[1]] +
    AS2P[Q2/rr, nf, 3] AS2P[Q2/rf, nf, 3] cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, 
    {eval}][[2]] . RAope[{j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[2]] +
    AS2P[Q2/rf, nf, 3]^2 cDVCS[{j, rf, rr}, {GPD2Sector[F], flavor}, {nf, scheme}, 
    {eval}][[1]] . RAope[{j, R}, {GPD2Sector[F], flavor}, {nf, scheme}, {eval}][[3]] ) . 
  Fcm[{j, xi, Delta2, Q02}, {{F, flavor}, argAns}, {NNLO, scheme, rf, rr}]
]


(* Compton form factor *)
(* ------------------- *)
 
CFF[{(xi_)?NumericQ, (Delta2_)?NumericQ, (Q2_)?NumericQ, {(Q02_)?NumericQ}}, 
     {{F_ /; GPD2Sector[F] == V, flavor_}, 
      argAns__}, {argAppr__}, {{c_:0.5, acc_, LAM_:10}, 
      eval_ /; eval == GEP || eval == GEPI || eval == GEPIfix}] := 
     (xi/2)^(-N[c] - 1) *(
    NIntegrate[(xi/2)^(- I*y)*(I + GPD2Trigonometric[F] @@ 
         {Pi*((N[c] + I*y)/2)})*(Gamma[5/2 + N[c] + I*y]/
        Gamma[3 + N[c] + I*y])*Fpw[{N[c] + I*y, xi, Delta2, Q2, {Q02}}, 
        {{F, flavor}, argAns}, {argAppr}, {N[c], eval}], {y, -LAM, -0.5, 0.5, LAM}, 
      AccuracyGoal -> acc, MinRecursion -> 2, MaxRecursion -> 10]
          )/(2*Gamma[3/2])

CFFmon[{(xi_)?NumericQ, (Delta2_)?NumericQ, (Q2_)?NumericQ, {(Q02_)?NumericQ}}, 
     {{F_ /; GPD2Sector[F] == V, flavor_}, 
      argAns__}, {argAppr__}, {{c_:0.5, acc_, LAM_:10}, 
      eval_ /; eval == GEP || eval == GEPI || eval == GEPIfix}] := Reap[
     (xi/2)^(-N[c] - 1) *(
    NIntegrate[(xi/2)^(- I*y)*(I + GPD2Trigonometric[F] @@ 
         {Pi*((N[c] + I*y)/2)})*(Gamma[5/2 + N[c] + I*y]/
        Gamma[3 + N[c] + I*y])*Fpw[{N[c] + I*y, xi, Delta2, Q2, {Q02}}, 
        {{F, flavor}, argAns}, {argAppr}, {N[c], eval}], {y, -LAM, -0.5, 0.5, LAM}, 
      AccuracyGoal -> acc, MinRecursion -> 2, MaxRecursion -> 10, EvaluationMonitor:>Sow[y]]
          )/(2*Gamma[3/2])][[2,1]]

CFFold[{(xi_)?NumericQ, (Delta2_)?NumericQ, (Q2_)?NumericQ, {(Q02_)?NumericQ}}, 
     {{F_ /; GPD2Sector[F] == V, flavor_}, 
      argAns__}, {argAppr__}, {{c_:0.5, acc_, LAM_:10}, 
      eval_ /; eval == GEP || eval == GEPI || eval == GEPIfix}] := 
     (xi/2)^(-N[c] - 1) *(
    NIntegrate[(xi/2)^(- I*y)*(I + GPD2Trigonometric[F] @@ 
         {Pi*((N[c] + I*y)/2)})*(Gamma[5/2 + N[c] + I*y]/
        Gamma[3 + N[c] + I*y])*Fpw[{N[c] + I*y, xi, Delta2, Q2, {Q02}}, 
        {{F, flavor}, argAns}, {argAppr}, {N[c], eval}], {y, -LAM, LAM}, 
      AccuracyGoal -> acc, MinRecursion -> 2, MaxRecursion -> 10]
          )/
     (2*Gamma[3/2])


CFFDM[{(xi_)?NumericQ, (Delta2_)?NumericQ, (Q2_)?NumericQ, {(Q02_)?NumericQ}}, 
     {{F_ /; GPD2Sector[F] == V, flavor_}, 
      argAns__}, {argAppr__}, {{c_:0.5, acc_, LAM_:10}, 
      eval_ /; eval == GEP || eval == GEPI || eval == GEPIfix}] := 
     (xi/2)^(-N[c] - 1) *(
    I* NIntegrate[Re[(xi/2)^(- I*y)*(Gamma[5/2 + N[c] + I*y]/
        Gamma[3 + N[c] + I*y])*Fpw[{N[c] + I*y, xi, Delta2, Q2, {Q02}}, 
        {{F, flavor}, argAns}, {argAppr}, {N[c], eval}]], {y, 0, 10, LAM}, 
      AccuracyGoal -> acc, MinRecursion -> 3, MaxRecursion -> 10] +
    NIntegrate[Re[(xi/2)^(- I*y)*(GPD2Trigonometric[F] @@ 
         {Pi*((N[c] + I*y)/2)})*(Gamma[5/2 + N[c] + I*y]/
        Gamma[3 + N[c] + I*y])*Fpw[{N[c] + I*y, xi, Delta2, Q2, {Q02}}, 
        {{F, flavor}, argAns}, {argAppr}, {N[c], eval}]], {y, 0, 0.5, 2, 10, LAM}, 
      AccuracyGoal -> acc, MinRecursion -> 2, MaxRecursion -> 10]
      )/(Gamma[3/2])

CFFphi[{(xi_)?NumericQ, (Delta2_)?NumericQ, (Q2_)?NumericQ, {(Q02_)?NumericQ}}, 
     {{F_ /; GPD2Sector[F] == V, flavor_}, 
      argAns__}, {argAppr__}, {{c_:0.5, acc_, LAM_:10}, 
      eval_ /; eval == GEP || eval == GEPI || eval == GEPIfix}] := Block[{eph=N[Exp[3 Pi I/4]]},
     (xi/2)^(-N[c] - 1) *(
    I* Im[NIntegrate[eph*(xi/2)^(- eph*y)*(Gamma[5/2 + N[c] + eph*y]/
        Gamma[3 + N[c] + eph*y])*Fpw[{N[c] + eph*y, xi, Delta2, Q2, {Q02}}, 
        {{F, flavor}, argAns}, {argAppr}, {N[c], eval}], {y, 0, 0.5, LAM}, 
      AccuracyGoal -> acc, MinRecursion -> 3, MaxRecursion -> 10]] +
    Im[NIntegrate[eph*(xi/2)^(- eph*y)*(GPD2Trigonometric[F] @@ 
         {Pi*((N[c] + eph*y)/2)})*(Gamma[5/2 + N[c] + eph*y]/
        Gamma[3 + N[c] + eph*y])*Fpw[{N[c] + eph*y, xi, Delta2, Q2, {Q02}}, 
        {{F, flavor}, argAns}, {argAppr}, {N[c], eval}], {y, 0, 0.5, LAM}, 
      AccuracyGoal -> acc, MinRecursion -> 2, MaxRecursion -> 10]]
      )/(Gamma[3/2]) ]


(* Interpolations of evolved Wilson coefficients *)
(* --------------------------------------------- *)

(* Following function creates array with cca (6 x np) points, divided in six areas: 
(-LAM,-1), (-1,-0.1), (-0.1,0), (0,0.1), (0.1, 1) and (1,LAM), each with equidistant np
points. Thus region around zero, where integrand is largest and changes itself 
more violently, is interpolated more accurately. *)

interpts[np_, LAM_:10.] := (
  Range[-LAM, -1.1, 10./np]~Join~Range[-1., -0.11, 1./np]~Join~
    Range[-0.1, 0.099, 0.1/np]~Join~Range[0.1, 0.99, 1./np]~Join~
    Range[1., LAM, 10./np] )

(* interpts[np_] := Join[Join[Join[Join[Range[-10., -1.1, 10./np], 
        Range[-1., -0.11, 1./np]], Range[-0.1, 0.099, 0.1/np]], 
      Range[0.1, 0.99, 1./np]], Range[1., 10., 10./np]] *)

(* 2D interpolation in j-R space *)
interpol[{points_, rf_, rr_, {Rstart_, Rend_, Rstep_}}, {V, SP}, 
     {nf_, scheme_}, {c_:0.5, LAM_:10.}] := Module[{aux, tmpQ, tmpG}, 
WriteString["stdout", "Interpolating ... "];
      Do[WriteString["stdout", O2NP[ord], " ... "]; 
        (*  Making list  {{j1, R1, CQ(j1), CG(j1)}, {j1, R2, ...} ... } *)
       aux = Flatten[Outer[
               Flatten[{#1, #2, common[c + I*#1 + 1, nf]; 
                        cDVCS[{c + I*#1, 1, 1, #2}, {V, SP}, {3, ord, bar2CS}, 
                                                {GEP}]}]& , 
                     interpts[points, LAM], Range[Rstart, Rend, Rstep]], 1];
        (* resorting above list into {{j1, R1, CQ(j1)}, {j1, R2,...} and interpolating *)
          tmpQ = Interpolation[Transpose[Drop[Transpose[aux], {4}]]]; 
        (* resorting above list into {{j1, R1, CG(j1)}, {j1, R2,...} and interpolating *)
          tmpG = Interpolation[Transpose[Drop[Transpose[aux], {3}]]]; 
        cDVCS[{(j_)?NumericQ, rf, rr, (R_)?NumericQ}, {V, SP}, 
          {nf, ord, bar2CS}, {N[c], GEPI}] = {tmpQ[Im[j], R], tmpG[Im[j], R]}, 
       {ord, 0, 2}]; 
WriteString["stdout", "Done."]]

 
(* interpolatin with fixed R, faster *)
interpolfix[{points_, rf_, rr_, Q2_, Q02_}, {V, SP}, {nf_, scheme_}, 
     {c_:0.5, LAM_:10.}] := Module[{aux, tmpQ, tmpG, R}, 
Print["Interpolating for fixed ratio: astrong[", Q2, "]/astrong[", Q02, "] ..."]; 
     Do[R = AS2P[Q2, nf, nappr]/AS2P[Q02, nf, nappr]; 
        WriteString["stdout", " ", O2NP[nappr], ":"]; 
        Do[WriteString["stdout", ".[", ord, "]."]; 
        (*  Making list  {{j1, CQ(j1), CG(j1)}}, {j2,...} ... } *)
          aux = Flatten[
            {#1, common[c + I*#1 + 1, nf]; 
             cDVCS[{c + I*#1, 1, 1, R}, {V, SP}, {3, ord, bar2CS}, {GEP}]}
                       ]&  /@ interpts[points, LAM]; 
        (* resorting above list into {{j1, CQ(j1)}, ...} and interpolating *)
          tmpQ = Interpolation[Transpose[Drop[Transpose[aux], {3}]]]; 
        (* resorting above list into {{j1, CG(j1)}, ...} and interpolating *)
          tmpG = Interpolation[Transpose[Drop[Transpose[aux], {2}]]]; 
          cDVCS[{(j_)?NumericQ, rf, rr, R}, {V, SP}, {nf, ord, bar2CS}, 
            {N[c], GEPIfix}] = {tmpQ[Im[j]], tmpG[Im[j]]}, {ord, 0, nappr}], 
       {nappr, 0, 2}];
WriteString["stdout", "... Done."]]

End[ ]

EndPackage[ ]
