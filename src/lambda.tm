
:Evaluate:      Wlambda[{G}, nf_, {V, sign_}, j_] := lo[nf, j+1][[Round[-sign/2+1.5]]]
:Begin:
:Function:      lambda
:Pattern:       lambda[nf_?IntegerQ /; 2<nf<7, n_?NumericQ]
:Arguments:     {nf, n}
:ArgumentTypes: {Manual}
:ReturnType:    Manual
:End:

