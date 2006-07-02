
:Evaluate:      ERfunc[{GEP}, nf_, {V, m_, R_}, j_?NumericQ] := erfunc[N[nf], N[j+1], N[R]]
:Begin:
:Function:      erfunc
:Pattern:       erfunc[nf_?NumericQ /; 2<nf<7, n_?NumericQ, r_?NumericQ]
:Arguments:     {nf, n, r}
:ArgumentTypes: {Manual}
:ReturnType:    Manual
:End:

