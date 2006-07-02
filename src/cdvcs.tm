:Evaluate:      cDVCS[{j_?NumericQ, rf_?NumericQ, rr_?NumericQ}, {V, SP}, {nf_?NumericQ, bar2CS}, {GEP}] := cdvcs[N[nf], N[j], N[rf], N[rr]]
:Evaluate:      cDVCS[{j_?NumericQ, rf_, rr_}, {V, SP}, {nf_, ord_, bar2CS}, {GEP}] := cdvcs[N[nf], N[j], N[rf], N[rr]][[ord+1]]
:Begin:
:Function:      cdvcs
:Pattern:       cdvcs[nf_?NumericQ, j_?NumericQ, rf2_?NumericQ, rr2_?NumericQ]
:Arguments:     {nf, j, rf2, rr2}
:ArgumentTypes: {Manual}
:ReturnType:    Manual
:End:

