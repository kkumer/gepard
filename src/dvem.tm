:Evaluate: int2f1::usage = "int2f1[mu, j, nu, k, ind, acc] does integral f(x) A^mu,j B_nu,k. ind should be 1-5 for f(x)=1, (2-z)z, (1-z)/z, z-2, z, respectively. acc=3-6 is accuracy. "
:Evaluate: int2f1::nan = "Error - Argument not in correct numerical form."
:Evaluate: errnan := Message[int2f1::nan]
:Begin:
:Function:      int2f1
:Pattern:       int2f1[mu_?NumericQ, j_?NumericQ, nu_?NumericQ, k_?NumericQ, ind_Integer, acc_Integer]
:Arguments:     {mu, j, nu, k, ind, acc}
:ArgumentTypes: {Manual}
:ReturnType:    Manual
:End:

:Evaluate: cdvem::usage = "cdvem[j, k, rgpdf2, rdaf2, rr2] returns three (unnormalized) DVEM Wilson coefficients C_jk: quark, pure singlet quark and gluon. rgpdf2, rdaf2 and rr2 are ratios of Q2 and squares of GPD factorization, DA factorization and renormalization scales squared, respectively"
:Evaluate: cdvem::nan = "Error - Argument not in correct numerical form."
:Evaluate: errnan := Message[cdvem::nan]
:Begin:
:Function:      cdvem
:Pattern:       cdvem[j_?NumericQ, k_?NumericQ, rgpdf2_?NumericQ, rdaf2_?NumericQ, rr2_?NumericQ]
:Arguments:     {j, k, rgpdf2, rdaf2, rr2}
:ArgumentTypes: {Manual}
:ReturnType:    Manual
:End:
