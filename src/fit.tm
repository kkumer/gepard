:Evaluate: MinuitInit::usage = "MinuitInit[a] performs initialization. a is integer which is presently ignored. It returns lists with tValues and jValues so it should be called like {tValues, jValues} = MinuitInit[1]"
:Begin:
:Function:      MinuitInit
:Pattern:       MinuitInit[a_Integer]
:Arguments:     {a}
:ArgumentTypes: {Integer}
:ReturnType:    Manual
:End:

:Evaluate: GepardInitInternal::usage = "GepardInitInternal[SPEED, P, SCHEME, ANSATZ] performs initialization of GeParD parameters. It reads GEPARD.INI for default parameters which will then be overridden by the arguments with same name if they are positive numbers or strings different then 'DFLT'. This function is not ment to be called by user who should use GepardInit[] instead."
:Begin:
:Function:      GepardInitInternal
:Pattern:       GepardInitInternal[speed_Integer, p_Integer, scheme_String, ansatz_String]
:Arguments:     {speed, p, scheme, ansatz}
:ArgumentTypes: {Integer, Integer, String, String}
:ReturnType:    String
:End:


:Evaluate: MinuitSetParameter::usage = "MinuitSetParameter[id, pname, vstart, step, lo, hi] declares MINUIT parameter by forwarding arguments to subroutine MNPARM(id, pname, vstart, step, lo, hi, iflag)."
:Begin:
:Function:      MinuitSetParameter
:Pattern:       MinuitSetParameter[id_Integer, pname_String, (vstart_)?NumericQ, (step_)?NumericQ, (lo_)?NumericQ, (hi_)?NumericQ]
:Arguments:     {id, pname, vstart, step, lo, hi}
:ArgumentTypes: {Integer, String, Real, Real, Real, Real}
:ReturnType:    Integer
:End:

:Evaluate: MinuitCommand::usage = "MinuitCommand[string] sends to MINUIT the command 'string' by forwarding it to subroutine MNCOMD."
:Begin:
:Function:      MinuitCommand
:Pattern:       MinuitCommand[cmd_String]
:Arguments:     {cmd}
:ArgumentTypes: {String}
:ReturnType:    Integer
:End:

:Evaluate: MinuitStatus::usage = "MinuitStatus[] returns current status of minimization in form of a list {chi-square, quality of covariance matrix: 0-3}."
:Begin:
:Function:      MinuitStatus
:Pattern:       MinuitStatus[]
:Arguments:     {}
:ArgumentTypes: {}
:ReturnType:    Manual
:End:

:Evaluate: MinuitGetParameter::usage = "MinuitGetParameter[id] returns the current value of a parameter number id by calling Minuit subroutine MNPOUT. It returns list {value, error, internal parameter number (or zero if constant, or negative if undefined)."
:Begin:
:Function:      MinuitGetParameter
:Pattern:       MinuitGetParameter[id_Integer]
:Arguments:     {id}
:ArgumentTypes: {Integer}
:ReturnType:    Manual
:End:

:Evaluate: MinuitCovarianceMatrix::usage = "MinuitCovarianceMatrix[npar] returns the current value of (npar x npar) covariance matrix by calling Minuit subroutine MNEMAT. npar should be equal to the number of variable parameters"
:Begin:
:Function:      MinuitCovarianceMatrix
:Pattern:       MinuitCovarianceMatrix[npar_Integer]
:Arguments:     {npar}
:ArgumentTypes: {Integer}
:ReturnType:    Manual
:End:


:Evaluate: cffHInternal::usage = "cffHInternal[xi, t, q2, q02, SPEED, P, SCHEME, ANSATZ] returns singlet CFF H(xi, t, q2, q02). This function is not ment to be called by user who should use cffH[] instead."
:Begin:
:Function:      cffHInternal
:Pattern:       cffHInternal[(xi_)?NumericQ, (t_)?NumericQ, (q2_)?NumericQ, (q02_)?NumericQ, speed_Integer, p_Integer, scheme_String, ansatz_String]
:Arguments:     {xi, t, q2, q02, speed, p, scheme, ansatz}
:ArgumentTypes: {Real, Real, Real, Real, Integer, Integer, String, String}
:ReturnType:    Manual
:End:
