:Evaluate: MinuitInit::usage = "MinuitInit[a] performs initialization. a is integer which is presently ignored. It returns lists with tValues and jValues so it should be called like {tValues, jValues} = MinuitInit[1]"
:Begin:
:Function:      MinuitInit
:Pattern:       MinuitInit[a_Integer]
:Arguments:     {a}
:ArgumentTypes: {Integer}
:ReturnType:    Manual
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

:Evaluate: cffH::usage = "cffH[xi, t, q2, q02, nf, P] is singlet CFF H(xi, t, q2, q02, nf, P). SPEED, SCHEME and ANSATZ have to be specified in GEPARD.INI"
:Begin:
:Function:      cffH
:Pattern:       cffH[(xi_)?NumericQ, (t_)?NumericQ, (q2_)?NumericQ, (q02_)?NumericQ, nf_Integer, p_Integer]
:Arguments:     {xi, t, q2, q02, nf, p}
:ArgumentTypes: {Manual}
:ReturnType:    Manual
:End:
