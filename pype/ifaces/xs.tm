:Evaluate: XSunp::usage = "XSunp[ModelID, charge, polarization, Ee, Ep, xB, Q2, t, phi] calls python routines and returns the value for cross section of lepton of energy Ee on unpolarized proton target of energy Ep. Charge=-1 is for electron. ModelID is one of  0 (debug, always returns 42), 1 (NPB fit without Hall A), 2 (NPB fit with Hall A), 3 (preliminary hybrid fit with gepard sea, from Trento presentation). "
:Begin:
:Function:      XSunp
:Pattern:       XSunp[id_Integer, Q_Integer, lam_Integer, (Ee_)?NumericQ, (Ep_)?NumericQ, (xB_)?NumericQ, (Q2_)?NumericQ, (t_)?NumericQ, (phi_)?NumericQ]
:Arguments:     {id, Q, lam, Ee, Ep, xB, Q2, t, phi}
:ArgumentTypes: {Integer, Integer, Integer, Real, Real, Real, Real, Real, Real}
:ReturnType:    Manual
:End:
