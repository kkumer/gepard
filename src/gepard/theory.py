"""Definition of theory frameworks."""

# from joblib import Parallel, delayed
from numpy import array, ndarray, sqrt
from scipy.stats import scoreatpercentile

from . import data, eff, gpd, model, quadrature

NCPU = 23  # how many CPUs to use in parallel


class Theory(object):
    """Class of theory frameworks for calculation of observables.

    This is a base class for a complete theory framework.
    It implements the generic routines for calculation of
    observables and uncertainty propagation.
    """

    def __init__(self, **kwargs) -> None:
        """This init is guaranteed to run.

        Args:
            name: short unique model name
            texname: TeX model name for (e.g. for plot annotation)
            description: longer description of the model
        """
        self.name = kwargs.setdefault('name', 'N/A')
        self.texname = kwargs.setdefault('texname', self.name)
        self.description = kwargs.setdefault('description', 'N/A')
        self.model = self       # to make old .m code work
        self.m = self       # to make old .m code work
        print('theory.Theory init done.')
        # We do NOT call super().__init__ here because object class will
        # not accapt **kwargs. This means that Theory has to be the last
        # class in mro.

    def chisq_single(self, points: data.DataSet, asym: bool = False,
                     **kwargs) -> float:
        """Return total chi-square.

        Args:
            points: measurements with uncertainties, observable is
                named in `yaxis` attribute of each DataPoint,
                value is `val` and uncertainty is `err`.
            asym: if measurements provide asymmetric uncertainties
                `errplus` and `errminus`, this enables their usage
            observable (str): overrides `yaxis` of DataPoints

        Notes:
            If the theory or model provide uncertainties, they are ignored -
            only experimental uncertainties are taken into account.
           """
        allpulls = []
        for pt in points:
            diff = (self.predict(pt, observable=pt.yaxis, **kwargs) - pt.val)
            if asym:
                if diff > 0:
                    allpulls.append(diff/pt.errplus)
                else:
                    allpulls.append(diff/pt.errminus)
            else:
                allpulls.append(diff/pt.err)
        chi = sum(p*p for p in allpulls)  # equal to m.fval if minuit fit is done
        return chi

    def pull(self, pt: data.DataPoint):
        """Return pull of a single Datapoint."""
        return (self.predict(pt, observable=pt.yaxis) - pt.val) / pt.err

#     def chisq_para(self, points: DataSet, asym: bool = False,
#                    **kwargs) -> float:
#         """Return total chi-square - parallel version.
# 
#         Warning:
#             Cannot be used until underlying global Fortran variables can change
#             during session. (Like different kinematics of data points.)
#         """
#         allpulls = Parallel(n_jobs=NCPU)(delayed(self.pull)(pt) for pt in points)
#         chi = sum(p*p for p in allpulls)  # equal to m.fval if minuit fit is done
#         return chi

    chisq = chisq_single

    def predict(self, pt, error=False, CL=False, **kwargs):
        """Give prediction for DataPoint pt.

        Args:
            pt: instance of DataPoint
            error: if available, produce tuple (mean, error)
            CL: (NNet only) error is not std.dev.,
                but 68% C.L. (mean, errplus, errminus)
            observable: string. Default is pt.yaxis. It is acceptable also
                        to pass CFF as observable, e.g., observable = 'ImH'
            parameters: dictionary which will temporarily update model's one
            orig_conventions: give prediction using original conventions of
                              the given DataPoint (e.g. for plotting)
        """
        if 'observable' in kwargs:
            obs = kwargs['observable']
        else:
            obs = pt.yaxis

        if 'parameters' in kwargs:
            old = m.parameters.copy()
            self.parameters.update(kwargs['parameters'])
        #elif isinstance(m, Model.ComptonNeuralNets):
        #    # It is not training (which always uses 'parameters'), and
        #    # we are not asked for particular net (call would again come
        #    # with 'parameters'), so we want mean of all nets
        #    self.parameters['nnet'] = 'ALL'
        #    result = getattr(self, obs)(pt)
        #    if error:
        #        return (result.mean(), result.std())
        #    else:
        #        return result.mean()

        # if obs in m.allCFFs or obs in m.allGPDs:
            # # we want GPD/CFF, which are model attributes
            # fun = getattr(m, obs)
        # else:
            # # we want a "real" observable, which is theory attribute
            # fun = getattr(self, obs)

        fun = getattr(self, obs)

        if error:
            try:
                # We now do standard propagation of error from m to observable
                pars = [p for p in self.parameters if not self.parameters_fix['p']]
                var = 0
                dfdp = {}
                for p in pars:
                    # calculating dfdp = derivative of observable w.r.t. parameter:
                    h=sqrt(self.covariance[p,p])
                    mem = self.parameters[p]
                    self.parameters[p] = mem+h/2.
                    up = fun(pt)
                    self.parameters[p] = mem-h/2.
                    down = fun(pt)
                    self.parameters[p] = mem
                    dfdp[p] = (up-down)/h
                for p1 in pars:
                    for p2 in pars:
                        var += dfdp[p1]*self.covariance[p1,p2]*dfdp[p2]
                result = (fun(pt), sqrt(var))
            except KeyError:
                # we have neural net
                allnets = fun(pt)
                if CL:
                    # 68% confidence level
                    m = allnets.mean()
                    try:
                        result = (m,
                              scoreatpercentile(allnets, 84)[0] - m,
                              m - scoreatpercentile(allnets, 16)[0])
                    except IndexError:
                        result = (m,
                              scoreatpercentile(allnets, 84) - m,
                              m - scoreatpercentile(allnets, 16))
                else:
                    # one sigma
                    result = (allnets.mean(), allnets.std())
        else:
            result = fun(pt)
            if isinstance(result, ndarray):
                # we have neural net
                result = result.mean()

        if 'parameters' in kwargs:
            # restore old values
            self.parameters.update(old)

        if kwargs.pop('orig_conventions', False):
            # express result in conventions of original datapoint
            try:
                result = (pt.orig_conventions(result[0]),) + result[1:]
            except (IndexError, TypeError):
                result = pt.orig_conventions(result)
        return result


    def _Xt4int(self, t, pt):
        """Same as _XDVCSt/_Xrhot but with additional variable t
        to facilitate integration over it.

        """
        aux = []
        for t_single in t:
            pt.t = t_single
            if hasattr(pt, 'process') and pt.process == 'gammastarp2rho0p':
                res = self._Xrhot(pt)
            else:
                res = self._XDVCSt(pt)
            del pt.t
            aux.append(res)
        return array(aux)


    def X(self, pt):
        """Total DVCS or DVMP cross section. """
        if 't' in pt or 'tm' in pt:
            # partial XS w.r.t momentum transfer t
            if hasattr(pt, 'process') and pt.process == 'gammastarp2rho0p':
                return self._Xrhot(pt)
            else:
                return self._XDVCSt(pt)

        else:
            # total XS
            if 'tmmax' in pt:
                tmmax = pt.tmmax
            else:
                tmmax = 1.  # default -t cuttoff in GeV^2
            res = quadrature.tquadrature(lambda t: self._Xt4int(t, pt), -tmmax, 0)
            return res
