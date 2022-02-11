"""Representation of experimental data.

DataPoint -- class for kinematical points, possibly representing measurements
DataSet   -- container for DataPoint instances

dset -- dictionary with public datasets
"""
from __future__ import annotations

import copy
import importlib_resources
import math
import os
import re
from typing import Union, List

import pandas as pd

from . import kinematics
from .constants import Mp, Mp2
from .datasets import (DIS, en2engamma, ep2epgamma, gammastarp2gammap,
                       gammastarp2Mp)

process_class_map = {
        'ep2epgamma': 'DVCS',
        'en2engamma': 'DVCS',
        'gammastarp2gammap': 'DVCS',
        'gammastarp2rho0p': 'DVMP',
        'gammastarp2phip': 'DVMP',
        'dis': 'DIS'
        }


def loaddata(resource):
    """Return dictionary {id : DataSet, ...}  out of files in resource package."""
    data = {}
    files = importlib_resources.files(resource).iterdir()
    for file in files:
        if file.suffix == '.dat':
            dataset = DataSet(datafile=file.read_text())
            for pt in dataset:
                pt.to_conventions()
            data[dataset.id] = dataset
    return data


class KinematicsError(Exception):
    """Exception thrown when kinematics is unphysical."""
    pass


class DataPoint(dict):
    """Kinematical point. May correspond to actual measurement.

    Args:
        xB (float): Bjorken x_B
        t (float): momentum transfer to target squared
        Q2 (float): Q^2
        phi (float): azimuthal angle
        FTn (int): harmonic of azimuthal angle, e.g. -1 for sin(phi)
        observable (str): name of the measured observable
        val (float): measurement value
        err (float): total uncertainty of `val`
        units (dict): pysical units of variables
        kindict (dict): for old, alternative passing of kinematics values,
                 like ``g.DataPoint({'xB': 0.1, 't': -0.2, 'Q2': 4.0})``

    There are other args as well.

    Examples:
            >>> pt = g.DataPoint(xB=0.1, t=-0.2, Q2=4.0)

    Todo:
        `kindict` is temporarily kept for backward compatibility.
        You should not rely on it.

    """

    def __init__(self, kindict: dict = None, **kwargs):
        # Just list some allowed attributes to help mypy:
        # -- Kinematics --
        self.xB = None
        self.xi = None
        self.j = None
        self.t = None
        self.Q2 = None
        # from https://stackoverflow.com/questions/4984647/
        super(DataPoint, self).__init__()
        self.__dict__ = self
        if kindict:
            self.update(kindict)
        self.update(kwargs)
        # calculate other determined kinematic variables:
        _fill_kinematics(self)

    def prepare(self):
        """Pre-calculate some functions of kinematics."""
        kinematics.prepare(self)
        return

    def copy(self):
        """Copy the DataPoint object."""
        # Do we need copy.deepcopy?
        new = copy.copy(self)
        new.__dict__ = new
        return new

    def update_from_grid(self, gridline, dataset):
        """Take data gridline, and update point atributes.

        Args:
            gridline (list): list constructed from one row of data grid in data file.
            dataset (:obj:`DataSet`): container that is to contain this `DataPoint`

        Note:
            It is assumed that data gridline is of the form::

               x1  x2 ....   y1  y1stat y1syst

            where `y1syst` need not be present, or can have one or two values,
            `syst+` and `syst-` (Further elements are ignored.)

        Todo:
            This passing of higher-level DataSet as argument sounds wrong!
            (See comment about aquisition in class docstring.)

        """
        # from https://stackoverflow.com/questions/4984647/
        # we want accessibility via both attributes and dict keys
        super(DataPoint, self).__init__()
        self.__dict__ = self
        self.errtypes = ['err', 'errminus', 'errplus', 'errstat', 'errsyst', 'errnorm']
        # 1. Put reference to container into attribute
        self.dataset = dataset
        # 2. Put data into attributes
        # 2a. first acquire also attributes of parent DataSet
        self.update(dataset.__dict__)
        # 2b. x-axes
        for name in self.xnames:
            nameindex = int(name[1:].split('name')[0])  # = 1, 0, 2, ...
            xname = getattr(self, name)  # = 't', 'xB', ...
            xval = getattr(self, 'x' + str(nameindex) + 'value')  # =  0.1
            if isinstance(xval, float) or isinstance(xval, int):
                # we have global instead of grid value
                setattr(self, xname, xval)    # pt.xB = 0.1, pt.FTn = 1, ...
            else:
                # take value from the grid
                columnindex = int(xval.split('column')[1])-1  # = 1, 0, 2, ...
                setattr(self, xname, gridline[columnindex])  # pt.xB = gridline[1]
        # 2c. y-axis
        self.val = gridline[int(self.y1value.split('column')[1])-1]
        # 2d. y-axis errors
        if 'y1error' in self:  # we are given total error already
            self.err = gridline[int(self.y1error.split('column')[1])-1]
            self.errplus = self.err
            self.errminus = self.err
        else:
            # we have to add various contributions. We do addition of variances:
            # stat = statsym + max(stat+,stat-)
            # syst_uncorr = systsym + max(syst+,syst-)
            # syst = syst_uncorr + norm
            # err = stat + syst   # This is used for fitting chisq
            # Following two are used for pulls:
            # err+ = statsym + stat+ + systsym + syst+ + norm
            # err- = statsym + stat- + systsym + syst- + norm
            varstat = 0
            varsyst = 0  # uncorrelated syst err
            varsym = 0
            varplus = 0
            varminus = 0
            varnorm = 0
            # 1. statistical error
            if 'y1errorstatistic' in self:
                es = gridline[int(self.y1errorstatistic.split('column')[1])-1]**2
                varstat += es
                varsym += es
            if 'y1errorstatisticplus' in self:
                ep = gridline[int(self.y1errorstatisticplus.split('column')[1])-1]**2
                em = gridline[int(self.y1errorstatisticminus.split('column')[1])-1]**2
                varstat += max(ep, em)
                varplus += ep
                varminus += em
            # 2. systematic error
            if 'y1errorsystematic' in self:
                es = gridline[int(self.y1errorsystematic.split('column')[1])-1]**2
                varsyst += es
                varsym += es
            if 'y1errorsystematicplus' in self:
                ep = gridline[int(self.y1errorsystematicplus.split('column')[1])-1]**2
                em = gridline[int(self.y1errorsystematicminus.split('column')[1])-1]**2
                varsyst += max(ep, em)
                varplus += ep
                varminus += em
        # 3. normalization error (specified as percentage)
            if 'y1errornormalization' in self:
                varnorm += (self.y1errornormalization * self.val)**2
            # 4. TOTAL errors
            self.errplus = math.sqrt(varsym + varplus + varnorm)
            self.errminus = math.sqrt(varsym + varminus + varnorm)
            self.errstat = math.sqrt(varstat)
            self.errsyst = math.sqrt(varsyst + varnorm)
            self.errnorm = math.sqrt(varnorm)
            # FIXME: One needs to make a choice here and we go conservative
            self.err = math.sqrt(varstat + varsyst + varnorm)
            # alternative:
            # self.err = (self.errplus+self.errminus)/2.
        # 2e. calculate standard kinematical variables
        _fill_kinematics(self)
        # 2f. polarizations
        # Unpolarized in1 particle
        if 'in1polarization' not in self:
            self.in1polarization = 0
        # For transversaly polarized target set, if needed and
        # if not already set, by default take dominant sine-varphi harmonic
        if ('in2polarizationvector' in self and self.in2polarizationvector == 'T' and
                'varFTn' not in self):
            self.varFTn = -1
        return

    def __repr__(self):
        """Print something useful."""
        try:
            return "DataPoint: " + self.observable + " = " + str(self.val)
        except AttributeError:
            return "DataPoint"

    def to_conventions(self):
        """Transform datapoint into gepard's conventions."""
        if hasattr(self, 'val'):
            self.origval = self.val  # to remember it for later convenience
        if hasattr(self, 'errtypes'):
            for errtype in self.errtypes:
                if hasattr(self, errtype):
                    setattr(self, 'orig'+errtype, getattr(self, errtype))
        # C1. azimutal angle phi should be in radians.
        if 'phi' in self and hasattr(self, 'units') and self.units['phi'][:3] == 'deg':
            self.phi = self.phi * math.pi / 180.
            self.newunits['phi'] = 'rad'
        # C2. phi_{Trento} -> (pi - phi_{BKM})
        if 'frame' in self and self.frame == 'Trento':
            if 'phi' in self:
                self.phi = math.pi - self.phi
            elif 'FTn' in self:
                if self.FTn == 1 or self.FTn == 3 or self.FTn == -2:
                    self.val = - self.val
        # C3. varphi_{Trento} -> (varphi_{BKM} + pi)
            if 'varphi' in self:
                self.varphi = self.varphi - math.pi
            elif 'varFTn' in self:
                if self.varFTn == 1 or self.varFTn == -1:
                    self.val = - self.val
                else:
                    raise ValueError('varFTn = {} not allowed. Only +/-1!'.format(
                                     self.varFTn))
            self.newframe = 'BMK'
        # C4. cross-sections should be in nb
        if hasattr(self, 'units') and self.units[self.observable] == 'pb/GeV^4':
            self.val = self.val/1000
            for errtype in self.errtypes:
                if hasattr(self, errtype):
                    err = getattr(self, errtype)
                    setattr(self, errtype, err/1000)
            self.newunits[self.observable] = 'nb/GeV^4'

    def from_conventions(self):
        """Transform stuff from Approach's conventions into original data's."""
        # C4. cross-sections should be in nb
        if hasattr(self, 'units') and self.units[self.observable] == 'pb/GeV^4':
            self.val = self.val*1000
            for errtype in self.errtypes:
                if hasattr(self, errtype):
                    err = getattr(self, errtype)
                    setattr(self, errtype, err*1000)
        # C2. phi_{BKM} -> (pi - phi_{Trento})
        if 'frame' in self and self.frame == 'Trento':
            if 'phi' in self:
                self.phi = math.pi - self.phi
            elif 'FTn' in self:
                if self.FTn == 1 or self.FTn == 3:
                    self.val = - self.val
        # C3. varphi_{Trento} -> (varphi_{BKM} + pi)
            if 'varphi' in self:
                self.varphi = self.varphi + math.pi
            elif 'varFTn' in self:
                if self.varFTn == 1 or self.varFTn == -1:
                    self.val = - self.val
            self.newframe = 'Trento'
        # C1. azimutal angle phi back to degrees
        if 'phi' in self and hasattr(self, 'units') and self.units['phi'][:3] == 'deg':
            self.phi = self.phi / math.pi * 180.

    def orig_conventions(self, val):
        """Like from_conventions, but for the prediction val."""
        # This doesn't touches self
        # C4. cross-sections nb --> pb
        if hasattr(self, 'units') and self.units[self.observable] == 'pb/GeV^4':
            val = val*1000
        # C2. phi_{BKM} --> (pi - phi_{Trento})
        if 'frame' in self and self.frame == 'Trento' and 'FTn' in self:
            if self.FTn == 1 or self.FTn == 3 or self.FTn == -2:
                val = - val
        if 'frame' in self and self.frame == 'Trento' and 'varFTn' in self:
            if self.varFTn == 1 or self.varFTn == -1:
                val = - val
        return val


class DataSet(list):
    """A container for `DataPoint` instances.

    Args:
        datapoints (list): list of :obj:`DataPoint` instances
        datafile (str): data file to be parsed, represented as string

    Either take explicit list of DataPoints or get them by parsing datafile.

    Information that is common to all data of a given dataset (i.e.
    which is contained in a preamble of datafile is accessible
    via attributes:

    Attributes:
        observable (str): name of observable measured
        collaboration (str): name of experimenatal collaboration
        units (dict): pysical units of variables
        newunits (dict): internal pysical units of variables

    There are other attributes as well.

    """
    def __init__(self, datapoints=None, datafile=None):
        if datapoints:
            list.__init__(self, datapoints)
        else:
            # we have datafile
            list.__init__(self)
            preamble, data = self.parse(datafile)

            # Create needed attributes before creating `DataPoint`s
            # Preamble stuff goes into attributes
            for key in preamble:
                try:  # try to convert to number everything that is number
                    setattr(self, key, _str2num(preamble[key]))
                except ValueError:  # rest stays as is
                    setattr(self, key, preamble[key])

            #  Extracting names of x-axes variables and observable
            #  xnames = ['x1name', 'x2name', ...], not necessarily sorted!
            #  xaxes = ['t', 'xB', ...]
            self.xnames = [key for key in preamble if re.match(r'^x\dname$', key)]
            self.xaxes = [preamble[key] for key in self.xnames]
            self.observable = preamble['y1name']
            # Good to have:
            self.filename = os.path.split(datafile)[-1]
            # Following dictionary will contain units for everything
            # i.e.  {'phi' : 'degrees', 't' : 'GeV^2', ...}
            self.units = dict((preamble[key], preamble[key[:2]+'unit'])
                              for key in self.xnames)
            self.units[self.observable] = preamble['y1unit']
            # Following dictionary will have units which are changed so that match
            # units used for internal theoretical formulas
            self.newunits = {}
            # charge of first particle FIXME: just electron treated
            if self.in1particle in ['e+', 'ep']:    # positron
                self.in1charge = +1
            elif self.in1particle in ['e', 'e-', 'em']:
                self.in1charge = -1
            # Mandelstam s, if specified
            try:
                if self.process in ['ep2epgamma', 'en2engamma']:
                    if self.exptype == 'fixed target':
                        self.s = 2 * Mp * self.in1energy + Mp2
                    elif self.exptype == 'collider':
                        self.s = 2 * self.in1energy * (self.in2energy + math.sqrt(
                            self.in2energy**2 - Mp2)) + Mp2
                    else:
                        pass  # FIXME: raise error
            except AttributeError:
                pass
            if hasattr(self, 'process'):
                self.process_class = process_class_map[self.process]

            for gridline in data:
                pt = DataPoint()
                pt.update_from_grid(gridline, self)
                self.append(pt)

    def __add__(self, rhs):
        """Add datasets.

        Args:
            rhs (:obj:`DataSet`): dataset to be appended

        Returns:
            joined dataset

        http://stackoverflow.com/questions/8180014/how-to-subclass-python-list-without-type-problems.

        """
        return DataSet(datapoints=list.__add__(self, rhs))

    def __repr__(self):
        """Pretty-print dataset."""
        return 'DataSet with {} points'.format(len(self))

    def __getitem__(self, key):
        """Get an element of dataset i.e. datapoint."""
        # From https://stackoverflow.com/questions/2936863/
        if isinstance(key, slice):
            lst = [self[k] for k in range(*key.indices(len(self)))]
            tmp = DataSet(lst)
            tmp.__dict__ = self.__dict__.copy()  # transfer the attributes
            return tmp
        elif isinstance(key, int):
            return list.__getitem__(self, key)
        else:
            raise TypeError("Invalid argument type.")

    def __getslice__(self, start, end):
        """Get a range of elements."""
        # This is called by [:], while [::] calls __getitem__()
        if start >= len(self):
            raise IndexError("""%s has only %d items and your slice
                starts at %d""" % (self, len(self), start))
        return DataSet(self[start:end:None])

    def parse(self, dataFile):
        """Parse `dataresource` string and return tuple (preamble, data).

        Args:
            dataFile (str): datafile to be parsed, as string

        Returns:
            (preamble, data) (tuple):

        `preamble` is dictionary obtained by converting datafile preamble
        items into dictionary items like this::

            y1 = AC from datafile goes into   {'y1' : 'AC', ...}

        `data` is actual numerical grid of experimental data converted
        into list of lists.

        """
        # [First] parsing the formatted ASCII file
        desc = {}   # description preamble (reference, kinematics, ...)
        data = []   # actual data grid  x1 x2  ... y1 dy1_stat dy1_syst ...
        for dataFileLine in dataFile.splitlines():
            # remove comments
            dataFileLine = dataFileLine.split('#')[0]
            # only lines with '=' (premble) or with numbers only (data grid) are parsed
            if re.search(r'=', dataFileLine):
                # converting preamble line into dictionary item
                desctpl = tuple([s.strip() for s in dataFileLine.split("=")])
                desc[desctpl[0]] = desctpl[1]
            if re.match(r'([ \t]*[-\.\d]+[ \t\r]+)+', dataFileLine):
                # FIXME: TAB-delimited columns are not handled! Only spaces are OK.
                snumbers = re.findall(r'[-\.\d]+', dataFileLine)
                numbers = []
                for s in snumbers:
                    f = float(s)
                    if (f - int(f)) == 0:  # we have integer
                        numbers.append(int(f))
                    else:
                        numbers.append(f)
                data.append(list(map(float, numbers)))

        return desc, data

    def df(self):
        """Return pandas DataFrame of a DataSet."""
        attrs = ['observable', 'collaboration', 'id', 'x', 'eta',
                 'xi', 'xB', 'Q2', 't', 'tm', 'in1energy', 'W',
                 'phi', 'FTn', 'varFTn', 'val',
                 'err', 'errminus', 'errplus', 'errstat', 'errsyst', 'errnorm']
        dat = []
        for pt in self:
            row = []
            for atr in attrs:
                if hasattr(pt, atr):
                    row.append(getattr(pt, atr))
                else:
                    row.append(None)
            row.append(pt)
            dat.append(row)
        return pd.DataFrame(dat, columns=attrs+['pt', ])


def _str2num(s: str) -> Union[int, float]:
    """Convert string to number, taking care if it should be int or float.

    Args:
        s: string to be converted

    Returns:
        number Union[int, float]

    http://mail.python.org/pipermail/tutor/2003-November/026136.html

    """
    if "." in s:
        return float(s)
    else:
        return int(s)


def _complete_xBWQ2(kin):
    """Make trio {xB, W, Q2} complete if two of them are given in 'kin'."""
    if 'W' in kin and 'Q2' in kin and 'xB' not in kin:
        kin.xB = kin.Q2 / (kin.W**2 + kin.Q2 - Mp2)
    elif 'xB' in kin and 'Q2' in kin and 'W' not in kin:
        kin.W = math.sqrt(kin.Q2 / kin.xB - kin.Q2 + Mp2)
    elif 'xB' in kin and 'W' in kin and 'Q2' not in kin:
        kin.Q2 = kin.xB * (kin.W**2 - Mp2) / (1. - kin.xB)
    else:
        raise KinematicsError('Exactly two of {xB, W, Q2} should be given.')
    return


def _complete_tmt(kin):
    """Make duo {t, tm} complete if one of them is given in 'kin'."""
    if 't' in kin and 'tm' not in kin:
        assert kin.t <= 0
        kin.tm = - kin.t
    elif 'tm' in kin and 't' not in kin:
        assert kin.tm >= 0
        kin.t = - kin.tm
    else:
        raise KinematicsError('Exactly one of {t, tm} should be given.')
    return


def _fill_kinematics(kin, old={}):
    """Update kinematics in place.

    Complete set of kinematical variables is {xB, t, Q2, W, s, xi, tm, phi}.
    Using standard identities, missing values are calculated, if possible, first
    solely from values given in 'kin', and then, second, using values in 'old',
    if provided.

    Args:
        kin (:obj:`DataPoint`): datapoint to be updated
        old (:obj:`DataPoint`): extra kinematical info

    """
    kkeys = set(kin.keys())
    trio = set(['xB', 'W', 'Q2'])
    if len(trio.intersection(kkeys)) == 3:
        raise KinematicsError('Overdetermined set {xB, W, Q2} given.')
    elif len(trio.intersection(kkeys)) == 2:
        _complete_xBWQ2(kin)
    elif len(trio.intersection(kkeys)) == 1 and old:
        given = trio.intersection(kkeys).pop()  # one variable given in 'kin'
        # We treat only the case when one of {xB, Q2} is given and second is
        # then taken from 'old'
        if given == 'xB':
            kin.Q2 = old.Q2
        elif given == 'Q2':
            kin.xB = old.xB
        _complete_xBWQ2(kin)
    else:
        # We have zero givens, so take all three from 'old'
        if old:
            for key in trio:
                kin.__setattr__(key, old.__getattribute__(key))
    # FIXME: xi is just fixed by xB - it cannot be given by user
    # There are t/Q2 corrections, cf. BMK Eq. (4), but they are
    # formally higher twist and it is maybe sensible to DEFINE xi,
    # the argument of CFF, as follows:
    if 'xB' in kin:
        kin.xi = kin.xB / (2. - kin.xB)
    duo = set(['t', 'tm'])
    if len(duo.intersection(kkeys)) == 2:
        raise KinematicsError('Overdetermined set {t, tm=-t} given.')
    elif len(duo.intersection(kkeys)) == 1:
        _complete_tmt(kin)
    else:
        # We have zero givens, so take both from 'old'
        if old:
            for key in duo:
                kin.__setattr__(key, old.__getattribute__(key))
    # s is just copied from old, if there is one
    if old and 's' in old:
        kin.s = old.s
    # phi and varphi are copied from old, if possible and necessary
    if 'phi' not in kin and 'phi' in old:
        kin.phi = old.phi
    if 'varphi' not in kin and 'varphi' in old:
        kin.varphi = old.varphi


def select(dataset: DataSet, criteria: List[str] = [], logic='AND'):
    """Filter point of dataset satisfying criteria.

    Args:
        dataset: DataSet or any collection of DataPoints
        criteria: list of strings describing selection criteria
        logic: 'AND' or 'OR': how to logically combine criteria

    Returns:
        DataSet of DataPoints satisfying criterion.

    Example:
        criteria=['xB > 0.1', 'observable == ALU']

    """
    selected = []
    for pt in dataset:
        if logic == 'OR':
            for criterion in criteria:
                if eval('pt.'+criterion):
                    selected.append(pt)
        elif logic == 'AND':
            ok = True
            for criterion in criteria:
                if not eval('pt.'+criterion):
                    ok = False
                    break
            if ok:
                selected.append(pt)
    # convert list to DataSet instance
    tmp = DataSet(selected)
    tmp.__dict__ = dataset.__dict__.copy()  # transfer the attributes
    return tmp


def list_data(ids):
    """List basic info about datasets specified by id numbers."""
    if not isinstance(ids, list):
        ids = [ids]
    for id in ids:
        try:
            dt = dset[id]
            ref = dt.reference.replace('arXiv:', '').replace('hep-ex', '').replace(
                  'nucl-ex', '').replace(
                  'from Morgan Murray, draft_90@hermes.desy.de, J. Burns and M. Murray',
                  'Morgan M.').replace('v1', '').replace(
                            'F. Ellinghaus, QCD02', 'Frank E.').replace(
                                    'PRELIMINARY', 'prelim.').strip('[]/ ')
            try:
                ref2 = dt.reference2
            except AttributeError:
                ref2 = ''
            print('[%3i] %8s %3i %9s %10s %s' % (
                  dt.id, dt.collaboration, len(dt), dt.observable, ref, ref2))
        except KeyError:
            pass


def describe_data(pts):
    """Print observables and where they come from."""
    all = []
    print("{:2s} x {:6s}  {:6s}  {:4s}   {:3s} {:12s}".format(
     'npt', 'obs', 'collab', 'FTn', 'id', 'ref.'))
    print(46*'-')
    for pt in pts:
        props = []
        for prop in ['observable', 'collaboration', 'FTn', 'id', 'reference']:
            if hasattr(pt, prop):
                props.append(str(getattr(pt, prop)))
            else:
                props.append('N/A')
        all.append(tuple(props))
    tot = len(all)
    uniqs = set(all)
    cc = 0
    for uniq in sorted(uniqs):
        n = all.count(uniq)
        cc += n
        print("{:2d} x {:6s}  {:6s}  {:4s}   {:3s} {:12s}".format(n, *uniq))
    assert cc == tot
    print(46*'-')
    print("TOTAL = {}".format(tot))


# Load all public datasets and make them globaly available
dset = loaddata(ep2epgamma)
dset.update(loaddata(gammastarp2Mp))
dset.update(loaddata(gammastarp2gammap))
dset.update(loaddata(en2engamma))
dset.update(loaddata(DIS))
