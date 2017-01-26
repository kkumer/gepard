"""
Classes for representing experimental data.

DataPoint -- object containing one experimental measurement
DataSet -- container holding `DataPoint` instances
DummyPoint -- class for points which have just few relevant attributes

"""

# FIXME: Should maybe use "properties" for setting attributes
# that should not be changed by user. See python help(property).
# To the same end, see about descriptors (properties are actually
# implemented using descriptor protocol).
# FIXME: Could use __slots__ to optimize DataPoint and the
# number and names of attributes

import os, re, math, copy
from numpy import pi

import pandas as pd

import utils
from constants import Mp, Mp2


class DummyPoint(object):
    """This is only used for creating simple DataPoint-like objects"""

    def __init__(self, init=None):
        if init:
            self.__dict__.update(init)

    def has_key(self, name):
        """Does point (or dataset) have attribute `name`?"""
        if self.__dict__.has_key(name):
            return True
        return False

    def keys(self):
        """Simulate dictionary."""
        return self.__dict__.keys()

    def to_conventions(self, approach):
        approach.to_conventions(self)
        return

    def prepare(self, approach):
        approach.prepare(self)
        return

    def copy(self):
        # Do we need copy.deepcopy?
        return copy.copy(self)

class DataPoint(DummyPoint):

    """One experimental measurement with all necessary information 
    about kinematics, all contained in attributes. E.g.

    `xB` -- x_Bjorken
    `Q2` -- squared momentum of virtual photon
    `val` -- measurement value
    `stat` -- statistical error of `val`
    `syst` -- systematic error of `val`
    `err` --  `stat` and `syst` added in quadrature
     ...

    Information that is common to all data of a given dataset (i.e.
    which is contained in a preamble of datafile is accessible 
    via `dataset` attribute:

    `dataset.yaxis` -- name of observable measured
    `dataset.xaxes` -- names of kinematical x-axes of data
    `dataset.collaboration` -- name of experimenatal collaboration
    `dataset.units` -- dictionary with pysical units of variables
    `dataset.newunits` -- dictionary with internal pysical units of variables
     ...

    For user's and programmer's convenience, these `dataset` attributes 
    are also inherited by `DataPoint` objects, so 
        point.dataset.yaxis == point.yaxis 
    (Is this type of ineritance, know also as "aquisition", good idea?)

    """ 

    def __init__(self, gridline, dataset):
        """Take data gridline, construct `DataPoint` object and append it to dataset

        `gridline` is a list constructed from one row of data grid in data file.
        It is assumed that data gridline is of the form:

              x1  x2 ....   y1  y1stat y1syst  

        where y1syst need not be present, or can have one or two values, syst+ and syst-

        (Further elements are ignored.)
        `dataset` is container `DataSet` that is to contain this `DataPoint`
        FIXME: this passing of higher-level DataSet as argument sounds wrong!
                (See comment about aquisition in class docstring.|

        FIXME: Should DataPoint be subclas of dict?
        """

        # 1. Put reference to container into attribute
        self.dataset = dataset
        # 2. Put data into attributes
        # 2a. first acquire also attributes of parent DataSet
        self.__dict__.update(dataset.__dict__)
        # 2b. x-axes
        for name in self.xnames:
            nameindex = int(name[1:].split('name')[0])  # = 1, 0, 2, ...
            xname = getattr(self, name)  # = 't', 'xB', ...
            xval = getattr(self, 'x' + str(nameindex) + 'value') #  =  0.1
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
        if self.has_key('y1error'):  # we are given total error already
            self.err = gridline[int(self.y1error.split('column')[1])-1]
            self.errplus = self.err
            self.errminus = self.err
        else:  # we have to add various contributions. We do addition of variances.
            self.varsym = 0     # symmetric contributions
            self.varplus = 0
            self.varminus = 0
            # 1. statistical error
            if self.has_key('y1errorstatistic'):
                self.varsym += gridline[int(self.y1errorstatistic.split('column')[1])-1]**2
            if self.has_key('y1errorstatisticplus'): 
                self.varplus += gridline[int(self.y1errorstatisticplus.split('column')[1])-1]**2
                self.varminus += gridline[int(self.y1errorstatisticminus.split('column')[1])-1]**2
            # 2. systematic error
            if self.has_key('y1errorsystematic'):
                self.varsym += gridline[int(self.y1errorsystematic.split('column')[1])-1]**2
            if self.has_key('y1errorsystematicplus'): 
                self.varplus += gridline[int(self.y1errorsystematicplus.split('column')[1])-1]**2
                self.varminus += gridline[int(self.y1errorsystematicminus.split('column')[1])-1]**2
	    # 3. normalization error
            if self.has_key('y1errornormalization'):
                self.varsym += (self.y1errornormalization * self.val)**2
            # 4. TOTAL variances and errors
            self.varplus += self.varsym  
            self.varminus += self.varsym
            self.errplus = math.sqrt(self.varplus)
            self.errminus = math.sqrt(self.varminus)
            self.err = (self.errplus+self.errminus)/2.
        # 2e. calculate standard kinematical variables
        utils.fill_kinematics(self)
        # 2f. polarizations
        # Unpolarized in1 particle
        if not self.has_key('in1polarization'):
            self.in1polarization = 0
        # For transversaly polarized target set, if needed and
        # if not already set, by default take dominant sine-varphi harmonic
        if (self.has_key('in2polarizationvector') and self.in2polarizationvector == 'T' and 
                not self.has_key('varFTn')):
            self.varFTn = -1
        return

    def __repr__(self):
        return "DataPoint: " + self.yaxis + " = " + str(self.val)


class DataSet(list):

    """A container for `DataPoint` instances.

    Information that is common to all data of a given dataset (i.e.
    which is contained in a preamble of datafile is accessible 
    via attributes:

    `yaxis` -- name of observable measured
    `xaxes` -- names of kinematical x-axes of data
    `collaboration` -- name of experimenatal collaboration
    `units` -- dictionary with pysical units of variables
    `newunits` -- dictionary with internal pysical units of variables
     ...

    """

    def __init__(self, datapoints=None, datafile=None):
        """Either take explicit list of `DataPoint`s or get them by parsing datafile."""

        if datapoints:
            list.__init__(self, datapoints)
        else:
            # we have datafile
            list.__init__(self)
            preamble, data = utils.parse(datafile)

            # Create needed attributes before creating `DataPoint`s
            # Preamble stuff goes into attributes
            for key in preamble:
                try: # try to convert to number everything that is number
                    setattr(self, key, utils.str2num(preamble[key]))
                except ValueError: # rest stays as is
                    setattr(self, key, preamble[key])

            #  Extracting names of x-axes variables 
            #  xnames = ['x1name', 'x2name', ...], not necessarily sorted!
            #  xaxes = ['t', 'xB', ...]
            self.xnames = [key for key in preamble if re.match('^x\dname$', key)]
            self.xaxes = [preamble[key] for key in self.xnames]

            # Good to have:
            self.yaxis = preamble['y1name']
            self.filename = os.path.split(datafile)[-1]
            # Following dictionary will contain units for everything
            # i.e.  {'phi' : 'degrees', 't' : 'GeV^2', ...}
            self.units = dict((preamble[key], preamble[key[:2]+'unit']) for key in self.xnames)
            self.units[self.yaxis] = preamble['y1unit']
            # Following dictionary will have units which are changed so that match
            # units used for internal theoretical formulas
            self.newunits = {}
            # charge of first particle FIXME: just electron treated
            if self.in1particle == 'e+' or self.in1particle == 'ep':    # positron
                self.in1charge = +1
            elif self.in1particle == 'e' or self.in1particle == 'e-' or self.in1particle == 'em':
                self.in1charge = -1
            # Mandelstam s
            if self.process == 'ep2epgamma':
                if self.exptype == 'fixed target':
                    self.s = 2 * Mp * self.in1energy + Mp2
                elif self.exptype == 'collider':
                    self.s = 2 * self.in1energy * (self.in2energy + math.sqrt(
                        self.in2energy**2 - Mp2)) + Mp2
                else:
                    pass # FIXME: raise error

            for gridline in data:
                self.append(DataPoint(gridline, self))
    
    def __add__(self, rhs):
        """http://stackoverflow.com/questions/8180014/how-to-subclass-python-list-without-type-problems"""
        return DataSet(datapoints=list.__add__(self,rhs))

    def __repr__(self):
        return "DataSet instance from '" + self.filename + "'"

    def __getslice__(self, start, end, step=None):
        if start >= len(self):
            raise IndexError, """%s has only %d items and your slice 
                starts at %d""" % (self, len(self), start)
        tmp = DataSet(self[start:end:step])
        tmp.__dict__ = self.__dict__.copy() # transfer the attributes
        return tmp

    def df(self):
        """Return pandas DataFrame of dataset."""
        attrs = ['y1name', 'collaboration', 'FTn', 'id',
                'xB', 'Q2', 't', 'tm', 'W', 'phi', 'FTn', 'varFTn', 'val',
                'err', 'stat', 'syst']#, 'systplus', 'systminus',
                #'statplus', 'statminus', 'normerr']
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
        return pd.DataFrame(dat, columns=attrs+['pt',])


