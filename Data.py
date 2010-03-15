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

import os, re 

import pylab as plt
from numpy import sqrt

import utils
from constants import Mp, Mp2


class DataPoint(object):

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

              x0  x1 ....   y0  y0stat y0syst 

        (Further elements are ignored.)
        `dataset` is container `DataSet` that is to contain this `DataPoint`
        FIXME: this passing of higher-level DataSet as argument sounds wrong!
                (See comment about aquisition in class docstring.|

        """

        # 1. Put reference to container into attribute
        self.dataset = dataset
        # 2. Put data into attributes
        # 2a. first acquire also attributes of parent DataSet
        self.__dict__.update(dataset.__dict__)
        # 2b. x-axes
        i = 0 # index of our position on gridline
        for xaxis in self.xaxes:
            setattr(self, xaxis, gridline[i])
            i = i + 1
        # 2c. y-axis 
        self.val = gridline[i]
        # 2d. y-axis errors
        try:
            self.stat = gridline[i+1]
            self.syst = gridline[i+2]
            self.err = sqrt( self.stat**2 + self.syst**2 )
        except IndexError: # we have just one error number
            self.err = gridline[i+1]
        # 2e. calculate standard kinematical variables
        if self.has('W') and self.has('Q2'):
            self.xB = self.Q2 / (self.W**2 + self.Q2 - Mp2)
        if self.has('xB') and self.has('Q2'):
            self.W = sqrt(self.Q2 / self.xB - self.Q2 + Mp2)
        if self.has('xB'):
            # There are t/Q2 corrections, cf. BMK Eq. (4), but they are 
            # formally higher twist and it is maybe sensible to DEFINE xi, 
            # the argument of CFF, as follows:
            self.xi = self.xB / (2. - self.xB)
        if self.has('t'):
            self.mt = - self.t
        elif self.has('mt'):
            self.t = - self.mt
        return

    def __repr__(self):
        return "DataPoint. Measurement: " + self.yaxis + " = " + str(self.val)

    def has(self, name):
        """Does point (or dataset) have attribute `name`?"""
        if self.dataset.__dict__.has_key(name) or self.__dict__.has_key(name):
            return True
        return False

    def to_conventions(self, approach):
        approach.to_conventions(self)
        return

    def prepare(self, approach):
        approach.prepare(self)
        return


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

            #  Extracting x-axes variables 
            xs = [key for key in preamble if re.match('^x\d$', key)]
            xs.sort()  # xs = ['x0', 'x1', ...]
            self.xaxes = [preamble[key] for key in xs]  # xaxes = ['t', 'xB', ...]

            # Good to have:
            self.yaxis = preamble['y0']
            self.filename = os.path.split(datafile)[-1]
            # Following dictionary will contain units for everything
            # i.e.  {'phi' : 'degrees', 't' : 'GeV^2', ...}
            self.units = dict((preamble[key], preamble[key+'unit']) for key in xs)
            self.units[self.yaxis] = preamble['y0unit']
            # Following dictionary will have units which are changed so that match
            # units used for internal theoretical formulas
            self.newunits = {}
            # charge of first particle FIXME: just electron treated
            if self.in1 == 'ep':                      # positron
                self.charge = +1
            elif self.in1 == 'e' or self.in1 == 'em':   # electron
                self.charge = -1
            # Mandelstam s
            if self.exptype == 'fixed target':
                self.s = 2 * Mp * self.in1energy + Mp2
            elif self.exptype == 'collider':
                self.s = 2 * self.in1energy * (self.in2energy + sqrt(
                    self.in2energy**2 - Mp2)) + Mp2
            else:
                pass # FIXME: raise error

            for gridline in data:
                self.append(DataPoint(gridline, self))
    
    def __repr__(self):
        return "DataSet instance from '" + self.filename + "'"

    def __getslice__(self, start, end, step=None):
        if start >= len(self):
            raise IndexError, """%s has only %d items and your slice 
                starts at %d""" % (self, len(self), start)
        tmp = DataSet(self[start:end:step])
        tmp.__dict__ = self.__dict__.copy() # transfer the attributes
        return tmp

    def plot(self, xaxis=None, kinlabels=[], fits=[], path=None, fmt='png'):
        """Plot the dataset, with fit lines if needed.

        Named arguments correspond to those of `utils.subplot`.
        If path to directory is given, figure is saved there in format `fmt`.

        """

        title = self.filename
        fig = plt.figure()
        fig.canvas.set_window_title(title)
        fig.suptitle(title)
        ax = fig.add_subplot(111)
        utils.subplot(ax, self, xaxis, kinlabels, fits)
        if path:
            fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
        else:
            fig.canvas.draw()
            fig.show()
        return 


class DummyPoint(object):
    """This is only used for creating simple DataPoint-like objects"""

    def __init__(self, init=None):
        if init:
            self.__dict__.update(init)

    def has(self, name):
        if self.__dict__.has_key(name):
            return True
        return False

    def prepare(self, approach):
        approach.prepare(self)
        return
