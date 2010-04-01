""" 
Some utility stuff needed all over the application 

AttrDict -- dictionary with attribute-style access
loaddata -- loads datafiles from directory
fill_kinematics --- calculates missing kinematical variables
parse -- parses datafiles
subplot -- creates subplot for matplotlib plot
npars -- number of free (not fixed) parameters of minuit object
str2num -- transforms string to float or int
prettyprint -- formatted printout of numbers
flatten -- flattens tuples
listFiles -- listfiles in subdirs matching pattern
"""

import os, re, string, fnmatch
import numpy as np

import Data
from constants import toTeX, Mp, Mp2

#from IPython.Debugger import Tracer; debug_here = Tracer()


class KinematicsError(Exception):
    pass

class AttrDict(dict):
    """A dictionary with attribute-style access. 
    
    It maps attribute access to the real dictionary.  
    By Keith Darth, http://code.activestate.com/recipes/473786/
    FIXME: Using this is bad for performance. __getitem__ is called 
    way to many times.
    """
    def __init__(self, init={}):
        dict.__init__(self, init)

    def __getstate__(self):
        return self.__dict__.items()

    def __setstate__(self, items):
        for key, val in items:
            self.__dict__[key] = val

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, dict.__repr__(self))

    def __str__(self):
        s = ""
        for key in self.keys():
            s += '%4s -> % .3g\n' % (key, self[key])
        s = s[:-1]
        return s

    def __setitem__(self, key, value):
        return super(AttrDict, self).__setitem__(key, value)

    def __getitem__(self, name):
        return super(AttrDict, self).__getitem__(name)

    def __delitem__(self, name):
        return super(AttrDict, self).__delitem__(name)


    __getattr__ = __getitem__
    __setattr__ = __setitem__

    def copy(self):
        ch = AttrDict(self)
        return ch


def loaddata(datadir='data', approach=None):
    """Return dictionary {id : `DataSet`, ...}  out of datadir/*dat files."""
    data = {}
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] == ".dat":
            dataset = Data.DataSet(datafile=os.path.join(datadir, file))
            if approach:
                [pt.to_conventions(approach) for pt in dataset]
                [pt.prepare(approach) for pt in dataset]
            data[dataset.id] = dataset
    return data

def _complete_xBWQ2(kin):
    """Make trio {xB, W, Q2} complete if two of them are given in 'kin'."""
    if kin.has_key('W') and kin.has_key('Q2') and not kin.has_key('xB'):
        kin.xB = kin.Q2 / (kin.W**2 + kin.Q2 - Mp2)
    elif kin.has_key('xB') and kin.has_key('Q2') and not kin.has_key('W'):
        kin.W = np.sqrt(kin.Q2 / kin.xB - kin.Q2 + Mp2)
    elif kin.has_key('xB') and kin.has_key('W') and not kin.has_key('Q2'):
        kin.Q2 = kin.xB * (kin.W**2 - Mp2) / (1. - kin.xB)
    else:
        raise KinematicsError, 'Exactly two of {xB, W, Q2} should be given.'
    return

def _complete_tmt(kin):
    """Make duo {t, tm} complete if one of them is given in 'kin'."""
    if kin.has_key('t') and not kin.has_key('tm'):
        assert kin.t <= 0
        kin.tm = - kin.t
    elif kin.has_key('tm') and not kin.has_key('t'):
        assert kin.tm >= 0
        kin.t = - kin.tm
    else:
        raise KinematicsError, 'Exactly one of {t, tm} should be given.'
    return

def fill_kinematics(kin, old={}):
    """Return complete up-to-date kinematical dictionary.
    
    Complete set of kinematical variables is {xB, t, Q2, W, s, xi, tm, phi}.
    Using standard identities, missing values are calculated, if possible, first
    solely from values given in 'kin', and then, second, using values in 'old',
    if provided.

    """
    if not (isinstance(kin, Data.DataPoint) or 
            isinstance(kin, Data.DummyPoint) or isinstance(kin, AttrDict)):
        kin = AttrDict(kin)  # fixing kin
    kkeys = set(kin.keys())
    trio = set(['xB', 'W', 'Q2'])
    if len(trio.intersection(kkeys)) == 3:
        raise KinematicsError, 'Overdetermined set {xB, W, Q2} given.'
    elif len(trio.intersection(kkeys)) == 2:
        _complete_xBWQ2(kin)
    elif len(trio.intersection(kkeys)) == 1 and old:
        given = trio.intersection(kkeys).pop() # one variable given in 'kin'
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
    kin.xi = kin.xB / (2. - kin.xB)
    duo = set(['t', 'tm'])
    if len(duo.intersection(kkeys)) == 2:
        raise KinematicsError, 'Overdetermined set {t, tm=-t} given.'
    elif len(duo.intersection(kkeys)) == 1:
        _complete_tmt(kin)
    else:
        # We have zero givens, so take both from 'old'
        if old:
            for key in duo:
                kin.__setattr__(key, old.__getattribute__(key))
    # s is just copied from old, if there is one
    if old and old.has_key('s'):
        kin.s = old.s
    # phi is copied from old, if possible and necessary
    if not kin.has_key('phi') and old.has_key('phi'):
        kin.phi = old.phi
    return kin

def parse(datafile):
    """Parse `datafile` and return tuple (preamble, data).

    `preamble` is dictionary obtained by converting datafile preamble
    items into dictionary items like this:
        y0 = BCA from datafile goes into   {'y0' : 'BCA', ...}

    `data` is actual numerical grid of experimental data converted 
    into list of lists
    
    """

    # [First] parsing the formatted ASCII file
    desc = {}   # description preamble (reference, kinematics, ...)
    data = []   # actual data grid  x1 x2  ... y1 dy1_stat dy1_syst ...
    dataFile = open(datafile, 'r')
    dataFileLine = dataFile.readline()
    while dataFileLine:
        # remove comments
        dataFileLine = dataFileLine.split('#')[0]
        # only lines with '=' (premble) or with numbers only (data grid) are parsed
        if re.search(r'=', dataFileLine):
            # converting preamble line into dictionary item
            desctpl = tuple(map(string.strip,string.split(dataFileLine,"=")))
            desc[desctpl[0]] = desctpl[1] 
        if re.match(r'([ ]*[-\.\d]+[ \r]+)+', dataFileLine):
            # FIXME: TAB-delimited columns are not handled! Only spaces are OK.
            snumbers = re.findall(r'[-\.\d]+', dataFileLine)
            numbers = []
            for s in snumbers:
                f = float(s)
                if (f - int(f)) == 0:  # we have integer
                    numbers.append(int(f))
                else:
                    numbers.append(f)
            data.append(map(float, numbers))
        dataFileLine = dataFile.readline()

    return desc, data


def subplot(ax, datasets, xaxis=None, kinlabels=[], fits=[]):
    """Plot datapoints together with fit/theory line(s).

    ax -- subplot of matplotlib's figure i.e. ax = figure.add_subplot(..)
    datasets --  list of `DataSet` instances to be plotted
    xaxis -- abscissa variable; if None, last of dataset.xaxes is taken
    kinlabels -- list of constant kinematic variables whose values will
                 be put on plot
    fits -- list of parameter sets/tuples describing fit curves for
            plotting. 

    """
    # first, fix the input if needed
    if not isinstance(kinlabels, list): kinlabels = [kinlabels]
    if not isinstance(fits, list): fits = [fits]
    if not xaxis: xaxis = dataset.xaxes[-1]
    # Data sets (or fits with errorbars)
    setshapes = ['o', 's']  # first circles, then squares ...
    setcolors = ['blue', 'black']  # circles are blue, squares are black, ...
    setn = 0
    for dataset in datasets:
        xval = []; yval = []; yerr = []
        for pt in dataset:
            xval.append(getattr(pt, xaxis)) 
            yval.append(pt.val)
            yerr.append(pt.err)
        ax.errorbar(xval, yval, yerr, linestyle='None', elinewidth=setn+1, 
                marker=setshapes[setn], color=setcolors[setn])
        setn += 1
    # Fit lines
    shapes = ['s', '^', 'd', 'h']  # first squares, then triangles, diamonds, hexagons
    colors = ['red', 'green', 'brown', 'purple']  # squares are red, etc.
    styles = ['-', '--', '-.', ':']
    fitn = 0
    for (approach, pars) in fits:
        # take abscissae from dataset
        line = [getattr(approach, pt.yaxis)(pt, pars) for pt 
                in datasets[0]]
        ## join the dots
        ax.plot(xval, line, color=colors[fitn], linestyle=styles[fitn], linewidth=2)
        ## put symbols on dots
        #ax.plot(xval, line, shapes[fitn], markersize=5,
        #        markerfacecolor=colors[fitn], markeredgecolor='black')
        fitn += 1
    # axes labels
    ax.set_xlabel(toTeX[xaxis], fontsize=15)
    ax.set_ylabel(toTeX[dataset[0].yaxis], fontsize=18)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    # constant kinematic variables positioning
    labx = min(0, min(xval)) + (max(xval) - min(0, min(xval))) * 0.5
    laby = min(0, min(yval)) + (max(yval) - min(0, min(yval))) * 0.05
    labtxt = ""
    for lab in kinlabels:
        try:
            labtxt += toTeX[lab] + ' = ' + str(getattr(dataset,lab)) + ', '
        except AttributeError:
            # If dataset doesn't have it, all points should have it 
            labtxt += toTeX[lab] + ' = ' + str(getattr(dataset[0],lab)) + ', '
    ax.text(labx, laby, labtxt[:-2])
    # ax.frame.set_linewidth(5) # ???
    return

def npars(m):
    """Return number of free (not fixed) parameters of MINUIT object m."""

    n = 0
    for key in m.fixed:
        if not m.fixed[key]:
            n += 1
    return n

def str2num(s):
    """Convert string to number, taking care if it should be int or float.
    
    http://mail.python.org/pipermail/tutor/2003-November/026136.html
    """

    return ("." in s and [float(s)] or [int(s)])[0]


def prettyprint(all_numbers):
    """Formatted printout of fit values. Using algorithm by Alex Martelli."""

    formatted_numbers = []
    space_on_the_left = []
    for number in all_numbers:
        formatted_string = "% 1.1f" % number
        formatted_numbers.append(formatted_string)
        on_the_left = formatted_string.find('.')
        if on_the_left < 0: on_the_left = len(formatted_string)
        space_on_the_left.append(on_the_left)

    left_total = max(space_on_the_left)

    for s, l in zip(formatted_numbers, space_on_the_left):
        padding = ' '*(left_total - l)
        print '%s%s' % (padding, s)
    return

def flatten(T):
    """Flatten the tuple."""
    if not isinstance(T, tuple): return (T,)
    elif len(T) == 0: return ()
    else: return flatten(T[0]) + flatten(T[1:]) 

def listFiles(root, patterns='*', recurse=1, return_folders=0):
    """Return the list of files matching patterns, recursively.

    From 4.19 Walking Directory Trees in Python Cookbook

    """
    # Expand patterns from semicolon-separated string to list
    pattern_list = patterns.split(';')
    # Collect input and output arguments into one bunch
    class Bunch:
        def __init__(self, **kwds):
            self.__dict__.update(kwds)

    arg = Bunch(recurse=recurse, pattern_list=pattern_list,
            return_folders=return_folders, results=[])

    def visit(arg, dirname, files):
        # Append to arg.results all relevant files (and perhaps folders)
        for name in files:
            fullname = os.path.normpath(os.path.join(dirname, name))
            if arg.return_folders or os.path.isfile(fullname):
                for pattern in arg.pattern_list:
                    if fnmatch.fnmatch(name, pattern):
                        arg.results.append(fullname)
                        break
        # Block recursion if recursion was disallowed
        if not arg.recurse: files[:]=[]

    os.path.walk(root, visit, arg)

    return arg.results


def select(dataset, criteria=[], logic='AND'):
    """Return list (not DataSet) of DataPoints satisfying criteria.

    logic='OR': select points satisfying any
    of the list of criteria.
    Example: criteria=['xB > 0.1', 'y0name == BSA']
    
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
    return selected

