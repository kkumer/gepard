""" 
Some utility stuff needed all over the application 

AttrDict -- dictionary with attribute-style access
loaddata -- loads datafiles from directory
fill_kinematics --- calculates missing kinematical variables
parse -- parses datafiles
npars -- number of free (not fixed) parameters of minuit object
str2num -- transforms string to float or int
prettyprint -- formatted printout of numbers
flatten -- flattens tuples
listFiles -- listfiles in subdirs matching pattern
select -- selecting DataPoints according to criteria
listdb --  listing the content of database of models
"""

import os, re, string, fnmatch
import numpy as np

import Data, Approach
from constants import Mp, Mp2

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

def loaddata(datadir='data', approach=Approach.hotfixedBMK):
    """Return dictionary {id : `DataSet`, ...}  out of datadir/*dat files.
    
    approach defines conventions for frame, kinematics etc. to which data
    is adapted. By default it is hotfixedBMK (Belitsky, Mueller).
    
    """
    data = {}
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] == ".dat":
            dataset = Data.DataSet(datafile=os.path.join(datadir, file))
            if approach and dataset.process == 'ep2epgamma':
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

def listdb(db):
    print "%-8s  |  %s" % ('name', 'description')
    print "%-8s--+--%s" % (8*'-', 60*'-')
    for key in db:
        print "%-8s  |  %s" % (key, db[key].description)
