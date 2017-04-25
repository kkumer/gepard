""" 
Some utility stuff needed all over the application 

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
hubDict -- merges two dicts, but not actually but by forwarding
stringcolor -- coloring string for output, if possible
FTF -- Fourier series fit
FTFMC -- Fourier series fit with MC error propagation
cvsets -- n-fold cross-validation sets
FTanalyse -- Determine the number of harmonics present in bins according to 
             n-fold cross-validation.
describe_data -- List observables in dataset.
compare_harmonics --  
"""

import os, re, string, fnmatch, itertools, logging
import numpy as np
import pandas as pd

import Data, Approach
from constants import Mp, Mp2

#from IPython.Debugger import Tracer; debug_here = Tracer()
_lg = logging.getLogger('p.%s' % __name__)
_lg.debug('Loading module %s' % __name__)


class KinematicsError(Exception):
    pass


def loaddata(datadir='data', approach=False):
    """Return dictionary {id : `DataSet`, ...}  out of datadir/*dat files.
    
    approach defines conventions for frame, kinematics etc. to which data
    is adapted.
    
    """
    data = {}
    for file in os.listdir(datadir):
        #_lg.debug('Loading datafile %s' % file)
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
    # phi and varphi are copied from old, if possible and necessary
    if not kin.has_key('phi') and old.has_key('phi'):
        kin.phi = old.phi
    if not kin.has_key('varphi') and old.has_key('varphi'):
        kin.varphi = old.varphi
    return kin

def parse(datafile):
    """Parse `datafile` and return tuple (preamble, data).

    `preamble` is dictionary obtained by converting datafile preamble
    items into dictionary items like this:
        y1 = BCA from datafile goes into   {'y1' : 'BCA', ...}

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
            data.append(map(float, numbers))
        dataFileLine = dataFile.readline()

    return desc, data

def npars(m):
    """Return number of free (not fixed) parameters of MINUIT object m."""

    n = 0
    for key in m.parameter_names:
        if m.parameters.has_key('fix_'+key) and not m.parameters['fix_'+key]:
            n += 1
    return n

def str2num(s):
    """Convert string to number, taking care if it should be int or float.
    
    http://mail.python.org/pipermail/tutor/2003-November/026136.html
    """

    return ("." in s and [float(s)] or [int(s)])[0]


def prettyprint(all_numbers):
    """Formatted printout of fit values. Using algorithm by Alex Martelli."""

    # Alternative:
    # http://stackoverflow.com/questions/1025379/decimal-alignment-formatting-in-python

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
    Example: criteria=['xB > 0.1', 'y1name == BSA']
    
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
    tmp = Data.DataSet(selected)
    tmp.__dict__ = dataset.__dict__.copy() # transfer the attributes
    return tmp

def listdb(db):
    print "%-17s--+--%s" % (17*'-', 60*'-')
    print "%-17s  |  %s" % ('name', 'description')
    print "%-17s--+--%s" % (17*'-', 60*'-')
    for key in db:
        print "%-17s  |  %s" % (key, db[key].description)
    #print "\n WARNING: gepard models are now likely broken. Reinitialize them!"

def listdata(ids, data):
    """List basic info about datasets specified by id numbers."""
    if not isinstance(ids, list): ids = [ids]
    for id in ids:
        try:
            dt = data[id]
            ref = dt.reference.replace('arXiv:', '').replace('hep-ex', '').replace('nucl-ex', '').replace('from Morgan Murray, draft_90@hermes.desy.de, J. Burns and M. Murray', 'Morgan M.').replace('v1', '').replace('F. Ellinghaus, QCD02', 'Frank E.').replace('PRELIMINARY', 'prelim.').strip('[]/ ')
            try:
                ref2 = dt.reference2
            except:
                ref2 =  ''
            print '[%3i] %8s %3i %9s %10s %s' % (dt.id, dt.collaboration, len(dt), dt.y1name, ref, ref2)
        except KeyError:
            pass

def listchis(ths, Q2cut=1., Q2max=1.e3, nsets=0, out='chis'):
    """Compare theories for subsets of data.
    
    What is printed out depends on 'out' keyword argument:
     'chis'  --  chisquares
    'probs'  --  probabilities of chisquares
    'pulls'  --  sum((th-exp)/err)/sqrt(size)

    """
    if not isinstance(ths, list): ths = [ths]
    from abbrevs import H1ZEUS, ALUIpts, BCApts, CLASpts, BSDwpoints, BSSwpoints,\
            AULpts, ALLpts, AUTIpts, CLAS14BSApts, CLAS14TSApts, CLAS14BTSApts,\
            BSACLAS_KKpoints, UNP5points, ALTGLO5points,\
            BSACLAS_DMpoints, CLASTSApts,\
            AUTICSpts, CLASKKpts, AUTDVCSpts, H_AULpts, C_AULpts,\
            H_BSDwpts, H_BSSw0pts, H_BSSw1pts,\
            H_BSDpts, H_BSS0pts, H_BSS1pts,\
            C_BSDwpts, C_BSSw0pts, C_BSSw1pts, H_BSD, H_BSS,\
            C_BSD, C_BSS
    #exps[0] = ['UNP5points', 'ALTGLO5', 'CLAS', 'CLASDM', 'BSDw', 'BSSw', 'TSA1', 'BTSA', 'TPpoints']
    #ptssets[0] = [UNP5points, ALTGLO5points, data[25], data[8], BSDwpoints, BSSwpoints, TSA1points, BTSApoints, TPpoints]
    sets = {}
    sets[0] = [('H1ZEUS', 'X_DVCS', H1ZEUS), ('HERMES', 'ALUI', ALUIpts),
            ('HERMES', 'BCA', BCApts), ('CLAS', 'BSA', CLASpts),
            ('Hall A', 'BSDw', BSDwpoints), ('Hall A', 'BSSw', BSSwpoints),
            ('HRM/CLS', 'AUL', AULpts), ('HERMES', 'ALL', ALLpts),
            ('HERMES', 'AUTI', AUTIpts)]
    sets[1] = [('CLAS07_KK', 'BSA', BSACLAS_KKpoints),
               ('CLAS14_KK', 'BSA', CLAS14BSApts),
               ('CLAS14_KK', 'TSA', CLAS14TSApts), 
               ('CLAS14_KK', 'BTSA', CLAS14BTSApts)]
    sets[2] = [('CLAS07_DM', 'BSA', BSACLAS_DMpoints),
               ('CLAS06', 'TSA', CLASTSApts),
               ('CLAS14_KK', 'BSA', CLAS14BSApts),
               ('CLAS14_KK', 'TSA', CLAS14TSApts), 
               ('CLAS14_KK', 'BTSA', CLAS14BTSApts)]
    sets[3] = [('CLAS0708', 'BSA', CLASKKpts),
            ('CLAS', 'AUL', C_AULpts), 
            ('CLAS14_KK', 'BSA', CLAS14BSApts),
            ('CLAS14_KK', 'TSA', CLAS14TSApts), 
            ('CLAS14_KK', 'BTSA', CLAS14BTSApts)]
    sets[4] = [ ('HERMES', 'BCA', BCApts), 
            ('HERMES', 'ALUI', ALUIpts),
            ('HERMES', 'AUL', H_AULpts), 
            ('HERMES', 'ALL', ALLpts),
            ('HERMES', 'AUTI', AUTIpts),
            ('HERMES', 'AUTICS', AUTICSpts),
            ('HERMES', 'AUTDVCS', AUTDVCSpts)]
    sets[5] = [
            ('H1ZEUS', 'X_DVCS', H1ZEUS), 
            ('HERMES', 'ALUI', ALUIpts),
            ('HERMES', 'BCA', BCApts),
            ('HRM/CLS', 'AUL', AULpts), ('HERMES', 'ALL', ALLpts),
            ('HERMES', 'AUTI', AUTIpts),
            ('CLAS', 'BSA', CLAS14BSApts),
            ('CLAS', 'TSA', CLAS14TSApts), 
            ('CLAS', 'BTSA', CLAS14BTSApts),
            ('CLAS', 'BSDw_s1', C_BSDwpts), ('CLAS', 'BSSw_c0', C_BSSw0pts),
            ('CLAS', 'BSSw_c1', C_BSSw1pts),
            ('Hall A', 'BSDw_s1', H_BSDwpts), ('Hall A', 'BSSw_c0', H_BSSw0pts),
            ('Hall A', 'BSSw_c1', H_BSSw1pts)
            ]
    sets[6] = [
            ('CLAS', 'BSDw_s1', C_BSDwpts), ('CLAS', 'BSSw_c0', C_BSSw0pts),
            ('CLAS', 'BSSw_c1', C_BSSw1pts),
            ('Hall A', 'BSDw_s1', H_BSDwpts), ('Hall A', 'BSSw_c0', H_BSSw0pts),
            ('Hall A', 'BSSw_c1', H_BSSw1pts)
            ]
    sets[7] = [
            ('CLAS', 'BSDw_s1', C_BSDwpts), ('CLAS', 'BSSw_c0', C_BSSw0pts),
            ('CLAS', 'BSSw_c1', C_BSSw1pts),
            ('Hall A', 'BSDw_s1', H_BSDwpts), ('Hall A', 'BSSw_c0', H_BSSw0pts),
            ('Hall A', 'BSSw_c1', H_BSSw1pts),
            ('Hall A', 'BSD_s1', H_BSDpts), ('Hall A', 'BSS_c0', H_BSS0pts),
            ('Hall A', 'BSS_c1', H_BSS1pts)
            ]
    sets[8] = [('HERMES', 'ALUI', ALUIpts),
            ('HERMES', 'BCA', BCApts), ('CLAS', 'BSA', CLASpts),
            ('Hall A', 'BSDw', BSDwpoints), ('Hall A', 'BSSw', BSSwpoints),
            ('HERMES', 'AUTI', AUTIpts)]
    names = [th.name[:10] for th in ths]
    sublines = ['------' for th in ths]
    if out == 'chis':
        # We want chisq/npts
        ftit = 21*' ' + len(names)*'{:^10s}'
        fstr = '{:9s} {:7s}: ' + len(names)*'{:10.2f}' + '   (np ={dof:3d})'
        chi_ind = 0
    else:    
        # We want probabilities
        ftit = 20*' ' + len(names)*'{:<10s}'
        fstr = '{:9s} {:7s}:  ' + len(names)*'{:<10.3g}' + '   (np ={dof:3d})'
        chi_ind = 2
    print ftit.format(*names)
    print ftit.format(*sublines)
    total_chis = np.array([0. for th in ths])
    total_npts = 0
    for collab, obs, pts in sets[nsets]:
        cutpts = select(pts, criteria=['Q2>=%f' % Q2cut, 'Q2<=%f' % Q2max])
        npts = len(cutpts)
        if out == 'pulls':
            quals = [th.chisq(cutpts, sigmas=True).sum()/np.sqrt(npts) for th in ths]
        else:
            chis = [th.chisq(cutpts)[chi_ind] for th in ths]
            total_chis += np.array(chis)
            total_npts += npts
            quals = [chi/npts for chi in chis]
        print fstr.format(collab, obs, *quals, dof=npts)
    if out == 'chis':
        # version with chisq/npts:
        total_chis = total_chis/total_npts
        print ftit.format(*sublines)
        print fstr.format('===', 'TOTAL', *total_chis.tolist(), dof=total_npts)



class hubDict(dict):
    """Merges two dictionaries, but not actually but just by forwarding."""

# most of the methods below are probably not necessary

    def __init__(self, da, db):
        self.d1 = da
        self.d2 = db

    def __getitem__(self, name):
        if self.d1.has_key(name):
            return self.d1[name]
        else:
            return self.d2[name]

    def __setitem__(self, name, value):
        if self.d1.has_key(name):
            self.d1[name] = value
        else:
            self.d2[name] = value

    def __iter__(self):
        return itertools.chain(self.d1.__iter__(), self.d2.__iter__())

    def has_key(self, name):
        if self.d1.has_key(name) or self.d2.has_key(name):
            return True
        else:
            return False

    def keys(self):
        return self.d1.keys() + self.d2.keys()

    def items(self):
        return self.d1.items() + self.d2.items()

    def iteritems(self):
        return itertools.chain(self.d1.iteritems(), self.d2.iteritems())

    def iterkeys(self):
        return itertools.chain(self.d1.iterkeys(), self.d2.iterkeys())

    def itervalues(self):
        return itertools.chain(self.d1.itervalues(), self.d2.itervalues())

    def copy(self):
        print "Can't copy hubDict yet!!!"

    def update(self, d):
        for key in d:
            self.__setitem__(key, d[key])

    def popitem(self):
        try:
            return self.d1.popitem()
        except KeyError:
            return self.d2.popitem()

    def __repr__(self):
        return 'First: %s\nSecond: %s' % (
                self.d1.__repr__(), self.d2.__repr__())


def _fakecolor(a, b):
    return a

try:
	from termcolor import colored
except ImportError:
	colored = _fakecolor

def stringcolor(a, c, colors=False):
    if colors:
        return colored(a, c)
    else:
        return _fakecolor(a, c)


def FTF(data, testdata=None, cosmax=None, sinmax=None, inverse=False):
    """Fourier series fit to data (takes and returns pandas dataframe).

    cosmax, sinmax  - index of highest harmonics
    inverse = True  - data values ARE harmonics, calculate function(phi)

    """
    N = len(data)
    # If length of series is not specified by the user, use maximal one
    #   which makes this fit equal to DFT
    if not cosmax and not (cosmax == 0):
        cosmax = int((N-1)/2.)  # maximal cos harmonic
    if not sinmax and not (sinmax == 0):
        sinmax = cosmax + ((N-1) % 2)  # maximal sin harmonic
    nharm = 1 + cosmax + sinmax  # number of harmonics
    if isinstance(testdata, pd.DataFrame):
        A = np.array([[np.cos(n*phi) for n in range(cosmax+1)] + [np.sin(n*phi) for n in range(1,sinmax+1)] for phi in testdata.phi.values])
    else:
        # take phis from first argument data
        A = np.array([[np.cos(n*phi) for n in range(cosmax+1)] + [np.sin(n*phi) for n in range(1,sinmax+1)] for phi in data.phi.values])
    B = data.val.values
    if not inverse:
        res = np.linalg.lstsq(A, B)
        #print "Sum of residuals = {}".format(res[1])
        z = np.zeros(N-len(res[0]), dtype=A.dtype)
        vals = np.concatenate((res[0], z), axis=0)
        df = pd.DataFrame({'phi': data.phi.values, 'val': vals})
        #print "Number of harmonics = {}".format(nharm)
    else:
        if isinstance(testdata, pd.DataFrame):
            res = np.dot(A, B[:nharm])
            df = pd.DataFrame({'phi': testdata.phi.values, 'val': testdata.val.values, 'pred': res})
        else:
            # we invert on phis of first argument
            res = np.dot(A, B[:nharm])
            df = pd.DataFrame({'phi': data.phi.values, 'val': res})
    return df


def FTFMC(data, nsamples=100, cosmax=None, sinmax=None, inverse=False):
    """Fourier series fit to data with MC error propagation (takes and returns pandas dataframe).
    
    cosmax, sinmax  - index of highest harmonics
    inverse = True  - data values ARE harmonics, calculate function(phi)
    
    """
    N = len(data)
    # If length of series is not specified by the user, use maximal one
    #   which makes this fit equal to DFT
    if not cosmax and not (cosmax == 0):
        cosmax = int((N-1)/2.)  # maximal cos harmonic
    if not sinmax and not (sinmax == 0):
        sinmax = cosmax + ((N-1) % 2)  # maximal sin harmonic
    nharm = 1 + cosmax + sinmax  # number of harmonics
    A = np.array([[np.cos(n*phi) for n in range(cosmax+1)] + [np.sin(n*phi) for n in range(1,sinmax+1)] for phi in data.phi.values])
    B = np.transpose((np.ones((nsamples,N))*data.val.values + np.random.randn(nsamples,N)*data.err.values))  # replicas
    if not inverse:
        res = np.linalg.lstsq(A, B)
        z = np.zeros((nsamples, N-len(res[0])), dtype=A.dtype)
        vals = np.concatenate((res[0].T, z), axis=1)
        df = pd.DataFrame({'phi': data.phi.values, 'val': vals.mean(axis=0), 'err': vals.std(axis=0)})
        #print "Number of harmonics = {}".format(nharm)
    else:
        res = np.dot(A, B[:nharm])
        df = pd.DataFrame({'phi': data.phi.values, 'val': res.mean(axis=1), 'err': res.std(axis=1)})
    return df


def cvsets(df_in, nfolds=3, shuffle=True):
        """Return n-fold cross-validation sets [(train1, valid1), (train2, valid2), ...].
        
        Copies data.
        
        """
        if shuffle:
            df = df_in.reindex(np.random.permutation(df_in.index))
        else:
            df = df_in.copy()
        # stolen from sklearn
        n = len(df)
        fold_sizes = (n // nfolds) * np.ones(nfolds, dtype=np.int)
        fold_sizes[:n % nfolds] += 1
        current = 0
        chunks = []
        for fold_size in fold_sizes:
            start, stop = current, current + fold_size
            #yield obj.idxs[start:stop]
            chunks.append(df[start:stop])
            current = stop
        sets = []
        # my ugly coding
        for f in range(nfolds):
            inds = range(nfolds)
            k = inds.pop(f)
            train = chunks[inds[0]]
            for i in inds[1:]:
                train = train.append(chunks[i])
            sets.append([train, chunks[f]])
        return sets


def FTanalyse(bins, HMAX=2, nf=4, Nrep=1):
    """Determine the number of harmonics present in bins according to n-fold cross-validation.

       bins - dictionary with bins
       HMAX - highest harmonic used in searches
         NS - number of replicas for MC error propagation
       Nrep - number of CV repetitions (probably wrong to make > 1)
         nf - number of folds for cross-validation

    """

    mins = []
    for nn in range(Nrep):
        for k in bins:
            df = bins[k]
            err_min = 100
            for CM in range(HMAX+1):
                for SM in range(HMAX+1):
                    errs = []
                    for cvtrain, cvtest in cvsets(df, nfolds=nf):
                        cvtest.reset_index(drop=True, inplace=True)
                        dfFT = FTF(cvtrain, cosmax=CM, sinmax=SM)
                        dftest = FTF(dfFT, testdata=cvtest, cosmax=CM, sinmax=SM, inverse=True)
                        errs.append(np.sum((dftest.val.values-dftest.pred.values)**2))
                    errs = np.array(errs)
                    # We divide by npts and not ndof because test set
                    # is NOT used for fitting, so there should be no penalization
                    # of model's complexity
                    err = errs.mean() / len(cvtest)
                    # However, experiments show that small penalty for complexity
                    # sometimes improves accuracy (but this is quite simple
                    # and creates problems with large number of folds)
                    #err = errs.mean() /(len(cvtest)-np.sqrt(CM+SM+1))
                    #print CM, SM, err
                    if err <= err_min:
                        err_min = err
                        h_min = (CM, SM)
            mins.append(h_min)
    nc, ns = np.array(mins).mean(axis=0)
    delc, dels = np.array(mins).std(axis=0)
    print "Highest extractable cos harmonic = {:.3f} +- {:.3f}".format(nc, delc)
    print "Highest extractable sin harmonic = {:.3f} +- {:.3f}\n".format(ns, dels)
    return int(np.round(nc)), int(np.round(ns))

def describe_data(pts):
    """Print observables and where they come from."""
    all = []
    print "{:2s} x {:5s}  {:6s}  {:4s}   {:3s} {:12s}".format(
     'npt', 'obs', 'collab', 'FTn', 'id', 'ref.')
    print 45*'-'
    #print "{:2s} x {:5s}  {:6s}  {:4s}".format(
     #'npt', 'obs', 'collab', 'FTn')
    #print 30*'-'
    for pt in pts:
        props = []
        for prop in ['y1name', 'collaboration', 'FTn', 'id', 'reference']:
        #for prop in ['y1name', 'collaboration', 'FTn']:
            if hasattr(pt, prop):
                if prop=='y1name' and pt.y1name=='X' and pt.has_key('t'):
                    props.append('Xt')
                else:
                    props.append(str(getattr(pt,prop)))
            else:
                props.append('N/A')
        all.append(tuple(props))
    tot = len(all)
    uniqs = set(all)
    cc = 0
    for uniq in sorted(uniqs):
        n = all.count(uniq)
        cc += n
        print "{:2d} x {:5s}  {:6s}  {:4s}   {:3s} {:12s}".format(n, *uniq)
        #print "{:2d} x {:5s}  {:6s}  {:4s}".format(n, *uniq)
    assert cc == tot
    print 45*'-'
    print "TOTAL = {}".format(tot)
    return tot

def _hmax(amp, trig):
    """Index of largest harmonic"""
    if amp == 'DVCS':
        return 2
    if amp == 'BH':
        if trig == 'c':
            return 2
        else:
            return 1
    return 3

def compare_harmonics(th, pt, pol='unp'):
    """Print relative sizes of nonzero harmonics."""
    amps = []
    for amp in ['BH', 'INT', 'DVCS']:
        prefac = getattr(th, 'PreFac'+amp)(pt)
        for trig in ['c','s']:
            for k in range(int(trig=='s'), _hmax(amp,trig)+1):
                name = trig+amp+str(k)+pol
                try:
                    res = prefac*getattr(th,name)(pt)
                    if abs(res) > 0:
                        amps.append((abs(res), np.sign(res), name))
                except AttributeError:
                    # amps.append(('N/A', name))
                    pass
    amps = sorted(amps, reverse=True)
    max = amps[0][0]
    for value, sign, name in amps:
        print "{:9s} = {: .3g}".format(name, sign*value/max)
