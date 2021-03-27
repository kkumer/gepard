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
stringcolor -- coloring string for output, if possible
describe_data -- List observables in dataset.
"""

import fnmatch
import os
import re

import numpy as np

import gepard as g
from gepard.constants import Mp, Mp2


class KinematicsError(Exception):
    pass


def loaddata(datadir='data', approach=False):
    """Return dictionary {id : `DataSet`, ...}  out of datadir/*dat files.
    
    approach defines conventions for frame, kinematics etc. to which data
    is adapted.
    
    """
    data = {}
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] == ".dat":
            dataset = g.data.DataSet(datafile=os.path.join(datadir, file))
            if approach and dataset.process in ['ep2epgamma', 'en2engamma']:
                [pt.to_conventions(approach) for pt in dataset]
                [pt.prepare(approach) for pt in dataset]
            data[dataset.id] = dataset
    return data

def _complete_xBWQ2(kin):
    """Make trio {xB, W, Q2} complete if two of them are given in 'kin'."""
    if 'W' in kin and 'Q2' in kin and 'xB' not in kin:
        kin.xB = kin.Q2 / (kin.W**2 + kin.Q2 - Mp2)
    elif 'xB' in kin and 'Q2' in kin and 'W' not in kin:
        kin.W = np.sqrt(kin.Q2 / kin.xB - kin.Q2 + Mp2)
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
        raise KinematicsError('Overdetermined set {xB, W, Q2} given.')
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
        dataFileLine = dataFile.readline()

    return desc, data

def npars(m):
    """Return number of free (not fixed) parameters of MINUIT object m."""

    n = 0
    for key in m.parameter_names:
        if 'fix_'+key in m.parameters.keys() and not m.parameters['fix_'+key]:
            n += 1
    return n

def str2num(s):
    """Convert string to number, taking care if it should be int or float.
    
    http://mail.python.org/pipermail/tutor/2003-November/026136.html
    """

    if "." in s:
        return float(s) 
    else:
        return int(s)


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
        print('%s%s' % (padding, s))
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
    """Return DataSet of DataPoints satisfying criteria.

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
    tmp = g.data.DataSet(selected)
    tmp.__dict__ = dataset.__dict__.copy() # transfer the attributes
    return tmp

def listdb(db):
    print("%-22s--+--%s" % (22*'-', 60*'-'))
    print("%-22s  |  %s" % ('name', 'description'))
    print("%-22s--+--%s" % (22*'-', 60*'-'))
    for key in db:
        print("%-22s  |  %s" % (key, db[key].description))
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
            print('[%3i] %8s %3i %9s %10s %s' % (dt.id, dt.collaboration, len(dt), dt.y1name, ref, ref2))
        except KeyError:
            pass

def listchis(ths, Q2cut=1., Q2max=1.e3, nsets=10, out='chis'):
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
            H17_BSDwpts, H17_BSSw0pts, H17_BSSw1pts,\
            H_BSDpts, H_BSS0pts, H_BSS1pts,\
            C_BSDwpts, C_BSSw0pts, C_BSSw1pts, H_BSD, H_BSS,\
            C_BSD, C_BSS, H20_nBSSw0pts, H20_nBSSw1pts
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
               ('CLAS15_KK', 'BTSA', CLAS14BTSApts)]
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
    sets[9] = [
            ('CLAS', 'BSA', CLAS14BSApts),
            ('CLAS', 'TSA', CLAS14TSApts), 
            ('CLAS', 'BTSA', CLAS14BTSApts),
            ('CLAS', 'BSDw_s1', C_BSDwpts), 
	    ('CLAS', 'BSSw_c0', C_BSSw0pts),
            ('CLAS', 'BSSw_c1', C_BSSw1pts),
            ('HallA 15', 'BSDw_s1', H_BSDwpts), 
	    ('HallA 15', 'BSSw_c0', H_BSSw0pts),
            ('HallA 15', 'BSSw_c1', H_BSSw1pts),
	    ('HallA 17', 'BSSw_c0', H17_BSSw0pts),
            ('HallA 17', 'BSSw_c1', H17_BSSw1pts)
            ]
    sets[10] = [
            ('H1ZEUS', 'X_DVCS', H1ZEUS), 
            ('HERMES', 'ALUI', ALUIpts),
            ('HERMES', 'BCA', BCApts),
            ('HRM/CLS', 'AUL', AULpts), 
	    ('HERMES', 'ALL', ALLpts),
            ('HERMES', 'AUTI', AUTIpts),
            ('CLAS', 'BSA', CLAS14BSApts),
            ('CLAS', 'TSA', CLAS14TSApts), 
            ('CLAS', 'BTSA', CLAS14BTSApts),
            ('CLAS', 'BSDw_s1', C_BSDwpts), 
	    ('CLAS', 'BSSw_c0', C_BSSw0pts),
            ('CLAS', 'BSSw_c1', C_BSSw1pts),
            ('HallA 15', 'BSDw_s1', H_BSDwpts), 
	    ('HallA 15', 'BSSw_c0', H_BSSw0pts),
            ('HallA 15', 'BSSw_c1', H_BSSw1pts),
            ('HallA 17', 'BSDw_s1', H17_BSDwpts), 
	    ('HallA 17', 'BSSw_c0', H17_BSSw0pts),
            ('HallA 17', 'BSSw_c1', H17_BSSw1pts)
            ]
    sets[11] = [
            ('CLAS', 'BSA', CLAS14BSApts),
            ('CLAS', 'TSA', CLAS14TSApts), 
            ('CLAS', 'BTSA', CLAS14BTSApts),
            ('CLAS', 'BSDw_s1', C_BSDwpts), 
	    ('CLAS', 'BSSw_c0', C_BSSw0pts),
            ('CLAS', 'BSSw_c1', C_BSSw1pts),
            ('HallA 15', 'BSDw_s1', H_BSDwpts), 
	    ('HallA 15', 'BSSw_c0', H_BSSw0pts),
            ('HallA 15', 'BSSw_c1', H_BSSw1pts),
	    ('HallA 17', 'BSSw_c0', H17_BSSw0pts),
            ('HallA 17', 'BSSw_c1', H17_BSSw1pts),
	    ('HallA 20', 'nBSSw_c0', H20_nBSSw0pts),
	    ('HallA 20', 'nBSSw_c1', H20_nBSSw1pts)
            ]
    names = [th.name[:10] for th in ths]
    sublines = ['------' for th in ths]
    if out == 'chis' or out == 'pulls':
        # We want chisq/npts
        ftit = 21*' ' + len(names)*'{:^10s}'
        fstr = '{:9s} {:7s}: ' + len(names)*'{:10.2f}' + '   (np ={dof:3d})'
        chi_ind = 0
    else:    
        # We want probabilities
        ftit = 20*' ' + len(names)*'{:<10s}'
        fstr = '{:9s} {:7s}:  ' + len(names)*'{:<10.3g}' + '   (np ={dof:3d})'
        chi_ind = 2
    print(ftit.format(*names))
    print(ftit.format(*sublines))
    total_chis = np.array([0. for th in ths])
    total_npts = 0
    for collab, obs, pts in sets[nsets]:
        cutpts = select(pts, criteria=['Q2>=%f' % Q2cut, 'Q2<=%f' % Q2max])
        npts = len(cutpts)
        if out == 'pulls':
            chis = [th.chisq(cutpts, pull=True)[chi_ind] for th in ths]
        else:
            chis = [th.chisq(cutpts)[chi_ind] for th in ths]
        total_chis += np.array(chis)
        total_npts += npts
        if out == 'pulls':
            quals = [chi/np.sqrt(npts) for chi in chis]
        else:
            quals = [chi/npts for chi in chis]
        print(fstr.format(collab, obs, *quals, dof=npts))
    if out == 'chis':
        # version with chisq/npts:
        total_chis = total_chis/total_npts
        print(ftit.format(*sublines))
        print(fstr.format('===', 'TOTAL', *total_chis.tolist(), dof=total_npts))


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


def describe_data(pts):
    """Print observables and where they come from."""
    all = []
    print("{:2s} x {:5s}  {:6s}  {:4s}   {:3s} {:12s}".format(
     'npt', 'obs', 'collab', 'FTn', 'id', 'ref.'))
    print(45*'-')
    #print "{:2s} x {:5s}  {:6s}  {:4s}".format(
     #'npt', 'obs', 'collab', 'FTn')
    #print 30*'-'
    for pt in pts:
        props = []
        for prop in ['y1name', 'collaboration', 'FTn', 'id', 'reference']:
        #for prop in ['y1name', 'collaboration', 'FTn']:
            if hasattr(pt, prop):
                if prop=='y1name' and pt.y1name=='X' and 't' in pt:
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
        print("{:2d} x {:5s}  {:6s}  {:4s}   {:3s} {:12s}".format(n, *uniq))
        #print "{:2d} x {:5s}  {:6s}  {:4s}".format(n, *uniq)
    assert cc == tot
    print(45*'-')
    print("TOTAL = {}".format(tot))
    return tot

