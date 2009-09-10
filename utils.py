""" 
Some utility functions needed all over the application 

loaddata -- loads datafiles from directory
parse -- parses datafiles
subplot -- creates subplot for matplotlib plot
npars -- number of free (not fixed) parameters of minuit object
str2num -- transforms string to float or int
prettyprint -- formatted printout of numbers
"""

import os, re, string

import Data
from constants import toTeX

def loaddata(datadir='data'):
    """Return dictionary {id : `DataSet`, ...}  out of datadir/*dat files."""

    data = {}
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] == ".dat":
            dataset = Data.DataSet(datafile=os.path.join(datadir, file))
            data[dataset.id] = dataset
    return data

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
            numbers = re.findall(r'[-\.\d]+', dataFileLine)
            data.append(map(float, numbers))
        dataFileLine = dataFile.readline()

    return desc, data


def subplot(ax, dataset, xaxis=None, kinlabels=[], fits=[]):
    """Plot datapoints together with fit/theory line(s).

    ax -- subplot of matplotlib's figure i.e. ax = figure.add_subplot(..)
    dataset -- `DataSet` instance to be plotted
    xaxis -- abscissa variable; if None, last of dataset.xaxes is taken
    kinlabels -- list of constant kinematic variables whose values will
                 be put on plott
    fits -- list of parameter sets/tuples describing fit curves for
            plotting. 

    TODO: - fits should take some sort of model specifications
          - fit should be real line i.e. not evaluated just at datapoints
          - coloring second, third etc. dataset and fit differently
          - different treatment of stat and syst errors
          - multiple datasets on the same panel

    """
    # first, fix the input if needed
    if not isinstance(kinlabels, list): kinlabels = [kinlabels]
    if not isinstance(fits, list): fits = [fits]
    if not xaxis: xaxis = dataset.xaxes[-1]
    xval = []; yval = []; yerr = []
    for pt in dataset:
        xval.append(getattr(pt, xaxis)) 
        yval.append(pt.val)
        yerr.append(pt.err)
    ax.errorbar(xval, yval, yerr, linestyle='None', marker='o', color='b')
    # Fit lines
    shapes = ['s', '^', 'd', 'h']  # first squares, then triangles, diamonds, hexagons
    colors = ['red', 'green', 'brown', 'purple']  # squares are red, etc.
    fitn = 0
    for (approach, pars) in fits:
        # take abscissae from dataset
        line = [getattr(approach, pt.yaxis)(pt, pars) for pt 
                in dataset]
        #ax.plot(xval, line, 'r-')  # join symbols by line
        ax.plot(xval, line, shapes[fitn], markersize=5, 
                markerfacecolor=colors[fitn], markeredgecolor='black')
        fitn += 1
    # axes labels
    ax.set_xlabel(toTeX[xaxis])
    ax.set_ylabel(toTeX[dataset.yaxis])
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

