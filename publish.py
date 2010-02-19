#!/usr/bin/env python

# This script processes all .dat files in the subdirectory
# 'datdir' of the current directory
# and creates corresponding .html files according to the
# datadir/expdata.tmpl and datadir/index.tmpl  cheetah templates
# as well as PNG figures with data.
# To see result, open datadir/index.html
#
# kkumer@phy.hr 2009-09-04 

import re, os, sys, shutil


from Cheetah.Template import Template
import matplotlib
matplotlib.use('AGG')

import Data 

datadir = 'data' # directory with exp data
wwwdir = 'www' # directory with templates and HTML files
outdatadir = 'outdata' # subdir of wwwdir for WWW-publishing of data
figdir = 'figs' # subdir of wwwdir for WWW-publishing of data
indexfilename = 'index.html'
version = '2009-09-04'

def cleanup(excludelist, dirname, fnames):
    """Remove all files (not directories) in dirname, apart
    from those in excludelist."""

    for item in excludelist:
        try:
            fnames.remove(item)
        except ValueError:
            pass
    for item in fnames:
        target = os.path.join(dirname, item)
        if not os.path.isdir(target): # don't remove dirs
            os.remove(target)

if len(sys.argv) > 1: 
    if sys.argv[1] == '-c':  #cleanup
        # Remove all files in 'www' dir, apart from excludelist
        excludelist = ['.svn', 'index.tmpl', 'expdata.tmpl']
        os.path.walk('www', cleanup, excludelist)
        sys.exit()
    else:
        sys.stderr.write('Unknown argument. Usage: publish.py [-c]\n')
        sys.exit(1)

# Select only *.dat files from the datadir for processing
datafiles = []
for file in os.listdir(datadir):
    if re.match(r'ep2epgamma.*.dat$', file):
        datafiles.append(file)

# Override by simpler sets for debugging:
#datafiles = datafiles[0:4]
#datafiles = ['ep2epgamma-BCA-HERMES-08-cos1_b.dat']

datasets = [] # list of DataSet instances

for file in datafiles:
    set = Data.DataSet(datafile = os.path.join(datadir, file))
    sys.stderr.write(set.filename + '\n')
    datasets.append(set)
    set.filepath = os.path.join(outdatadir, set.filename)
    set.figfilename = set.filename + '.html'
    set.figpath = os.path.join(figdir, set.filename + '.png')

    # apply to template and print to *.html
    f = open(os.path.join(wwwdir, set.figfilename), 'w')
    f.write(str(Template(file = os.path.join(wwwdir, 'expdata.tmpl'), searchList=[set])))
    f.close()

    # put actual data u data subdir
    shutil.copy(os.path.join(datadir, file), os.path.join(wwwdir, outdatadir))

    # plot fig in fig subdir
    set.plot(set.xaxes[-1], path=os.path.join(wwwdir, figdir))


# do index page
f = open(os.path.join(wwwdir, indexfilename), 'w')
datasets.sort(lambda s1, s2: int(s1.id)-int(s2.id))
indexpage = {'datasets' : datasets, 'version' : version}
f.write(str(Template(file = os.path.join(wwwdir, 'index.tmpl'), searchList=[indexpage])))
f.close()

