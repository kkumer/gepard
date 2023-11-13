#!/usr/bin/env python

import os
import re
import string
import sys


def procdir(dirname, datafiles):
    """Processing single directory and its files."""
    ids = []
    for file in datafiles:
        dataFile = open(os.path.join(dirname, file), 'r')
        dataFileLine = dataFile.readline()
        while dataFileLine:
            # remove comments
            dataFileLine = dataFileLine.split('#')[0]
            # look for line with 'id ='
            if re.search(r'id( )*=', dataFileLine):
                # converting preamble line into dictionary item
                id = dataFileLine.split("=")[-1].strip()
                try:
                    if ids.count(int(id)):
                        if int(id)>30:
                            sys.stderr.write(
                             'Duplicate id (%d) of file %s.\n' % (int(id), file))
                    elif int(id) > 1000:
                        sys.stderr.write(
                          'Mock/pseudo data (id=%d) in file %s.\n' % (int(id), file))
                    else:
                        ids.append(int(id))
                    break
                except ValueError:
                    sys.stderr.write('File %s has funny id.\n' % file)
                break
            dataFileLine = dataFile.readline()

    ids.sort()
    if len(ids) > 0:
        print("-----------\n[{}] used ids: {}".format(dirname, ids))
        print("[{}] first unused: {}".format(dirname, ids[-1]+1))


for root, d_names, f_names in os.walk(os.path.curdir):
    datafiles = [f for f in f_names if f.split('.')[-1]=="dat"]
    procdir(root, datafiles)
