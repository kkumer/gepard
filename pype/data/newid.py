#!/usr/bin/env python2

import os, re, string, sys
import utils

# files with ".dat" extension
#datafiles = [f for f in os.listdir(os.path.curdir) if f.split('.')[-1]=="dat"]
datafiles = utils.listFiles('data', '*dat') 


ids = []
for file in datafiles:
    dataFile = open(file, 'r')
    dataFileLine = dataFile.readline()
    while dataFileLine:
        # remove comments
        dataFileLine = dataFileLine.split('#')[0]
        # look for line with 'id ='
        if re.search(r'id( )*=', dataFileLine):
            # converting preamble line into dictionary item
            id = string.split(dataFileLine,"=")[-1].strip()
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

print("Used ids are: ")
print(ids)
print("So first unused is: [ %d ]" % (ids[-1]+1,))
