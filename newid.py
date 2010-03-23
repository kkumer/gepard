#!/usr/bin/env python

import os, re, string
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
            ids.append(int(id))
            break
        dataFileLine = dataFile.readline()

ids.sort()

print "Used ids are: "
print ids
print "So first unused is: [ %d ]" % (ids[-1]+1,)
