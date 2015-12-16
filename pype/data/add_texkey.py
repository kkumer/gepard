#!/usr/bin/env python2
'''Go over files specified on command line and add some content.
This particular example adds texkeys after reference key.
'''

import fileinput, re, subprocess

for line in fileinput.input(inplace=1): #, backup='.bak'):
    print line,
    m = re.search(r'(([A-Za-z]+:[12][0-9][0-9][0-9][a-z]+)|((hep|gr|nucl|quant)-(ph|th|ex|lat|qc)/[0-9\.]{7})|([0-9]{4}\.[0-9]{4,5}))', line)
    # If regex is found add next line:
    if m:
        # getbib is external program returning texkey given arXiv id
        key = subprocess.check_output(['getbib', m.groups()[0]]).strip()
        print 'inspiretex = {}'.format(key)

