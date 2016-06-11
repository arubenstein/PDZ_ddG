#!/usr/bin/env python

import sys, shutil

pdblist = open(sys.argv[1]).read().split()

for pdb in pdblist:
    print("copying", pdb)
    shutil.copy(pdb + '/mmgbsa/ddG.csv', 'ddG/' + pdb + '.csv')
