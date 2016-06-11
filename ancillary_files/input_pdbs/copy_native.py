#!/usr/bin/env python

import os, shutil

pdblist = open("../../results/refine/amber/pdblist_extra.txt").read().split()

for pdb in pdblist:
    fn = pdb[:6]

    src = fn + '_0001_complex.pdb'
    dest = fn + '_complex.pdb'
    if os.path.exists(src):
        print(src, dest)
        shutil.copy(src, dest)
