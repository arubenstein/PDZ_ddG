#!/usr/bin/env python

import os
import argparse

def read_csv_dict(csv_file_name):
    with open(csv_file_name) as c:
        lines = c.readlines()

    csv_dict = dict( (l.strip().split(',')[0], l.strip().split(',')[1:] ) for l in lines )

    return csv_dict

def read_csv_list(csv_file_name, protocol=None):
    with open(csv_file_name) as c:
        lines = c.readlines()

    print lines

    csv_list = [ l.strip().split(',') for l in lines ]
    print csv_list
    if protocol == "rosetta":
        csv_list = [ [l[0][0:-5]].extend(l[1:]) for l in csv_list ]
    elif protocol == "amber":
        csv_list = [ [l[0][38:-4]].extend(l[1:]) for l in csv_list ]

    print csv_list
    
    return csv_list

def main(rec_corr, amber_csv_path, rosetta_csv_path, amber_pdb_path, rosetta_pdb_path, out_csv_path):
    #root, ext = os.path.splitext(in_pdb)	 
    #pdb_path, pdb_id = os.path.split(root)

    rec_corr_list = read_csv_list(rec_corr)
    print rec_corr_list
    amber_csv_dict = {}

    #read in all amber csvs that correspond to column 3 in rec_corr file and rosetta ones too
    for record_id, prefix, filename in rec_corr_list:
        amber_csv_dict[filename] = read_csv_list(os.path.join(amber_csv_path,filename+".csv"), protocol="amber") 
        rosetta_csv_dict[filename] = read_csv_list(os.path.join(rosetta_csv_path,filename+".csv"), protocol="rosetta")

    print amber_csv_dict
    print rosetta_csv_dict

    #find wt csv that correspond to wt row in rec_corr_file (column 2)

    #find average of ddgs in wt csv as well as lowest rmsd and lowest energy (second and third columns)

    #loops thru other records in rec_corr_dict
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument ('--rec_corr', help="record correspondence file")
    parser.add_argument ('--amber_csv_path', help="amber csv dg file")
    parser.add_argument ('--rosetta_csv_path', help="rosetta csv dg file")
    parser.add_argument ('--amber_pdb_path', help="In path for Amber pdbs")
    parser.add_argument ('--rosetta_pdb_path', help="In path for Rosetta pdbs")
    parser.add_argument ('--out_csv_path', help="Out path for final csv files")
    parser.add_argument ('--out_pdb_path', help="Out path for final pdb files")

    args = parser.parse_args()

    main(args.rec_corr, args.amber_csv_path, args.rosetta_csv_path, args.amber_pdb_path, args.rosetta_pdb_path, args.out_csv_path) 
