import sys
import pytraj as pt

native_dir = '../../../../ancillary_files/input_pdbs/'

def main():
    pdblist = open(sys.argv[1]).read().split()
    csv_files = [pdb + '.csv' for pdb in pdblist]
    
    for fn in csv_files:
        print("processing {}".format(fn))
        native_pdb = native_dir + fn[:6] + '_complex.pdb'
        print(native_pdb)
        ref = pt.iterload(native_pdb)
    
        root = '../{}/'.format(fn.split('.')[0])
        with open(fn) as fh:
            rst7_files = [root + line.split(',')[0].strip('../') for line in fh.readlines()]
    
        traj = pt.iterload(rst7_files, root + 'prmtop')
        ca_rmsd = pt.rmsd(traj, ref=ref, mask='@CA')
    
        with open(fn) as fh:
            lines = [','.join((line.strip('\n'), str(rmsd))) for line, rmsd in zip(fh.readlines(), ca_rmsd)]
    
        with open('output/' + fn, 'w') as fw:
            fw.write('\n'.join(lines))

main()
