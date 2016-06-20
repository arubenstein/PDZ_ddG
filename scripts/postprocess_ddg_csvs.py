#!/usr/bin/env python

import glob
import os
import argparse
import numpy as np
from plot import scatterplot
from plot import conv
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

def read_csv_dict(csv_file_name):
    with open(csv_file_name) as c:
        lines = c.readlines()

    csv_dict = dict( (l.strip().split(',')[0], l.strip().split(',')[1:] ) for l in lines )

    return csv_dict

def read_csv_list(csv_file_name, protocol=None):
    with open(csv_file_name) as c:
        lines = c.readlines()


    csv_list = [ l.strip().split(',') for l in lines ]
    if protocol == "rosetta":
        csv_list = [ [l[0][0:-5]] + l[1:] for l in csv_list ]
    elif protocol == "amber":
        csv_list = [ [l[0][38:-4]] + [float(l[1])/2.0] + l[2:]  for l in csv_list ]

    
    return csv_list

def get_stats(csv_dict):
    mean_dg = np.mean([ float(dg) for key, dg, rmsd in csv_dict ])
    min_dg_dg = min([ float(dg) for key, dg, rmsd in csv_dict ])
    min_rmsd_dg = float(sorted(csv_dict, key=lambda x: float(x[2]))[0][1])
    bottom3_dg = np.mean(sorted([float(dg) for key,dg,rmsd in csv_dict])[0:3])
    return mean_dg, min_dg_dg, min_rmsd_dg, bottom3_dg

def main(rec_corr_path, amber_csv_path, rosetta_csv_path, amber_pdb_path, rosetta_pdb_path, out_csv_path):
    #root, ext = os.path.splitext(in_pdb)	 
    #pdb_path, pdb_id = os.path.split(root)

    fig_all, axarr_all = conv.create_ax(4, 3)
    fig_all_corr, axarr_all_corr = conv.create_ax(1, 2)

    list_rec_corr_names = glob.glob(rec_corr_path + "*.rc")

    corr_values_dict = { "Mean" : [ [] for _ in range(0,6) ], "Min DG" : [ [] for _ in range(0,6) ], "Min RMSD" : [ [] for _ in range(0,6) ], "Bottom 3" : [ [] for _ in range(0,6) ] }
    all_amean_ddg = []
    all_rmean_ddg = []
    k_ddg = []

    for rec_corr in list_rec_corr_names:
	print rec_corr
        rec_corr_list = read_csv_list(rec_corr)
        if len(rec_corr_list[0]) == 3:
            continue
	amber_csv_dict = {}
        rosetta_csv_dict = {}

        #read in all amber csvs that correspond to column 3 in rec_corr file and rosetta ones too
        for record_id, prefix, filename, known_ddg in rec_corr_list:
            amber_csv_dict[filename] = read_csv_list(os.path.join(amber_csv_path,filename+".csv"), protocol="amber") 
            rosetta_csv_dict[filename] = read_csv_list(os.path.join(rosetta_csv_path,filename+".csv"), protocol="rosetta")

        #find wt csv that correspond to wt row in rec_corr_file (column 2)
        wt_csv_name = [ rec[2] for rec in rec_corr_list if "wt" in rec[1] ][0]
        wt_a_csv_dict = amber_csv_dict[wt_csv_name]
        wt_r_csv_dict = rosetta_csv_dict[wt_csv_name]

        #find average of ddgs in wt csv as well as lowest rmsd and lowest energy (second and third columns)
        amean_dg, amin_dg_dg, amin_rmsd_dg, abottom3_dg = get_stats(wt_a_csv_dict)
        rmean_dg, rmin_dg_dg, rmin_rmsd_dg, rbottom3_dg = get_stats(wt_r_csv_dict)

        known_ddg = []
        a_mean_ddg = []
        a_min_ddg = []
        a_minr_ddg = []
        a_bottom3_ddg = []
        r_mean_ddg = []
        r_min_ddg = []
        r_minr_ddg = []
        r_bottom3_ddg = []

        #loops thru other records in rec_corr_dict
        for rec, prefix, filename, k in rec_corr_list:
            if "wt" not in prefix:
                a_amean_dg, a_amin_dg_dg, a_amin_rmsd_dg, a_abottom3_dg = get_stats(amber_csv_dict[filename])
                r_rmean_dg, r_rmin_dg_dg, r_rmin_rmsd_dg, r_rbottom3_dg = get_stats(rosetta_csv_dict[filename])

                known_ddg.append(float(k))
	        a_mean_ddg.append(a_amean_dg-amean_dg)
 	        a_min_ddg.append(a_amin_dg_dg-amin_dg_dg)
                a_minr_ddg.append(a_amin_rmsd_dg-amin_rmsd_dg)
                a_bottom3_ddg.append(a_abottom3_dg-abottom3_dg)
                r_mean_ddg.append(r_rmean_dg-rmean_dg) 
	        r_min_ddg.append(r_rmin_dg_dg-rmin_dg_dg)
                r_minr_ddg.append(r_rmin_rmsd_dg-rmin_rmsd_dg)
                r_bottom3_ddg.append(r_rbottom3_dg-rbottom3_dg)

        all_amean_ddg.extend(a_mean_ddg)
        all_rmean_ddg.extend(r_mean_ddg)
        k_ddg.extend(known_ddg)

        fig, axarr = conv.create_ax(4, 2, shx=True, shy=True)

        a_vals = draw_plot(axarr[0,0], known_ddg, a_mean_ddg, 'b', "Mean", "Known DDG", "Amber DDG")    
        r_vals = draw_plot(axarr[1,0], known_ddg, r_mean_ddg, 'b', "Mean", "Known DDG", "Rosetta DDG")
        append_values_to_corr_dict( corr_values_dict, a_vals, r_vals, "Mean")

        a_vals = draw_plot(axarr[0,1], known_ddg, a_min_ddg, 'b', "Min DG", "Known DDG", "Amber DDG")    
	r_vals = draw_plot(axarr[1,1], known_ddg, r_min_ddg, 'b', "Min DG", "Known DDG", "Rosetta DDG")
        append_values_to_corr_dict( corr_values_dict, a_vals, r_vals, "Min DG")

        a_vals = draw_plot(axarr[0,2], known_ddg, a_minr_ddg, 'b', "Min RMSD", "Known DDG", "Amber DDG")    
        r_vals = draw_plot(axarr[1,2], known_ddg, r_minr_ddg, 'b', "Min RMSD", "Known DDG", "Rosetta DDG")
        append_values_to_corr_dict( corr_values_dict, a_vals, r_vals, "Min RMSD")

        a_vals = draw_plot(axarr[0,3], known_ddg, a_bottom3_ddg, 'b', "Bottom 3", "Known DDG", "Amber DDG")    
        r_vals = draw_plot(axarr[1,3], known_ddg, r_bottom3_ddg, 'b', "Bottom 3", "Known DDG", "Rosetta DDG")
        append_values_to_corr_dict( corr_values_dict, a_vals, r_vals, "Bottom 3")

        conv.save_fig(fig, out_csv_path + "/" + os.path.splitext(os.path.basename(rec_corr))[0] + ".txt", "ddg", 16, 8) 
    labels_vals=["-PCC","-Rho","-Mae"]
    labels=["Mean", "Min DG", "Min RMSD", "Bottom 3"]
    for x_ind in range(0,3):
        for y_ind,label in enumerate(labels):
	    vals_list = corr_values_dict[label]
            scatterplot.draw_actual_plot(axarr_all[x_ind,y_ind], vals_list[x_ind], vals_list[x_ind+3],'b',label+labels_vals[x_ind],"Amber","Rosetta", size=40)
            axarr_all[x_ind,y_ind].relim()
            # update ax.viewLim using the new dataLim
            axarr_all[x_ind,y_ind].autoscale_view()
 	    if x_ind == 2:   
 	        axarr_all[x_ind,y_ind].set_xlim([-0.2,10.0])
                axarr_all[x_ind,y_ind].set_ylim([-0.2,10.0])
    	    	scatterplot.add_x_y_line(axarr_all[x_ind,y_ind],0.0,10.0)
	    else:
		axarr_all[x_ind,y_ind].set_xlim([-1.2,1.2])
                axarr_all[x_ind,y_ind].set_ylim([-1.2,1.2])
                scatterplot.add_x_y_line(axarr_all[x_ind,y_ind],-1.0,1.0)
	    
    conv.save_fig(fig_all, out_csv_path + "/all.txt", "ddg", 16, 12)

    draw_plot(axarr_all_corr[0,0], k_ddg, all_amean_ddg, 'b', "All Mean", "Known DDG", "Amber DDG")
    draw_plot(axarr_all_corr[1,0], k_ddg, all_rmean_ddg, 'b', "All Mean", "Known DDG", "Rosetta DDG")

    conv.save_fig(fig_all_corr, out_csv_path + "/all_corr.txt", "ddg", 4, 8)

def append_values_to_corr_dict( corr_values_dict, a_vals, r_vals, label ):
    corr_values_dict[label][0].append(a_vals[0])
    corr_values_dict[label][1].append(a_vals[1])
    corr_values_dict[label][2].append(a_vals[2])
    corr_values_dict[label][3].append(r_vals[0])
    corr_values_dict[label][4].append(r_vals[1])
    corr_values_dict[label][5].append(r_vals[2])


def mean_abs_error(known, pred):
    sum_err=0.0
    for k,p in zip(known, pred):
        sum_err += np.abs(k-p)
    return sum_err/len(known)   

def draw_plot(ax, x, y, color, x_axis, y_axis, title):
    scatterplot.draw_actual_plot(ax, x, y, color, x_axis, y_axis, title, size=40)
    coeff, pval = pearsonr(x, y)
    rho, pval = spearmanr(x, y)
    mae = mean_abs_error(x, y)
    conv.add_text_dict(ax, { "PCC" : coeff, "Rho" : rho, "MAE" : mae })
        
    scatterplot.add_x_y_line(ax, min_val=min(x), max_val=max(x))

    return [coeff, rho, mae]

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
