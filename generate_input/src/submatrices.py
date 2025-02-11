# -*- coding: utf-8 -*-
"""
This program generates submatrices from an original matrix. 
v0.2 : 2023/09/04
"""
import sys

print("Importing pcSAC libraries...")

# pcSAC module
from . import utils
from . import globals
from . import reconstruct

import numpy as np
import scipy as scp
import pandas as pd
import hicstraw
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.pyplot as plt 
import seaborn as sns
from math import trunc
import os
import glob
import subprocess
np.set_printoptions(threshold=sys.maxsize)

"""
Testing
"""

SAVE = True
hic_file = globals.hic_file
clamp_n = True
norm_q = False
diag_n = False


'''Variables for mapping_submat function
    matrix: original matrix whole
    submatrix: submatrix to be mapped
    o_rang: original range of the submatrix'''

def mapping_submat(matrix,submatrix,o_rang):
    sub_coord = reconstruct.mat2coord(submatrix,o_rang) # converting my submatrix to coordinates
    org_coord = reconstruct.mat2coord(matrix,0) #converting my main big matrix to coordinates
    sub_df = pd.DataFrame(sub_coord, columns=['X', 'Y', 'Value'],dtype=object) #making data frames for handeling data
    org_df = pd.DataFrame(org_coord, columns=['X', 'Y', 'Value'],dtype=object)

    #mapping submatrix values to original coordinates in main matrix
    sub_df.set_index(['X', 'Y'], inplace=True)
    org_df.set_index(['X', 'Y'], inplace=True)
    org_df.update(sub_df) # Update the values in org_df with the matching values from submatrix coordinates
    org_df.reset_index(inplace=True)

    #re-doing original matrix with updated values 
    updt_org = org_df.values.tolist()
    updt_org_mat = reconstruct.coord2mat(updt_org,0)

    return updt_org_mat
    

def map_submatrix_to_main(main_matrix, submatrix, start_coord):
    """
    Map the submatrix values onto the main matrix based on the given starting coordinate.

    :param main_matrix: The main matrix where submatrix values will be mapped.
    :param submatrix: The submatrix to be mapped onto the main matrix.
    :param start_coord: A tuple indicating the starting coordinate (row, col) of the submatrix relative to the main matrix.
    :return: Modified main matrix with the submatrix values.
    """
    
    # Extract the dimensions of the submatrix
    sub_rows, sub_cols = submatrix.shape

    # Extract the starting coordinates
    start_row, start_col = (start_coord,start_coord)

    # Map the submatrix onto the main matrix
    main_matrix[start_row:start_row+sub_rows, start_col:start_col+sub_cols] = submatrix

    return main_matrix

''' Variables for function merge_submatrices
    coords: vector of coordinates for the submatrices, start:end
    number_chains: number of chains to be reconstructed'''

def merge_submatrices(coords,number_chains,downsample_f,hic_res,lrf):
    str_coords=[coord.split(":")[0] for coord in coords]
    submat_files = globals.submat_main  
    reconstructed_subm = []
    for i in str_coords:
        base_path = os.path.join(submat_files, "input_matrices")
        search_path = os.path.join(base_path, f"*{hic_res}*_{i}_start*{downsample_f}*{lrf}*")
        submatrix_i = glob.glob(search_path)
        extract_parameters = os.path.basename(submatrix_i[0])
        parameters=extract_parameters.split("_")
        #Need to check this section because it will take all resolutions and not only the one we are working with
        cell=parameters[0]
        hic_res=parameters[2]
        df=parameters[8]
        LRF=parameters[10]   
        for dir_path in submatrix_i:
            reconstructed = glob.glob(os.path.join(dir_path, f'reconstructed*{number_chains}_chains*'))
            O = np.loadtxt(reconstructed[0])
            reconstructed_subm.append(O)
    print("Scanning for original whole matrix, if not generated, run the pipeline with --run_option whole")
    if clamp_n:
        original_mat = globals.restrictions_path
        HR_mat = globals.input_matrices
    else:
        original_mat = globals.restrictions_path_nc
        HR_mat = globals.input_matrices_nc

    filename_results = os.path.join(original_mat, 
                                    f"{df}_{LRF}_kb_LR")
    HR_pattern = f"{cell}*{hic_res}*{df}*{LRF}*"
    HR_mat_path = os.path.join(HR_mat,HR_pattern)

    pattern = f"{cell}*{hic_res}*{df}*{LRF}*{number_chains}_chains*"
    whole_reconstructed_mat = glob.glob(os.path.join(filename_results, pattern))

    if not whole_reconstructed_mat:
        raise ValueError("No whole reconstructed matrix found. Run the pipeline with --run_option whole and then submatrices.")
    else:
        #print(whole_reconstructed_mat[0])
        collision=int(hic_res) / 100 #bp to amstrongs
        interaction_distance=int((collision * 2) + 5)
        O = reconstruct.read_results(whole_reconstructed_mat[0],number_chains,interaction_distance)
        
        M_N = glob.glob(HR_mat_path +"/"+ 'whole*')
        M_N = np.loadtxt(M_N[0])
        
    updt_org_mat = O.copy()
    for i, (submatrix, str_coord) in enumerate(zip(reconstructed_subm, str_coords)):
        updt_org_mat = map_submatrix_to_main(updt_org_mat,submatrix,int(str_coord))
    
    print("Plotting subtmatrices merged into main matrix...")
    plot_path = os.path.dirname(HR_mat)
    plot_path = os.path.dirname(plot_path)
    p_output_dir = os.path.join(plot_path, "submatrices/plots")
    reconstruct.merged_submatrices_plot(M_N,updt_org_mat,p_output_dir,LRF,number_chains,df,hic_res)

    return updt_org_mat,M_N

''' Reconstruction section '''
submat_main = globals.submat_main

''' Variables for function run_reconstruction

    lrf: low resolution factor i.e if HR is 5kb and lrf is 4, LR is 20kb
    downsampling_factor: downsampling factor for the HiC matrix
    run_option: submatrices or whole
    number_chains: number of chains to be reconstructed
    region_size: size of the region to be reconstructed
    hic_res: resolution of the HiC matrix
    chromosome : chromosome where the region to be reconstructed is located
    start: start coordinate of the region to be reconstructed

    '''
#coords = [(0,180,0,180),(130,310,130,310),(220,400,220,400)]
def run_reconstruction(downsampling_factor,lrf,run_option,number_chains,region_size,hic_res,chromosome,start,cell_t="GM12878"):
    #Each submatrix is ran with that interaction_distance, therefore main matrix should be ran with the same interaction_distance
    smooth_interactions = True

    if hic_res == 5000:
        collision=int(((hic_res) / 100)) #bp to amstrongs
        interaction_distance=int((collision * 2) - 1)
    else:
        collision=int(((hic_res) / 100) - 10) #bp to amstrongs
        interaction_distance=int((collision * 2) - 10)
        if smooth_interactions:
            interaction_distance = int(interaction_distance - 10)
    cmd = [
        "sbatch", "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/03.compare_matrices.sh",
        "-d", str(downsampling_factor),
        "-z", str(region_size),
        "-h", str(hic_res),
        "-l", str(lrf),
        "-o", str(run_option),
        "-n", str(number_chains),
        "-c", str(chromosome),
        "-s", str(start),
        "-i", str(interaction_distance),
        '-t', str(cell_t)
    ]
    print(cmd)
    subprocess.run(cmd)


# At the moment is doing it separately for each submatrix but we need to gather all of them 
def submatrices_f(directory,reconstruction_path,number_chains,plots_main,input_matrices,dc):
    print("Reconstructing submatrices...")
    print("Reconstruction path:",reconstruction_path)
    tmp = "_"+ str(number_chains) + "_chains"
    plots_specific = reconstruction_path.split(tmp)[0]
    output_dir=os.path.join(plots_main, plots_specific, 
                 str(number_chains)+"_chains"+ 
                 (reconstruction_path.split(tmp)[1])+"_" 
                 + str(dc)+ "_dc")
    os.makedirs(output_dir, exist_ok=True)
    submatrix_i = os.path.join(input_matrices + plots_specific)
    print("Submatrix path:",submatrix_i)
    interaction =  glob.glob(submatrix_i +"/"+ 'interaction_*')

    M_N = glob.glob(submatrix_i +"/"+ 'whole*')
    M_N = np.loadtxt(M_N[0])

    P_N = glob.glob(submatrix_i + "/"+'down*')
    P_N = np.loadtxt(P_N[0])

    P_LR = glob.glob(submatrix_i+"/"+ 'lr*')
    P_LR = np.loadtxt(P_LR[0])


    O = reconstruct.read_results(directory,number_chains,dc)
    reconstruct_name = "reconstructed_matrix_"+str(number_chains)+"_chains"+".txt" #whole matrix is normalized
    reconstruct_path = os.path.join(submatrix_i,reconstruct_name)
    n = np.asarray(O)
    np.savetxt(reconstruct_path, n, delimiter=" ")

    D= abs(O-M_N)

    return (M_N,P_N,P_LR,O,D,output_dir)
    
    #reconstruct.

def reconstruct_matrix(norm_option, run_option,number_chains, directory, reconstruction_path,dc):
    if (norm_option == "clamp" and run_option == "submatrices"):
        submat_main = globals.submat_main
        plots_main = submat_main + "plots/"
        input_matrices = submat_main + "input_matrices/"
    elif (norm_option != "clamp" and run_option == "submatrices"):
        submat_main = globals.submat_main_nc
        plots_main = submat_main + "plots/"
        input_matrices = submat_main + "input_matrices/"

    if(run_option == "submatrices"):
        main_matrix = globals.reconstructed_matrix
        sub_reconstruction_matrix = globals.sub_reconstruction
        tmp = "_"+ str(number_chains) + "_chains"
        M_N,P_N,P_LR,O,D,output_dir = submatrices_f(directory,reconstruction_path,number_chains,plots_main,input_matrices,dc)

    return(M_N,P_N,P_LR,O,D,output_dir)




''' Variables for function ensambles_3D
    lrf: low resolution factor i.e if HR is 5kb and lrf is 4, LR is 20kb
    downsampling_factor: downsampling factor for the HiC matrix
    run_option: submatrix or whole
    number_chains: number of chains to be generated
    region_size: size of the region to be reconstructed
    hic_res: resolution of the HiC matrix
    '''

def ensambles_3D(downsampling_factor,lrf,run_option,number_chains,region_size,hic_res,split_chains,cell):
    cmd = [
        "sbatch", "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/02.restrictions.sh",
        "-d", str(downsampling_factor),
        "-r", str(region_size),
        "-h", str(hic_res),
        "-l", str(lrf),
        "-o", str(run_option),
        "-n", str(number_chains),
        "-s", str(split_chains),
        "-c", str(cell)
    ]
    subprocess.run(cmd)

''' Creates HiC plots for each of the submatrices'''
def plot_sm(dense_matrix, maxcolor, suffix,ax,start,end):
    LR_F = int((hic_res * l_res) / 1000)
    im =  plt.matshow(dense_matrix, cmap="rocket_r", vmin=0, vmax=maxcolor)
    cbar = plt.colorbar(im, shrink=0.8)
    cbar.ax.set_ylabel(ylabel='', rotation=270, labelpad=11)
    filename = os.path.splitext(os.path.basename(hic_file))[0]
    current_directory = os.getcwd()
    pwd = os.path.dirname(current_directory)

    #output_dir = pwd + "/" + "results" + "/" + "plots" + "/" + filename + "_data_" + str(hic_res) +"_resolution_"+ str(size) + "_Mb_" + str(d_factor) + "downSample_" + str(LR_F) + "_kb_LR"
    if(clamp_n != True):
        output_dir = os.path.join(pwd, 
                    "results",
                    "submatrices", 
                    "plots", 
                    "{}_data_{}_resolution_{}_start_{}_end_{}_downSample_{}_kb_LR".format(
                        filename, 
                        hic_res, 
                        start, end, 
                        d_factor, 
                        LR_F ))
        
    if(clamp_n == True):
        norm_type = "clamp"
        output_dir = os.path.join(pwd, 
                    "results",norm_type, 
                    "submatrices",
                    "plots", 
                    "{}_data_{}_resolution_{}_start_{}_end_{}_downSample_{}_kb_LR".format(
                        filename, 
                        hic_res, 
                        start, end,
                        d_factor, 
                        LR_F ))
        
    os.makedirs(output_dir, exist_ok=True)
    print(output_dir)
    pln = suffix + '_HiCMatrix.png'
    plp = os.path.join(output_dir, pln)
    plt.savefig(plp)


''' Must be implemented in the generate data section, add options for whole or submatrices'''

''' Creates interaction files for each of the submatrices'''

def interaction_sm(interact_matrix,start_row,end_row,start_col,end_col,suffix):
    LR_F = int((hic_res * l_res) / 1000)
    rm = (int(hic_res)/100)/2 #radius in amstrongs
    filename = os.path.splitext(os.path.basename(hic_file))[0]
    current_directory = os.getcwd()
    pwd = os.path.dirname(current_directory)

    if(clamp_n != True):
        output_dir = os.path.join(pwd, 
                    "results",
                    "submatrices", 
                    "input_matrices", 
                    "{}_data_{}_resolution_{}_start_{}_end_{}_downSample_{}_kb_LR".format(
                        filename, 
                        hic_res, 
                        start_row, end_row, 
                        d_factor, 
                        LR_F ))
        
    if(clamp_n == True):
        norm_type = "clamp"
        output_dir = os.path.join(pwd, 
                    "results",norm_type, 
                    "submatrices",
                     "input_matrices",  
                    "{}_data_{}_resolution_{}_start_{}_end_{}_downSample_{}_kb_LR".format(
                        filename, 
                        hic_res, 
                        start_row, end_row,
                        d_factor, 
                        LR_F ))
        
    os.makedirs(output_dir, exist_ok=True)
    
    if(suffix == 'M_N'):
        p2n = "int_mat_seg_len_"+filename+".txt"
        p2 = os.path.join(output_dir,p2n)
        Path(p2).touch(exist_ok=True)

        p3n = "whole_matrix_"+filename+".txt" #whole matrix is normalized
        p3 = os.path.join(output_dir,p3n)
        n = np.asarray(interact_matrix)
        np.savetxt(p3, n, delimiter=" ")
        
        m = (interact_matrix.shape[0]-1)*[str(2*rm)+"\n"]
        with open(p2, "w+") as g:
            g.writelines(m)
        

    if(suffix == 'P_N'):
        p4n = "downsample_matrix_"+filename+".txt" #P_N
        p4 = os.path.join(output_dir,p4n)
        p_arr = np.asarray(interact_matrix)
        np.savetxt(p4, p_arr, delimiter=" ")


    if(suffix == 'P_LR'):
        pn = "interaction_matrix_"+filename+".txt"
        p = os.path.join(output_dir,pn)
        Path(p).touch(exist_ok=True)
        #l can be modified to LR if we want to reduce the number of restrictions
        l = [str(i*l_res)+" "+str((i+1)*l_res-1)+" "+str(j*l_res)+" "+str((j+1)*l_res-1)+" "+str((interact_matrix)[i][j])+"\n" for j in range(1,(interact_matrix).shape[0]) for i in range(j)]
        
        with open(p, "w+") as f:
            f.writelines(l)

        p5n = "lr_downsample_matrix_"+filename+".txt" #P_LR
        p5 = os.path.join(output_dir,p5n)
        lr = np.asarray(interact_matrix)
        np.savetxt(p5, lr, delimiter=" ") #Saving LR matrixes
  
''' Binning main matrix into required submatrices'''

def bin_matrix(matrix, suffix, coordinates):
    ncols = len(coordinates)  # The number of submatrices is the number of coordinates
    nrows = 1  # Assuming you want one column of subplots
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(18, 18), sharex=True, sharey=True)
    axes = axes.ravel()  # Flatten axes to make it easier to iterate over
    
    for idx, (start_row, end_row, start_col, end_col) in enumerate(coordinates):
        ax = axes[idx]
        if suffix == "P_LR":
            submatrix = matrix[start_row//l_res:end_row//l_res, start_col//l_res:end_col//l_res]
        else :
            submatrix = matrix[start_row:end_row, start_col:end_col]
        
        if suffix == "P_LR":
            print(f"Creating matrix {suffix} plot from row {start_row} to {end_row} and from column {start_col} to {end_col}")
            plot_sm(submatrix, 1, suffix, ax, figsize=(4, 4))
            if SAVE:
                print(f"Creating matrix {suffix} file from row {start_row} to {end_row} and from column {start_col} to {end_col}")
                utils.interaction_file(submatrix, start_row, end_row, start_col, end_col, suffix)
        else:
            print(f"Creating matrix {suffix} plot from row {start_row} to {end_row} and from column {start_col} to {end_col}")
            plot_sm(submatrix, 1, suffix, ax, figsize=(4, 4))
            if SAVE:
                print(f"Creating matrix {suffix} file from row {start_row} to {end_row} and from column {start_col} to {end_col}")
                utils.interaction_file(submatrix, start_row, end_row, start_col, end_col, suffix)

    plt.tight_layout()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--func', type=str, help='Function name to run')
    parser.add_argument('--number_chains', type=int, help='Number of chains for merge_submatrices')
    parser.add_argument('--downsample_f', type=int, help='Downsample factor of main and submatrices (should match)')
    parser.add_argument('--hic_res', type=int, help='Resolution of the HiC matrix and submatrices (should match)')
    parser.add_argument('--lrf', type=int, help='Low resolution factor of the HiC matrix and submatrices (should match)')

    args = parser.parse_args()

    if args.func == 'merge_submatrices':
        print("Individual submatrices are being merged with main matrix...")
        R_M,M_N = merge_submatrices(globals.coords,args.number_chains,args.downsample_f,args.hic_res,args.lrf)
