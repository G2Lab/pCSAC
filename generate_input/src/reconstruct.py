# pcSAC modules
import sys
sys.path.append("/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/src/")

import globals

# External modules
## Hi-C modules
import hicstraw
#import fanc

## Correlation metrics
from skimage.metrics import structural_similarity
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics.pairwise import cosine_similarity
from math import trunc
from scipy.spatial.distance import cdist


## General modules
import sys
import numpy as np
import scipy as scp
import pandas as pd
from pathlib import Path
import random
import os
import glob
import re

## Plotting modules
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.pyplot as plt 
import seaborn as sns



def ssim(x, y):
    print("SSIM")
    data_range = y.max()-y.min()
    if (x.flatten() == y.flatten()).all():
        return -100
    return structural_similarity(x, y, data_range=data_range)

def pearsons(x, y):
    print("PCC")
    r, p_value = pearsonr(x.flatten(), y.flatten())
    return r

def mse(x, y):
    print("MSE")
    return mean_squared_error(x, y)

def mae(x, y):
    print("MAE")
    return mean_absolute_error(x, y)

def psnr(x, y):
    print("PSNR")
    data_range = np.max(x) - np.min(x)
    err = mean_squared_error(x, y)
    if err == 0:
        return 100
    else:
        return 10*np.log10((data_range**2)/err)

def spearmans(x, y):
    print("SCC")
    r, p_value = spearmanr(x.flatten(), y.flatten())
    return r

def cosine(x, y):
    print("COSINE")
    x = x.flatten()
    y = y.flatten()
    return cosine_similarity(x.reshape(1,-1), y.reshape(1,-1)).item()

def mse_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        r = mean_squared_error(vector1, vector2)
        result[key] = r

    return result

def mae_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        r = mean_absolute_error(vector1, vector2)
        result[key] = r

    return result

def pearson_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        r,_ = pearsonr(vector1, vector2)
        result[key] = r

    return result
        

def spearman_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        r,_ = spearmanr(vector1, vector2)
        result[key] = r

    return result

def cosine_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        r = cosine_similarity(vector1, vector2)
        result[key] = r 
        
    return result

def GenomicDistance_correlation(reconstructed_matrix, whole_M,resolution,genomic_metric):
    df_rec = matrix2df(reconstructed_matrix)
    df_whole = matrix2df(whole_M)

    gdist_rec = genomicDistance(df_rec, resolution)
    gdist_org = genomicDistance(df_whole, resolution)


    cs_gd = genomic_metric(gdist_rec, gdist_org)
    min_val = resolution * 2
    max_val = resolution * 150
    cs_gd_sel = {k: v for k, v in cs_gd.items() if min_val <= k <= max_val}

    return cs_gd_sel

def jaccard_index(vector1, vector2):
    intersection = np.sum(vector1 & vector2)
    union = np.sum(vector1 | vector2)

    return intersection / union if union != 0 else 0

def jaccard_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = np.array(dict1[key]).astype(bool)
        vector2 = np.array(dict2[key]).astype(bool)
        r = jaccard_index(vector1, vector2)
        result[key] = r 
        
    return result

def GenomicDistance_jaccard(reconstructed_matrix, whole_M,resolution):
    print("Jaccard only works for binary matrices")
    df_rec = matrix2df(reconstructed_matrix)
    df_whole = matrix2df(whole_M)

    gdist_rec = genomicDistance(df_rec, resolution)
    gdist_org = genomicDistance(df_whole, resolution)


    cs_gd = jaccard_gd(gdist_rec, gdist_org)
    min_val = resolution
    max_val = resolution * 200
    cs_gd_sel = {k: v for k, v in cs_gd.items() if min_val <= k <= max_val}

    return cs_gd_sel


def plot_GDC(reconstructed, original, resolution, save_path,save=True):
    metric_fun = [mse_gd, mae_gd, pearson_gd, spearman_gd]
    metric_res = []
    metric_name = ("MSE", "MAE", "PCC", "SCC")

    for metric, name in zip(metric_fun, metric_name):
        gd = GenomicDistance_correlation(reconstructed, original, resolution, metric)
        metric_res.append(gd)
        lists = sorted(gd.items())
        x, y = zip(*lists)

        # Plotting
        plt.figure(figsize=(4.5, 3.5))
        plt.plot(x, y, color='#AD343E', linewidth=1.5, label='pcSAC')

        plt.xlabel('Genomic Distance')
        plt.ylabel(name) 
        plt.grid(True, which="both", ls="--", c='0.85')
        plt.ylim(0,None)
        plt.xlim(resolution*2, resolution * 100)
        plt.legend(loc='best')

        ax = plt.gca()
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

        if save:
            save_filename = os.path.join(save_path, f"{name}.png")
            plt.savefig(save_filename, bbox_inches='tight')

        plt.show()

# Modified function to plot in two columns and two rows

def plot_GDC_compare(pcSAC_matrix, matdhic_matrix, original_matrix, resolution):
    metric_fun = [mse_gd, mae_gd, pearson_gd, spearman_gd]
    metric_name = ("MSE", "MAE", "PCC", "SCC")

    # Creating a figure with 2 rows and 2 columns of subplots
    fig, axs = plt.subplots(2, 2, figsize=(9, 7))  # Adjust the figsize as needed
    axs = axs.flatten()  # Flattening the 2D array of axes for easy iteration

    for i, (metric, name) in enumerate(zip(metric_fun, metric_name)):
        gd = GenomicDistance_correlation(pcSAC_matrix, original_matrix, resolution, metric)
        lists = sorted(gd.items())
        x, y = zip(*lists)

        gd_another = GenomicDistance_correlation(matdhic_matrix, original_matrix, resolution, metric)
        lists_another = sorted(gd_another.items())
        x_another, y_another = zip(*lists_another)

        # Plotting on the respective subplot
        axs[i].plot(x, y, color='#AD343E', linewidth=1.5, label='pcSAC')
        axs[i].plot(x_another, y_another, color='#202C59', linewidth=1.5, label='DeepHiC')

        axs[i].set_xlabel('Genomic Distance')
        axs[i].set_ylabel(name)
        axs[i].grid(True, which="both", ls="--", c='0.85')
        axs[i].set_xlim(resolution * 2, resolution * 100)
        # axs[i].set_ylim(0, 0.3)
        axs[i].legend(loc='best')

        for spine in axs[i].spines.values():
            spine.set_visible(False)
        axs[i].spines['bottom'].set_visible(True)
        axs[i].spines['left'].set_visible(True)

    plt.tight_layout()
    plt.show()

# Modified function to plot just one Jaccard Index graph

def plot_jaccard(pcSAC_matrix, matdhic_matrix, original_matrix, resolution, start_bin=1 ,end_bin=100):
    metric_fun = jaccard_gd
    metric_name = "Jaccard Index"

    plt.figure(figsize=(5, 4))

    gd = GenomicDistance_jaccard(pcSAC_matrix, original_matrix, resolution)
    lists = sorted(gd.items())
    x, y = zip(*lists)

    gd_another = GenomicDistance_jaccard(matdhic_matrix, original_matrix, resolution)
    lists_another = sorted(gd_another.items())
    x_another, y_another = zip(*lists_another)

    # Plotting the Jaccard Index
    plt.plot(x, y, color='#DF8088', linewidth=1.5, linestyle='dashed', label='pcSAC')
    plt.plot(x_another, y_another, color='#546191', linewidth=1.5, linestyle='dashed', label='DeepHiC')

    plt.xticks(rotation=45, fontsize=9)
    plt.xlabel('Genomic Distance')
    plt.ylabel(metric_name)
    plt.grid(True, which="both", ls="--", c='0.85')
    plt.xlim(resolution * start_bin, resolution * end_bin)
    plt.legend(loc='best')

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    plt.tight_layout()
    plt.show()


def create_mask(matrix, n_off_diagonals=2):
    mask = np.ones_like(matrix, dtype=bool)
    for k in range(-n_off_diagonals, n_off_diagonals + 1):
        np.fill_diagonal(mask, 0, wrap=True)
        mask = np.logical_and(mask, np.eye(*matrix.shape, k=k, dtype=bool) == 0)
    return mask

def upper_mask(reconstructed_matrix):
    mask = create_mask(reconstructed_matrix)
    upper_reconstructed = np.triu(reconstructed_matrix) * mask
    return upper_reconstructed

def all_metrics(matrix1,matrix2):
    metrics = {
        "SSIM": ssim(matrix1, matrix2),
        "PCC": pearsons(matrix1, matrix2),
        "SCC": spearmans(matrix1, matrix2),
        "MSE": mse(matrix1, matrix2),
        "MAE": mae(matrix1, matrix2),
        "PSNR": psnr(matrix1, matrix2),
        "COSINE": cosine(matrix1, matrix2)
    
    }
    return metrics

def gaussian_smooth(matrix,sigma):
    from scipy.ndimage.filters import gaussian_filter
    return gaussian_filter(matrix, sigma=sigma)


def compare_metrics(metrics_upper_reconstructed, metrics_upper_reconstructed2):
    metrics = list(metrics_upper_reconstructed.keys())
    num_metrics = len(metrics)

    fig, axs = plt.subplots(1, num_metrics, figsize=(num_metrics * 3, 4))

    for i, metric in enumerate(metrics):
        values = [metrics_upper_reconstructed[metric], metrics_upper_reconstructed2[metric]]
        labels = ['pcSAC', 'DeepHiC']

        axs[i].plot(labels, values, marker='o', color='#3256a8', linewidth=1.5)
        axs[i].grid(True, which="both", ls="--", c='0.85')
        axs[i].set_title(metric)
        axs[i].set_xticks(labels)
        axs[i].tick_params(axis='x', rotation=45)

        if metric == 'MAE' or metric == 'MSE':
            axs[i].set_ylim(min(values)-0.001, max(values)+0.001) 
        else:  
            axs[i].set_ylim(min(values)-0.01, max(values) + 0.01)

        for spine in axs[i].spines.values():
            spine.set_visible(False)

    plt.tight_layout()
    plt.show()

''' Parameters:
    <array> @matrix: Hi-C matrix to be converted
    <int> @resolution: Resolution of the matrix
    <string> @chromosome: Chromosome name
    <int> @start: Start position of the matrix on the chromosome
    <int> @end: End position of the matrix on the chromosome'''

def matrix2coordinates(matrix, resolution, chromosome, start, end):
    
    # Get the upper triangle indices
    rows, cols = np.triu_indices(matrix.shape[0])
    #non_zero_indices = matrix[rows, cols] != 0
    #rows = rows[non_zero_indices]
    #cols = cols[non_zero_indices]
    values = matrix[rows, cols]

    rows_transformed = [(x * resolution) + start for x in rows]
    cols_transformed = [(x * resolution) + start for x in cols]

    df = pd.DataFrame({
        'bin1_id': rows_transformed,
        'bin2_id': cols_transformed,
        'count': values,
    })

    # Reorder columns and sort
    df = df[['bin1_id','bin2_id', 'count']]
    #df = df[df['bin1_id'] != df['bin2_id']]
    df = df.sort_values(by=['bin1_id', 'bin2_id'])
    
    return df


def matrix2coordinates_noVals(matrix, resolution, chromosome, start, end):
    
    # Get the upper triangle indices
    rows, cols = np.triu_indices(matrix.shape[0])
    #non_zero_indices = matrix[rows, cols] != 0
    #rows = rows[non_zero_indices]
    #cols = cols[non_zero_indices]
    values = matrix[rows, cols]

    rows_transformed = [(x * resolution) + start for x in rows]
    cols_transformed = [(x * resolution) + start for x in cols]

    df = pd.DataFrame({
        'bin1_id': rows_transformed,
        'bin2_id': cols_transformed,
    })

    # Reorder columns and sort
    df = df[['bin1_id','bin2_id']]

    #df = df[df['bin1_id'] != df['bin2_id']]
    df = df.sort_values(by=['bin1_id', 'bin2_id'])
    
    return df

def matrix2coords_for_cool_for_cmd(matrix, resolution, chromosome, start_pos, output_file):
    """
    Convert a Hi-C matrix to a list of interaction coordinates suitable for cooler.

    Args:
    - matrix: 2D numpy array of Hi-C interaction counts.
    - resolution: The size of each bin (base pairs).
    - chromosome: The name of the chromosome (e.g., 'chr1').
    - start_pos: The genomic start position of the first bin (base pairs).
    - output_file: Path to the output file to save the interaction pairs.

    Returns:
    - Saves the interaction pairs to a file in a format suitable for cooler.
    """
    # Get the upper triangle indices
    rows, cols = np.triu_indices(matrix.shape[0])
    values = matrix[rows, cols]


    start1 = rows * resolution + start_pos
    start2 = cols * resolution + start_pos

    df = pd.DataFrame({
        'chrom1': chromosome,
        'start1': start1,
        'end1': start1 + resolution,
        'chrom2': chromosome,
        'start2': start2,
        'end2': start2 + resolution,
        'count': values
    })


    df.sort_values(by=['start1', 'start2'], inplace=True)

    

    return df

def matrix_to_pixels(matrix, resolution):
    """
    Convert a Hi-C interaction matrix to a DataFrame of pixel values with bin indices.
    
    Args:
    - matrix: 2D numpy array of Hi-C interaction counts.
    - resolution: The size of each bin in base pairs.
    
    Returns:
    - A DataFrame with 'bin1_id', 'bin2_id', and 'count' columns.
    """
    num_bins = matrix.shape[0] -1
    print(num_bins)
    bin1_id, bin2_id = np.triu_indices(num_bins, k=0)
    counts = (matrix[bin1_id, bin2_id])*10000
    
    # Filter out zero counts to reduce file size and processing time
    nonzero = counts > -1
    bin1_id = bin1_id[nonzero]
    bin2_id = bin2_id[nonzero]
    counts = counts[nonzero]
    
    return pd.DataFrame({'bin1_id': bin1_id, 'bin2_id': bin2_id, 'count': counts})

import math
def cool_generate_bins(chromosome, start_pos, end_pos, resolution):
    """Generate genomic bins for a single chromosome, ensuring full coverage."""
    # Use math.ceil to ensure rounding up, thus including any partial bin at the end
    num_bins = math.ceil((end_pos - start_pos) / resolution) 
    bins = pd.DataFrame({
        'chrom': [chromosome] * num_bins,
        'start': [start_pos + i * resolution for i in range(num_bins)],
        'end': [min(start_pos + (i + 1) * resolution, end_pos) for i in range(num_bins)]
    })
    print(len(bins))
    return bins

def coordinates2matrix(df,resolution,start,end):
        max_index = int((end - start) / resolution)
        matrix = np.zeros((max_index + 1, max_index + 1))

        for index, row in df.iterrows():
               
                row_idx = int((row[0] -  start) / resolution)
                col_idx = int((row[1] -  start) / resolution)

                if 0 <= row_idx <= max_index and 0 <= col_idx <= max_index:
                    matrix[row_idx, col_idx] = row[2]
                    matrix[col_idx, row_idx] = row[2]

        return matrix

''' Parameters for matrix2cool:
    <array> @matrix: Hi-C matrix to be converted
    <int> @bsize: Resolution of the matrix
   <string> @chr: Chromosome name 'chrN' (N is the number of the chromosome)
    <int> @start: Start position of the matrix on the chromosome
    <int> @end: End position of the matrix on the chromosome
    <int> @n_chains: Number of chains used to reconstruct the matri
    <string> @outpath: Path to save the cooler file
    <string> @name: Name for the cooler file
    return: Cooler file
'''


def matrix2cool(matrix,outpath,name,bsize,chr,start,end,n_chains):
    from cooler.create import ArrayLoader
    import cooler as clr     
    
    start = int(start) # retrieve from globals temporal
    end = int(end)
    max_val = n_chains #counts are normalized by the number of chains
    valid_chr = str(chr) # temporal only valid 
    chromsizes = clr.util.fetch_chromsizes('hg19').loc[valid_chr:valid_chr]
    bins_tut = clr.binnify(chromsizes, bsize) #all bins in the chromosome
    bins = bins_tut[(bins_tut['start'] >= start) & (bins_tut['end'] <= end)] #bins in the region of interest

    n_trial = np.triu(matrix[:-1, :-1]*(max_val)) #original counts matrix upper triangle
    

    filename=f"{outpath}/{name}.cool"

    clr_1 = clr.create.ArrayLoader(bins, n_trial, chunksize=100000)
    clr.create_cooler(filename, bins, clr_1, assembly='hg19')

    cool_file = clr.Cooler(filename)
    
    return cool_file

    
''' Parameters for calculate_HiCRep:
    <cool> @cool_1: Original Hi-C matrix
    <cool> @cool_2: Reconstructed Hi-C matrix
    <int> @bsize: Resolution of the matrix

'''
def calculate_HiCRep(cool_1,cool_2,bsize):
    from hicrep import hicrepSCC
    # maximal genomic distance to include in the calculation
    binSize=bsize
    dBPMax=100000 #2Mb
    h=10#smoothing parameter
    bDownSample=True
    hiC_rep = hicrepSCC(cool_1, cool_2, h, dBPMax, bDownSample)

    return hiC_rep[0]


def chain_from_file(df): 
    #only takes the first 3 columns of the df because contain the coordinates of the monomers
    x = df[df.columns[0]] #takes first column of the df
    y = df[df.columns[1]] #similarly for the second column
    z = df[df.columns[2]] #similarly for the third column
    a = [np.array([x[i], y[i], z[i]]) for i in range(len(x))] #get all the coordinates of the monomers and store them in a list
    return a

def interactions_from_chains(chain, dist):
    df = pd.DataFrame(chain)
    distances = cdist(df, df)
    distances[distances < dist] = 1
    distances[distances > dist] = 0
    
    return distances

def matrix2df(matrix):
    df = pd.DataFrame([(i+1, j+1, matrix[i, j]) for i in range(matrix.shape[0]) for j in range(matrix.shape[1])])
    return df

def genomicDistance(df,resolution):
    df[0] = (df[0]) * resolution
    df[1] = (df[1]) * resolution

    df['distance'] = abs(df[1] - df[0])

    result_dict = {}
    for _, row in df.iterrows():
        diff = row['distance']
        if diff not in result_dict:
            result_dict[diff] = []
        result_dict[diff].append(row[2])
        
    return result_dict

def cosine_similarity_genomicD(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}

    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]

        dot_product = sum([v1*v2 for v1, v2 in zip(vector1, vector2)])
        magnitude1 = sum([v**2 for v in vector1]) ** 0.5
        magnitude2 = sum([v**2 for v in vector2]) ** 0.5
        
        if magnitude1 * magnitude2 == 0:
            similarity = 0 
        else:
            similarity = dot_product / (magnitude1 * magnitude2)
        
        result[key] = similarity
        
    return result

def run_reconstruction(main_path,downsampling_factor,lrf,run_option,number_chains,region_size,hic_res,chromosome,start,cell_t="GM12878"):
    import subprocess
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
        "-m", str(main_path),
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


def run_metrics(original,reconstructed,path,hic_r):
    import statistics
    import sys
    

    original=upper_mask(original)
    reconstructed=upper_mask(reconstructed)

    pattern = r'data_(\d+)_resolution'
    match = re.search(pattern, path)
    r = hic_r
    #Saving metrics
    original_stdout = sys.stdout
    cfile = 'correlation_metrics.txt'
    cpath = os.path.join(path +"/"+"correlation/")
    if not os.path.exists(cpath):
        os.makedirs(cpath)
    
    plot_GDC(reconstructed,original,r,cpath)

    output_file  = os.path.join(cpath, cfile)
    with open(output_file, 'w') as file:
        # Redirect the standard output to the file
        sys.stdout = file
        print("Correlation Metrics")
        print(ssim(np.triu(reconstructed) , np.triu(original)))
        print(pearsons(np.triu(reconstructed) , np.triu(original)))
        print(spearmans(np.triu(reconstructed) , np.triu(original)))
        print(mse(np.triu(reconstructed) , np.triu(original)))
        print(mae(np.triu(reconstructed) , np.triu(original)))
        print(psnr(np.triu(reconstructed) , np.triu(original)))
        print(cosine(np.triu(reconstructed) , np.triu(original)))
    sys.stdout = original_stdout

def annotation(path,lrf):
    pattern = r'data_(\d+)_resolution'
    match = re.search(pattern, path)
    resolution = match.group(1)

    celltype = path.split("plots/")[1].split("/")[0].split("-")[0]

    pattern = r'(\d+)_downSample'
    match = re.search(pattern, path)
    df = match.group(1)

    pattern = r'(\d+)_chains'
    match = re.search(pattern, path)
    chain_sample = match.group(1)
    lrf_str = str(lrf)

    table = {"Parameter":["Cell Type","Downsampling factor","HiC resolution", "Reconstructed chains","LRF Matrix (in kb)"],
            "Value":[celltype,df,resolution,chain_sample,lrf_str]}

    df_table = pd.DataFrame(table)

    return df_table

def main_plot(M_N,P_N,P_LR,O,D,path,LRF):
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(22, 4))
    # Plot the first matrix M_N
    im0 = axes[0].matshow(M_N, cmap="rocket_r", vmin=0, vmax=1)
    cbar0 = fig.colorbar(im0, ax=axes[0], shrink=0.8)
    cbar0.ax.set_ylabel(ylabel='', rotation=270, labelpad=11)

    # Plot the third matrix P_LR
    im2 = axes[1].matshow(P_LR, cmap="rocket_r", vmin=0, vmax=1)
    cbar2 = fig.colorbar(im2, ax=axes[1], shrink=0.8)
    cbar2.ax.set_ylabel(ylabel='', rotation=270, labelpad=11)

    # Plot the 4 matrix O
    im3 = axes[2].matshow(O, cmap="rocket_r", vmin=0, vmax=1)
    cbar3 = fig.colorbar(im3, ax=axes[2], shrink=0.8)
    cbar3.ax.set_ylabel(ylabel='', rotation=270, labelpad=11)

    # Plot the 5matrix D
    im4 = axes[3].matshow(D, cmap='RdGy_r', interpolation='none', vmax=1, vmin=-1)
    cbar4 = fig.colorbar(im4, ax=axes[3], shrink=0.8)
    cbar4.ax.set_ylabel(ylabel='', rotation=270, labelpad=11)
    avg_error = 100*np.average(D)
    axes[3].text(0.5, -0.1, f"Total error: {avg_error:.2f}%", ha="center", transform=axes[3].transAxes)

    axes[4].axis('off')
    lowfactor = LRF
    print(path)
    print(LRF)
    df_table = annotation(path,str(lowfactor))
    table = axes[4].table(cellText=df_table.values, colLabels=df_table.columns, cellLoc='center', loc='center')
    table.scale(1.4, 3)

    for cell in table.properties()["celld"].values():
        cell.set_fontsize(60) 

    plt.tight_layout()

    plt.savefig(os.path.join(path, "all_matrices.pdf"))

def read_results(chains,chain_v,dc,verbose=False):
    files = sorted([f for f in os.listdir(chains) if os.path.isfile(os.path.join(chains, f))])
    dist=dc
    n_chains=chain_v
    cumulative_matrix = None
    chains_rec = []
    if chain_v == chain_v:
        print("Chains were successfully generated. Calculating the reconstructed matrix...")
        for i in range(int(chain_v)):
            if verbose:
                print(f"Processing chain {i+1}/{n_chains}")
                
            chain = chain_from_file(pd.read_csv((chains + "/" + files[i]), sep="\t", comment="#",header=None))
            interactions = interactions_from_chains(chain,dist)
            if cumulative_matrix is None:
                cumulative_matrix = np.zeros_like(interactions)
            cumulative_matrix += interactions
        
        return cumulative_matrix/n_chains
    
    else:
        
        print("Chains were not successfully generated. \nAborting...")
        sys.exit(1)

def main_compare(chr,reconstruction_path,chain_v,directory,dc,norm):  
    if (norm == "clamp"):
        d_all = globals.all
    else:
    #Setting up directories to plots
        d_all = globals.all_nc

    tmp = "_"+ str(chain_v) + "_chains"
    #GM12878-HRC_data_5000_resolution_2000000_Mb_16downSample_40_kb_LR_5000_chains_1_rho2 = reconstruction_path
    plots_specific = reconstruction_path.split(tmp)[0]

    p_output_dir = os.path.join(d_all,f"chr{chr}", plots_specific, str(chain_v)+"_chains"+ (reconstruction_path.split(tmp)[1])+"_" + str(dc)+ "_dc")
    print(p_output_dir)
    
    os.makedirs(p_output_dir, exist_ok=True)

    default_sigma=0.6

    #Retrieving the original files to compare
    

    if (norm == "clamp"):
        whole_M = globals.input_matrices
    else:
        whole_M = globals.input_matrices_nc
        
    whole_M = os.path.join(whole_M +f"/chr{chr}/"+ plots_specific)
    interaction = glob.glob(whole_M +"/"+ 'interaction_*')

    M_N = glob.glob(whole_M +"/"+ 'whole*')
    M_N = np.loadtxt(M_N[0])

    P_N = glob.glob(whole_M + "/"+'down*')
    P_LR = glob.glob(whole_M +"/"+ 'lr*')
    P_N = np.loadtxt(P_N[0])
    P_LR = np.loadtxt(P_LR[0])
    O = read_results(directory,chain_v,dc)
    
    O = gaussian_smooth(O,0.6)
    o_name = "reconstructed_matrix_"+str(chain_v)+"chains_"+str(dc)+"_dc.txt"
    o_path = os.path.join(whole_M, o_name)
    o_arr = np.asarray(O)
    np.savetxt(o_path, o_arr, delimiter=" " )

    D=abs(O-M_N)

    return (M_N,P_N,P_LR,O,D,p_output_dir)

def merged_submatrices_plot(M_N,R_M,path,LRF,number_chains,df,hic_res):
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

    # Plot the first matrix M_N
    im0 = axes[0].matshow(M_N, cmap="rocket_r", vmin=0, vmax=1)
    cbar0 = fig.colorbar(im0, ax=axes[0], shrink=0.8)
    cbar0.ax.set_ylabel(ylabel='', rotation=270, labelpad=11)

    # Plot the reconstructed matrix R_M
    im2 = axes[1].matshow(R_M, cmap="rocket_r", vmin=0, vmax=1)
    cbar2 = fig.colorbar(im2, ax=axes[1], shrink=0.8)
    cbar2.ax.set_ylabel(ylabel='', rotation=270, labelpad=11)

    # Adding a table with the parameters on the third subplot
    parameters = {'LRF': LRF,
                  'number_chains': number_chains,
                  'downsampling_factor': df,
                  'hic_res': hic_res}
    
    table_data = list(parameters.items())
    ax_table = axes[2].axis('tight')
    ax_table = axes[2].axis('off')
    axes[2].table(cellText=table_data, colLabels=['Parameter', 'Value'], cellLoc='center', loc='center')
    
    plt.tight_layout()

    filename = f"reconstructed_matrices_chains_{number_chains}_downSample_{df}_resolution_{hic_res}_{LRF}_kb.png"
    print("Plotting was successful. Saving the plot as:", filename)
    plt.savefig(os.path.join(path, filename))


