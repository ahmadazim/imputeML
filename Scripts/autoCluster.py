## Write script to automatically cluster cells

import scanpy as sc
import numpy as np
import pandas as pd
from math import log
from statistics import median

def cluster(inputMat, modelName):
    print('Working on {} cells and {} genes'.format(*inputMat.shape))
    dataPath = str('/home/ahmadazim/data/modelImputations' + modelName + '.txt')
    # Output data as txt file
    np.savetxt(dataPath, inputMat)
    # Import data (export and then import to keep record of data)
    data = sc.read_text(dataPath)

    print("Data imported.")

    data.var_names_make_unique()
    data.var['mt'] = data.var_names.str.startswith('MT-')  
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

    maxGene = input("Filter out all cells with n_genes_by_counts greater than: ")
    maxMT = input("Filter out all cells with pct_counts_mt greater than (input \"NA\" to ignore): ")

    data = data[data.obs.n_genes_by_counts < int(maxGene), :]
    if maxMT != "NA":
        data = data[data.obs.pct_counts_mt < int(maxMT), :]

    sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)

    data = data[:, data.var.highly_variable]
    sc.pp.regress_out(data, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(data, max_value=10)

    print("QC steps done.")

    sc.tl.pca(data, svd_solver='arpack')
    sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
    sc.tl.umap(data)
    sc.tl.leiden(data)

    print("Plotting PCA and UMAP...")    

    sc.pl.pca(data, color = 'leiden', save= str(modelName + '.png')) 
    sc.pl.umap(data, color = 'leiden', save= str(modelName + '.png'))


# def evaluateAccuracy(inputMat, pubAssign):
