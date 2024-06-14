#!/usr/bin/env python
import sys
import os
import hdbscan
import numpy as np
def hdbscan_clustering(file_name,mat_file,minpts,minsamps,cores):
    path = os.getcwd()
    os.chdir(path)
    d = np.genfromtxt(mat_file, delimiter=',')
    clusterer = hdbscan.HDBSCAN(min_cluster_size = minpts, min_samples = minsamps, core_dist_n_jobs = cores).fit(d)
    np.savetxt('clus.%s.%s.%s.txt' % (file_name, minpts, minsamps), clusterer.labels_.astype(int), fmt = '%i')

