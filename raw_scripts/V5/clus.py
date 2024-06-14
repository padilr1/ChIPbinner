#!/usr/bin/env python
import sys
import hdbscan
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import seaborn as sns
# import argparse
# import logging
import os
def hdbscan_clustering(file_name,mat_file,minpts,minsamps,baseline_samp_label,treated_samp_label,cores):
    path = os.getcwd()
    os.chdir(path)
    d = np.genfromtxt(mat_file, delimiter=',')
    clusterer = hdbscan.HDBSCAN(min_cluster_size = minpts, min_samples = minsamps, core_dist_n_jobs = cores).fit(d)
    # color_palette = sns.color_palette('deep', 8)
    # cluster_colors = [color_palette[x] if x >= 0
    #                 else (0.5, 0.5, 0.5)
    #                 for x in clusterer.labels_]
    # cluster_member_colors = [sns.desaturate(x, p) for x, p in
    #                         zip(cluster_colors, clusterer.probabilities_)]
    np.savetxt('clus.%s.%s.%s.txt' % (file_name, minpts, minsamps), clusterer.labels_.astype(int), fmt = '%i')
    # plt.scatter(*d.T, s=50, linewidth=0, c=cluster_member_colors, alpha=0.25)
    # plt.xlabel(baseline_samp_label)
    # plt.ylabel(treated_samp_label)
    # plt.savefig('clus.%s.%s.%s.png' % (file_name, minpts, minsamps), transparent=True)
