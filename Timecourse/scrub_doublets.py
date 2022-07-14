#!/usr/bin/env python
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import random

print(sys.path)

def predict_doublets(matrix):
    counts_matrix = scipy.io.mmread(matrix).T.tocsc()
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.10)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    return scrub.doublet_scores_obs_
    
