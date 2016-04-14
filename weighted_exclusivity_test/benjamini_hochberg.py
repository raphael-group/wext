#!/usr/bin/env python

import numpy as np

def benjamini_hochberg_rejection(unchanged_p_values, alpha, independent=True):
    """
    Compute the Benjamini-Hochberg multiple-hypothesis rejections for p-values p and false discovery rate alpha.

    Example (from Benjamini and Hochberg (1995)):
    In [1]: p_values = [0.0201, 0.7590, 0.0344, 0.0278, 0.0095,
                        0.0298, 0.6719, 1.0000, 0.0001, 0.4262,
                        0.0004, 0.6528, 0.3240, 0.0459, 0.0019]
    In [2]: alpha = 0.05
    In [3]: benjamini_hochberg_rejection(p_values, alpha)
    Out[3]: [False False False False  True
             False False False  True False
             True False False False  True]
    """
    p_values = np.asarray(unchanged_p_values)
    m = np.size(p_values)
    c = 1.0 if independent else np.sum(1.0/np.arange(1, m+1, dtype=np.float64))

    indices = np.argsort(p_values)
    sorted_p_values = p_values[indices]

    k = 0
    for i, p in enumerate(sorted_p_values):
        if p<=alpha*c*(i+1)/m:
            k = i+1

    rejections = np.zeros(m, dtype=np.bool)
    rejections[:k] = True
    rejections[k:] = False

    return rejections[np.argsort(indices)]

def benjamini_hochberg_correction(unchanged_p_values, independent=True):
    """
    Compute the Benjamini-Hochberg multiple-hypothesis correction for p-values p.

    Example (from Benjamini and Hochberg (1995)):
    In [1]: p_values = [0.0201, 0.7590, 0.0344, 0.0278, 0.0095,
                        0.0298, 0.6719, 1.0000, 0.0001, 0.4262,
                        0.0004, 0.6528, 0.3240, 0.0459, 0.0019]
    In [2]: benjamini_hochberg_correction(p_values)
    Out[2]: [ 0.0603  0.8132  0.0645  0.0639  0.0356
              0.0639  0.7753  1.      0.0015  0.5812
              0.003   0.7753  0.486   0.0765  0.0095]
    """
    p_values = np.asarray(unchanged_p_values)
    m = np.size(p_values)
    c = 1.0 if independent else np.sum(1.0/np.arange(1, m+1, dtype=np.float64))

    indices = np.argsort(p_values)
    sorted_p_values = p_values[indices]

    sorted_q_values = sorted_p_values/(c*np.arange(1, m+1, dtype=np.float64)/m)
    adjusted_p_values = np.zeros(m)
    adjusted_p_values[m-1] = min(sorted_q_values[m-1], 1.0)
    for i in reversed(range(m-1)):
        adjusted_p_values[i] = min(sorted_q_values[i], adjusted_p_values[i+1])

    return adjusted_p_values[np.argsort(indices)]
