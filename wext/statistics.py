#!/usr/bin/env python

import numpy as np

def multiple_hypothesis_correction(p_values_, method='BH'):
    """
    Compute a multiple-hypothesis correction using one of multiple
    methods: Bonferroni (bonferroni), Benjamini-Hochberg (BH), or
    Benjamini-Hochberg-Yekutieli (BY).

    Example (from Benjamini and Hochberg (1995)):
    In [1]: p_values = [0.0201, 0.7590, 0.0344, 0.0278, 0.0095,
                        0.0298, 0.6719, 1.0000, 0.0001, 0.4262,
                        0.0004, 0.6528, 0.3240, 0.0459, 0.0019]
    In [2]: q_values = multiple_hypothesis_correction(p_values, method='BH')
    In [3]: q_values
    Out[3]: array([ 0.0603    ,  0.81321429,  0.0645    ,  0.06385714,  0.035625  ,
                    0.06385714,  0.77526923,  1.        ,  0.0015    ,  0.58118182,
                    0.003     ,  0.77526923,  0.486     ,  0.0765    ,  0.0095    ])
    """
    if method not in ['bonferroni', 'BH', 'BY']:
        raise NotImplementedError('{} method not implemented'.format(method))

    valid_indices = [i for i, p_value in enumerate(p_values_) if 0.0<=p_value<=1.0]
    invalid_indices = [i for i, p_value in enumerate(p_values_) if not 0.0<=p_value<=1.0]

    p_values = np.asarray(p_values_)[valid_indices]
    n = len(p_values)

    if method=='bonferroni':
        q_values = np.minimum(n*p_values, 1)

    elif method=='BH':
        sorted_indices = np.argsort(p_values)
        sorted_p_values = p_values[sorted_indices]

        sorted_q_values = np.zeros(n)
        sorted_q_values[n-1] = min(sorted_p_values[n-1], 1.0)
        for i in reversed(range(n-1)):
            sorted_q_values[i] = min(float(n)/float(i+1)*sorted_p_values[i], sorted_q_values[i+1])

        q_values = np.zeros(n)
        q_values[sorted_indices] = sorted_q_values

    elif method=='BY':
        sorted_indices = np.argsort(p_values)
        sorted_p_values = p_values[sorted_indices]

        c = np.sum(1.0/np.arange(1, n+1, dtype=np.float64))
        sorted_q_values = np.zeros(n)
        sorted_q_values[n-1] = min(c*sorted_p_values[n-1], 1.0)
        for i in reversed(range(n-1)):
            sorted_q_values[i] = min(c*(float(n)/float(i+1))*sorted_p_values[i], sorted_q_values[i+1])

        q_values = np.zeros(n)
        q_values[sorted_indices] = sorted_q_values

    q_values_ = np.zeros(len(p_values_))
    q_values_[valid_indices] = q_values
    q_values_[invalid_indices] = float('nan')

    return q_values_
