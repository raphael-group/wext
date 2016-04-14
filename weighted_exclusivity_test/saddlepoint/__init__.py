import numpy as np
from ..constants import *
from saddlepoint_exclusive_2_fast import saddlepoint_exclusive_2
from saddlepoint_exclusive_3_fast import saddlepoint_exclusive_3

def saddlepoint(t, x, p, tail=ONE_GREATER, verbose=False):
    if tail!=ONE_GREATER:
        raise NotImplementedError('Saddlepoint approximation not implemented for given tail')

    k = len(x)
    q = [np.asarray(a, dtype=np.float64) for a in p]

    if k==2:
        saddlepoint_p_value = saddlepoint_exclusive_2(t, x[0], x[1], q[0], q[1])
    elif k==3:
        saddlepoint_p_value = saddlepoint_exclusive_3(t, x[0], x[1], x[2], q[0], q[1], q[2])
    else:
        raise NotImplementedError('Saddlepoint approximation not implemented for k={}'.format(k))

    return saddlepoint_p_value

def unweighted_saddlepoint(t, x, p, tail=ONE_GREATER, verbose=False):
    k = len(p)
    n = len(p[0])
    unweighted_p = np.ones((k, n), dtype=np.float64)

    return saddlepoint(t, x, unweighted_p, tail, verbose)
