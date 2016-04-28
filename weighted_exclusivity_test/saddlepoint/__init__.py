import numpy as np
from ..constants import *
from saddlepoint import saddlepoint as saddlepoint_approximation

def saddlepoint(t, x, p, tail=ONE_GREATER, verbose=False):
    if tail!=ONE_GREATER:
        raise NotImplementedError('Saddlepoint approximation not implemented for given tail')

    return saddlepoint_approximation(t, x, p)
