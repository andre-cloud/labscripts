import numpy as np
from scipy.constants import R


def boltzmann_perc(en, T):
    return np.exp(-(en-np.min(en))/R*T)/np.sum(np.exp(-(en-np.min(en))/R*T))

