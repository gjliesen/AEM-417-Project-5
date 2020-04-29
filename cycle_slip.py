import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import linalg as la
import navpy
import constants as cn


def carrier_phase_range_lambda_phi(r_i, n_i, k, K_i, del_t, del_T_i,
                                   d_eph_i, d_iono_i, d_tropo_i, pseudorange):
    lambda_phi = r_i - n_i * cn.wavelength + (k - K_i) + cn.c \
                 * (del_t - del_T_i) + d_eph_i + d_iono_i + d_tropo_i \
                 + pseudorange
