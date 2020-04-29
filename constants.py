import numpy as np


SATS = [2, 4, 5, 9, 10, 12, 17, 23, 25]
sats = [2, 4, 5, 9, 10, 12, 17, 23, 25]

gps_ephm_columns = ['prn', 'M0', 'DeltaN', 'e', 'sqrtA', 'Omega0',
                    'Io', 'omega', 'OmegaDot', 'IDOT', 'Cuc',
                    'Cus', 'Crc', 'Crs', 'Cic', 'Cis', 'Toe',
                    'IODE', 'GPS_week', 'Toc', 'Af0', 'Af1', 'Af2',
                    '0']

receiver_col_start = ['Time', 'week_num', 'x', 'y',
                      'z', 'lat', 'lon', 'alt',
                      'num_sat', 'v_n', 'v_e', ' v_d']
it = 0
sat_cols = []

for sat in sats:
    temp = ['prn ' + str(sat),
            'snr ' + str(sat),
            'csc ' + str(sat),
            'pseudorange ' + str(sat),
            'cp ' + str(sat)]
    sat_cols.append(temp)
    receiver_col_start += temp

l_1_freq = 1575.42e6

SIGMA = 0.025

c = 2.99792458e8

WAVELENGTH = c / l_1_freq

BASELINE_TRUE = 0.36

index = [353, 354, 355, 347, 359, 361, 365, 371, 372]