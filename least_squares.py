import math as m
import numpy as np
import pandas as pd
from scipy import linalg as la
import navpy
import constants as cn

import sys

sys.path.append('../..utils')
sys.path.append('./LAMBDA_Python_1.0')


def rounding(N):
    N_round = N.round(0)
    print()


def geometry_free(rho, cp):
    g_free = (rho * (1 / cn.WAVELENGTH) - cp)
    n_g_free = g_free.mean(2).round(0)
    return n_g_free


def brute_force():
    true_baseline = [36, 0.36]
    print()


def ecef_to_ned(x_ecef, R):
    x_temp = R.T @ x_ecef
    return x_temp


def get_ref_vector(x_ref):
    vector = x_ref / np.linalg.norm(x_ref)
    return vector


def calc_least_squares_float(data, H, cp, lat, long, h):
    float_solution = pd.DataFrame(index=data.index,
                                  columns=['x', 'y', 'z', 'lat', 'long', 'alt',
                                           'p_inv'])
    for i in range(len(H)):
        p_inv = la.pinv(H[i].T @ H[i])
        x_hat = p_inv @ H[i].T @ cp[i]
        ned = [float(x_hat[0]), float(x_hat[1]), float(x_hat[2])]
        float_solution.iloc[i, 0] = float(x_hat[0])
        float_solution.iloc[i, 1] = float(x_hat[1])
        float_solution.iloc[i, 2] = float(x_hat[2])
        coord = navpy.ned2lla(ned, lat * 180 / m.pi, long * 180 / m.pi, h,
                              latlon_unit='deg', alt_unit='m',
                              model='wgs84')
        float_solution.iloc[i, 3] = coord[0]
        float_solution.iloc[i, 4] = coord[1]
        float_solution.iloc[i, 5] = coord[2]
        float_solution.iloc[i, 6] = p_inv.flatten()
    return float_solution


def p_range_multi(base_df, sat_pos, rover_df, R):
    H_cand = []
    cp_list = []
    HH = []
    for Time in base_df.index:
        sats_temp = cn.sats.copy()
        cp = np.array([[]])
        H = np.array([[]])
        flag = True
        ref = base_df.Ref.loc[Time]

        x_ecef = sat_pos[ref].loc[Time, ['x', 'y', 'z']]
        x_ecef = x_ecef.to_numpy().reshape((-1, 1))

        del sats_temp[ref]

        x_ref = ecef_to_ned(x_ecef, R)
        vec = get_ref_vector(x_ref)

        for i in sats_temp:
            if sat_pos[cn.sats.index(i)].at[Time, 'x'] != 0:
                x_ecef_temp = sat_pos[cn.sats.index(i)].loc[Time,
                                                            ['x', 'y', 'z']]
                x_ecef_temp = x_ecef_temp.to_numpy().reshape((-1, 1))
                x_cur_temp = ecef_to_ned(x_ecef_temp, R)
                vec_temp = get_ref_vector(x_cur_temp)
                H_temp = (vec - vec_temp).T
                H_temp *= 1 / cn.WAVELENGTH
                cp_rover_1 = rover_df.loc[
                    Time, 'cp ' + str(i)]
                cp_base_1 = base_df.loc[
                    Time, 'cp ' + str(i)]
                cp_rover_2 = rover_df.loc[
                    Time, 'cp ' + str(cn.sats[ref])]
                cp_base_2 = base_df.loc[
                    Time, 'cp ' + str(cn.sats[ref])]
                cp_temp = (cp_rover_1 - cp_base_1) - (
                        cp_rover_2 - cp_base_2)

                if flag:
                    H = H_temp
                    cp = np.array([[cp_temp]])
                    flag = False
                else:
                    H = np.vstack((H, H_temp))
                    cp = np.vstack((cp, cp_temp))
        H_cand.append(H)
        ident = -1 * np.eye(len(sats_temp))
        HH.append(np.hstack((H, ident)))
        cp_list.append(cp)
    return [H_cand, HH, cp_list]


def least_squares_pos_solution(base_df, sat_pos, rover_df, R, lat, long, h):
    [H_cand, HH, cp] = p_range_multi(base_df, sat_pos, rover_df, R)
    temp = calc_least_squares_float(base_df, HH, cp, lat, long, h)
    N_float = temp['p_inv']
    return temp
    print()
