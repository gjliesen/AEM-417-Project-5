import math as m
import numpy as np
import pandas as pd
from scipy import linalg as la
import navpy
import constants as cn

import sys

sys.path.append('../..utils')
sys.path.append('./LAMBDA_Python_1.0')


def get_static_solution(HHH, YYY):
    sol = []
    for i in range(len(HHH)):
        p_inv = la.pinv(HHH[i].T @ HHH[i]) @ HHH[i].T @ YYY[i]
        sol.append(p_inv)
    return sol


def get_static_x(static_sol):
    static_x = []
    for i in range(len(static_sol)):
        static_x.append(static_sol[i][0:3])
    return static_x


def get_N_float(static_sol):
    N_float = []
    for i in range(len(static_sol)):
        N_float.append(static_sol[i][3::])
    return N_float


def get_N_round(N):
    N_round = []
    for i in range(len(N)):
        N_round.append(N[i].round(0))
    return N_round


def get_geometry_free(rho, YYY):
    N_gfree = []
    for i in range(len(rho)):
        g_free = (rho[i] * (1 / cn.WAVELENGTH) - YYY[i])
        n_g_free = g_free.mean().round(0)
        N_gfree.append(n_g_free)
    return n_g_free


def create_R_matrix(YYY, num_sats):
    R = np.zeros((len(YYY), num_sats - 1))
    print()


def brute_force():
    true_baseline = [36, 0.36]
    N_cand
    print()


def ecef_to_ned(x_ecef, R):
    x_temp = R.T @ x_ecef
    return x_temp


def get_ref_vector(x_ref):
    vector = x_ref / np.linalg.norm(x_ref)
    return vector


def static_x_to_NED(data, static_x, lat, long, h):
    static_NED = pd.DataFrame(index=data.index,
                              columns=['x', 'y', 'z',
                                       'lat', 'long', 'alt'])
    for i in range(len(static_x)):
        ned = [float(static_x[i][0]), float(static_x[i][1]),
               float(static_x[i][2])]
        static_NED.iloc[i, 0] = float(static_x[i][0])
        static_NED.iloc[i, 1] = float(static_x[i][1])
        static_NED.iloc[i, 2] = float(static_x[i][2])
        coord = navpy.ned2lla(ned, lat * 180 / m.pi, long * 180 / m.pi, h,
                              latlon_unit='deg', alt_unit='m',
                              model='wgs84')
        static_NED.iloc[i, 3] = coord[0]
        static_NED.iloc[i, 4] = coord[1]
        static_NED.iloc[i, 5] = coord[2]
    return static_NED


def p_range_multi(base_df, sat_pos, rover_df, R):
    HHH_cand = []
    HHH = []
    YYY = []
    for Time in base_df.index:
        YY = []
        HH = []
        flag = True
        ref = base_df.Ref.loc[Time]

        # ecef_conversions
        x_ref_ecef = sat_pos[ref].loc[Time, ['x', 'y', 'z']]
        x_ref_ecef = x_ref_ecef.to_numpy().reshape((-1, 1))

        x_ref_ned = ecef_to_ned(x_ref_ecef, R)
        x_ref = get_ref_vector(x_ref_ned)

        for sat in cn.sats:
            if sat != cn.sats[ref]:
                # ecef_conversions
                x_ecef_sat = sat_pos[cn.sats.index(sat)].loc[Time,
                                                             ['x', 'y', 'z']]
                x_ecef_sat = x_ecef_sat.to_numpy().reshape((-1, 1))
                x_sat_ned = ecef_to_ned(x_ecef_sat, R)
                x_sat = get_ref_vector(x_sat_ned)

                # Measurement Matrix
                H = (1 / cn.WAVELENGTH) * (x_ref - x_sat).T

                # Psuedorange Difference Measurement Vector
                cp_rover_1 = rover_df.loc[Time, 'cp ' + str(sat)]
                cp_base_1 = base_df.loc[Time, 'cp ' + str(sat)]
                cp_rover_2 = rover_df.loc[Time, 'cp ' + str(cn.sats[ref])]
                cp_base_2 = base_df.loc[Time, 'cp ' + str(cn.sats[ref])]
                Y = (cp_rover_1 - cp_base_1) - (cp_rover_2 - cp_base_2)

                # Stack Measurements
                if flag:
                    HH = H
                    YY = np.array([Y])
                    flag = False
                else:
                    HH = np.vstack((HH, H))
                    YY = np.vstack((YY, Y))
        HHH_cand.append(HH)
        ident = -1 * np.eye(6)
        temp = np.hstack((HH, ident))
        HHH.append(temp)
        YYY.append(YY)

    return [HHH_cand, HHH, YYY]


def least_squares_pos_solution(base_df, sat_pos, rover_df, R, lat, long, h, rho):
    [HHH_cand, HHH, YYY] = p_range_multi(base_df, sat_pos, rover_df, R)
    static_sol = get_static_solution(HHH, YYY)
    static_x = get_static_x(static_sol)
    static_NED = static_x_to_NED(base_df, static_x, lat, long, h)
    N_float = get_N_float(static_sol)
    N_Round = get_N_round(N_float)
    N_gfree = get_geometry_free(rho, YYY)
    R = create_R_matrix(YYY, 7)
    return(static_NED)
