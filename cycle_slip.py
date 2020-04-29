import constants as cn


def get(rover_df, base_df, ephem_df_list):
    drops = []
    for idx, sat in enumerate(cn.SATS):
        column = 'csc ' + str(sat)
        if rover_df[column].iat[0] != rover_df[column].iat[-1]:
            rover_df = rover_df.drop(cn.sat_cols[idx], axis=1)
            base_df = base_df.drop(cn.sat_cols[idx], axis=1)
            cn.sats.remove(sat)
            drops.append(idx)
    del ephem_df_list[drops[0]]
    del ephem_df_list[drops[1]-1]
    return [rover_df, base_df, ephem_df_list]
