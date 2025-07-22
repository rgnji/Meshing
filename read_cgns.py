import CGNS.MAP as cgnsmap

tree, links, paths = cgnsmap.load('practice.cgns')

cgns_base = tree[2][1]

zone_conn = []

for item in cgns_base[2]:
    if item[3] == 'Zone_t':
        zone_grid_conn = item[2][4]
        grid_conn_1to1 = zone_grid_conn[2][0]
        transform = grid_conn_1to1[2][0][1]
        point_range = grid_conn_1to1[2][1][1]
        point_range_donor = grid_conn_1to1[2][2][1]
        zone_conn.append({'transform': transform,
                          'point_range': point_range,
                          'point_range_donor': point_range_donor})
