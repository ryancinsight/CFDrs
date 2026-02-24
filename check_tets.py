import itertools

def check_jacobian(v_idx):
    v = [
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    ]
    p0, p1, p2, p3 = v[v_idx[0]], v[v_idx[1]], v[v_idx[2]], v[v_idx[3]]
    
    # (p1 - p0) . ((p2 - p0) x (p3 - p0))
    v1 = (p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2])
    v2 = (p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2])
    v3 = (p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2])
    
    # cross product v2 x v3
    cp = (
        v2[1]*v3[2] - v2[2]*v3[1],
        v2[2]*v3[0] - v2[0]*v3[2],
        v2[0]*v3[1] - v2[1]*v3[0]
    )
    
    # dot product v1 . cp
    dot = v1[0]*cp[0] + v1[1]*cp[1] + v1[2]*cp[2]
    return dot

tets_a = [
    [0, 1, 3, 4],
    [1, 2, 3, 6],
    [4, 5, 6, 1],
    [4, 7, 6, 3],
    [1, 3, 4, 6],
]

tets_b = [
    [1, 0, 2, 5],
    [3, 0, 2, 7],
    [4, 0, 5, 7],
    [6, 2, 5, 7],
    [0, 2, 5, 7],
]

for name, tets in [("tets_a", tets_a), ("tets_b", tets_b)]:
    print(f"--- {name} ---")
    for tet in tets:
        val = check_jacobian(tet)
        if val < 0:
            # swap last two to flip sign
            corr = [tet[0], tet[1], tet[3], tet[2]]
            print(f"Fixing {tet} -> {corr} (val={check_jacobian(corr)})")
        else:
            print(f"OK {tet}")

