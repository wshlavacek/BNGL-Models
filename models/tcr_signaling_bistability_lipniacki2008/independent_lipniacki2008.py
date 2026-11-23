"""The independent 37-state right-hand side of Lipniacki et al. (2008), Appendix A.

The block between the two marker lines is copied byte for byte from the cell of
`verify_lipniacki2008.ipynb` that defines it, where it is described: it is written from the
paper's Appendix A equations and does not read BioNetGen's generated network.
`verify_lipniacki2008_steady_state.ipynb` imports it from here and asserts that the block is
still identical to that cell, so the two notebooks cannot drift apart.

States y1..y37 (index 0..36): y1 free agonist pMHC; y2..y10 agonist-TCR complexes (no LCK;
LCK with (Y, S59) = (u, u), (u, p); then (p, u), (p, p) at ITAM u, p, pp); y11..y20 the same
for antagonist; y21 free TCR; y22 TCR.pSHP1; y23 agonist-TCR.pSHP1; y24 agonist-TCR.LCK.pSHP1;
y25, y26 the antagonist pair; y27 free LCK; y28 SHP1; y29 pSHP1; y30, y31 ZAP70 u, p; y32..y34
MEK u, p, pp; y35..y37 ERK u, p, pp.
"""

import numpy as np

# --- begin: copied verbatim from verify_lipniacki2008.ipynb ---
p = {
    "b1": 0.3 / 30000,
    "b2": 0.3 / 30000,
    "d1": 0.05,
    "d2": 1.0,
    "lb": 0.3 / 100000,
    "ly1": 5.0 / 300000,
    "ly2": 0.3,
    "ls1": 0.1,
    "ls2": 0.5 / 300000,
    "tp": 0.05,
    "s0": 1e-5,
    "s1": 30.0 / 300000,
    "s2": 0.0006,
    "s3": 0.05,
    "z0": 2e-6,
    "z1": 5.0 / 100000,
    "z2": 0.02,
    "m1": 5.0 / 100000,
    "m2": 0.02,
    "e1": 5.0 / 300000,
    "e2": 0.02,
}


def rhs(_time, y):
    (
        y1, y2, y3, y4, y5, y6, y7, y8, y9, y10,
        y11, y12, y13, y14, y15, y16, y17, y18, y19, y20,
        y21, y22, y23, y24, y25, y26, y27, y28, y29, y30,
        y31, y32, y33, y34, y35, y36, y37,
    ) = y
    q = p
    dy = np.zeros(37)

    dy[0] = q["d1"] * sum((y2, y3, y4, y5, y6, y7, y8, y9, y10, y23, y24))
    dy[0] -= q["b1"] * y1 * y21
    dy[1] = q["b1"] * y1 * y21 + (q["s2"] + q["s3"]) * y23
    dy[1] -= (q["d1"] + q["lb"] * y27 + q["ly1"] * y29) * y2
    dy[2] = q["lb"] * y27 * y2 + q["ls1"] * y4
    dy[2] += (q["s2"] + q["s3"]) * y24
    dy[2] -= (q["d1"] + q["ly2"] + q["ls2"] * y37 + q["ly1"] * y29) * y3
    dy[3] = q["ls2"] * y37 * y3 - (q["d1"] + q["ly2"] + q["ls1"]) * y4
    dy[4] = q["ly2"] * y3 + q["ls1"] * y6
    dy[4] -= (q["d1"] + q["tp"] + q["ls2"] * y37 + q["ly1"] * y29) * y5
    dy[5] = q["ly2"] * y4 + q["ls2"] * y37 * y5
    dy[5] -= (q["d1"] + q["tp"] + q["ls1"]) * y6
    dy[6] = q["tp"] * y5 + q["ls1"] * y8
    dy[6] -= (q["d1"] + q["tp"] + q["ls2"] * y37 + q["ly1"] * y29) * y7
    dy[7] = q["tp"] * y6 + q["ls2"] * y37 * y7
    dy[7] -= (q["d1"] + q["tp"] + q["ls1"]) * y8
    dy[8] = q["tp"] * y7 + q["ls1"] * y10
    dy[8] -= (q["d1"] + q["ls2"] * y37 + q["ly1"] * y29) * y9
    dy[9] = q["tp"] * y8 + q["ls2"] * y37 * y9
    dy[9] -= (q["d1"] + q["ls1"]) * y10

    dy[10] = q["d2"] * sum(
        (y12, y13, y14, y15, y16, y17, y18, y19, y20, y25, y26)
    )
    dy[10] -= q["b2"] * y11 * y21
    dy[11] = q["b2"] * y11 * y21 + (q["s2"] + q["s3"]) * y25
    dy[11] -= (q["d2"] + q["lb"] * y27 + q["ly1"] * y29) * y12
    dy[12] = q["lb"] * y27 * y12 + q["ls1"] * y14
    dy[12] += (q["s2"] + q["s3"]) * y26
    dy[12] -= (q["d2"] + q["ly2"] + q["ls2"] * y37 + q["ly1"] * y29) * y13
    dy[13] = q["ls2"] * y37 * y13 - (q["d2"] + q["ly2"] + q["ls1"]) * y14
    dy[14] = q["ly2"] * y13 + q["ls1"] * y16
    dy[14] -= (q["d2"] + q["tp"] + q["ls2"] * y37 + q["ly1"] * y29) * y15
    dy[15] = q["ly2"] * y14 + q["ls2"] * y37 * y15
    dy[15] -= (q["d2"] + q["tp"] + q["ls1"]) * y16
    dy[16] = q["tp"] * y15 + q["ls1"] * y18
    dy[16] -= (q["d2"] + q["tp"] + q["ls2"] * y37 + q["ly1"] * y29) * y17
    dy[17] = q["tp"] * y16 + q["ls2"] * y37 * y17
    dy[17] -= (q["d2"] + q["tp"] + q["ls1"]) * y18
    dy[18] = q["tp"] * y17 + q["ls1"] * y20
    dy[18] -= (q["d2"] + q["ls2"] * y37 + q["ly1"] * y29) * y19
    dy[19] = q["tp"] * y18 + q["ls2"] * y37 * y19
    dy[19] -= (q["d2"] + q["ls1"]) * y20

    dy[20] = q["d1"] * sum((y2, y3, y4, y5, y6, y7, y8, y9, y10))
    dy[20] += q["d2"] * sum((y12, y13, y14, y15, y16, y17, y18, y19, y20))
    dy[20] += (q["s2"] + q["s3"]) * y22
    dy[20] -= q["b1"] * y1 * y21 + q["b2"] * y11 * y21 + q["ly1"] * y29 * y21
    dy[21] = q["ly1"] * y29 * y21
    dy[21] += q["d1"] * (y23 + y24) + q["d2"] * (y25 + y26)
    dy[21] -= (q["s2"] + q["s3"]) * y22
    dy[22] = q["ly1"] * y29 * y2 - (q["s2"] + q["s3"] + q["d1"]) * y23
    dy[23] = q["ly1"] * y29 * (y3 + y5 + y7 + y9)
    dy[23] -= (q["s2"] + q["s3"] + q["d1"]) * y24
    dy[24] = q["ly1"] * y29 * y12 - (q["s2"] + q["s3"] + q["d2"]) * y25
    dy[25] = q["ly1"] * y29 * (y13 + y15 + y17 + y19)
    dy[25] -= (q["s2"] + q["s3"] + q["d2"]) * y26
    dy[26] = q["d1"] * sum((y3, y4, y5, y6, y7, y8, y9, y10, y24))
    dy[26] += q["d2"] * sum((y13, y14, y15, y16, y17, y18, y19, y20, y26))
    dy[26] -= q["lb"] * (y2 + y12) * y27
    active_lck_1 = y5 + y7 + y9
    active_lck_2 = y15 + y17 + y19
    bound_shp = y22 + y23 + y24 + y25 + y26
    dy[27] = q["s2"] * (y29 + bound_shp)
    dy[27] -= q["s1"] * (active_lck_1 + active_lck_2) * y28 + q["s0"] * y28
    recruitable_tcr = y2 + y3 + y5 + y7 + y9 + y12 + y13 + y15 + y17 + y19 + y21
    dy[28] = q["s1"] * (active_lck_1 + active_lck_2) * y28
    dy[28] += q["s3"] * bound_shp + q["s0"] * y28
    dy[28] -= q["s2"] * y29 + q["ly1"] * recruitable_tcr * y29

    active_tcr = y9 + y10 + y19 + y20
    dy[29] = q["z2"] * y31 - q["z1"] * active_tcr * y30 - q["z0"] * y30
    dy[30] = q["z1"] * active_tcr * y30 + q["z0"] * y30 - q["z2"] * y31
    dy[31] = q["m2"] * (y33 + y34) - q["m1"] * y31 * y32
    dy[32] = q["m1"] * y31 * (y32 - y33) - q["m2"] * y33
    dy[33] = q["m1"] * y31 * y33 - q["m2"] * y34
    dy[34] = q["e2"] * (y36 + y37) - q["e1"] * y34 * y35
    dy[35] = q["e1"] * y34 * (y35 - y36) - q["e2"] * y36
    dy[36] = q["e1"] * y34 * y36 - q["e2"] * y37
    return dy
# --- end: copied verbatim from verify_lipniacki2008.ipynb ---
