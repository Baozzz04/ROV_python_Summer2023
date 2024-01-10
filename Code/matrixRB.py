import numpy as np
import globals

def MatrixRB(RB, Cdc, r1, l1, r2, l2, r3, l3, r4, l4, r5, l5, b, b1, b2, H5, x1, x2, Ca, Ca1):
    for i in range(6):
        for j in range(6):
            RB[i, j] = 0

    RB[0, 0] = globals.m
    RB[0, 4] = globals.m * globals.zg
    RB[0, 5] = -globals.m * globals.yg

    RB[1, 1] = globals.m
    RB[1, 3] = -globals.m * globals.zg
    RB[1, 5] = globals.m * globals.xg

    RB[2, 2] = globals.m
    RB[2, 3] = globals.m * globals.yg
    RB[2, 4] = -globals.m * globals.xg

    RB[3, 1] = -globals.m * globals.zg
    RB[3, 2] = globals.m * globals.yg
    RB[3, 3] = globals.Ixx
    RB[3, 4] = globals.Ixy
    RB[3, 5] = globals.Ixz

    RB[4, 0] = globals.m * globals.zg
    RB[4, 2] = -globals.m * globals.xg
    globals.Iyx = globals.Ixy
    globals.Izy = globals.Iyz
    globals.Izx = globals.Ixz
    RB[4, 3] = globals.Iyx
    RB[4, 4] = globals.Iyy
    RB[4, 5] = globals.Iyz

    RB[5, 0] = -globals.m * globals.yg
    RB[5, 1] = globals.m * globals.xg
    RB[5, 3] = globals.Izx
    RB[5, 4] = globals.Izy
    RB[5, 5] = globals.Izz

    return (
        RB, Cdc, r1, l1, r2, l2, r3, l3, r4, l4, r5, l5, b, b1, b2, H5, x1, x2, Ca, Ca1
    )