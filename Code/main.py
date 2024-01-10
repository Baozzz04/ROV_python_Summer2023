import numpy as np
import globals
import TSDV
import Lucchanvit_0
import matran_damping
import addedmass
import matrixRB
import matranB
import chaycham_rov

if __name__ == "__main__": 
    Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1, Cdc = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2
    Xuu,Yvv,Zww,Kpp,Mqq,Nrr,Kvv,Muu,m11_3D,m22,m33,m44,m55,m66,m11 = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    add = np.zeros((6, 6))
    DP = np.zeros((6, 6))
    added = np.zeros((6, 6))
    RB = np.zeros((6, 6))
    matran_B = np.zeros((6, 6))
    cc = np.zeros((6, 6))

    i = 0
    j = 0
    n = 0

    globals.initialize()

    #Goi thong so dau vao
    (Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1, globals.F_tb1, globals.F_tb2, globals.phi_tb1, globals.si_tb1, globals.teta_tb1, globals.phi_tb2, globals.si_tb2, globals.teta_tb2, globals.rx_tb1, globals.ry_tb1, globals.rz_tb1, globals.rx_tb2, globals.ry_tb2, globals.rz_tb2) = TSDV.TSDV()
    #Goi luc chan vit
    (globals.F_tb1,globals.F_tb2,globals.phi_tb1,globals.si_tb1,globals.teta_tb1,globals.phi_tb2,globals.si_tb2,globals.teta_tb2,globals.rx_tb1,globals.ry_tb1,globals.rz_tb1,globals.rx_tb2,globals.ry_tb2,globals.rz_tb2) = Lucchanvit_0.Lucchanvit_0(globals.F_tb1,globals.F_tb2,globals.phi_tb1,globals.si_tb1,globals.teta_tb1,globals.phi_tb2,globals.si_tb2,globals.teta_tb2,globals.rx_tb1,globals.ry_tb1,globals.rz_tb1,globals.rx_tb2,globals.ry_tb2,globals.rz_tb2)
    #Goi ma tran suy giam
    (DP, Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1) = matran_damping.matran_damping(DP, Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1)
    #matran damping
    globals.damping = DP
    #Goi ma tran khoi luong nuoc kem
    (m11_3D, m22, m33, m44, m55, m66, Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1, add) = addedmass.addedmass(m11_3D, m22, m33, m44, m55, m66, Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1, add)
    globals.added_mass = add

    print('Add')
    print(add)
    print('dp')
    print(DP)

    # Goi ma tran RB
    (RB, Cdc, globals.r1, globals.l1, r2, l2, r3, l3, globals.r4, l4, r5, l5, b, b1, globals.b2, H5, x1, x2, Ca, Ca1) = matrixRB.MatrixRB(RB,Cdc, globals.r1, globals.l1, r2, l2, r3, l3, globals.r4, l4, r5, l5, b, b1, globals.b2, H5, x1, x2, Ca, Ca1)
    print('MatrixRB')
    print(RB)

    # Goi ma tran B
    (add, RB, matran_B) = matranB.matranB(add, RB, matran_B)
    print('MatrixB')
    print(matran_B)

    for i in range(6):
        for j in range(6):
            globals.bb[i, j] = matran_B[i, j]

    globals.bb = np.linalg.inv(globals.bb)
    print('MatrixRB nghichdao')
    print(globals.bb)

    # Goi ham chaychamrov
    chaycham_rov.chaycham_rov()

globals.fileID.close()
globals.fileID2.close()
