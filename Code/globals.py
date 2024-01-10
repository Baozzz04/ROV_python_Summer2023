import numpy as np

def initialize():
    global bb
    bb = np.zeros((6, 6))
    global damping
    damping = np.zeros((6, 6))
    global added_mass
    added_mass = np.zeros((6, 6))

    global Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Iyx, Izx, Izy
    Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Iyx,Izx,Izy = 0,0,0,0,0,0,0,0,0 

    global m, Css
    m = 0
    Css = 0

    global V_rov, xg, yg, zg, xb, yb, zb, ro 
    V_rov,xg,yg,zg,xb,yb,zb,ro = 0,0,0,0,0,0,0,0

    global F_tb1, F_tb2, phi_tb1, si_tb1, teta_tb1, phi_tb2, si_tb2, teta_tb2, rx_tb1, ry_tb1, rz_tb1, rx_tb2, ry_tb2, rz_tb2
    F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2 = 0,0,0,0,0,0
    si_tb2,teta_tb2 = 0,0
    rx_tb1,ry_tb1,rz_tb1,rx_tb2,ry_tb2,rz_tb2 = 0,0,0,0,0,0

    global Apxy, Apxz, Apyz
    Apxy = 0
    Apxz = 0
    Apyz = 0

    global fileID, fileID2
    fileID = open('lucCanTongHop.txt', 'w')
    fileID2 = open('ketquara_diaphuong.txt', 'w')

    global Xprop, Yprop, Zprop, Kprop, Mprop, Nprop
    Xprop, Yprop, Zprop, Kprop, Mprop, Nprop = 0,0,0,0,0,0

    global fmt_2
    fmt_2 = np.zeros((3, 3))

    global r1, l1, r4, b2
    r1, l1, r4, b2 = 0,0,0,0

    global kq, kq2
    kq, kq2 = 0,0

    global r4_m, r1_m, b2_m, c2_m, l1_m
    r4_m, r1_m, b2_m, c2_m, l1_m = 0,0,0,0,0

    global X_HS, Y_HS, Z_HS, fK_HS, fM_HS, fN_HS
    X_HS, Y_HS, Z_HS, fK_HS, fM_HS, fN_HS = 0,0,0,0,0,0

    global Time_run, u_dauvao, deltaT
    Time_run, u_dauvao, deltaT = 0,0,0
    global cd_luc, vv_f
    vv_f = np.zeros(3)

    global fmt
    fmt = np.zeros((3, 3))

    global phi_teta_si
    phi_teta_si = np.zeros((3, 1))