import globals

#TSDV
def TSDV():
    with open('thongsodauvao1.txt', 'r') as fid:
        global Ixx, Iyy, Izz, Ixy, Ixz
        global Iyz, Iyx, Izx, Izy
        global m
        global Time_run, u_dauvao
        global deltaT
        global V_rov, xg, yg, zg, xb, yb, zb, ro
        global Css
        global Apxy, Apxz, Apyz
        global r1, b2, r4, l1

        lines = fid.readlines()

        globals.ro, Cdc, globals.m, globals.r1, globals.l1, r2, l2, r3, l3, globals.r4, l4, r5, l5, \
        globals.Ixx, globals.Iyy, globals.Izz, globals.Ixy, globals.Ixz, globals.Iyz, globals.Apxy, globals.Apxz, globals.Apyz, \
        b, b1, globals.b2, H5, x1, x2, globals.xg, globals.yg, globals.zg, Ca, Ca1, \
        globals.deltaT, globals.V_rov, globals.Time_run, globals.u_dauvao, globals.F_tb1, globals.F_tb2, \
        globals.phi_tb1, globals.si_tb1, globals.teta_tb1, globals.phi_tb2, globals.si_tb2, globals.teta_tb2, \
        globals.rx_tb1, globals.ry_tb1, globals.rz_tb1, globals.rx_tb2, globals.ry_tb2, globals.rz_tb2, globals.Css = map(float, lines)

    globals.Iyx = globals.Ixy
    globals.Izx = globals.Ixz
    globals.Izy = globals.Iyz

    return (
        Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1, 
        globals.F_tb1, globals.F_tb2, globals.phi_tb1, globals.si_tb1, globals.teta_tb1, globals.phi_tb2, globals.si_tb2, 
        globals.teta_tb2, globals.rx_tb1, globals.ry_tb1, globals.rz_tb1, globals.rx_tb2, globals.ry_tb2, globals.rz_tb2
    )