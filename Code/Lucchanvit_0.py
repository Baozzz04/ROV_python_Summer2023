import numpy as np
import globals

#Lucchanvit_0
def Lucchanvit_0(F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2,si_tb2,teta_tb2,rx_tb1,ry_tb1,rz_tb1,rx_tb2,ry_tb2,rz_tb2):
    fmt = np.zeros((3, 3))
    c1 = np.zeros((3, 1))
    c2 = np.zeros((3, 1))
    FM_tb = np.zeros((6, 1))

    (F_tb1,phi_tb1,si_tb1,teta_tb1,rx_tb1,ry_tb1,rz_tb1,FM_tb) = luctuocbin(F_tb1,phi_tb1,si_tb1,teta_tb1,rx_tb1,ry_tb1,rz_tb1,FM_tb)
    globals.Xprop = FM_tb[0]
    globals.Yprop = FM_tb[1]
    globals.Zprop = FM_tb[2]
    globals.Kprop = FM_tb[3]
    globals.Mprop = FM_tb[4]
    globals.Nprop = FM_tb[5]

    (F_tb2,phi_tb2,si_tb2,teta_tb2,rx_tb2, ry_tb2,rz_tb2,FM_tb) = luctuocbin(F_tb2,phi_tb2,si_tb2,teta_tb2,rx_tb2,ry_tb2,rz_tb2,FM_tb)
    globals.Nprop = -globals.Nprop
    globals.Xprop = globals.Xprop + FM_tb[0]
    globals.Yprop = globals.Yprop + FM_tb[1]
    globals.Zprop = globals.Zprop + FM_tb[2]
    globals.Kprop = globals.Kprop + FM_tb[3]
    globals.Mprop = globals.Zprop + FM_tb[4]
    globals.Nprop = globals.Nprop + FM_tb[5]

    #viet vao file
    with open('Luc chan vit.txt', 'w') as fid:
        fid.write(str(globals.Xprop))
        fid.write(str(globals.Yprop))
        fid.write(str(globals.Zprop))
        fid.write(str(globals.Kprop))
        fid.write(str(globals.Mprop))
        fid.write(str(globals.Nprop))

    return (
        F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2,si_tb2,teta_tb2,rx_tb1,ry_tb1,rz_tb1,rx_tb2,ry_tb2,rz_tb2
    )
#tinh luc tuoc bin
def luctuocbin(F_tb, phi_tb, si_tb, teta_tb, x, y, z, F):
    c1 = np.zeros((3, 1))
    c1[0] = F_tb
    c1[1] = 0
    c1[2] = 0

    fmt = np.zeros((3, 3))
    (phi_tb, si_tb, teta_tb, fmt) = matranquanhe1(phi_tb, si_tb, teta_tb, fmt)
    c2 = np.dot(fmt, c1)

    F[0] = c2[0]
    F[1] = c2[1]
    F[2] = c2[2]
    F[3] = y * F[2] - z * F[1]
    F[4] = -x * F[2] + z * F[0]
    F[5] = x * F[1] - y * F[0]

    return F_tb, phi_tb, si_tb, teta_tb, x, y, z, F
#tinh ma tran quan he
def matranquanhe1(phi,si,teta,fmt):
    globals.fmt_2 = np.zeros((3, 3))
    globals.fmt_2 = fmt[0:3, 0:3]

    fmt[0, 0] = np.cos(si) * np.cos(teta)
    fmt[0, 1] = -np.sin(si) * np.cos(phi) + np.cos(si) * np.sin(teta) * np.sin(phi)
    fmt[0, 2] = np.sin(si) * np.sin(phi) + np.cos(si) * np.sin(teta) * np.cos(phi)
    fmt[1, 0] = np.sin(si) * np.cos(teta)
    fmt[1, 1] = np.cos(si) * np.cos(phi) + np.sin(si) * np.sin(teta) * np.sin(phi)
    fmt[1, 2] = -np.cos(si) * np.sin(phi) + np.sin(si) * np.sin(teta) * np.cos(phi)
    fmt[2, 0] = -np.sin(teta)
    fmt[2, 1] = np.cos(teta) * np.sin(phi)
    fmt[2, 2] = np.cos(teta) * np.cos(phi)
    globals.fmt_2[0:3, 0:3] = fmt[0:3, 0:3]

    return phi, si, teta, fmt