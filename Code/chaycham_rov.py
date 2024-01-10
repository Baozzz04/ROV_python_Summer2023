import numpy as np
import globals
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#tinh chay cham
def chaycham_rov():
    nstep = 5000

    k = 0
    u = globals.u_dauvao
    v, w, p, q, r, rr, x_i, y_i, z_i = 0,0,0,0,0,0,0,0,0
    p0, q0, r0, h, u_i, v_i, w_i = 0,0,0,0,0,0,0
    c2 = np.zeros((3, 1))
    c4 = np.zeros((3, 1))
    vantoc_phi, vantoc_teta, vantoc_si = 0,0,0
    record_x = np.zeros(nstep + 1)
    record_y = np.zeros(nstep + 1)
    biendem = 0
    phi = 0
    si = 0
    tetay = 0

    while (k <= nstep):
        (phi, tetay, si) = luchoiphuc_diaphuong(phi, tetay, si)
        xxx = 0
        yyy = 0
        zzz = 0
        p0 = 0
        q0 = 0
        r0 = 0
        # Giai pt trong toa do dia phuong
        (u, v, w, p, q, rr, xxx, yyy, zzz, p0, q0, r0) = giaipt_uwqteta1_6(u, v, w, p, q, rr, xxx, yyy, zzz, p0, q0, r0)
        # Doi sang toan cuc
        (phi, tetay, si, u, v, w, u_i, v_i, w_i, xxx, yyy, zzz, p, q, rr, vantoc_phi, vantoc_teta, vantoc_si, p0, q0, r0, c2, c4) = tinhvantoctoancuc(phi, tetay, si, u, v, w, u_i, v_i, w_i, xxx, yyy, zzz, p, q, rr, vantoc_phi, vantoc_teta, vantoc_si, p0, q0, r0, c2, c4)
        # lay gia tri cua cd_luc va vv_f, cho vao mang record_x, record_y
        record_x[biendem] = globals.cd_luc
        record_y[biendem] = globals.vv_f[0]
        biendem += 1
        
        x_i += xxx
        z_i += zzz  # theo do sau
        y_i += yyy  # theo chieu ngang
        
        h += zzz

        with open('ketquara.txt', 'w') as fid1:
            fid1.write(str(x_i) + '\n')
            fid1.write(str(y_i) + '\n')
            fid1.write(str(z_i) + '\n')
            fid1.write(str(u_i) + '\n')
            fid1.write(str(v_i) + '\n')
            fid1.write(str(w_i) + '\n')
            fid1.write(str(phi) + '\n')
            fid1.write(str(tetay) + '\n')
            fid1.write(str(si) + '\n')
        
        with open('ketquara_Z.txt', 'w') as fid2:
            fid2.write(str(z_i) + '\n')

        print(k, x_i, y_i, z_i, u_i, v_i, w_i, phi, tetay, si)
        
        (u, v, w, p, q, rr, c2, c4, phi, si, tetay) = vantocdiaphuongmoi(u, v, w, p, q, rr, c2, c4, phi, si, tetay)
        
        k += 1

    # ve do thi su dung 2 mang da luu
    plt.plot(record_y, record_x)
    plt.show()
    # Dong file
    fid1.close()
    fid2.close()
#Luc thuy tinh or luc phuc hoi
def luchoiphuc_diaphuong(phi, teta, si):
    v = globals.V_rov
    g = 9.81
    
    globals.X_HS = -(globals.m * g - v * globals.ro * g) * np.sin(teta)
    globals.Y_HS = (globals.m * g - v * globals.ro * g) * np.cos(teta) * np.sin(phi)
    globals.Z_HS = (globals.m * g - v * globals.ro * g) * np.cos(teta) * np.cos(phi)
    
    globals.fK_HS = (globals.yg * globals.m * g - globals.yb * v * globals.ro * g) * np.cos(teta) * np.cos(phi) - \
        (globals.zg * globals.m * g - globals.zb * v * globals.ro * g) * np.cos(teta) * np.sin(phi)
    globals.fM_HS = -(globals.zg * globals.m * g - globals.zb * v * globals.ro * g) * np.sin(teta) - (globals.xg * \
        globals.m * g - globals.xb * v * globals.ro * g) * np.cos(teta) * np.cos(phi)
    globals.fN_HS = (globals.xg * globals.m * g - globals.xb * v * globals.ro * g) * np.cos(teta) * np.sin(phi) - \
        (globals.yg * globals.m * g - globals.yb * v * globals.ro * g) * np.sin(teta)

    return phi, teta, si
#giai vi phan
def giaipt_uwqteta1_6(u, v, w, p, q, rr, x, y, z, p0, q0, r0):
    n = 12
    yy1 = np.zeros(n)
    y_out = np.zeros(n)
    YPRIME = np.zeros(n)

    # gan gia tri x..rr vao ma tran yy1
    yy1[0] = x
    yy1[1] = y
    yy1[2] = z
    yy1[3] = p0
    yy1[4] = q0
    yy1[5] = r0
    yy1[6] = u
    yy1[7] = v
    yy1[8] = w
    yy1[9] = p
    yy1[10] = q
    yy1[11] = rr

    tend = globals.deltaT
    T = 0
    lenw = n * n + 11 * n + 300
    #ham dao ham
    def uwqteta_cham_6_Luc(Yy, t):
        g = 9.81
        c = np.zeros(6)
        cc = np.zeros(6)
        CRB = np.zeros((6, 6))
        CA = np.zeros((6, 6))
        DV = np.zeros((6, 6))
        CCD = np.zeros((6, 6))
        VV = np.zeros(6)
        VVV = np.zeros(6)
        fluc1 = np.zeros(3)
        dydt = np.zeros(12)

        u = Yy[6]
        v = Yy[7]
        w = Yy[8]
        p = Yy[9]
        q = Yy[10]
        rr = Yy[11]
        
        globals.m11 = globals.added_mass[0, 0]
        globals.m22 = globals.added_mass[1, 1]
        globals.m33 = globals.added_mass[2, 2]
        globals.m44 = globals.added_mass[3, 3]
        globals.m55 = globals.added_mass[4, 4]
        globals.m66 = globals.added_mass[5, 5]

        Xuu = globals.damping[0, 0]
        Yvv = globals.damping[1, 1]
        Zww = globals.damping[2, 2]
        Kpp = globals.damping[3, 3]
        Mqq = globals.damping[4, 4]
        Nrr = globals.damping[5, 5]
        Muu = globals.damping[4, 0]
        Kvv = globals.damping[3, 1]

        CRB[0, 1] = -globals.m * rr
        CRB[0, 2] = globals.m * q
        CRB[1, 0] = globals.m * rr
        CRB[1, 2] = -globals.m * p
        CRB[2, 0] = -globals.m * q
        CRB[2, 1] = globals.m * p
        CRB[3, 4] = globals.Izz * rr
        CRB[3, 5] = -globals.Iyy * q
        CRB[4, 3] = -globals.Izz * rr
        CRB[4, 5] = globals.Ixx * p
        CRB[5, 3] = globals.Iyy * q
        CRB[5, 4] = -globals.Ixx * p

        CA[0, 4] = globals.m33 * w
        CA[0, 5] = -globals.m22 * v
        CA[1, 3] = -globals.m33 * w
        CA[1, 5] = globals.m11 * u
        CA[2, 3] = globals.m22 * v
        CA[2, 4] = -globals.m11 * u
        CA[3, 1] = globals.m33 * w
        CA[3, 2] = -globals.m22 * v
        CA[3, 4] = globals.m66 * rr
        CA[3, 5] = -globals.m55 * q
        CA[4, 0] = -globals.m33 * w
        CA[4, 2] = globals.m11 * u
        CA[4, 3] = -globals.m66 * rr
        CA[4, 5] = globals.m44 * p
        CA[5, 0] = globals.m22 * v
        CA[5, 1] = -globals.m11 * u
        CA[5, 3] = globals.m55 * q
        CA[5, 4] = -globals.m44 * p
        CA = -CA

        DV[0, 0] = Xuu * abs(u)
        DV[1, 1] = Yvv * abs(v)
        DV[2, 2] = Zww * abs(w)
        DV[3, 3] = Kpp * abs(p)
        DV[4, 4] = Mqq * abs(q)
        DV[5, 5] = Nrr * abs(rr)
        DV[4, 0] = Muu * abs(u)
        DV[3, 1] = Kvv * abs(v)

        CCD = CRB + CA - DV

        VV[0] = u
        VV[1] = v
        VV[2] = w
        VV[3] = p
        VV[4] = q
        VV[5] = rr

        VVV = np.dot(CCD, VV)
        c[0] = -(VVV[0]) + globals.Xprop + globals.X_HS
        c[1] = -(VVV[1]) + globals.Yprop + globals.Y_HS
        c[2] = -(VVV[2]) + globals.Zprop + globals.Z_HS
        c[3] = -(VVV[3]) + globals.Kprop + globals.fK_HS
        c[4] = -(VVV[4]) + globals.Mprop + globals.fM_HS
        c[5] = -(VVV[5]) + globals.Nprop + globals.fN_HS
        
        add_dp = globals.damping

        fluc = np.dot(add_dp, VV) #trong he toa do dia phuong

        fluc1[0] = fluc[0]
        fluc1[1] = fluc[1]
        fluc1[2] = fluc[2]
        globals.fmt[0:3, 0:3] = globals.fmt_2[0:3, 0:3]
        fluc1 = np.dot(globals.fmt, fluc1)  #trong he toa do toan cuc
        
        globals.vv_f[0:3] = VV[0:3]
        globals.vv_f = np.dot(globals.fmt, globals.vv_f)
        globals.cd_luc = abs(fluc1[0]) / (0.5 * globals.Apyz * globals.ro * globals.vv_f[0] * globals.vv_f[0])
        #tinh ve phai
        
        cc = np.dot(globals.bb, c)
        #dydt vi phan
        dydt[0] = u
        dydt[1] = v
        dydt[2] = w
        dydt[3] = p
        dydt[4] = q
        dydt[5] = rr
        dydt[6] = cc[0]
        dydt[7] = cc[1]
        dydt[8] = cc[2]
        dydt[9] = cc[3]
        dydt[10] = cc[4]
        dydt[11] = cc[5]

        return dydt
    #giai vi phan, dau ra la ma tran y_out
    Yy = odeint(uwqteta_cham_6_Luc, yy1, [0, globals.deltaT])

    x = Yy[1, 0]
    y = Yy[1, 1]
    z = Yy[1, 2]
    p0 = Yy[1, 3]
    q0 = Yy[1, 4]
    r0 = Yy[1, 5]
    u = Yy[1, 6]
    v = Yy[1, 7]
    w = Yy[1, 8]
    p = Yy[1, 9]
    q = Yy[1, 10]
    rr = Yy[1, 11]

    # Viet ra file
    #fileID, fileID2

    with open('lucCanTongHop.txt', 'w') as fid1:
        fid1.write(str(globals.vv_f[0]) + '\n')
        fid1.write(str(globals.cd_luc) + '\n')
    
    with open('ketquara_diaphuong.txt', 'w') as fid2:
        fid2.write(str(x) + '\n')
        fid2.write(str(y) + '\n')
        fid2.write(str(z) + '\n')
        fid2.write(str(p0) + '\n')
        fid2.write(str(q0) + '\n')
        fid2.write(str(r0) + '\n')
        fid2.write(str(u) + '\n')
        fid2.write(str(v) + '\n')
        fid2.write(str(w) + '\n')
        fid2.write(str(p) + '\n')
        fid2.write(str(q) + '\n')
        fid2.write(str(rr) + '\n')

    return (
        u, v, w, p, q, rr, x, y, z, p0, q0, r0
    )
#tinh van toc toan cuc
def tinhvantoctoancuc(phi, teta, si, u, v, w, u_i, v_i, w_i, xxx, yyy, zzz, p, q, rr, vantoc_phi, vantoc_teta, vantoc_si, p0, 
q0, r0, c2, c4):
    fmt = np.zeros((3, 3))

    c1 = np.zeros((3, 1))
    c1[0, 0] = u
    c1[1, 0] = v
    c1[2, 0] = w
    (phi, si, teta, fmt) = matranquanhe1_giaipt(phi, si, teta, fmt) #matran J1(eta2)  cong thuc 3.2
    c2 = np.dot(fmt, c1) #vecto van toc trong toa do toan cuc
    u_i = c2[0, 0]
    v_i = c2[1, 0]
    w_i = c2[2, 0]

    c1[0, 0] = xxx #dich chuyen trong dia phuong
    c1[1, 0] = yyy
    c1[2, 0] = zzz
    c1 = np.dot(fmt, c1) 
    xxx = c1[0, 0] #dich chuyen trong toan cuc
    yyy = c1[1, 0]
    zzz = c1[2, 0]

    c1[0, 0] = p
    c1[1, 0] = q
    c1[2, 0] = rr
    (phi, si, teta, fmt) = matranquanhe2_giaipt(phi, si, teta, fmt)
    c4 = np.dot(fmt, c1)
    vantoc_phi = c4[0, 0] #van toc goc trong toa do toan cuc moi
    vantoc_teta = c4[1, 0]
    vantoc_si = c4[2, 0]

    globals.phi_teta_si[0, 0] = c4[0, 0]
    globals.phi_teta_si[1, 0] = c4[1, 0]
    globals.phi_teta_si[2, 0] = c4[2, 0]

    c1[0, 0] = p0
    c1[1, 0] = q0
    c1[2, 0] = r0
    c3 = np.dot(fmt, c1) #delta goc tren he toa do toan cuc
    phi += c3[0, 0]
    teta += c3[1, 0]
    si += c3[2, 0]

    return (
        phi, teta, si, u, v, w, u_i, v_i, w_i, xxx, yyy, zzz, p, q, rr, vantoc_phi, vantoc_teta, vantoc_si, p0, q0, r0, c2, c4
    )
#van toc dia phuong moi
def vantocdiaphuongmoi(u, v, w, p, q, rr, c2, c4, phi, si, teta):
    fmt = np.zeros((3, 3))
    (phi, si, teta, fmt) = matranquanhe1_giaipt(phi, si, teta, fmt) #ct 3.2
    fmt = np.transpose(fmt)
    c3 = np.dot(fmt, c2)
    u = c3[0, 0] #u,v,w o toa do dia phuong moi
    v = c3[1, 0]
    w = c3[2, 0]

    (phi, si, teta, fmt) = matranquanhe2_giaipt(phi, si, teta, fmt) #ct 3.5
    fmt = np.linalg.inv(fmt)
    c1 = np.dot(fmt, c4) #p,q,r trong toa do dia phuong moi
    p = c1[0, 0]
    q = c1[1, 0]
    rr = c1[2, 0]

    return (
        u, v, w, p, q, rr, c2, c4, phi, si, teta
    )

def matranquanhe2_giaipt(phi, si, teta, fmt):
    fmt[0, 0] = 1
    fmt[0, 1] = np.sin(phi) * np.tan(teta)
    fmt[0, 2] = np.cos(phi) * np.tan(teta)
    fmt[1, 0] = 0
    fmt[1, 1] = np.cos(phi)
    fmt[1, 2] = -np.sin(phi)
    fmt[2, 0] = 0
    fmt[2, 1] = np.sin(phi) / np.cos(teta)
    fmt[2, 2] = np.cos(phi) / np.cos(teta)

    return phi, si, teta, fmt

def matranquanhe1_giaipt(phi, si, teta, fmt):
    global fmt_2_v2
    fmt_2_v2 = fmt[0:3, 0:3]

    fmt[0, 0] = np.cos(si) * np.cos(teta)
    fmt[0, 1] = -np.sin(si) * np.cos(phi) + np.cos(si) * np.sin(teta) * np.sin(phi)
    fmt[0, 2] = np.sin(si) * np.sin(phi) + np.cos(si) * np.sin(teta) * np.cos(phi)
    fmt[1, 0] = np.sin(si) * np.cos(teta)
    fmt[1, 1] = np.cos(si) * np.cos(phi) + np.sin(si) * np.sin(teta) * np.sin(phi)
    fmt[1, 2] = -np.cos(si) * np.sin(phi) + np.sin(si) * np.sin(teta) * np.cos(phi)
    fmt[2, 0] = -np.sin(teta)
    fmt[2, 1] = np.cos(teta) * np.sin(phi)
    fmt[2, 2] = np.cos(teta) * np.cos(phi)
    fmt_2_v2[0:3, 0:3] = fmt[0:3, 0:3]

    return phi, si, teta, fmt