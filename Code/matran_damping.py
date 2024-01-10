import numpy as np
import globals
import Tinhtichphan

#Ch??ng tr�nh t�nh c�c h? s? l?c c?n
def matran_damping(DP, Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1):
    a, c, c1, c2 = 0, 0, 0, 0
    Cd1, Cd2, Cd3, Cd4, Cd5, aa, bb, aa2, bb2 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    cd4, cdc = 0, 0
    hsmmc4, hsmmc1 = 0, 0
    Afx1, Afx2, Afx3, Afx4, Afx5 = 0, 0, 0, 0, 0
    Afy1, Afy2, Afy3, Afy4, Afy5 = 0, 0, 0, 0, 0
    Afz1, Afz2, Afz3, Afz4, Afz5 = 0, 0, 0, 0, 0
    Ap1, Ap2, Ap3, Ap4, Ap5 = 0, 0, 0, 0, 0
    m11_3D = 0
    m22, m33, m44, m55, m66, m11, add = 0, 0, 0, 0, 0, 0, 0
    dp, i, j = 0, 0, 0

    [Cd1, Cd2, Cd3, Cd4, Cd5, Afx1, Afx2, Afx3, Afx4, Afx5, Ap4, Cdc, r2, l2, r3, l3, l4, r5, l5, globals.Ixx, globals.Iyy, globals.Izz, globals.Ixy, globals.Ixz, globals.Iyz, globals.Apxy, globals.Apxz, globals.Apyz, b, b1, H5, x1, x2, Ca, Ca1] = Hesocandoc(Cd1, Cd2, Cd3, Cd4, Cd5, Afx1, Afx2, Afx3, Afx4, Afx5, Ap4, Cdc, r2, l2, r3, l3, l4, r5, l5, globals.Ixx, globals.Iyy, globals.Izz, globals.Ixy, globals.Ixz, globals.Iyz, globals.Apxy, globals.Apxz, globals.Apyz, b, b1, H5, x1, x2, Ca, Ca1)
    #T�nh di?n t�ch chi?u theo ph??ng OY
    Afy1 = np.pi * globals.r1**2 + 2 * globals.r1 * (globals.l1 - 2 * globals.r1)
    Afy2 = 2 * r2 * l2
    Afy3 = 2 * r3 * l3
    Afy4 = 2 * globals.r4 * l4
    Afy5 = 2 * r5 * l5
    #T�nh di?n t�ch chi?u theo ph??ng OZ
    Afz1 = np.pi * globals.r1**2 + 2 * globals.r1 * (globals.l1 - 2 * globals.r1)
    Afz2 = 2 * r2 * l2
    Afz3 = 2 * r3 * l3
    Afz4 = 3.13 * globals.r4**2
    Afz5 = 2 * r5 * l5
    
    a = globals.r1
    c = 2 * r2
    c1 = 2 * r3
    c2 = 2 * globals.r4
    aa = a + globals.b2
    bb = a + globals.b2 + c2
    aa2 = -globals.l1 / 2
    bb2 = globals.l1 / 2
    #T�nh h? s? l?c c?n Xuu
    Xuu1 = -globals.ro * Cd1 * Afx1 / 2
    Xuu2 = -globals.ro * Cd2 * Afx2 / 2
    Xuu3 = -globals.ro * Cd3 * Afx3 / 2
    Xuu4 = -globals.ro * Cdc * Ap4 / 2
    Xuu5 = -globals.ro * Cd5 * Afx5 / 2
    Xuu = Xuu1 + 2 * Xuu2 + 2 * Xuu3 + 2 * Xuu4 + 2 * Xuu5
    #T�nh h? s? l?c c?n Yvv
    Yvv1 = -globals.ro * Cdc * Afy1 / 2
    Yvv2 = -globals.ro * Cdc * Afy2 / 2
    Yvv3 = -globals.ro * Cdc * Afy3 / 2
    Yvv4 = -globals.ro * Cdc * Afy4 / 2
    Yvv5 = -globals.ro * Cdc * Afy5 / 2
    Yvv = Yvv1 + 2 * Yvv2 + 2 * Yvv3 + 2 * Yvv4 + 2 * Yvv5
    #T�nh h? s? l?c c?n Zww
    Zww1 = -globals.ro * Cdc * Afz1 / 2
    Zww2 = -globals.ro * Cdc * Afz2 / 2
    Zww3 = -globals.ro * Cdc * Afz3 / 2
    Zww4 = -globals.ro * Cd4 * Afz4 / 2
    Zww5 = -globals.ro * Cdc * Afz5 / 2
    Zww = Zww1 + 2 * Zww2 + 2 * Zww3 + 2 * Zww4 + 2 * Zww5
    #T�nh m� men c?n quay quanh tr?c OX: Kpp
    Kpp1 = 0
    Kpp2 = -globals.ro * Cdc * l2 * ((a + b + c)**4 - (a + b)**4) / 8
    Kpp3 = -globals.ro * Cdc * l3 * ((a + b1 + c1)**4 - (a + b1)**4) / 8
    #T�nh Kpp4 do ??ng c? th?ng ??ng g�y ra
    (globals.ro, Cd4, Cdc, aa, bb, aa2, bb2, hsmmc4, hsmmc1) = Tinhtichphan.Tinhtichphan(globals.ro, cd4, Cdc, aa, bb, aa2, bb2, hsmmc4, 
    hsmmc1)
    Kpp4 = hsmmc4
    #T�nh Kpp5 do thanh tr?ng l??ng g�y ra
    Kpp5 = Yvv5 * H5**3
    Kpp = Kpp1 + 2 * Kpp2 + 2 * Kpp3 + 2 * Kpp4 + 2 * Kpp5
    #T�nh m� men c?n quay quanh tr?c OY: Mqq
    Mqq1 = -(globals.ro * Cdc * globals.r1 * (globals.l1 / 2)**4) / 2
    Mqq2 = -(globals.ro * Cdc * r2 * ((x1 + l2)**4 + x1**4)) / 4
    Mqq3 = -(globals.ro * Cdc * r3 * ((x2 + l3)**4 + x2**4)) / 4
    Mqq4 = -(globals.ro * Cdc * globals.r4 * (l4 / 2)**4) / 2
    Mqq5 = -(globals.ro * Cdc * r5 * (l5 / 2)**4) / 2
    Mqq = Mqq1 + 2 * Mqq2 + 2 * Mqq3 + 2 * Mqq4 + 2 * Mqq5
    #T�nh m� men c?n quay quanh tr?c OZ: Nrr
    Nrr1 = hsmmc1
    Nrr2 = -globals.ro * Cdc * r2 * ((x1 + l2)**4 - x1**4) / 4
    Nrr3 = -globals.ro * Cdc * r3 * ((x2 + l3)**4 - x2**4) / 4
    Nrr4 = -globals.ro * Cdc * l4 * ((a + b + 2 * globals.r4)**4 - (a + b)**4) / 8
    Nrr5 = -globals.ro * Cdc * r5 * l5**4 / 32
    Nrr = Nrr1 + 2 * Nrr2 + 2 * Nrr3 + 2 * Nrr4 + 2 * Nrr5
    #T�nh h? s? m� men c?n Kvv quanh truc OX do v g�y ra
    Kvv = -2 * globals.ro * Cdc * r5 * l5 * H5
    #T�nh h? s? m� men c?n Muu quay quanh OY do u g�y ra
    Muu = -globals.ro * Cd5 * np.pi * r5**2 * H5
    #viet file vao Hesoluccan
    with open('Hesoluccan.txt', 'w') as fid:
        fid.write(str(Xuu) + '\n')
        fid.write(str(Yvv) + '\n')
        fid.write(str(Zww) + '\n')
        fid.write(str(Kpp) + '\n')
        fid.write(str(Mqq) + '\n')
        fid.write(str(Nrr) + '\n')
        fid.write(str(Kvv) + '\n')
        fid.write(str(Muu) + '\n')
    
    DP = np.zeros((6, 6))
    DP[0, 0] = Xuu
    DP[1, 1] = Yvv
    DP[2, 2] = Zww
    DP[3, 3] = Kpp
    DP[4, 4] = Mqq
    DP[5, 5] = Nrr
    DP[3, 1] = Kvv
    DP[4, 0] = Muu
    
    return DP, Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1
#He so can doc
def Hesocandoc(Cd1, Cd2, Cd3, Cd4, Cd5, Afx1, Afx2, Afx3, Afx4, Afx5, Ap4, Cdc, r2, l2, r3, l3, l4, r5, l5, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Apxy, Apxz, Apyz, b, b1, H5, x1, x2, Ca, Ca1):
    #T�nh di?n t�ch chi?u theo ph??ng OX
    Afx1 = np.pi * globals.r1**2
    Afx2 = np.pi * r2**2
    Afx3 = np.pi * r3**2
    Afx4 = 2 * globals.r4 * l4
    Afx5 = np.pi * r5**2
    #T�nh di?n t�ch m?t ph?ng chi?u Ap
    Ap1 = 2 * globals.r1 * globals.l1
    Ap2 = 2 * r2 * l2
    Ap3 = 2 * r3 * l3
    Ap4 = 2 * globals.r4 * l4
    Ap5 = 2 * r5 * l5
    
    Cd1 = globals.Css * np.pi * Ap1 * (1 + 60 * ((2 * globals.r1) / globals.l1)**3 + 0.0025 * (globals.l1 / (2 * globals.r1))) / Afx1
    Cd2 = globals.Css * np.pi * Ap2 * (1 + 60 * ((2 * r2) / l2)**3 + 0.0025 * (l2 / (2 * r2))) / Afx2
    Cd3 = globals.Css * np.pi * Ap3 * (1 + 60 * ((2 * r3) / l3)**3 + 0.0025 * (l3 / (2 * r3))) / Afx3
    Cd4 = 0.0
    Cd5 = globals.Css * np.pi * Ap5 * (1 + 60 * ((2 * r5) / l5)**3 + 0.0025 * (l5 / (2 * r5))) / Afx5
    #mo file hesocandoc1
    with open('hesocandoc1.txt', 'w') as fid:
        fid.write(str(Cd1) + '\n')
        fid.write(str(Cd2) + '\n')
        fid.write(str(Cd3) + '\n')
        fid.write(str(Cd4) + '\n')
        fid.write(str(Cd5) + '\n')
    
    return Cd1, Cd2, Cd3, Cd4, Cd5, Afx1, Afx2, Afx3, Afx4, Afx5, Ap4, Cdc, r2, l2, r3, l3, l4, r5, l5, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Apxy, Apxz, Apyz, b, b1, H5, x1, x2, Ca, Ca1