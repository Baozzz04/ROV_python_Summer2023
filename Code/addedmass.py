import numpy as np
import globals
#�nh h? s? kh?i l??ng th�m Added mass
def addedmass(m11_3D, m22, m33, m44, m55, m66, Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1, add):
    m11_3DS = 0
    landa = 0
    m11 = 0
    
    Axy = np.pi * (globals.l1**2) / 4 #Di?n t�ch tham chi?u c?a v?t th? r?n l�n XOY
    Axz = np.pi * (globals.l1**2) / 4 #Di?n t�ch tham chi?u c?a v?t th? r?n l�n XOZ
    Ayz = np.pi * (globals.r1**2) #Di?n t�ch tham chi?u c?a v?t th? r?n l�n YOZ
    Vr = 2 * np.pi * globals.r1**2 * globals.l1 / 3 #Th? t�ch tham chi?u c?a v?t th? r?n
    Cpxy = globals.Apxy / Axy
    Cpxz = globals.Apxz / Axz
    Cpyz = globals.Apyz / Ayz
    #T�nh h? s? m11 trong kh�ng gian 3 chi?u
    m11_3D = Ca1 * 2 * np.pi * globals.r1**2 * globals.l1 * globals.ro * (Cpyz**2) * Cpxz * Cpxy / 3
    m11 = m11_3D
    #T�nh h? s? m11_3DS theo ph??ng ph�p l� thuy?t m?nh
    m11_3DS = 4 * np.pi * globals.r1**3 * globals.ro * Cpyz / 3
    #T�nh h? s? t? l? sai kh�c gi?a 3D v� Strip
    landa = m11_3D / m11_3DS
    #T�nh h? s? m44
    m44 = 0
    #T�nh h? s? m22
    m22 = np.pi * globals.l1**2 * globals.r1 * globals.ro * Ca * Cpxz * landa / 3
    #T�nh h? s? m55
    m55 = globals.ro * np.pi * Cpxz * globals.r1 * (globals.l1**2 - globals.r1**2) * landa / 24
    #T�nh h? s? m33
    m33 = globals.ro * np.pi * globals.l1**2 * globals.r1 * Cpxy * landa / 3
    #T�nh h? s? m66
    m66 = globals.ro * np.pi * Cpxy * globals.r1 * (globals.l1**2 - globals.r1**2) * landa / 24
    
    add[0, 0] = m11
    add[1, 1] = m22
    add[2, 2] = m33
    add[3, 3] = m44
    add[4, 4] = m55
    add[5, 5] = m66
    
    # mo file He so khoi luong them
    with open('He so khoi luong them.txt', 'w') as fid:
        fid.write(str(m11) + '\n')
        fid.write(str(m22) + '\n')
        fid.write(str(m33) + '\n')
        fid.write(str(m44) + '\n')
        fid.write(str(m55) + '\n')
        fid.write(str(m66) + '\n')

    return (
        m11_3D, m22, m33, m44, m55, m66, Cdc, r2, l2, r3, l3, l4, r5, l5, b, b1, H5, x1, x2, Ca, Ca1, add
    )