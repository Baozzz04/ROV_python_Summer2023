import sympy as sp
import numpy as np
import globals

def Tinhtichphan(ro, Cd4, Cdc, aa, bb, aa2, bb2, hsmmc4, hsmmc1):
    #global kq, r1, b2, r4, kq2, l1, c2
    #global r4_m, r1_m, b2_m, c2_m, l1_m

    xx = sp.symbols('xx')
    
    globals.r4_m = globals.r4
    globals.b2_m = globals.b2
    globals.r1_m = globals.r1
    globals.l1_m = globals.l1
    #globals.c2_m = globals.c2
    kq2 = int(0)
    
    x_tmp = sp.symbols('x_tmp')
    f3 = sp.sqrt(sp.Abs(globals.r4_m**2 - sp.sqrt(sp.Abs(globals.r1_m + globals.b2_m + globals.r4_m - x_tmp))**2)) * x_tmp**3
    kq = sp.integrate(f3, (x_tmp, aa, bb))
    kq = -kq * ro * Cd4
    hsmmc4 = kq
    
    s = globals.l1_m / 2 - globals.r1_m

    f4 = sp.Piecewise(
        (sp.sqrt(sp.Abs((xx - s)**2 - globals.r1_m**2)) * xx**2 * sp.Abs(xx), (xx > s) & (xx < globals.l1_m/2)),
        (globals.r1_m * xx**2 * sp.Abs(xx), (xx < s) & (xx > -s)),
        (sp.sqrt(sp.Abs((xx + s)**2 - globals.r1_m**2)) * xx**2 * sp.Abs(xx), True) 
    )

    hsmmc1 = (-1) * sp.integrate(f4, (xx, aa2, bb2)) * ro * Cdc
    
    return ro, Cd4, Cdc, aa, bb, aa2, bb2, hsmmc4, hsmmc1
