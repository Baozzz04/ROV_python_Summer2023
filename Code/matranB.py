import numpy as np

def matranB(add, rb, matran_B):
    for i in range(6):
        for j in range(6):
            matran_B[i, j] = add[i, j] + rb[i, j]
    
    return add, rb, matran_B