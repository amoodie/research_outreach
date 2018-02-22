# channel functions
import numpy as np


def set_B(B0, mou, thet, nx, dx):
    # set B for channel width
    B = np.zeros(nx + 1)
    chanLen = int(mou * nx)
    B[:chanLen] = B0
    B[chanLen:] = B0 + 2 * ((np.arange(1, (nx - chanLen+2))) \
        * (1600e3 / 400) * np.tan(np.radians(2)))
    # print(B)
    return B

def get_slope(eta, nx, dx):
    # return slope of input bed (eta)
    S = np.zeros(nx + 1)
    S[0] = (eta[0] - eta[1]) / dx
    S[1:nx] = (eta[0:nx-1] - eta[2:nx+1]) / (2 * dx)
    S[nx:] = (eta[nx-1] - eta[nx]) / dx
    return S