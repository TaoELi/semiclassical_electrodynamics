'''
Given a vector field in real space, this script can calculate the transverse component
of this field by doing 2D FFT twice.
'''

import numpy as np
import math

# Given the expression of Px(R), Py(R) and Pz(R), I need to calculate the Fourier component
def calc_transverse_components(size = [5.0, 5.0], dx = 0.5, sigma = 1.0, mu12 = 0.10, local_size=100):

    sigma0 = sigma
    sigma = 2.5 * sigma
    x = np.arange(-local_size/2.0, local_size/2.0, dx)
    y = np.arange(-local_size/2.0, local_size/2.0, dx)
    [X, Y] = np.meshgrid(x, y)

    R2 = X * X + Y * Y

    Pz = (mu12 / (2.0 * np.pi)**(1.5) / sigma**2) * (2*np.pi)**(0.5) * np.exp(-R2 / 2.0 /sigma**2)
    Pz0 = (mu12 / (2.0 * np.pi)**(1.5) / sigma0**2) * (2*np.pi)**(0.5) * np.exp(-R2 / 2.0 /sigma0**2)

    start_idx_x, end_idx_x = int((local_size - size[0])/2/dx) , int((local_size + size[0])/2/dx)
    start_idx_y, end_idx_y = int((local_size - size[1])/2/dx) , int((local_size + size[1])/2/dx)

    Pz =   Pz[start_idx_x:end_idx_x, start_idx_y:end_idx_y]
    Pz0 =   Pz0[start_idx_x:end_idx_x, start_idx_y:end_idx_y]
    Pz = np.reshape(Pz0, (math.ceil(size[0]/dx), math.ceil(size[1]/dx), 1))

    Pz = np.copy(Pz.astype(np.complex128), order='C')

    return Pz