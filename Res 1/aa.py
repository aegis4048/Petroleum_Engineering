# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:18:51 2018

@author: TK
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

data = pd.read_excel('LonghornPVT1.xlsx')

# create regression lines
p = np.polyfit(data['P(Psia)'][:13], data['Bo (RB/STB)'][:13], deg=2)
j = np.polyfit(data['P(Psia)'][14:], data['Bo (RB/STB)'][14:], deg=2)

plt.scatter(data['P(Psia)'], data['Bo (RB/STB)'])
x = data['P(Psia)'][:13]
X = data['P(Psia)'][14:]
plt.plot(x, p[0] * x ** 2 + p[1] * x + p[2])
plt.plot(X, j[0] * X ** 2 + j[1] * X + j[2])
plt.xlabel('P(Psia)')
plt.ylabel('Bo (RB/STB)')
plt.title('Pressure vs. Bo')
plt.show()

wells = pd.read_excel('LonghornReservoir.xlsx')
press = wells['Well Pressure (psi)']
wells['Bo (RB/STB)'] = j[0] * press ** 2 + j[1] * press + j[2]

plt.scatter(wells["Well-X (ft)"], wells["Well-Y (ft)"], c=wells["Porosity"])
plt.show()
'''
%  The function is to interpolate netpay, porosity and Bo by kriging method

%  Input: 
% - well location; 
% - well data;
% Output:
% - Contour Plot of Depth

% % Kriging Method:
% - Spherical Model
% - lambda = inv(K)*M2, where lambda is the weight vector

% % Kriging Parameters: 
% Nugget - variance at zero distance
% Range - the distance at which the semivariogram levels off and beyond
%         which the semivariance is constant
% Sill - the constant semivariance value beyond the range
'''


def gamma_comp(C0, a, sill, h):
    # Variogram Function

    # Spherical Model
    if h == 0:
        y = 0
    elif h > a:
        y = C0 + sill
    else:
        y = C0 + sill * (1.5 * h / a - 0.5 * (h / a) ** 3)
    return y


def krige(wells):
    ReservoirGeometry = wells

    # Well Location
    well_x = wells['Well-X (ft)']
    well_y = wells['Well-Y (ft)']
    # Net Pay, ft
    well_netpay = wells['Well-Y (ft)']
    # Porosity
    well_phi = wells['Porosity']
    # Formation Volume Factor
    well_Bo = wells['Bo (RB/STB)']

    num_well = len(well_x)  # The Number of Wells
    # Grid Size
    dx = 55
    dy = 55  # ft

    # Number of grid blocks
    num_grid_x = 106
    num_grid_y = 69

    # Variagram Function
    C0 = 2  # Nugget
    a = 2000  # Range
    sill = 20  # Sill

    # Matrix Construction
    # K-Matrix
    K = np.zeros((num_well + 1, num_well + 1))  # Initialization
    for i in range(1, num_well - 1):
        for j in range((i + 1), num_well):
            h = math.sqrt((well_x[i - 1] - well_x[j - 1]) ** 2 + (
                        well_y[i - 1] - well_y[j - 1]) ** 2)  # Distance between point-i to point-j
            K[i - 1, j - 1] = C0 + sill - gamma_comp(C0, a, sill, h)  # Cij
            K[j - 1, i - 1] = K[i - 1, j - 1]  # Cji = Cij

    for i in range(1, num_well):
        K[i - 1, i - 1] = C0 + sill
        K[i - 1, num_well] = 1
        K[num_well, i - 1] = 1

    K[num_well - 1, num_well - 1] = 0;
    # Grids
    X = np.linspace(dx / 2, dx, num=num_grid_x * dx - dx / 2)
    Y = np.linspace(dy / 2, dy, num=num_grid_y * dy - dy / 2)
    [x, y] = np.meshgrid(X, Y)

    # Compute the netpay, porosity and Bo at each grid block
    for i in range(1, np.size((x, 1))):
        for j in range(1, np.size((x, 2))):
            if ReservoirGeometry[i - 1, j - 1] == 1:
                for k in range(1, num_well):
                    h = math.sqrt((x[i - 1, j - 1] - well_x[k - 1]) ** 2 + (
                                y[i - 1, j - 1] - well_y[k - 1]) ** 2)  # Distance between point-0 to point-k
                    # M2-Matrix
                    M2[k - 1, 0] = C0 + sill - gamma_comp(C0, a, sill, h)

                M2[num_well, 0] = 1
                lambd = numpy.linalg.lstsq(k.T, M2.T)[0].T
                netpay[i - 1, j - 1] = well_netpay * lambd[1, end - 1]  # Net Pay at Each Grid
                phi[i - 1, j - 1] = well_phi * lambd[1, end - 1]  # Porosity at Each Grid
                Bo[i - 1, j - 1] = well_Bo * lambd[1, end - 1]  # Bo at Each Grid
            else:
                netpay[i - 1, j - 1] = 0
                phi[i - 1, j - 1] = 0
                Bo[i - 1, j - 1] = 0

    # Transpose
    netpay = np.transpose(netpay)
    phi = np.transpose(phi)
    Bo = np.transpose(Bo)
    return (netpay, phi, Bo)