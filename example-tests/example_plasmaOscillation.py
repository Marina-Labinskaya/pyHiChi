import sys
sys.path.append("../bin/")
import pyHiChi as hichi

import numpy as np
import math as ma

import matplotlib.pyplot as plt
import matplotlib.animation as animation

nx = 64
ny = 8
nz = 8
L = 1.0
dx = L / nx
dy = dx
dz = dx
Nip = 256
Np = 2
N = Nip * Np
T0 = 0.01 * hichi.ELECTRON_MASS * hichi.c * hichi.c
D = T0 / (8 * hichi.pi * hichi.ELECTRON_CHARGE * hichi.ELECTRON_CHARGE * 0.25 * dx * dx)
wp = ma.sqrt(4.0 * hichi.pi * hichi.ELECTRON_CHARGE * hichi.ELECTRON_CHARGE * D / hichi.ELECTRON_MASS)
dt = 2 * (hichi.pi / wp) / Nip
Nc = 30
w = D * dx * dx * dx / Nc
A = 0.05

field_size = hichi.Vector3d(nx, ny, nz)
min_coords = hichi.Vector3d(0.0, 0.0, 0.0)
field_step = hichi.Vector3d(dx, dy, dz)

field = hichi.YeeField(field_size, min_coords, field_step, dt)

x = np.arange(min_coords.x, L, dx)

def valueEx(x, y, z):
    v1 = hichi.Vector3d(x, y, z)
    v2 = hichi.Vector3d(min_coords.x + v1.x, min_coords.y + v1.y, min_coords.z + v1.z)
    Ex = -2 * L * D * A * hichi.ELECTRON_CHARGE * ma.cos(2 * hichi.pi * v2.x / L)
    return Ex

def valueEy(x, y, z):
    Ey = 0
    return Ey

def valueEz(x, y, z):
    Ez = 0 
    return Ez

def valueBx(x, y, z):
    Bx = 0
    return Bx

def valueBy(x, y, z):
    By = 0  
    return By

def valueBz(x, y, z):
    Bz = 0
    return Bz

def DensityFunc(x, y, z):
    T0 = 0.01 * hichi.ELECTRON_MASS * hichi.c * hichi.c
    L = 1.0
    dx = L / 64.0
    D = T0 / (8 * hichi.pi * hichi.ELECTRON_CHARGE * hichi.ELECTRON_CHARGE * 0.25 * dx * dx)
    return D * (1 + 0.05 * ma.sin(2 * hichi.pi * x / L))

def InitialTemperatureFunc(x, y, z):
    T0 = 0.01 * hichi.ELECTRON_MASS * hichi.c * hichi.c
    return T0

def InitialMomentumFunc(x, y, z):
    return hichi.Vector3d(0.0, 0.0, 0.0)

# def ExArray(i, Ex):
#     if i % 16 == 0:
#         res = []
#         for j in range(nx):
#             #idx = hichi.Vector3d(0, 0, 0)
#             #internalCoords = hichi.Vector3d(0.0, 0.0, 0.0)
#             res.append(Ex)
#     return res

def get_fields():
    global field, x, nx
    #print(field)
    Ex = np.zeros(shape=(nx))
    Ey = np.zeros(shape=(nx))
    Ez = np.zeros(shape=(nx))
    Bx = np.zeros(shape=(nx))
    By = np.zeros(shape=(nx))
    Bz = np.zeros(shape=(nx))
    for ix in range(nx):
        coord_x = hichi.Vector3d(x[ix], 0.0, 0.0)
        E = field.get_E(coord_x)
        Ex[ix] = E.x
        Ey[ix] = E.y
        Ez[ix] = E.z
        B = field.get_B(coord_x)
        Bx[ix] = B.x
        By[ix] = B.y
        Bz[ix] = B.z
    return Ex, Ey, Ez, Bx, By, Bz

# def update_data():
#     for i in range(1000):
#         field.update_fields()

field.set_E(valueEx, valueEy, valueEz)
field.set_B(valueBx, valueBy, valueBz)
field.set_periodical_BC()

pusher = hichi.BorisPusher()
currentDeposition = hichi.FirstOrderCurrentDepositionYee()

J_BC = hichi.periodical_J_BC()
particle_BC = hichi.periodical_particle_BC()

p_array = hichi.ParticleArray()

particle_generator = hichi.ParticleGenerator()
particle_generator(p_array, field, DensityFunc, InitialTemperatureFunc, InitialMomentumFunc, w, hichi.ELECTRON) # field != Grid*

for i in range (N + 1):
    (Ex, Ey, Ez, Bx, By, Bz) = get_fields()
    # fdtd
    field.update_fields()
    # interpolate
    fields = []
    for j in range(p_array.size()):
        E = hichi.Vector3d()
        B = hichi.Vector3d()
        #field.getFieldsCIC(particleArray[i].getPosition(), E, B) ?
        fields.append(hichi.FieldValue(E, B))
    # pusher
    pusher(p_array, fields, dt)
    # periodical particle BC
    particle_BC.updateParticlePosition(field, p_array, dt)
    # current deposition
    # fdtd.updateDims() ?
    currentDeposition(field, p_array) # field != Grid*
    J_BC.updateCurrentBondaries()

    


# for (int i = 0; i <= N; ++i) {
#     if (i % 16 == 0)
#         for (int j = 0; j < grid.numInternalCells.x; ++j) {
#             Int3 idx; FP3 internalCoords;
#             fout << i << " " << minCoords.x + j * grid.steps.x << " " << grid.Ex(j + grid.getNumExternalLeftCells().x, grid.getNumExternalLeftCells().y, grid.getNumExternalLeftCells().z) << std::endl;
#             //fout << i << " " << minCoords.x + j * grid.steps.x << " " << grid.Ex(j, grid.getNumExternalLeftCells().y, grid.getNumExternalLeftCells().z) << std::endl;
#         };
#     if (i % 16 == 0) {
#         for (int k = 0; k < grid.numInternalCells.x; ++k) {
#             int electronCount = 0;
#             FP minCellCoords = minCoords.x + k * grid.steps.x;
#             FP maxCellCoords = minCoords.x + (k + 1) * grid.steps.x;
#             //FP minCellCoords = grid.origin.x + k * grid.steps.x;
#             //FP maxCellCoords = grid.origin.x + (k + 1) * grid.steps.x;

#             for (int j = 0; j < particleArray.size(); ++j) {
#                 if ((particleArray[j].getPosition().x >= minCellCoords) && (particleArray[j].getPosition().x <= maxCellCoords))
#                     electronCount++;
#             }
#             // to sinchronize the output with Picador
#             FP electronDensity = electronCount * w;  // density through plane
#             if (i % 16 == 0)
#                 fout2 << i << " " << minCellCoords << " " << electronDensity << std::endl;
#         }
#     }
# }
