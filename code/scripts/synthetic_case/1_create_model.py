import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pji import PetMod

# create and save sensors
snum = 72; sspc = 3; sstartx = 20
sensors = np.arange(sstartx, sstartx + snum * sspc, sspc, dtype="float")
print(sensors)
print("Sensors", len(sensors))
sensors.dump("sensors.npy")

# Model creation
world = mt.createWorld([0, -40], [253, 0], layers=[-1 , -25], worldMarker=False)

# landfill body
lfb = mt.createPolygon([[35, -1], [218, -1], [203, -18], [50, -18]],
                         isClosed=True)
                         
# preferential flow path
pfp = mt.createPolygon([[90, -1], [95, -1], [85, -18], [70, -18]],
                         isClosed=True)
                         
geom = mt.mergePLC([world, lfb, pfp])
geom.addRegionMarker((10, -.5), 0)
geom.addRegionMarker((10, -10), 1)
geom.addRegionMarker((10, -30), 2)
geom.addRegionMarker((50, -10), 3)
geom.addRegionMarker((150, -10), 3)
geom.addRegionMarker((75, -15), 4)
geom.save("geom.bms")

# ~ ax, _ = pg.show(geom, markers=True)
# ~ ax.plot(sensors, np.zeros_like(sensors), 'ko')

# ~ pg.plt.show()

mesh = mt.createMesh(geom, area=1.0); pg.plt.show()

#~ pg.show(mesh, markers=True)

# Model creation based on pore fractions
philayers = np.array([0.9, 0.25, 0.05, 0.45, 0.5])
frlayers = 1 - philayers
fwlayers = np.array([0.4, 0.25, 0.03, 0.2, 0.3])
falayers = philayers - fwlayers

falayers[np.isclose(falayers, 0.0)] = 0.0
print(falayers)

# ~ # Save for covariance calculations
# ~ Fsyn = np.vstack((fwlayers, filayers, falayers, frlayers))
# ~ np.savetxt("syn_model.dat", Fsyn)

pm = PetMod(phi=philayers, vw=1500., va=330., vr=5000, a=1.2, n=1.5, m=2,
            rhow=100.)

print(falayers + fwlayers + frlayers)
rholayers = pm.rho(fwlayers, falayers, frlayers)
vellayers = 1. / pm.slowness(fwlayers, falayers, frlayers)

print(rholayers)
print(vellayers)


def to_mesh(data):
    return data[mesh.cellMarkers()]


rhotrue = to_mesh(rholayers)
veltrue = to_mesh(vellayers)

# %%
# Save true model and mesh
fa = to_mesh(falayers)
fw = to_mesh(fwlayers)
fr = to_mesh(frlayers)

pm.fr = fr
pm.phi = 1 - fr
#~ pm.show(mesh, rhotrue, veltrue)

assert np.allclose(fa + fw + fr, 1)

np.savez("true_model.npz", rho=rhotrue, vel=veltrue, fa=fa, fw=fw,
         fr=fr)

mesh.save("mesh.bms")
np.savetxt("rhotrue.dat", rhotrue)
np.savetxt("veltrue.dat", veltrue)

#~ pg.plt.show()
