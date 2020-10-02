import sys
import numpy as np
import pygimli as pg
from pji import PetMod

mesh = pg.load("mesh.bms")

if len(sys.argv) > 1:
    pd = pg.load("paraDomain.bms")
    resinv = np.loadtxt("res_conventional.dat")
    vest = np.loadtxt("vel_conventional.dat")
    scenario = "phom"
    pg.boxprint(scenario)
    phi = 0.3  # Porosity assumed to calculate fi, fa, fw with 4PM
else:
    pd = pg.load("paraDomain.bms")
    resinv = np.loadtxt("res_conventional.dat")
    vest = np.loadtxt("vel_conventional.dat")
    scenario = "phet"
    pg.boxprint(scenario)
    frtrue = np.load("true_model.npz")["fr"]

    phi = []
    for cell in pd.cells():
        idx = mesh.findCell(cell.center()).id()
        phi.append(1 - frtrue[idx])
    phi = np.array(phi)

# Save some stuff
fpm = PetMod(phi=phi, vw=1500., va=330., vr=5000, a=1.2, n=1.5, m=2,
            rhow=100.)
fae, fwe, maske = fpm.all(resinv, vest)
print(np.min(fwe), np.max(fwe))
np.savez("conventional_%s.npz" % scenario, vel=np.array(vest),
         rho=np.array(resinv), fa=fae, fw=fwe, mask=maske)
