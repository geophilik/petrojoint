import sys

import numpy as np

import pygimli as pg
from pji import PetMod, JointInv, JointMod
from pybert.manager import ERTManager
from pygimli.physics import Refraction

# Settings
if len(sys.argv) > 1:
    scenario = "pest"
    poro = 0.3  # startmodel if poro is estimated
    fix_poro = False
    poro_min = 0.15
    poro_max = 0.45
    lam = 10
    case = 2
else:
    scenario = "pfix"
    fix_poro = True
    poro_min = 0
    poro_max = 1
    lam = 20
    case = 1

############
# Settings
maxIter = 15
############

# Poro to rock content (inversion parameter)
fr_min = 1 - poro_max
fr_max = 1 - poro_min

# Load meshes and data
mesh = pg.load("mesh.bms")
true = np.load("true_model.npz")
sensors = np.load("sensors.npy", allow_pickle=True)

veltrue, rhotrue, fatrue, fwtrue = true["vel"], true["rho"], true[
    "fa"], true["fw"]

ertScheme = pg.DataContainerERT("erttrue.dat")

meshRST = pg.load("paraDomain.bms")
meshERT = pg.load("meshERT.bms")

if fix_poro:
    frtrue = np.load("true_model.npz")["fr"]
    phi = []
    for cell in meshRST.cells():
        idx = mesh.findCell(cell.center()).id()
        phi.append(1 - frtrue[idx])
    phi = np.array(phi)
    fr_min = 0
    fr_max = 1
else:
    phi = poro

pm = PetMod(phi=phi)

# Setup managers and equip with meshes
ert = ERTManager()
ert.setMesh(meshERT)
ert.setData(ertScheme)
ert.fop.createRefinedForwardMesh()

rst = Refraction("tttrue.dat", verbose=True)
ttData = rst.dataContainer
rst.setMesh(meshRST, secNodes=3)

# Setup joint modeling and inverse operators
JM = JointMod(meshRST, ert, rst, pm, fix_poro=fix_poro, zWeight=0.25)

data = pg.cat(ttData("t"), ertScheme("rhoa"))
error = pg.cat(rst.relErrorVals(ttData), ertScheme("err"))

# Set gradient starting model of f_ice, f_water, f_air = phi/3
velstart = np.loadtxt("rst_startmodel.dat")
rhostart = np.ones_like(velstart) * np.mean(ertScheme("rhoa"))
fas, fws, _ = pm.all(rhostart, velstart)
frs = np.ones_like(fas) - pm.phi
if not fix_poro:
    frs[frs <= fr_min] = fr_min + 0.01
    frs[frs >= fr_max] = fr_max - 0.01
startmodel = np.concatenate((fws, fas, frs))
# Fix small values to avoid problems in first iteration
startmodel[startmodel <= 0.01] = 0.01

inv = JointInv(JM, data, error, startmodel, lam=lam, frmin=fr_min,
               frmax=fr_max, maxIter=maxIter)
inv.setModel(startmodel)

# Run inversion
model = inv.run()
pg.boxprint(("Chi squared fit:", inv.getChi2()), sym="+")

# Save results
fwe, fae, fre = JM.fractions(model)
fsum = fwe + fae + fre

print("Min/Max sum:", min(fsum), max(fsum))

rhoest = JM.pm.rho(fwe, fae, fre)
velest = 1. / JM.pm.slowness(fwe, fae, fre)

array_mask = np.array(((fae < 0) | (fae > 1 - fre))
                      | ((fwe < 0) | (fwe > 1 - fre))
                      | ((fre < 0) | (fre > 1))
                      | (fsum > 1.01))

np.savez("joint_inversion_%s.npz" % scenario, vel=np.array(velest),
         rho=np.array(rhoest), fa=fae, fw=fwe, fr=fre, mask=array_mask)

print("#" * 80)
ertchi, _ = JM.ERTchi2(model, error)
rstchi, _ = JM.RSTchi2(model, error, ttData("t"))
print("ERT chi^2", ertchi)
print("RST chi^2", rstchi)
print("#" * 80)

fig = JM.showFit(model)
title = "Overall chi^2 = %.2f" % inv.getChi2()
title += "\nERT chi^2 = %.2f" % ertchi
title += "\nRST chi^2 = %.2f" % rstchi
fig.suptitle(title)
fig.savefig("datafit_%s.png" % case, dpi=150)
