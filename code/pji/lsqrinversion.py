from math import sqrt

import numpy as np

import pygimli as pg
from pygimli.utils import boxprint

from .mylsqr import lsqr


class LSQRInversion(pg.RInversion):
    """LSQR solver based inversion"""

    def __init__(self, *args, **kwargs):
        """Init."""
        pg.RInversion.__init__(self, *args, **kwargs)
        self.G = None
        self.c = None
        self.my = 1.0

    def setParameterConstraints(self, G, c, my=1.0):
        """Set parameter constraints G*p=c."""
        self.G = G
        self.c = c
        self.my = my

    def run(self, **kwargs):
        """Run."""
        print("model", min(self.model()), max(self.model()))
        for i in range(self.maxIter()):
            boxprint("Iteration #%d" % i, width=80, sym="+")
            self.oneStep()
            boxprint("Iteration: %d | Chi^2: %.2f | RMS: %.2f%%" % (i,
                     self.chi2(), self.relrms()))
            if self.chi2() <= 1.0:
                print("Done. Reached target data misfit of chi^2 <= 1.")
                break
            # ~ phi = self.getPhi()
            model = self.model()
            phi = self.getPhi(model[:self.fop().cellCount * 3],
                              self.fop().response(model))
            if i > 2:
                print("Phi / oldphi", phi / oldphi)
            if (i > 10) and (phi / oldphi >
                             (1 - self.deltaPhiAbortPercent() / 100)):
                print("Done. Reached data fit criteria of delta phi < %.2f%%."
                      % self.deltaPhiAbortPercent())
                break
            if i + 1 == self.maxIter():
                print("Done. Maximum number of iterations reached.")
            oldphi = phi
        return self.model()

    def oneStep(self):
        """One inversion step."""
        model = self.model()
        
        # update cec
        model = self.fop().updateCEC(model)
        self.setModel(model)
        
        # ~ print('#' * 30)
        # ~ print('begin of oneStep()')
        # ~ print(model)
        # ~ print('#' * 30)
        if len(self.response()) != len(self.data()):
            self.setResponse(self.forwardOperator().response(model))

        self.forwardOperator().createJacobian(model)
        self.checkTransFunctions()
        tD = self.transData()
        tM = self.transModel()
        nData = self.data().size()
        #        nModel = len(model)
        self.A = pg.BlockMatrix()  # to be filled with scaled J and C matrices
        # part 1: data part
        J = self.forwardOperator().jacobian()
        # self.dScale = 1.0 / pg.log(self.error()+1.0)
        self.dScale = 1.0 / (
            tD.deriv(self.data()) * self.error() * self.data())
        self.leftJ = tD.deriv(self.response()) * self.dScale
        #        self.leftJ = self.dScale / tD.deriv(self.response())
        self.rightJ = 1.0 / tM.deriv(model[:self.fop().cellCount * 3])
        # ~ print("#" * 30)
        # ~ print("J",J.cols(), J.rows())
        # ~ print("leftJ",len(self.leftJ))
        # ~ print("rightJ",len(self.rightJ))
        # ~ print("#" * 30)
        self.JJ = pg.matrix.MultLeftRightMatrix(J, self.leftJ, self.rightJ)
        #        self.A.addMatrix(self.JJ, 0, 0)
        self.mat1 = self.A.addMatrix(self.JJ)
        self.A.addMatrixEntry(self.mat1, 0, 0)
        # part 2: normal constraints
        self.checkConstraints()
        self.C = self.forwardOperator().constraints()
        self.leftC = pg.RVector(self.C.rows(), 1.0)
        self.rightC = pg.RVector(self.C.cols(), 1.0)
        self.CC = pg.matrix.MultLeftRightMatrix(self.C, self.leftC,
                                                self.rightC)
        self.mat2 = self.A.addMatrix(self.CC)
        lam = self.getLambda()
        self.A.addMatrixEntry(self.mat2, nData, 0, sqrt(lam))
        # % part 3: parameter constraints
        if self.G is not None:
            self.rightG = 1.0 / tM.deriv(model[:self.fop().cellCount * 3])
            # ~ tmp = 1.0 / tM.deriv(model)
            # ~ tmp[self.fop().cellCount*2:self.fop().cellCount*3] = 1.
            # ~ self.rightG = tmp
            # ~ print(self.rightG[self.fop().cellCount*2:self.fop().cellCount*3])
            self.GG = pg.matrix.MultRightMatrix(self.G, self.rightG)
            self.mat3 = self.A.addMatrix(self.GG)
            nConst = self.C.rows()
            self.A.addMatrixEntry(self.mat3, nData + nConst, 0, sqrt(self.my))
        self.A.recalcMatrixSize()
        # right-hand side vector
        deltaD = (tD.fwd(self.data()) - tD.fwd(self.response())) * self.dScale
        deltaC = -(self.CC * tM.fwd(model[:self.fop().cellCount * 3]) * sqrt(lam))
        deltaC *= 1.0 - self.localRegularization()  # operates on DeltaM only
        rhs = pg.cat(deltaD, deltaC)
        if self.G is not None:
            deltaG = (self.c - self.G * model[:self.fop().cellCount * 3]) * sqrt(self.my)
            rhs = pg.cat(pg.cat(deltaD, deltaC), deltaG)

        dM = lsqr(self.A, rhs)
        # ~ dM = pg.cat(dM, np.ones(self.fop().cellCount) * model[self.fop().cellCount * 3:] * -1)
        dM = pg.cat(dM, np.ones(self.fop().cellCount))
        tau, responseLS = self.lineSearchInter(dM, model)
        if tau < 0.1:  # did not work out
            tau = self.lineSearchQuad(dM, responseLS)
        if tau > 0.9:  # save time and take 1
            tau = 1.0
        else:
            # ~ self.forwardOperator().response(self.model())
            self.forwardOperator().response(model)

        if tau < 0.1:  # still not working
            tau = 0.1  # try a small value

        # ~ self.setModel(tM.update(self.model(), dM * tau))
        umodel = tM.update(model, dM * tau)
        umodel = pg.cat(umodel[:self.fop().cellCount * 3], model[self.fop().cellCount * 3:])
        
        # ~ # update cec
        # ~ umodel = self.fop().updateCEC(umodel)
        
        # set updated model
        self.setModel(umodel)
        # ~ self.setModel(tM.update(model, dM * tau))
        # ~ print('#' * 30)
        # ~ print('after update')
        # ~ print(self.model())
        # ~ print('#' * 30)        
        # print("model", min(self.model()), max(self.model()))
        if tau == 1.0:
            self.setResponse(responseLS)
        else:  # compute new response
            self.setResponse(self.forwardOperator().response(self.model()))

        self.setLambda(self.getLambda() * self.lambdaFactor())
        return True

    def lineSearchInter(self, dM, model, nTau=100):
        """Optimizes line search parameter by linear response interpolation."""
        tD = self.transData()
        tM = self.transModel()
        # ~ model = self.model()
        response = self.response()
        # ~ print(model)
        # ~ print(dM)
        modelLS = tM.update(model, dM)
        modelLS = pg.cat(modelLS[:self.fop().cellCount * 3], model[self.fop().cellCount * 3:])
        responseLS = self.forwardOperator().response(modelLS)
        taus = np.linspace(0.0, 1.0, nTau)
        phi = np.ones_like(taus) * self.getPhi()
        # ~ print('#' * 30)
        # ~ print(modelLS)
        # ~ print(responseLS)
        # ~ print('#' * 30)
        phi[-1] = self.getPhi(modelLS[:self.fop().cellCount * 3], responseLS)
        t0 = tD.fwd(response)
        t1 = tD.fwd(responseLS)
        for i in range(1, len(taus) - 1):
            tau = taus[i]
            modelI = tM.update(model, dM * tau)
            responseI = tD.inv(t1 * tau + t0 * (1.0 - tau))
            phi[i] = self.getPhi(modelI[:self.fop().cellCount * 3], responseI)

        return taus[np.argmin(phi)], responseLS

    def lineSearchQuad(self, dM, responseLS):
        """Optimize line search by fitting parabola by Phi(tau) curve."""
        return 0.1


if __name__ == '__main__':
    nlay = 4  # number of layers
    lam = 200.  # (initial) regularization parameter
    errPerc = 3.  # relative error of 3 percent
    ab2 = np.logspace(-1, 2, 50)  # AB/2 distance (current electrodes)
    mn2 = ab2 / 3.  # MN/2 distance (potential electrodes)
    f = pg.DC1dModelling(nlay, ab2, mn2)
    synres = [100., 500., 20., 800.]  # synthetic resistivity
    synthk = [0.5, 3.5, 6.]  # synthetic thickness (nlay-th layer is infinite)
    rhoa = f(synthk + synres)
    rhoa = rhoa * (pg.randn(len(rhoa)) * errPerc / 100. + 1.)
    transLog = pg.RTransLog()
    inv = LSQRInversion(rhoa, f, transLog, transLog, True)
    inv.setRelativeError(errPerc / 100)
    startModel = pg.cat(
        pg.Vector(nlay - 1, 5), pg.Vector(nlay, pg.median(rhoa)))
    print(inv.response())
    inv.setModel(startModel)
    inv.setMarquardtScheme()
    inv.setLambda(1000)
    G = pg.RMatrix(rows=1, cols=len(startModel))
    for i in range(3):
        G[0][i] = 1
    c = pg.Vector(1, pg.sum(synthk))
    inv.setParameterConstraints(G, c, 100)
    # print("Start", inv.chi2(), inv.relrms(), pg.sum(inv.model()(0, nlay-1)))
    if 0:
        for i in range(10):
            inv.oneStep()
            print(i, inv.chi2(), inv.relrms(), pg.sum(inv.model()(0,
                                                                  nlay - 1)))
    else:
        inv.run()
        print(inv.chi2(), inv.relrms(), pg.sum(inv.model()(0, nlay - 1)))
