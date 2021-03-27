import matplotlib.pyplot as plt
import numpy as np

import pybert as pb
import pygimli as pg


class JointMod(pg.ModellingBase):
    def __init__(self, mesh, ertlofop, erthifop, srtfop, petromodel, fix_poro=True,
                 zWeight=1, verbose=True, corr_l=None, fix_water=False,
                 fix_air=False, fix_cec=False):
        """Joint petrophysical modeling operator.

        Parameters
        ----------
        mesh : pyGIMLi mesh
        ertlofop : ERT low forward operator
        erthifop : ERT high forward operator
        srtfop : SRT forward operator
        petromodel : Petrophysical four-phase model
        zWeight : zWeight for more or less layering
        verbose : Be more verbose
        corr_l : tuple
            Horizontal and vertical correlation lengths. If provided,
            geostatistical regularization will be used and classical smoothing
            with zWeight will be ignored.
        fix_poro|water|air : boolean or vector
            Fix to starting model or provide weight vector for particular cells.
        """
        pg.ModellingBase.__init__(self, verbose)
        self.mesh = pg.Mesh(mesh)
        self.ERTlo = ertlofop
        self.ERThi = erthifop
        self.SRT = srtfop
        self.fops = [self.SRT, self.ERThi, self.ERTlo]
        self.pm = petromodel
        self.cellCount = self.mesh.cellCount()
        self.fix_water = fix_water
        self.fix_air = fix_air
        # ~ self.fix_cec = fix_cec
        self.fix_poro = fix_poro
        self.zWeight = zWeight
        # self.fix_cells = fix_cells
        self.corr_l = corr_l
        self.createConstraints()

    def fractions(self, model):
        """Split model vector into individual distributions"""
        return np.reshape(model, (4, self.cellCount))

    def createJacobian(self, model):
        #~ print('*' * 30)
        #~ print('here')
        #~ print('*' * 30)
        # ~ fw, fa, cec, fr = self.fractions(model)
        fw, fa, fr, cec = self.fractions(model)
        
        # ~ print(fw, fa, cec, fr)

        rholo = self.pm.rholo(fw, fa, cec, fr)
        rhohi = self.pm.rhohi(fw, fa, cec, fr)
        s = self.pm.slowness(fw, fa, fr)
        
        # ~ print(rholo, rhohi, s)

        self.ERTlo.fop.createJacobian(rholo)
        #~ print('*' * 30)
        #~ print('ERTlo done')
        #~ print('*' * 30)
        self.ERThi.fop.createJacobian(rhohi)
        #~ print('*' * 30)
        #~ print('ERThi done')
        #~ print('*' * 30)
        self.SRT.fop.createJacobian(s)
        #~ print('*' * 30)
        #~ print('SRT done')
        #~ print('*' * 30)
        
        jacERTlo = self.ERTlo.fop.jacobian()
        jacERThi = self.ERThi.fop.jacobian()
        jacSRT = self.SRT.fop.jacobian()

        # Setting inner derivatives
        self.jacSRTW = pg.MultRightMatrix(jacSRT, r=1. / self.pm.vw)
        self.jacSRTA = pg.MultRightMatrix(jacSRT, r=1. / self.pm.va)
        # ~ self.jacSRTC = pg.MultRightMatrix(jacSRT, r=0)
        self.jacSRTR = pg.MultRightMatrix(jacSRT, r=1. / self.pm.vr)

        self.jacERTloW = pg.MultRightMatrix(
            jacERTlo, r=self.pm.rholo_deriv_fw(fw, fa, cec, fr))
        self.jacERTloA = pg.MultRightMatrix(
            jacERTlo, r=self.pm.rholo_deriv_fa(fw, fa, cec, fr))
        # ~ self.jacERTloC = pg.MultRightMatrix(
            # ~ jacERTlo, r=self.pm.rholo_deriv_cec(fw, fa, cec, fr))
        self.jacERTloR = pg.MultRightMatrix(
            jacERTlo, r=self.pm.rholo_deriv_fr(fw, fa, cec, fr))

        self.jacERThiW = pg.MultRightMatrix(
            jacERThi, r=self.pm.rhohi_deriv_fw(fw, fa, cec, fr))
        self.jacERThiA = pg.MultRightMatrix(
            jacERThi, r=self.pm.rhohi_deriv_fa(fw, fa, cec, fr))
        # ~ self.jacERThiC = pg.MultRightMatrix(
            # ~ jacERThi, r=self.pm.rhohi_deriv_cec(fw, fa, cec, fr))
        self.jacERThiR = pg.MultRightMatrix(
            jacERThi, r=self.pm.rhohi_deriv_fr(fw, fa, cec, fr))

        # Putting subjacobians together in block matrix
        self.jac = pg.BlockMatrix()
        nData = 0
        self.jac.addMatrix(self.jacSRTW, nData, 0)
        self.jac.addMatrix(self.jacSRTA, nData, self.cellCount)
        # ~ self.jac.addMatrix(self.jacSRTC, nData, self.cellCount * 2)
        self.jac.addMatrix(self.jacSRTR, nData, self.cellCount * 2)
        nData += self.SRT.fop.data().size()  # update total vector length
        self.jac.addMatrix(self.jacERTloW, nData, 0)
        self.jac.addMatrix(self.jacERTloA, nData, self.cellCount)
        # ~ self.jac.addMatrix(self.jacERTloC, nData, self.cellCount * 2)
        self.jac.addMatrix(self.jacERTloR, nData, self.cellCount * 2)
        nData += self.ERTlo.fop.data().size()
        self.jac.addMatrix(self.jacERThiW, nData, 0)
        self.jac.addMatrix(self.jacERThiA, nData, self.cellCount)
        # ~ self.jac.addMatrix(self.jacERThiC, nData, self.cellCount * 2)
        self.jac.addMatrix(self.jacERThiR, nData, self.cellCount * 2)
        
        #~ print('*' * 30)
        #~ print('and then here')
        #~ print('*' * 30)
        
        self.setJacobian(self.jac)

    def createConstraints(self):
        # First order smoothness matrix
        self._Ctmp = pg.RSparseMapMatrix()

        if self.corr_l is None:
            pg.info("Using smoothing with zWeight = %.2f." % self.zWeight)
            rm = self.SRT.fop.regionManager()
            rm.fillConstraints(self._Ctmp)

            # Set zWeight
            rm.setZWeight(self.zWeight)
            self.cWeight = pg.RVector()
            rm.fillConstraintsWeight(self.cWeight)
            self._CW = pg.LMultRMatrix(self._Ctmp, self.cWeight)
        else:
            pg.info("Using geostatistical constraints with " + str(self.corr_l))
            # Geostatistical constraints by Jordi et al., GJI, 2018
            CM = pg.utils.geostatistics.covarianceMatrix(self.mesh, I=self.corr_l)
            self._Ctmp = pg.matrix.Cm05Matrix(CM)
            self._CW = self._Ctmp

        # Putting together in block matrix
        self._C = pg.RBlockMatrix()
        cid = self._C.addMatrix(self._CW)
        self._C.addMatrixEntry(cid, 0, 0)
        self._C.addMatrixEntry(cid, self._Ctmp.rows(), self.cellCount)
        # ~ self._C.addMatrixEntry(cid, self._Ctmp.rows() * 2, self.cellCount * 2)
        self._C.addMatrixEntry(cid, self._Ctmp.rows() * 3, self.cellCount * 2)
        self.setConstraints(self._C)

        # Identity matrix for interparameter regularization
        self._I = pg.IdentityMatrix(self.cellCount)

        self._G = pg.RBlockMatrix()
        iid = self._G.addMatrix(self._I)
        self._G.addMatrixEntry(iid, 0, 0)
        self._G.addMatrixEntry(iid, 0, self.cellCount)
        # ~ self._G.addMatrixEntry(iid, 0, self.cellCount * 2)
        self._G.addMatrixEntry(iid, 0, self.cellCount * 2)

        self.fix_val_matrices = {}
        # Optionally fix phases to starting model globally or in selected cells
        phases = ["water", "air", "cec", "rock matrix"]
        for i, phase in enumerate([self.fix_water, self.fix_air,# self.fix_cec,
                                   self.fix_poro]):
            name = phases[i]
            vec = pg.RVector(self.cellCount)
            if phase is True:
                pg.info("Fixing %s content globally." % name)
                vec += 1.0
            elif hasattr(phase, "__len__"):
                pg.info("Fixing %s content at selected cells." % name)
                phase = np.asarray(phase, dtype="int")
                vec[phase] = 1.0
            self.fix_val_matrices[name] = pg.matrix.DiagonalMatrix(vec)
            self._G.addMatrix(self.fix_val_matrices[name],
                              self._G.rows(), self.cellCount * i)

    def showModel(self, model):
        # ~ fw, fa, cec, fr = self.fractions(model)
        fw, fa, fr, cec = self.fractions(model)

        rholo = self.pm.rholo(fw, fa, cec, fr)
        rhohi = self.pm.rhohi(fw, fa, cec, fr)
        s = self.pm.slowness(fw, fa, fr)

        _, axs = plt.subplots(3, 3)
        pg.show(self.mesh, fw, ax=axs[0, 0], label="Water content", hold=True,
                logScale=False, cMap="Blues")
        pg.show(self.mesh, fa, ax=axs[2, 0], label="Air content", hold=True,
                logScale=False, cMap="Greens")
        pg.show(self.mesh, fr, ax=axs[2, 1], label="Rock matrix content",
                hold=True, logScale=False, cMap="Oranges")
        pg.show(self.mesh, rholo, ax=axs[0, 1], label="Rho low", hold=True,
                cMap="Spectral_r")
        pg.show(self.mesh, rhohi, ax=axs[0, 2], label="Rho high", hold=True,
                cMap="Spectral_r")
        pg.show(self.mesh, 1 / s, ax=axs[1, 1], label="Velocity")

    def showFit(self, model):
        resp = self.response(model)

        fig, axs = plt.subplots(2, 3, figsize=(10, 10))
        t_resp = resp[:self.SRT.dataContainer.size()]
        rholoa_resp = resp[self.SRT.dataContainer.size():self.SRT.dataContainer.size()+self.ERTlo.data.size()]
        rhohia_resp = resp[self.SRT.dataContainer.size() + self.ERTlo.data.size():]
        self.SRT.showData(response=t_resp, ax=axs[0, 0])

        t_fit = t_resp - self.SRT.dataContainer("t")
        lim = np.max(np.abs(t_fit))
        axs[0, 0].set_title("Traveltime curves with fit")
        axs[1, 0].set_title("Deviation between traveltimes")
        self.SRT.showVA(vals=t_fit, ax=axs[1, 0], cMin=-lim, cMax=lim,
                        cmap="RdBu_r")

        rholoa_fit = (self.ERTlo.data("rhoa") - rholoa_resp) / rholoa_resp * 100
        lim = np.max(np.abs(rholoa_fit))
        pb.show(self.ERTlo.data, ax=axs[0, 1], label=r"Measured data $\rho_a$")
        pb.show(self.ERTlo.data, vals=rhoa_fit, cMin=-lim, cMax=lim,
                label="Relative fit (%%)", cMap="RdBu_r", ax=axs[1, 1])

        rhohia_fit = (self.ERThi.data("rhoa") - rhohia_resp) / rhohia_resp * 100
        lim = np.max(np.abs(rhohia_fit))
        pb.show(self.ERThi.data, ax=axs[0, 2], label=r"Measured data $\rho_a$")
        pb.show(self.ERThi.data, vals=rhoa_fit, cMin=-lim, cMax=lim,
                label="Relative fit (%%)", cMap="RdBu_r", ax=axs[1, 2])
        fig.show()
        return fig

    def ERTlochi2(self, model, error):  # chi2 and relative rms for the rhoa data
        resp = self.response(model)
        resprholoa = resp[self.SRT.dataContainer.size():self.SRT.dataContainer.size()+self.ERTlo.data.size()]
        rholoaerr = error[self.SRT.dataContainer.size():self.SRT.dataContainer.size()+self.ERTlo.data.size()]
        chi2rholoa = pg.utils.chi2(self.ERTlo.data("rhoa"), resprholoa, rholoaerr)
        rmsrholoa = pg.rrms(self.ERTlo.data("rhoa"), resprholoa)
        return chi2rholoa, rmsrholoa

    def ERThichi2(self, model, error):  # chi2 and relative rms for the rhoa data
        resp = self.response(model)
        resprhohia = resp[self.SRT.dataContainer.size() + self.ERTlo.data.size():]
        rhohiaerr = error[self.SRT.dataContainer.size() + self.ERTlo.data.size():]
        chi2rhohia = pg.utils.chi2(self.ERThi.data("rhoa"), resprhohia, rhohiaerr)
        rmsrhohia = pg.rrms(self.ERThi.data("rhoa"), resprhohia)
        return chi2rhohia, rmsrhohia

    def SRTchi2(self, model, error,
                data):  # chi2 and relative rms for the travel time data
        resp = self.response(model)
        resptt = resp[:self.SRT.dataContainer.size()]
        tterr = error[:self.SRT.dataContainer.size()]
        chi2tt = pg.utils.chi2(data, resptt, tterr)
        rmstt = np.sqrt(np.mean((resptt - data)**2))
        return chi2tt, rmstt

    def response(self, model):
        return self.response_mt(model)

    def response_mt(self, model, i=0):
        model = np.nan_to_num(model)
        fw, fa, fr, cec = self.fractions(model)

        rholo = self.pm.rholo(fw, fa, cec, fr)
        rhohi = self.pm.rhohi(fw, fa, cec, fr)
        s = self.pm.slowness(fw, fa, fr)

        print("=" * 30)
        print("        Min. | Max.")
        print("-" * 30)
        print(" Water: %.2f | %.2f" % (np.min(fw), np.max(fw)))
        print(" Air:   %.2f | %.2f" % (np.min(fa), np.max(fa)))
        print(" Rock:  %.2f | %.2f" % (np.min(fr), np.max(fr)))
        print("-" * 30)
        print(" SUM:   %.2f | %.2f" % (np.min(fa + fw + fr),
                                       np.max(fa + fw + fr)))
        print("-" * 30)
        print(" CEC: %.2f | %.2f" % (np.min(cec), np.max(cec)))
        print("-" * 30)
        print("=" * 30)
        print(" Rho:   %.2e | %.2e" % (np.min(rholo), np.max(rholo)))
        print(" Rho:   %.2e | %.2e" % (np.min(rhohi), np.max(rhohi)))
        print(" Vel:   %d | %d" % (np.min(1 / s), np.max(1 / s)))

        t = self.SRT.fop.response(s)
        rholoa = self.ERTlo.fop.response(rholo)
        rhohia = self.ERThi.fop.response(rhohi)

        return pg.cat(t, pg.cat(rholoa, rhohia))
