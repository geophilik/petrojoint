import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg


class PPMod():

    def __init__(self, vw=1500., va=330., vr=5500, a=1., n=2., m=1.3,
                 phi=0.4, rhow=150., R=.1, B=3.1*1e-9, l=3.0*1e-10, rhog=2650):
        """Petrophysical model. Estimates fraction
        of air and water from electrical bulk resistivity and seismic
        velocity.

        Parameters
        ----------
        vw : float or array type
            Velocity of water in m/s (the default is 1500.).
        va : float or array type
            Velocity of air in m/s (the default is 330.).
        vr : float or array type
            Velocity of rock in m/s (the default is 5000).
        a : float or array type
            Archie parameter `a` (the default is 1).
        n : float or array type
            Archie parameter `n` (the default is 2).
        m : float or array type
            Archie parameter `m` (the default is 1.3).
        phi : float or array type
            Porosity `phi` (the default is 0.4).
        rhow : float or array type
            Water resistivity `rhow` (the default is 150).
        R : float
            Dimensionless number `R` relating normalized chargeability
            and surface conductivity (the default is 0.1).
        B : float
            Apparent mobility of counterions for surface conduction `B`
            in m**-2 s**-1 V**-1 (the default is 3.1 * 1e-9)
        l : float
            apparent mobility of counterions for the polarization
            `lambda` in m**-2 s**-1 V**-1 (the default is 3.0 * 1e-10)
        rhog : float (or array type)
            grain density `rhog` in kg/m**3 (the default is 2650)
        """

        # Velocities of water, air, ice and rock (m/s)
        self.vw = vw
        self.va = va
        self.vr = vr

        # Archie parameter
        self.a = a
        self.m = m
        self.n = n
        self.phi = phi
        self.fr = 1 - self.phi  # fraction of rock
        self.rhow = rhow
        
        # surface conductivity
        self.R = R
        self.B = B
        self.l = l
        self.rhog = rhog

    def water_old(self, rho):
        fw = ((self.a * self.rhow * self.phi**self.n) /
              (rho * self.phi**self.m))**(1. / self.n)
        fw[np.isclose(fw, 0)] = 0
        return fw

    def water(self, rholo, rhohi):
        sigmahi = 1. / rhohi
        sigmalo = 1. / rholo
        
        mn = sigmahi - sigmalo
        # ~ if np.mean(sigmahi) > np.mean(sigmalo):
            # ~ mn = sigmahi - sigmalo
        # ~ else:
            # ~ mn = np.abs(sigmahi - sigmalo)
        # ~ mn[mn < 0] = 0
        
        # ~ fw = (self.rhow * self.phi**self.n / self.phi**self.m * \
              # ~ (sigmahi - (mn) / self.R))**(1. / self.n)
        fw = (self.rhow * (self.phi**(self.n) / self.phi**(self.m)) * \
              (sigmahi - mn / self.R))**(1. / self.n)
        fw[np.isclose(fw, 0)] = 0
        return fw

    # ~ def cec(self, rholo, rhohi, fw):
        # ~ sigmahi = 1. / rhohi
        # ~ sigmalo = 1. / rholo
        
        # ~ mn = sigmahi - sigmalo
        # ~ # mn[mn < 0] = 0
        
        # ~ # cec = (self.phi / fw)**(self.n-1) * self.phi**(self.m-1) * (mn) / (self.rhog * self.l)
        # ~ cec = self.phi**(self.n - self.m) * mn / (fw**(self.n-1) * self.rhog * self.l)
        
        # ~ return cec

    def cec(self, rholo, rhohi):
        sigmahi = 1. / rhohi
        sigmalo = 1. / rholo
        
        mn = sigmahi - sigmalo
        print(mn)
        # ~ if np.mean(sigmahi) > np.mean(sigmalo):
            # ~ mn = sigmahi - sigmalo
        # ~ else:
            # ~ mn = np.abs(sigmahi - sigmalo)
        # ~ mn[mn < 0] = 0
        
        # ~ cec = (self.phi / fw)**(self.n-1) * self.phi**(self.m-1) * (mn) / (self.rhog * self.l)
        cec = self.phi**(self.n - self.m) * mn / (self.water(rholo, rhohi)**(self.n-1) * self.rhog * self.l)
        
        cec[np.isclose(cec, 0)] = 0
        
        return cec

    # ~ def air(self, rho, v):
        # ~ fa = ((self.vi * self.va / (self.vi - self.va) * (
            # ~ 1. / v - self.fr / self.vr - self.phi / self.vi - self.water(rho) *
            # ~ (1. / self.vw - 1. / self.vi))))
        # ~ fa[np.isclose(fa, 0)] = 0
        # ~ return fa
        
    def air(self, rholo, rhohi, v):
        fa = self.va * (1. / v - self.fr / self.vr - self.water(rholo, rhohi) / self.vw)
        fa[np.isclose(fa, 0)] = 0
        return fa
        
    # ~ def air(self, rho, v):
        # ~ fa = (2 - self.water(rho) + (1 + v / self.vw) - self.fr * (1 + v / self.vr)) / (1 + v / self.va)
        # ~ fa[np.isclose(fa, 0)] = 0
        # ~ return fa
        
    # ~ def air(self, rho, v):
        # ~ fa = (1 / 2) * (1 - self.water(rho) - self.fr + self.va * (1 / v - self.water(rho) / self.vw - self.fr / self.vr))
        # ~ fa[np.isclose(fa, 0)] = 0
        # ~ return fa

    def rho_old(self, fw, fa, fr=None):
        """Return electrical resistivity based on fraction of water `fw`."""
        if fr is None:
            phi = fw + fa
        else:
            phi = 1 - fr

        rho = self.a * self.rhow * phi**(-self.m) * (fw / phi)**(-self.n)
        if (rho <= 0).any():
            pg.warn(
                "Found negative resistivity, setting to nearest above zero.")
            rho[rho <= 0] = np.min(rho[rho > 0])
        return rho

    # ~ def rholo(self, fw, fa, cec, fr=None):
        # ~ if fr is None:
            # ~ phi = fw + fa
        # ~ else:
            # ~ phi = 1 - fr
        
        # ~ #~ print('+' * 30)
        # ~ #~ print(fw.mean())
        # ~ #~ print('+' * 30)
        # ~ #~ print('+' * 30)
        # ~ #~ print(fa.mean())
        # ~ #~ print('+' * 30)
        # ~ #~ print('+' * 30)
        # ~ #~ print(cec.mean())
        # ~ #~ print('+' * 30)
        # ~ #~ print('+' * 30)
        # ~ #~ print(fr.mean())
        # ~ #~ print('+' * 30)
        # ~ rho = phi**(self.n - self.m) * (fw**self.n * self.rhow**(-1) + \
                # ~ fw**(self.n-1) * self.rhog * (self.B - self.l) * cec)
        
        # ~ #~ print('+' * 30)
        # ~ #~ print(rho.mean())
        # ~ #~ print('+' * 30)
        
        # ~ if (rho <= 0).any():
            # ~ pg.warn(
                # ~ "Found negative resistivity, setting to nearest above zero.")
            # ~ rho[rho <= 0] = np.min(rho[rho > 0])
        # ~ return rho

    # ~ def rhohi(self, fw, fa, cec, fr=None):
        # ~ if fr is None:
            # ~ phi = fw + fa
        # ~ else:
            # ~ phi = 1 - fr
            
        # ~ rho = phi**(self.n - self.m) * (fw**self.n * self.rhow**(-1) + \
                # ~ fw**(self.n-1) * self.rhog * self.B * cec)
        
        # ~ if (rho <= 0).any():
            # ~ pg.warn(
                # ~ "Found negative resistivity, setting to nearest above zero.")
            # ~ rho[rho <= 0] = np.min(rho[rho > 0])
        # ~ return rho

    def rhohi(self, fw, fa, cec, fr=None):
        """Return high frequency electrical resistivity based on 
        fraction of water `fw` and cation exchange capacity `CEC`."""

        rho = 1. / self.sigmahi(fw, fa, cec, fr)
        if (rho <= 0).any():
            pg.warn(
                "Found negative resistivity, setting to nearest above zero.")
            rho[rho <= 0] = np.min(rho[rho > 0])
        return rho
    
    def sigmahi(self, fw, fa, cec, fr):
        if fr is None:
            phi = fw + fa
        else:
            phi = 1 - fr
            
        return (fw / phi)**self.n * phi**self.m * (1 / self.rhow) + \
               (fw / phi)**(self.n-1) * (phi**self.m / phi) * self.rhog * self.B * cec
    
    def rholo(self, fw, fa, cec, fr=None):
        """Return low frequency electrical resistivity based on 
        fraction of water `fw` and cation exchange capacity `CEC`."""

        rho = 1. / self.sigmalo(fw, fa, cec, fr)
        if (rho <= 0).any():
            pg.warn(
                "Found negative resistivity, setting to nearest above zero.")
            rho[rho <= 0] = np.min(rho[rho > 0])
        return rho
    
    def sigmalo(self, fw, fa, cec, fr):
        if fr is None:
            phi = fw + fa
        else:
            phi = 1 - fr

        return (fw / phi)**self.n * phi**self.m * (1 / self.rhow) + \
               (fw / phi)**(self.n-1) * (phi**self.m / phi) * self.rhog * (self.B - self.l) * cec
    
    # ~ def rho_deriv_fw(self, fw, fa, fr):
        # ~ return self.rho(fw, fa, fr) * -self.n / fw
    
    # ~ def rholo_deriv_fw(self, fw, fa, cec, fr):
        # ~ return self.rholo(fw, fa, cec, fr) * \
            # ~ (-self.n * fw**(self.n-1) * self.rhow**(-1) - (self.n-1) * fw**(self.n-2) * self.rhog * (self.B - self.l) * cec) / \
            # ~ (fw**self.n * self.rhow**(-1) + fw**(self.n-1) * self.rhog * (self.B - self.l) * cec)
    
    # ~ def rhohi_deriv_fw(self, fw, fa, cec, fr):
        # ~ return self.rhohi(fw, fa, cec, fr) * \
            # ~ (-self.n * fw**(self.n-1) * self.rhow**(-1) - (self.n-1) * fw**(self.n-2) * self.rhog * self.B * cec) / \
            # ~ (fw**self.n * self.rhow**(-1) + fw**(self.n-1) * self.rhog * self.B * cec)
    
    def rholo_deriv_fw(self, fw, fa, cec, fr):
        return ((1 / fw) * (-self.n * (fw / (1-fr))**self.n * (1-fr)**self.m * (1/self.rhog) - \
            (self.n-1) * (fw / (1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * (self.B - self.l) * cec)) / \
            self.sigmahi(fw, fa, cec, fr)**2
    
    def rhohi_deriv_fw(self, fw, fa, cec, fr):
        return ((1 / fw) * (-self.n * (fw / (1-fr))**self.n * (1-fr)**self.m * (1/self.rhog) - \
            (self.n-1) * (fw / (1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * self.B * cec)) / \
            self.sigmahi(fw, fa, cec, fr)**2

    # ~ def rho_deriv_fr(self, fw, fa, fr):
        # ~ return self.rho(fw, fa, fr) * (self.n - self.m) / (fr - 1)

    # ~ def rholo_deriv_fr(self, fw, fa, cec, fr):
        # ~ return self.rholo(fw, fa, cec, fr) * \
            # ~ (self.m - self.n) / (1 - fr)

    # ~ def rhohi_deriv_fr(self, fw, fa, cec, fr):
        # ~ return self.rhohi(fw, fa, cec, fr) * \
            # ~ (self.m - self.n) / (1 - fr)

    def rholo_deriv_fr(self, fw, fa, cec, fr):
        return ((1/(1-fr)) * (-(1-self.m) * (fw/(1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * (self.B - self.l) * cec - \
            (self.n-1)*(fw/(1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * (self.B - self.l) * cec + \
            self.m * (fw/(1-fr))**self.n * (1-fr)**self.m * (1/self.rhow)) - \
            self.n * (fw/(1-fr))**self.n * (1-fr)**self.m * (1/self.rhow)) / \
            self.sigmahi(fw, fa, cec, fr)**2

    def rhohi_deriv_fr(self, fw, fa, cec, fr):
        return ((1/(1-fr)) * (-(1-self.m) * (fw/(1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * self.B * cec - \
            (self.n-1)*(fw/(1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * self.B * cec + \
            self.m * (fw/(1-fr))**self.n * (1-fr)**self.m * (1/self.rhow)) - \
            self.n * (fw/(1-fr))**self.n * (1-fr)**self.m * (1/self.rhow)) / \
            self.sigmahi(fw, fa, cec, fr)**2
        
    # ~ def rho_deriv_fa(self, fw, fa, fr):
        # ~ return 0

    def rholo_deriv_fa(self, fw, fa, cec, fr):
        return 0

    def rhohi_deriv_fa(self, fw, fa, cec, fr):
        return 0

    # ~ def rholo_deriv_cec(self, fw, fa, cec, fr):
        # ~ return self.rholo(fw, fa, cec, fr) * \
            # ~ (-fw**(self.n-1) * self.rhog * (self.B - self.l)) / \
            # ~ (fw**self.n * self.rhow**(-1) + fw**(self.n-1) * self.rhog * (self.B - self.l) * cec)

    # ~ def rhohi_deriv_cec(self, fw, fa, cec, fr):
        # ~ return self.rhohi(fw, fa, cec, fr) * \
            # ~ (-fw**(self.n-1) * self.rhog * (self.B - self.l)) / \
            # ~ (fw**self.n * self.rhow**(-1) + fw**(self.n-1) * self.rhog * (self.B - self.l) * cec)

    def slowness(self, fw, fa, fr=None):
        """Return slowness based on fraction of water `fw` and air `fa`."""
        if fr is None:
            fr = 1 - (fw + fa)

        s = fw / self.vw + fr / self.vr + fa / self.va
        if (s <= 0).any():
            pg.warn("Found negative slowness, setting to nearest above zero.")
            s[s <= 0] = np.min(s[s > 0])
        return s

    def all(self, rholo, rhohi, v, mask=False):
        """Syntatic sugar for all fractions including a mask for unrealistic
        values."""

        # RVectors sometimes cause segfaults
        rholo = np.array(rholo)
        rhohi = np.array(rhohi)
        v = np.array(v)

        fa = self.air(rholo, rhohi, v)
        fw = self.water(rholo, rhohi)
        cec = self.cec(rholo, rhohi)

        # Check that fractions are between 0 and 1
        array_mask = np.array(((fa < 0) | (fa > 1 - self.fr))
                              | ((fw < 0) | (fw > 1 - self.fr))
                              | ((self.fr < 0) | (self.fr > 1)))
        if array_mask.sum() > 1:
            print("WARNING: %d of %d fraction values are unphysical." % (int(
                array_mask.sum()), len(array_mask.ravel())))

        if mask:
            fa = np.ma.array(fa, mask=(fa < 0) | (fa > 1 - self.fr))
            fw = np.ma.array(fw, mask=(fw < 0) | (fw > 1 - self.fr))

        return fa, fw, cec, array_mask

    def show(self, mesh, rho, vel, mask=True, **kwargs):
        fa, fw, mask = self.all(rho, vel, mask=mask)

        fig, axs = plt.subplots(3, 2, figsize=(16, 10))
        pg.show(mesh, fw, ax=axs[0, 0], label="Water content", hold=True,
                logScale=False, cMap="Blues", **kwargs)
        pg.show(mesh, fa, ax=axs[1, 0], label="Air content", hold=True,
                logScale=False, cMap="Greens", **kwargs)
        pg.show(mesh, rho, ax=axs[0, 1], label="Rho", hold=True,
                cMap="Spectral_r", logScale=True, **kwargs)
        pg.show(mesh, vel, ax=axs[1, 1], label="Velocity", logScale=False,
                hold=True, **kwargs)
        pg.show(mesh, self.phi, ax=axs[2, 0], label="Porosity", logScale=False,
                cMap='Oranges', hold=True, **kwargs)
        pg.show(mesh, 1 - self.phi, ax=axs[2, 1], label="Rock content", logScale=False,
                cMap='Oranges', hold=True, **kwargs)
        for ax in axs.flat:
            ax.set_facecolor("0.5")
        return fig, axs

class PetroMod():

    def __init__(self, vw=1500., va=330., vr=5500, vi=3750,
                 n=2., m=2.,
                 phi=0.4,
                 rhow=150., rhow_T0=21,
                 rhog=2650, 
                 # ~ cec=18.8, 
                 t=-4,
                 alpha_t=.02, tf=-1, tc=4,
                 R=.1, 
                 B=3.1*1e-9, B_T0=25,
                 l=3.0*1e-10, l_T0=25):
        """Petrophysical model. Estimates fraction
        of air and water from electrical bulk resistivity and seismic
        velocity.

        Parameters
        ----------
        vw : float or array type
            Velocity of water in m/s (the default is 1500.).
        va : float or array type
            Velocity of air in m/s (the default is 330.).
        vr : float or array type
            Velocity of rock in m/s (the default is 5000).
        n : float or array type
            Archie parameter `n` (the default is 2).
        m : float or array type
            Archie parameter `m` (the default is 1.3).
        phi : float or array type
            Porosity `phi` (the default is 0.4).
        rhow : float or array type
            Water resistivity `rhow` (the default is 150).
        R : float
            Dimensionless number `R` relating normalized chargeability
            and surface conductivity (the default is 0.1).
        B : float
            Apparent mobility of counterions for surface conduction `B`
            in m**-2 s**-1 V**-1 (the default is 3.1 * 1e-9)
        l : float
            apparent mobility of counterions for the polarization
            `lambda` in m**-2 s**-1 V**-1 (the default is 3.0 * 1e-10)
        rhog : float (or array type)
            grain density `rhog` in kg/m**3 (the default is 2650)
        """

        # Velocities of water, air, ice and rock (m/s)
        self.vw = vw
        self.va = va
        self.vr = vr
        self.vi = vi

        # Archie parameter
        self.m = m
        self.n = n
        self.phi = phi
        self.fr = 1 - self.phi  # fraction of rock
        
        # surface conductivity
        self.rhog = rhog
        # ~ self.cec = cec
        
        # additional new parameters (related to temperature)
        self.t = t
        self.alpha_t = alpha_t
        self.tf = tf
        self.tc = tc
        
        # temperature dependent parameters
        # ~ self.B = B
        # ~ self.l = l
        # ~ self.rhow = rhow
        self.B = self._compute_temp_dep(B, B_T0)
        self.l = self._compute_temp_dep(l, l_T0)
        self.rhow = self._compute_temp_dep(rhow, rhow_T0)
        self.R = self.l/self.B
        
        # ~ print("B", self.B)

    def _compute_temp_dep(self, val_T0, T0):
        return val_T0*(1 + self.alpha_t * (self.t - T0))
                
    def water_hauck2011(self, rho):
        fw = ((self.a * self.rhow * self.phi**self.n) /
              (rho * self.phi**self.m))**(1. / self.n)
        fw[np.isclose(fw, 0)] = 0
        return fw

    def water(self, rholo, rhohi):
        sigmahi = 1. / rhohi
        sigmalo = 1. / rholo
        
        mn = sigmahi - sigmalo

        fw = (self.rhow * self.phi**(self.n-self.m) * \
              (sigmahi - mn / self.R))**(1. / self.n)
        # ~ print("<<< fw before >>>", fw.min(), fw.max())
        # ~ fw[self.t<=self.tf] = fw[self.t<=self.tf] * np.exp(-((self.t[self.t<=self.tf]-self.tf)/self.tc)**2)
        # ~ print("<<< fw after >>>", fw.min(), fw.max())
              
        # ~ fw = (self.rhow * (sigmahi - mn/self.R))**(1/self.m)
        # ~ fw[np.isclose(fw, 0, atol=1e-2)] = 1e-2
        fw[np.isclose(fw, 0)] = 0
        return fw
        
    def water_coperey2019(self, rholo, rhohi):
        '''After Coperey et al. (2019)'''
        sigmahi = 1. / rhohi
        sigmalo = 1. / rholo
        
        mn = sigmahi - sigmalo

        fw = mn/(self.rhog*self.l*self.cec)
        fw[np.isclose(fw, 0)] = 0
        return fw

    def cec(self, rholo, rhohi):
        sigmahi = 1. / rhohi
        sigmalo = 1. / rholo
        
        mn = sigmahi - sigmalo
        
        # ~ print('<<< phi >>>')
        # ~ print(self.phi)
        cec = self.phi**(self.n - self.m) * mn / (self.water(rholo, rhohi)**(self.n-1) * self.rhog * self.l)
        # ~ cec = (self.water(rholo, rhohi)/self.phi)**(1-self.n)*self.phi**(1-self.m)*(mn/(self.rhog*self.l))
        # ~ cec = mn / (self.water(rholo, rhohi)**(self.m-1) * self.rhog * self.l)
        cec[np.isclose(cec, 0)] = 0
        # ~ print(cec)
        return cec

    def air_hauck2011(self, rho, v):
        fa = ((self.vi * self.va / (self.vi - self.va) * (
            1. / v - self.fr / self.vr - self.phi / self.vi - self.water_hauck2011(rho) *
            (1. / self.vw - 1. / self.vi))))
        fa[np.isclose(fa, 0)] = 0
        return fa
        
    def air(self, rholo, rhohi, v):
        fa = (self.vi * self.va) / (self.vi - self.va) * (
            1 / v - (1 - self.phi) / self.vr - self.phi / self.vi - self.water(rholo, rhohi) * 
            (1 / self.vw - 1 / self.vi))
        fa[np.isclose(fa, 0)] = 0
        return fa
        
    def ice_hauck2011(self, rho, v):
        fi = (self.va * self.vi) / (self.va - self.vi) * (
            1 / v - (1 - self.phi) / self.vr - self.phi / self.va - self.water_hauck2011(rho) *
            (1 / self.vw - 1 / self.va))
        fi[np.isclose(fi, 0)] = 0
        return fi        
        
    def ice(self, rholo, rhohi, v):
        fi = (self.va * self.vi) / (self.va - self.vi) * (
            1 / v - (1 - self.phi) / self.vr - self.phi / self.va - self.water(rholo, rhohi) * 
            (1 / self.vw - 1 / self.va))
        fi[np.isclose(fi, 0)] = 0
        return fi

    def rho_archie1942(self, fw, fa, fr=None):
        """Return electrical resistivity based on fraction of water `fw`."""
        if fr is None:
            phi = fw + fa
        else:
            phi = 1 - fr

        rho = self.a * self.rhow * phi**(-self.m) * (fw / phi)**(-self.n)
        if (rho <= 0).any():
            pg.warn(
                "Found negative resistivity, setting to nearest above zero.")
            rho[rho <= 0] = np.min(rho[rho > 0])
        return rho

    def rhohi(self, fw, fi, fa, cec, fr=None):
        """Return high frequency electrical resistivity based on 
        fraction of water `fw` and cation exchange capacity `CEC`."""

        rho = 1. / self.sigmahi(fw, fi, fa, cec, fr)
        if (rho <= 0).any():
            pg.warn(
                "Found negative resistivity, setting to nearest above zero.")
            rho[rho <= 0] = np.min(rho[rho > 0])
        return rho
    
    def sigmahi(self, fw, fi, fa, cec, fr):
        if fr is None:
            phi = fw + fa + fi
        else:
            phi = 1 - fr
        
        # ~ print("B", self.B)
        
        ret = (fw / phi)**self.n * phi**self.m * (1 / self.rhow) + \
               (fw / phi)**(self.n-1) * phi**(self.m-1) * self.rhog * self.B * cec
               
        # ~ print(ret.shape)
        
        return ret
        # ~ return fw**self.m * (1 / self.rhow) + \
               # ~ fw**(self.m-1) * self.rhog * self.B * cec
    
    def rholo(self, fw, fi, fa, cec, fr=None):
        """Return low frequency electrical resistivity based on 
        fraction of water `fw` and cation exchange capacity `CEC`."""

        rho = 1. / self.sigmalo(fw, fi, fa, cec, fr)
        if (rho <= 0).any():
            pg.warn(
                "Found negative resistivity, setting to nearest above zero.")
            rho[rho <= 0] = np.min(rho[rho > 0])
        return rho
    
    def sigmalo(self, fw, fi, fa, cec, fr):
        if fr is None:
            phi = fw + fa + fi
        else:
            phi = 1 - fr

        return (fw / phi)**self.n * phi**self.m * (1 / self.rhow) + \
               (fw / phi)**(self.n-1) * phi**(self.m-1) * self.rhog * (self.B - self.l) * cec
        # ~ return fw**self.m * (1 / self.rhow) + \
               # ~ fw**(self.m-1) * self.rhog * (self.B - self.l) * cec
    
    def rholo_deriv_fw(self, fw, fi, fa, cec, fr):
        return ((-1 / fw) * (self.n * (fw / (1-fr))**self.n * (1-fr)**self.m * (1/self.rhow) + \
            (self.n-1) * (fw / (1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * (self.B - self.l) * cec)) / \
            self.sigmalo(fw, fi, fa, cec, fr)**2
    
    # ~ def rholo_deriv_fw(self, fw, fi, fa, cec, fr):
        # ~ return (-((self.B - self.l) * cec * fw**(self.m-1) * self.rhog * (self.m-1))/(fw)-(fw**self.m * self.m * (1/self.rhow))/(fw)) / \
            # ~ self.sigmalo(fw, fi, fa, cec, fr)**2
    
    def rhohi_deriv_fw(self, fw, fi, fa, cec, fr):
        return ((-1 / fw) * (self.n * (fw / (1-fr))**self.n * (1-fr)**self.m * (1/self.rhow) + \
            (self.n-1) * (fw / (1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * self.B * cec)) / \
            self.sigmahi(fw, fi, fa, cec, fr)**2
    
    # ~ def rhohi_deriv_fw(self, fw, fi, fa, cec, fr):
        # ~ return (-(self.B * cec * fw**(self.m-1) * self.rhog * (self.m-1))/(fw)-(fw**self.m * self.m * (1/self.rhow))/(fw)) / \
            # ~ self.sigmahi(fw, fi, fa, cec, fr)**2

    # ~ def rholo_deriv_fr(self, fw, fi, fa, cec, fr):
        # ~ return 0
        
    def rholo_deriv_fr(self, fw, fi, fa, cec, fr):
        return ((-1/(1-fr)) * ((1-self.m) * (fw/(1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * (self.B - self.l) * cec + \
            (self.n-1)*(fw/(1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * (self.B - self.l) * cec - \
            self.m * (fw/(1-fr))**self.n * (1-fr)**self.m * (1/self.rhow)) + \
            self.n * (fw/(1-fr))**self.n * (1-fr)**self.m * (1/self.rhow)) / \
            self.sigmalo(fw, fi, fa, cec, fr)**2

    # ~ def rhohi_deriv_fr(self, fw, fi, fa, cec, fr):
        # ~ return 0
        
    def rhohi_deriv_fr(self, fw, fi, fa, cec, fr):
        return ((-1/(1-fr)) * (+(1-self.m) * (fw/(1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * self.B * cec + \
            (self.n-1)*(fw/(1-fr))**(self.n-1) * (1-fr)**(self.m-1) * self.rhog * self.B * cec - \
            self.m * (fw/(1-fr))**self.n * (1-fr)**self.m * (1/self.rhow)) + \
            self.n * (fw/(1-fr))**self.n * (1-fr)**self.m * (1/self.rhow)) / \
            self.sigmahi(fw, fi, fa, cec, fr)**2

    def rholo_deriv_fa(self, fw, fi, fa, cec, fr):
        return 0

    def rhohi_deriv_fa(self, fw, fi, fa, cec, fr):
        return 0

    def rholo_deriv_fi(self, fw, fi, fa, cec, fr):
        return 0

    def rhohi_deriv_fi(self, fw, fi, fa, cec, fr):
        return 0

    def rholo_deriv_t(self, fw, fi, fa, cec, fr):
        return 0

    def rhohi_deriv_t(self, fw, fi, fa, cec, fr):
        return 0
        
    def slowness(self, fw, fi, fa, fr=None):
        """Return slowness based on fraction of water `fw` and air `fa`."""
        if fr is None:
            fr = 1 - (fw + fa + fi)
        
        # ~ s = np.zeros_like(fw)
        
        # ~ if self.t > self.tf:
            # ~ mask = self.t > self.tf
            # ~ s[mask] = fw[mask] / self.vw + fr[mask] / self.vr + fa[mask] / self.va
        # ~ else:
            # ~ mask = self.t <= self.t
            # ~ s[mask] = fw[mask] / self.vw + fr[mask] / self.vr + fa[mask] / self.va + (fi[mask] / self.vi)*np.exp(-((self.t[mask]-self.tf)/self.tc))
        
        s = fw / self.vw + fi / self.vi + fr / self.vr + fa / self.va
        
        if (s <= 0).any():
            pg.warn("Found negative slowness, setting to nearest above zero.")
            s[s <= 0] = np.min(s[s > 0])
        return s
        
    def slowness_deriv_fw(self):
        return 1./self.vw
        
    def slowness_deriv_fa(self):
        return 1./self.va
        
    def slowness_deriv_fi(self):
        
        # ~ s_d_fi = np.zeros_like(self.t)
        
        # ~ s_d_fi[self.t > self.tf] = 1./self.vi
        # ~ mask = self.t <= self.tf
        # ~ s_d_fi[mask] = 1./self.vi*np.exp(-(self.t[mask]-self.tf)/self.tc)
        
        # ~ return s_d_fi
        
        return 1./self.vi
        
    def slowness_deriv_fr(self):
        return 1./self.vr
        
    def slowness_deriv_t(self):
        return 0

    def all(self, rholo, rhohi, v, mask=False):
        """Syntatic sugar for all fractions including a mask for unrealistic
        values."""

        # RVectors sometimes cause segfaults
        rholo = np.array(rholo)
        rhohi = np.array(rhohi)
        v = np.array(v)

        fw = self.water(rholo, rhohi)
        fa = self.air(rholo, rhohi, v)
        fi = self.ice(rholo, rhohi, v)
        cec = self.cec(rholo, rhohi)

        # Check that fractions are between 0 and 1
        array_mask = np.array(((fa < 0) | (fa > 1 - self.fr))
                              | ((fi < 0) | (fi > 1 - self.fr))
                              | ((fw < 0) | (fw > 1 - self.fr))
                              | ((self.fr < 0) | (self.fr > 1)))
        if array_mask.sum() > 1:
            print("WARNING: %d of %d fraction values are unphysical." % (int(
                array_mask.sum()), len(array_mask.ravel())))

        if mask:
            fa = np.ma.array(fa, mask=(fa < 0) | (fa > 1 - self.fr))
            fi = np.ma.array(fi, mask=(fi < 0) | (fi > 1 - self.fr))
            fw = np.ma.array(fw, mask=(fw < 0) | (fw > 1 - self.fr))

        return fa, fi, fw, cec, array_mask

    def show(self, mesh, rholo, rhohi, vel, mask=True, **kwargs):
        fa, fi, fw, mask = self.all(rholo, rhohi, vel, mask=mask)

        fig, axs = plt.subplots(4, 2, figsize=(16, 10))
        pg.show(mesh, fw, ax=axs[0, 0], label="Water content", hold=True,
                logScale=False, cMap="Blues", **kwargs)
        pg.show(mesh, fi, ax=axs[1, 0], label="Ice content", hold=True,
                logScale=False, cMap="Greens", **kwargs)
        pg.show(mesh, fa, ax=axs[2, 0], label="Air content", hold=True,
                logScale=False, cMap="Greens", **kwargs)
        pg.show(mesh, rholo, ax=axs[0, 1], label="Rho_0", hold=True,
                cMap="Spectral_r", logScale=True, **kwargs)
        pg.show(mesh, rhohi, ax=axs[1, 1], label="Rho_inf", hold=True,
                cMap="Spectral_r", logScale=True, **kwargs)
        pg.show(mesh, vel, ax=axs[2, 1], label="Velocity", logScale=False,
                hold=True, **kwargs)
        pg.show(mesh, self.phi, ax=axs[2, 1], label="Porosity", logScale=False,
                cMap='Oranges', hold=True, **kwargs)
        pg.show(mesh, 1 - self.phi, ax=axs[2, 0], label="Rock content", logScale=False,
                cMap='Oranges', hold=True, **kwargs)
        for ax in axs.flat:
            ax.set_facecolor("0.5")
        return fig, axs

def testPetroMod():
    # Based on parameters proposed by Hauck et al. (2011)
    pm = PetroMod(vw=1500, va=300, vi=3500, vr=6000, phi=0.5, 
                 n=2., m=2.,
                 rhow=150., rhow_T0=21,
                 rhog=2650, cec=18.8, 
                 t=-4, 
                 alpha_t=.02, tf=-1, tc=4,
                 R=.1, 
                 B=3.1*1e-9, B_T0=25,
                 l=3.0*1e-10, l_T0=25)
    # ~ assert pm.water(10.0) == 10.0
    v = np.linspace(500, 6000, 1000)
    rholo = np.logspace(2, 7, 1000)
    rhohi = np.logspace(2.01, 7.01, 1000)
    vvec, rholovec = np.meshgrid(v, rholo)
    vvec, rhohivec = np.meshgrid(v, rhohi)

    fa, fi, fw, cec, mask = pm.all(rholovec, rhohivec, vvec, mask=True)

    cmap = plt.cm.get_cmap('Spectral_r', 41)
    fig, axs = plt.subplots(3, figsize=(6, 4.5), sharex=True)
    labels = ["Air content", "Ice content", "Water content"]
    for data, ax, label in zip([fa, fi, fw], axs, labels):
        im = ax.imshow(
            data[::-1], cmap=cmap, extent=[
                v.min(),
                v.max(),
                np.log10(rholo.min()),
                np.log10(rholo.max())
            ], aspect="auto", vmin=0, vmax=0.5)
        plt.colorbar(im, ax=ax, label=label)

    axs[1].set_ylabel(r"Log resistivity ($\Omega$m)")
    axs[-1].set_xlabel("Velocity (m/s)")

    fig.tight_layout()

    plt.figure()
    im = plt.imshow(fa + fi + fw, vmin=0, vmax=0.5)
    plt.colorbar(im)

    return fig


if __name__ == '__main__':
    import seaborn
    seaborn.set(font="Fira Sans", style="ticks")
    plt.rcParams["image.cmap"] = "viridis"
    fig = testPetroMod()
    fig.savefig("4PM_value_range.pdf")
