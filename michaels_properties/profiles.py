from tangos.properties import HaloProperties, TimeChunkedProperty
from tangos.properties.spherical_region import SphericalRegionHaloProperties
import numpy as np
import pynbody
from tangos.properties.centring import centred_calculation

kb = pynbody.array.SimArray(1.380658e-16, 'erg K**-1')
mh = pynbody.array.SimArray(1.6726219e-24, 'g')

def emissivity(rho, T, mu, tcool):
    return 3. / 2. * rho * kb * T / (mu * mh * tcool)


def tcool(rho, T, mu):
    # taken from https://arxiv.org/abs/astro-ph/9809159 eq 12
    fm = 1.0  # metalicity dependent factor, 1.0 for solar, 0.03 for pristine
    C1 = 3.88e11
    C2 = 5e7
    return C1 * mu * mh * T ** 0.5 / (rho * (1 + C2 * fm / T))


@pynbody.analysis.profile.Profile.profile_property
def Tew(self):
    temp = np.zeros(self.nbins)
    for i in range(self.nbins):
        subs = self.sim[self.binind[i]]
        #use = np.where(subs.g['temp'] > temp_cut)[0]
        mu = 0.58
        tc = tcool(subs.g['rho'].in_units('g cm**-3'), subs.g['temp'], mu)
        em = emissivity(subs.g['rho'].in_units('g cm**-3'), subs.g['temp'], mu, tc)
        temp[i] = np.sum(em * subs.g['temp']) / np.sum(em)

    return kb.in_units('keV K**-1') * temp


@pynbody.analysis.profile.Profile.profile_property
def Tmw(self):
    temp = np.zeros(self.nbins)
    for i in range(self.nbins):
        subs = self.sim[self.binind[i]]
        #use = np.where(subs.g['temp'] > temp_cut)[0]
        temp[i] = np.sum(subs.g['mass'] * subs.g['temp']) / np.sum(subs.g['mass'])
    return kb.in_units('keV K**-1') * temp

@pynbody.analysis.profile.Profile.profile_property
def rho_e(self):
    n_e = np.zeros(self.nbins)
    for i in range(self.nbins):
        subs = self.sim[self.binind[i]]
        #use = np.where(subs.g['temp'] > temp_cut)[0]
        n_e[i] = np.sum(subs.g['ne'] * subs.g['mass'].in_units('m_p'))/self._binsize.in_units('cm**'+str(int(self.ndim)))[i]
    return n_e


class GasProfiles(HaloProperties):
    _temp_cut = 1.26e6
    @classmethod
    def name(self):
        return "Tew_profile", "Tmw_profile", "rho_e_profile"

    def plot_x0(cls):
        return 0.05

    @classmethod
    def plot_xdelta(cls):
        return 0.1

    @classmethod
    def plot_xlabel(cls):
        return "R/kpc"

    @classmethod
    def plot_ylog(cls):
        return False

    @classmethod
    def plot_xlog(cls):
        return False

    @staticmethod
    def plot_ylabel():
        return "T$_{ew}$ keV", "T$_{mw}$ keV"

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    @centred_calculation
    def calculate(self, halo, existing_properties):
        #halo['pos'] -= existing_properties['SSC']
        #halo.wrap()
        delta = self.plot_xdelta()
        nbins = int(existing_properties['max_radius']/ delta)
        maxrad = delta * (nbins + 1)
        ps = pynbody.analysis.profile.Profile(halo.g[pynbody.filt.HighPass('temp',self._temp_cut)], type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)
        Tew = ps['Tew']
        Tmw = ps['Tmw']
        rho_e = ps['rho_e']
        return Tew, Tmw, rho_e


class OldHaloDensityProfile(SphericalRegionHaloProperties):

    @classmethod
    def name(self):
        return "dm_density_profile", "dm_mass_profile", "tot_density_profile", "tot_mass_profile", "gas_density_profile", "gas_mass_profile", "star_density_profile", "star_mass_profile"

    @classmethod
    def plot_x0(cls):
        return 0.05

    @classmethod
    def plot_xdelta(cls):
        return 0.1

    @classmethod
    def plot_xlabel(cls):
        return "r/kpc"

    @staticmethod
    def plot_ylabel():
        return r"$\rho/M_{\odot}\,kpc^{-3}$", r"$M/M_{\odot}$", r"$\rho/M_{\odot}\,kpc^{-3}$", r"$M/M_{\odot}$", r"$\rho/M_{\odot}\,kpc^{-3}$", r"$M/M_{\odot}$", r"$\rho/M_{\odot}\,kpc^{-3}$", r"$M/M_{\odot}$"

    def rstat(self, halo, maxrad, cen, delta=0.1):
        mass_a = []
        rho_a = []

        mass_x = 0

        V_x = 0
        halo['pos'] -= cen
        halo.wrap()

        nbins = int(maxrad / delta)
        maxrad = delta * (nbins + 1)

        pro = pynbody.analysis.profile.Profile(halo, type='lin', ndim=3,
                                               min=0, max=maxrad, nbins=nbins)

        rho_a = pro['density']
        mass_a = pro['mass_enc']

        halo['pos'] += cen
        halo.wrap()
        rho_a = np.array(rho_a)
        mass_a = np.array(mass_a)

        return rho_a, mass_a

    def calculate(self, halo, existing_properties):

        halo.dm['mass']
        try:
            halo.gas['mass']
            halo.star['mass']
        except:
            pass
        delta = existing_properties.get('delta',0.1)
        self.mark_timer('dm')
        dm_a, dm_b = self.rstat(
            halo.dm, existing_properties["Rvir"], existing_properties["SSC"],delta)
        self.mark_timer('tot')
        tot_a, tot_b = self.rstat(
            halo, existing_properties["Rvir"], existing_properties["SSC"],delta)
        self.mark_timer('gas')
        gas_a, gas_b = self.rstat(
            halo.gas, existing_properties["Rvir"], existing_properties["SSC"],delta)
        self.mark_timer('star')
        star_a, star_b = self.rstat(
            halo.star, existing_properties["Rvir"], existing_properties["SSC"],delta)
        return dm_a, dm_b, tot_a, tot_b, gas_a, gas_b, star_a, star_b

    def requires_property(self):
        return ["SSC", "Rvir"]