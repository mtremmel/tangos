#from __future__ import absolute_import

from tangos.properties import HaloProperties, TimeChunkedProperty
import numpy as np
import pynbody

class InnerStarFormHistogram(TimeChunkedProperty):

    _maxr_frac = 0.1
    _ncrit = 200

    requires_particle_data = True
    names = "inner_SFR_histogram"

    def requires_property(self):
        return ["shrink_center","tot_mass_profile"]

    def radius(self, halo, existing_properties):
        a = halo.properties['a']
        H_z = pynbody.analysis.cosmology._a_dot(a, halo.properties['h'], halo.properties['omegaM0'], halo.properties['omegaL0']) / a
        H_z = pynbody.units.Unit("100 km s^-1 Mpc^-1") * H_z
        rho_crit = (3 * H_z ** 2) / (8 * np.pi * pynbody.units.G)
        rho_mean = existing_properties['tot_mass_profile']/(4./3. * np.pi * ((np.arange(len(existing_properties['tot_mass_profile']))+1)*0.1)**3)
        return np.where(rho_mean>rho_crit.in_units('Msol kpc**-3')*self._ncrit)[0][-1]*0.1

    def calculate(self, halo, existing_properties):
        filter = pynbody.filt.Sphere(self._maxr_frac*self.radius(halo,existing_properties),cen=existing_properties['shrink_center'])
        M,_ = np.histogram(halo.st[filter]['tform'].in_units("Gyr"),weights=halo.st[filter]['massform'].in_units("Msol"),bins=self.nbins,range=(0,self.tmax_Gyr))
        t_now = halo.properties['time'].in_units("Gyr")
        M/=self.delta_t
        M = M[self.store_slice(t_now)]

        return M

    @classmethod
    def reassemble(cls, *options):
        reassembled = TimeChunkedProperty.reassemble(*options)
        return reassembled/1e9 # Msol per Gyr -> Msol per yr
