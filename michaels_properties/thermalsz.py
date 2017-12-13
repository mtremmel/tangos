import pynbody
from tangos.properties import HaloProperties


@pynbody.snapshot.SimSnap.derived_quantity
def electron_density(sim):
    sim['HII'] = pynbody.snapshot.tipsy.HII(sim)
    sim['HeIII'] = pynbody.snapshot.tipsy.HeIII(sim)
    number_density = (sim['rho']).in_units("m_p cm^-3")*pynbody.snapshot.tipsy.ne(sim)#sim['ne']
    number_density.units="cm^-3"
    number_density[number_density<0]=0 # negative ne due to missing metals
    return number_density

@pynbody.snapshot.SimSnap.derived_quantity
def y_integrand(sim):
    ne = sim['electron_density']
    kT_by_me = (sim['temp']*pynbody.units.k/(pynbody.units.m_e*pynbody.units.c**2)).in_units("1")
    thomson = pynbody.units.Unit("6.652e-29 m**2")
    return (thomson*ne*kT_by_me).in_units("kpc^-1")


class ComptonImage(HaloProperties):
    @classmethod
    def name(self):
        return "compton_y_map"

    @classmethod
    def plot_extent(self):
        return 500.0

    def requires_property(self):
        return ["shrink_center"]

    def region_specification(self, db_data):
        import pynbody
        return pynbody.filt.Sphere("1 Mpc", db_data['shrink_center'])

    def calculate(self, particle_data, properties):
        particle_data['pos']-=properties['shrink_center']
        try:
            size = self.plot_extent()
            g = self._render_projected(particle_data.gas, size)
        finally:
            particle_data['pos']+=properties['shrink_center']
        return g

    def _render_projected(self, f, size):
        import pynbody.plot
        im = pynbody.plot.sph.image(f[pynbody.filt.BandPass(
            'z', -size / 2, size / 2)], 'y_integrand', size, units="1",
            noplot=True, resolution=1000)
        return im
