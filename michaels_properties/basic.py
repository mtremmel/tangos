from . import LiveHaloProperties
import numpy as np
import pynbody

class Radius(LiveHaloProperties):
	def __init__(self, simulation, n_crit=200):
		super(Radius, self).__init__(simulation)
		self._omegaM0 = 0.3086
		self._omegaL0 = 0.6914
		self._h0 = 0.67
		self._ncrit = n_crit

	@classmethod
	def name(cls):
		return 'radius'

	def requires_property(self):
		return['tot_mass_profile']

	def live_calculate(self, halo,*args):
		ts = halo.timestep
		z = ts.redshift
		a = 1.0 / (1.0 + z)
		H_z = pynbody.analysis.cosmology._a_dot(a, self._h0, self._omegaM0, self._omegaL0) / a
		H_z = pynbody.units.Unit("100 km s^-1 Mpc^-1") * H_z
		rho_crit = (3 * H_z ** 2) / (8 * np.pi * pynbody.units.G)
		rho_mean = halo['tot_mass_profile']/(4./3. * np.pi * ((np.arange(len(halo['tot_mass_profile']))+1)*0.1)**3)
		return np.where(rho_mean>rho_crit.in_units('Msol kpc**-3')*self._ncrit)[0][-1]*0.1
