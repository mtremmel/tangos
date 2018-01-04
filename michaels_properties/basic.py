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
		return ['tot_mass_profile']

	def live_calculate(self, halo,*args):
		ts = halo.timestep
		z = ts.redshift
		a = 1.0 / (1.0 + z)
		H_z = pynbody.analysis.cosmology._a_dot(a, self._h0, self._omegaM0, self._omegaL0) / a
		H_z = pynbody.units.Unit("100 km s^-1 Mpc^-1") * H_z
		rho_crit = (3 * H_z ** 2) / (8 * np.pi * pynbody.units.G)
		rho_mean = halo['tot_mass_profile']/(4./3. * np.pi * ((np.arange(len(halo['tot_mass_profile']))+1)*0.1)**3)
		return np.where(rho_mean>rho_crit.in_units('Msol kpc**-3')*self._ncrit)[0][-1]*0.1

class EscapeEnergy(LiveHaloProperties):

	names='escape_energy_profile'

	def requires_property(self):
		return ['tot_mass_profile']

	def plot_x0(cls):
		return 0.05

	def plot_xdelta(cls):
		return 0.1

	def live_calculate(self,halo):
		mass = np.asarray(halo['tot_mass_profile'])

		max_x = len(mass) * 0.1
		pot_dx = 0.01

		r = np.arange(0.05, 0.04 + 0.1 * len(mass), 0.1)
		rnew = np.arange(pot_dx, max_x, pot_dx)

		min_r_i = int(0.35 / pot_dx)

		mass_new = np.ones(len(rnew))

		import scipy
		mass_new[min_r_i:] = scipy.interp(rnew[min_r_i:], r, mass)
		mass_new[:min_r_i] = mass_new[min_r_i] * (rnew[:min_r_i] / rnew[min_r_i]) ** 3

		force = mass_new / (rnew ** 2) * pynbody.units.G.in_units('cm**2 s**-2 Msol**-1 kpc')


		pot = np.cumsum(pot_dx * force)
		radius = halo.calculate('radius()')
		pot_diff = pot[int(radius/0.01)] - pot[::10]
		#in cgs units
		return pot_diff

