
from tangos.properties import HaloProperties, TimeChunkedProperty, LiveHaloProperties
import numpy as np
import pynbody
from tangos_nbodyshop_properties.BH import BHShortenedLog
import scipy

class BHAccAveHistogram(TimeChunkedProperty):
    @classmethod
    def name(self):
        return "BH_mdot_histogram_ave"

    def requires_property(self):
        return []


    def preloop(self, f, db_timestep):
        self.log = BHShortenedLog.get_existing_or_new(db_timestep.filename)

    @classmethod
    def no_proxies(self):
        return True

    def calculate(self, halo, properties):

        halo = halo.s

        if len(halo)!=1:
            raise RuntimeError("Not a BH!")

        if halo['tform'][0]>0:
            raise RuntimeError("Not a BH!")

        mask = self.log.vars['bhid']==halo['iord']
        if(mask.sum()==0):
            raise RuntimeError("Can't find BH in .orbit file")

        t_orbit = self.log.vars['time'][mask]
        Mdot_orbit = self.log.vars['mdotmean'][mask]
        order = np.argsort(t_orbit)

        t_max = properties.timestep.time_gyr
        dt = self.tmax_Gyr/self.nbins
        t_grid = np.arange(0, self.tmax_Gyr+dt, dt)
        #print t_grid
        mdot_grid_n, _ = np.histogram(t_orbit[order],bins=t_grid)
        mdot_grid_sum, _ = np.histogram(t_orbit[order],bins=t_grid,weights=Mdot_orbit[order])
        mdot_grid_ave = mdot_grid_sum/mdot_grid_n.astype(np.float)


        #Mdot_grid = scipy.interpolate.interp1d(t_orbit[order], Mdot_orbit[order], bounds_error=False)(t_grid)


        #print t_max
        #print Mdot_grid

        return mdot_grid_ave[self.store_slice(t_max)]

class BHGalAccHistogram(LiveHaloProperties,TimeChunkedProperty):

    def __init__(self, simulation=None, choose='BH_mdot_ave', minmax='max', bhtype='BH_central'):
        super(BHGalAccHistogram, self).__init__(simulation)
        self._maxmin = minmax
        self._choicep = choose
        self._bhtype = bhtype

    @classmethod
    def name(cls):
        return "bh_mdot_histogram_galaxy"

    def requires_property(self):
        return []

    def live_calculate(self, halo, *args):
        try:
            bh = halo.calculate('bh('+self._choicep+', '+self._maxmin+','+self._bhtype+')')
        except:
            return None
        mdot = bh.calculate('raw(BH_mdot_histogram)')
        return mdot