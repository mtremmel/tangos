
from tangos.properties import HaloProperties, TimeChunkedProperty
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
        t_grid = np.linspace(0, self.tmax_Gyr+dt, self.nbins)
        mdot_grid_n, = np.histogram(t_orbit[order],bins=t_grid)
        mdot_grid_sum, = np.histogram(t_orbit[order],bins=t_grid,weights=Mdot_orbit[order])
        mdot_grid_ave = mdot_grid_sum/mdot_grid_n.astype(np.float)


        #Mdot_grid = scipy.interpolate.interp1d(t_orbit[order], Mdot_orbit[order], bounds_error=False)(t_grid)


        #print t_max
        #print Mdot_grid

        return mdot_grid_ave[self.store_slice(t_max)]