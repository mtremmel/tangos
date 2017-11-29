
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

class BHGalHistogram(LiveHaloProperties,TimeChunkedProperty):

    def __init__(self, simulation=None, property='BH_mdot_histogram', operation='max', bhtype='BH_central'):
        super(BHGalHistogram, self).__init__(simulation)
        self._operation=operation
        self._property = property
        self._bhtype = bhtype

    @classmethod
    def name(cls):
        return "bh_histogram_galaxy"

    def requires_property(self):
        return []

    def live_calculate(self, halo, *args):
        if halo.object_typecode != 0:
            return None

        if self._bhtype not in list(halo.keys()):
            mdot = np.zeros(self.nbins)[self.store_slice(halo.timestep.time_gyr)]
            return mdot

        if type(halo[self._bhtype]) is list:
            all_hists = []
            for bh in halo[self._bhtype]:
                if self._property not in list(bh.keys()):
                    mdot_part = np.zeros(self.nbins)[self.store_slice(halo.timestep.time_gyr)]
                else:
                    mdot_part = bh.calculate('raw('+self._property+')')
                all_hists.append(mdot_part)
            if len(all_hists) != len(halo[self._bhtype]):
                raise RuntimeError, "bad size! "+str(halo)
            if self._operation=='sum':
                print "in sum", self._operation
                mdot =  np.sum(all_hists,axis=0)
            if self._operation=='max':
                print "in max", self._operation
                mdot =  np.max(all_hists,axis=0)
        else:
            bh = halo[self._bhtype]
            if self._property not in list(bh.keys()):
                mdot = np.zeros(self.nbins)[self.store_slice(halo.timestep.time_gyr)]
            else:
                mdot = bh.calculate('raw('+self._property+')')
        return mdot
