from tangos.properties import HaloProperties, TimeChunkedProperty
from tangos.properties.pynbody.spherical_region import SphericalRegionHaloProperties
import numpy as np
import pynbody
from tangos.properties.pynbody.centring import centred_calculation

class NewFlowProfile(SphericalRegionHaloProperties):
    #_xmax = 100.0
    _threshold_vel = 20.0

    def region_specification(self, db_data):
        return pynbody.filt.Sphere(db_data['max_radius'], db_data['shrink_center']) & \
               (pynbody.filt.FamilyFilter(pynbody.family.gas)|pynbody.filt.FamilyFilter(pynbody.family.star))

    @classmethod
    def name(cls):
        return "gas_inflow_Mdot", "gas_outflow_Mdot", \
               "gas_inflow_vel", "gas_outflow_vel", \
               "gas_inflow_vel2", "gas_outflow_vel2", \
               "gas_inflow_temp", "gas_outflow_temp"
#               "inflow_Mdot_dm", "outflow_Mdot_dm", \
#               "inflow_vel_dm", "outflow_vel_dm", \
#               "inflow_vel2_dm", "outflow_vel2_dm"


    def plot_x0(cls):
        return 0.5

    @classmethod
    def plot_xdelta(cls):
        return 1.0

    def profile_calculation(self, f_gas, vr_cut, rvir):
        f_gas['Mdot'] = f_gas['mass'] * f_gas['vr_mean'] / (pynbody.units.Unit("kpc")*self.plot_xdelta())
        if vr_cut<0:
            f_gas['Mdot']*=(f_gas['vr_mean']<vr_cut).view(np.ndarray)
            f_gas['Mdot'] *= -1
        else:
            f_gas['Mdot']*=(f_gas['vr_mean']>vr_cut).view(np.ndarray)



        pro = pynbody.analysis.profile.Profile(f_gas, min=0.0,max=rvir,
                                               nbins=int(rvir/self.plot_xdelta()),
                                               ndim=3, weight_by='Mdot')
        return pro['weight_fn'].in_units("Msol yr^-1"), pro['vr_mean'], pro['vr_mean2'], pro['temp']

    #@centred_calculation
    def calculate(self, halo, properties):
        original_positions = np.array(halo.ancestor['pos'])
        halo.ancestor['pos'] -= properties['shrink_center']
        halo.ancestor.wrap()
        with pynbody.analysis.halo.vel_center(halo):
            halo.gas['vr_mean'] = np.sum(halo.gas['v_mean']*halo.gas['pos'],axis=1)/halo.g['r']
            halo.gas['vr_mean2'] = halo.gas['vr_mean'] ** 2
            inflow_Mdot, inflow_vel, inflow_vel2, inflow_temp = self.profile_calculation(halo.ancestor.gas, -self._threshold_vel, properties['max_radius'])
            outflow_Mdot, outflow_vel, outflow_vel2, outflow_temp = self.profile_calculation(halo.ancestor.gas, self._threshold_vel,properties['max_radius'])
            #inflow_Mdot_dm, inflow_vel_dm, inflow_vel2_dm, _ = self.profile_calculation(halo.ancestor.dm, -self._threshold_vel,properties['max_radius'])
            #outflow_Mdot_dm, outflow_vel_dm, outflow_vel2_dm, _ = self.profile_calculation(halo.ancestor.dm, self._threshold_vel,properties['max_radius'])
        halo.ancestor['pos'] = original_positions
        return inflow_Mdot, outflow_Mdot, -inflow_vel, outflow_vel,  inflow_vel2, outflow_vel2, inflow_temp, outflow_temp


def get_outflow_particles(fgas):
    threshold_vel = 20.0
    med_v = np.median(fgas['vr'])
    mad_v = np.mean(np.abs(fgas['vr']-np.median(fgas['vr'])))
    if med_v > mad_v:
        return (fgas['vr']>threshold_vel).view(np.ndarray)
    else:
        return (fgas['vr']>max(threshold_vel,med_v+mad_v)).view(np.ndarray)


@pynbody.analysis.profile.Profile.profile_property
def outflow(self,data='mdot'):
    out = np.zeros(self.nbins)
    dr = self['bin_edges'][1:] - self['bin_edges'][:-1]
    for i in range(self.nbins):
        subs = self.sim[self.binind[i]]
        f_gas = subs.g
        mdot = f_gas['mass'] * f_gas['vr'] / dr[i]
        mdot *= get_outflow_particles(f_gas)
        if data=='mdot':
            out[i] = np.sum(mdot.in_units("Msol yr^-1"))
        else:
            out[i] = np.sum(mdot*f_gas[data])/np.sum(mdot)
    return out


class WindProfile(SphericalRegionHaloProperties):

    names = "winds_mdot", "winds_v", "winds_v2","winds_temp",
    def region_specification(self, db_data):
        return pynbody.filt.Sphere(db_data['max_radius'], db_data['shrink_center']) & \
               (pynbody.filt.FamilyFilter(pynbody.family.gas)|pynbody.filt.FamilyFilter(pynbody.family.star))

    def plot_x0(cls):
        return 0.5

    @classmethod
    def plot_xdelta(cls):
        return 1.0

    def profile_calculation(self, f_gas, rvir):
        pro = pynbody.analysis.profile.Profile(f_gas, min=0.0,max=rvir,
                                               nbins=int(rvir/self.plot_xdelta()), ndim=3)
        return pro['outflow, mdot'], pro['outflow, vr'], pro['outflow,vr2'], pro['outflow,temp']

    def calculate(self, halo, properties):
        original_positions = np.array(halo.ancestor['pos'])
        original_velocities = np.array(halo.ancestor['vel'])
        halo.ancestor['pos'] -= properties['shrink_center']
        halo.ancestor.wrap()
        with pynbody.analysis.halo.vel_center(halo):
            halo.gas['vr2'] = halo.gas['vr'] ** 2
            outflow_Mdot, outflow_vel, outflow_vel2, outflow_temp = self.profile_calculation(halo.ancestor.gas,properties['max_radius'])
        halo.ancestor['pos'] = original_positions
        halo.ancestor['vel'] = original_velocities
        return outflow_Mdot, outflow_vel, outflow_vel2, outflow_temp
