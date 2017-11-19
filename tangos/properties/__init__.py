from __future__ import absolute_import
import numpy as np
from tangos.util import timing_monitor
import six
from six.moves import zip
import importlib
import warnings
from .. import live_calculation

class HaloPropertiesMetaClass(type):
    # Present to register new subclasses of HaloProperties, so that subclasses can be dynamically
    # instantiated when required based on their cls.name() values. Stored as a dictionary so that
    # reloaded classes overwrite their old versions.
    def __init__(cls, name, bases, dict):
        type.__init__(cls, name, bases, dict)
        cls._all_classes[name] = cls

class HaloProperties(six.with_metaclass(HaloPropertiesMetaClass,object)):
    _all_classes = {}

    @classmethod
    def all_classes(cls):
        return list(cls._all_classes.values())

    def __init__(self, simulation):
        """Initialise a HaloProperties calculation object

        :param simulation: The simulation from which the properties will be derived
        :type simulation: tangos.core.simulation.Simulation
        """
        self._simulation = simulation
        self.timing_monitor = timing_monitor.TimingMonitor()

    def requires_array(self):
        """Returns a list of loaded arrays required to
        calculate this property"""
        return []

    @classmethod
    def requires_particle_data(self):
        """If this returns false, the class can do its
        calculation without any raw simulation data loaded
        (i.e. derived from other properties)"""
        return True

    @classmethod
    def name(self):
        """Returns either the name or a list of names of
        properties that will be calculated by this class"""
        return "undefined"

    @classmethod
    def index_of_name(cls, name):
        """Returns the index of the named property in the
        results returned from calculate().

        For example, for a BasicHaloProperties object X,
        X.calculate(..)[X.index_of_name("SSC")] returns the SSC.
        """
        name = name.split("(")[0]
        return cls.name().index(name)

    @classmethod
    def no_proxies(self):
        """Returns True if the properties MUST be supplied
        as an actual Halo object rather than the normal
        dictionary-like proxy that is supplied. Only use if
        absolutely necessary; has adverse consequences e.g.
        because uncommitted updates are not reflected into
        your calculate function"""
        return False

    def requires_property(self):
        """Returns a list of existing properties
        required to calculate this property"""
        return []

    def preloop(self, sim, db_timestep):
        """Perform one-per-snapshot calculations, given the loaded simulation data and TimeStep object"""
        pass

    def region_specification(self, db_data):
        """Returns an abstract specification of the region that this halo property is to be calculated on,
        or None if we want the halo particles as defined by the finder.

        See spherical_region.SphericalRegionHaloProperties for an example and useful base class for returning
        everything within the virial radius."""
        return None


    def mark_timer(self, label=None):
        """Called by subfunctions to mark a time"""
        self.timing_monitor.mark(label)


    def accept(self, db_entry):
        for x in self.requires_property():
            if db_entry.get(x, None) is None:
                return False
        return True

    def calculate(self, pynbody_halo_data, halo_entry):
        """Calculate the properties using the given data

        :param pynbody_halo_data: The halo data, if available
        :type pynbody_halo_data: pynbody.snapshot.SimSnap (when the pynbody backend is in use, otherwise could be a yt snapshot etc)

        :param halo_entry: The database object associated with the halo, if available
        :type halo_entry: tangos.core.halo.Halo
        :return: All properties as named by names()
        """
        raise NotImplementedError

    def live_calculate(self, halo_entry, *input_values):
        """Calculate the result of a function, using the existing data in the database alone

        :param halo_entry: The database object associated with the halo
        :type halo_entry: tangos.core.halo.Halo

        :param input_values: Input values for the function
        :return: All function values as named by self.name()
        """
        if self.requires_particle_data():
            raise(RuntimeError("Cannot live-calculate a property that requires particle data"))
        return self.calculate(None, halo_entry)

    def live_calculate_named(self, name, halo_entry, *input_values):
        """Calculate the result of a function, using the existing data in the database alone

        :param name: The name of the one property to return (which must be one of the values specified by self.name())
        :param halo_entry: The database object associated with the halo
        :type halo_entry: tangos.core.halo.Halo

        :param input_values: Input values for the function
        :return: The single named value
        """
        values = self.live_calculate(halo_entry, *input_values)
        names = self.name()
        if isinstance(names, six.string_types):
            return values
        else:
            return values[self.name().index(name)]

    def calculate_from_db(self, db):
        if self.requires_particle_data():
            region_spec =  self.region_specification(self)
            if region_spec:
                halo_particles = db.timestep.load_region(region_spec)
            else:
                halo_particles = db.load()

            preloops_done = getattr(halo_particles.ancestor, "_did_preloop", [])
            halo_particles.ancestor._did_preloop = preloops_done

            if str(self.__class__) not in preloops_done:
                self.preloop(halo_particles.ancestor, db.timestep)
                preloops_done.append(str(self.__class__))
        else:
            halo_particles = None
        return self.calculate(halo_particles, db)

    def plot_x_values(self, for_data):
        """Return a suitable array of x values to match the
        given y values"""
        return np.arange(self.plot_x0(), self.plot_x0() + self.plot_xdelta() * (len(for_data) - 0.5), self.plot_xdelta())

    def plot_x_extent(self):
        return None

    def plot_extent(self):
        return None

    def plot_x0(self):
        return 0

    def plot_xdelta(self):
        return 1.0

    def plot_xlabel(self):
        return None

    def plot_ylabel(self):
        return None

    def plot_yrange(self):
        return None

    def plot_xlog(self):
        return True

    def plot_ylog(self):
        return True

    def plot_clabel(self):
        return None

class TimeChunkedProperty(HaloProperties):
    """TimeChunkedProperty implements a special type of halo property where chunks of a histogram are stored
    at each time step, then appropriately reassembled when the histogram is retrieved."""

    nbins = 2000
    tmax_Gyr = 20.0
    minimum_store_Gyr = 1.0

    @property
    def delta_t(self):
        return self.tmax_Gyr/self.nbins

    @classmethod
    def bin_index(self, time):
        """Convert a time (Gyr) to a bin index in the histogram"""
        index = int(self.nbins*time/self.tmax_Gyr)
        if index<0:
            index = 0
        return index

    @classmethod
    def store_slice(self, time):
        """Tells subclasses which have generated a histogram over all time which slice of that histogram
        they should store."""
        return slice(self.bin_index(time-self.minimum_store_Gyr), self.bin_index(time))

    @classmethod
    def reassemble(cls, property, halo, reassembly_type='major'):
        """Reassemble a histogram by suitable treatment of the merger tree leading up to the current halo.

        This function is normally called by the framework (see Halo.get_data_with_reassembly_options) and you would
        rarely call it directly yourself. From within a live-calculation, it can be accessed using the
        reassemble(halo_property, options...) function. See live_calculation.md for more information.

        :param: property - the halo property for which the reassembly should occur

        :param: reassembly_type - if 'major' (default), return the histogram for the major progenitor branch
                                - if 'sum', return the histogram summed over all progenitors
                                (e.g. this can be used to return SFR histograms that count infalling material as well
                                as the major progenitor)
                                - if 'place', return only the histogram stored at this step but place it within
                                a correctly zero-padded array
                                - if 'raw', return the raw data
        """

        print 'here in reassemble', property, type(property)
        from tangos import relation_finding as rfs
        if reassembly_type=='major':
            return cls._reassemble_using_finding_strategy(property, halo, strategy = rfs.MultiHopMajorProgenitorsStrategy)
        elif reassembly_type=='major_across_simulations':
            return cls._reassemble_using_finding_strategy(property, halo, strategy = rfs.MultiHopMajorProgenitorsStrategy,
                                                          strategy_kwargs = {'target': None})
        elif reassembly_type=='sum':
            return cls._reassemble_using_finding_strategy(property, halo, strategy = rfs.MultiHopAllProgenitorsStrategy)
        elif reassembly_type=='place':
            return cls._place_data(property.halo.timestep.time_gyr, property.data_raw)
        elif reassembly_type=='raw':
            return property.data_raw
        else:
            raise ValueError("Unknown reassembly type")

    @classmethod
    def _place_data(cls, time, raw_data):
        final = np.zeros(cls.bin_index(time))
        end = len(final)
        start = end - len(raw_data)
        final[start:] = raw_data
        return final

    @classmethod
    def _reassemble_using_finding_strategy(cls, property, halo, strategy, strategy_kwargs={}):
        #print "in reassemble finding strat", property, isinstance(property, live_calculation.Calculation)
        if not isinstance(property, live_calculation.Calculation):
            name = property.name.text
            t, stack = halo.calculate_for_descendants("t()", "raw(" + name + ")", strategy=strategy,strategy_kwargs=strategy_kwargs)
        else:
            name = property.name()
            t, stack = halo.calculate_for_descendants("t()", "raw(" + name + "())", strategy=strategy, strategy_kwargs=strategy_kwargs)
        final = np.zeros(cls.bin_index(t[0]))
        previous_time = -1
        for t_i, hist_i in zip(t, stack):
            end = cls.bin_index(t_i)
            start = end - len(hist_i)
            valid = hist_i == hist_i
            if t_i != previous_time:
                # new timestep; overwrite what was there previously
                final[start:end][valid] = hist_i[valid]
            else:
                # same timestep, multiple halos; accumulate
                final[start:end][valid] += hist_i[valid]
            previous_time = t_i
        return final


    def plot_xdelta(cls):
        return cls.tmax_Gyr/cls.nbins

    def plot_xlog(cls):
        return False

    def plot_ylog(cls):
        return False



class LiveHaloProperties(HaloProperties):
    def __init__(self, simulation, *args):
        super(LiveHaloProperties, self).__init__(simulation)
        self._nargs = len(args)

    @classmethod
    def requires_particle_data(self):
        return False

    def calculate(self, _, halo):
        return self.live_calculate(halo, *([None]*self._nargs))


class LiveHaloPropertiesInheritingMetaProperties(LiveHaloProperties):
    """LiveHaloProperties which inherit the meta-data (i.e. x0, delta_x values etc) from
    one of the input arguments"""
    def __init__(self, simulation, inherits_from, *args):
        """
        :param simulation: The simulation DB entry for this instance
        :param inherits_from: The HaloProperties description from which the metadata should be inherited
        :type inherits_from: HaloProperties
        """
        super(LiveHaloPropertiesInheritingMetaProperties, self).__init__(simulation)
        self._inherits_from = inherits_from(simulation)

    def plot_x0(self):
        return self._inherits_from.plot_x0()

    def plot_xdelta(self):
        return self._inherits_from.plot_xdelta()

class ProxyHalo(object):

    """Used to return pointers to halos within this snapshot to the database"""

    def __init__(self, value):
        self.value = value

    def __int__(self):
        return int(self.value)



##############################################################################
# UTILITY FUNCTIONS
##############################################################################

def all_property_classes():
    """Return list of all classes derived from HaloProperties"""

    return HaloProperties.all_classes()



def _check_class_provided_name(name):
    if "(" in name or ")" in name:
        raise ValueError("Property names must not include brackets; %s not suitable"%name)

def all_properties():
    """Return list of all properties which can be calculated using
    classes derived from HaloProperties"""
    classes = all_property_classes()
    pr = []
    for c in classes:
        try:
            i = c(None)
        except TypeError:
            continue
        name = i.name()
        if type(name) == str:
            _check_class_provided_name(name)
            pr.append(name)
        else:
            for name_j in name:
                _check_class_provided_name(name_j)
                pr.append(name_j)

    return pr


def providing_class(property_name, silent_fail=False):
    """Return providing class for given property name"""
    classes = all_property_classes()
    property_name = property_name.lower().split("(")[0]
    for c in classes:
        name = c.name()
        if isinstance(name, tuple) or isinstance(name, list):
            for name_j in name:
                if name_j.lower() == property_name:
                    return c
        elif name.lower() == property_name:
            return c
    if silent_fail:
        return None
    raise NameError("No providing class for property " + property_name)


def providing_classes(property_name_list, silent_fail=False):
    """Return providing classes for given list of property names"""
    classes = []
    for property_name in property_name_list:
        cl = providing_class(property_name, silent_fail)
        if cl not in classes and cl != None:
            classes.append(cl)

    return classes

def instantiate_classes(simulation, property_name_list, silent_fail=False):
    instances = []
    for property_identifier in property_name_list:
        instances.append(providing_class(property_identifier, silent_fail)(simulation))

    return instances

def instantiate_class(simulation, property_name, silent_fail=False):
    instance = instantiate_classes(simulation, [property_name],silent_fail)
    if len(instance)==0:
        return None
    else:
        return instance[0]

def get_required_properties(property_name):
    return providing_class(property_name).requires_property()

def import_configured_property_modules():
    from ..config import property_modules
    for pm in property_modules:
        if pm=="": continue
        try:
            importlib.import_module(pm)
        except ImportError:
            warnings.warn("Failed to import requested property module %r. Some properties may be unavailable."%pm,
                          ImportWarning)

import_configured_property_modules()
from . import live_profiles, intrinsic, BH, zoom, centring, profile, images
