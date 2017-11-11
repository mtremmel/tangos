from __future__ import absolute_import

import glob
import os
import os.path
import time
import weakref
import re
import numpy as np
import pynbody

from . import halo_stat_files, finding
from . import SimulationOutputSetHandler
from .. import config
from ..log import logger
from ..parallel_tasks import pynbody_server as ps
from six.moves import range


_loaded_halocats = {}

class DummyTimeStep(object):
    def __init__(self, filename):
        self.filename = filename


    def __repr__(self):
        return self.filename

    pass

class PynbodyOutputSetHandler(SimulationOutputSetHandler):
    patterns = [] # should be specified by child class

    @classmethod
    def best_matching_handler(cls, basename):
        handler_names = []
        handler_timestep_lengths = []
        base = os.path.join(config.base, basename)
        if len(cls.__subclasses__())==0:
            return cls
        for possible_handler in cls.__subclasses__():
            timesteps_detected = finding.find(basename = base+"/", patterns = possible_handler.patterns)
            handler_names.append(possible_handler)
            handler_timestep_lengths.append(len(timesteps_detected))
        return handler_names[np.argmax(handler_timestep_lengths)]



    def enumerate_timestep_extensions(self):
        base = os.path.join(config.base, self.basename)
        extensions = finding.find(basename=base + "/", patterns=self.patterns)
        for e in extensions:
            if self._pynbody_can_load_halos_for(e):
                yield e[len(base)+1:]

    def _pynbody_can_load_halos_for(self, filepath):
        try:
            f = pynbody.load(filepath)
            h = f.halos()
            return True
        except (IOError, RuntimeError):
            return False

    def get_timestep_properties(self, ts_extension):
        ts_filename =  self._extension_to_filename(ts_extension)
        f = pynbody.load(ts_filename)
        try:
            time_gyr = f.properties['time'].in_units("Gyr",**f.conversion_context())
        except:
            time_gyr = -1

        results = {'time_gyr': time_gyr, 'redshift': float(f.properties['z']),
                   'available': True}
        return results

    def load_timestep_without_caching(self, ts_extension, mode=None):
        if mode=='partial' or mode is None:
            f = pynbody.load(self._extension_to_filename(ts_extension))
            f.physical_units()
            return f
        elif mode=='server' or mode=='server-partial':
            return ps.RemoteSnapshotConnection(self._extension_to_filename(ts_extension))
        else:
            raise NotImplementedError("Load mode %r is not implemented"%mode)

    def load_region(self, ts_extension, region_specification, mode=None):
        if mode is None:
            timestep = self.load_timestep(ts_extension)
            return timestep[region_specification]
        elif mode=='server':
            timestep = self.load_timestep(ts_extension)
            return timestep.get_view(region_specification)
        elif mode=='server-partial':
            timestep = self.load_timestep(ts_extension)
            view = timestep.get_view(region_specification)
            load_index = view['remote-index-list']
            logger.info("Partial load %r, taking %d particles",ts_extension,len(load_index))
            f = pynbody.load(self._extension_to_filename(ts_extension), take=load_index)
            f.physical_units()
            return f
        elif mode=='partial':
            raise NotImplementedError("For partial loading to work with custom regions, you need load-mode=server-partial (instead of load-mode=partial)")
        else:
            raise NotImplementedError("Load mode %r is not implemented"%mode)

    def load_object(self, ts_extension, halo_number, object_typetag='halo', mode=None):
        if mode=='partial':
            h = self._construct_halo_cat(ts_extension, object_typetag)
            h_file = h.load_copy(halo_number)
            h_file.physical_units()
            return h_file
        elif mode=='server':
            timestep = self.load_timestep(ts_extension)
            return timestep.get_view(halo_number)
        elif mode=='server-partial':
            timestep = self.load_timestep(ts_extension)
            view = timestep.get_view(halo_number)
            load_index = view['remote-index-list']
            logger.info("Partial load %r, taking %d particles", ts_extension, len(load_index))
            f = pynbody.load(self._extension_to_filename(ts_extension), take=load_index)
            f.physical_units()
            return f
        elif mode is None:
            h = self._construct_halo_cat(ts_extension, object_typetag)
            return h[halo_number]
        else:
            raise NotImplementedError("Load mode %r is not implemented"%mode)

    def load_tracked_region(self, ts_extension, track_data, mode=None):
        f = self.load_timestep(ts_extension)
        indices = self._get_indices_for_snapshot(f, track_data)
        if mode=='partial':
            return pynbody.load(f.filename, take=indices)
        elif mode is None:
            return f[indices]
        else:
            raise NotImplementedError("Load mode %r is not implemented"%mode)


    def _get_indices_for_snapshot(self, f, track_data):
        pt = track_data.particles
        if track_data.use_iord is True:

            dm_part = f.dm[np.in1d(f.dm['iord'], pt)]

            try:
                star_part = f.star[np.in1d(f.star['iord'], pt)]
            except KeyError:
                star_part = f[0:0]

            try:
                gas_part = f.gas[np.in1d(f.gas['iord'], pt)]
            except KeyError:
                gas_part = f[0:0]

            # fx = dm_part.union(star_part)
            # fx = fx.union(gas_part)
            # return fx
            ilist = np.hstack((dm_part.get_index_list(f),
                               star_part.get_index_list(f),
                               gas_part.get_index_list(f)))
            ilist = np.sort(ilist)
            return ilist
        else:
            return pt



    def _construct_halo_cat(self, ts_extension, object_typetag):
        f = self.load_timestep(ts_extension)
        loaded_cats = getattr(f, '_tangos_cached_cats', None)

        if loaded_cats is None:
            loaded_cats = f._cached_cats = {}

        h = loaded_cats.get(object_typetag, None)
        if h is None:
            h = self._construct_halo_cat_without_caching(f, object_typetag)
            loaded_cats[object_typetag] = h # keeps alive for lifetime of simulation
        return h

    def _construct_halo_cat_without_caching(self, f, object_typetag):
        if object_typetag!= 'halo':
            raise ValueError("Unknown object type %r" % object_typetag)
        h = f.halos()
        if isinstance(h, pynbody.halo.SubfindCatalogue):
            h = f.halos(subs=True)
        return h

    def match_halos(self, ts1, ts2, halo_min, halo_max,
                    dm_only=False, threshold=0.005, object_typetag='halo'):
        if dm_only:
            only_family='dm'
        else:
            only_family=None

        f1 = self.load_timestep(ts1)
        f2 = self.load_timestep(ts2)

        h1 = self._construct_halo_cat(ts1, object_typetag)
        h2 = self._construct_halo_cat(ts2, object_typetag)

        return f1.bridge(f2).fuzzy_match_catalog(halo_min, halo_max, threshold=threshold,
                                                 only_family=only_family, groups_1=h1, groups_2=h2)

    def enumerate_objects(self, ts_extension, object_typetag="halo", min_halo_particles=config.min_halo_particles):
        ts = DummyTimeStep(self._extension_to_filename(ts_extension))
        ts.redshift = self.get_timestep_properties(ts_extension)['redshift']
        try:
            statfile = halo_stat_files.HaloStatFile(ts)
            if object_typetag != "halo":
                raise StopIteration
            logger.info("Reading halos for timestep %r using a stat file",ts)
            for X in statfile.iter_rows("n_dm", "n_star", "n_gas"):
                yield X
        except IOError:
            logger.warn("No halo statistics file found for timestep %r",ts)

            try:
                h = self._construct_halo_cat(ts_extension, object_typetag)
            except:
                logger.warn("Unable to read %ss using pynbody, skipping step", object_typetag)
                raise StopIteration

            logger.info("Reading %ss directly using pynbody", object_typetag)

            istart = 1

            if isinstance(h, pynbody.halo.SubfindCatalogue):
                istart = 0 # subfind indexes from zero

            if hasattr(h, 'precalculate'):
                h.precalculate()


            for i in range(istart, len(h)+istart):
                try:
                    hi = h[i]
                    if len(hi.dm) > min_halo_particles:
                        yield i, len(hi.dm), len(hi.star), len(hi.gas)
                except (ValueError, KeyError) as e:
                    pass

    def get_properties(self):
        timesteps = list(self.enumerate_timestep_extensions())
        if len(timesteps)>0:
            f = self.load_timestep_without_caching(timesteps[0])
            return {'approx_resolution_kpc': self._estimate_resolution(f)}
        else:
            return {}

    def _estimate_resolution(self, f):
        f.physical_units()
        if "eps" in f.dm.loadable_keys():
            # Interpret the eps array as a softening, and assume that it is not a comoving softening (as the
            # pynbody units system might naively tell us) but actually already a physical softening. Note that
            # whether or not this is a correct interpretation depends on the code in use, and the flags passed
            # to that code.
            return float(f.dm['eps'].in_units('kpc a').min())
        else:
            # There is no softening information available, so take a best guess as to what a reasonable
            # softening might be as 1/100th of the mean interparticle distance (in the deepest zoom level)
            tot_box_mass = f.dm['mass'].sum()
            min_mass = f.dm['mass'].min()
            frac_mass = min_mass/tot_box_mass
            frac_length = frac_mass ** (1. / 3)
            estimated_eps = 0.01 * frac_length * f.properties['boxsize'].in_units('kpc a', **f.conversion_context())
            return float(estimated_eps)


class RamsesHOPOutputSetHandler(PynbodyOutputSetHandler):
    patterns = ["output_0????"]



class GadgetSubfindOutputSetHandler(PynbodyOutputSetHandler):
    patterns = ["snapshot_???"]

    def load_object(self, ts_extension, halo_number, object_typetag='halo', mode=None):
        if mode=='subfind_properties':
            h = self._construct_halo_cat(ts_extension, object_typetag)
            return h.get_halo_properties(halo_number,with_unit=False)
        else:
            return super(GadgetSubfindOutputSetHandler, self).load_object(ts_extension, halo_number, object_typetag, mode)

    def enumerate_timestep_extensions(self):
        base = os.path.join(config.base, self.basename)
        extensions = finding.find(basename=base + "/", patterns=["snapshot_???"])
        for e in extensions:
            if self._pynbody_can_load_halos_for(e):
                yield e[len(base)+1:]

    def _construct_halo_cat_without_caching(self, f, object_typetag):
        if object_typetag== 'halo':
            h = f.halos(subs=True)
        elif object_typetag== 'group':
            h = f.halos(subs=False)
        else:
            raise ValueError("Unknown halo type %r" % object_typetag)

        assert isinstance(h, pynbody.halo.SubfindCatalogue)
        return h

class EagleLikeHDFOutputSetHandler(PynbodyOutputSetHandler):
    patterns = ["snap_???_*.*.hdf5"]

    def enumerate_timestep_extensions(self):
        base = os.path.join(config.base, self.basename)
        extensions = finding.find(basename=base + "/", patterns=self.patterns)
        # join together different CPUs
        extensions_no_cpu = set()
        for e in extensions:
            matched_name = re.match("(.*)\.[0-9]+\.hdf5$", e)
            if matched_name:
                extensions_no_cpu.add(matched_name.group(1))
            else:
                logger.warn("Could not determine the snapshot that file %r belongs to", e)

        for e in extensions_no_cpu:
            if self._pynbody_can_load_halos_for(e):
                yield e[len(base)+1:]

    def _construct_halo_cat_without_caching(self, f, object_typetag):
        if object_typetag == 'halo':
            h = f.halos(subs=True)
        elif object_typetag == 'group':
            h = f.halos(subs=False)
        else:
            raise ValueError("Unknown halo type %r" % object_typetag)

        assert isinstance(h, pynbody.halo.GrpCatalogue)
        return h



class ChangaOutputSetHandler(PynbodyOutputSetHandler):
    flags_include = ["dPhysDenMin", "dCStar", "dTempMax",
                     "dESN", "bLowTCool", "bSelfShield", "dExtraCoolShutoff"]

    patterns = ["*.00???","*.00????"]


    def get_properties(self):
        parent_prop_dict = super(ChangaOutputSetHandler, self).get_properties()

        pfile = self._get_paramfile_path()

        if pfile is None:
            logger.warn("Param file cannot be found - no simulation properties will be available")
            return {}
        else:
            logger.info("Param file is %s", pfile)

        pfile_dict = self._param_file_to_dict(pfile)
        log_path, prop_dict = self._get_log_path(pfile, pfile_dict)

        if log_path:
            prop_dict.update(self._get_properties_from_log(log_path))

        prop_dict.update(self._filter_paramfile_properties(pfile_dict))
        prop_dict.update(parent_prop_dict)

        return prop_dict

    def _get_paramfile_path(self):
        try:
            pfile = self._get_param_file_for_output(self._extension_to_filename(""))
        except RuntimeError:
            pfile = None
        return pfile

    def _get_log_path(self, paramfile_name, paramfile_dict):
        prop_dict = {}
        log_fn = paramfile_dict.get("achOutName","") + ".log"
        log_path = paramfile_name.split("/")[:-1]
        log_path.append(log_fn)
        log_path = "/".join(log_path)
        if os.path.exists(log_path):
            logger.info("Log file is %s", log_path)
        else:
            logger.warn("Cannot find log file (%s)", log_path)
            log_path = None
        return log_path, prop_dict

    def _filter_paramfile_properties(self, pfile_dict):
        filtered_pfile_dict = {}
        for f in self.flags_include:
            if f in pfile_dict:
                filtered_pfile_dict[f] = pfile_dict[f]
        return filtered_pfile_dict

    def _get_properties_from_log(self, log_path):
        prop_dict = {}
        with open(log_path, 'r') as f:
            for l in f:
                if "# Code compiled:" in l:
                    prop_dict["compiled"] = time.strptime(
                        l.split(": ")[1].strip(), "%b %d %Y %H:%M:%S")
                if "# Preprocessor macros: " in l:
                    prop_dict["macros"] = l.split(": ")[1].strip()
                    break
        return prop_dict

    @staticmethod
    def _get_param_file_for_output(output_file):
        """Work out the param file corresponding to the
        specified output"""

        q = "/".join(output_file.split("/")[:-1])
        if len(q) != 0:
            path = "/".join(output_file.split("/")[:-1]) + "/"
        else:
            path = ""

        candidates = glob.glob(path + "*.param")

        if len(candidates) == 0:
            candidates = glob.glob(path + "../*.param")

        if len(candidates) == 0:
            raise RuntimeError("No .param file in " + path + \
                                " (or parent) -- please supply or create tipsy.info manually")

        candidates = [x for x in candidates if "direct" not in x and "mpeg_encode" not in x]

        if len(candidates) > 1:
            raise RuntimeError("Can't resolve ambiguity -- too many param files matching " + \
                                path)

        return candidates[0]

    @staticmethod
    def _param_file_to_dict(param_file):
        f = open(param_file)
        out = {}

        for line in f:
            try:
                s = line.split()
                if s[1] == "=" and "#" not in s[0]:
                    key = s[0]
                    v = s[2]

                    if key[0] == "d":
                        v = float(v)
                    elif key[0] == "i" or key[0] == "n" or key[0] == "b":
                        v = int(v)

                    out[key] = v
            except (IndexError, ValueError):
                pass
        return out
