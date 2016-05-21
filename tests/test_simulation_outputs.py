import halo_db as db
import halo_db.simulation_output_handlers.pynbody as pynbody_outputs
import halo_db.tools.add_simulation as add
import os
import numpy.testing as npt
import pynbody
import gc

def setup():
    global output_manager
    db.init_db("sqlite://")
    db.config.base = os.path.join(os.path.dirname(__name__), "test_simulations")
    output_manager = pynbody_outputs.ChangaOutputSetHandler("test_tipsy")


def test_get_handler():
    assert db.simulation_output_handlers.get_named_handler_class('pynbody.ChangaOutputSetHandler') == pynbody_outputs.ChangaOutputSetHandler

def test_handler_name():
    assert pynbody_outputs.ChangaOutputSetHandler.handler_class_name()=="pynbody.ChangaOutputSetHandler"

def test_enumerate():
    assert set(output_manager.enumerate_timestep_extensions())==set(["tiny.000640","tiny.000832"])

def test_timestep_properties():
    props = output_manager.get_timestep_properties("tiny.000640")
    npt.assert_allclose(props['time_gyr'],2.17328504831)
    npt.assert_allclose(props['redshift'], 2.96382819878)

def test_enumerate_halos():
    halos = list(output_manager.enumerate_halos("tiny.000640"))
    assert len(halos)==9
    assert halos[0]==[1,2041986, 364232, 198355]
    assert halos[1]==[2, 421027, 30282, 57684]

def test_properties():
    props = output_manager.get_properties()
    assert props['dPhysDenMin']==0.2 # from param file
    assert props['macros'].startswith("CHANGESOFT COOLING_COSMO") # from log file
    assert "dBHSinkAlpha" not in props # in the param file but not in the list of parameters we want to expose



def test_load_timestep():
    add_test_simulation_to_db()
    pynbody_f = db.get_timestep("test_tipsy/tiny.000640").load()
    assert isinstance(pynbody_f, pynbody.snapshot.SimSnap)
    assert pynbody_f.filename.endswith("tiny.000640")

def test_load_halo():
    add_test_simulation_to_db()
    pynbody_h = db.get_halo("test_tipsy/tiny.000640/1").load()
    assert isinstance(pynbody_h, pynbody.snapshot.SubSnap)
    assert len(pynbody_h)==200
    assert_is_subview_of_full_file(pynbody_h)

def test_load_tracker_halo():
    add_test_simulation_to_db()
    pynbody_h = db.get_halo("test_tipsy/tiny.000640/1.1").load()
    assert len(pynbody_h)==4

    # test that we have a subview of the whole file
    assert_is_subview_of_full_file(pynbody_h)

def test_partial_load_tracker_halo():
    add_test_simulation_to_db()
    pynbody_h = db.get_halo("test_tipsy/tiny.000640/1.1").load(partial=True)
    assert len(pynbody_h)==4
    assert pynbody_h.ancestor is pynbody_h

def test_load_persistence():
    f = db.get_timestep("test_tipsy/tiny.000640").load()
    f2 = db.get_timestep("test_tipsy/tiny.000640").load()
    h = db.get_halo("test_tipsy/tiny.000640/1").load()
    h_tracker = db.get_halo("test_tipsy/tiny.000640/1.1").load()
    assert id(f)==id(f2)
    assert id(h.ancestor)==id(f)
    assert id(h_tracker.ancestor)==id(f)

    old_id = id(f)

    del f, f2, h, h_tracker
    gc.collect()

    f3 = db.get_timestep("test_tipsy/tiny.000640").load()
    assert id(f3)!=old_id

def test_load_tracker_iord_halo():
    add_test_simulation_to_db()
    h_direct = db.get_halo("test_tipsy/tiny.000640/1.1").load(partial=True)
    h_iord = db.get_halo("test_tipsy/tiny.000640/1.2").load(partial=True)
    assert (h_direct['iord']==h_iord['iord']).all()


_added_to_db = False
tracked_particles = [2, 4, 6, 8]
tracked_iord = [20000,40000,60000,80000]

def add_test_simulation_to_db():
    global _added_to_db

    if not _added_to_db:
        add.SimulationAdderUpdater(output_manager).scan_simulation_and_add_all_descendants()
        tx = db.core.tracking.TrackData(db.get_simulation("test_tipsy"))
        tx.particles = tracked_particles
        tx.use_iord = False
        tx = db.core.get_default_session().merge(tx)
        tx.create_halos()

        tx = db.core.tracking.TrackData(db.get_simulation("test_tipsy"))
        tx.particles = tracked_iord
        tx.use_iord = True
        tx = db.core.get_default_session().merge(tx)
        tx.create_halos()
        _added_to_db=True

def assert_is_subview_of_full_file(pynbody_h):
    assert len(pynbody_h.ancestor) == len(db.get_timestep("test_tipsy/tiny.000640").load())

