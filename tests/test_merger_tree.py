from __future__ import absolute_import
from __future__ import print_function
import tangos.core.dictionary
import tangos.core.halo
import tangos.core.halo_data
import tangos.core.simulation
import tangos.core.timestep
import tangos
import tangos.testing.simulation_generator

from tangos.relation_finding import tree
import tangos.testing as testing
from nose.tools import assert_raises

def setup():
    testing.init_blank_db_for_testing()

    generator = tangos.testing.simulation_generator.TestSimulationGenerator()
    generator.add_timestep() # ts1
    generator.add_objects_to_timestep(7)

    generator.add_timestep() # ts2
    generator.add_objects_to_timestep(7)
    generator.link_last_halos()

    generator.add_timestep() # ts3
    generator.add_objects_to_timestep(6)
    generator.link_last_halos_using_mapping({1:1, 2:1, 3:2, 4:3, 5:4, 6:5, 7:6}) # ts2->ts3: merger of halos 1 & 2

    generator.add_timestep() # ts4
    generator.add_objects_to_timestep(6)
    generator.link_last_halos()

    generator.add_timestep() # ts5
    generator.add_objects_to_timestep(5)
    generator.link_last_halos_using_mapping({1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 1})  # ts4->ts5: merger of halos 1 & 6

    generator.add_timestep() # ts6
    generator.add_objects_to_timestep(5)
    generator.link_last_halos()

def test_default_tree_has_correct_structure():
    mt = tree.MergerTree(tangos.get_halo("%/ts6/1"))
    mt.construct()
    assert mt.summarise()=="1(1(1(1(1(1),2(2))),6(6(7(7)))))"

    mt = tree.MergerTree(tangos.get_halo("%/ts6/2"))
    mt.construct()
    assert mt.summarise() == "2(2(2(2(3(3)))))"

def test_filter_tree_by_minweight():
    old = tree.mergertree_min_fractional_weight
    try:
        tree.mergertree_min_fractional_weight = 0.8
        mt = tree.MergerTree(tangos.get_halo("%/ts6/1"))
        mt.construct()
        assert mt.summarise() == "1(1(1(1(1(1),2(2)))))"
    finally:
        tree.mergertree_min_fractional_weight = old

def test_filter_tree_by_NDM():
    tree.mergertree_min_fractional_NDM = 0.2
    mt = tree.MergerTree(tangos.get_halo("%/ts6/1"))
    mt.construct()
    assert mt.summarise()=="1(1(1(1(1(1),2(2))),6))"