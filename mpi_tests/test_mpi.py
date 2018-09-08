#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import sys
import tangos.parallel_tasks as pt
import time
from six.moves import range


def test_function():
    lock = pt.ExclusiveLock("hello")

    print("Hello from rank",pt.backend.rank())
    for i in pt.distributed(range(10)):
        with lock:
            print("Task",i)
            time.sleep(0.1)

    if pt.backend.rank()==1:
        print()
        print("OK")

if len(sys.argv)!=2:
    print("Syntax: test_mpi.py [backend name]")
else:
    pt.use(sys.argv[1])
    pt.launch(test_function, 8)




