TANGOS - The Amazing Numerical Galaxy Organisation System
---------------------------------------------------------

[![Build Status](https://travis-ci.com/N-BodyShop/halo_database.svg?token=Kwgna3AKWpdHTHRrmaYX&branch=master)](https://travis-ci.com/N-BodyShop/halo_database)

This repository contains the complete code for _tangos_, a system for building and querying databases summarising the
results of cosmological simulations. _Tangos_:

 - is written for python 2.7 and 3.5 or later;
 - can be customised to work with multiple python modules such as pynbody or yt to process raw simulation data;
 - uses sqlalchemy to store the resulting data and therefore can connect to many different flavours of database;
 - provides a web interface to the data;
 - allows users to construct efficient, science-focussed queries from python without typing a line of SQL.
 
 
Before you start
----------------
**This is a beta version of tangos**. When _tangos_ is released there will be an accompanying paper. 
_Tangos_ is GPL-licenced but good scientific practice requires you to acknowledge its use. Until the paper is
available please use the following acknowledgement or equivalent:

> This work made use of the _tangos_ analysis stack (Pontzen et al in prep); see www.github.com/pynbody/tangos.


Installation
------------

To install _tangos_ first clone the repository, then use the standard setuptools `install` command:

```
git clone git@github.com:pynbody/tangos.git 
cd tangos
python setup.py install
```

Alternatively if you intend to 

This should check for and install the _minimum_ prerequisites, but doesn't install _pynbody_. That's because _tangos_ is
written to be agnostic about how the underlying simulation snapshots are read so in principle you could use e.g. _yt_.
For all current tutorials, _pynbody_ is the preferred reading system and so for an easy life you should install it:

```
pip install git+ssh://git@github.com/pynbody/pynbody.git
```

Once installed, you should check that _tangos_ is functioning correctly by entering the `tests` folder and 
typing `nosetests`. You should see a bunch of text scrolling by, ultimately finishing with the simple message `OK`. 
If you get a failure message instead of `OK`, report it (with as much detail of your setup as possible) in the 
github issue tracker.

Setting up paths
----------------

By default tangos will look for raw simulation data in your home folder and create its database file there as well. 
If you don't want it to do this, you can set the environment variables `TANGOS_SIMULATION_FOLDER` (for the simulation folder) 
and `TANGOS_DB_CONNECTION` (for the database file). For example, in bash:

```
export TANGOS_SIMULATION_FOLDER=/path/to/simulation/data/
export TANGOS_DB_CONNECTION=~/scratch/Romulus/DatabaseFiles/cosmo25/data_romulus25.db 
```
or, in cshell:
```
setenv TANGOS_SIMULATION_FOLDER /nobackupp8/mtremmel/Romulus/
setenv TANGOS_DB_CONNECTION /nobackupp8/mtremmel/DataBaseFiles/romulus8/data_romulus8.db
```
The top line in each example points to the parent directory for all of your simulation data directories. 
If you don't have any simulations (i.e. you are just using a database object already created) then you 
should not have to worry about this variable. The second line points to the database object you wish to analyze;
by default this will be a sqlite file but you can also specify a 
[sqlalchemy URL](http://docs.sqlalchemy.org/en/latest/core/engines.html#database-urls).

Remember, you will need to set these environment variables *every* time you start a new session on your computer prior 
to booting up the database, either with the webserver or the python interface (see below).

Where next?
-----------

Now that you've set up the basics, you can either [make your first _tangos_ database](tutorials/first_steps.md) 
using some tutorial data or [download an existing database to perform data analysis](tutorials/data_exploration.md).