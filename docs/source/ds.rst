Data Store (ds) Module
=======================

The Data Store refers to an original source for raw data from an
experiment.  In the NDK context, data stores are usually local field
potential recordings, either in native format (e.g., Spike2, Plexon,
etc. data files) or in Cassandra, in the form of an apache Cassandra
database containing the time-series data recorded during an
experiment.  The fastest available data store for NDK is the NBF
format, which is the preferred way to preserve recording data and
metadata.

NDK's Data Store API can read files using python neo (hence any file
that neo supports) and WAV formats, as well as nbf.

.. toctree::
   :maxdepth: 3

   nbf
   neo_in
   probe
   wav
   edf_in



