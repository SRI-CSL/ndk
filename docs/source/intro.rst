Introduction
============

The Neurogram Deconvolution Kit (ndk) is a software suite for
warehousing, compression, visualization, and analysis of peripheral
nerve neurograms (spike trains).  The use of Python affords us a
certain degree of portability, and a large degree of interoperability
with other languages or systems (e.g., Matlab).  Much of NDK's core
concepts are derived from earlier efforts to organize and analyze
multiunit recording data from operant conditioning experiments (Jog,
et al.).  These data included recordings from cortex and striatum, as
well as collection of behavioral markers.  The NDK can handle such
data, but also contains analysis tools suited for analysis of signals
taken from peripheral nerves.

The NDK itself is a Python package, rewritten to accomodate a variety
of signal types in both the peripheral and central nervous systems.
NDK appeals to classical function approximation theory to represent
discrete signals, such as neuronal spikes or EEG events, in terms of
basis functions that are specific to the phenomena encountered in a
data set.  Events and their corresponding basis sets may occur at
different time scales, or require different temporal support.  For
example, spikes emanating from vagus nerve A fibers have a waveform
shape and temporal extent that is distinct from C fiber spikes.
Moreover, NDK needs to accomodate an arbitrary channel count, as well
as a variety of behavioral and stimulus markers that may be present in
a wide variety of experimental contexts.  In some cases, multiple
channel configurations provide information about the direction and
speed of action potential propagation.  NDK's design (especially the
basis set representation) is intended to assist directionality
analysis.  All of these issues motivate the development within NDK of
relational data models that must support our needs, and the consequent
use of modern database technology.

NDK separates the concepts of data store and event store.  The former
is responsible for accessing and serving dense time-series data, while
the latter is responsible for access, analysis, and visualization of
important events that are discovered in the dense time-series data.

Data Store
----------

A "data store" implements storage and fast direct access for
high-volume time series data.  In this context, a data store could be
as simple as a recording file (e.g., a local field potential recording
in Plexon or Spike2 formats), or as flexible and accessible as a
public Apache Cassandra
(http://cassandra.apache.org)
database.  The NDK exploits the python "neo" package
(https://pypi.python.org/pypi/neo)
for data store I/O using files, but generally, the NDK should export
such data into a Cassandra data store to provide a uniform interface
to time series data.  Cassandra allows for a standard, SQL-like direct
access to recorded samples..

The current preferred data store type for NDK is NBF - NDK Binary
Format.  This is a memory-mapped binary format meant for fast local
access to recording data and metadata.  Think of this as an internal
representation through which data can be transcoded for various
purposes.  NBF supports multiple channels, revisioning and a
pyramid-based storage scheme.


Event Store
-----------

An "event store", on the other hand, is an efficient representation
for discrete events found in the data store, characterized by a time
stamp and a "recipe" for reconstruction of the original signal.  Each
event class is represented as a linear combination of a basis set.
Different basis sets can be constructed for different phenomena (e.g.,
spikes emanating from A, B, and C fibers).  The event store is very
much like a classical database, allowing efficient access, analysis,
and visualization of relational data derived from the denser
time-series data.  In addition, the event store holds all the metadata
associated with the experiment from which the events were derived,
including links back to the raw data sources.

The events of greatest interest in neurogram analysis are neuronal
spikes, so the NDK has been developed in this context.  We anticipate,
however, that other events in other neuronal data stores could be
represented using NDK.  This includes phenomena like slow waves and
spindles in EEG datasets.  Nonetheless, our primary focus is currently
on peripheral nerve bundles.

The current preferred event store type is a db3 - sqlite file.  The
combination of sqlite event store with NBF data store is meant to be
portable (modulo endian-ness) and locally re-rootable (i.e., pathnames
are relative).


Events
------

One of the driving philosophies behind NDK is that each "event" (e.g.,
an action potential), is manifested in the time series recording as a
recognizable feature that can be represented using classical function
approximation techniques.  That is, punctate events, whether they are
spikes from diverse neuronal sources, or spindles (EEG), or slow waves
(also EEG), can be represented using classical approximation theory,
and marked with the function basis ID used for reconstruction.  Each
basis is a finite set of N orthonormal functions.  The inner product
of the raw signal and each basis function provides up to N
coefficients that represent the event in that basis.  These
coefficients are then used to represent the event for most analysis
(e.g., spike clustering), and provide a least-squares approximation of
the original signal with respect to the basis.  For spike events, by
default the Gaussian and its derivatives are used to generate each
basis, parameterized by variance and window width in milliseconds.


References
----------

Jog, M. S., Y. Kubota, C. Connolly, V. Hillegaart, A. M. Graybiel (1999) Building Neural Representations of Habits. *Science*, vol. 286, no. 5445, Nov 1999.
