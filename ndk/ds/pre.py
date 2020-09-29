import numpy as np
import glob
import math
import ndk.features
from uritools import urisplit, uriunsplit

#
# This file contains code to read in data from .prexx files.  This is
# an old one-off format used by the Graybiel lab (borrowed from a
# similar format used by the Wilson lab) at MIT circa 1995.  These are
# not LFPs, but rather threshold-triggered windows of spikes as
# recorded from tetrodes.
#
# Returns a list of things that can be written out to a db3 event
# store.  This function also reads the data sample window, and it's
# not clear how that should be handled.  A few possibilities:
#
# 1) Represent each spike event in basis 0 (the default basis).
#
# 2) Try to partition (cluster) these spikes and create other bases
# that provide better fits to the data.
#
# 3) Export the data in some "reasonable" way to Cassandra and do more
# research into representations.

def read_tuple(f):
    """Read and interpret a tuple from a .pre?? file header."""
    pair = f.readline().split(': ')
    key = pair[0]
    line = pair[1].split(' ')
    val = []
    for x in line:
        val.append(int(x))
    return key.lower(), val


def parse_waveforms(list, gain=10000.0, zero_mean=True):
    """Parse a set of waveform data points that appeared in a tetrode .pre file.  Returns a 4x32 Numpy array of ints."""
    w = []
    w.append(np.zeros(32))
    w.append(np.zeros(32))
    w.append(np.zeros(32))
    w.append(np.zeros(32))
    k = 0
    for j in range(32):
        for i in range(4):
            w[i][j] = float(list[k]) / gain
            k += 1
    if zero_mean:
        w[0] = ndk.features.zero_mean(w[0])
        w[1] = ndk.features.zero_mean(w[1])
        w[2] = ndk.features.zero_mean(w[2])
        w[3] = ndk.features.zero_mean(w[3])
    return w


def read_event_file(filename):
    """Reads a .event file *filename* and returns a dict containing the header and records.  Each record contains a timestamp and an external event name."""
    master = []
    i = 0
    hdr = []
    start = False
    with open(filename, 'r') as f:
        for line in f:
            if not start:
                hdr.append(line)
                start = 'Series' in line
            else:
                rec = line.split(' ')
                ts = int(rec[0])
                tag = rec[1]
                master.append((ts,tag))
            i += 1
    return master, hdr
        

def read_pre_file(filename):
    """Reads a .pre file *filename* and returns a dict containing the header and records.  Each record contains a timestamp, tetrode number, and 128 waveform data points."""
    master = {}
    with open(filename, 'r') as f:
        id = f.readline()
        master['file_id'] = id

        tetrode_count = int(f.readline().split(': ')[1])
        master['tetrode_count'] = tetrode_count

        ts_line = f.readline().split(': ')[1]
        ts = ts_line.split(' ')
        t0 = int(ts[0])
        t1 = int(ts[1])
        master['limits'] = [t0, t1]

        recs = (f.readline().split(': ')[1]).split(' ')
        master['records'] = [int(recs[0]), int(recs[1])]

        mins = (f.readline().split(': ')[1]).split(' ')
        master['mins'] = [int(mins[0]), int(mins[1])]

        maxes = (f.readline().split(': ')[1]).split(' ')
        master['maxes'] = [int(maxes[0]), int(maxes[1])]

        rate = float(f.readline().split(': ')[1])
        master['rate'] = rate / (tetrode_count * 4)

        # Electrodes:
        key, val = read_tuple(f)
        master[key] = val

        # Gains:
        key, val = read_tuple(f)
        master[key] = val

        # Impedances:
        key, val = read_tuple(f)
        master[key] = val

        # Low cutoff:
        key, val = read_tuple(f)
        master[key] = val

        # High cutoff:
        key, val = read_tuple(f)
        master[key] = val

        series = f.readline().split(' ; ')
        l = [ series[0].split(' ')[1], series[1] ]
        master['series'] = l

        k = 0
        for line in f:
            vals = line.split(' ')
            ts = vals[0]
            tetrode = vals[1]
            waveforms = parse_waveforms(vals[2:])
            master[k] = [ts, tetrode, waveforms]
            k += 1
            # add_to_es(ts, tetrode, waveforms)

    return master
            
def mean_event_times(elist):
    # The event files don't always even have all of the events.
    # Clearly, there were bugs in the event collection code, so absent
    # any other criteria, just use the first event file:
    return elist[0]

# The correct way to use this is to pass either a filename string or a
# URI up to, but not including the suffixes, e.g.:
#         'pre://localhost/opt/spikes/060597aa'
# or      '/opt/spikes/060597aa'

class pre_ds():
    """Data store class for .pre files.  Initialized by passing in the
name of the .pre file without the tetrode numbers.  The init routine
will populate the data store with data from all available tetrodes."""
    def __init__(self, uri, sigma=0.15, window=1.024):
        if isinstance(uri, str):
            self.desc = urisplit(uri)
        else:
            self.desc = uri
        basename = self.desc.path
        print(basename)
        self.sigma = sigma
        self.window = window
        self.pre_files = glob.glob(basename+".pre*")
        self.event_files = glob.glob(basename+".events*")
        print("Pre files: {}".format(self.pre_files))
        self.info = []
        events = []
        self.eventhdr = []
        for fname in self.pre_files:
            m = read_pre_file(fname)
            self.info.append(m)
        for fname in self.event_files:
            m, hdr = read_event_file(fname)
            events.append(m)
            self.eventhdr.append(hdr)
        self.events = mean_event_times(events)

        
    def get_signal(self, channel, start=None, end=None, into_vec=None):
        print("get_signal is not yet available for .pre data sources.")

    def get_chunk(self, channel, start, end):
        """Returns the portion of the time series for the specified channel that is between start and end inclusive (time)."""
        which = int(channel / 4)
        m = self.info[which]
        # print("Interval: [{}, {}]".format(start,end))
        for key, value in m.items():
            ts = int(value[0]) - 16
            # print("Looking: ts={}".format(ts))
            if ts >= start and ts <= end:
                waveforms = value[2]
                ic = channel % 4
                return waveforms[ic]
        
    def filter_data(self, low_cutoff, high_cutoff):
        print('filter_data() for pre stores is not yet implemented!')
        

    # We also need to define a get_chunk function that will return to
    # the caller a waveform vector.
    def fill_event_store(self, es, tetrode=0):
        """Given an event store object 'es', populate it from the information
contained in this pre_ds data store object.  Only handles one tetrode
at a time."""
        samprate = self.info[0]['rate']
        nchannels = len(self.info) * 8
        # If we separate the tetrodes, then only use 4 channels:
        es.make_meta_table(samprate, 4, uriunsplit(self.desc))
        # See if the default sigma (0.32 ms) and window width (2ms)
        # will work here:

        es.make_basis_table()
        es.make_event_table()
        es.make_spiketimes_table()
        es.make_spikebasis_table()
        # es.make_channel_table()
        es.make_spikecoefs_table()
        es.make_spikelabels_table()
        es.make_default_labels()

        es.add_basis(0, self.sigma, self.window, nfuncs=7, name='pre_basis')
        seq = es.get_basis(0)
        spike_id = 0
        which = int(tetrode / 2)
        tnum = tetrode % 2
        m = self.info[which]
        # Does not handle behavioral events:
        for key, value in m.items():
            # Only consider integer keys - these are data records:
            if isinstance(key, int):
                ts, tetrode, w = value
                t = int(tetrode)
                if t == tnum:
                    es.set_spike_time(spike_id, ts)
                    es.set_spike_basisid(spike_id, 0)
                    es.set_spike_label(spike_id, 0)
                    for j in range(4):
                        chan = j
                        signal = w[j]
                        coefs = ndk.features.waveform_coefs(signal, seq)
                        es.set_spike_coefs(spike_id, chan, coefs)
                    spike_id += 1
        for ts, tag in self.events:
            es.add_event(tag, ts)

        
