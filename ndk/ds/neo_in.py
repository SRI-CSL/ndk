from __future__ import print_function
from uritools import urisplit, uriunsplit, urijoin
from ndk.features import butter_low, butter_high, butter_bandpass
from scipy.signal import convolve, butter, lfilter, filtfilt
import numpy as np
import neo.io
import os
from .dsfilter import *

def spike2_to_floats(filename):
    """Convert the time series in a segment of a Spike2 file into an array of floats.  Also returns the rate, start time, and end time for the data."""
    print('Converting {}'.format(filename))
    r = neo.io.Spike2IO(filename=filename)
    # print('neo.io.Spike2IO returns {}'.format(r))
    hdr = r.header
    sig_channels = hdr['signal_channels']
    print(sig_channels)
    seg = r.read_segment()
    signals = seg.analogsignals
    print("signals[0] shape: {}".format(signals[0].shape))
    nchannels = len(signals)
    print('nchannels = {}'.format(nchannels))
    samprate = float(signals[0].sampling_rate)
    # len(signals) is nchannels.  len(signals[0]) is the number of
    # samples:
    print('spike2 data interval: [ {}, {} ]'.format(signals[0].t_start, signals[0].t_stop))
    t0 = int(signals[0].t_start * samprate)
    t1 = int(signals[0].t_stop * samprate)
    print('spike2 sample num interval: [ {}, {} ]'.format(t0, t1))
    print('spike2 signal array length: {}'.format(len(signals[0])))
    # Here, signals[i] represents an entire time series for channel i:
    return signals, samprate, t0, t1, seg


class spike2():
    def __init__(self, uri):
        if isinstance(uri, str):
            self.desc = urisplit(uri)
        else:
            self.desc = uri
            # self.filename = filename
        self.filename = self.desc.path
        print('# Filename: {}'.format(self.filename))
        sig, samprate, start, end, seg = spike2_to_floats(self.filename)
        self.seg = seg
        self.signal = sig
        self.samprate = samprate
        self.t0 = start
        self.t1 = end
        print('# Dataset interval: [{}, {}]'.format(self.t0, self.t1))
        self.dsf = dsfilter(samprate)
                
    def filter_data(self, low_cutoff, high_cutoff, order=5):
        self.dsf.filter_data(low_cutoff, high_cutoff, order)

    def get_signal(self, channel):
        return self.signal[channel]

    
    def get_chunk(self, channel, start, end):
        """Returns the portion of the time series for the specified channel that is between start and end inclusive (time)."""
        x = self.signal[channel]
        i0 = max(0, int(start - self.t0))
        #i1 = max(0, int(end+1 - self.t0))
        i1 = max(0, int(end - self.t0))
        # We should probably call numpy.asarray to coerce this to a pure array:
        r = np.asarray(x[i0:i1])
        result = self.dsf.maybe_filter(r)
        result = result.reshape(len(r))
        return result

    

    def fill_event_store(self, es, sigma=0.32, window=2.0, nfuncs=5):
        """Given an event store object 'es', populate it from the information
contained in this pre_ds data store object.  Only handles one tetrode
at a time."""
        samprate = self.samprate
        nchannels = len(self.signal)
        samp_sigma = sigma * 0.001 * samprate
        samp_window = window * 0.001 * samprate
        es.make_meta_table(samprate, nchannels, uriunsplit(self.desc))
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

        es.add_basis(0, samp_sigma, samp_window)

        seq = es.get_basis(0)

        zc = find_feature_peaks(self.signal[0], seq)
        last = len(zc)

        i = 0
        iprev = 0
        spike_id = -1
        prev_id = 0
        while i < last:
            spikep = (zc[i] != 0)
            if spikep:
                ts = t0 + i*dt
                if ( (i - iprev) > width ):
                    spike_id += 1
                    es.set_spike_time( spike_id, ts )
                    es.set_spike_basisid( spike_id, 0 )
                    es.set_spike_label(spike_id, 0)
                    for chan in range(nchannels):
                        x = ndk.features.get_window(signal[chan], ts, winsize)
                        coefs = ndk.features.waveform_coefs(w, seq)
                        es.set_spike_coefs(spike_id, chan, coefs)
                    iprev = i
            i = i + 1
            if spike_id - prev_id >= 100:
                if spikeid-previd >= 100:
                    print( '{} '.format(spikeid), end="")
                    prev_id = spike_id

                    for j in range(4):
                        chan = j
                        signal = w[j]
                        coefs = ndk.features.waveform_coefs(signal, seq)
                        es.set_spike_coefs(spike_id, chan, coefs)
                    spike_id += 1
        for ts, tag in self.events:
            es.add_event(tag, ts)



def neuralynx_to_floats(filename):
    """Convert the time series in a segment of a Neuralynx file into an array of floats.  Also returns the rate, start time, and end time for the data."""
    dir = os.path.dirname(filename)
    # r = neo.io.NeuralynxIO(dir, use_cache='never')
    r = neo.io.NeuralynxIO(dir, use_cache=True, cache_path='home')
    seg = r.read_segment()
    signals = seg.analogsignals  # nchan x nsamp - opposite of spike2.
    sigout = signals[0].transpose()
    nchannels = sigout.shape[0]
    print('Number of channels: {} ; Sample rate: {}'.format(nchannels, signals[0].sampling_rate))
    samprate = float(signals[0].sampling_rate.simplified)
    # len(signals) is nchannels.  len(signals[0]) is the number of
    # samples:
    print('neuralynx data interval: [ {}, {} ]'.format(signals[0].t_start, signals[0].t_stop))
    t0 = int(signals[0].t_start.simplified * samprate)
    t1 = int(signals[0].t_stop.simplified * samprate)
    print('neuralynx sample num interval: [ {}, {} ]'.format(t0, t1))
    print('neuralynx signal array length: {}'.format(len(sigout[0])))
    # Not good enough - misses start / end recording events:
    # events = r.get_event_timestamps()
    events = []
    for x in seg.events:
        times = x.times
        labels = x.labels
        for i in range(len(times)):
            ts = times[i]
            label = labels[i]
            events.append((int(float(ts)*samprate), label))
    print(events)
    return sigout, samprate, t0, t1, seg, events


def plexon_to_floats(filename):
    """Convert the time series in a segment of a Plexon file into an array of floats.  Also returns the rate, start time, and end time for the data."""
    dir = os.path.dirname(filename)
    print("Calling PlexonIO('{}')".format(filename))
    r = neo.io.PlexonIO(filename)
    print(r)
    seg = r.read_segment()
    signals = seg.analogsignals  # nchan x nsamp - opposite of spike2.

    # this is one difference between neo neuralynx and neo plexon:
    nchannels = len(signals)
    nsamples = signals[0].shape[0]
    samprate = float(signals[0].sampling_rate.simplified)
    print('plexon_to_floats - Number of channels: {} ; Sample rate: {}; nsamples: {}'.format(nchannels, samprate, nsamples))

    out = np.zeros((nchannels, nsamples))
    for k in range(nchannels):
        out[k,:] = signals[k].transpose()

    # Plexon uses sample numbers!!!!!
    print('plexon data interval: [ {}, {} ]'.format(signals[0].t_start, signals[0].t_stop))
    # t0 = int(signals[0].t_start.simplified * samprate)
    # t1 = int(signals[0].t_stop.simplified * samprate)
    t0 = int(signals[0].t_start.simplified)
    t1 = int(signals[0].t_stop.simplified)
    print('plexon sample num interval: [ {}, {} ]'.format(t0, t1))

    # Not good enough - misses start / end recording events:
    # events = r.get_event_timestamps()
    events = []
    for x in seg.events:
        times = x.times
        labels = x.labels
        for i in range(len(times)):
            ts = times[i]
            label = labels[i]
            events.append((int(float(ts)*samprate), label))
    print(events)
    return out, samprate, t0, t1, seg, events

