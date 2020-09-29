import mmap
import numpy as np

class mmap_meta:
    def __init__(self, filename, npts=0, nchan=1, samprate=32000.0):
        self.sample_rate = samprate
        self.nchannels = nchan
        self.npts = npts
        self.filename = filename


class mmap_ds():
    """Primary class for direct I/O to temporary memory-mapped float arrays files."""
    def __init__(self, filename, npts=0, nchan=1, samprate=32000.0):
        # Caller MUST provide npts and nchan.  Use put_signal or
        # put_chunk to add data, and set_rate / get_rate for rate
        # info.
        self.do_filter = False
        self.butter_info = None
        self.order = 5
        self.hi = 0.0
        self.lo = 0.0
        self.dsf = dsfilter(samprate)

        # r+ is read-write:
        data_shape = (nchan, npts+1)
        self.data = np.memmap(filename, np.dtype('f4'), 'r+', shape=data_shape)
        # Super-simple format: The initial element is the data rate.
        # The size of the file dictates the number of elements.
        self.rate = self.data[0,0]

    # This should be part of a base class that can filter any chunk of
    # a dense time series:
    def filter_data(self, low_cutoff, high_cutoff, order=5):
        self.dsf.filter_data(low_cutoff, high_cutoff, order)


    def set_rate(self,rate):
        self.data[0,0] = rate
        self.rate = rate

    def get_rate(self):
        self.rate = self.data[0,0]
        return self.rate

    #
    def write_data(self, filename):
        self.data.flush()

    def read_data(self, filename):
        return self.data


    
    # It is always meaningful to call get_signal:
    def get_signal(self, channel, from_time=None, to_time=None):
        if from_time == None:
            from_time = self.start
        if to_time == None:
            to_time = self.end
        # Should this be to_time or to_time-1 ??
        return self.data[channel,from_time:to_time]

    def put_signal(self, data, channel, from_time=None, to_time=None):
        if from_time == None:
            from_time = 0
        if to_time == None:
            to_time = len(data)
        # Should this be to_time or to_time-1 ??
        self.data[channel,from_time:to_time] = data[:]

    def get_chunk(self, channel, start, end):
        """Returns the portion of the time series for the specified channel that is between start and end inclusive (time)."""
        x = self.data[channel,:]
        i0 = max(0, int(start - self.t0))
        i1 = max(0, int(end+1 - self.t0))
        # We should probably call numpy.asarray to coerce this to a pure array:
        r = np.asarray(x[i0:i1])
        if self.do_filter:
            b, a = self.butter_info
            r = filtfilt(b, a, r, padlen=0)
        return r

