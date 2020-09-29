import wave
import numpy as np

class wav_ds():
    """Primary class for handling data stores that use WAV files."""
    def __init__(self, filename, data=None, rate=32000.0):
        self.filename = filename
        self.nchannels = 1
        # This data source has two possibilities: If 'data' is
        # provided, we assume that the caller wants to save the data
        # in the named file.  Otherwise, we interpret the file as a
        # WAV file with data.
        if data is None:
            self.direction = 'input'
            self.data, self.rate, self.start, self.end, junk = wav_to_floats(filename)
        else:
            self.direction = 'output'
            self.data = data
            self.rate = rate
            self.start = 0
            self.end = len(data)

        self.t0 = self.start
        self.t1 = self.end
        print('# Dataset interval: [{}, {}]'.format(self.t0, self.t1))
        self.samprate = self.rate
        self.do_filter = False
        self.butter_info = None
        self.order = 5
        self.hi = 0.0
        self.lo = 0.0

    #
    def write_data(self, filename):
        """Save the data points into a WAV file in the given filename."""
        s = wave.open(filename,'w')
        rate = s.setframerate(self.rate)
        # For wav files, we assume 16-bit signed.
        s.setsampwidth(2)
        s.setnchannels(1)
        out = np.zeros(len(self.data), dtype='short')
        xmin = self.data.min()
        xmax = self.data.max()
        out[:] = (32767.0 * self.data[:]) / (xmax-xmin)
        s.writeframes(out)
        s.close()
        return True

        
    # It is always meaningful to call get_signal:
    def get_signal(self, channel, from_time=None, to_time=None):
        if from_time == None:
            from_time = self.start
        if to_time == None:
            to_time = self.end
        # Should this be to_time or to_time-1 ??
        return self.data[from_time:to_time]


    # In theory, it is always meaningful to call put_signal, but if
    # direction is 'input', the result is not persistent, since it
    # will not be written (yet).  
    def put_signal(self, channel, data, from_time=None, to_time=None):
        if from_time == None:
            from_time = self.start
        if to_time == None:
            to_time = self.end

        n = to_time - from_time
        if n > len(data):
            n = len(data)
        
        # Should this be to_time or to_time-1 ??
        self.data[from_time:to_time] = data[0:n]




    def get_chunk(self, channel, start, end):
        """Returns the portion of the time series for the specified channel that is between start and end inclusive (time)."""
        x = self.data[channel]
        i0 = max(0, int(start - self.t0))
        i1 = max(0, int(end - self.t0))
        # We should probably call numpy.asarray to coerce this to a pure array:
        r = np.asarray(x[i0:i1])
        if self.do_filter:
            b, a = self.butter_info
            r = filtfilt(b, a, r, padlen=0)
        return r

        

    # This should be part of a base class that can filter any chunk of
    # a dense time series:
    def filter_data(self, low_cutoff, high_cutoff, order=5):
        self.lo = low_cutoff
        self.hi = high_cutoff
        self.order = order

        nyquist = 0.5 * self.samprate
        print('# Nyquist freq. is {}.'.format(nyquist))
        if self.lo is None:
            if self.hi is None:
                self.do_filter = False
                print('# No filtering.')
            else:
                f = self.hi / nyquist
                self.butter_info = butter(self.order, f, btype='low')
                print('# Low-pass filtering (f<{} Hz).'.format(self.hi))
                self.do_filter = True
        else:
            if self.hi is None:
                f = self.lo / nyquist
                self.butter_info = butter(self.order, f, btype='high')
                print('# High-pass filtering (f>{} Hz).'.format(self.lo))
                self.do_filter = True
            else:
                f0 = self.lo / nyquist
                f1 = self.hi / nyquist
                self.butter_info = butter(self.order, [f0,f1], btype='band')
                print('# Bandpass filtering [{}, {}].'.format(self.lo, self.hi))
                self.do_filter = True


# At present, we can only handle mono WAV files.
def wav_to_floats(filename):
    """Convert the WAV file in the given filename into an array of floats.  Also returns the rate, start time, and end time for the data."""
    s = wave.open(filename,'r')
    rate = s.getframerate()
    w = s.getsampwidth()
    strsig = s.readframes(s.getnframes())
    t0 = 0
    t1 = 0
    result = []
    if w == 2:
        x = np.fromstring(strsig, np.short)
        xmin = float(x.min())
        xmax = float(x.max())
        print( xmin, xmax )
        y = x / (xmax-xmin)
        data = y.astype('float32')
        out = np.reshape(data, (1,len(y)))
        t1 = len(y)
    else:
        print( "don't know how to interpret sampwidth of {}".format(w) )
    s.close()
    return out, rate, t0, t1, []
