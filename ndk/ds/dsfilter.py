"""Module supplies filtering to data source classes."""
import numpy as np
from ndk.features import butter_low, butter_high, butter_bandpass
from scipy.signal import convolve, butter, lfilter, filtfilt
from ndk.ui import iprint,wprint,eprint

#
# This class implements bandpass filtering for data.  It could be made
# cleaner.  Nonetheless, it's an attempt to isolate the filtering
# operations that we may or may not want to perform during raw data
# delivery.  The desire for filtering has to be countered by the fact
# that the basis set representation that we use for events is itself
# effectively a bandpass filter, so we should probably avoid filtering
# as much as possible.
#

class dsfilter():
    def __init__(self, samprate, lazy=True):
        self.samprate = samprate
        self.do_filter = False
        self.butter_info = None
        self.order = 5
        self.hi = 0.0
        self.lo = 0.0
        self.lazy = lazy
        self.filtered_data = None

    # axis=0 appears to be essential here!
    def maybe_filter(self, r):
        if self.do_filter and len(r) > 0:
            b, a = self.butter_info
            r = filtfilt(b, a, r, padlen=0, axis=0)
        return r

    # This should be part of a base class that can filter any chunk of
    # a dense time series:
    def filter_data(self, low_cutoff, high_cutoff, order=6):
        # Only sets the filter parameters.  Does not actually filter data yet.
        self.lo = low_cutoff
        self.hi = high_cutoff
        self.order = order

        # nyquist = 0.5 * self.samprate
        iprint('Nyquist freq. is {}.'.format(self.samprate/2))
        if self.lo is None:
            if self.hi is None:
                self.do_filter = False
                iprint('No filtering.')
            else:
                f = 2 * self.hi / self.samprate
                self.butter_info = butter(self.order, f, btype='lowpass')
                iprint('Low-pass filtering (f<{} Hz).'.format(self.hi))
                self.do_filter = True
        else:
            if self.hi is None:
                f = 2 * self.lo / self.samprate
                self.butter_info = butter(self.order, f, btype='highpass')
                iprint('High-pass filtering (f>{} Hz).'.format(self.lo))
                self.do_filter = True
            else:
                f0 = 2 * self.lo / self.samprate
                f1 = 2 * self.hi / self.samprate
                self.butter_info = butter(self.order, (f0,f1), btype='bandpass')
                iprint('Bandpass filtering [{}, {}].'.format(self.lo, self.hi))
                self.do_filter = True
                
