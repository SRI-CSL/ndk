"""Module reads and writes header and data for EDF data."""
import numpy as np
from uritools import urisplit, uriunsplit, urijoin
from logging import getLogger
from .dsfilter import *

lg = getLogger(__name__)

#Poor man's version of:
#'https://github.com/breuderink/eegtools/blob/master/eegtools/io/edfplus.py'

#Values are slightly different from those computed by FieldTrip, however they
#are identical to those computed by Biosig and EDFBrowser. The difference is due
#to the calibration.

#Brazenly lifted from phypno by C. Connolly because of module incompatibilities. 

#Author sequence: 
#Boris Reuderink (http://nl.linkedin.com/in/breuderink/)
#Gio Piantoni (https://github.com/gpiantoni)

from datetime import datetime, timedelta
from math import floor
from re import findall
from struct import pack

from numpy import empty, asarray, fromstring, iinfo, abs, max

EDF_FORMAT = 'int16'  # by definition
edf_iinfo = iinfo(EDF_FORMAT)
DIGITAL_MAX = edf_iinfo.max
DIGITAL_MIN = -1 * edf_iinfo.max  # so that digital 0 = physical 0


def _assert_all_the_same(items):
    """Check that all the items in a list are the same.

    """
    assert all(items[0] == x for x in items)


class Edf:
    """Provide class EDF, which can be used to read the header and the data."""

#    Parameters
#    ----------
#    edffile : str
#        Full path for the EDF file

#    Attributes
#    ----------
#    h : dict
#        disorganized header information

    def __init__(self, uri, use_file=None):
        print(uri)
        if isinstance(uri, (str, unicode)):
            check = uri.split(':')
            if check[0] == 'edf':
                self.desc = urisplit(uri)
            else:
                self.desc = urisplit('edf://localhost' + uri)
                self.desc = urisplit(uri)
        else:
            print('Non-string uri: ', uri)
            self.desc = uri
        self.filename = self.desc.path
        print('# Filename: {}'.format(self.filename))
        self._read_hdr()

    def _read_hdr(self):
        """Read header from EDF file.  It only reads the header for internal
        purposes and adds a hdr."""

        with open(self.filename, 'rb') as f:

            hdr = {}
            assert f.tell() == 0
            assert f.read(8) == b'0       '

            # recording info
            hdr['subject_id'] = f.read(80).decode('utf-8').strip()
            hdr['recording_id'] = f.read(80).decode('utf-8').strip()

            # parse timestamp
            (day, month, year) = [int(x) for x in findall('(\d+)',
                                  f.read(8).decode('utf-8'))]
            (hour, minute, sec) = [int(x) for x in findall('(\d+)',
                                   f.read(8).decode('utf-8'))]
            hdr['start_time'] = datetime(year + 2000, month, day, hour, minute,
                                         sec)

            # misc
            hdr['header_n_bytes'] = int(f.read(8))
            f.seek(44, 1)  # reserved for EDF+
            hdr['n_records'] = int(f.read(8))
            hdr['record_length'] = float(f.read(8))  # in seconds
            nchannels = hdr['n_channels'] = int(f.read(4))

            # read channel info
            channels = range(hdr['n_channels'])
            hdr['label'] = [f.read(16).decode('utf-8').strip() for n in
                            channels]
            hdr['transducer'] = [f.read(80).decode('utf-8').strip()
                                 for n in channels]
            hdr['physical_dim'] = [f.read(8).decode('utf-8').strip() for n in
                                   channels]
            hdr['physical_min'] = asarray([float(f.read(8))
                                           for n in channels])
            hdr['physical_max'] = asarray([float(f.read(8))
                                           for n in channels])
            hdr['digital_min'] = asarray([float(f.read(8)) for n in channels])
            hdr['digital_max'] = asarray([float(f.read(8)) for n in channels])
            hdr['prefiltering'] = [f.read(80).decode('utf-8').strip()
                                   for n in channels]
            hdr['n_samples_per_record'] = [int(f.read(8)) for n in channels]
            f.seek(32 * nchannels, 1)  # reserved

            assert f.tell() == hdr['header_n_bytes']

            self.hdr = hdr

    def return_hdr(self):
        """Return the header for further use."""

#        Returns
#        -------
#        subj_id : str
#            subject identification code
#        start_time : datetime
#            start time of the dataset
#        s_freq : float
#            sampling frequency
#        chan_name : list of str
#            list of all the channels
#        n_samples : int
#            number of samples in the dataset
#        orig : dict
#            additional information taken directly from the header
#


        subj_id = self.hdr['subject_id']
        start_time = self.hdr['start_time']
        _assert_all_the_same(self.hdr['n_samples_per_record'])
        s_freq = (self.hdr['n_samples_per_record'][0] /
                  self.hdr['record_length'])
        chan_name = self.hdr['label']
        n_samples = (self.hdr['n_samples_per_record'][0] *
                     self.hdr['n_records'])

        return subj_id, start_time, s_freq, chan_name, n_samples, self.hdr

    def _read_dat(self, i_chan, begsam, endsam):
        """Read raw data from a single EDF channel.  Reads only one channel at
the time. Very initial implementation, very simple and
probably not very fast"""

#        Parameters
#        ----------
#        i_chan : int
#            index of the channel to read
#        begsam : int
#            index of the first sample
#        endsam : int
#            index of the last sample
#
#        Returns
#        -------
#        numpy.ndarray
#            A vector with the data as written on file, in 16-bit precision

        assert begsam < endsam

        # print(i_chan)
        begsam = float(begsam)
        endsam = float(endsam)

        n_sam_rec = self.hdr['n_samples_per_record']

        begrec = int(floor(begsam / n_sam_rec[i_chan]))
        begsam_rec = int(begsam % n_sam_rec[i_chan])

        endrec = int(floor(endsam / n_sam_rec[i_chan]))
        endsam_rec = int(endsam % n_sam_rec[i_chan])

        dat = empty(shape=(int(endsam) - int(begsam)), dtype='int16')
        i_dat = 0

        with open(self.filename, 'rb') as f:
            for rec in range(begrec, endrec + 1):
                if rec == begrec:
                    begpos_rec = begsam_rec
                else:
                    begpos_rec = 0

                if rec == endrec:
                    endpos_rec = endsam_rec
                else:
                    endpos_rec = n_sam_rec[i_chan]

                begpos = begpos_rec + sum(n_sam_rec) * rec + sum(
                    n_sam_rec[:i_chan])
                endpos = endpos_rec + sum(n_sam_rec) * rec + sum(
                    n_sam_rec[:i_chan])

                f.seek(begpos * 2 + self.hdr['header_n_bytes'])
                samples = f.read(2 * (endpos - begpos))

                i_dat_end = i_dat + endpos - begpos
                dat[i_dat:i_dat_end] = fromstring(samples, dtype='<i2')
                i_dat = i_dat_end

        return dat

    def return_dat(self, chan, begsam, endsam):
        """Read data from an EDF file.  Reads channel by channel, and adjusts the values by calibration."""
#        Parameters
#        ----------
#        chan : list of str
#            index (indices) of the channels to read
#        begsam : int
#            index of the first sample
#        endsam : int
#            index of the last sample

#        Returns
#        -------
#        numpy.ndarray
#            A 2d matrix, where the first dimension is the channels and the
#            second dimension are the samples.

        hdr = self.hdr
        dig_min = hdr['digital_min']
        phys_min = hdr['physical_min']
        phys_range = hdr['physical_max'] - hdr['physical_min']
        dig_range = hdr['digital_max'] - hdr['digital_min']
        assert all(phys_range > 0)
        assert all(dig_range > 0)
        gain = phys_range / dig_range

        dat = empty(shape=(len(chan), endsam - begsam), dtype='float64')

        for i, i_chan in enumerate(chan):
            d = self._read_dat(i, begsam, endsam).astype('float64')
            dat[i, :] = (d - dig_min[i]) * gain[i] + phys_min[i]

        return dat

    def return_markers(self):
        """"""
        return []


def write_edf(data, filename, physical_max=1000):
    """Export data to FieldTrip."""

#    Parameters
#    ----------
#    data : instance of ChanTime
#        data with only one trial
#    filename : path to file
#        file to export to (include '.mat')
#    physical_max : int
#        values above this parameter will be considered saturated (and also
#        those that are too negative). This parameter defines the precision.

#    Notes
#    -----
#    Data is always recorded as 2 Byte int (which is 'int16'), so precision is
#    limited. You can control the precision with physical_max. To get the
#    precision:

#    >>> precision = physical_max / DIGITAL_MAX

#    where DIGITAL_MAX is 32767.
    if data.start_time is None:
        raise ValueError('Data should contain a valid start_time (as datetime)')
    start_time = data.start_time + timedelta(seconds=data.axis['time'][0][0])

    if physical_max is None:
        physical_max = max(abs(data.data[0]))

    precision = physical_max / DIGITAL_MAX
    lg.info('Data exported to EDF will have precision ' + str(precision))

    physical_min = -1 * physical_max
    dat = data.data[0] / physical_max * DIGITAL_MAX
    dat = dat.astype(EDF_FORMAT)
    dat[dat > DIGITAL_MAX] = DIGITAL_MAX
    dat[dat < DIGITAL_MIN] = DIGITAL_MIN

    with open(filename, 'wb') as f:
        f.write('{:<8}'.format(0).encode('ascii'))
        f.write('{:<80}'.format('X X X X').encode('ascii'))  # subject_id
        f.write('{:<80}'.format('Startdate X X X X').encode('ascii'))
        f.write(start_time.strftime('%d.%m.%y').encode('ascii'))
        f.write(start_time.strftime('%H.%M.%S').encode('ascii'))

        n_smp = data.data[0].shape[1]
        s_freq = int(data.s_freq)
        n_records = n_smp // s_freq  # floor
        record_length = 1
        n_channels = data.number_of('chan')[0]

        header_n_bytes = 256 + 256 * n_channels
        f.write('{:<8d}'.format(header_n_bytes).encode('ascii'))
        f.write((' ' * 44).encode('ascii'))  # reserved for EDF+

        f.write('{:<8}'.format(n_records).encode('ascii'))
        f.write('{:<8d}'.format(record_length).encode('ascii'))
        f.write('{:<4}'.format(n_channels).encode('ascii'))

        for chan in data.axis['chan'][0]:
            f.write('{:<16}'.format(chan).encode('ascii'))  # label
        for _ in range(n_channels):
            f.write(('{:<80}').format('').encode('ascii'))  # tranducer
        for _ in range(n_channels):
            f.write('{:<8}'.format('uV').encode('ascii'))  # physical_dim
        for _ in range(n_channels):
            f.write('{:<8}'.format(physical_min).encode('ascii'))
        for _ in range(n_channels):
            f.write('{:<8}'.format(physical_max).encode('ascii'))
        for _ in range(n_channels):
            f.write('{:<8}'.format(DIGITAL_MIN).encode('ascii'))
        for _ in range(n_channels):
            f.write('{:<8}'.format(DIGITAL_MAX).encode('ascii'))
        for _ in range(n_channels):
            f.write('{:<80}'.format('').encode('ascii'))  # prefiltering
        for _ in range(n_channels):
            f.write('{:<8d}'.format(s_freq).encode('ascii'))  # n_smp in record
        for _ in range(n_channels):
            f.write((' ' * 32).encode('ascii'))

        l = s_freq * n_channels  # length of one record
        for i in range(n_records):
            i0 = i * s_freq
            i1 = i0 + s_freq
            x = dat[:, i0:i1].flatten(order='C')  # assumes it's ChanTimeData
            f.write(pack('<' + 'h' * l, *x))


#
# This case underscores the need for a better abstraction for data
# stores:
#
def edf_to_floats(filename):
    """Convert a multichannel EDF file into an array of floats.  Also returns the rate, start time, and end time for the data."""
    ds = Edf(filename, use_file=True)
    hdr = ds.return_hdr()
    samprate = hdr[2] / ds.hdr['record_length']
    channels = ds.hdr['label']
    nchan = len(channels)
    nsamp = ds.hdr['n_records']
    print('Channels: {}\nSample Rate: {}\nSample Count: {}\n'.format(nchan, samprate, nsamp))
    data = ds.return_dat(channels, 0, nsamp)
    return data, samprate, 0.0, nsamp / samprate, data


class edf_ds():
    """Primary class for handling data stores that use EDF files."""
    def __init__(self, filename):
        self.filename = filename
        self.edf = Edf(filename)
        self.hdr = self.edf.return_hdr()
        self.samprate = self.hdr[2] / self.edf.hdr['record_length']
        self.channels = self.edf.hdr['label']
        self.nchannels = len(self.channels)
        self.nsamp = self.edf.hdr['n_records']
        self.data = self.edf.return_dat(range(self.nchannels), 0, self.nsamp)
        self.start = 0.0
        self.end = self.nsamp
        self.dsf = dsfilter(self.samprate)

    def filter_data(self, low_cutoff, high_cutoff, order=5):
        self.dsf.filter_data(low_cutoff, high_cutoff, order)

    def get_signal(self, channel, from_time=None, to_time=None, into_array=None):
        if from_time == None:
            from_time = self.start
        if to_time == None:
            to_time = self.end
        if into_array != None:
            result = into_array
            result[:] = 0
            to_time = from_time + len(result)
            if to_time >= self.end:
                to_time = self.end - 1
            n = len(result)
        else:
            n = (to_time - from_time)
            result = numpy.zeros(n)
        result[:] = self.data[channel, from_time:to_time]
        return result


    def get_chunk(self, channel, start, end):
        x = self.data[channel, int(start):int(end)]
        r = np.asarray(x)
        return self.dsf.maybe_filter(r)
                    
    
    def get_dataset_interval(self):
        return self.start, self.end



def read_edf_file(filename):
    return edf_to_floats(filename)
