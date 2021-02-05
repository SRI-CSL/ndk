import wfdb
import numpy as np
from uritools import urisplit, uriunsplit, urijoin
import os
import pickle
from .dsfilter import *
from ndk.ui import iprint,wprint,eprint
import ndk.ds

NEW_NBM = True
NEW_FILTERING = True
#
# nbf = NDK Binary Format: This is a pyramid-style "native"
# representation for NDK.
#
# Version 1: Uses mmap to map a raw data file into memory and use that
# as the primary data structure.  Modern architectures will support
# large-ish recordings with ~500M of data.  If we need to scale up
# from there, we will need to get creative about adding windowing to
# this scheme.  For now, define the basic data structure and let the
# file itself be a raw memory map of the structure.


class nb_metadata:
    # Problem: this initialization will, by default, create a bunch of
    # absolute pathnames, which we don't want.
    """Primary class for handling NDK Binary Format data stores."""
    def __init__(self, dirname, sample_rate=32000.0, nchannels=1, start=0, end=0, relative=False):
        self.dirname = os.path.abspath(dirname)
        self.name = os.path.basename(self.dirname)
        self.filename = os.path.join(self.dirname, '{}.nbm'.format(self.name))
        self.relative = relative    # Dictates whether the nbm file has an absolute or relative path in it.
        # Don't hard-code these:
        #self.rawfile = os.path.join(dirname, '{}.raw'.format(name))
        #self.meanfile = os.path.join(dirname, '{}-mean.pyr'.format(name))
        #self.minfile = os.path.join(dirname, '{}-min.pyr'.format(name))
        #self.maxfile = os.path.join(dirname, '{}-max.pyr'.format(name))

        # Always in units of samples per second:
        self.samprate = sample_rate
        self.nchannels = nchannels
        # Always in units of samples:
        self.start = start
        self.end = end

    def rawfile(self):
        """Returns the absolute pathname for the raw data file."""
        return(os.path.join(self.dirname, '{}.raw'.format(self.name)))

    def meanfile(self):
        return(os.path.join(self.dirname, '{}-mean.pyr'.format(self.name)))

    def minfile(self):
        return(os.path.join(self.dirname, '{}-min.pyr'.format(self.name)))

    def maxfile(self):
        return(os.path.join(self.dirname, '{}-max.pyr'.format(self.name)))

               
    def npts(self):
        """Returns the number of samples in an NBF."""
        return self.end-self.start

    # Probably superfluous now:
    def basedir(self):
        return os.path.dirname(self.filename)

    # We should not need this anymore, at least in its current form:
    def maybe_fix_filenames(self, dirname):
        old = os.path.dirname(self.filename)
        iprint('maybe_fix_filenames: dirname={}'.format(dirname))
        iprint('maybe_fix_filenames: old dirname={}'.format(old))
        if dirname != old:
            wprint('Directory was moved.  Resetting...')
            name = os.path.basename(dirname)
            self.filename = os.path.join(dirname, '{}.nbm'.format(name))
            self.rawfile = os.path.join(dirname, '{}.raw'.format(name))
            self.meanfile = os.path.join(dirname, '{}-mean.pyr'.format(name))
            self.minfile = os.path.join(dirname, '{}-min.pyr'.format(name))
            self.maxfile = os.path.join(dirname, '{}-max.pyr'.format(name))

    # This still saves one absolute pathname.  Once we debug this,
    # maybe we can allow dirname to be '.' by default.
    def save(self, source=None):
        """Saves NBF metadata."""
        filename = self.filename
        if os.path.isfile(filename):
            i = 1
            while os.path.isfile("{}.{}".format(filename, i)):
                i += 1
            os.rename(filename, "{}.{}".format(filename, i))
        with open(filename, 'w') as f:
            p1 = os.path.dirname(self.filename)
            if self.relative:
                f.write("dirname=.\n")
            else:
                f.write("dirname={}\n".format(p1))
            f.write("sample_rate={}\n".format(self.samprate))
            f.write("nchannels={}\n".format(self.nchannels))
            f.write("start={}\n".format(self.start))
            f.write("end={}\n".format(self.end))
            # This is a hack that assumes that the nbm file is created
            # "near" the neo-compatible source data files.  We can
            # override this by providing a source argument:
            if source is None:
                f.write("source={}\n".format(p1))
            else:
                f.write("source={}\n".format(source))
                


    
def read_nbm(filename):
    """Reads NBF metadata, which is usually saved in an 'nbm' file.  This
file contains basic information about the recording."""
    attrs = {}
    iprint('Reading nbm file {}.'.format(filename))
    with open(filename) as f:
        for line in f:
            x = line[0:-1].split('=')   # No spaces allowed please!!
            if len(x) == 2:
                attrs[x[0]] = x[1]
                iprint('Attribute: {} = {}'.format(x[0],x[1]))
            else:
                wprint("Don't know what to do with this line: {}".format(line))
    if attrs['dirname'][0] == '.':
        new_dirname = os.path.dirname(os.path.abspath(filename))
        attrs['dirname'] = new_dirname
    metadata = nb_metadata(attrs['dirname'], 
                           sample_rate=float(attrs['sample_rate']),
                           nchannels=int(attrs['nchannels']), 
                           start=int(attrs['start']), 
                           end=int(attrs['end']))
    return metadata


        
# Note that 'data' must be an array of 2 dimensions: nchan x npts

def floats_to_nbf(dirname, data, sample_rate, start, end, events=None, relative=True, source=None, permute=None):
    """Converts an array of floats (data) into NBF form in the specified
dirname, with the sample rate, start and end times given by the
arguments.  If supplied, 'events' is a list of time,event pairs."""
    name = os.path.basename(dirname)
    data_shape = data.shape
    nchan = data_shape[0]
    npts = data_shape[1]
        
    md = nb_metadata(dirname, sample_rate, nchan, start, end, relative)
    filename = md.filename

    # If the directory doesn't exist, create it:
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    if events:
        eventfile = os.path.join(dirname, 'events.dat')
        print("Events found:")
        with open(eventfile, 'w') as f:
            for (ts, label) in events:
                print("{}: {}".format(label, ts))
                if label[-1] != '\n':
                    f.write('{}:{}\n'.format(ts, label))
                else:
                    f.write('{}:{}'.format(ts, label))

    md.save(source)

    # Initialize the raw data file:
    rawfile = md.rawfile()
    iprint('Mapping raw file {} using shape {}.'.format(rawfile, data_shape))
    raw = np.memmap(rawfile, np.dtype('f4'), 'w+', shape=data_shape)
    # print(raw)
    iprint('Initializing raw file: {}'.format(rawfile))
    if permute is None:
        rawchan = [k for k in range(nchan)]
    else:
        rawchan = permute

    for j in range(nchan):
        k = rawchan[j]
        raw[j][:] = data[k][:]
        iprint('Wrote {} data elements'.format(len(data[k])))
    # print(raw)
    del raw


# Other physiological data:
def wfdb_to_nbf(filename, to_dir, channel_names=['V1', 'V2', 'V3', 'V4']):
    iprint("Creating nbf from file: {}".format(filename))

    data,fields = wfdb.rdsamp(filename)

    print("Available fields: {}".format(fields))
    # end sample number:
    start = 0
    end = fields['sig_len']

    # I think this is the sample rate:
    sample_rate = fields['fs']

    # This will barf if the channel_names don't exist in the file, but
    # collect the indices of the channels that we are interested in
    # (the Vn leads on the chest, for now):
    indices = [fields['sig_name'].index(x) for x in channel_names]
    nchan = len(indices)
    
    md = nb_metadata(to_dir, sample_rate, nchan, start, end, True)

    nbf_filename = md.filename

    # If the directory doesn't exist, create it:
    if not os.path.isdir(to_dir):
        os.makedirs(to_dir)

    # See floats_to_nbf for hints on Events.  For now, we aren't
    # looking for events, but this is where we would create an
    # events.dat file.

    # Save the metadata file:
    md.save(source=filename)

    # Create the raw data file:
    rawfile = md.rawfile()

    iprint('Mapping raw file {} using shape {}.'.format(rawfile, ( nchan, end )))
    raw = np.memmap(rawfile, np.dtype('f4'), 'w+', shape=(nchan, end) )
    # print(raw)
    iprint('Initializing raw file: {}'.format(rawfile))

    for j in range(nchan):
        k = indices[j]
        raw[j][:] = data[:,k]
        iprint('Wrote {} data elements'.format(len(data[:,k])))
    # print(raw)
    del raw



    
def old(rawfile, data, nchan):
    with open(rawfile, 'wb') as f:
        for k in range(nchan):
            f.write(data[k])
            iprint('Wrote {} data elements'.format(len(data[k])))


def nbf_filename(uri, type):
    desc = urisplit(uri)
    fname = desc.path
    dirname = os.path.dirname(fname)
    basename = os.path.basename(fname)
    name = basename.split('.')[0]
    result = os.path.join(dirname, '{}.{}'.format(name, type))
    return result
    


def get_version_count(basename):
    l = os.listdir(os.path.dirname(basename))
    version = 0
    for name in l:
        comps = name.split('.')
        if len(comps) > 2 and comps[-2] == 'raw':
            this_version = int(comps[-1])
            if this_version > version:
                version = this_version
    return version


class nbf:
    """Class that holds information about a wfdb ECG recording (data and metadata for a cardiac recording)."""
    def __init__(self, filename, start=0, end=0, sample_rate=32000.0, nchannels=1, version=None):
        
        iprint("Creating nbf from file: {}".format(filename))

        self.filename = filename
        self.filtered_data = {}
        iprint("Set up filtered_data list: {}".format(self.filtered_data))
        self.lazy = False
    
        dirname = os.path.dirname(self.filename)
        self.dirname = dirname
        iprint('Using directory {}'.format(dirname))

        # These are defaults, can be overridden later:
        name = os.path.basename(self.filename)
        name = name.split('.')[0]
        new = False
        # If the directory doesn't exist, create it:
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
            
        # Compute the constituent filenames:
        self.mdfile = os.path.join(dirname, '{}.nbm'.format(name))
        efile = os.path.join(dirname, 'events.dat')
        def getfirst(y):
            return y[0]
        # If there are any events, generate the list:
        self.events = []
        if os.path.exists(efile):
            print("Events found:")
            with open(efile, 'r') as f:
                for line in f:
                    x = line.split(':')
                    elabel = x[1].split('\n')[0]
                    print("{}: {}".format(elabel, x[0]))
                    self.events.append( (int(x[0]), elabel) )
                self.events.sort(key=getfirst)
        # How?  Each channel is a "strip" of consecutive samples at a rate of
        # sample_rate * 2^(-i) samples per second.

        # Thus, the first index is the channel, the second is sample.

        # This will load metadata, but otherwise will create it with
        # defaults:
        iprint('Checking for metadata file {}'.format(self.mdfile))
        if not os.path.exists(self.mdfile):
            dirname = os.path.dirname(os.path.abspath(self.mdfile))
            self.md = nb_metadata(dirname, sample_rate=sample_rate, start=start, end=end, nchannels=nchannels)
        elif NEW_NBM:
            self.md = read_nbm(self.mdfile)
        else:
            with open(self.mdfile, 'r') as f:
                self.md = pickle.load(f)
                self.md.maybe_fix_filenames(dirname)
                
        self.dsf = dsfilter(self.md.samprate)
        self.data_shape = (self.md.nchannels, self.md.npts())
        iprint('Initializing nbf object with data shape {}.'.format(self.data_shape))
        self.t0 = self.md.start
        self.t1 = self.md.end
        self.version = None
        self.map_raw_data(version)
        # print(self.raw)


    def map_raw_data(self, version):
        """Maps raw data files into memory for direct access to data."""
        rawfile = self.md.rawfile()
        vmax = get_version_count(rawfile)
        iprint('Versions found: {}'.format(vmax))
        if version is not None and version <= vmax:
            self.version = version
            iprint('Selecting version {}.'.format(self.version))
            rawfile = rawfile + '.{}'.format(version)
            
        iprint('Checking for raw file: {}'.format(rawfile))
        
        if not os.path.exists(rawfile):
            iprint('Creating raw file...')
            arr = np.zeros(self.data_shape, np.dtype('f4'))
            with open(rawfile, 'wb') as f:
                f.write(arr)
        
        iprint('Mapping raw file with shape {}...'.format(self.data_shape))
        self.raw = np.memmap(rawfile, np.dtype('f4'), 'r+', shape=self.data_shape)


    def make_filtered_data(self, force=False):
        """Makes a cache of the filtered data, with filter parameters specified by the dsf object."""
        dsf = self.dsf
        dir = self.dirname
        fdatafile = os.path.join(dir, 'filtered_{}_{}.raw'.format(dsf.lo, dsf.hi) )
        if force or not os.path.exists(fdatafile):
            iprint('Creating raw file {} using shape {}.'.format(fdatafile, self.data_shape))
            fdata = np.memmap(fdatafile, np.dtype('f4'), 'w+', shape=self.data_shape)
            iprint('Initializing filtered data file: {}'.format(fdatafile))
            for k in range(self.data_shape[0]):
                data = dsf.maybe_filter(self.raw[k])
                fdata[k,:] = data[:]
            # print(raw)
            del fdata


    def map_filtered_data(self):
        """Maps filtered data files into memory for direct access to filtered data."""
        dsf = self.dsf
        dir = self.dirname
        if len(self.filtered_data) == 0:
            fdatafile = os.path.join(dir, 'filtered_{}_{}.raw'.format(dsf.lo, dsf.hi) )
            iprint('Checking for filtered data file: {}'.format(fdatafile))
        
            if not os.path.exists(fdatafile):
                self.make_filtered_data()
        
            iprint('Mapping filtered data file with shape {}...'.format(self.data_shape))
            self.filtered_data = np.memmap(fdatafile, np.dtype('f4'), 'r+', shape=self.data_shape)




    def get_filtered_data(self, channel):
        if not self.dsf.do_filter:
            return self.raw[channel]
        else:
            self.map_filtered_data()
            return self.filtered_data[channel]
            
        
    def get_dataset_interval(self):
        """Returns the time interval corresponding to this dataset, as sample numbers."""
        return self.t0, self.t1
    
    def filter_data(self, low_cutoff, high_cutoff, order=5):
        """Sets the desired filter parameters for accessing the raw data."""
        self.dsf.filter_data(low_cutoff, high_cutoff, order)
            
    def get_signal(self, channel):
        """Gets the entire signal for a given channel."""
        return self.raw[channel]

    # Rethink this - we might want to filter the WHOLE dataset and
    # just extract from that.  Edge effects seem to creep in when we
    # ask for chunks.
    def get_chunk(self, channel, start, end):
        """Returns the portion of the time series for the specified channel that is between start and end inclusive (time)."""
        i0 = max(0, int(start - self.t0))
        i1 = max(0, int(end - self.t0))

        # We should probably call numpy.asarray to coerce this to a pure array:
        # r = np.asarray(x[i0:i1])
        if self.lazy:
            x = self.raw[channel]
            r = x[i0:i1].reshape(i1-i0)
            result = self.dsf.maybe_filter(r)
        elif NEW_FILTERING:
            x = self.get_filtered_data(channel)
        else:
            # If we miss with filtered_data, run the filter and add it
            # to the filtered_data dictionary.  This is where we could
            # try filesystem caching:
            try:
                x = self.filtered_data[channel]
            except KeyError:
                x = self.raw[channel]
                iprint("In get_chunk({},{},{},{}): Initializing filtered copy".format(self, channel, start, end))
                #self.get_filtered_data(x, channel)
                self.filtered_data[channel] = self.dsf.maybe_filter(x)
                x = self.filtered_data[channel]
        try:
            result = x[i0:i1].reshape(i1-i0)
        except:
            # Most likely, we reached an edge, so just return a 0 vector:
            result = np.zeros((i1-i0))
        # result = result.reshape(len(r))
                
        return result

    def get_chunk1(self, channel, start, end):
        """Returns the portion of the time series for the specified channel that is between start and end inclusive (time)."""
        print("Please don't use get_chunk1!")
        x = self.raw[channel]
        i0 = max(0, int(start - self.t0))
        i1 = max(0, int(end - self.t0))
        di = i1-i0

        # We should probably call numpy.asarray to coerce this to a pure array:
        r = np.ndarray(buffer=x, offset=i0, dtype=np.dtype('f4'), shape=(di))
        result = self.dsf.maybe_filter(r)
        # result = result.reshape(len(r))
        return result
    

    def save_metadata(self):
        self.md.save()
            

#    def close(self):
#        np.memunmap
