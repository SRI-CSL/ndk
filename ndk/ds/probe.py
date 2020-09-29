import os
import os.path
import ndk.ds
from uritools import urisplit, uriunsplit
from ndk.ui import iprint,wprint,eprint

def parse_uri(uri, filename=None):
    true_uri, rel = resolve_relative_paths(uri, filename)
    if rel and filename is None:
        wprint('URI path is relative, but no base filename is specified.')
    else:
        uri = true_uri
        iprint('URI: {}'.format(uri))
    desc = urisplit(uri)
    return desc, rel

#
# ONLY works with filenames, not URIs:
#
def read_raw_data(filename):
    """Read the time series (LFP) data in the given filename and return an array of floats."""
    suffix = filename.split('.')[-1]
    if (suffix == 'wav'):
        return ndk.ds.wav_to_floats(filename)
    elif (suffix == 'smr'):
        return ndk.ds.spike2_to_floats(filename)
    elif (suffix == 'plx'):
        return ndk.ds.plexon_to_floats(filename)
    elif (suffix == 'ncs'):
        return ndk.ds.neuralynx_to_floats(filename)
    elif (suffix == 'edf'):
        return ndk.ds.read_edf_file(filename)
    elif (suffix == 'ds'):
        return ndk.ds.memmap(filename)
    elif (suffix == 'nbm'):
        nbf = ndk.ds.nbf(filename)
        #nbf = ndk.ds.open(filename)
        return nbf.raw, nbf.md.samprate, nbf.md.start, nbf.md.end, nbf, nbf.events
    else:
        print( "Can't interpret {} as a WAV, SMR, ncs, plx, EDF, or DS file.".format(filename) )


def find_ndk_data():
    """Try to find and return the search path for NDK time series data."""
    try:
        datapath = os.environ['NDKDATA']
        dirlist = datapath.split(':')
    except KeyError:
        dirlist = []
    return dirlist

#
# Only a step.  Ultimately, we will want to devise a URI that works
# for Cassandra as well as local data files:
#
def fullname(file):
    """Return the full (absolute) pathname of the given file."""
    # If the filename is absolute and exists, then return it:
    if (os.path.isabs(file) and os.path.isfile(file)):
        return file
    else:
        # If the NDKDATA environment variable is set, assume this is a
        # colon-separated list of directories that can contain raw
        # data files in some neo-friendly format.
        dirlist = find_ndk_data()
        for dir in dirlist:
            x = os.path.join(dir, file)
            if os.path.isfile(x):
                return x
        return file


def resolve_relative_paths(uri, filename):
    r = uri.split('/')
    uri2 = uri
    has_rel_path = False

    if len(r) > 3:
        has_rel_path = (r[3] == '.')

    if has_rel_path and filename is not None:
        dir = os.path.dirname(os.path.abspath(filename))
        rnew = r[0:3] + [ dir[1:] ] + r[4:]
        uri2 = '/'.join(rnew)
    return uri2, has_rel_path


#
# Resource descriptors look like Universal Resource Indicators (URIs).
#
# Grammar:
#
# <uri>      ::= <scheme>://<host>/<pathname>
#
# <scheme>   ::= pre | smr | cass | wav | ds | edf | nbm
# <host>     ::= Any legal IP hostname, including 'localhost'
# <pathname> ::= <rel_path> | <abs_path>
# <rel_path> ::= './'<path>
# <abs_path> ::= Any legal pathname not beginning with './'
# <path>     ::= Any legal pathname delimited by '/'
#
# A bit loose - to be clear, URI pathnames are nearly identical to URL
# paths, except for the case where the path begins with './'.  The
# latter is reserved for paths in an implied filesystem that are
# relative with respect to the location of the event store.  For
# example, if the event store is accessed as a file, and it's metadata
# table refers to a URI like "nbm://localhost/./uplps09/uplps09.nbm",
# then the '.' is replaced by the os.path.dirname for the event store
# to obtain the path to the data store.  This is primarily intended to
# allow event-stores to be arbitrarily "re-rooted" with their data
# stores.
#  
def open(uri, filename=None):
    """Given a URI for a data source, open it and return the appropriate data source object."""
    if True:
        desc, rel = parse_uri(uri, filename)
    else:
        new_uri, rel = resolve_relative_paths(uri, filename)
        if rel and filename is None:
            wprint('Data store path is relative, but no event store is specified.')
        else:
            uri = new_uri
        desc = urisplit(uri)

    iprint('Opening data store at {}'.format(uriunsplit(desc)))
    s = desc.scheme
    if s is None:
        return ndk.ds.neo_in.spike2(desc)
    elif s == 'pre':
        return ndk.ds.pre.pre_ds(desc)
    elif s == 'smr':
        return ndk.ds.neo_in.spike2(desc)
    # elif s == 'file':
    #    return ndk.ds.neo.neo_in(desc)
    elif s == 'cass':
        return ndk.ds.cass.cdb(desc)
    elif s == 'wav':
        return ndk.ds.wav.wav_ds(desc.path)
    elif s == 'ds':
        return ndk.ds.mmap.mmap_ds(desc.path)
    elif s == 'edf':
        return ndk.ds.edf_ds(desc)
    elif s == 'nbm':
        iprint('NBF path: {}'.format(desc.path))
        if rel and desc.path[0] == '/' and desc.path[1] == '.':
            return ndk.ds.nbf(desc.path[1:]) # Hack!
        else:
            return ndk.ds.nbf(desc.path)
    else:
        print("Don't know what do with this URI scheme: {} (in {})".format(s, uri))


# Parses pathnames that Harold has recorded.  Create a list of
# attributes from the filename.
def parse_path(pathname):
    parts = pathname.split('/')
    group = parts[-2]
    info = parts[-1]
    md = info.split('_')
    date = md[0]
    time = md[1]
    desc = ' '.join(md[2:])
    return [group, date, time, desc]

    
