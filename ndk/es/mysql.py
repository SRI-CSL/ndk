from __future__ import print_function
#
# sqlite3 submodule for ndk - an sqlite3-friendly way to store spike metadata.  TBD...
#

import sys
import numpy
import os.path
import pymysql as s3
#import sqlite3 as s3
import ndk.features
import ndk.ds
import math
from ndk.es.schema import *

#
# This is a template for "event level" data store objects, i.e.,
# parametric compression of features in the raw data stream using a
# database.
#
# This class should be decomposed into a parent class that is
# "SQL-like", so that the data model and queries can be exposed and
# modified there.  Although we use sqlite3 now, we might opt for MySQL
# or some other DB substrate in the future.
#
class mysqldb():
    """Primary class for handling event stores via sqlite3."""
    def __init__(self, db, host='localhost', port=3306):
        self.schema = schema()
        self.conn = s3.connect(host='localhost', port=3306, user='root', passwd='', db=db)
        self.cur = self.conn.cursor()
        self.sample_rate = 0
        self.nspikes = 0
        self.nchannels = 0
        self.recording = ''
        self.bases = {}  # Cache the basis functions as needed
        # Waveform min / max - these next three are examples of
        # properties that should be cached within the database itself,
        # somehow:
        self.modified = False
        self.min = None
        self.max = None
        try:
            self.cur.execute('select * from meta;')
            result = self.cur.fetchone()
            self.sample_rate = result[0]
            self.nchannels = result[1]
            self.recording = result[2]
        except:
            print ( "Uninitialized db." )

    # Very low-level - this simply reflects the SQL to the caller.  We
    # need to abstract this up a bit for a generic data store:
    def execute(self, string):
        return self.cur.execute(string)

    def fetchone(self):
        return self.cur.fetchone()

    def fetchall(self):
        return self.cur.fetchall()

    def waveform_min_max(self, channel):
        """Returns the min and max values for the waveform on the specified channel."""
        if self.min == None:
            self.cur.execute('select max(spikeID) from spiketimes;')
            x = self.cur.fetchone()
            if x != None:
                end_id = x[0]

            w = self.get_spike_waveform(0, channel)
            self.min = w.min()
            self.max = w.max()
            canFlush = True
            try:
                sys.stdout.flush()
            except:
                canFlush = False

            print ( "recomputing min,max: ", end="")
            for id in range(end_id):
                w = self.get_spike_waveform(id+1, channel)
                if id % 1000 == 0:
                    print ( id, end="" )
                    if canFlush:
                        sys.stdout.flush()
                if w.min() < self.min:
                    self.min = w.min()
                if w.max() > self.max:
                    self.max = w.max()
        return self.min,self.max

    def get_spike_amplitudes(self, channel):
        """Returns an array of relative spike amplitudes on the specified channel."""
        self.cur.execute('select spikeID from spiketimes order by samplenum;')
        r = self.cur.fetchall()
        amp = numpy.zeros(len(r))
        for k in range(len(r)):
            w = self.get_spike_waveform(r[k][0], channel)
            v1 = w.min()
            v2 = w.max()
            amp[k] = v2 - 0.5
        return amp

    def get_dataset_interval(self):
        """Returns the start and end times, in samples, for this dataset."""
        self.cur.execute('select min(samplenum) from spiketimes;')
        t0 = self.cur.fetchone()
        if t0 == None:
            return None
        self.cur.execute('select max(samplenum) from spiketimes;')
        t1 = self.cur.fetchone()
        return t0[0], t1[0]

    def get_spike_count(self):
        """Returns the number of spike events in this dataset."""
        self.cur.execute('select count(*) from spiketimes;')
        r = self.cur.fetchone()
        return r[0]

    # If we haven't computed coefficients yet, then this is 0:
    def get_waveform_count(self):
        """Returns the number of waveforms in this dataset.  Mainly for internal use."""
        self.cur.execute('select count(*) from spikecoefs where cindex == 0 and channel == 0;')
        r = self.cur.fetchone()
        return r[0]



    # 
    def get_spike_waveform(self, spike_id, channel):
        """Returns the waveform associated with the spike_id on the specified channel as a vector of values."""
        # Basis ID is really an index:
        k = self.get_spike_basisid(spike_id)
        try:
            # If we already computed the functions, use them:
            seq = self.bases[k]
        except:
            # We did not compute this basis yet, so get the sequence and
            # save it:
            seq = self.get_basis(k)
            self.bases[k] = seq
        a = self.get_spike_coefs(spike_id, channel)
        x = ndk.features.gen_waveform(a, seq)
        return x



    # 'width' used to be a discrete sample count.  However, to
    # accomodate different sampling rates, it needs to be a
    # floating-point time in milliseconds.  Window width in samples is
    # implied.
    def add_basis(self, basis_id, sigma, width, nfuncs=5, name="gaussian_basis"):
        """Adds a Gaussian basis with the specified sigma and width in
milliseconds to the database, with the desired ID."""
        ret = True
        self.cur.execute('SELECT * from basis where basisID == {};'.format(basis_id))
        result = self.cur.fetchall()
        print ( result )
        if result == None or len(result) == 0:
            self.cur.execute('INSERT INTO basis values ({}, {}, {}, {}, "{}");'.format(basis_id, sigma, width, nfuncs, name))
            self.conn.commit()
        else:
            result = result[0]
            if math.fabs(result[1] - sigma) > 1.0e-06 or result[2] != width or result[3] != nfuncs:
                print ( "Basis ID {} is already in use with these parameters: sigma={}, window width={}, n={}".format(result[0], result[1], result[2], result[3]) )
                print ( "Here are the basis ID's in use in this file:" )
                self.cur.execute('SELECT * from basis;')
                result = self.cur.fetchall()
                print ( "BasisID  sigma   window   nfuncs   name" )
                for tuple in result:
                    print ( "{}        {}    {}     {}      {}".format(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4]) )
                ret = False
            else:
                print ( "Your basis set is already available in {}.".format(self) )
                ret = True
        return ret


    def get_basis(self, basis_id):
        """Returns the set of basis functions corresponding to 'basis_id'."""
        try:
            seq = self.bases[basis_id]
        except:
            self.cur.execute('SELECT * from basis where basisID == {};'.format(basis_id))
            result = self.cur.fetchall()
            if len(result) < 1:
                return None
            rate = self.sample_rate
            # Convert to seconds, then to sample numbers:
            sigma = result[0][1] * 0.001 * rate
            width = int(result[0][2] * 0.001 * rate)
            nfuncs = int(result[0][3])
            print ( "Window width = ", width )
            gen_fn = self.schema.basis_generator(result[0][4])
            seq = gen_fn(sigma, width, nfuncs)
            self.bases[basis_id] = seq
        return seq



    def set_spike_time(self, spike_id, spike_timestamp):
        """Forces the specified spike_id to have the given spike_timestamp (in samples)."""
        self.cur.execute('INSERT or REPLACE INTO spiketimes VALUES({}, {});'.format(spike_id, spike_timestamp))

    def get_spike_time(self, spike_id):
        """Returns the timestamp (sample number) associated with the given spike_id."""
        self.cur.execute('SELECT samplenum FROM spiketimes where spikeID == {};'.format(spike_id))
        result = self.cur.fetchone()
        if result != None:
            return result[0]
        else:
            return result


    def set_spike_basisid(self, spike_id, basis_id):
        """Sets the basis_id associated with the given spike_id."""
        self.cur.execute('INSERT or REPLACE INTO spikebasis VALUES({}, {});'.format(spike_id, basis_id))

    def get_spike_basisid(self, spike_id):
        """Returns the basis_id associated with the given spike_id."""
        self.cur.execute('SELECT basisID FROM spikebasis where spikeID == {};'.format(spike_id))
        result = self.cur.fetchone()
        return result[0]


    def set_spike_label(self, spike_id, label):
        """Sets the spike label for the given spike_id.  This is usually a cluster number."""
        self.cur.execute('INSERT or REPLACE into spikelabels values ({}, {});'.format(spike_id, label))

    def get_spike_label(self, spike_id):
        """Returns the spike label for the given spike_id.  This is usually a cluster number."""
        self.cur.execute('select label from spikelabels where spikeID == {};'.format(spike_id))
        result = self.cur.fetchone()
        return result[0]

    def clear_spike_labels(self):
        """Clear all of the spike labels (cluster numbers)."""
        self.cur.execute('delete from spikelabels;')


    def set_spike_coefs(self, spike_id, channel, coefs):
        """Set the spike coefficients for spike_id on channel 'channel' to be the values of coefs (a 5-element vector)."""
        for i in range(len(coefs)):
            self.cur.execute('INSERT or REPLACE INTO spikecoefs VALUES({}, {}, {}, {});'.format
                         (spike_id, channel, i, coefs[i]))

    def get_spike_coefs(self, spike_id, channel):
        """Returns a 5-vector of coefficients for the given spike_id."""
        self.cur.execute('SELECT cindex,coef FROM spikecoefs where spikeID == {} and channel == {};'.format(spike_id, channel))
        result = self.cur.fetchall()
        a = numpy.zeros(len(result))
        # (id, channel, c0, c1, c2, c3, c4) = result[channel]
        # print ( result )
        for i in range(len(result)):
            r = result[i]
            k = r[0]
            a[k] = r[1]
        return a



    # The function ndk.features.parse_coef_spec returns pairs of
    # (coef, chan).  Use this to extract a numpy array of coefficient
    # values from a database:
    def get_coef_array(self, pairs, basis=0):
        """Get the array of coefficients specified by 'pairs' for the specified basis (default 0)."""
        # Gather the data into vectors:
        pts = None
        ncols = len(pairs)
        # Count the number of spike coefficient vectors were found for the desired basis:
        self.cur.execute('select count(*) from spikecoefs join spikebasis where\
        spikecoefs.spikeID=spikebasis.spikeID and\
        spikecoefs.channel=0 and spikebasis.basisID={} and spikecoefs.cindex=0;'.format(basis))
        x = self.cur.fetchone()
        if len(x) > 0:
            npts = x[0]
            if npts > 0:
                pts = numpy.zeros((ncols, npts))
                k = 0
                # We should modify this to be able to obtain / compute
                # coefficient vector magnitudes.
                for tuple in pairs:
                    i = tuple[0]      # coef
                    chan = tuple[1]   # channel
                    self.cur.execute('select coef from spikecoefs join spikebasis where\
                    spikecoefs.cindex = {} and\
                    spikecoefs.spikeID=spikebasis.spikeID and\
                    spikecoefs.channel={} and spikebasis.basisID={}\
                    order by spikecoefs.spikeID;'.format(i, chan, basis))
                    r = self.cur.fetchall()
                    for j in range(len(r)):
                        pts[k][j] = r[j][0]
                    pts[k] -= pts[k].min()
                    if pts[k].max() > 0:
                        scale = 1.0 / pts[k].max()
                    else:
                        scale = 1.0
                    pts[k] = pts[k]*scale - 0.5
                    k += 1
                return pts

    # Now, you MUST supply a basis, although the default is basis 0:
    def get_label_array(self, basis=0):
        """Get the array of labels (sorted by spikeID) for the given basis (default 0)."""
        self.cur.execute('select label from spikelabels join spikebasis where \
        spikelabels.spikeID=spikebasis.spikeID and \
        spikebasis.basisID={} order by spikelabels.spikeID;'.format(basis))
        r = self.cur.fetchall()
        result = []
        for j in range(len(r)):
            result.append(r[j][0])
        return result



    def add_event(self, eventlabel, event_time):
        """Adds a generic event with label 'eventlabel' and timestamp 'event_time' (in samples)."""
        self.cur.execute("insert into event values ('{}', {})".format(eventlabel, event_time))

    def clear_events(self):
        """Clears all generic events."""
        self.cur.execute('delete from event;')

    def get_events(self):
        """Returns a list of tuples for behavioral or external events, where the first elements are timestamps, and the second elements are event labels."""
        self.cur.execute('SELECT * from event')
        r = self.cur.fetchall()
        events = []
        for j in range(len(r)):
            label = r[j][0]
            ts = int( r[j][1] )
            events.append((ts, label))
        return events


    def get_metadata(self, verbose=False):
        """Returns the sample rate, the number of channels, and the full name of the raw data source for this event store."""
        try:
            self.cur.execute('SELECT * from meta;')
        except:
            return None

        # Try to return the full name of the data resource:
        result = self.cur.fetchone()
        sample_rate = result[0]
        nchannels = result[1]
        return sample_rate, nchannels, ndk.ds.fullname(result[2])




    # So far, we haven't run into trouble, but throughout this file,
    # we should be mindful of the need to commit changes:
    def commit(self):
        """Commit all changes to the db object."""
        self.conn.commit()

    def close(self):
        """Close the event store, committing all changes."""
        self.conn.commit()
        self.conn.close()
        self.conn = None
        self.cur = None

    def label_histogram(self, basis=0):
        """Return a vector that contains the number of spikes in each label, for the given basis set (default 0)."""
        self.cur.execute('select max(label) from spikelabels join spikebasis where spikelabels.spikeID=spikebasis.spikeID and spikebasis.basisID={};'.format(basis))
        ml = self.cur.fetchone()
        n = ml[0] + 1
        h = []
        for i in range(n):
            self.cur.execute('select count(*) from spikelabels join spikebasis where spikelabels.spikeID=spikebasis.spikeID and spikebasis.basisID={} and spikelabels.label == {};'.format(basis,i))
            r = self.cur.fetchone()
            h.append(r[0])
        return h

    # Philosophy: Don't represent waveforms as actual samples.  Let
    # Cassandra (or equivalent low-end) mechanism do that.  Represent
    # spike events as a timestamp + coefficients for a set of basis
    # functions that allow us to reconstruct the waveform.  The
    # duration of the spike (or other event) is stored in the basis
    # table.  How far can we go with the Gaussian basis?

    ################################################################
    #
    # This file is the SOLE SOURCE for "standard" table definitions.
    # Users can create other tables at will, but this is the minimal set.

    def make_all_tables(self):
        for key in self.schema.defs:
            self.schema.create_table(self, key)

    def make_basis_table(self, verbose=False):
        """Make the basis table for this event store."""
        self.schema.create_table(self, 'basis')


    # On the one hand, it's not clear whether we really need a primary key
    # on this table.  On the other hand, there's a risk that we would
    # duplicate events:

    def make_event_table(self, verbose=False):
        """Make the generic event table for this event store.  Note that generic events are not signal-derived events.  They are usually external, e.g., behavioral events or stimulus delivery."""
        self.schema.create_table(self, 'event')


    # The spiketimes table holds pairs of Spike ID + sample number:
    
    def make_spiketimes_table(self, verbose=False):
        """Make the spiketimes table, to hold spike IDs and their timestamps (in samples)."""
        self.schema.create_table(self, 'spiketimes')

    # The spikebasis table holds pairs of Spike ID + basis ID:
    
    def make_spikebasis_table(self, verbose=False):
        """Make the spike basis table.  This associates each spike ID with a reconstruction basis set."""
        self.schema.create_table(self, 'spikebasis')

    #
    # Channels are often given names (strings).  This table uses an
    # integer channel ID.  There should be an optional table that maps
    # channel IDs to names:
    #

    def make_channel_table(self, channel_info):
        """Make the channel attribute table for this event store, specifying information line gain, bandpass, impedance, quality."""
        self.schema.create_table(self, 'channel')


    def make_meta_table(self, samprate, nchannels, fname, verbose=False):
        """Make the metadata table for this event store, specifying the sampling rate, number of channels, and full name of the raw data source."""
        ret = True
        self.sample_rate = samprate
        self.nchannels = nchannels
        self.recording = fname
        self.schema.create_table(self, 'meta')
        self.execute('SELECT * from meta;')
        result = self.fetchone()
        print ( result )
        if result == None:
            self.execute('INSERT INTO meta values ({}, {}, "{}");'.format(samprate, nchannels, fname))
            ret = True
        else:
            if result[0] != samprate or result[1] != nchannels or result[2] != fname:
                print ( "You supplied arguments that are inconsistent with this db file!" )
                print ( "db file says: {}".format(result) )
                print ( "args say: {} {} {}".format(samprate, nchannels, fname) )
                ret = False
            else:
                if verbose:
                    print ( "Valid db file for these arguments." )
                    ret = True
        return ret

    # This version of the coefficient table adds a coefficient index,
    # allowing us to represent an arbitrary number of coefficients to
    # represent an event.  We use a different name so that we can pull
    # the switch when we believe it works.
    def make_spikecoefs_table(self, verbose=False):
        """Make the spike coefficient table for this event store.  This saves
the reconstruction coefficients associated with a spike ID and one
channel of the recording.  Each coefficient also has an index that
tells us which basis component it belongs to. """
        self.schema.create_table(self, 'spikecoefs')


    # Old schoole coefs table: stores 5 coefficients in each row.   Deprecated.
    def make_spikecoefs2_table(self, verbose=False):
        """Old style of coefficient storage: Make the spike coefficient table
for this event store.  This saves the reconstruction coefficients
associated with a spike ID and one channel of the recording.  This
table currently stores 5 coefficients per spike / channel
combination."""
        self.schema.create_table(self, 'spikecoefs2')


    def make_spikelabels_table(self, verbose=False):
        """Make the spike label table for this event store.  This associates a label (usually a cluster number) with each spike ID.  Note that clusters are INTRA-BASIS ONLY.  That is, the label is qualified by the basis ID, so that cluster 0 in basis ID 0 is completely disjoint from cluster 0 in basis ID 1.  Different bases are assumed to represent disctinct phenomena, and are therefore clustered independently."""
        self.schema.create_table(self, 'spikelabels')


    ################################################################


    def make_default_labels(self):
        """Sets all spikes to have label 0."""
        self.execute('insert into spikelabels select spikeID,0 from spiketimes;')


    # We represent each spike as a vector of coefficients for a particular
    # orthonormal basis.  The basis is represented in a table 'basis' that
    # contains a basisID, a sigma, and a width in milliseconds.  The
    # latter two quantities are used to generate an orthonormal basis (see
    # ndk.features.basis) *expressed in terms of samples*.  Spikes
    # themselves need only refer to the basisID.


    def create_feature(c, name, units=None):
        """Experimental:  Creates a new feature table that can be used to store user-defined features."""
        try:
            c.execute('CREATE TABLE ' + name + ' (timestamp integer primary key, value real)')
        except s3.OperationalError as e:
            print ( e )
        finally:
            return True


    def add_value(c, name, timestamp, value):
        """Do not use."""
        c.execute('INSERT INTO '+ name + 'values ({},{});'.format(timestamp, value))



    def add_label(self, spike_id, label):
        """Set the given spike ID to have the given 'label'."""
        self.execute('replace into spikelabels values({}, {});'.format(spike_id, label))


    #
    # Old - unify this with get_coef_array below:
    def get_db_coefs(self, which, into=None):
        """Deprecated."""
        meta = self.get_metadata()
        nchannels = meta[1]
        print ( "get_db_coefs has been DEPRECATED!" )
        if into != None:
            r = into
        else:
            self.execute('select count(*) from spikecoefs where channel == 0;')
            result = self.fetchone()
            nspikes = result[0]
            r = numpy.zeros((nspikes, nchannels))
        for k in range(nchannels):
            self.execute('SELECT coef from spikecoefs where cindex={} and channel=={} order by spikeID;'.format(which,k))
            result = self.fetchall()
            for j in range(nspikes):
                r[j][k] = result[j][0]
        return r



    def get_signal(self, channel, from_time, to_time, into_array=None):
        """Reconstruct the signal on *channel* between *from_time* and *to_time* (sample numbers) and return the resulting array of values."""
        if into_array != None:
            result = into_array
            result[:] = 0
            to_time = from_time + len(result)
            n = len(result)
        else:
            n = (to_time - from_time)
            result = numpy.zeros(n)
        self.execute('select spikeID, samplenum from spiketimes where samplenum between {} and {} order by samplenum'.format(from_time, to_time))
        qr = self.fetchall()
        for i in range(len(qr)):
            spike_id = qr[i][0]
            sampnum = qr[i][1]
            w = self.get_spike_waveform(spike_id, channel)
            offset = (sampnum - from_time) 
            for k in range(0, len(w)):
                j = offset + k
                if j >= 0 and j < n:
                    result[j] = w[k]
        return result
        
    # Convention for tetrode data: If T is the tetrode number, and C
    # is the channel within the tetrode, then the es channel number is
    # N = 4*T + C.
    def fill_from_pre(self, pre):
        """Given a dict *pre* containing information derived from a .pre
file, populate the event store with relevant information from
*pre_dict*."""
        samprate = pre['rate']
        tetrodes = pre['electrodes']
        nchan = 6 * 4 # number of tetrodes * number of channels per tetrode.
        fname = 'multiple'
        self.make_meta_table(pre_dict['rate'], nchannels, fname)
        self.make_basis_table()
        self.add_basis(0, 8.0, 32)
        self.make_spikebasis_table()
        self.make_spikelabels_table()
        self.make_spiketimes_table()
        
        # Fill the channel attributes table here
        for k, entry in pre:
            if isinstance(k, int):
                spikeID = k
                spike_ts = entry[0]
                spike_tetrode = entry[1]
                waveforms = entry[2]
                self.set_spike_time(spike_id, spike_ts)

#def open_event_store(filename, create=False):
#    """The approved way of opening and returning an independent event store file."""
#    if create or os.path.isfile(filename):
#        return mysqldb(filename)


