import numpy as np
from cassandra.cluster import Cluster
from cassandra.protocol import NumpyProtocolHandler
from cassandra.query import tuple_factory, SimpleStatement
#
# At present, Cassandra service on vegas.csl.sri.com must be accessed
# by ssh tunneling.  Use this command to log into vegas and tunnel
# port 9042 to the local host:
#
# ssh vegas.csl.sri.com -L 9042:vegas.csl.sri.com:9042
#
# Remember to 'pip install cqlsh' - should this be part of setup.py?

#
# This is a template for "event level" data store objects:
#
class cdb():
    """Primary class for handling data stores via Apache Cassandra (to be implemented)."""
    def __init__(self, dataset, addr='127.0.0.1'):
        self.dataset = dataset
        # print("Creating Cassandra datastore object, descriptor = {}".format(dataset))
        self.datafile = self.dataset.path.split('/')[-1]
        self.cluster = Cluster([addr])
        self.session = self.cluster.connect('ndk')
        self.session.row_factory = tuple_factory
        self.session.client_protocol_handler=NumpyProtocolHandler
        row = self.session.execute("select samplingrate from sequentialdata limit 1")
        self.samprate = row[0]['samplingrate'][0]

    # One of the main functions for any data store:
    def get_signal(self, channel):
        """Return the entire time-series for the specified channel."""
        cnum = 'CSC{}'.format(channel)
        fname = self.dataset
        gnum = '{}_{}'.format(self.dataset.split('.')[0], cnum)
        stmt =  "SELECT value from sequentialdata where \
                 namedatafile='{}' and channelnumber='{}' \
                 and groupnumber='{}' and \
                 samplingrate = {}".format(fname, cnum, gnum, self.samprate)
        # print( query )
        rows = self.session.execute(query)
        return rows

    def get_chunk(self, channel, start, end):
        """Returns the portion of the time series for the specified channel that is between start and end inclusive (time)."""
        cnum = 'CSC{}'.format(channel+1)
        fname = self.datafile
        gnum = '{}_{}'.format(fname.split('.')[0], cnum)
        query = "SELECT value from sequentialdata where \
                 namedatafile='{}' and channelnumber='{}' \
                 and groupnumber='{}' and \
                 samplingrate = {} and \
                 samplen >= {} and samplen <= {}".format(fname, cnum, gnum, self.samprate, int(start), int(end))
        npts = int(end - start)
        q = SimpleStatement(query, fetch_size=npts)
        rows = self.session.execute(q)
        out = np.zeros(npts+1)
        k = 0
        for r in rows:
            next = r['value']
            n = len(next)
            out[k:n] = next[:]
        return out


    def filter_data(self, low_cutoff, high_cutoff):
        print('filter_data() not yet implemented for Apache Cassandra data store objects.')
        
    def put_signal(self, channel, data, start=None, end=None):
        print('put_signal() not yet implemented for Apache Cassandra data store objects.')


    def fill_event_store(self, es, sigma=0.32, window=2.0, nfuncs=5):
        print('fill_event_store not yet implemented for Apache Cassandra data store objects.')

# c = cass.cdb('CervicalVagusRecording7.smr')
# x = c.get_channel(0)

